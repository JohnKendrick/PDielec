#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA
#
r"""Layered Raman scattering calculator using GTMcore for optical field propagation.

The Raman amplitude for phonon mode m is computed as a source-overlap integral
over the Raman-active layers in a multilayer stack::

    A_m = sum_l integral E_S(z)^T R_{l,m}^lab E_L(z) dz

where E_L is the incident laser field (evaluated at the laser frequency ν_L)
and E_S is the reciprocal scattered field (launched back from the detector
channel, evaluated at the scattered frequency ν_S = ν_L − ν_m).  R_{ℓ,m}^lab
is the Raman tensor of mode m in layer ℓ rotated into the laboratory frame.

Phase 2 features (all enabled by default):

- E_S evaluated at ν_S per mode (set ``approximate_es=True`` for the Phase 1
  approximation E_S ≈ E_L).
- Forward scattering: set ``collection_side='substrate'`` to launch E_S from
  the substrate side using a reversed stack (GTMcore.System.reversed_system).
- Independent collection angle: ``collection_angle_rad`` for E_S (defaults to
  the incident angle for backscattering).
- Layer combination: set ``coherent_layers=True`` to sum amplitudes
  across Raman-active layers before squaring (GUI setting: "Coherent amplitudes").

Phase 3c Jones-vector polarisation:

- Incident and detected channels are internally represented as Jones vectors
  ``[cp, cs]`` (complex coefficients for p- and s-pol field components).  The
  string shortcuts ``incident_pol='p'`` → ``[1, 0]``, ``'s'`` → ``[0, 1]``,
  and ``detected_pol='unpolarised'`` → incoherent sum of \|A_p\|² + \|A_s\|² are
  converted automatically.  Direct Jones vector control can be added later via
  a subclass or optional parameter.

The integral is evaluated numerically using Gauss-Legendre quadrature within
each Raman-active layer.
"""

import dataclasses
import inspect
import logging
from functools import partial

import numpy as np

from PDielec.Constants import boltzmann_si, planck_si, speed_light_si
from PDielec.OpticalChannelResolver import OpticalChannelResolver, compute_modal_fields_at_points
from PDielec.PhononFinalStateResolver import (
    BULK_PHASE_MATCHED,
    COHERENT_EXTERNAL_CHANNEL,
    COHERENT_FILM,
    DISCARDED_INTERNAL_COMPONENT,
    RESOLVED_FINAL_STATE,
    PhononFinalStateResolver,
)
from PDielec.RamanAmplitudeAccumulator import RamanAmplitudeAccumulator

logger = logging.getLogger(__name__)

DEPTH_INTEGRATION_COHERENT = "Coherent amplitude"
DEPTH_INTEGRATION_INCOHERENT = "Incoherent intensity"
DEPTH_INTEGRATION_OPTIONS = (DEPTH_INTEGRATION_COHERENT, DEPTH_INTEGRATION_INCOHERENT)

# Modal-pair summation policies (only relevant when modal_pairs=True).
# GROUP_Q  : coherently sum amplitudes that share the same phonon q-vector and
#            detector channel, across all layers.  Physically recommended default.
# INCOHERENT: square each (layer, i_L, j_S, det_pol) amplitude independently —
#            diagnostic / pre-q-grouping-fix compatibility mode.
# COHERENT_ALL: sum ALL modal-pair amplitudes into a single group before squaring —
#            debug only; mixes distinct phonon-momentum final states.
@dataclasses.dataclass
class RamanContribution:
    """A single coherent-group amplitude contribution to one phonon mode.

    Returned by ``LayeredRamanCalculator.calculate_mode_intensities`` when
    ``return_contributions=True``.  Each instance represents one coherent
    group — the amplitudes inside the group are combined before squaring;
    groups are summed incoherently (as intensities).

    Attributes
    ----------
    mode_idx : int
        Index into ``raman_layers[0].phonon_frequencies_cm1``.
    frequency : float
        Phonon frequency in cm⁻¹.
    sigma : float
        Phonon linewidth in cm⁻¹.
    group_key : object
        Hashable identifier for the coherent group.  Meaning depends on
        the calculation path:

        * Modal-pairs ``GROUP_Q``       — ``(q_key_tuple, det_pol_int)``
        * Modal-pairs ``INCOHERENT``    — ``(layer_idx, i_ch, det_pol, j_ch)``
        * Modal-pairs ``COHERENT_ALL``  — the string ``"all"``
        * Standard coherent layers     — ``("all", det_pol_int)``
        * Standard incoherent layers   — ``(layer_idx, det_pol_int)``
        * Standard incoherent depth    — ``(layer_idx, "depth")``
    amplitude : complex
        Depth-integrated complex Raman amplitude for this group.
        Zero for the incoherent-depth path (use ``local_intensity`` instead).
    intensity : float
        This group's contribution to the final mode intensity ``I_m``
        (Bose-Einstein factor already applied).
    local_intensity : float
        Non-zero only for the incoherent-depth integration path, where
        ``amplitude`` is not meaningful.  Equal to the raw depth-integrated
        local intensity before the Bose factor.
    """
    mode_idx: int
    frequency: float
    sigma: float
    group_key: object
    amplitude: complex
    intensity: float
    local_intensity: float = 0.0


MODAL_PAIR_GROUP_Q      = "Group q channels"
MODAL_PAIR_INCOHERENT   = "Incoherent pairs"
MODAL_PAIR_COHERENT_ALL = "Coherent all pairs"
MODAL_PAIR_OPTIONS = (MODAL_PAIR_GROUP_Q, MODAL_PAIR_INCOHERENT, MODAL_PAIR_COHERENT_ALL)

# Modes below this frequency (cm⁻¹) are treated as acoustic and skipped.
_ACOUSTIC_THRESHOLD_CM1 = 10.0

# Target z-points per pool chunk for the parallel E-field computation.
# Larger values reduce dispatch overhead; smaller values increase parallelism.
_EFIELD_CHUNK_SIZE = 1000

# Maximum order for one Gauss-Legendre panel.  Larger requested point counts
# are implemented as a composite rule made from repeated panels of this order.
_MAX_GAUSS_LEGENDRE_PANEL_ORDER = 64

# Tolerance for grouping Berreman modes that are numerically split members of
# the same propagation-q subspace.  At exact normal incidence on a high-symmetry
# surface the two in-plane Berreman modes should be degenerate, but DFT numerical
# noise in the permittivity tensor (ε_∞ and ionic oscillator strengths) typically
# produces |Δn| ~ 1e-8 to 1e-6 even after point-group symmetrisation of ε_∞.
# Physically distinct birefringence (e.g. ZnO ordinary vs extraordinary at 10°
# incidence) gives |Δn| ~ 0.06, so a tolerance of 1e-4 safely groups the
# numerically-split degenerate modes while leaving genuine birefringence intact.
_MODAL_Q_GROUP_TOL = 1.0e-4


def _gauss_legendre_nodes_weights(n):
    """Return Gauss-Legendre quadrature nodes and weights on [-1, 1].

    For moderate point counts this is one ordinary Gauss-Legendre rule.  For
    large point counts it becomes a composite Gauss-Legendre rule with repeated
    fixed-order panels.  The composite form keeps setup cost and numerical
    conditioning practical for mm/cm-thick layers where the GUI can request
    hundreds of thousands of integration points.

    Parameters
    ----------
    n : int
        Number of quadrature points.  Must be at least 1.

    Returns
    -------
    nodes : ndarray, shape (n,)
        Gauss-Legendre abscissas on [-1, 1].
    weights : ndarray, shape (n,)
        Gauss-Legendre weights.  Their sum equals 2.0 (= length of [-1, 1]).

    """
    n = max(1, int(n))
    if n <= _MAX_GAUSS_LEGENDRE_PANEL_ORDER:
        return np.polynomial.legendre.leggauss(n)

    order = _MAX_GAUSS_LEGENDRE_PANEL_ORDER
    n_panels = int(np.ceil(n / order))
    panel_nodes, panel_weights = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(-1.0, 1.0, n_panels + 1)
    nodes = np.empty(n_panels * order, dtype=float)
    weights = np.empty(n_panels * order, dtype=float)

    for panel in range(n_panels):
        left = edges[panel]
        right = edges[panel + 1]
        half_width = 0.5 * (right - left)
        midpoint = 0.5 * (right + left)
        start = panel * order
        end = start + order
        nodes[start:end] = midpoint + half_width * panel_nodes
        weights[start:end] = half_width * panel_weights

    return nodes, weights


def _compute_efield_chunk_worker(shared, z_chunk):
    """Compute the GTM electric field for a contiguous chunk of z-positions.

    Used to parallelise the E_L (and approximate E_S) field computation across
    pool workers.  Each worker handles an independent, contiguous sub-array of
    the full GL z-grid.

    Parameters
    ----------
    shared : tuple
        ``(system, freq_cm1, angle_rad)`` — the optical system, evaluation
        frequency in cm⁻¹, and angle of incidence in radians.  The system is
        pickled once per ``functools.partial`` binding.
    z_chunk : ndarray, shape (M,)
        Sorted z-coordinates (metres) for this chunk.

    Returns
    -------
    E_chunk : ndarray, shape (6, M)
        Electric field at the requested z-positions.

    """
    system, freq_cm1, angle_rad = shared
    freq_hz = freq_cm1 * speed_light_si * 1e2
    system.initialize_sys(freq_hz)
    zeta = np.sin(angle_rad) * np.sqrt(system.superstrate.epsilon[0, 0])
    _, E_chunk, _ = system.calculate_Efield(freq_hz, zeta, z_vect=z_chunk)
    return E_chunk


def _compute_raman_mode_worker(shared, mode_args):
    """Compute the Raman intensity for a single phonon mode.

    Designed to be called via ``pool.imap(partial(_compute_raman_mode_worker, shared),
    mode_args_list)`` so that the large shared arrays (E_L, GL grid) are bound once
    to the partial function and only the small per-mode data is iterated.

    Parameters
    ----------
    shared : tuple
        ``(es_system, collection_angle_rad, z_s_arr, E_L_out, E_S_fixed,
        gl_phys_weights, gl_layer_slices, rotation_matrices,
        incident_jones, detected_jones, coherent_layers, depth_integration,
        temperature_K)``

        ``es_system`` is ``None`` when ``E_S_fixed`` is provided (approximate-ES
        or pre-computed forward-scatter field).  ``E_S_fixed`` is ``None`` when
        ``es_system`` should be used to compute E_S fresh at ``nu_S``.
    mode_args : tuple
        ``(mode_idx, nu_m, sigma, nu_S, mode_raman_tensors, mode_selected)``

        ``nu_S`` is the scattered frequency in cm⁻¹ (ignored when ``E_S_fixed``
        is not ``None``).  ``mode_raman_tensors`` is a list of (3, 3) complex
        arrays, one per Raman-active layer, for this mode.

    Returns
    -------
    tuple or None
        ``(mode_idx, nu_m, I_m, sigma)`` for an active mode, or ``None`` if
        the mode is below the acoustic threshold or the scattered frequency is
        non-positive.

    """
    (es_system, collection_angle_rad, z_s_arr,
     E_L_out, E_S_fixed,
     gl_phys_weights, gl_layer_slices,
     rotation_matrices,
     incident_jones, detected_jones,
     coherent_layers, depth_integration, temperature_K) = shared

    (mode_idx, nu_m, sigma, nu_S, mode_raman_tensors, mode_selected) = mode_args

    if not mode_selected:
        return None
    if abs(nu_m) < _ACOUSTIC_THRESHOLD_CM1:
        return None

    # Obtain E_S: either precomputed or freshly evaluated at nu_S
    if E_S_fixed is not None:
        E_S_out = E_S_fixed
    else:
        if nu_S <= 0.0:
            return None
        freq_hz = nu_S * speed_light_si * 1e2
        es_system.initialize_sys(freq_hz)
        zeta = np.sin(collection_angle_rad) * np.sqrt(es_system.superstrate.epsilon[0, 0])
        _, E_S_out, _ = es_system.calculate_Efield(freq_hz, zeta, z_vect=z_s_arr)

    cp_L, cs_L = incident_jones

    if depth_integration == DEPTH_INTEGRATION_INCOHERENT:
        I_m = 0.0
        for sl, R_crystal, G in zip(gl_layer_slices, mode_raman_tensors, rotation_matrices):
            R_lab = G @ np.asarray(R_crystal, dtype=complex) @ G.T
            E_L = cp_L * E_L_out[0:3, sl] + cs_L * E_L_out[3:6, sl]
            R_E_L = R_lab @ E_L
            w = gl_phys_weights[sl]
            if detected_jones is None:
                integrand_p = np.einsum("ij,ij->j", E_S_out[0:3, sl], R_E_L)
                integrand_s = np.einsum("ij,ij->j", E_S_out[3:6, sl], R_E_L)
                I_m += np.dot(w, np.abs(integrand_p) ** 2 + np.abs(integrand_s) ** 2)
            else:
                cp_S, cs_S = detected_jones
                E_S = cp_S * E_S_out[0:3, sl] + cs_S * E_S_out[3:6, sl]
                integrand = np.einsum("ij,ij->j", E_S, R_E_L)
                I_m += np.dot(w, np.abs(integrand) ** 2)
    elif coherent_layers:
        total_amp_p = 0.0 + 0.0j
        total_amp_s = 0.0 + 0.0j
        for sl, R_crystal, G in zip(gl_layer_slices, mode_raman_tensors, rotation_matrices):
            R_lab = G @ np.asarray(R_crystal, dtype=complex) @ G.T
            E_L = cp_L * E_L_out[0:3, sl] + cs_L * E_L_out[3:6, sl]
            R_E_L = R_lab @ E_L
            w = gl_phys_weights[sl]
            if detected_jones is None:
                total_amp_p += np.dot(w, np.einsum("ij,ij->j", E_S_out[0:3, sl], R_E_L))
                total_amp_s += np.dot(w, np.einsum("ij,ij->j", E_S_out[3:6, sl], R_E_L))
            else:
                cp_S, cs_S = detected_jones
                E_S = cp_S * E_S_out[0:3, sl] + cs_S * E_S_out[3:6, sl]
                total_amp_p += np.dot(w, np.einsum("ij,ij->j", E_S, R_E_L))
        I_m = abs(total_amp_p) ** 2 + abs(total_amp_s) ** 2 if detected_jones is None else abs(total_amp_p) ** 2
    else:
        I_m = 0.0
        for sl, R_crystal, G in zip(gl_layer_slices, mode_raman_tensors, rotation_matrices):
            R_lab = G @ np.asarray(R_crystal, dtype=complex) @ G.T
            E_L = cp_L * E_L_out[0:3, sl] + cs_L * E_L_out[3:6, sl]
            R_E_L = R_lab @ E_L
            w = gl_phys_weights[sl]
            if detected_jones is None:
                amp_p = np.dot(w, np.einsum("ij,ij->j", E_S_out[0:3, sl], R_E_L))
                amp_s = np.dot(w, np.einsum("ij,ij->j", E_S_out[3:6, sl], R_E_L))
                I_m += abs(amp_p) ** 2 + abs(amp_s) ** 2
            else:
                cp_S, cs_S = detected_jones
                E_S = cp_S * E_S_out[0:3, sl] + cs_S * E_S_out[3:6, sl]
                I_m += abs(np.dot(w, np.einsum("ij,ij->j", E_S, R_E_L))) ** 2

    I_m *= bose_factor(nu_m, temperature_K)
    return (mode_idx, nu_m, I_m, sigma)


def _modal_pair_layer_worker_data(raman_layers):
    """Return pickle-friendly Raman layer data for modal-pair workers."""
    return [
        (
            rl.layer_index,
            rl.phonon_frequencies_cm1,
            rl.raman_tensors,
            rl.rotation_matrix,
            rl.modes_selected,
        )
        for rl in raman_layers
    ]


def _raman_layer_mode_is_selected(modes_selected, mode_idx):
    """Return whether a mode is selected for a RamanLayer in worker code."""
    if modes_selected is None:
        return True
    return mode_idx < len(modes_selected) and bool(modes_selected[mode_idx])


def _compute_modal_pair_mode_worker(shared, mode_idx):
    """Compute modal-pair spectral rows for one phonon mode."""
    (
        gl_layer_slices,
        gl_phys_weights,
        channels_L,
        channels_S_base,
        channels_S_by_frequency,
        nac_cache,
        modal_pair_q_vectors,
        raman_layer_data,
        linewidths_cm1,
        modal_pair_use_nac,
        detected_pol_indices,
        detected_pol,
        modal_pair_combination,
        coherent_layers,
        depth_integration,
        approximate_es,
        resolver,
        temperature_K,
    ) = shared

    _ref_layer_index, ref_freqs, _ref_tensors, _ref_rotation, _ref_selected = raman_layer_data[0]
    fallback_nu_m = ref_freqs[mode_idx]
    sigma = linewidths_cm1[mode_idx] if mode_idx < len(linewidths_cm1) else 5.0
    use_q_resolved_nac = (
        modal_pair_use_nac is None
        or mode_idx >= len(modal_pair_use_nac)
        or bool(modal_pair_use_nac[mode_idx])
    )
    accumulator = RamanAmplitudeAccumulator()
    local_integrands = {}
    line_metadata = {}

    for layer_data, sl in zip(raman_layer_data, gl_layer_slices, strict=False):
        layer_index, _phonon_frequencies_cm1, raman_tensors, rotation_matrix, modes_selected = layer_data
        layer_channels_L = channels_L.get(layer_index, [])
        for i_channel, channel_L in enumerate(layer_channels_L):
            for det_pol_idx in detected_pol_indices:
                layer_channels_S_base = channels_S_base[det_pol_idx].get(layer_index, [])
                for j_channel, _channel_S_base in enumerate(layer_channels_S_base):
                    cache_key = (layer_index, i_channel, det_pol_idx, j_channel)
                    cache_val = nac_cache.get(cache_key) if use_q_resolved_nac else None
                    fallback_selected = _raman_layer_mode_is_selected(modes_selected, mode_idx)
                    nu_m, nac_tensors, line_sigma, mode_selected = _modal_pair_phonon_data(
                        mode_idx, fallback_nu_m, raman_tensors, sigma, fallback_selected, cache_val
                    )
                    if not mode_selected or abs(nu_m) < _ACOUSTIC_THRESHOLD_CM1:
                        continue
                    if mode_idx >= len(nac_tensors):
                        continue

                    if approximate_es:
                        channels_S = channels_S_base
                    else:
                        freq_key = round(float(nu_m), 10)
                        channels_S = channels_S_by_frequency.get(freq_key)
                        if channels_S is None:
                            continue

                    layer_channels_S = channels_S[det_pol_idx].get(layer_index, [])
                    if j_channel >= len(layer_channels_S):
                        continue
                    channel_S = layer_channels_S[j_channel]

                    R_crystal = nac_tensors[mode_idx]
                    R_lab = rotation_matrix @ R_crystal @ rotation_matrix.T

                    w = gl_phys_weights[sl]
                    R_E_L = R_lab @ channel_L["field"]
                    integrand = np.einsum("ij,ij->j", channel_S["field"], R_E_L)
                    amp_ij = np.dot(w, integrand)

                    if modal_pair_combination == MODAL_PAIR_GROUP_Q:
                        q_ph = modal_pair_q_vectors.get(cache_key, np.zeros(3, dtype=float))
                        classification = resolver.resolve_pair(
                            mode_idx,
                            q_ph,
                            is_polar=use_q_resolved_nac,
                            frequency_cm1=nu_m,
                            raman_tensor=R_lab,
                            linewidth_cm1=line_sigma,
                            selected=mode_selected,
                            reference_k=max(np.linalg.norm(resolver.q_ext), 1.0),
                        )
                        if classification.kind == DISCARDED_INTERNAL_COMPONENT:
                            continue
                        if not use_q_resolved_nac:
                            q_key = ("q_independent",)
                        elif classification.kind == RESOLVED_FINAL_STATE:
                            q_key = classification.final_state.key[-1]
                        elif classification.kind == COHERENT_EXTERNAL_CHANNEL:
                            q_key = resolver.q_class_key(q_ph)
                        else:
                            q_key = resolver.q_class_key(q_ph)
                        det_key = det_pol_idx if detected_pol == "unpolarised" else 0
                        final_state_key = (
                            "modal_pair",
                            det_key,
                            mode_idx,
                            _line_key(nu_m, line_sigma),
                            q_key,
                        )
                        group_key = final_state_key if coherent_layers else (layer_index, final_state_key)
                    elif modal_pair_combination == MODAL_PAIR_INCOHERENT:
                        if not use_q_resolved_nac:
                            # Non-polar mode: all pairs are q-independent — group coherently
                            # so the combination option does not affect non-polar intensities.
                            final_state_key = ("modal_pair_incoherent_np", mode_idx, _line_key(nu_m, line_sigma))
                            group_key = "all" if coherent_layers else (layer_index, "all")
                        else:
                            final_state_key = ("modal_pair_incoherent", cache_key, _line_key(nu_m, line_sigma))
                            group_key = cache_key
                    else:
                        final_state_key = ("modal_pair_all", mode_idx, _line_key(nu_m, line_sigma))
                        group_key = "all" if coherent_layers else (layer_index, "all")

                    line_metadata[final_state_key] = (nu_m, line_sigma)
                    if depth_integration == DEPTH_INTEGRATION_INCOHERENT:
                        local_key = (layer_index, sl.start, sl.stop)
                        integrand_key = (final_state_key, group_key, local_key)
                        if integrand_key not in local_integrands:
                            local_integrands[integrand_key] = [np.zeros_like(integrand), w]
                        local_integrands[integrand_key][0] += integrand
                    else:
                        accumulator.add_amplitude(final_state_key, amp_ij, coherence_group=group_key)

    for (final_state_key, _group_key, _local_key), (integrand, weights) in local_integrands.items():
        accumulator.add_incoherent_intensity(final_state_key, float(np.dot(weights, np.abs(integrand) ** 2)))

    line_accumulator = {}
    for final_state_key, raw_intensity in accumulator.intensities().items():
        nu_m, sigma = line_metadata[final_state_key]
        intensity = raw_intensity * bose_factor(nu_m, temperature_K)
        lkey = _line_key(nu_m, sigma)
        if lkey not in line_accumulator:
            line_accumulator[lkey] = [nu_m, intensity, sigma]
        else:
            line_accumulator[lkey][1] += intensity

    return mode_idx, list(line_accumulator.values())


def _modal_q_channels_from_fields(modal_fields, qs_by_layer, layer_indices, layer_slices, pol_idx):
    """Group modal field contributions into qz subspace channels."""
    channels_by_layer = {}
    for layer_index, sl in zip(layer_indices, layer_slices, strict=True):
        qs_arr = qs_by_layer.get(layer_index, np.zeros(4, dtype=complex))
        groups = []
        for mode_idx, qz in enumerate(qs_arr):
            qz_re = float(np.real(qz))
            group = None
            for candidate in groups:
                if abs(qz_re - candidate["qz"]) <= _MODAL_Q_GROUP_TOL:
                    group = candidate
                    break
            if group is None:
                group = {
                    "qz": qz_re,
                    "field": np.zeros((3, sl.stop - sl.start), dtype=np.complex128),
                }
                groups.append(group)
            group["field"] += modal_fields[mode_idx, pol_idx, :, sl]

        channels_by_layer[layer_index] = [
            {"qz": data["qz"], "field": data["field"]}
            for data in sorted(groups, key=lambda item: item["qz"], reverse=True)
        ]
    return channels_by_layer


def _modal_pair_stokes_channels(
    es_system,
    collection_angle_rad,
    z_s,
    layer_indices,
    gl_layer_slices,
    detected_system_layer_indices,
    detected_pol_indices,
    detected_jones,
    laser_frequency_cm1,
    mode_frequency_cm1,
):
    """Return scattered modal channels for a Stokes-shifted phonon frequency."""
    nu_S = laser_frequency_cm1 - mode_frequency_cm1
    if nu_S <= 0.0:
        logger.warning("modal_pairs: mode %.1f cm⁻¹ exceeds laser freq; skipping.", mode_frequency_cm1)
        return None

    freq_hz = nu_S * speed_light_si * 1e2
    es_system.initialize_sys(freq_hz)
    zeta = np.sin(collection_angle_rad) * np.sqrt(es_system.superstrate.epsilon[0, 0])
    modal_fields_S, qs_S_mode = compute_modal_fields_at_points(
        es_system,
        freq_hz,
        zeta,
        np.asarray(z_s, dtype=float),
        layer_indices,
        gl_layer_slices,
        system_layer_indices=detected_system_layer_indices,
    )
    channels_S = {
        pol_idx: _modal_q_channels_from_fields(
            modal_fields_S,
            qs_S_mode,
            layer_indices,
            gl_layer_slices,
            pol_idx,
        )
        for pol_idx in detected_pol_indices
    }
    if detected_jones is not None:
        for pol_idx, pol_channels in channels_S.items():
            coeff = detected_jones[pol_idx]
            for channels in pol_channels.values():
                for channel in channels:
                    channel["field"] = coeff * channel["field"]
    return channels_S


def _modal_pair_phonon_data(mode_idx, fallback_frequency, fallback_tensors,
                            fallback_sigma, fallback_selected, cache_val):
    """Return q-resolved phonon data and the selection flag for one pair."""
    if cache_val is not None and cache_val[0] is not None:
        nac_freqs = cache_val[0]
        nac_tensors = cache_val[1]
        nac_sigmas = cache_val[2] if len(cache_val) > 2 else None
        nac_selected = cache_val[3] if len(cache_val) > 3 else None
        dominant_to_by_nac = cache_val[4] if len(cache_val) > 4 else None
        nac_idx = mode_idx
        if dominant_to_by_nac is not None:
            candidates = np.where(np.asarray(dominant_to_by_nac, dtype=int) == mode_idx)[0]
            if len(candidates) == 0:
                return float(fallback_frequency), fallback_tensors, fallback_sigma, False
            if len(candidates) == 1:
                nac_idx = int(candidates[0])
            else:
                fallback_frequency = float(fallback_frequency)
                nac_idx = int(candidates[np.argmin(np.abs(np.asarray(nac_freqs)[candidates] - fallback_frequency))])
        if nac_idx < len(nac_freqs) and nac_idx < len(nac_tensors):
            selected = fallback_selected
            if nac_selected is not None and nac_idx < len(nac_selected):
                selected = bool(nac_selected[nac_idx])
            sigma = fallback_sigma
            if nac_sigmas is not None and nac_idx < len(nac_sigmas):
                sigma = float(nac_sigmas[nac_idx])
            mapped_tensors = list(fallback_tensors)
            if mode_idx < len(mapped_tensors):
                mapped_tensors[mode_idx] = nac_tensors[nac_idx]
            return float(nac_freqs[nac_idx]), mapped_tensors, sigma, selected
    return float(fallback_frequency), fallback_tensors, fallback_sigma, fallback_selected


def _line_key(frequency_cm1, sigma_cm1):
    """Hashable key for merging numerically identical spectral lines."""
    return (round(float(frequency_cm1), 10), round(float(sigma_cm1), 10))


def _resolver_coherence_regime(depth_integration):
    """Return the final-state resolver regime matching the depth treatment."""
    if depth_integration == DEPTH_INTEGRATION_INCOHERENT:
        return BULK_PHASE_MATCHED
    return COHERENT_FILM


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

class RamanLayer:
    """A Raman-active layer descriptor for LayeredRamanCalculator.

    Parameters
    ----------
    layer_index : int
        Zero-based index into ``System.layers`` (the list of finite-thickness
        layers, excluding the semi-infinite superstrate and substrate).
    phonon_frequencies_cm1 : array_like, shape (n_modes,)
        All phonon frequencies in cm⁻¹, including acoustic branches.
    raman_tensors : list of ndarray, shape (n_modes, 3, 3)
        3×3 Raman susceptibility tensors in the DFT crystal frame, one per
        mode, using the physical bulk ``R_epsilon`` convention from
        ``Raman-Theory.tex``.
    rotation_matrix : ndarray, shape (3, 3)
        Combined rotation G that maps crystal-frame vectors to the lab frame:
        ``v_lab = G @ v_crystal``.  Typically this is the product of the
        global azimuthal rotation and the HKL surface-normal rotation, i.e.
        ``G = R_z(psi_global) @ layer.euler`` where ``layer.euler`` is the
        rotation computed by ``SingleCrystalLayer.calculate_euler_matrix()``.

    Notes
    -----
    For a rank-2 tensor the lab-frame rotation is:
        ``R_lab = G @ R_crystal @ G.T``

    """

    def __init__(self, layer_index, phonon_frequencies_cm1, raman_tensors, rotation_matrix,
                 nac_function=None, modes_selected=None):
        """Initialise a RamanLayer descriptor."""
        self.layer_index = layer_index
        self.phonon_frequencies_cm1 = np.asarray(phonon_frequencies_cm1, dtype=float)
        self.raman_tensors = [np.asarray(R, dtype=complex) for R in raman_tensors]
        self.rotation_matrix = np.asarray(rotation_matrix, dtype=float)
        # Optional closure for on-demand NAC computation (Level 3 modal_pairs).
        # Contract: f(q_hat_lab: ndarray[3]) -> (nac_freqs, nac_tensors, nac_sigmas[, nac_selected])
        #           or (None, None, None[, None]) when |q_ph| is negligible.
        self.nac_function = nac_function
        self.modes_selected = None if modes_selected is None else np.asarray(modes_selected, dtype=bool)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def bose_factor(nu_cm1, temperature_K):
    """Return the Stokes Bose-Einstein thermal prefactor (n+1)/ν.

    Parameters
    ----------
    nu_cm1 : float
        Phonon frequency in cm⁻¹.  Must be positive.
    temperature_K : float
        Temperature in Kelvin.

    Returns
    -------
    float
        ``(n(ν) + 1) / ν`` where ``n(ν) = 1 / (exp(hcν / kT) - 1)``.
        At T = 0 the limit ``1/ν`` is returned.

    """
    if temperature_K <= 0.0:
        return 1.0 / nu_cm1
    hcnu_over_kT = (planck_si * speed_light_si * 100.0 * nu_cm1) / (boltzmann_si * temperature_K)
    if hcnu_over_kT > 100.0:
        # Effectively zero-temperature limit: n ≈ 0
        return 1.0 / nu_cm1
    n = 1.0 / (np.exp(hcnu_over_kT) - 1.0)
    return (n + 1.0) / nu_cm1


def lorentzian_broaden(mode_frequencies_cm1, mode_intensities, linewidths_cm1, freq_axis_cm1):
    """Broaden a set of delta-function Raman intensities with Lorentzian lineshapes.

    Parameters
    ----------
    mode_frequencies_cm1 : array_like, shape (M,)
        Mode frequencies in cm⁻¹ (only active modes, after acoustic filtering).
    mode_intensities : array_like, shape (M,)
        Integrated mode intensities before broadening.
    linewidths_cm1 : array_like, shape (M,)
        Lorentzian half-widths at half-maximum σ_m in cm⁻¹.
    freq_axis_cm1 : array_like, shape (N,)
        Output Raman-shift axis in cm⁻¹.

    Returns
    -------
    spectrum : ndarray, shape (N,)
        Broadened Raman spectrum I(Δν).

    """
    freq_axis_cm1 = np.asarray(freq_axis_cm1)
    spectrum = np.zeros(len(freq_axis_cm1))
    for nu_m, I_m, sigma in zip(mode_frequencies_cm1, mode_intensities, linewidths_cm1, strict=False):
        if sigma <= 0.0:
            continue
        spectrum += I_m * sigma / ((freq_axis_cm1 - nu_m) ** 2 + sigma ** 2)
    return spectrum


def _polarisation_to_jones(pol, name, allow_unpolarised=False):
    """Return a normalised p/s Jones vector, or None for unpolarised detection."""
    if isinstance(pol, str):
        if pol == "p":
            return np.array([1.0 + 0j, 0.0 + 0j])
        if pol == "s":
            return np.array([0.0 + 0j, 1.0 + 0j])
        if allow_unpolarised and pol == "unpolarised":
            return None
        allowed = "'p' or 's'" if not allow_unpolarised else "'p', 's', or 'unpolarised'"
        raise ValueError(f"{name} must be {allowed}, got '{pol}'")

    arr = np.asarray(pol, dtype=complex)
    if arr.shape != (2,):
        raise ValueError(f"{name} Jones vector must have shape (2,), got {arr.shape}")
    norm = np.linalg.norm(arr)
    if norm == 0.0:
        raise ValueError(f"{name} Jones vector must be non-zero")
    return arr / norm


# ---------------------------------------------------------------------------
# Main calculator class
# ---------------------------------------------------------------------------

class LayeredRamanCalculator:
    r"""Calculate the Raman spectrum of a multilayer system via depth integration.

    The calculation proceeds as follows:

    1. Build Gauss-Legendre quadrature grids within each Raman-active layer.
    2. Compute the electric field E_L at the laser frequency for both p and s
       incidence in a single call to ``system.calculate_Efield``.
    3. For each phonon mode m (skipping acoustic modes):

       a. Compute E_S at the scattered frequency ν_S = ν_L − ν_m on the
          appropriate system (original for backscattering, reversed for
          forward scattering).  When ``approximate_es=True`` the Phase 1
          approximation E_S ≈ E_L is used instead (one field call total).
       b. For each Raman-active layer ℓ, evaluate the integrand

              a_{ℓ,m}(z_j) = E_S(z_j)^T R_{ℓ,m}^lab E_L(z_j)

          and integrate with the GL weights to obtain A_{ℓ,m}.
       c. Combine layer contributions:
          - ``coherent_layers=False`` (default): I_m ∝ Σ_ℓ \|A_{ℓ,m}\|²
          - ``coherent_layers=True``:            I_m ∝ \|Σ_ℓ A_{ℓ,m}\|²
       d. Apply the Bose-Einstein prefactor: I_m ← (n(ν_m)+1)/ν_m × I_m.

    4. Broaden all mode intensities with Lorentzian lineshapes and sum.

    Parameters
    ----------
    system : GTMcore.System or GTMcore.ScatteringMatrixSystem
        Fully constructed GTM multilayer system.  The layers must already have
        their Euler orientations set (via ``layer.set_euler`` or equivalent)
        before this class is used.
    raman_layers : list of RamanLayer
        Descriptors for Raman-active layers, in any order.  Each
        ``RamanLayer.layer_index`` must be a valid index into ``system.layers``.
    laser_frequency_cm1 : float
        Laser frequency in cm⁻¹.
    incident_angle_rad : float
        Angle of incidence in radians, measured from the surface normal.
    incident_pol : {'p', 's'} or array_like, shape (2,)
        Incident polarisation channel, or a Jones vector ``[cp, cs]``.
    detected_pol : {'p', 's', 'unpolarised'} or array_like, shape (2,)
        Detected polarisation channel.  'unpolarised' sums the p and s
        intensities incoherently (no analyser on the detector side).  A Jones
        vector ``[cp, cs]`` describes a linear analyser channel.
    temperature_K : float
        Sample temperature in Kelvin for the Bose-Einstein factor.
    linewidths_cm1 : array_like, shape (n_modes,)
        Lorentzian half-widths in cm⁻¹, one entry per phonon mode (including
        acoustic modes; those are skipped internally).
    n_gauss : int, optional
        Number of Gauss-Legendre quadrature points per Raman-active layer.
        Default is 21.
    collection_side : {'superstrate', 'substrate'}, optional
        Which side the detector is on.  ``'superstrate'`` (default) gives
        backscattering geometry; ``'substrate'`` gives forward scattering and
        causes E_S to be computed on the reversed stack.
    collection_angle_rad : float or None, optional
        Collection (detector) angle in radians from the surface normal.
        Defaults to ``incident_angle_rad`` (appropriate for backscattering).
    coherent_layers : bool, optional
        If ``True``, amplitudes are summed across Raman-active layers before
        squaring (coherent combination).  Default is ``False`` (incoherent:
        intensities are summed).
    approximate_es : bool, optional
        If ``True``, use the Phase 1 approximation E_S ≈ E_L (both fields at
        the laser frequency).  Default is ``False``.
    depth_integration : {'Coherent amplitude', 'Incoherent intensity'}, optional
        How to combine the source along the depth of each Raman-active layer.
        ``'Coherent amplitude'`` integrates the complex amplitude before
        squaring and is appropriate for thin coherent films.  ``'Incoherent
        intensity'`` integrates the local intensity and is more stable for
        thick or bulk samples where long-range phase coherence is not physical.
    Notes
    -----
    The E-field returned by ``calculate_Efield`` has shape ``(6, N)`` where
    rows 0–2 are the x, y, z components for p-polarised incidence and rows 3–5
    are the components for s-polarised incidence.  The laboratory frame uses the
    GTMcore convention: x in the plane of incidence, y along the slit (s-pol
    direction), z normal to the layers (depth axis).

    The Raman tensor contraction uses an ordinary (non-conjugate) transpose for
    E_S, consistent with the optical-reciprocity formulation of layered Raman:

        A = Σ_j w_j  E_S(z_j)^T  R_lab  E_L(z_j)

    For forward scattering (``collection_side='substrate'``), the z-coordinates
    for E_S are remapped as ``z_rev = total_thickness − z_orig`` before being
    passed to ``calculate_Efield`` on the reversed system; the quadrature
    weights are unchanged.

    """

    def __init__(
        self,
        system,
        raman_layers,
        laser_frequency_cm1,
        incident_angle_rad,
        incident_pol,
        detected_pol,
        temperature_K,
        linewidths_cm1,
        n_gauss=21,
        collection_side="superstrate",
        collection_angle_rad=None,
        coherent_layers=False,
        approximate_es=False,
        depth_integration=DEPTH_INTEGRATION_COHERENT,
        modal_pairs=False,
        modal_pair_combination=MODAL_PAIR_GROUP_Q,
        modal_pair_use_nac=None,
        q_tol_deg=0.0,
        modes_selected=None,
    ):
        """Initialise LayeredRamanCalculator with system, layers and calculation parameters."""
        incident_jones = _polarisation_to_jones(incident_pol, "incident_pol", allow_unpolarised=False)
        detected_jones = _polarisation_to_jones(detected_pol, "detected_pol", allow_unpolarised=True)
        if collection_side not in ("superstrate", "substrate"):
            raise ValueError(f"collection_side must be 'superstrate' or 'substrate', got '{collection_side}'")
        if depth_integration not in DEPTH_INTEGRATION_OPTIONS:
            raise ValueError(
                f"depth_integration must be one of {DEPTH_INTEGRATION_OPTIONS}, got {depth_integration}"
            )
        if modal_pair_combination not in MODAL_PAIR_OPTIONS:
            raise ValueError(
                f"modal_pair_combination must be one of {MODAL_PAIR_OPTIONS}, got '{modal_pair_combination}'"
            )

        self.system = system
        self.raman_layers = list(raman_layers)
        self.laser_frequency_cm1 = float(laser_frequency_cm1)
        self.incident_angle_rad = float(incident_angle_rad)
        self.incident_pol = incident_pol
        self.detected_pol = detected_pol
        self.temperature_K = float(temperature_K)
        self.linewidths_cm1 = np.asarray(linewidths_cm1, dtype=float)
        self.n_gauss = max(1, int(n_gauss))
        self.collection_side = collection_side
        self.collection_angle_rad = float(collection_angle_rad) if collection_angle_rad is not None else float(incident_angle_rad)
        self.coherent_layers = bool(coherent_layers)
        self.approximate_es = bool(approximate_es)
        self.depth_integration = depth_integration
        self.modal_pairs = bool(modal_pairs)
        self.modal_pair_combination = modal_pair_combination
        self.modal_pair_use_nac = None if modal_pair_use_nac is None else np.asarray(modal_pair_use_nac, dtype=bool)
        self.q_tol_deg = float(q_tol_deg)
        self.modes_selected = None if modes_selected is None else np.asarray(modes_selected, dtype=bool)
        self._nac_cache = {}
        self._modal_pair_q_vectors = {}

        # Phase 3c: internal Jones vectors for incident and detected channels.
        # _incident_jones : (cp, cs) complex pair — determines the linear combination
        #     of p- and s-pol incident fields used in the overlap integral.
        # _detected_jones : (cp, cs) or None — None signals incoherent (no-analyser)
        #     detection: both p and s amplitudes are squared and summed.
        self._incident_jones = incident_jones
        self._detected_jones = detected_jones

        # Gauss-Legendre nodes and weights on [-1, 1].  Large point counts use
        # a composite fixed-order Gauss rule so thick layers remain practical.
        self._quad_nodes, self._quad_weights = _gauss_legendre_nodes_weights(self.n_gauss)
        self.n_gauss = len(self._quad_nodes)

        # Built by _build_gl_grid()
        self._gl_z = None               # ndarray: concatenated z points (m) for E_L
        self._gl_phys_weights = None    # ndarray: corresponding physical integration weights (m)
        self._gl_layer_slices = None    # list of slice: one per RamanLayer

        self._build_gl_grid()

    def _mode_is_selected(self, mode_idx):
        """Return True when a phonon mode is enabled in the caller's mode mask."""
        if self.modes_selected is None or mode_idx >= len(self.modes_selected):
            return True
        return bool(self.modes_selected[mode_idx])

    def _layer_mode_is_selected(self, raman_layer, mode_idx):
        """Return True when a mode is enabled for a specific Raman layer."""
        layer_selected = getattr(raman_layer, "modes_selected", None)
        if layer_selected is not None and mode_idx < len(layer_selected):
            return bool(layer_selected[mode_idx])
        return self._mode_is_selected(mode_idx)

    def _any_layer_mode_is_selected(self, mode_idx):
        """Return True when at least one Raman-active layer enables this mode."""
        return any(self._layer_mode_is_selected(rl, mode_idx) for rl in self.raman_layers)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_gl_grid(self):
        """Map Gauss-Legendre quadrature points onto each Raman-active layer interval.

        Uses ``system.get_layers_boundaries()`` to find each layer's z extent.
        Boundaries are in metres with z = 0 at the superstrate/first-layer
        interface (the convention used by ``calculate_Efield``).

        Populates
        ---------
        _gl_z : ndarray
            Concatenated quadrature z coordinates for all Raman-active layers, in metres
            (the unit ``calculate_Efield`` expects).
        _gl_phys_weights : ndarray
            Corresponding integration weights in Ångström (quadrature weight × half-thickness × 1e10).
            The Å unit cancels the Å² in the Raman tensor more cleanly than metres,
            keeping intensity magnitudes in a readable range.
        _gl_layer_slices : list of slice
            One slice per ``RamanLayer`` indexing into ``_gl_z`` / ``_gl_phys_weights``.
        """
        boundaries = self.system.get_layers_boundaries()
        # boundaries[0]   = -superstrate.thick
        # boundaries[1]   = 0  (incident surface)
        # boundaries[k+1] = top of system.layers[k] (k = 0 … n_layers-1)
        # boundaries[k+2] = bottom of system.layers[k]

        z_parts = []
        w_parts = []
        slices = []
        offset = 0
        for rl in self.raman_layers:
            idx = rl.layer_index
            z_start = boundaries[idx + 1]
            z_end   = boundaries[idx + 2]
            half    = 0.5 * (z_end - z_start)
            mid     = 0.5 * (z_start + z_end)
            z_j = mid + half * self._quad_nodes             # metres (for calculate_Efield)
            w_j = half * self._quad_weights * 1.0e10       # Å (for the intensity integral)
            z_parts.append(z_j)
            w_parts.append(w_j)
            slices.append(slice(offset, offset + self.n_gauss))
            offset += self.n_gauss

        self._gl_z = np.concatenate(z_parts)
        self._gl_phys_weights = np.concatenate(w_parts)
        self._gl_layer_slices = slices

    def _get_field_at_gl_points(self, freq_cm1, system, angle_rad, z_arr):
        """Compute the electric field at the given z positions.

        Parameters
        ----------
        freq_cm1 : float
            Frequency in cm⁻¹ at which to evaluate the field.
        system : GTMcore.System
            The optical system to use for the field calculation.
        angle_rad : float
            Angle of incidence on ``system``'s superstrate, in radians.
        z_arr : ndarray
            z-coordinates (metres) at which to evaluate the field.  Must be
            in the coordinate system of ``system`` (z = 0 at superstrate/
            first-layer interface).

        Returns
        -------
        E_out : ndarray, shape (6, len(z_arr))
            Electric field: rows 0–2 for p-pol incidence, rows 3–5 for s-pol.

        """
        freq_hz = freq_cm1 * speed_light_si * 1e2
        system.initialize_sys(freq_hz)
        zeta = np.sin(angle_rad) * np.sqrt(system.superstrate.epsilon[0, 0])
        _, E_out, _ = system.calculate_Efield(freq_hz, zeta, z_vect=z_arr)
        return E_out

    @staticmethod
    def _rotate_raman_tensor(R_crystal, G):
        """Rotate a Raman tensor from crystal frame to lab frame.

        Parameters
        ----------
        R_crystal : ndarray, shape (3, 3)
            Raman tensor in the DFT crystal frame.
        G : ndarray, shape (3, 3)
            Rotation matrix satisfying ``v_lab = G @ v_crystal``.

        Returns
        -------
        ndarray, shape (3, 3)
            Raman tensor in the lab frame: ``G @ R_crystal @ G.T``.

        """
        return G @ R_crystal @ G.T

    def _layer_amplitude(self, E_L_out, E_S_out, R_lab, sl):
        """Compute the Raman amplitude for one mode in one layer.

        Evaluates the numerical quadrature:

            A = Σ_j w_j  E_S(z_j)^T  R_lab  E_L(z_j)

        Parameters
        ----------
        E_L_out : ndarray, shape (6, N_total)
            Incident laser field array from ``calculate_Efield`` at ν_L.
        E_S_out : ndarray, shape (6, N_total)
            Reciprocal scattered field array from ``calculate_Efield`` at ν_S.
            May be equal to ``E_L_out`` when ``approximate_es=True`` in
            backscattering geometry.
        R_lab : ndarray, shape (3, 3)
            Raman tensor in the lab frame.
        sl : slice
            Slice selecting the GL points for this layer.

        Returns
        -------
        amplitude_p : complex
            Amplitude for the p-detected channel.
        amplitude_s : complex or None
            Amplitude for the s-detected channel when ``detected_pol ==
            'unpolarised'``; ``None`` otherwise.

        """
        # Incident field: Jones vector combination of p and s incidence.
        # E_L_out rows 0:3 = p-pol incidence, rows 3:6 = s-pol incidence.
        cp_L, cs_L = self._incident_jones
        E_L = cp_L * E_L_out[0:3, sl] + cs_L * E_L_out[3:6, sl]  # (3, n_gauss)

        # R_lab @ E_L  →  shape (3, n_gauss)
        R_E_L = R_lab @ E_L

        # Physical quadrature weights for this layer
        w = self._gl_phys_weights[sl]  # shape (n_gauss,)

        if self._detected_jones is None:
            # Unpolarised (no analyser): incoherent sum of p and s detected channels.
            E_S_p = E_S_out[0:3, sl]
            E_S_s = E_S_out[3:6, sl]
            integrand_p = np.einsum("ij,ij->j", E_S_p, R_E_L)  # shape (n_gauss,)
            integrand_s = np.einsum("ij,ij->j", E_S_s, R_E_L)
            return np.dot(w, integrand_p), np.dot(w, integrand_s)
        # Polarised detection: Jones vector combination of p and s detected fields.
        cp_S, cs_S = self._detected_jones
        E_S = cp_S * E_S_out[0:3, sl] + cs_S * E_S_out[3:6, sl]
        integrand = np.einsum("ij,ij->j", E_S, R_E_L)
        return np.dot(w, integrand), None

    def _layer_depth_intensity(self, E_L_out, E_S_out, R_lab, sl):
        """Compute the incoherent depth-integrated Raman intensity for one layer."""
        cp_L, cs_L = self._incident_jones
        E_L = cp_L * E_L_out[0:3, sl] + cs_L * E_L_out[3:6, sl]
        R_E_L = R_lab @ E_L
        w = self._gl_phys_weights[sl]

        if self._detected_jones is None:
            integrand_p = np.einsum("ij,ij->j", E_S_out[0:3, sl], R_E_L)
            integrand_s = np.einsum("ij,ij->j", E_S_out[3:6, sl], R_E_L)
            return np.dot(w, np.abs(integrand_p) ** 2 + np.abs(integrand_s) ** 2)

        cp_S, cs_S = self._detected_jones
        E_S = cp_S * E_S_out[0:3, sl] + cs_S * E_S_out[3:6, sl]
        integrand = np.einsum("ij,ij->j", E_S, R_E_L)
        return np.dot(w, np.abs(integrand) ** 2)

    def _get_modal_fields_at_gl_points(self, freq_cm1, system, angle_rad, z_arr, system_layer_indices=None):
        """Compute per-Berreman-mode electric field contributions at GL quadrature points.

        For each Raman-active layer and each of its GL z-points, the field
        contribution from each of the 4 Berreman modes is evaluated separately.
        The total field at any z-point equals the sum of all 4 modal contributions,
        which reproduces ``calculate_Efield`` exactly.

        Parameters
        ----------
        freq_cm1 : float
            Frequency in cm⁻¹ at which to solve the GTM system.
        system : GTMcore.System
            Optical system (TransferMatrixSystem or ScatteringMatrixSystem).
        angle_rad : float
            Angle of incidence on ``system``'s superstrate, in radians.
        z_arr : ndarray
            GL z-coordinates in metres (in system's coordinate frame).
        system_layer_indices : list of int or None
            Optional mapping from Raman-layer physical indices to indices in
            ``system.layers``.  This is required for reversed systems used by
            forward-scattering reciprocal fields.

        Returns
        -------
        modal_fields : ndarray, shape (4, 2, 3, N_gl_total)
            ``modal_fields[mode_n, pol_idx, xyz, z_pt]``
            mode_n: 0..3 (trans-p, trans-s, refl-p, refl-s).
            pol_idx: 0 = p-pol incidence on system, 1 = s-pol incidence.
            xyz: Ex, Ey, Ez components.
        qs_by_layer : dict of {int: ndarray(4,)}
            Maps physical ``rl.layer_index`` to layer.qs (4 kz eigenvalues)
            from the corresponding optical-system layer.
        zeta : complex
            In-plane wavevector used for this system/frequency/angle.
        """
        freq_hz = freq_cm1 * speed_light_si * 1e2
        system.initialize_sys(freq_hz)
        zeta = np.sin(angle_rad) * np.sqrt(system.superstrate.epsilon[0, 0])
        layer_indices = [rl.layer_index for rl in self.raman_layers]
        modal_fields, qs_by_layer = compute_modal_fields_at_points(
            system,
            freq_hz,
            zeta,
            np.asarray(z_arr, dtype=float),
            layer_indices,
            self._gl_layer_slices,
            system_layer_indices=system_layer_indices,
        )
        return modal_fields, qs_by_layer, zeta

    def _call_get_modal_fields_at_gl_points(self, freq_cm1, system, angle_rad, z_arr, system_layer_indices=None):
        """Call the modal-field helper while tolerating legacy test doubles."""
        helper = self._get_modal_fields_at_gl_points
        helper_identity = getattr(helper, "__func__", helper)
        cached_identity = getattr(self, "_modal_fields_helper_identity", None)
        if cached_identity is not helper_identity:
            self._modal_fields_helper_identity = helper_identity
            self._modal_fields_helper_accepts_system_layer_indices = (
                "system_layer_indices" in inspect.signature(helper).parameters
            )
        if self._modal_fields_helper_accepts_system_layer_indices:
            return helper(freq_cm1, system, angle_rad, z_arr, system_layer_indices=system_layer_indices)
        return helper(freq_cm1, system, angle_rad, z_arr)

    def _get_modal_q_channels(self, modal_fields, qs_by_layer, pol_idx):
        """Group modal fields into propagation-q subspaces for one polarisation.

        Berreman eigenvectors inside an exactly or nearly degenerate qz
        subspace are not unique.  At normal incidence in a uniaxial layer, for
        example, the two forward in-plane optical modes have the same qz and
        may rotate within that subspace as the crystal azimuth changes.  The
        physical p/s field contribution is the coherent sum over the whole qz
        subspace, not an individual eigenvector selected by its modal index.
        """
        layer_indices = [rl.layer_index for rl in self.raman_layers]
        return _modal_q_channels_from_fields(
            modal_fields,
            qs_by_layer,
            layer_indices,
            self._gl_layer_slices,
            pol_idx,
        )

    def _calculate_mode_intensities_modal_pairs(self, progress_callback=None, pool=None, return_contributions=False):
        """Level 3 modal-pairs Raman intensity calculation.

        For each phonon mode, sums over incident and detected propagation-q
        channels.  Each channel is the coherent field contribution from all
        Berreman eigenmodes in the same qz subspace for the selected
        polarisation.  This avoids basis dependence when the optical eigenmodes
        are degenerate or nearly degenerate, as at normal incidence in a
        uniaxial layer.

        The phonon wavevector for each channel pair is q_ph = k_L − k_S,
        giving a distinct NAC-corrected phonon frequency and Raman tensor per
        pair.  NAC results are cached per layer/channel pair and reused across
        all phonon modes.  The returned arrays are q-resolved: one input
        ``mode_idx`` can produce multiple spectral rows when different q-groups
        have different NAC frequencies.  Groups with numerically identical
        frequencies and linewidths are merged as incoherent intensities.

        NAC closures are evaluated only while building the serial cache.  When a
        worker pool is supplied, the per-mode field and intensity calculations
        are dispatched with picklable arrays and GTM systems.

        Returns
        -------
        active_frequencies : ndarray
        active_intensities : ndarray
        active_linewidths : ndarray
        """
        # --- Set up ES system and z-array (mirrors calculate_mode_intensities) ---
        if self.collection_side == "substrate":
            total_thick = sum(layer.thick for layer in self.system.layers)
            z_s = total_thick - self._gl_z
            es_system = self.system.reversed_system()
            detected_system_layer_indices = [len(self.system.layers) - 1 - rl.layer_index for rl in self.raman_layers]
        else:
            z_s = self._gl_z
            es_system = self.system
            detected_system_layer_indices = [rl.layer_index for rl in self.raman_layers]

        ref_layer = self.raman_layers[0]
        n_modes = len(ref_layer.phonon_frequencies_cm1)

        incident_pol_indices = [idx for idx, value in enumerate(self._incident_jones) if abs(value) > 1.0e-14]
        detected_pol_indices = (
            [0, 1]
            if self._detected_jones is None
            else [idx for idx, value in enumerate(self._detected_jones) if abs(value) > 1.0e-14]
        )

        q_ext = OpticalChannelResolver().external_momentum_transfer(
            self.system,
            self.laser_frequency_cm1,
            self.incident_angle_rad,
            es_system,
            self.laser_frequency_cm1,
            self.collection_angle_rad,
            self.collection_side,
        )
        use_external_q_for_modal_nac = self.depth_integration == DEPTH_INTEGRATION_INCOHERENT

        # Laser modal fields at laser frequency + incident angle
        modal_fields_L, qs_L_dict, zeta_L = self._call_get_modal_fields_at_gl_points(
            self.laser_frequency_cm1, self.system, self.incident_angle_rad, self._gl_z
        )
        zeta_L_re = float(np.real(zeta_L))
        channels_L_by_pol = {
            pol_idx: self._get_modal_q_channels(modal_fields_L, qs_L_dict, pol_idx)
            for pol_idx in incident_pol_indices
        }
        channels_L = {}
        for pol_idx, pol_channels in channels_L_by_pol.items():
            coeff = self._incident_jones[pol_idx]
            for layer_index, channels in pol_channels.items():
                merged = channels_L.setdefault(layer_index, [])
                for channel in channels:
                    merged.append({"qz": channel["qz"], "field": coeff * channel["field"]})

        # ES modal fields at laser frequency (used for q_ph / NAC computation,
        # and also for the field integral when approximate_es=True)
        modal_fields_S_base, qs_S_dict, zeta_S = self._call_get_modal_fields_at_gl_points(
            self.laser_frequency_cm1,
            es_system,
            self.collection_angle_rad,
            z_s,
            system_layer_indices=detected_system_layer_indices,
        )
        zeta_S_re = float(np.real(zeta_S))
        channels_S_base = {
            pol_idx: self._get_modal_q_channels(modal_fields_S_base, qs_S_dict, pol_idx)
            for pol_idx in detected_pol_indices
        }
        if self._detected_jones is not None:
            for pol_idx, pol_channels in channels_S_base.items():
                coeff = self._detected_jones[pol_idx]
                for channels in pol_channels.values():
                    for channel in channels:
                        channel["field"] = coeff * channel["field"]

        # Build NAC cache upfront (one call per q-channel pair per layer,
        # independent of phonon mode).  Keys use channel ordinal positions rather
        # than qz values so the cache remains valid when E_S is recomputed at
        # the Stokes frequency and qz shifts slightly.
        self._nac_cache = {}
        self._modal_pair_q_vectors = {}
        for rl in self.raman_layers:
            layer_channels_L = channels_L.get(rl.layer_index, [])
            for i_channel, channel_L in enumerate(layer_channels_L):
                for det_pol_idx in detected_pol_indices:
                    layer_channels_S = channels_S_base[det_pol_idx].get(rl.layer_index, [])
                    for j_channel, channel_S in enumerate(layer_channels_S):
                        cache_key = (rl.layer_index, i_channel, det_pol_idx, j_channel)

                        # q_ph = k_laser − k_scattered.  Channel qz values are
                        # coherent sums over degenerate Berreman eigenmodes.
                        q_pair = np.array([
                            zeta_L_re - zeta_S_re,
                            0.0,
                            channel_L["qz"] - channel_S["qz"],
                        ])
                        q_ph = q_ext if use_external_q_for_modal_nac else q_pair
                        self._modal_pair_q_vectors[cache_key] = q_ph

                        q_ph_norm = np.linalg.norm(q_ph)
                        if q_ph_norm < 1e-8 or rl.nac_function is None:
                            self._nac_cache[cache_key] = None  # use TO baseline
                        else:
                            q_hat_lab = q_ph / q_ph_norm
                            self._nac_cache[cache_key] = rl.nac_function(q_hat_lab)

        modal_pair_use_nac = self.modal_pair_use_nac
        resolver = PhononFinalStateResolver(
            coherence_regime=_resolver_coherence_regime(self.depth_integration),
            q_ext=q_ext,
            angular_tolerance_deg=max(float(self.q_tol_deg), 1.0e-10),
            q_zero_tol=1.0e-3,
        )

        active_freqs = []
        active_intensities = []
        active_sigmas = []
        contributions: list = []

        if pool is not None and not return_contributions:
            layer_indices = [rl.layer_index for rl in self.raman_layers]
            raman_layer_data = _modal_pair_layer_worker_data(self.raman_layers)
            channels_S_by_frequency = {}
            if not self.approximate_es:
                required_stokes_channels = {}
                for mode_idx in range(n_modes):
                    fallback_nu_m = ref_layer.phonon_frequencies_cm1[mode_idx]
                    sigma = self.linewidths_cm1[mode_idx] if mode_idx < len(self.linewidths_cm1) else 5.0
                    use_q_resolved_nac = (
                        modal_pair_use_nac is None
                        or mode_idx >= len(modal_pair_use_nac)
                        or bool(modal_pair_use_nac[mode_idx])
                    )
                    for layer_data in raman_layer_data:
                        layer_index, _freqs, raman_tensors, _rotation_matrix, modes_selected = layer_data
                        layer_channels_L = channels_L.get(layer_index, [])
                        for i_channel, _channel_L in enumerate(layer_channels_L):
                            for det_pol_idx in detected_pol_indices:
                                layer_channels_S_base = channels_S_base[det_pol_idx].get(layer_index, [])
                                for j_channel, _channel_S_base in enumerate(layer_channels_S_base):
                                    cache_key = (layer_index, i_channel, det_pol_idx, j_channel)
                                    cache_val = self._nac_cache.get(cache_key) if use_q_resolved_nac else None
                                    fallback_selected = _raman_layer_mode_is_selected(modes_selected, mode_idx)
                                    nu_m, nac_tensors, _line_sigma, mode_selected = _modal_pair_phonon_data(
                                        mode_idx,
                                        fallback_nu_m,
                                        raman_tensors,
                                        sigma,
                                        fallback_selected,
                                        cache_val,
                                    )
                                    if not mode_selected or abs(nu_m) < _ACOUSTIC_THRESHOLD_CM1:
                                        continue
                                    if mode_idx >= len(nac_tensors):
                                        continue
                                    required_stokes_channels.setdefault(round(float(nu_m), 10), nu_m)
                for freq_key, nu_m in required_stokes_channels.items():
                    channels_S = _modal_pair_stokes_channels(
                        es_system,
                        self.collection_angle_rad,
                        z_s,
                        layer_indices,
                        self._gl_layer_slices,
                        detected_system_layer_indices,
                        detected_pol_indices,
                        self._detected_jones,
                        self.laser_frequency_cm1,
                        nu_m,
                    )
                    if channels_S is not None:
                        channels_S_by_frequency[freq_key] = channels_S

            shared = (
                self._gl_layer_slices,
                self._gl_phys_weights,
                channels_L,
                channels_S_base,
                channels_S_by_frequency,
                self._nac_cache,
                self._modal_pair_q_vectors,
                raman_layer_data,
                self.linewidths_cm1,
                modal_pair_use_nac,
                detected_pol_indices,
                self.detected_pol,
                self.modal_pair_combination,
                self.coherent_layers,
                self.depth_integration,
                self.approximate_es,
                resolver,
                self.temperature_K,
            )
            worker_fn = partial(_compute_modal_pair_mode_worker, shared)
            for _mode_idx, rows in pool.imap(worker_fn, range(n_modes), chunksize=1):
                if progress_callback is not None:
                    progress_callback()
                for nu_m, intensity, sigma in rows:
                    active_freqs.append(nu_m)
                    active_intensities.append(intensity)
                    active_sigmas.append(sigma)
            return (
                np.array(active_freqs),
                np.array(active_intensities),
                np.array(active_sigmas),
            )

        for mode_idx in range(n_modes):
            if progress_callback is not None:
                progress_callback()

            fallback_nu_m = ref_layer.phonon_frequencies_cm1[mode_idx]
            sigma = self.linewidths_cm1[mode_idx] if mode_idx < len(self.linewidths_cm1) else 5.0
            use_q_resolved_nac = (
                modal_pair_use_nac is None
                or mode_idx >= len(modal_pair_use_nac)
                or bool(modal_pair_use_nac[mode_idx])
            )
            accumulator = RamanAmplitudeAccumulator()
            local_integrands = {}
            line_metadata = {}
            group_metadata = {}
            channels_s_cache = {}

            for rl, sl in zip(self.raman_layers, self._gl_layer_slices):
                layer_channels_L = channels_L.get(rl.layer_index, [])
                for i_channel, channel_L in enumerate(layer_channels_L):
                    for det_pol_idx in detected_pol_indices:
                        layer_channels_S_base = channels_S_base[det_pol_idx].get(rl.layer_index, [])
                        for j_channel, _channel_S_base in enumerate(layer_channels_S_base):
                            cache_key = (rl.layer_index, i_channel, det_pol_idx, j_channel)
                            cache_val = self._nac_cache.get(cache_key) if use_q_resolved_nac else None
                            fallback_selected = self._layer_mode_is_selected(rl, mode_idx)
                            nu_m, nac_tensors, line_sigma, mode_selected = _modal_pair_phonon_data(
                                mode_idx, fallback_nu_m, rl.raman_tensors, sigma, fallback_selected, cache_val
                            )
                            if not mode_selected:
                                continue

                            if abs(nu_m) < _ACOUSTIC_THRESHOLD_CM1:
                                continue

                            if mode_idx >= len(nac_tensors):
                                continue

                            if self.approximate_es:
                                channels_S = channels_S_base
                            else:
                                nu_S = self.laser_frequency_cm1 - nu_m
                                if nu_S <= 0.0:
                                    logger.warning("modal_pairs: mode %.1f cm⁻¹ exceeds laser freq; skipping.", nu_m)
                                    continue
                                freq_key = round(float(nu_m), 10)
                                if freq_key not in channels_s_cache:
                                    modal_fields_S, qs_S_mode, _ = self._call_get_modal_fields_at_gl_points(
                                        nu_S,
                                        es_system,
                                        self.collection_angle_rad,
                                        z_s,
                                        system_layer_indices=detected_system_layer_indices,
                                    )
                                    channels_s_cache[freq_key] = {
                                        pol_idx: self._get_modal_q_channels(modal_fields_S, qs_S_mode, pol_idx)
                                        for pol_idx in detected_pol_indices
                                    }
                                    if self._detected_jones is not None:
                                        for _pol, _pol_channels in channels_s_cache[freq_key].items():
                                            _coeff = self._detected_jones[_pol]
                                            for _channels in _pol_channels.values():
                                                for _channel in _channels:
                                                    _channel["field"] = _coeff * _channel["field"]
                                    for _rl in self.raman_layers:
                                        for _pol in detected_pol_indices:
                                            _n_base = len(channels_S_base[_pol].get(_rl.layer_index, []))
                                            _n_mode = len(channels_s_cache[freq_key][_pol].get(_rl.layer_index, []))
                                            if _n_base != _n_mode:
                                                logger.warning(
                                                    "modal_pairs approximate_es=False: layer %d pol %d has "
                                                    "%d qz channels at nu_S=%.1f cm-1 but the NAC cache was "
                                                    "built with %d channels at the laser frequency; NAC "
                                                    "tensors may be misassigned for this mode.",
                                                    _rl.layer_index, _pol, _n_mode, nu_S, _n_base,
                                                )
                                channels_S = channels_s_cache[freq_key]

                            layer_channels_S = channels_S[det_pol_idx].get(rl.layer_index, [])
                            if j_channel >= len(layer_channels_S):
                                continue
                            channel_S = layer_channels_S[j_channel]

                            R_crystal = nac_tensors[mode_idx]
                            R_lab = self._rotate_raman_tensor(R_crystal, rl.rotation_matrix)

                            w = self._gl_phys_weights[sl]
                            R_E_L = R_lab @ channel_L["field"]
                            integrand = np.einsum("ij,ij->j", channel_S["field"], R_E_L)
                            amp_ij = np.dot(w, integrand)

                            if self.modal_pair_combination == MODAL_PAIR_GROUP_Q:
                                q_ph = self._modal_pair_q_vectors.get(cache_key, np.zeros(3, dtype=float))
                                classification = resolver.resolve_pair(
                                    mode_idx,
                                    q_ph,
                                    is_polar=use_q_resolved_nac,
                                    frequency_cm1=nu_m,
                                    raman_tensor=R_lab,
                                    linewidth_cm1=line_sigma,
                                    selected=mode_selected,
                                    reference_k=max(np.linalg.norm(q_ext), 1.0),
                                )
                                if classification.kind == DISCARDED_INTERNAL_COMPONENT:
                                    continue
                                if not use_q_resolved_nac:
                                    q_key = ("q_independent",)
                                elif classification.kind == RESOLVED_FINAL_STATE:
                                    q_key = classification.final_state.key[-1]
                                elif classification.kind == COHERENT_EXTERNAL_CHANNEL:
                                    q_key = resolver.q_class_key(q_ph)
                                else:
                                    q_key = resolver.q_class_key(q_ph)
                                det_key = det_pol_idx if self.detected_pol == "unpolarised" else 0
                                final_state_key = (
                                    "modal_pair",
                                    det_key,
                                    mode_idx,
                                    _line_key(nu_m, line_sigma),
                                    q_key,
                                )
                                group_key = final_state_key if self.coherent_layers else (rl.layer_index, final_state_key)
                            elif self.modal_pair_combination == MODAL_PAIR_INCOHERENT:
                                if not use_q_resolved_nac:
                                    # Non-polar mode: all pairs are q-independent — group coherently
                                    # so the combination option does not affect non-polar intensities.
                                    final_state_key = ("modal_pair_incoherent_np", mode_idx, _line_key(nu_m, line_sigma))
                                    group_key = "all" if self.coherent_layers else (rl.layer_index, "all")
                                else:
                                    final_state_key = ("modal_pair_incoherent", cache_key, _line_key(nu_m, line_sigma))
                                    group_key = cache_key  # already per-layer via layer_index in cache_key
                            else:  # MODAL_PAIR_COHERENT_ALL
                                final_state_key = ("modal_pair_all", mode_idx, _line_key(nu_m, line_sigma))
                                group_key = "all" if self.coherent_layers else (rl.layer_index, "all")

                            local_key = (rl.layer_index, sl.start, sl.stop)
                            line_metadata[final_state_key] = (nu_m, line_sigma)
                            group_metadata[(final_state_key, group_key)] = (nu_m, line_sigma)
                            if self.depth_integration == DEPTH_INTEGRATION_INCOHERENT:
                                integrand_key = (final_state_key, group_key, local_key)
                                if integrand_key not in local_integrands:
                                    local_integrands[integrand_key] = [np.zeros_like(integrand), w]
                                local_integrands[integrand_key][0] += integrand
                            else:
                                accumulator.add_amplitude(final_state_key, amp_ij, coherence_group=group_key)

            for (final_state_key, _group_key, _local_key), (integrand, weights) in local_integrands.items():
                accumulator.add_incoherent_intensity(final_state_key, float(np.dot(weights, np.abs(integrand) ** 2)))

            line_accumulator = {}
            for final_state_key, raw_intensity in accumulator.intensities().items():
                nu_m, sigma = line_metadata[final_state_key]
                intensity = raw_intensity * bose_factor(nu_m, self.temperature_K)
                lkey = _line_key(nu_m, sigma)
                if lkey not in line_accumulator:
                    line_accumulator[lkey] = [nu_m, intensity, sigma]
                else:
                    line_accumulator[lkey][1] += intensity

            if return_contributions:
                for (final_state_key, group_key), amplitude in accumulator.grouped_amplitudes().items():
                    nu_m, sigma = group_metadata[(final_state_key, group_key)]
                    intensity = abs(amplitude) ** 2 * bose_factor(nu_m, self.temperature_K)
                    contributions.append(RamanContribution(
                        mode_idx=mode_idx, frequency=nu_m, sigma=sigma,
                        group_key=group_key, amplitude=amplitude, intensity=intensity,
                    ))
                for (final_state_key, group_key, _local_key), (integrand, weights) in local_integrands.items():
                    nu_m, sigma = group_metadata[(final_state_key, group_key)]
                    local_intensity = float(np.dot(weights, np.abs(integrand) ** 2))
                    contributions.append(RamanContribution(
                        mode_idx=mode_idx, frequency=nu_m, sigma=sigma,
                        group_key=group_key, amplitude=0.0 + 0.0j,
                        intensity=local_intensity * bose_factor(nu_m, self.temperature_K),
                        local_intensity=local_intensity,
                    ))

            for nu_m, intensity, sigma in line_accumulator.values():
                active_freqs.append(nu_m)
                active_intensities.append(intensity)
                active_sigmas.append(sigma)

        result = (
            np.array(active_freqs),
            np.array(active_intensities),
            np.array(active_sigmas),
        )
        return (*result, contributions) if return_contributions else result

    def resolve_optical_channels(self, scattered_frequency_cm1=None):
        """Return Phase-3 optical channel diagnostics for the configured stack.

        This parallel diagnostic path exposes the incident laser field,
        reciprocal detector field, q-resolved Berreman channel groups, and
        external momentum transfer used by the final-state refactor.  It does
        not alter the existing standard or modal-pairs Raman calculations.
        """
        if self.collection_side == "substrate":
            total_thick = sum(layer.thick for layer in self.system.layers)
            z_s = total_thick - self._gl_z
            es_system = self.system.reversed_system()
        else:
            z_s = self._gl_z
            es_system = self.system

        nu_s = self.laser_frequency_cm1 if scattered_frequency_cm1 is None else float(scattered_frequency_cm1)
        incident_pol_indices = [idx for idx, value in enumerate(self._incident_jones) if abs(value) > 1.0e-14]
        detected_pol_indices = (
            [0, 1]
            if self._detected_jones is None
            else [idx for idx, value in enumerate(self._detected_jones) if abs(value) > 1.0e-14]
        )
        layer_indices = [rl.layer_index for rl in self.raman_layers]
        if self.collection_side == "substrate":
            detected_system_layer_indices = [len(self.system.layers) - 1 - layer_index for layer_index in layer_indices]
        else:
            detected_system_layer_indices = layer_indices
        resolver = OpticalChannelResolver(optical_subspace_tolerance=_MODAL_Q_GROUP_TOL)

        incident = resolver.solve(
            self.system,
            self.laser_frequency_cm1,
            self.incident_angle_rad,
            self._gl_z,
            layer_indices,
            self._gl_layer_slices,
            incident_pol_indices,
        )
        detected = resolver.solve(
            es_system,
            nu_s,
            self.collection_angle_rad,
            z_s,
            layer_indices,
            self._gl_layer_slices,
            detected_pol_indices,
            system_layer_indices=detected_system_layer_indices,
        )
        q_ext = resolver.external_momentum_transfer(
            self.system,
            self.laser_frequency_cm1,
            self.incident_angle_rad,
            es_system,
            nu_s,
            self.collection_angle_rad,
            self.collection_side,
        )
        return {
            "incident": incident,
            "detected": detected,
            "q_ext": q_ext,
            "z_detected": z_s,
            "collection_system": es_system,
        }

    def calculate_phase3_modal_pair_diagnostic_intensities(self):
        """Run a Phase-3 q-independent modal diagnostic path.

        This uses ``OpticalChannelResolver``, ``PhononFinalStateResolver`` and
        ``RamanAmplitudeAccumulator`` together.  It is intentionally diagnostic:
        it covers q-independent/TO tensor accumulation and is meant to compare
        against simple existing modal-pairs cases before Phase 4 replacement.
        """
        if not self.approximate_es:
            msg = (
                "Phase-3 modal-pair diagnostics currently require approximate_es=True; "
                "per-mode Stokes-frequency reciprocal channels will be added before Phase 4."
            )
            raise NotImplementedError(msg)
        diagnostics = self.resolve_optical_channels(scattered_frequency_cm1=self.laser_frequency_cm1)
        incident = diagnostics["incident"]
        detected = diagnostics["detected"]
        phase_resolver = PhononFinalStateResolver()
        active_freqs = []
        active_intensities = []
        active_sigmas = []

        ref_layer = self.raman_layers[0]
        n_modes = len(ref_layer.phonon_frequencies_cm1)
        incident_pol_indices = [idx for idx, value in enumerate(self._incident_jones) if abs(value) > 1.0e-14]
        detected_pol_indices = (
            [0, 1]
            if self._detected_jones is None
            else [idx for idx, value in enumerate(self._detected_jones) if abs(value) > 1.0e-14]
        )
        if self.modal_pair_combination == MODAL_PAIR_INCOHERENT:
            self._raise_if_phase3_incoherent_has_degenerate_channels(incident, detected, incident_pol_indices, detected_pol_indices)

        for mode_idx in range(n_modes):
            nu_m = float(ref_layer.phonon_frequencies_cm1[mode_idx])
            if abs(nu_m) < _ACOUSTIC_THRESHOLD_CM1 or not self._any_layer_mode_is_selected(mode_idx):
                continue
            sigma = float(self.linewidths_cm1[mode_idx]) if mode_idx < len(self.linewidths_cm1) else 5.0
            accumulator = RamanAmplitudeAccumulator()

            for rl, layer_slice in zip(self.raman_layers, self._gl_layer_slices, strict=True):
                if not self._layer_mode_is_selected(rl, mode_idx):
                    continue
                r_lab = self._rotate_raman_tensor(rl.raman_tensors[mode_idx], rl.rotation_matrix)
                weights = self._gl_phys_weights[layer_slice]
                layer_channels_l = []
                for pol_idx in incident_pol_indices:
                    coeff = self._incident_jones[pol_idx]
                    for channel_idx, channel in enumerate(incident.channels_by_pol.get(pol_idx, {}).get(rl.layer_index, [])):
                        layer_channels_l.append((pol_idx, channel_idx, coeff, channel))

                layer_channels_s = []
                for pol_idx in detected_pol_indices:
                    coeff = 1.0 + 0.0j if self._detected_jones is None else self._detected_jones[pol_idx]
                    for channel_idx, channel in enumerate(detected.channels_by_pol.get(pol_idx, {}).get(rl.layer_index, [])):
                        layer_channels_s.append((pol_idx, channel_idx, coeff, channel))

                for _pol_l, channel_l_idx, coeff_l, channel_l in layer_channels_l:
                    for det_pol_idx, channel_s_idx, coeff_s, channel_s in layer_channels_s:
                        q_ph = np.array([
                            np.real(incident.zeta - detected.zeta),
                            0.0,
                            np.real(channel_l.qz - channel_s.qz),
                        ])
                        classification = phase_resolver.resolve_pair(
                            mode_idx,
                            q_ph,
                            is_polar=False,
                            frequency_cm1=nu_m,
                            raman_tensor=r_lab,
                            linewidth_cm1=sigma,
                            selected=True,
                        )
                        r_e_l = r_lab @ (coeff_l * channel_l.field[:, layer_slice])
                        integrand = np.einsum("ij,ij->j", coeff_s * channel_s.field[:, layer_slice], r_e_l)
                        amplitude = np.dot(weights, integrand)
                        if self.modal_pair_combination == MODAL_PAIR_INCOHERENT:
                            group = (rl.layer_index, channel_l_idx, det_pol_idx, channel_s_idx)
                        elif self.coherent_layers:
                            group = ("all", det_pol_idx)
                        else:
                            group = (rl.layer_index, det_pol_idx)
                        accumulator.add_classified(classification, amplitude, coherence_group=group)

            intensity = accumulator.total_intensity() * bose_factor(nu_m, self.temperature_K)
            if intensity > 0.0:
                active_freqs.append(nu_m)
                active_intensities.append(intensity)
                active_sigmas.append(sigma)

        return np.array(active_freqs), np.array(active_intensities), np.array(active_sigmas)

    def _raise_if_phase3_incoherent_has_degenerate_channels(
        self,
        incident,
        detected,
        incident_pol_indices,
        detected_pol_indices,
    ):
        """Reject incoherent diagnostics when q-channels contain modal subspaces."""
        for rl in self.raman_layers:
            for pol_idx in incident_pol_indices:
                for channel in incident.channels_by_pol.get(pol_idx, {}).get(rl.layer_index, []):
                    if len(channel.mode_indices) > 1:
                        msg = (
                            "Phase-3 MODAL_PAIR_INCOHERENT diagnostics require one Berreman mode per optical "
                            "channel; degenerate subspace channels cannot reproduce raw per-mode incoherent sums."
                        )
                        raise NotImplementedError(msg)
            for pol_idx in detected_pol_indices:
                for channel in detected.channels_by_pol.get(pol_idx, {}).get(rl.layer_index, []):
                    if len(channel.mode_indices) > 1:
                        msg = (
                            "Phase-3 MODAL_PAIR_INCOHERENT diagnostics require one Berreman mode per optical "
                            "channel; degenerate subspace channels cannot reproduce raw per-mode incoherent sums."
                        )
                        raise NotImplementedError(msg)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def calculate_mode_intensities(self, progress_callback=None, pool=None, return_contributions=False):
        """Compute the per-mode Raman intensities (before broadening).

        Modes below ``_ACOUSTIC_THRESHOLD_CM1`` are excluded.  Layer
        contributions are combined incoherently or coherently depending on
        ``self.coherent_layers``.

        Parameters
        ----------
        progress_callback : callable, optional
            If supplied, called once per phonon-mode iteration (including modes
            that are skipped as acoustic), so the caller can drive a progress
            bar.  The callback takes no arguments.
        pool : multiprocessing.Pool or multiprocessing.dummy.Pool, optional
            Worker pool for parallel execution.  When provided, each phonon mode
            (i.e. each E_S field evaluation at its own scattered frequency) is
            dispatched to a pool worker.  E_L is computed once on the calling
            process and broadcast to all workers via ``functools.partial``.
            When ``None`` (default), the calculation runs sequentially.

        return_contributions : bool, optional
            When True, return a fourth element — a list of
            :class:`RamanContribution` objects exposing the per-group complex
            amplitudes before squaring.  Only available in the sequential path
            (``pool`` must be ``None``).

        Returns
        -------
        active_frequencies : ndarray, shape (M,)
            Frequencies of the active (non-acoustic) spectral contributions in
            cm⁻¹.  In ``modal_pairs`` mode, one original phonon index may
            produce multiple q-resolved entries.
        active_intensities : ndarray, shape (M,)
            Corresponding Raman intensities I_m (real, non-negative).
        active_linewidths : ndarray, shape (M,)
            Linewidths of the active modes in cm⁻¹.
        contributions : list of RamanContribution
            Only present when ``return_contributions=True``.  One entry per
            (mode, coherent-group) pair, ordered by mode then group.

        Notes
        -----
        This method assumes all ``RamanLayer`` objects in ``self.raman_layers``
        have the same set of phonon modes (same DFT calculation).  If layers
        have different materials this method must be extended.

        When ``approximate_es=False`` (default), ``calculate_Efield`` is called
        once per active mode for E_S at ν_S = ν_L − ν_m, plus once for E_L.
        For systems with many modes this can be slow; set ``approximate_es=True``
        to recover the single-field-call behaviour at the cost of accuracy.

        """
        if return_contributions and pool is not None:
            raise ValueError(
                "return_contributions=True is only supported in the sequential path "
                "(pool must be None)."
            )

        if not self.raman_layers:
            logger.warning("calculate_mode_intensities: no Raman-active layers defined")
            empty = (np.array([]), np.array([]), np.array([]))
            return (*empty, []) if return_contributions else empty

        # Level 3 dispatch: build NAC cache serially, then use the pool for
        # per-mode field/intensity work when one is available.
        if self.modal_pairs:
            return self._calculate_mode_intensities_modal_pairs(
                progress_callback, pool=pool, return_contributions=return_contributions
            )

        # --- Set up the system and z-array for E_S (common to both paths) ---
        if self.collection_side == "substrate":
            # Forward scattering: launch E_S from the substrate side
            total_thick = sum(layer.thick for layer in self.system.layers)
            z_s = total_thick - self._gl_z   # z in the reversed system
            es_system = self.system.reversed_system()
        else:
            # Backscattering: same system, same z-array
            z_s = self._gl_z
            es_system = self.system

        # Use frequencies from the first RamanLayer (all layers share the same DFT modes)
        ref_layer = self.raman_layers[0]
        n_modes = len(ref_layer.phonon_frequencies_cm1)

        # Precompute per-layer rotation matrices (shared across modes)
        rotation_matrices = [rl.rotation_matrix for rl in self.raman_layers]

        if pool is not None:
            # Parallel path: E_L (and E_S_fixed when needed) are computed
            # inside _calculate_modes_parallel using the pool.
            return self._calculate_modes_parallel(
                pool, ref_layer, n_modes,
                es_system, z_s,
                rotation_matrices,
                progress_callback,
            )

        # Sequential path -------------------------------------------------
        # Compute E_L once on the calling process.
        E_L_out = self._get_field_at_gl_points(
            self.laser_frequency_cm1, self.system, self.incident_angle_rad, self._gl_z
        )

        # If approximate_es, compute E_S once at the laser frequency
        E_S_out_fixed = None
        if self.approximate_es:
            if self.collection_side == "substrate":
                E_S_out_fixed = self._get_field_at_gl_points(
                    self.laser_frequency_cm1, es_system, self.collection_angle_rad, z_s
                )
            else:
                # Pure Phase-1 approximation: E_S = E_L (no extra call)
                E_S_out_fixed = E_L_out

        active_freqs = []
        active_intensities = []
        active_sigmas = []
        contributions: list = []

        for mode_idx in range(n_modes):
            if progress_callback is not None:
                progress_callback()

            if not self._any_layer_mode_is_selected(mode_idx):
                continue

            nu_m = ref_layer.phonon_frequencies_cm1[mode_idx]

            # Skip acoustic modes
            if abs(nu_m) < _ACOUSTIC_THRESHOLD_CM1:
                continue

            sigma = self.linewidths_cm1[mode_idx] if mode_idx < len(self.linewidths_cm1) else 5.0

            # Compute or reuse E_S
            if self.approximate_es:
                E_S_out = E_S_out_fixed
            else:
                nu_S = self.laser_frequency_cm1 - nu_m
                if nu_S <= 0.0:
                    logger.warning(
                        "Mode at %.1f cm⁻¹ exceeds laser frequency; skipping.", nu_m
                    )
                    continue
                E_S_out = self._get_field_at_gl_points(
                    nu_S, es_system, self.collection_angle_rad, z_s
                )

            # Accumulate contributions from all Raman-active layers.
            # group_amps maps group_key → complex amplitude; groups are summed
            # incoherently (as |amp|²).  For the incoherent-depth path, local
            # intensities (not amplitudes) are collected in local_intens instead.
            group_amps: dict = {}
            local_intens: dict = {}

            if self.depth_integration == DEPTH_INTEGRATION_INCOHERENT:
                for rl, sl in zip(self.raman_layers, self._gl_layer_slices, strict=True):
                    if not self._layer_mode_is_selected(rl, mode_idx):
                        continue
                    R_crystal = rl.raman_tensors[mode_idx]
                    R_lab = self._rotate_raman_tensor(R_crystal, rl.rotation_matrix)
                    local_intens[(rl.layer_index, "depth")] = (
                        self._layer_depth_intensity(E_L_out, E_S_out, R_lab, sl)
                    )
                I_m = sum(local_intens.values())
            else:
                # layer_key "all" → coherent across layers; int → incoherent per layer
                for rl, sl in zip(self.raman_layers, self._gl_layer_slices, strict=True):
                    if not self._layer_mode_is_selected(rl, mode_idx):
                        continue
                    R_crystal = rl.raman_tensors[mode_idx]
                    R_lab = self._rotate_raman_tensor(R_crystal, rl.rotation_matrix)
                    amp_p, amp_s = self._layer_amplitude(E_L_out, E_S_out, R_lab, sl)
                    lkey = "all" if self.coherent_layers else rl.layer_index
                    group_amps[(lkey, 0)] = group_amps.get((lkey, 0), 0.0 + 0.0j) + amp_p
                    if self._detected_jones is None:
                        group_amps[(lkey, 1)] = group_amps.get((lkey, 1), 0.0 + 0.0j) + amp_s
                I_m = sum(abs(amp) ** 2 for amp in group_amps.values())

            # Apply Bose-Einstein thermal prefactor
            bose = bose_factor(nu_m, self.temperature_K)
            I_m *= bose

            active_freqs.append(nu_m)
            active_intensities.append(I_m)
            active_sigmas.append(sigma)

            if return_contributions:
                if local_intens:
                    for gkey, local_I in local_intens.items():
                        contributions.append(RamanContribution(
                            mode_idx=mode_idx, frequency=nu_m, sigma=sigma,
                            group_key=gkey, amplitude=0.0 + 0.0j,
                            intensity=local_I * bose, local_intensity=local_I,
                        ))
                else:
                    for gkey, amp in group_amps.items():
                        contributions.append(RamanContribution(
                            mode_idx=mode_idx, frequency=nu_m, sigma=sigma,
                            group_key=gkey, amplitude=amp,
                            intensity=abs(amp) ** 2 * bose,
                        ))

        result = (
            np.array(active_freqs),
            np.array(active_intensities),
            np.array(active_sigmas),
        )
        return (*result, contributions) if return_contributions else result

    def _get_field_parallel(self, pool, freq_cm1, system, angle_rad, z_arr):
        """Compute the GTM E-field at *z_arr* in parallel using *pool*.

        Splits ``z_arr`` into contiguous chunks of at most ``_EFIELD_CHUNK_SIZE``
        z-positions, dispatches each chunk to a pool worker, then reassembles
        the full ``(6, N)`` field array.  This avoids the long single-threaded
        Python loop inside ``calculate_Efield`` by distributing the z-iteration
        across pool workers.

        Parameters
        ----------
        pool : multiprocessing.Pool or compatible
            Worker pool used to distribute field calculations over chunks of
            z-coordinates.
        freq_cm1 : float
            Evaluation frequency in cm⁻¹.
        system : GTMcore.System
            Optical system (copied to each worker via pickling).
        angle_rad : float
            Angle of incidence in radians.
        z_arr : ndarray, shape (N,)
            Sorted z-coordinates in metres.

        Returns
        -------
        E_out : ndarray, shape (6, N)

        """
        n_chunks = max(1, len(z_arr) // _EFIELD_CHUNK_SIZE)
        chunks = np.array_split(z_arr, n_chunks)
        shared_ef = (system, freq_cm1, angle_rad)
        worker_fn = partial(_compute_efield_chunk_worker, shared_ef)
        E_parts = list(pool.imap(worker_fn, chunks, chunksize=1))
        return np.concatenate(E_parts, axis=1)

    def _calculate_modes_parallel(
        self, pool, ref_layer, n_modes,
        es_system, z_s,
        rotation_matrices,
        progress_callback,
    ):
        """Parallel implementation of mode-intensity calculation via pool.imap.

        E_L is computed in parallel by splitting the GL z-grid into chunks and
        dispatching each chunk to a pool worker.  Each phonon-mode E_S call is
        also dispatched to a pool worker.  This eliminates both the pre-pool
        single-threaded bottleneck (E_L chunk computation is parallelised)
        and the per-mode E_S bottleneck.

        Parameters
        ----------
        pool : multiprocessing.Pool or compatible
            Active worker pool.
        ref_layer : RamanLayer
            First Raman-active layer (source of phonon frequencies).
        n_modes : int
            Total number of phonon modes (including acoustic).
        es_system : GTMcore.System
            System used to evaluate E_S (original or reversed).
        z_s : ndarray
            z-coordinates (m) for the E_S field evaluation.
        rotation_matrices : list of ndarray
            Crystal-to-lab rotation matrices, one per RamanLayer.
        progress_callback : callable or None
            Called once per mode (in result-collection order) if not ``None``.

        Returns
        -------
        tuple of ndarray
            ``(active_frequencies, active_intensities, active_linewidths)``

        """
        # --- E_L: compute in parallel by splitting z-array into chunks ----------
        E_L_out = self._get_field_parallel(
            pool, self.laser_frequency_cm1, self.system, self.incident_angle_rad,
            self._gl_z,
        )

        # --- E_S_out_fixed: computed in parallel when approximate_es=True -------
        if self.approximate_es:
            if self.collection_side == "substrate":
                E_S_out_fixed = self._get_field_parallel(
                    pool, self.laser_frequency_cm1, es_system, self.collection_angle_rad,
                    z_s,
                )
            else:
                # Pure Phase-1 approximation: E_S = E_L
                E_S_out_fixed = E_L_out
        else:
            E_S_out_fixed = None

        # Shared data bound via partial — pickled once per pool worker process.
        # When E_S_out_fixed is already available, es_system is never used in
        # the worker, so pass None to avoid unnecessary serialisation.
        shared = (
            None if E_S_out_fixed is not None else es_system,
            self.collection_angle_rad,
            z_s,
            E_L_out,
            E_S_out_fixed,
            self._gl_phys_weights,
            self._gl_layer_slices,
            rotation_matrices,
            self._incident_jones,
            self._detected_jones,
            self.coherent_layers,
            self.depth_integration,
            self.temperature_K,
        )
        worker_fn = partial(_compute_raman_mode_worker, shared)

        # Per-mode args: only small data iterated per task
        mode_args_list = []
        for mode_idx in range(n_modes):
            nu_m = float(ref_layer.phonon_frequencies_cm1[mode_idx])
            sigma = float(self.linewidths_cm1[mode_idx]) if mode_idx < len(self.linewidths_cm1) else 5.0
            nu_S = self.laser_frequency_cm1 - nu_m if E_S_out_fixed is None else 0.0
            # Raman tensors for this mode across all layers
            mode_raman_tensors = [
                rl.raman_tensors[mode_idx]
                if self._layer_mode_is_selected(rl, mode_idx)
                else np.zeros_like(rl.raman_tensors[mode_idx])
                for rl in self.raman_layers
            ]
            mode_args_list.append((
                mode_idx, nu_m, sigma, nu_S, mode_raman_tensors,
                self._any_layer_mode_is_selected(mode_idx),
            ))

        active_freqs = []
        active_intensities = []
        active_sigmas = []

        for result in pool.imap(worker_fn, mode_args_list, chunksize=1):
            if progress_callback is not None:
                progress_callback()
            if result is None:
                continue
            _mode_idx, nu_m, I_m, sigma = result
            active_freqs.append(nu_m)
            active_intensities.append(I_m)
            active_sigmas.append(sigma)

        return (
            np.array(active_freqs),
            np.array(active_intensities),
            np.array(active_sigmas),
        )

    def calculate_spectrum(self, freq_axis_cm1, progress_callback=None):
        """Compute the broadened Raman spectrum on the given frequency axis.

        Parameters
        ----------
        freq_axis_cm1 : array_like, shape (N,)
            Raman-shift axis in cm⁻¹ on which to evaluate the spectrum.
        progress_callback : callable, optional
            Forwarded to :meth:`calculate_mode_intensities` so a caller can
            drive a progress bar.  Called once per phonon mode.

        Returns
        -------
        spectrum : ndarray, shape (N,)
            Raman intensity I(Δν) on the requested axis.

        """
        freq_axis_cm1 = np.asarray(freq_axis_cm1)
        active_freqs, active_intensities, active_sigmas = self.calculate_mode_intensities(
            progress_callback=progress_callback
        )
        if len(active_freqs) == 0:
            return np.zeros(len(freq_axis_cm1))
        return lorentzian_broaden(active_freqs, active_intensities, active_sigmas, freq_axis_cm1)
