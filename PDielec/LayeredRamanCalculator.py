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
import logging
from functools import partial

import numpy as np

from PDielec.Constants import boltzmann_si, planck_si, speed_light_si

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
# the same propagation-q subspace.  Normal-incidence uniaxial modes can differ
# by ~1e-10 in qz from eigensolver round-off, while physically distinct forward
# and backward subspaces differ by O(1).
_MODAL_Q_GROUP_TOL = 1.0e-8


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


def _compute_modal_fields_and_channels(freq_cm1, system, angle_rad, z_arr, layer_indices, layer_slices, pol_indices):
    """Compute modal fields and qz channels for a set of polarisations."""
    freq_hz = freq_cm1 * speed_light_si * 1e2
    system.initialize_sys(freq_hz)
    zeta = np.sin(angle_rad) * np.sqrt(system.superstrate.epsilon[0, 0])

    modal_amps, _zn = system.calculate_modal_amplitudes(freq_hz, zeta)
    boundaries = system.get_layers_boundaries()

    n_total = len(z_arr)
    modal_fields = np.zeros((4, 2, 3, n_total), dtype=np.complex128)
    qs_by_layer = {}

    for layer_index, sl in zip(layer_indices, layer_slices, strict=True):
        if layer_index not in modal_amps or layer_index >= len(system.layers):
            continue
        layer = system.layers[layer_index]
        z_reference = boundaries[layer_index + 2]
        amps = modal_amps[layer_index]
        qs_by_layer[layer_index] = layer.qs.copy()

        z_pts = z_arr[sl]
        for pt_offset, z_j in enumerate(z_pts):
            global_idx = sl.start + pt_offset
            dKiz = np.array([
                np.exp(layer.propagation_exponents[n] * (z_j - z_reference) / layer.thick)
                for n in range(4)
            ], dtype=np.complex128)

            for n in range(4):
                Eprop_p = dKiz[n] * amps[n]
                Eprop_s = dKiz[n] * amps[n + 4]
                modal_fields[n, 0, :, global_idx] = Eprop_p * layer.gamma[n, :]
                modal_fields[n, 1, :, global_idx] = Eprop_s * layer.gamma[n, :]

    channels = {
        pol_idx: _modal_q_channels_from_fields(modal_fields, qs_by_layer, layer_indices, layer_slices, pol_idx)
        for pol_idx in pol_indices
    }
    return modal_fields, qs_by_layer, zeta, channels


def _modal_pair_phonon_data(mode_idx, fallback_frequency, fallback_tensors,
                            fallback_sigma, fallback_selected, cache_val):
    """Return q-resolved phonon data and the selection flag for one pair."""
    if cache_val is not None and cache_val[0] is not None:
        nac_freqs = cache_val[0]
        nac_tensors = cache_val[1]
        nac_sigmas = cache_val[2] if len(cache_val) > 2 else None
        nac_selected = cache_val[3] if len(cache_val) > 3 else None
        if mode_idx < len(nac_freqs) and mode_idx < len(nac_tensors):
            selected = fallback_selected
            if nac_selected is not None and mode_idx < len(nac_selected):
                selected = bool(nac_selected[mode_idx])
            sigma = fallback_sigma
            if nac_sigmas is not None and mode_idx < len(nac_sigmas):
                sigma = float(nac_sigmas[mode_idx])
            return float(nac_freqs[mode_idx]), nac_tensors, sigma, selected
    return float(fallback_frequency), fallback_tensors, fallback_sigma, fallback_selected


def _line_key(frequency_cm1, sigma_cm1):
    """Hashable key for merging numerically identical spectral lines."""
    return (round(float(frequency_cm1), 10), round(float(sigma_cm1), 10))


def _compute_modal_pair_mode_worker(shared, mode_args):
    """Compute one modal-pair Raman mode using picklable cached NAC data."""
    (
        es_system,
        collection_angle_rad,
        z_s,
        channels_L,
        channels_S_base,
        detected_pol_indices,
        layer_indices,
        gl_layer_slices,
        gl_phys_weights,
        rotation_matrices,
        raman_tensors_by_layer,
        modes_selected_by_layer,
        nac_cache,
        modal_pair_q_keys,
        modal_pair_combination,
        approximate_es,
        temperature_K,
    ) = shared

    mode_idx, fallback_nu_m, sigma, laser_frequency_cm1, fallback_mode_selected = mode_args

    # Accumulate amplitudes into groups according to the summation policy.
    # The group key determines which amplitudes are coherently combined
    # before squaring.
    pair_groups = {}  # group_key -> {"amp", "frequency", "sigma"}
    channels_s_cache = {}

    for layer_pos, (layer_index, sl) in enumerate(zip(layer_indices, gl_layer_slices, strict=True)):
        layer_channels_L = channels_L.get(layer_index, [])
        for i_channel, channel_L in enumerate(layer_channels_L):
            for det_pol_idx in detected_pol_indices:
                layer_channels_S_base = channels_S_base[det_pol_idx].get(layer_index, [])
                for j_channel, _channel_S_base in enumerate(layer_channels_S_base):
                    cache_key = (layer_index, i_channel, det_pol_idx, j_channel)
                    cache_val = nac_cache.get(cache_key)

                    fallback_selected = fallback_mode_selected
                    layer_selected = modes_selected_by_layer[layer_pos]
                    if layer_selected is not None and mode_idx < len(layer_selected):
                        fallback_selected = bool(layer_selected[mode_idx])

                    nu_m, nac_tensors, line_sigma, mode_selected = _modal_pair_phonon_data(
                        mode_idx,
                        fallback_nu_m,
                        raman_tensors_by_layer[layer_pos],
                        sigma,
                        fallback_selected,
                        cache_val,
                    )
                    if not mode_selected:
                        continue
                    if abs(nu_m) < _ACOUSTIC_THRESHOLD_CM1:
                        continue

                    if mode_idx >= len(nac_tensors):
                        continue

                    if approximate_es:
                        channels_S = channels_S_base
                    else:
                        nu_S_cm1 = laser_frequency_cm1 - nu_m
                        if nu_S_cm1 <= 0.0:
                            continue
                        freq_key = round(float(nu_m), 10)
                        if freq_key not in channels_s_cache:
                            _modal_fields_S, _qs_S, _zeta_S, channels_s_cache[freq_key] = (
                                _compute_modal_fields_and_channels(
                                    nu_S_cm1,
                                    es_system,
                                    collection_angle_rad,
                                    z_s,
                                    layer_indices,
                                    gl_layer_slices,
                                    detected_pol_indices,
                                )
                            )
                        channels_S = channels_s_cache[freq_key]

                    layer_channels_S = channels_S[det_pol_idx].get(layer_index, [])
                    if j_channel >= len(layer_channels_S):
                        continue
                    channel_S = layer_channels_S[j_channel]

                    R_crystal = nac_tensors[mode_idx]
                    R_lab = rotation_matrices[layer_pos] @ R_crystal @ rotation_matrices[layer_pos].T

                    w = gl_phys_weights[sl]
                    R_E_L = R_lab @ channel_L["field"]
                    integrand = np.einsum("ij,ij->j", channel_S["field"], R_E_L)
                    amp_ij = np.dot(w, integrand)

                    if modal_pair_combination == MODAL_PAIR_GROUP_Q:
                        # Coherently sum all pairs that share the same phonon
                        # q-vector and detected-pol channel, across all layers.
                        group_key = modal_pair_q_keys.get(cache_key, ("unknown", det_pol_idx))
                    elif modal_pair_combination == MODAL_PAIR_INCOHERENT:
                        # Each (layer, i_L, j_S, det_pol) pair is squared
                        # independently — diagnostic / compatibility mode.
                        group_key = cache_key
                    else:  # MODAL_PAIR_COHERENT_ALL
                        # Single group: sum all amplitudes then square.
                        # Debug only — mixes distinct phonon-q final states.
                        group_key = "all"

                    existing = pair_groups.get(group_key)
                    if existing is None:
                        pair_groups[group_key] = {
                            "amp": amp_ij,
                            "frequency": nu_m,
                            "sigma": line_sigma,
                        }
                    elif abs(existing["frequency"] - nu_m) <= 1.0e-10:
                        existing["amp"] += amp_ij
                    else:
                        split_key = (group_key, _line_key(nu_m, line_sigma))
                        split_existing = pair_groups.get(split_key)
                        if split_existing is None:
                            pair_groups[split_key] = {
                                "amp": amp_ij,
                                "frequency": nu_m,
                                "sigma": line_sigma,
                            }
                        else:
                            split_existing["amp"] += amp_ij

    line_accumulator = {}
    for group_data in pair_groups.values():
        nu_m = group_data["frequency"]
        sigma = group_data["sigma"]
        intensity = abs(group_data["amp"]) ** 2 * bose_factor(nu_m, temperature_K)
        key = _line_key(nu_m, sigma)
        if key not in line_accumulator:
            line_accumulator[key] = [nu_m, intensity, sigma]
        else:
            line_accumulator[key][1] += intensity

    return [(mode_idx, nu_m, intensity, sigma) for nu_m, intensity, sigma in line_accumulator.values()]


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
    incident_pol : {'p', 's'}
        Incident polarisation channel.
    detected_pol : {'p', 's', 'unpolarised'}
        Detected polarisation channel.  'unpolarised' sums the p and s
        intensities incoherently (no analyser on the detector side).
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
        modes_selected=None,
    ):
        """Initialise LayeredRamanCalculator with system, layers and calculation parameters."""
        if incident_pol not in ("p", "s"):
            raise ValueError(f"incident_pol must be 'p' or 's', got '{incident_pol}'")
        if detected_pol not in ("p", "s", "unpolarised"):
            raise ValueError(f"detected_pol must be 'p', 's', or 'unpolarised', got '{detected_pol}'")
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
        self.modes_selected = None if modes_selected is None else np.asarray(modes_selected, dtype=bool)
        # NAC cache for Level 3 modal_pairs: keyed by (layer_index, i_mode, j_mode)
        self._nac_cache = {}
        self._modal_pair_q_keys = {}

        # Phase 3c: internal Jones vectors for incident and detected channels.
        # _incident_jones : (cp, cs) complex pair — determines the linear combination
        #     of p- and s-pol incident fields used in the overlap integral.
        # _detected_jones : (cp, cs) or None — None signals incoherent (no-analyser)
        #     detection: both p and s amplitudes are squared and summed.
        _pol_to_jones = {"p": np.array([1.0 + 0j, 0.0]), "s": np.array([0.0, 1.0 + 0j])}
        self._incident_jones = _pol_to_jones[incident_pol]
        self._detected_jones = None if detected_pol == "unpolarised" else _pol_to_jones[detected_pol]

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

    # ------------------------------------------------------------------
    # Level 3 (modal_pairs) helpers
    # ------------------------------------------------------------------

    def _get_active_modal_pairs(self):
        """Return the list of (i_L, j_S) Berreman mode-index pairs.

        Each Berreman eigenmode has both a forward (kz > 0) and a backward
        (kz < 0) variant per polarisation:

          Mode 0: forward p-pol  (kz > 0)
          Mode 1: forward s-pol  (kz > 0)
          Mode 2: backward p-pol (kz < 0)
          Mode 3: backward s-pol (kz < 0)

        Both propagation directions contribute to the total field, so the
        physically complete Level-3 sum includes all forward/backward pairs
        for the appropriate polarisation.

        For p-pol incidence: i_modes = {0, 2}  (both p modes, fwd and bwd).
        For s-pol incidence: i_modes = {1, 3}.
        Detected modes follow the same convention.

        In backscattering the dominant amplitude comes from pair (0, 2) because
        ``qs_L[0] + qs_S[2] ≈ 0`` gives a non-oscillating integrand.  Pairs
        with both modes forward or both backward produce rapidly oscillating
        integrands that are suppressed for thick layers, so the incoherent sum
        ``Σ|A^{ij}|²`` is still dominated by the correct backscattering terms.

        Returns
        -------
        list of (int, int)
            Cartesian product of incident × detected mode sets.
        """
        # incident_pol is restricted to 'p' or 's' by __init__ validation.
        i_modes = [0, 2] if self.incident_pol == "p" else [1, 3]
        # Detected modes from detected_pol ('p', 's', or 'unpolarised').
        if self.detected_pol == "p":
            j_modes = [0, 2]
        elif self.detected_pol == "s":
            j_modes = [1, 3]
        else:
            j_modes = [0, 1, 2, 3]
        return [(i, j) for i in i_modes for j in j_modes]

    def _get_modal_fields_at_gl_points(self, freq_cm1, system, angle_rad, z_arr):
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

        Returns
        -------
        modal_fields : ndarray, shape (4, 2, 3, N_gl_total)
            ``modal_fields[mode_n, pol_idx, xyz, z_pt]``
            mode_n: 0..3 (trans-p, trans-s, refl-p, refl-s).
            pol_idx: 0 = p-pol incidence on system, 1 = s-pol incidence.
            xyz: Ex, Ey, Ez components.
        qs_by_layer : dict of {int: ndarray(4,)}
            Maps rl.layer_index → layer.qs (4 kz eigenvalues) from ``system``.
            Note: for forward scatter (reversed system), rl.layer_index is
            used directly as an index into system.layers, which is an
            approximation; the physical-layer correspondence is exact only
            for backscattering (es_system == self.system).
        zeta : complex
            In-plane wavevector used for this system/frequency/angle.
        """
        from PDielec.Constants import speed_light_si  # already imported at module level

        freq_hz = freq_cm1 * speed_light_si * 1e2
        system.initialize_sys(freq_hz)
        zeta = np.sin(angle_rad) * np.sqrt(system.superstrate.epsilon[0, 0])

        modal_amps, zn = system.calculate_modal_amplitudes(freq_hz, zeta)

        # Layer boundaries in system's coordinate frame
        boundaries = system.get_layers_boundaries()

        N_total = len(z_arr)
        modal_fields = np.zeros((4, 2, 3, N_total), dtype=np.complex128)
        qs_by_layer = {}

        for rl_idx, (rl, sl) in enumerate(zip(self.raman_layers, self._gl_layer_slices)):
            sys_k = rl.layer_index
            if sys_k not in modal_amps or sys_k >= len(system.layers):
                continue
            layer = system.layers[sys_k]
            # calculate_modal_amplitudes returns the F_bk vector indexed as
            # k+1 in GTMcore.calculate_Efield.  That vector is propagated from
            # the back interface of finite layer k, so the partial-depth
            # exponent must use boundaries[sys_k + 2] as its reference.
            z_reference = boundaries[sys_k + 2]

            amps = modal_amps[sys_k]  # shape (8,): [0:4] p-pol, [4:8] s-pol
            qs_by_layer[rl.layer_index] = layer.qs.copy()

            z_pts = z_arr[sl]
            for pt_offset, z_j in enumerate(z_pts):
                global_idx = sl.start + pt_offset
                # Partial-depth propagation factors (same formula as calculate_Efield)
                dKiz = np.array([
                    np.exp(layer.propagation_exponents[n] * (z_j - z_reference) / layer.thick)
                    for n in range(4)
                ], dtype=np.complex128)

                for n in range(4):
                    Eprop_p = dKiz[n] * amps[n]      # p-pol incidence, mode n
                    Eprop_s = dKiz[n] * amps[n + 4]  # s-pol incidence, mode n
                    modal_fields[n, 0, :, global_idx] = Eprop_p * layer.gamma[n, :]
                    modal_fields[n, 1, :, global_idx] = Eprop_s * layer.gamma[n, :]

        return modal_fields, qs_by_layer, zeta

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
        else:
            z_s = self._gl_z
            es_system = self.system

        ref_layer = self.raman_layers[0]
        n_modes = len(ref_layer.phonon_frequencies_cm1)

        incident_pol_idx = 0 if self.incident_pol == "p" else 1
        detected_pol_indices = [0, 1] if self.detected_pol == "unpolarised" else [
            0 if self.detected_pol == "p" else 1
        ]

        # Laser modal fields at laser frequency + incident angle
        modal_fields_L, qs_L_dict, zeta_L = self._get_modal_fields_at_gl_points(
            self.laser_frequency_cm1, self.system, self.incident_angle_rad, self._gl_z
        )
        zeta_L_re = float(np.real(zeta_L))
        channels_L = self._get_modal_q_channels(modal_fields_L, qs_L_dict, incident_pol_idx)

        # ES modal fields at laser frequency (used for q_ph / NAC computation,
        # and also for the field integral when approximate_es=True)
        modal_fields_S_base, qs_S_dict, zeta_S = self._get_modal_fields_at_gl_points(
            self.laser_frequency_cm1, es_system, self.collection_angle_rad, z_s
        )
        zeta_S_re = float(np.real(zeta_S))
        channels_S_base = {
            pol_idx: self._get_modal_q_channels(modal_fields_S_base, qs_S_dict, pol_idx)
            for pol_idx in detected_pol_indices
        }

        # Build NAC cache upfront (one call per q-channel pair per layer,
        # independent of phonon mode).  Keys use channel ordinal positions rather
        # than qz values so the cache remains valid when E_S is recomputed at
        # the Stokes frequency and qz shifts slightly.
        self._nac_cache = {}
        self._modal_pair_q_keys = {}
        for rl in self.raman_layers:
            layer_channels_L = channels_L.get(rl.layer_index, [])
            for i_channel, channel_L in enumerate(layer_channels_L):
                for det_pol_idx in detected_pol_indices:
                    layer_channels_S = channels_S_base[det_pol_idx].get(rl.layer_index, [])
                    for j_channel, channel_S in enumerate(layer_channels_S):
                        cache_key = (rl.layer_index, i_channel, det_pol_idx, j_channel)

                        # q_ph = k_laser − k_scattered.  Channel qz values are
                        # coherent sums over degenerate Berreman eigenmodes.
                        q_ph = np.array([
                            zeta_L_re - zeta_S_re,
                            0.0,
                            channel_L["qz"] - channel_S["qz"],
                        ])

                        q_ph_norm = np.linalg.norm(q_ph)
                        det_key = det_pol_idx if self.detected_pol == "unpolarised" else 0
                        if q_ph_norm < 1e-8 or rl.nac_function is None:
                            self._modal_pair_q_keys[cache_key] = ("to", det_key)
                            self._nac_cache[cache_key] = None  # use TO baseline
                        else:
                            q_hat_lab = q_ph / q_ph_norm
                            q_key = tuple(float(x) for x in np.round(q_ph, decimals=10))
                            self._modal_pair_q_keys[cache_key] = (q_key, det_key)
                            self._nac_cache[cache_key] = rl.nac_function(q_hat_lab)

        layer_indices = [rl.layer_index for rl in self.raman_layers]
        rotation_matrices = [rl.rotation_matrix for rl in self.raman_layers]
        raman_tensors_by_layer = [rl.raman_tensors for rl in self.raman_layers]
        modes_selected_by_layer = [rl.modes_selected for rl in self.raman_layers]

        if pool is not None:
            shared = (
                es_system,
                self.collection_angle_rad,
                z_s,
                channels_L,
                channels_S_base,
                detected_pol_indices,
                layer_indices,
                self._gl_layer_slices,
                self._gl_phys_weights,
                rotation_matrices,
                raman_tensors_by_layer,
                modes_selected_by_layer,
                self._nac_cache,
                self._modal_pair_q_keys,
                self.modal_pair_combination,
                self.approximate_es,
                self.temperature_K,
            )
            mode_args_list = []
            for mode_idx in range(n_modes):
                fallback_nu_m = float(ref_layer.phonon_frequencies_cm1[mode_idx])
                sigma = float(self.linewidths_cm1[mode_idx]) if mode_idx < len(self.linewidths_cm1) else 5.0
                mode_args_list.append((
                    mode_idx, fallback_nu_m, sigma, self.laser_frequency_cm1,
                    self._mode_is_selected(mode_idx),
                ))

            active_freqs = []
            active_intensities = []
            active_sigmas = []
            worker_fn = partial(_compute_modal_pair_mode_worker, shared)
            for result in pool.imap(worker_fn, mode_args_list, chunksize=1):
                if progress_callback is not None:
                    progress_callback()
                if result is None:
                    continue
                for _mode_idx, nu_m, I_m, sigma in result:
                    active_freqs.append(nu_m)
                    active_intensities.append(I_m)
                    active_sigmas.append(sigma)

            return (
                np.array(active_freqs),
                np.array(active_intensities),
                np.array(active_sigmas),
            )

        active_freqs = []
        active_intensities = []
        active_sigmas = []
        contributions: list = []

        for mode_idx in range(n_modes):
            if progress_callback is not None:
                progress_callback()

            fallback_nu_m = ref_layer.phonon_frequencies_cm1[mode_idx]
            sigma = self.linewidths_cm1[mode_idx] if mode_idx < len(self.linewidths_cm1) else 5.0
            pair_groups = {}  # group_key -> {"amp", "frequency", "sigma"}
            channels_s_cache = {}

            for rl, sl in zip(self.raman_layers, self._gl_layer_slices):
                layer_channels_L = channels_L.get(rl.layer_index, [])
                for i_channel, channel_L in enumerate(layer_channels_L):
                    for det_pol_idx in detected_pol_indices:
                        layer_channels_S_base = channels_S_base[det_pol_idx].get(rl.layer_index, [])
                        for j_channel, _channel_S_base in enumerate(layer_channels_S_base):
                            cache_key = (rl.layer_index, i_channel, det_pol_idx, j_channel)
                            cache_val = self._nac_cache.get(cache_key)

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
                                    modal_fields_S, qs_S_mode, _ = self._get_modal_fields_at_gl_points(
                                        nu_S, es_system, self.collection_angle_rad, z_s
                                    )
                                    channels_s_cache[freq_key] = {
                                        pol_idx: self._get_modal_q_channels(modal_fields_S, qs_S_mode, pol_idx)
                                        for pol_idx in detected_pol_indices
                                    }
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
                                group_key = self._modal_pair_q_keys.get(cache_key, ("unknown", det_pol_idx))
                            elif self.modal_pair_combination == MODAL_PAIR_INCOHERENT:
                                group_key = cache_key
                            else:  # MODAL_PAIR_COHERENT_ALL
                                group_key = "all"

                            existing = pair_groups.get(group_key)
                            if existing is None:
                                pair_groups[group_key] = {"amp": amp_ij, "frequency": nu_m, "sigma": line_sigma}
                            elif abs(existing["frequency"] - nu_m) <= 1.0e-10:
                                existing["amp"] += amp_ij
                            else:
                                split_key = (group_key, _line_key(nu_m, line_sigma))
                                split_existing = pair_groups.get(split_key)
                                if split_existing is None:
                                    pair_groups[split_key] = {"amp": amp_ij, "frequency": nu_m, "sigma": line_sigma}
                                else:
                                    split_existing["amp"] += amp_ij

            line_accumulator = {}
            for gkey, group_data in pair_groups.items():
                nu_m = group_data["frequency"]
                sigma = group_data["sigma"]
                amp = group_data["amp"]
                bose = bose_factor(nu_m, self.temperature_K)
                intensity = abs(amp) ** 2 * bose
                lkey = _line_key(nu_m, sigma)
                if lkey not in line_accumulator:
                    line_accumulator[lkey] = [nu_m, intensity, sigma]
                else:
                    line_accumulator[lkey][1] += intensity
                if return_contributions:
                    contributions.append(RamanContribution(
                        mode_idx=mode_idx, frequency=nu_m, sigma=sigma,
                        group_key=gkey, amplitude=amp, intensity=intensity,
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
