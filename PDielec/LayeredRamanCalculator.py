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
"""Layered Raman scattering calculator using GTMcore for optical field propagation.

The Raman amplitude for phonon mode m is computed as a source-overlap integral
over the Raman-active layers in a multilayer stack:

    A_m = Σ_ℓ ∫ E_S(z)^T R_{ℓ,m}^lab E_L(z) dz

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
- Coherent layer summation: set ``coherent_layers=True`` to sum amplitudes
  across Raman-active layers before squaring.

Phase 3c Jones-vector polarisation:
- Incident and detected channels are internally represented as Jones vectors
  ``[cp, cs]`` (complex coefficients for p- and s-pol field components).  The
  string shortcuts ``incident_pol='p'`` → ``[1, 0]``, ``'s'`` → ``[0, 1]``,
  and ``detected_pol='unpolarised'`` → incoherent sum of |A_p|² + |A_s|² are
  converted automatically.  Direct Jones vector control can be added later via
  a subclass or optional parameter.

The integral is evaluated numerically using Gauss-Legendre quadrature within
each Raman-active layer.
"""

import logging

import numpy as np
from numpy.polynomial.legendre import leggauss

from PDielec.Constants import boltzmann_si, planck_si, speed_light_si

logger = logging.getLogger(__name__)

# Modes below this frequency (cm⁻¹) are treated as acoustic and skipped.
_ACOUSTIC_THRESHOLD_CM1 = 10.0


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
        mode.  The stored reader tensors have units (Å/amu)^0.5 (i.e. they
        are divided by √V_cell); multiply by ``np.sqrt(volume_ang3)`` before
        passing here so that intensities are proportional to Å² × V / amu.
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

    def __init__(self, layer_index, phonon_frequencies_cm1, raman_tensors, rotation_matrix):
        """Initialise a RamanLayer descriptor."""
        self.layer_index = layer_index
        self.phonon_frequencies_cm1 = np.asarray(phonon_frequencies_cm1, dtype=float)
        self.raman_tensors = [np.asarray(R, dtype=float) for R in raman_tensors]
        self.rotation_matrix = np.asarray(rotation_matrix, dtype=float)


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
    """Calculate the Raman spectrum of a multilayer system via depth integration.

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
          - ``coherent_layers=False`` (default): I_m ∝ Σ_ℓ |A_{ℓ,m}|²
          - ``coherent_layers=True``:            I_m ∝ |Σ_ℓ A_{ℓ,m}|²
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
        Default is 20, which is accurate for smooth fields in thin layers.
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
    passed to ``calculate_Efield`` on the reversed system; the GL quadrature
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
        n_gauss=20,
        collection_side="superstrate",
        collection_angle_rad=None,
        coherent_layers=False,
        approximate_es=False,
    ):
        """Initialise LayeredRamanCalculator with system, layers and calculation parameters."""
        if incident_pol not in ("p", "s"):
            raise ValueError(f"incident_pol must be 'p' or 's', got '{incident_pol}'")
        if detected_pol not in ("p", "s", "unpolarised"):
            raise ValueError(f"detected_pol must be 'p', 's', or 'unpolarised', got '{detected_pol}'")
        if collection_side not in ("superstrate", "substrate"):
            raise ValueError(f"collection_side must be 'superstrate' or 'substrate', got '{collection_side}'")

        self.system = system
        self.raman_layers = list(raman_layers)
        self.laser_frequency_cm1 = float(laser_frequency_cm1)
        self.incident_angle_rad = float(incident_angle_rad)
        self.incident_pol = incident_pol
        self.detected_pol = detected_pol
        self.temperature_K = float(temperature_K)
        self.linewidths_cm1 = np.asarray(linewidths_cm1, dtype=float)
        self.n_gauss = int(n_gauss)
        self.collection_side = collection_side
        self.collection_angle_rad = float(collection_angle_rad) if collection_angle_rad is not None else float(incident_angle_rad)
        self.coherent_layers = bool(coherent_layers)
        self.approximate_es = bool(approximate_es)

        # Phase 3c: internal Jones vectors for incident and detected channels.
        # _incident_jones : (cp, cs) complex pair — determines the linear combination
        #     of p- and s-pol incident fields used in the overlap integral.
        # _detected_jones : (cp, cs) or None — None signals incoherent (no-analyser)
        #     detection: both p and s amplitudes are squared and summed.
        _pol_to_jones = {"p": np.array([1.0 + 0j, 0.0]), "s": np.array([0.0, 1.0 + 0j])}
        self._incident_jones = _pol_to_jones[incident_pol]
        self._detected_jones = None if detected_pol == "unpolarised" else _pol_to_jones[detected_pol]

        # Gauss-Legendre nodes and weights on [-1, 1]
        self._gl_nodes, self._gl_weights = leggauss(self.n_gauss)

        # Built by _build_gl_grid()
        self._gl_z = None               # ndarray: concatenated z points (m) for E_L
        self._gl_phys_weights = None    # ndarray: corresponding physical integration weights (m)
        self._gl_layer_slices = None    # list of slice: one per RamanLayer

        self._build_gl_grid()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_gl_grid(self):
        """Map GL quadrature points onto each Raman-active layer interval.

        Uses ``system.get_layers_boundaries()`` to find each layer's z extent.
        Boundaries are in metres with z = 0 at the superstrate/first-layer
        interface (the convention used by ``calculate_Efield``).

        Populates
        ---------
        _gl_z : ndarray
            Concatenated GL z coordinates for all Raman-active layers, in metres
            (the unit ``calculate_Efield`` expects).
        _gl_phys_weights : ndarray
            Corresponding integration weights in Ångström (GL weight × half-thickness × 1e10).
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
            z_j = mid + half * self._gl_nodes              # metres (for calculate_Efield)
            w_j = half * self._gl_weights * 1.0e10         # Å (for the intensity integral)
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

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def calculate_mode_intensities(self):
        """Compute the per-mode Raman intensities (before broadening).

        Modes below ``_ACOUSTIC_THRESHOLD_CM1`` are excluded.  Layer
        contributions are combined incoherently or coherently depending on
        ``self.coherent_layers``.

        Returns
        -------
        active_frequencies : ndarray, shape (M,)
            Frequencies of the active (non-acoustic) modes in cm⁻¹.
        active_intensities : ndarray, shape (M,)
            Corresponding Raman intensities I_m (real, non-negative).
        active_linewidths : ndarray, shape (M,)
            Linewidths of the active modes in cm⁻¹.

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
        if not self.raman_layers:
            logger.warning("calculate_mode_intensities: no Raman-active layers defined")
            return np.array([]), np.array([]), np.array([])

        # --- E_L: incident laser field on the original system ---
        E_L_out = self._get_field_at_gl_points(
            self.laser_frequency_cm1, self.system, self.incident_angle_rad, self._gl_z
        )

        # --- Set up the system and z-array for E_S ---
        if self.collection_side == "substrate":
            # Forward scattering: launch E_S from the substrate side
            total_thick = sum(layer.thick for layer in self.system.layers)
            z_s = total_thick - self._gl_z   # z in the reversed system
            es_system = self.system.reversed_system()
        else:
            # Backscattering: same system, same z-array
            z_s = self._gl_z
            es_system = self.system

        # If approximate_es, compute E_S once at the laser frequency
        if self.approximate_es:
            if self.collection_side == "substrate":
                E_S_out_fixed = self._get_field_at_gl_points(
                    self.laser_frequency_cm1, es_system, self.collection_angle_rad, z_s
                )
            else:
                # Pure Phase-1 approximation: E_S = E_L (no extra call)
                E_S_out_fixed = E_L_out

        # Use frequencies from the first RamanLayer (all layers share the same DFT modes)
        ref_layer = self.raman_layers[0]
        n_modes = len(ref_layer.phonon_frequencies_cm1)

        active_freqs = []
        active_intensities = []
        active_sigmas = []

        for mode_idx in range(n_modes):
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

            # Accumulate amplitude contributions from all Raman-active layers
            if self.coherent_layers:
                # Coherent: sum amplitudes first
                total_amp_p = 0.0 + 0.0j
                total_amp_s = 0.0 + 0.0j
                for rl, sl in zip(self.raman_layers, self._gl_layer_slices, strict=True):
                    R_crystal = rl.raman_tensors[mode_idx]
                    R_lab = self._rotate_raman_tensor(R_crystal, rl.rotation_matrix)
                    amp_p, amp_s = self._layer_amplitude(E_L_out, E_S_out, R_lab, sl)
                    total_amp_p += amp_p
                    if self._detected_jones is None:
                        total_amp_s += amp_s
                if self._detected_jones is None:
                    I_m = abs(total_amp_p) ** 2 + abs(total_amp_s) ** 2
                else:
                    I_m = abs(total_amp_p) ** 2
            else:
                # Incoherent: sum intensities
                I_m = 0.0
                for rl, sl in zip(self.raman_layers, self._gl_layer_slices, strict=True):
                    R_crystal = rl.raman_tensors[mode_idx]
                    R_lab = self._rotate_raman_tensor(R_crystal, rl.rotation_matrix)
                    amp_p, amp_s = self._layer_amplitude(E_L_out, E_S_out, R_lab, sl)
                    if self._detected_jones is None:
                        I_m += abs(amp_p) ** 2 + abs(amp_s) ** 2
                    else:
                        I_m += abs(amp_p) ** 2

            # Apply Bose-Einstein thermal prefactor
            I_m *= bose_factor(nu_m, self.temperature_K)

            active_freqs.append(nu_m)
            active_intensities.append(I_m)
            active_sigmas.append(sigma)

        return (
            np.array(active_freqs),
            np.array(active_intensities),
            np.array(active_sigmas),
        )

    def calculate_spectrum(self, freq_axis_cm1):
        """Compute the broadened Raman spectrum on the given frequency axis.

        Parameters
        ----------
        freq_axis_cm1 : array_like, shape (N,)
            Raman-shift axis in cm⁻¹ on which to evaluate the spectrum.

        Returns
        -------
        spectrum : ndarray, shape (N,)
            Raman intensity I(Δν) on the requested axis.

        """
        freq_axis_cm1 = np.asarray(freq_axis_cm1)
        active_freqs, active_intensities, active_sigmas = self.calculate_mode_intensities()
        if len(active_freqs) == 0:
            return np.zeros(len(freq_axis_cm1))
        return lorentzian_broaden(active_freqs, active_intensities, active_sigmas, freq_axis_cm1)
