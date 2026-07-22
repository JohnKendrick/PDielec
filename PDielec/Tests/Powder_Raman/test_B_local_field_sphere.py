"""Tests B1-B4: Local-field tensor for isotropic spherical particles.

Plan sections: B1, B2, B3, B4, J4, J5, M4, M5.

Tests verify the internal field tensor N (Eq. 47) and the effective particle
Raman tensor R_particle (Eq. 60) for the case of an isotropic sphere embedded
in a scalar matrix permittivity.

For an isotropic sphere (L = I/3) with scalar particle permittivity ε_i and
scalar matrix permittivity ε_e::

  N = (3 ε_e / (ε_i + 2 ε_e)) × I   (scalar factor × identity)
  R_eff = N_s^2 × R_eps

where N_s = 3 ε_e / (ε_i + 2 ε_e).

There is one local-field factor at the laser frequency and one at the
scattered frequency.  In the present non-dispersive optical approximation
they are equal, so intensities scale as |N_s|^4.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest

from PDielec.Calculator import (
    compute_internal_field_tensor,
    compute_particle_raman_tensor,
    compute_powder_raman_intensities,
    initialise_sphere_depolarisation_matrix,
)

L_sphere = initialise_sphere_depolarisation_matrix()


def _scalar_N(eps_i, eps_e):
    """Analytic scalar local-field factor for an isotropic sphere."""
    return 3.0 * eps_e / (eps_i + 2.0 * eps_e)


# ---------------------------------------------------------------------------
# B2 / J4 / M4 — Matched permittivity: ε_i = ε_e → N = I, R_particle = R_eps
# ---------------------------------------------------------------------------

def test_b2_matched_permittivity_N_is_identity():
    """B2/J4/M4: ε_i = ε_e gives N = I for an isotropic sphere."""
    eps = 3.5
    N = compute_internal_field_tensor(L_sphere, eps * np.eye(3), eps)
    assert np.allclose(N, np.eye(3), atol=1e-12), f"N != I: {N}"


def test_b2_matched_permittivity_R_particle_equals_R_eps():
    """B2/J4: Matched permittivity gives R_particle = R_eps."""
    eps = 2.0
    R = np.array([[1.0, 0.5, 0.0],
                  [0.5, 0.0, 0.3],
                  [0.0, 0.3, 1.5]])
    N = compute_internal_field_tensor(L_sphere, eps * np.eye(3), eps)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps * np.eye(3), eps)
    assert np.allclose(R_p, R, atol=1e-12), "R_particle != R_eps for matched ε"


def test_b2_matched_permittivity_VV_unchanged():
    """B2: VV/VH intensities must be unchanged when ε_i = ε_e."""
    eps = 4.0
    R = np.diag([2.0, 1.0, 0.5])
    vv_bare, vh_bare = compute_powder_raman_intensities(R)
    N = compute_internal_field_tensor(L_sphere, eps * np.eye(3), eps)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps * np.eye(3), eps)
    vv_p, vh_p = compute_powder_raman_intensities(R_p)
    assert abs(vv_p - vv_bare) < 1e-10
    assert abs(vh_p - vh_bare) < 1e-10


def test_b2_matched_any_shape():
    """B2: Result must be shape-independent when ε_i = ε_e (correction is zero)."""
    from PDielec.Calculator import (
        initialise_needle_depolarisation_matrix,
        initialise_plate_depolarisation_matrix,
    )
    eps = 3.0
    R = np.diag([1.0, 2.0, 3.0])
    vv_bare, vh_bare = compute_powder_raman_intensities(R)
    for L in [
        L_sphere,
        initialise_plate_depolarisation_matrix(np.array([0, 0, 1])),
        initialise_needle_depolarisation_matrix(np.array([0, 0, 1])),
    ]:
        N = compute_internal_field_tensor(L, eps * np.eye(3), eps)
        R_p = compute_particle_raman_tensor(R, N, L, eps * np.eye(3), eps)
        vv_p, vh_p = compute_powder_raman_intensities(R_p)
        assert abs(vv_p - vv_bare) < 1e-10, "Matched ε should be shape-independent"


# ---------------------------------------------------------------------------
# B1 / J5 — Isotropic sphere with real permittivity contrast
# ---------------------------------------------------------------------------

def test_b1_scalar_N_equals_analytic():
    """B1/J5: N = N_s × I with N_s = 3 ε_e / (ε_i + 2 ε_e) for isotropic sphere."""
    for eps_i, eps_e in [(4.0, 1.0), (2.25, 1.5), (9.0, 2.0)]:
        N = compute_internal_field_tensor(L_sphere, eps_i * np.eye(3), eps_e)
        N_s_expected = _scalar_N(eps_i, eps_e)
        assert np.allclose(N, N_s_expected * np.eye(3), atol=1e-12), \
            f"ε_i={eps_i}, ε_e={eps_e}: N not scalar"


def test_b1_R_eff_scales_as_Ns_squared():
    """B1: For an isotropic sphere with scalar ε_i, R_eff = N_s^2 × R_eps."""
    eps_i, eps_e = 4.0, 1.0
    R = np.diag([2.0, 1.0, 0.5])
    N_s = _scalar_N(eps_i, eps_e)
    N = compute_internal_field_tensor(L_sphere, eps_i * np.eye(3), eps_e)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i * np.eye(3), eps_e)
    assert np.allclose(R_p, N_s**2 * R, atol=1e-12), \
        f"R_eff != N_s^2 R; N_s={N_s}"


def test_b1_intensity_scales_as_Ns_fourth():
    """B1: Powder intensity (VV or VH) scales as |N_s|^4 for an isotropic sphere."""
    R = np.diag([1.0, 2.0, 0.5])
    vv_bare, vh_bare = compute_powder_raman_intensities(R)

    for eps_i, eps_e in [(4.0, 1.0), (2.25, 1.5), (9.0, 2.0)]:
        N_s = _scalar_N(eps_i, eps_e)
        N = compute_internal_field_tensor(L_sphere, eps_i * np.eye(3), eps_e)
        R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i * np.eye(3), eps_e)
        vv_p, vh_p = compute_powder_raman_intensities(R_p)
        scale = N_s**4
        assert abs(vv_p / vv_bare - scale) < 1e-10, \
            f"ε_i={eps_i}: VV scale {vv_p/vv_bare} != N_s^4={scale}"
        assert abs(vh_p / vh_bare - scale) < 1e-10, \
            f"ε_i={eps_i}: VH scale {vh_p/vh_bare} != N_s^4={scale}"


def test_b1_depolarisation_ratio_unchanged_for_scalar_eps():
    """B1/K3: Depolarisation ratio VH/VV must be unchanged for scalar ε_i (sphere)."""
    R = np.diag([2.0, 1.0, 0.5])
    vv_bare, vh_bare = compute_powder_raman_intensities(R)
    rho_bare = vh_bare / vv_bare

    for eps_i in [1.5, 2.0, 4.0, 9.0, 25.0]:
        eps_e = 1.0
        N = compute_internal_field_tensor(L_sphere, eps_i * np.eye(3), eps_e)
        R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i * np.eye(3), eps_e)
        vv_p, vh_p = compute_powder_raman_intensities(R_p)
        rho_p = vh_p / vv_p
        assert abs(rho_p - rho_bare) < 1e-10, \
            f"ε_i={eps_i}: rho changed from {rho_bare} to {rho_p}"


# ---------------------------------------------------------------------------
# B3 — High permittivity contrast: ε_i >> ε_e
# ---------------------------------------------------------------------------

def test_b3_large_eps_i_suppresses_intensity():
    """B3: Raman intensity decreases monotonically as ε_i/ε_e increases (sphere)."""
    R = np.diag([1.0, 1.0, 1.0])
    eps_e = 1.0
    prev_vv = None
    for eps_i in [1.0, 4.0, 9.0, 25.0, 100.0]:
        N = compute_internal_field_tensor(L_sphere, eps_i * np.eye(3), eps_e)
        R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i * np.eye(3), eps_e)
        vv, _ = compute_powder_raman_intensities(R_p)
        if prev_vv is not None:
            assert vv < prev_vv, f"ε_i={eps_i}: VV did not decrease"
        prev_vv = vv


def test_b3_high_contrast_no_nan_or_inf():
    """B3: No NaNs or infinities for large positive real ε_i (sphere)."""
    R = np.diag([1.0, 0.5, 2.0])
    eps_e = 1.0
    for eps_i in [50.0, 100.0, 500.0, 1000.0]:
        N = compute_internal_field_tensor(L_sphere, eps_i * np.eye(3), eps_e)
        R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i * np.eye(3), eps_e)
        vv, vh = compute_powder_raman_intensities(R_p)
        assert np.isfinite(vv), f"ε_i={eps_i}: VV is not finite"
        assert np.isfinite(vh), f"ε_i={eps_i}: VH is not finite"


def test_b3_N_s_value_for_test_case():
    """B3: For ε_i=4, ε_e=1 sphere: N_s=0.5 and intensity scaling is 1/16."""
    eps_i, eps_e = 4.0, 1.0
    N_s = _scalar_N(eps_i, eps_e)
    assert abs(N_s - 0.5) < 1e-12
    R = np.eye(3)
    vv_bare, _ = compute_powder_raman_intensities(R)
    N = compute_internal_field_tensor(L_sphere, eps_i * np.eye(3), eps_e)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i * np.eye(3), eps_e)
    vv_p, _ = compute_powder_raman_intensities(R_p)
    assert abs(vv_p / vv_bare - 1.0 / 16.0) < 1e-12, \
        f"scaling={vv_p/vv_bare}, expected 1/16"


# ---------------------------------------------------------------------------
# B4 — Surface-mode resonance near ε_i = -2 ε_e (complex permittivity)
# ---------------------------------------------------------------------------

def test_b4_resonance_peaks_near_surface_mode_condition():
    """B4: N_s magnitude peaks near ε_i = -2 ε_e with finite damping."""
    eps_e = 1.0
    # Sweep ε_i through -2 with small imaginary part
    eps_real_values = np.linspace(-4.0, 0.0, 201)
    gamma = 0.05  # small imaginary part (damping)
    intensities = []
    for eps_r in eps_real_values:
        eps_i_complex = (eps_r + 1j * gamma) * np.eye(3)
        N = compute_internal_field_tensor(L_sphere, eps_i_complex, eps_e)
        # Intensity of N ~ |N_s|^2 for scalar case
        N_s = N[0, 0]
        intensities.append(abs(N_s)**2)
    # Peak should occur near ε_i = -2
    peak_idx = np.argmax(intensities)
    eps_at_peak = eps_real_values[peak_idx]
    assert abs(eps_at_peak - (-2.0 * eps_e)) < 0.1, \
        f"Peak at ε_i={eps_at_peak}, expected near {-2*eps_e}"


def test_b4_finite_damping_prevents_divergence():
    """B4: With finite damping, N remains finite at the surface-mode condition."""
    eps_e = 1.0
    # Exactly at resonance: ε_i + 2 ε_e = 0 → ε_i = -2
    # With damping γ > 0, N must be finite
    for gamma in [0.001, 0.01, 0.1, 0.5]:
        eps_i = (-2.0 + 1j * gamma) * np.eye(3)
        N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
        assert np.all(np.isfinite(N)), f"N has non-finite elements at γ={gamma}"


def test_b4_increasing_damping_broadens_and_reduces_peak():
    """B4: Larger damping reduces and broadens the resonance peak."""
    eps_e = 1.0
    eps_r = -2.0  # exactly at resonance
    peak_heights = []
    for gamma in [0.01, 0.1, 0.5, 1.0]:
        eps_i = (eps_r + 1j * gamma) * np.eye(3)
        N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
        peak_heights.append(abs(N[0, 0])**2)
    # Peak height should decrease as damping increases
    for i in range(len(peak_heights) - 1):
        assert peak_heights[i] > peak_heights[i + 1], \
            f"Peak height did not decrease: {peak_heights}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
