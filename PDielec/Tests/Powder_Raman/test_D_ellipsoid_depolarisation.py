"""Tests D1-D4: Ellipsoidal particle depolarisation tensors.

Plan sections: D1, D2, D3, D4, J7.

Tests verify:
- the spherical limit (a = b = c) gives L = I/3 (D1);
- the sum rule L_x + L_y + L_z = 1 (D2 / J7);
- limiting behaviour for needles and plates (D2);
- isotropic particle in ellipsoid gives axis-by-axis local-field factors (D3);
- rotation of particle shape relative to crystal axes (D4).

Uses PDielec.Calculator depolarisation functions directly.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import (
    compute_internal_field_tensor,
    compute_particle_raman_tensor,
    compute_powder_raman_intensities,
    initialise_sphere_depolarisation_matrix,
    initialise_ellipsoid_depolarisation_matrix,
    initialise_plate_depolarisation_matrix,
    initialise_needle_depolarisation_matrix,
)


def _random_rotation(rng=None):
    if rng is None:
        rng = np.random.default_rng(12345)
    A = rng.standard_normal((3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


# ---------------------------------------------------------------------------
# D1 — Spherical limit of ellipsoid (a = b = c)
# ---------------------------------------------------------------------------

def test_d1_ellipsoid_sphere_limit_L_is_identity_over_3():
    """D1: Ellipsoid with a/b = 1 gives L = I/3 (sphere depolarisation)."""
    unique = np.array([0.0, 0.0, 1.0])
    L_ellipsoid = initialise_ellipsoid_depolarisation_matrix(unique, 1.0)
    L_sphere = initialise_sphere_depolarisation_matrix()
    assert np.allclose(L_ellipsoid, L_sphere, atol=1e-6), \
        f"Ellipsoid a/b=1 != sphere L:\n{L_ellipsoid}\nvs\n{L_sphere}"


def test_d1_ellipsoid_sphere_limit_all_eigenvalues_equal():
    """D1: All depolarisation factors should equal 1/3 in the sphere limit."""
    unique = np.array([1.0, 0.0, 0.0])
    L = initialise_ellipsoid_depolarisation_matrix(unique, 1.0)
    eigenvalues = np.linalg.eigvalsh(np.real(L))
    assert np.allclose(eigenvalues, 1.0 / 3.0, atol=1e-6)


def test_d1_sphere_limit_recovers_sphere_Raman():
    """D1: Ellipsoid a/b=1 Raman result equals sphere Raman result."""
    unique = np.array([0.0, 0.0, 1.0])
    L_ell = initialise_ellipsoid_depolarisation_matrix(unique, 1.0)
    L_sph = initialise_sphere_depolarisation_matrix()
    eps_i = np.diag([3.0, 3.0, 3.0]).astype(complex)
    eps_e = 1.0
    R = np.diag([2.0, 1.0, 0.5])

    N_ell = compute_internal_field_tensor(L_ell, eps_i, eps_e)
    N_sph = compute_internal_field_tensor(L_sph, eps_i, eps_e)
    R_ell = compute_particle_raman_tensor(R, N_ell, L_ell, eps_i, eps_e)
    R_sph = compute_particle_raman_tensor(R, N_sph, L_sph, eps_i, eps_e)

    vv_ell, vh_ell = compute_powder_raman_intensities(R_ell)
    vv_sph, vh_sph = compute_powder_raman_intensities(R_sph)
    assert abs(vv_ell - vv_sph) < 1e-6
    assert abs(vh_ell - vh_sph) < 1e-6


# ---------------------------------------------------------------------------
# D2 / J7 — Sum rule L_x + L_y + L_z = 1
# ---------------------------------------------------------------------------

def _depolarisation_principal_values(L):
    """Extract principal depolarisation factors as eigenvalues of the real part of L."""
    return np.sort(np.linalg.eigvalsh(np.real(L)))


def test_d2_sphere_sum_rule():
    """D2/J7: Sphere depolarisation factors sum to 1."""
    L = initialise_sphere_depolarisation_matrix()
    total = np.trace(np.real(L))
    assert abs(total - 1.0) < 1e-12


def test_d2_plate_sum_rule():
    """D2/J7: Plate depolarisation factors sum to 1."""
    L = initialise_plate_depolarisation_matrix(np.array([0.0, 0.0, 1.0]))
    total = np.trace(np.real(L))
    assert abs(total - 1.0) < 1e-12


def test_d2_needle_sum_rule():
    """D2/J7: Needle depolarisation factors sum to 1."""
    L = initialise_needle_depolarisation_matrix(np.array([0.0, 0.0, 1.0]))
    total = np.trace(np.real(L))
    assert abs(total - 1.0) < 1e-12


def test_d2_ellipsoid_sum_rule_various_aspect_ratios():
    """D2/J7: Ellipsoid depolarisation factors sum to 1 for various aspect ratios."""
    unique = np.array([0.0, 0.0, 1.0])
    for aoverb in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        L = initialise_ellipsoid_depolarisation_matrix(unique, aoverb)
        total = np.trace(np.real(L))
        assert abs(total - 1.0) < 1e-10, f"a/b={aoverb}: trace={total}"


def test_d2_ellipsoid_principal_values_between_0_and_1():
    """D2: All depolarisation factors must lie in [0, 1]."""
    unique = np.array([0.0, 0.0, 1.0])
    for aoverb in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        L = initialise_ellipsoid_depolarisation_matrix(unique, aoverb)
        vals = _depolarisation_principal_values(L)
        assert np.all(vals >= -1e-10), f"a/b={aoverb}: negative eigenvalue {vals}"
        assert np.all(vals <= 1.0 + 1e-10), f"a/b={aoverb}: eigenvalue > 1: {vals}"


def test_d2_prolate_needle_limit():
    """D2: For very elongated prolate ellipsoid (a >> b), L along long axis → 0."""
    unique = np.array([0.0, 0.0, 1.0])
    L_needle = initialise_needle_depolarisation_matrix(unique)
    # Needle: along unique axis L ≈ 0; perpendicular L ≈ 0.5
    L_along = float(np.real((unique @ L_needle @ unique)))
    assert L_along < 0.1, f"Needle L_along_axis={L_along}, expected ~0"


def test_d2_oblate_plate_limit():
    """D2: For very flat oblate ellipsoid (a << b), L along flat axis → 1."""
    unique = np.array([0.0, 0.0, 1.0])
    L_plate = initialise_plate_depolarisation_matrix(unique)
    L_normal = float(np.real((unique @ L_plate @ unique)))
    assert L_normal > 0.9, f"Plate L_normal={L_normal}, expected ~1"


def test_d2_ellipsoid_prolate_long_axis_small():
    """D2: Prolate ellipsoid (a/b > 1) has smaller L along the unique axis."""
    unique = np.array([0.0, 0.0, 1.0])
    L_prolate = initialise_ellipsoid_depolarisation_matrix(unique, 3.0)
    L_along = float(np.real((unique @ L_prolate @ unique)))
    # Should be less than 1/3 for prolate
    assert L_along < 1.0 / 3.0, \
        f"Prolate a/b=3: L_along={L_along}, expected < 1/3"


def test_d2_ellipsoid_oblate_normal_axis_large():
    """D2: Oblate ellipsoid (a/b < 1) has larger L along the unique axis."""
    unique = np.array([0.0, 0.0, 1.0])
    L_oblate = initialise_ellipsoid_depolarisation_matrix(unique, 0.3)
    L_along = float(np.real((unique @ L_oblate @ unique)))
    # Should be greater than 1/3 for oblate
    assert L_along > 1.0 / 3.0, \
        f"Oblate a/b=0.3: L_along={L_along}, expected > 1/3"


# ---------------------------------------------------------------------------
# D3 — Isotropic particle in ellipsoid: axis-by-axis local-field factors
# ---------------------------------------------------------------------------

def test_d3_isotropic_particle_ellipsoid_N_diagonal_in_principal_frame():
    """D3: Isotropic particle in ellipsoid gives diagonal N in the ellipsoid principal frame."""
    unique = np.array([0.0, 0.0, 1.0])
    aoverb = 2.0
    L = initialise_ellipsoid_depolarisation_matrix(unique, aoverb)
    eps_i_scalar = 4.0
    eps_e = 1.0
    eps_i = eps_i_scalar * np.eye(3, dtype=complex)
    N = compute_internal_field_tensor(L, eps_i, eps_e)
    # L is already diagonal when unique = z-axis; check N is also diagonal
    off_diag = N - np.diag(np.diag(N))
    assert np.allclose(off_diag, 0, atol=1e-12)


def test_d3_N_axis_formula_for_scalar_eps_and_ellipsoid():
    """D3: Each N_j satisfies N_j = ε_e / (ε_e + L_j*(ε_i - ε_e)) for scalar ε_i."""
    unique = np.array([0.0, 0.0, 1.0])
    aoverb = 2.0
    L = initialise_ellipsoid_depolarisation_matrix(unique, aoverb)
    eps_i_s = 4.0
    eps_e = 1.0
    eps_i = eps_i_s * np.eye(3, dtype=complex)
    N = compute_internal_field_tensor(L, eps_i, eps_e)
    L_diag = np.diag(np.real(L))
    for j in range(3):
        N_j_expected = eps_e / (eps_e + L_diag[j] * (eps_i_s - eps_e))
        assert abs(N[j, j] - N_j_expected) < 1e-12, \
            f"j={j}: N[{j},{j}]={N[j,j]}, expected {N_j_expected}"


def test_d3_sphere_limit_matches_sphere_N():
    """D3: Ellipsoid a/b=1 local field matches sphere N = (3ε_e/(ε_i+2ε_e))*I."""
    unique = np.array([0.0, 0.0, 1.0])
    L_ell = initialise_ellipsoid_depolarisation_matrix(unique, 1.0)
    L_sph = initialise_sphere_depolarisation_matrix()
    eps_i_s = 4.0
    eps_e = 1.0
    eps_i = eps_i_s * np.eye(3, dtype=complex)
    N_ell = compute_internal_field_tensor(L_ell, eps_i, eps_e)
    N_sph = compute_internal_field_tensor(L_sph, eps_i, eps_e)
    assert np.allclose(N_ell, N_sph, atol=1e-6)


# ---------------------------------------------------------------------------
# D4 — Rotation of particle shape relative to crystal axes
# ---------------------------------------------------------------------------

def test_d4_rotating_all_tensors_together_leaves_powder_average_unchanged():
    """D4: Co-rotating R, ε_i, and L leaves the powder-averaged intensities unchanged."""
    unique = np.array([0.0, 0.0, 1.0])
    L_crystal = initialise_ellipsoid_depolarisation_matrix(unique, 2.0)
    eps_i_crystal = np.diag([2.0, 2.0, 6.0]).astype(complex)
    R_crystal = np.diag([1.0, 1.5, 0.5])

    N_crystal = compute_internal_field_tensor(L_crystal, eps_i_crystal, 1.0)
    R_p_crystal = compute_particle_raman_tensor(R_crystal, N_crystal, L_crystal, eps_i_crystal, 1.0)
    vv_ref, vh_ref = compute_powder_raman_intensities(R_p_crystal)

    rng = np.random.default_rng(55)
    for _ in range(6):
        G0 = _random_rotation(rng)
        L_rot = G0 @ np.real(L_crystal) @ G0.T
        eps_i_rot = G0 @ eps_i_crystal @ G0.T
        R_rot = G0 @ R_crystal @ G0.T
        N_rot = compute_internal_field_tensor(L_rot, eps_i_rot, 1.0)
        R_p_rot = compute_particle_raman_tensor(R_rot, N_rot, L_rot, eps_i_rot, 1.0)
        vv, vh = compute_powder_raman_intensities(R_p_rot)
        assert abs(vv - vv_ref) < 1e-8, f"VV changed after co-rotation: {vv} vs {vv_ref}"
        assert abs(vh - vh_ref) < 1e-8, f"VH changed after co-rotation: {vh} vs {vh_ref}"


def test_d4_rotating_only_R_changes_result_for_anisotropic_case():
    """D4: Rotating only R (not L or ε_i) changes the single-orientation result."""
    unique = np.array([0.0, 0.0, 1.0])
    L = initialise_ellipsoid_depolarisation_matrix(unique, 2.0)
    eps_i = np.diag([2.0, 2.0, 6.0]).astype(complex)
    R = np.diag([1.0, 1.5, 0.5])
    eps_e = 1.0

    N = compute_internal_field_tensor(L, eps_i, eps_e)
    R_p_orig = compute_particle_raman_tensor(R, N, L, eps_i, eps_e)
    vv_ref, vh_ref = compute_powder_raman_intensities(R_p_orig)

    rng = np.random.default_rng(321)
    G0 = _random_rotation(rng)
    R_rotated = G0 @ R @ G0.T
    R_p_rotated = compute_particle_raman_tensor(R_rotated, N, L, eps_i, eps_e)
    vv, vh = compute_powder_raman_intensities(R_p_rotated)
    # For anisotropic L and ε_i, rotating only R should change the single-orientation result
    assert not (abs(vv - vv_ref) < 1e-10 and abs(vh - vh_ref) < 1e-10), \
        "Rotating only R did not change the result (expected a change for anisotropic case)"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
