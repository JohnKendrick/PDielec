"""Tests C1-C3: Anisotropic permittivity in a spherical particle.

Plan sections: C1, C2, C3, J6, J8, K3, K4, M6.

For a sphere (L = I/3) with diagonal anisotropic permittivity
    ε_i = diag(ε_x, ε_y, ε_z)
the local-field tensor N is diagonal in the same frame with components
    N_j = 3 ε_e / (ε_j + 2 ε_e).

Tests verify:
- diagonal structure of N for diagonal ε_i (C1, J6);
- rotation invariance of powder averages when all tensors are co-rotated (C2, J8, M6);
- complex symmetric off-diagonal ε_i handled correctly (C3);
- isotropic R with anisotropic ε_i can make the effective tensor anisotropic (K4).
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
    sobol_rotations,
)


def _random_rotation(rng=None):
    if rng is None:
        rng = np.random.default_rng(12345)
    A = rng.standard_normal((3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q

L_sphere = initialise_sphere_depolarisation_matrix()


# ---------------------------------------------------------------------------
# C1 / J6 — Diagonal anisotropic ε_i aligned with Raman tensor
# ---------------------------------------------------------------------------

def _analytic_N_diagonal_sphere(eps_x, eps_y, eps_z, eps_e):
    """Analytic N tensor for sphere with diagonal ε_i = diag(ε_x, ε_y, ε_z)."""
    Nx = 3.0 * eps_e / (eps_x + 2.0 * eps_e)
    Ny = 3.0 * eps_e / (eps_y + 2.0 * eps_e)
    Nz = 3.0 * eps_e / (eps_z + 2.0 * eps_e)
    return np.diag([Nx, Ny, Nz])


def test_c1_N_is_diagonal_for_diagonal_eps_i():
    """C1/J6: N is diagonal for sphere with diagonal ε_i aligned with Raman tensor."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    eps_i = np.diag([eps_x, eps_y, eps_z]).astype(complex)
    N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
    off_diag = N - np.diag(np.diag(N))
    assert np.allclose(off_diag, 0, atol=1e-12), f"N has off-diagonal elements: {off_diag}"


def test_c1_diagonal_N_matches_analytic():
    """C1/J6: Diagonal elements of N match the per-axis analytic formula."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    eps_i = np.diag([eps_x, eps_y, eps_z]).astype(complex)
    N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
    N_analytic = _analytic_N_diagonal_sphere(eps_x, eps_y, eps_z, eps_e)
    assert np.allclose(N, N_analytic, atol=1e-12), \
        f"N={np.diag(N)} != analytic={np.diag(N_analytic)}"


def test_c1_R_particle_off_diagonal_remains_zero():
    """C1: Off-diagonal R components remain zero when ε_i and R are both diagonal."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    eps_i = np.diag([eps_x, eps_y, eps_z]).astype(complex)
    R = np.diag([1.0, 2.0, 0.5])
    N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i, eps_e)
    off_diag = R_p - np.diag(np.diag(R_p))
    assert np.allclose(off_diag, 0, atol=1e-12)


def test_c1_diagonal_R_particle_components():
    """C1: Each diagonal element of R_particle satisfies R_p_jj = N_j^2*(1-...) R_jj."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    eps_i_diag = np.array([eps_x, eps_y, eps_z])
    eps_i = np.diag(eps_i_diag).astype(complex)
    R = np.diag([1.0, 2.0, 0.5])
    N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i, eps_e)

    for j, eps_j in enumerate(eps_i_diag):
        N_j = 3.0 * eps_e / (eps_j + 2.0 * eps_e)
        # For diagonal case, R_p_jj = N_j^3 * R_jj
        expected = N_j**3 * R[j, j]
        assert abs(R_p[j, j] - expected) < 1e-12, \
            f"j={j}: R_p[{j},{j}]={R_p[j,j]}, expected {expected}"


def test_c1_isotropic_eps_recovers_sphere_result():
    """C1: ε_x = ε_y = ε_z should recover the isotropic sphere result."""
    eps_scalar = 3.0
    eps_e = 1.0
    eps_i_aniso = np.diag([eps_scalar, eps_scalar, eps_scalar]).astype(complex)
    eps_i_iso = eps_scalar * np.eye(3, dtype=complex)
    R = np.diag([2.0, 1.0, 0.5])

    N_aniso = compute_internal_field_tensor(L_sphere, eps_i_aniso, eps_e)
    N_iso = compute_internal_field_tensor(L_sphere, eps_i_iso, eps_e)
    assert np.allclose(N_aniso, N_iso, atol=1e-12)

    R_p_aniso = compute_particle_raman_tensor(R, N_aniso, L_sphere, eps_i_aniso, eps_e)
    R_p_iso = compute_particle_raman_tensor(R, N_iso, L_sphere, eps_i_iso, eps_e)
    assert np.allclose(R_p_aniso, R_p_iso, atol=1e-12)


# ---------------------------------------------------------------------------
# C2 / J8 / M6 — Rotation invariance of powder average
# ---------------------------------------------------------------------------

def test_c2_powder_average_invariant_under_global_rotation():
    """C2/J8/M6: Powder-averaged VV and VH are unchanged when R and ε_i are co-rotated."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    eps_i_crystal = np.diag([eps_x, eps_y, eps_z]).astype(complex)
    R_crystal = np.diag([2.0, 1.0, 0.5])

    # Powder average in crystal frame (using analytical invariants for sphere)
    N_crystal = compute_internal_field_tensor(L_sphere, eps_i_crystal, eps_e)
    R_p_crystal = compute_particle_raman_tensor(R_crystal, N_crystal, L_sphere, eps_i_crystal, eps_e)
    vv_ref, vh_ref = compute_powder_raman_intensities(R_p_crystal)

    # Apply a fixed arbitrary rotation to both tensors
    rng = np.random.default_rng(99)
    for _ in range(8):
        G0 = _random_rotation(rng)
        eps_i_rot = G0 @ eps_i_crystal @ G0.T
        R_rot = G0 @ R_crystal @ G0.T
        N_rot = compute_internal_field_tensor(L_sphere, eps_i_rot, eps_e)
        R_p_rot = compute_particle_raman_tensor(R_rot, N_rot, L_sphere, eps_i_rot, eps_e)
        vv, vh = compute_powder_raman_intensities(R_p_rot)
        assert abs(vv - vv_ref) < 1e-8, f"VV changed after rotation: {vv} vs {vv_ref}"
        assert abs(vh - vh_ref) < 1e-8, f"VH changed after rotation: {vh} vs {vh_ref}"


def test_c2_rotation_invariance_off_diagonal_R():
    """C2/J8: Rotation invariance holds for an off-diagonal Raman tensor with anisotropic ε."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    eps_i_crystal = np.diag([eps_x, eps_y, eps_z]).astype(complex)
    R_crystal = np.array([[0.0, 1.0, 0.0],
                          [1.0, 0.0, 0.3],
                          [0.0, 0.3, 0.0]])
    N_crystal = compute_internal_field_tensor(L_sphere, eps_i_crystal, eps_e)
    R_p_crystal = compute_particle_raman_tensor(R_crystal, N_crystal, L_sphere, eps_i_crystal, eps_e)
    vv_ref, vh_ref = compute_powder_raman_intensities(R_p_crystal)

    rng = np.random.default_rng(77)
    for _ in range(8):
        G0 = _random_rotation(rng)
        eps_i_rot = G0 @ eps_i_crystal @ G0.T
        R_rot = G0 @ R_crystal @ G0.T
        N_rot = compute_internal_field_tensor(L_sphere, eps_i_rot, eps_e)
        R_p_rot = compute_particle_raman_tensor(R_rot, N_rot, L_sphere, eps_i_rot, eps_e)
        vv, vh = compute_powder_raman_intensities(R_p_rot)
        assert abs(vv - vv_ref) < 1e-8, "VV not rotation-invariant for off-diagonal R"
        assert abs(vh - vh_ref) < 1e-8, "VH not rotation-invariant for off-diagonal R"


# ---------------------------------------------------------------------------
# C3 — Complex symmetric off-diagonal ε_i
# ---------------------------------------------------------------------------

def test_c3_complex_offdiagonal_eps_reduces_to_diagonal_when_delta_zero():
    """C3: Off-diagonal ε reduces to diagonal result when the off-diagonal element is zero."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    delta = 0.0 + 0.0j
    eps_i_diag = np.diag([eps_x, eps_y, eps_z]).astype(complex)
    eps_i_offdiag = np.array([[eps_x, delta, 0.0],
                              [delta, eps_y, 0.0],
                              [0.0, 0.0, eps_z]], dtype=complex)
    N_diag = compute_internal_field_tensor(L_sphere, eps_i_diag, eps_e)
    N_offdiag = compute_internal_field_tensor(L_sphere, eps_i_offdiag, eps_e)
    assert np.allclose(N_diag, N_offdiag, atol=1e-12)


def test_c3_complex_offdiagonal_eps_varies_smoothly():
    """C3: N varies smoothly as the off-diagonal element increases from zero."""
    eps_x, eps_y, eps_z, eps_e = 2.0, 4.0, 6.0, 1.0
    N_prev = None
    for delta in [0.0, 0.01, 0.05, 0.1, 0.2]:
        eps_i = np.array([[eps_x, delta, 0.0],
                          [delta, eps_y, 0.0],
                          [0.0, 0.0, eps_z]], dtype=complex)
        N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
        assert np.all(np.isfinite(N)), f"N has non-finite elements at delta={delta}"
        if N_prev is not None:
            diff = np.max(np.abs(N - N_prev))
            assert diff < 0.1, f"N changed abruptly: delta={delta}, diff={diff}"
        N_prev = N.copy()


def test_c3_complex_eps_full_tensor_inversion_matches_direct():
    """C3: N via tensor inversion agrees with direct computation for complex symmetric ε_i."""
    eps_x, eps_y, eps_z, eps_e = 2.0 + 0.1j, 4.0 + 0.2j, 6.0 + 0.3j, 1.0
    delta = 0.5 + 0.1j
    eps_i = np.array([[eps_x, delta, 0.0],
                      [delta, eps_y, 0.0],
                      [0.0, 0.0, eps_z]], dtype=complex)
    # compute_N uses matrix inversion; verify against explicit formula
    I3 = np.eye(3, dtype=complex)
    N_direct = np.linalg.inv(I3 + (1.0 / eps_e) * L_sphere @ (eps_i - eps_e * I3))
    N_helper = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
    assert np.allclose(N_direct, N_helper, atol=1e-12)


# ---------------------------------------------------------------------------
# K4 — Isotropic R with anisotropic ε_i can produce anisotropic effective tensor
# ---------------------------------------------------------------------------

def test_k4_isotropic_R_anisotropic_eps_makes_anisotropic_R_particle():
    """K4: R = a*I with anisotropic ε_i gives non-zero VH in the particle Raman tensor."""
    a = 1.0
    R = a * np.eye(3)
    eps_x, eps_y, eps_z, eps_e = 2.0, 6.0, 12.0, 1.0
    eps_i = np.diag([eps_x, eps_y, eps_z]).astype(complex)
    N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i, eps_e)
    vv, vh = compute_powder_raman_intensities(R_p)
    # With strongly anisotropic ε, VH should be non-zero
    assert vh > 0.0, "Expected non-zero VH for isotropic R but anisotropic ε"


def test_k4_isotropic_eps_gives_zero_VH_for_isotropic_R():
    """K4: R = a*I with isotropic ε_i gives zero VH (no anisotropy)."""
    a = 2.0
    R = a * np.eye(3)
    eps_scalar = 4.0
    eps_i = eps_scalar * np.eye(3, dtype=complex)
    eps_e = 1.0
    N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
    R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i, eps_e)
    _, vh = compute_powder_raman_intensities(R_p)
    assert abs(vh) < 1e-12, f"VH={vh}, expected 0 for isotropic ε and R"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
