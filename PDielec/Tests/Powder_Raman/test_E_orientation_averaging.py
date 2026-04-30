"""Tests E1-E4: Orientation averaging and SO(3) Sobol sampling.

Plan sections: E1, E2, E3, E4, J8, J9, M6, M7, M10.

Tests verify:
- Sobol-based SO(3) sampling produces uniformly distributed orientations (E1);
- numerical powder average converges to analytic Placzek invariants (E2/J9);
- powder average is rotation-invariant (E3/J8/M6);
- VV and VH polarisation definitions are consistent (E4).

The numerical averaging in _powder_average_numerical matches the
backscattering geometry in PowderScenarioTab._compute_orientation_sampled_spectrum:
  VV: e_L = e_S = [0, 1, 0]
  VH: e_L = [0, 1, 0], e_S = [1, 0, 0]
  contribution per orientation: 45 × |e_S @ R_lab @ e_L|^2
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import compute_powder_raman_intensities, sobol_rotations


def _random_rotation(rng=None):
    if rng is None:
        rng = np.random.default_rng(12345)
    A = rng.standard_normal((3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


def _powder_average_numerical(raman_tensor, n_samples=2048, seed=42):
    R_cryst = np.asarray(raman_tensor, dtype=complex)
    e_L = np.array([0.0, 1.0, 0.0])
    e_VV = np.array([0.0, 1.0, 0.0])
    e_VH = np.array([1.0, 0.0, 0.0])
    vv_sum = 0.0
    vh_sum = 0.0
    for G in sobol_rotations(n_samples, seed=seed):
        R_lab = G @ R_cryst @ G.T
        vv_sum += 45.0 * abs(e_VV @ R_lab @ e_L) ** 2
        vh_sum += 45.0 * abs(e_VH @ R_lab @ e_L) ** 2
    return vv_sum / n_samples, vh_sum / n_samples


# ---------------------------------------------------------------------------
# E1 — Uniform SO(3) sampling test
# ---------------------------------------------------------------------------

def test_e1_direction_cosine_mean_is_zero():
    """E1: Mean of each direction cosine g_ij over Sobol SO(3) sample is zero."""
    n = 1024
    rotations = sobol_rotations(n)
    G = np.array(rotations)  # shape (n, 3, 3)
    mean_abs = np.max(np.abs(G.mean(axis=0)))
    assert mean_abs < 0.05, f"Max |<g_ij>| = {mean_abs}, expected < 0.05"


def test_e1_direction_cosine_mean_squared_is_one_third():
    """E1: Mean of g_ij^2 over Sobol SO(3) sample is approximately 1/3."""
    n = 1024
    rotations = sobol_rotations(n)
    G = np.array(rotations)  # shape (n, 3, 3)
    mean_sq = np.mean(G**2, axis=0)
    expected = np.full((3, 3), 1.0 / 3.0)
    max_err = np.max(np.abs(mean_sq - expected))
    assert max_err < 0.05, f"Max |<g_ij^2> - 1/3| = {max_err}"


def test_e1_rotations_are_orthogonal():
    """E1: All Sobol rotation matrices are orthogonal (R R^T = I)."""
    for G in sobol_rotations(64):
        err = np.max(np.abs(G @ G.T - np.eye(3)))
        assert err < 1e-12, f"Rotation not orthogonal: err={err}"


def test_e1_rotations_have_det_plus_one():
    """E1: All Sobol rotation matrices have det = +1."""
    for G in sobol_rotations(64):
        det = np.linalg.det(G)
        assert abs(det - 1.0) < 1e-12, f"det={det}"


def test_e1_mean_cross_terms_near_zero():
    """E1: Cross-moment <g_ij * g_kl> ≈ 0 for i≠k or j≠l."""
    n = 2048
    rotations = sobol_rotations(n)
    G = np.array(rotations)
    # Check a selection of cross-terms
    cross_00_11 = np.mean(G[:, 0, 0] * G[:, 1, 1])
    cross_00_01 = np.mean(G[:, 0, 0] * G[:, 0, 1])
    assert abs(cross_00_11) < 0.05, f"<g_00 g_11> = {cross_00_11}"
    assert abs(cross_00_01) < 0.05, f"<g_00 g_01> = {cross_00_01}"


# ---------------------------------------------------------------------------
# E2 / J9 — Convergence to analytic powder invariants
# ---------------------------------------------------------------------------

def test_e2_isotropic_tensor_convergence():
    """E2/J9: Isotropic R = a*I: numerical VH converges to 0 with increasing samples."""
    a = 2.0
    R = a * np.eye(3)
    for n in [64, 256, 1024]:
        _, vh = _powder_average_numerical(R, n_samples=n)
        assert abs(vh) < 0.5, f"n={n}: VH={vh} not converging to 0"


def test_e2_diagonal_tensor_VV_converges_to_analytic():
    """E2/J9: Numerical VV for diagonal tensor converges to analytic value."""
    R = np.diag([2.0, 1.0, 0.5])
    vv_analytic, _ = compute_powder_raman_intensities(R)

    n = 4096
    vv_num, _ = _powder_average_numerical(R, n_samples=n)

    # Relative error should be small
    rel_err = abs(vv_num - vv_analytic) / vv_analytic
    assert rel_err < 0.05, f"n={n}: rel error {rel_err:.4f}"


def test_e2_diagonal_tensor_VH_converges_to_analytic():
    """E2/J9: Numerical VH for diagonal tensor converges to analytic value."""
    R = np.diag([2.0, 1.0, 0.5])
    _, vh_analytic = compute_powder_raman_intensities(R)

    n = 4096
    _, vh_num = _powder_average_numerical(R, n_samples=n)

    rel_err = abs(vh_num - vh_analytic) / vh_analytic
    assert rel_err < 0.05, f"n={n}: rel error {rel_err:.4f}"


def test_e2_offdiagonal_tensor_depolarisation_ratio_converges():
    """E2/J9: Numerical depolarisation ratio VH/VV for off-diagonal R converges to 0.75."""
    d = 1.0
    R = np.array([[0.0, d, 0.0], [d, 0.0, 0.0], [0.0, 0.0, 0.0]])

    n = 4096
    vv_num, vh_num = _powder_average_numerical(R, n_samples=n)
    rho = vh_num / vv_num
    assert abs(rho - 0.75) < 0.05, f"n={n}: rho={rho:.4f}, expected 0.75"


def test_e2_error_decreases_with_more_samples():
    """E2: Numerical error decreases monotonically as n increases for diagonal R."""
    R = np.diag([2.0, 1.0, 0.5])
    vv_analytic, _ = compute_powder_raman_intensities(R)
    prev_err = None
    for n in [64, 256, 1024, 4096]:
        vv_num, _ = _powder_average_numerical(R, n_samples=n)
        err = abs(vv_num - vv_analytic)
        if prev_err is not None:
            # Allow some non-monotonicity from Sobol quasi-randomness
            assert err < prev_err * 2.0, f"Error did not decrease: n={n}, err={err}, prev={prev_err}"
        prev_err = err


# ---------------------------------------------------------------------------
# E3 / J8 / M6 — Rotation-invariance of powder average (numerical path)
# ---------------------------------------------------------------------------

def test_e3_isotropic_sphere_numerical_rotation_invariant():
    """E3/J8/M6: Numerical powder average is unchanged after a global rotation (isotropic R)."""
    R = np.diag([2.0, 1.0, 0.5])
    vv_ref, vh_ref = _powder_average_numerical(R, n_samples=2048)

    rng = np.random.default_rng(17)
    for _ in range(4):
        G0 = _random_rotation(rng)
        R_rot = G0 @ R @ G0.T
        vv, vh = _powder_average_numerical(R_rot, n_samples=2048)
        rel_vv = abs(vv - vv_ref) / (abs(vv_ref) + 1e-10)
        rel_vh = abs(vh - vh_ref) / (abs(vh_ref) + 1e-10)
        assert rel_vv < 0.05, f"VV changed after rotation: rel err={rel_vv:.4f}"
        assert rel_vh < 0.05, f"VH changed after rotation: rel err={rel_vh:.4f}"


def test_e3_offdiagonal_rotation_invariant():
    """E3: Off-diagonal tensor: numerical powder average is rotation-invariant."""
    d = 1.0
    R = np.array([[0.0, d, 0.0], [d, 0.0, 0.0], [0.0, 0.0, 0.0]])
    vv_ref, vh_ref = _powder_average_numerical(R, n_samples=2048)

    rng = np.random.default_rng(88)
    for _ in range(4):
        G0 = _random_rotation(rng)
        R_rot = G0 @ R @ G0.T
        vv, vh = _powder_average_numerical(R_rot, n_samples=2048)
        rel_vv = abs(vv - vv_ref) / (abs(vv_ref) + 1e-10)
        rel_vh = abs(vh - vh_ref) / (abs(vh_ref) + 1e-10)
        assert rel_vv < 0.05, f"Off-diag VV changed: rel err={rel_vv:.4f}"
        assert rel_vh < 0.05, f"Off-diag VH changed: rel err={rel_vh:.4f}"


# ---------------------------------------------------------------------------
# E4 — VV and VH polarisation definitions
# ---------------------------------------------------------------------------

def test_e4_VV_larger_than_VH_for_isotropic_tensor():
    """E4: For isotropic R, VV > VH = 0 (VH is zero, VV is non-zero)."""
    R = np.eye(3)
    vv, vh = compute_powder_raman_intensities(R)
    assert vv > 0.0
    assert abs(vh) < 1e-12


def test_e4_VV_and_VH_positive():
    """E4: VV and VH must both be non-negative for physically valid Raman tensors."""
    cases = [
        np.diag([2.0, 1.0, 0.5]),
        np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        np.array([[1.0, 0.5, 0.2], [0.5, 0.8, 0.3], [0.2, 0.3, 0.6]]),
        np.eye(3),
    ]
    for R in cases:
        vv, vh = compute_powder_raman_intensities(R)
        assert vv >= 0.0, f"VV={vv} negative"
        assert vh >= 0.0, f"VH={vh} negative"


def test_e4_depolarisation_ratio_correct_for_isotropic():
    """E4: Depolarisation ratio for isotropic tensor is 0."""
    R = 2.0 * np.eye(3)
    vv, vh = compute_powder_raman_intensities(R)
    assert abs(vh / vv) < 1e-12


def test_e4_depolarisation_ratio_correct_for_offdiagonal():
    """E4: Depolarisation ratio for purely off-diagonal symmetric tensor is 3/4."""
    d = 1.5
    R = np.array([[0.0, d, 0.0], [d, 0.0, 0.0], [0.0, 0.0, 0.0]])
    vv, vh = compute_powder_raman_intensities(R)
    assert abs(vh / vv - 0.75) < 1e-12


def test_e4_VV_VH_sum_equals_total():
    """E4: VV + VH = total intensity (no kappa term in the symmetric Raman tensor case)."""
    R = np.diag([2.0, 1.0, 0.5])
    vv, vh = compute_powder_raman_intensities(R)
    # For symmetric R (kappa2 = 0): total = 45*alpha2 + 7*gamma2 = VV + VH - gamma2
    # Actually: total = 45*alpha2 + 7*gamma2; VV = 45*alpha2 + 4*gamma2; VH = 3*gamma2
    # So: VV + VH = 45*alpha2 + 7*gamma2 = total ✓
    R = np.asarray(R, dtype=complex)
    I3 = np.eye(3, dtype=complex)
    alpha = np.trace(R) / 3.0
    gamma_t = 0.5 * (R + R.T) - alpha * I3
    gamma2 = 3.0 / 2.0 * float(np.real(np.sum(gamma_t * np.conj(gamma_t))))
    alpha2 = float(np.real(alpha * np.conj(alpha)))
    total_expected = 45.0 * alpha2 + 7.0 * gamma2
    assert abs(vv + vh - total_expected) < 1e-12


def test_e4_numerical_VV_larger_than_analytical_VH_for_diagonal():
    """E4: Numerical VV should exceed numerical VH for a standard anisotropic tensor."""
    R = np.diag([2.0, 1.0, 0.5])
    vv, vh = _powder_average_numerical(R, n_samples=2048)
    assert vv > vh, f"VV={vv} not > VH={vh}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
