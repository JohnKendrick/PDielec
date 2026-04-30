"""Tests A1-A4: Raman tensor algebra and powder invariants without local-field effects.

Plan sections: A1, A2, A3, A4, J1, J2, J3, M1, M2, M3.

All tests use N = I (no local-field correction), so R_particle = R_eps.
The analytical invariant formulas tested here are those in Eqs. 90-99 of
Raman-Theory.pdf, as implemented in PowderScenarioTab._calculate_raman.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import compute_powder_raman_intensities, sobol_rotations


# ---------------------------------------------------------------------------
# A1 / J1 / M1 — Isotropic tensor R = a*I
# ---------------------------------------------------------------------------

def test_a1_isotropic_vv_equals_45_a_squared():
    """A1/J1: VV intensity for R = a*I must equal 45 a^2 (Eqs. 90-99)."""
    a = 3.7
    R = a * np.eye(3)
    vv, vh = compute_powder_raman_intensities(R)
    assert abs(vv - 45.0 * a**2) < 1e-10, f"VV={vv}, expected {45*a**2}"


def test_a1_isotropic_vh_is_zero():
    """A1/J1: VH intensity for R = a*I must be zero."""
    a = 2.5
    R = a * np.eye(3)
    vv, vh = compute_powder_raman_intensities(R)
    assert abs(vh) < 1e-12, f"VH={vh}, expected 0"


def test_a1_isotropic_rotation_invariant_vv():
    """A1: VV must be unchanged when R = a*I is rotated by any orientation."""
    a = 1.5
    R = a * np.eye(3)
    vv_ref, _ = compute_powder_raman_intensities(R)
    for G in sobol_rotations(64):
        R_rot = G @ R @ G.T
        vv, vh = compute_powder_raman_intensities(R_rot)
        assert abs(vv - vv_ref) < 1e-10, "VV changed under rotation for isotropic R"
        assert abs(vh) < 1e-11, "VH non-zero under rotation for isotropic R"


def test_a1_isotropic_scaling():
    """A1: VV scales as a^2 (doubling a gives 4x VV)."""
    vv1, _ = compute_powder_raman_intensities(1.0 * np.eye(3))
    vv2, _ = compute_powder_raman_intensities(2.0 * np.eye(3))
    assert abs(vv2 / vv1 - 4.0) < 1e-10


def test_a1_negative_isotropic():
    """A1: Negative a gives same VV as positive a (intensity is \|a\|^2-based)."""
    a = 2.0
    vv_pos, _ = compute_powder_raman_intensities(a * np.eye(3))
    vv_neg, _ = compute_powder_raman_intensities(-a * np.eye(3))
    assert abs(vv_pos - vv_neg) < 1e-10


# ---------------------------------------------------------------------------
# A2 / J2 / M2 — Diagonal anisotropic tensor R = diag(a, b, c)
# ---------------------------------------------------------------------------

def _placzek_invariants_diagonal(a, b, c):
    """Exact Placzek invariants for a real diagonal tensor diag(a, b, c).

    Returns (vv_analytic, vh_analytic) using the code's convention.

    Code convention (Eq. 90-99):
        alpha  = (a + b + c) / 3
        gamma2 = 3/2(a-alpha)^2 + (b-alpha)^2 + (c-alpha)^2  (traceless-sym sum-of-squares)
        VV     = 45*alpha^2 + 4*gamma2
        VH     = 3*gamma2
    """
    alpha = (a + b + c) / 3.0
    gamma2 = 3.0 / 2.0 * ((a - alpha)**2 + (b - alpha)**2 + (c - alpha)**2)
    vv = 45.0 * alpha**2 + 4.0 * gamma2
    vh = 3.0 * gamma2
    return vv, vh


def test_a2_diagonal_tensor_vv_matches_analytic():
    """A2/J2: Diagonal R = diag(2, 1, 0.5) VV matches analytic invariant formula."""
    R = np.diag([2.0, 1.0, 0.5])
    vv, vh = compute_powder_raman_intensities(R)
    vv_ref, vh_ref = _placzek_invariants_diagonal(2.0, 1.0, 0.5)
    assert abs(vv - vv_ref) < 1e-10, f"VV={vv}, expected {vv_ref}"


def test_a2_diagonal_tensor_vh_matches_analytic():
    """A2/J2: Diagonal R = diag(2, 1, 0.5) VH matches analytic invariant formula."""
    R = np.diag([2.0, 1.0, 0.5])
    vv, vh = compute_powder_raman_intensities(R)
    vv_ref, vh_ref = _placzek_invariants_diagonal(2.0, 1.0, 0.5)
    assert abs(vh - vh_ref) < 1e-10, f"VH={vh}, expected {vh_ref}"


def test_a2_isotropic_special_case():
    """A2: a=b=c should give zero VH (isotropic limit)."""
    R = np.diag([2.0, 2.0, 2.0])
    vv, vh = compute_powder_raman_intensities(R)
    vv_ref, vh_ref = _placzek_invariants_diagonal(2.0, 2.0, 2.0)
    assert abs(vh) < 1e-12
    assert abs(vv - vv_ref) < 1e-10


def test_a2_uniaxial_special_case():
    """A2: a=b != c (uniaxial tensor) VH is non-zero and matches analytic."""
    R = np.diag([1.0, 1.0, 3.0])
    vv, vh = compute_powder_raman_intensities(R)
    vv_ref, vh_ref = _placzek_invariants_diagonal(1.0, 1.0, 3.0)
    assert abs(vv - vv_ref) < 1e-10
    assert abs(vh - vh_ref) < 1e-10
    assert vh > 0.0


def test_a2_several_diagonal_tensors():
    """A2: Multiple diagonal tensors all match analytic formula."""
    cases = [
        (1.0, 2.0, 3.0),
        (0.5, 0.5, 2.0),
        (1.0, -0.5, 0.3),
        (0.0, 0.0, 1.0),
    ]
    for a, b, c in cases:
        R = np.diag([a, b, c])
        vv, vh = compute_powder_raman_intensities(R)
        vv_ref, vh_ref = _placzek_invariants_diagonal(a, b, c)
        assert abs(vv - vv_ref) < 1e-10, f"diag({a},{b},{c}) VV mismatch"
        assert abs(vh - vh_ref) < 1e-10, f"diag({a},{b},{c}) VH mismatch"


def test_a2_rotation_invariance_diagonal():
    """A2: Rotating a diagonal tensor before powder averaging must not change VV/VH."""
    R = np.diag([2.0, 1.0, 0.5])
    vv_ref, vh_ref = compute_powder_raman_intensities(R)
    for G in sobol_rotations(32):
        R_rot = G @ R @ G.T
        vv, vh = compute_powder_raman_intensities(R_rot)
        assert abs(vv - vv_ref) < 1e-8, "VV changed under rotation"
        assert abs(vh - vh_ref) < 1e-8, "VH changed under rotation"


# ---------------------------------------------------------------------------
# A3 / J3 / M3 — Off-diagonal symmetric tensor
# ---------------------------------------------------------------------------

def test_a3_offdiagonal_vh_over_vv_equals_0p75():
    """A3/J3/M3: Purely off-diagonal R_xy: depolarisation ratio VH/VV = 0.75."""
    d = 2.0
    R = np.array([[0.0, d, 0.0],
                  [d,  0.0, 0.0],
                  [0.0, 0.0, 0.0]])
    vv, vh = compute_powder_raman_intensities(R)
    rho = vh / vv
    assert abs(rho - 0.75) < 1e-12, f"rho={rho}, expected 0.75"


def test_a3_offdiagonal_alpha_is_zero():
    """A3: Isotropic invariant alpha must be zero for a traceless off-diagonal tensor."""
    d = 1.5
    R = np.array([[0.0, d, 0.0],
                  [d,  0.0, 0.0],
                  [0.0, 0.0, 0.0]])
    alpha = np.trace(R) / 3.0
    assert abs(alpha) < 1e-12


def test_a3_offdiagonal_rotation_invariant():
    """A3: Rotating the off-diagonal tensor before powder averaging must not change VV/VH."""
    d = 1.0
    R = np.array([[0.0, d, 0.0],
                  [d,  0.0, 0.0],
                  [0.0, 0.0, 0.0]])
    vv_ref, vh_ref = compute_powder_raman_intensities(R)
    for G in sobol_rotations(32):
        R_rot = G @ R @ G.T
        vv, vh = compute_powder_raman_intensities(R_rot)
        assert abs(vv - vv_ref) < 1e-8, "VV changed under rotation for off-diagonal R"
        assert abs(vh - vh_ref) < 1e-8, "VH changed under rotation for off-diagonal R"


def test_a3_depolarisation_ratio_scale_invariant():
    """A3: Depolarisation ratio VH/VV is independent of tensor amplitude d."""
    for d in [0.5, 1.0, 2.0, 5.0]:
        R = np.array([[0.0, d, 0.0], [d, 0.0, 0.0], [0.0, 0.0, 0.0]])
        vv, vh = compute_powder_raman_intensities(R)
        rho = vh / vv
        assert abs(rho - 0.75) < 1e-12, f"d={d}: rho={rho}"


def test_a3_depolarisation_ratio_all_off_diagonal():
    """A3: All three off-diagonal symmetric tensors give depolarisation ratio 0.75."""
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        R = np.zeros((3, 3))
        R[i, j] = 1.0
        R[j, i] = 1.0
        vv, vh = compute_powder_raman_intensities(R)
        rho = vh / vv
        assert abs(rho - 0.75) < 1e-12, f"R[{i},{j}]+R[{j},{i}]: rho={rho}"


# ---------------------------------------------------------------------------
# A4 — Antisymmetric tensor (kappa term)
# ---------------------------------------------------------------------------

def test_a4_antisymmetric_vh_nonzero():
    """A4: Antisymmetric R_xy tensor contributes to both VV and VH via kappa^2 term."""
    d = 1.0
    R = np.array([[ 0.0, d,  0.0],
                  [-d,  0.0, 0.0],
                  [ 0.0, 0.0, 0.0]])
    vv, vh = compute_powder_raman_intensities(R)
    assert vv > 0.0, "VV should be non-zero for antisymmetric R"
    assert vh > 0.0, "VH should be non-zero for antisymmetric R"


def test_a4_antisymmetric_alpha_is_zero():
    """A4: Antisymmetric tensor has zero trace, so alpha = 0."""
    d = 2.0
    R = np.array([[ 0.0, d,  0.0],
                  [-d,  0.0, 0.0],
                  [ 0.0, 0.0, 0.0]])
    alpha = np.trace(R) / 3.0
    assert abs(alpha) < 1e-12


def test_a4_antisymmetric_kappa_contribution():
    """A4: Antisymmetric R contributes via the 5*kappa^2 term, not gamma^2."""
    d = 1.0
    R_anti = np.array([[ 0.0, d,  0.0],
                       [-d,  0.0, 0.0],
                       [ 0.0, 0.0, 0.0]])
    R_sym  = np.array([[0.0, d, 0.0],
                       [d,  0.0, 0.0],
                       [0.0, 0.0, 0.0]])
    vv_anti, vh_anti = compute_powder_raman_intensities(R_anti)
    vv_sym, vh_sym   = compute_powder_raman_intensities(R_sym)
    # Antisymmetric: gamma2=0, kappa2=(3/2)*2=3 => VV=VH=5*kappa2=15*d^2
    # Symmetric off-diagonal: kappa2=0, gamma2=(3/2)*2=3 => VV=4*gamma2=12*d^2, VH=3*gamma2=9*d^2
    assert abs(vv_anti - 15.0 * d**2) < 1e-12, f"antisymmetric VV={vv_anti}"
    assert abs(vh_anti - 15.0 * d**2) < 1e-12, f"antisymmetric VH={vh_anti}"
    assert abs(vv_sym - 12.0 * d**2) < 1e-12, f"symmetric VV={vv_sym}"
    assert abs(vh_sym - 9.0 * d**2) < 1e-12, f"symmetric VH={vh_sym}"


# ---------------------------------------------------------------------------
# Complex phase invariance (M-plan: phase of R does not change intensity)
# ---------------------------------------------------------------------------

def test_global_complex_phase_invariance():
    """Plan M: Multiplying all Raman tensors by a global complex phase leaves intensity unchanged."""
    R = np.array([[2.0, 0.5, 0.1],
                  [0.5, 1.0, 0.3],
                  [0.1, 0.3, 0.5]])
    vv_ref, vh_ref = compute_powder_raman_intensities(R)
    for phase in [1j, np.exp(1j * 0.7), np.exp(1j * 1.23), -1.0 + 0.0j]:
        R_phased = phase * R
        vv, vh = compute_powder_raman_intensities(R_phased)
        assert abs(vv - vv_ref) < 1e-10, f"phase={phase}: VV changed"
        assert abs(vh - vh_ref) < 1e-10, f"phase={phase}: VH changed"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
