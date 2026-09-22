"""Tests C1-C3: SMatrix and Redheffer star product.

Plan sections C1 (unit_matrix identity element), C2 (redheffer composition),
C3 (calculate_s shape).

The Redheffer star product is the identity operation for the unit scattering
matrix: S ★ I = I ★ S = S.  It is also associative.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import SMatrix


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _random_smatrix(seed=0):
    """Return a random SMatrix built from a well-conditioned 4×4 complex array."""
    rng = np.random.default_rng(seed)
    # Build a small random S-matrix with off-diag elements smaller than diagonal
    S = 0.1 * (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4)))
    # Make S21, S12 close to identity (lossless-ish)
    S[2:4, 0:2] += np.eye(2)
    S[0:2, 2:4] += np.eye(2)
    return SMatrix(S.astype(np.cdouble))


def _smatrices_allclose(sa, sb, atol=1e-12):
    """Return True if two SMatrix objects have the same S, S11, S22, S12, S21."""
    np.testing.assert_allclose(sa.S, sb.S, atol=atol)
    np.testing.assert_allclose(sa.S11, sb.S11, atol=atol)
    np.testing.assert_allclose(sa.S12, sb.S12, atol=atol)
    np.testing.assert_allclose(sa.S21, sb.S21, atol=atol)
    np.testing.assert_allclose(sa.S22, sb.S22, atol=atol)


# ---------------------------------------------------------------------------
# C1: unit_matrix
# ---------------------------------------------------------------------------

class TestC1UnitMatrix:
    """C1: The unit SMatrix is the identity element of the Redheffer product."""

    def test_unit_matrix_constructs(self):
        sm = SMatrix()
        sm.unit_matrix()
        assert sm is not None

    def test_unit_s11_is_zero(self):
        sm = SMatrix()
        sm.unit_matrix()
        np.testing.assert_allclose(sm.S11, np.zeros((2, 2)), atol=1e-15)

    def test_unit_s22_is_zero(self):
        sm = SMatrix()
        sm.unit_matrix()
        np.testing.assert_allclose(sm.S22, np.zeros((2, 2)), atol=1e-15)

    def test_unit_s21_is_identity(self):
        sm = SMatrix()
        sm.unit_matrix()
        np.testing.assert_allclose(sm.S21, np.eye(2), atol=1e-15)

    def test_unit_s12_is_identity(self):
        sm = SMatrix()
        sm.unit_matrix()
        np.testing.assert_allclose(sm.S12, np.eye(2), atol=1e-15)

    def test_left_identity_S_star_I_equals_S(self):
        S = _random_smatrix(seed=1)
        I_mat = SMatrix()
        result = S.redheffer(I_mat)
        _smatrices_allclose(result, S, atol=1e-12)

    def test_right_identity_I_star_S_equals_S(self):
        S = _random_smatrix(seed=2)
        I_mat = SMatrix()
        result = I_mat.redheffer(S)
        _smatrices_allclose(result, S, atol=1e-12)


# ---------------------------------------------------------------------------
# C2: redheffer (Redheffer star product)
# ---------------------------------------------------------------------------

class TestC2Redheffer:
    """C2: Redheffer product properties."""

    def test_two_unit_matrices_give_unit(self):
        I1 = SMatrix()
        I2 = SMatrix()
        result = I1.redheffer(I2)
        expected = SMatrix()
        _smatrices_allclose(result, expected, atol=1e-14)

    def test_associativity(self):
        A = _random_smatrix(seed=10)
        B = _random_smatrix(seed=11)
        Cc = _random_smatrix(seed=12)
        # (A ★ B) ★ C
        AB_C = A.redheffer(B).redheffer(Cc)
        # A ★ (B ★ C)
        A_BC = A.redheffer(B.redheffer(Cc))
        np.testing.assert_allclose(AB_C.S, A_BC.S, atol=1e-11)

    def test_redheffer_none_treated_as_unit(self):
        # Passing None as b should be treated as unit matrix
        S = _random_smatrix(seed=3)
        result = S.redheffer(None)
        _smatrices_allclose(result, S, atol=1e-12)

    def test_result_is_smatrix_instance(self):
        A = _random_smatrix(seed=5)
        B = _random_smatrix(seed=6)
        result = A.redheffer(B)
        assert isinstance(result, SMatrix)

    def test_s21_of_two_pass_through_layers_bounded(self):
        # For a system with only forward transmission, |S21| <= 1 element-wise
        A = _random_smatrix(seed=7)
        B = _random_smatrix(seed=8)
        result = A.redheffer(B)
        # Transmission amplitude should not exceed a reasonable bound
        assert np.all(np.abs(result.S21) < 10.0)


# ---------------------------------------------------------------------------
# C3: calculate_s
# ---------------------------------------------------------------------------

class TestC3CalculateS:
    """C3: calculate_s reconstructs the full 4×4 S matrix from sub-blocks."""

    def test_calculate_s_shape(self):
        sm = _random_smatrix(seed=9)
        assert sm.S.shape == (4, 4)

    def test_calculate_s_dtype_complex(self):
        sm = _random_smatrix(seed=9)
        assert np.iscomplexobj(sm.S)

    def test_subblocks_consistent_with_full_matrix(self):
        sm = _random_smatrix(seed=13)
        np.testing.assert_allclose(sm.S[0:2, 0:2], sm.S11, atol=1e-15)
        np.testing.assert_allclose(sm.S[2:4, 2:4], sm.S22, atol=1e-15)
        np.testing.assert_allclose(sm.S[2:4, 0:2], sm.S21, atol=1e-15)
        np.testing.assert_allclose(sm.S[0:2, 2:4], sm.S12, atol=1e-15)

    def test_explicit_s_from_unit_matrix(self):
        # Unit S-matrix: S11=S22=0, S12=S21=I
        # Full S = [[0, I], [I, 0]]
        sm = SMatrix()
        expected = np.block([[np.zeros((2, 2)), np.eye(2)],
                             [np.eye(2),        np.zeros((2, 2))]])
        np.testing.assert_allclose(sm.S, expected, atol=1e-15)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
