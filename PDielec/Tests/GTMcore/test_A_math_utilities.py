"""Tests A1-A4: Mathematical utility functions in GTMcore.

Plan sections A1 (vacuum_eps), A2 (exact_inv_2x2), A3 (exact_inv_3x3),
A4 (exact_inv_4x4).

Pure algebra — no Layer or System objects required.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import vacuum_eps, exact_inv_2x2, exact_inv_3x3, exact_inv_4x4


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

rng = np.random.default_rng(42)


def _random_complex_invertible(n):
    """Return a random well-conditioned complex n×n matrix."""
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    # Shift diagonal to improve conditioning
    A += n * np.eye(n, dtype=complex)
    return A.astype(np.cdouble)


# ---------------------------------------------------------------------------
# A1: vacuum_eps
# ---------------------------------------------------------------------------

class TestA1VacuumEps:
    """A1: vacuum_eps returns the vacuum permittivity (1+0j) at any frequency."""

    def test_scalar_frequency_returns_complex(self):
        result = vacuum_eps(1e12)
        assert isinstance(result, complex), f"Expected complex, got {type(result)}"

    def test_scalar_value_is_one(self):
        for f in [0.0, 100.0, 1e13, 3e14]:
            val = vacuum_eps(f)
            assert abs(val - (1.0 + 0j)) < 1e-15, f"vacuum_eps({f}) = {val} != 1+0j"

    def test_array_frequency_returns_array(self):
        freqs = np.array([1e12, 2e12, 3e12])
        result = vacuum_eps(freqs)
        assert hasattr(result, "__len__"), "Expected array result for array input"
        assert len(result) == 3

    def test_array_all_ones(self):
        freqs = np.linspace(1e11, 1e15, 10)
        result = vacuum_eps(freqs)
        np.testing.assert_allclose(result, np.ones(10), atol=1e-15)


# ---------------------------------------------------------------------------
# A2: exact_inv_2x2
# ---------------------------------------------------------------------------

class TestA2ExactInv2x2:
    """A2: exact_inv_2x2 computes the correct inverse of a 2×2 matrix."""

    def test_identity_inverts_to_identity(self):
        I = np.eye(2, dtype=np.cdouble)
        result = exact_inv_2x2(I)
        np.testing.assert_allclose(result, I, atol=1e-14)

    def test_M_times_inv_M_is_identity(self):
        M = _random_complex_invertible(2)
        Minv = exact_inv_2x2(M)
        product = M @ Minv
        np.testing.assert_allclose(product, np.eye(2), atol=1e-12)

    def test_known_inverse(self):
        # [[1,2],[3,4]]^-1 = [[-2,1],[1.5,-0.5]]
        M = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.cdouble)
        expected = np.array([[-2.0, 1.0], [1.5, -0.5]], dtype=np.cdouble)
        result = exact_inv_2x2(M)
        np.testing.assert_allclose(result, expected, atol=1e-14)

    def test_result_is_complex(self):
        M = _random_complex_invertible(2)
        result = exact_inv_2x2(M)
        assert result.dtype in (np.cdouble, np.complex128, np.complex64)


# ---------------------------------------------------------------------------
# A3: exact_inv_3x3
# ---------------------------------------------------------------------------

class TestA3ExactInv3x3:
    """A3: exact_inv_3x3 computes the correct inverse of a 3×3 matrix."""

    def test_identity_inverts_to_identity(self):
        I = np.eye(3, dtype=np.cdouble)
        result = exact_inv_3x3(I)
        np.testing.assert_allclose(result, I, atol=1e-14)

    def test_M_times_inv_M_is_identity(self):
        M = _random_complex_invertible(3)
        Minv = exact_inv_3x3(M)
        product = M @ Minv
        np.testing.assert_allclose(product, np.eye(3), atol=1e-12)

    def test_diagonal_matrix_inverse(self):
        d = np.array([2.0 + 0j, 3.0 + 0j, 4.0 + 0j])
        M = np.diag(d)
        result = exact_inv_3x3(M)
        expected = np.diag(1.0 / d)
        np.testing.assert_allclose(result, expected, atol=1e-14)

    def test_agrees_with_numpy_inv(self):
        M = _random_complex_invertible(3)
        result = exact_inv_3x3(M)
        expected = np.linalg.inv(M)
        np.testing.assert_allclose(result, expected, atol=1e-10)


# ---------------------------------------------------------------------------
# A4: exact_inv_4x4
# ---------------------------------------------------------------------------

class TestA4ExactInv4x4:
    """A4: exact_inv_4x4 computes the correct inverse of a 4×4 matrix."""

    def test_identity_inverts_to_identity(self):
        I = np.eye(4, dtype=np.cdouble)
        result = exact_inv_4x4(I)
        np.testing.assert_allclose(result, I, atol=1e-14)

    def test_M_times_inv_M_is_identity(self):
        M = _random_complex_invertible(4)
        Minv = exact_inv_4x4(M)
        product = M @ Minv
        np.testing.assert_allclose(product, np.eye(4), atol=1e-11)

    def test_diagonal_matrix_inverse(self):
        d = np.array([1.0 + 0j, 2.0 + 0j, 3.0 + 0j, 5.0 + 0j])
        M = np.diag(d)
        result = exact_inv_4x4(M)
        expected = np.diag(1.0 / d)
        np.testing.assert_allclose(result, expected, atol=1e-14)

    def test_agrees_with_numpy_inv(self):
        M = _random_complex_invertible(4)
        result = exact_inv_4x4(M)
        expected = np.linalg.inv(M)
        np.testing.assert_allclose(result, expected, atol=1e-10)

    def test_complex_matrix_matches_numpy(self):
        M = (rng.standard_normal((4, 4)) + 5.0 * np.eye(4)
             + 1j * rng.standard_normal((4, 4))).astype(np.cdouble)
        result = exact_inv_4x4(M)
        expected = np.linalg.inv(M)
        np.testing.assert_allclose(result, expected, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
