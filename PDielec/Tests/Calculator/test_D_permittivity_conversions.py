"""Tests D: Permittivity / refractive-index conversion functions in Calculator.py.

Covers: calculate_permittivity, calculate_refractive_index,
        calculate_refractive_index_scalar.

Notes
-----
- calculate_permittivity(n) returns n * n  (works with scalars or tensors).
- calculate_refractive_index(eps_3x3) traces the 3x3 tensor and returns a
  complex scalar.
- calculate_refractive_index_scalar(eps_scalar) works on a scalar complex.
- The refractive index is chosen so that Im(n) >= 0 (causal convention).
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import (
    calculate_permittivity,
    calculate_refractive_index,
    calculate_refractive_index_scalar,
)


# ---------------------------------------------------------------------------
# D1: calculate_permittivity
# ---------------------------------------------------------------------------

class TestD1CalculatePermittivity:
    """D1: calculate_permittivity(n) == n * n."""

    def test_n_equals_one_gives_one(self):
        result = calculate_permittivity(1.0 + 0j)
        assert result == pytest.approx(1.0 + 0j)

    def test_n_sqrt2_gives_2(self):
        n = np.sqrt(2.0) + 0j
        result = calculate_permittivity(n)
        assert result.real == pytest.approx(2.0, rel=1e-12)

    def test_complex_n(self):
        n = 1.0 + 1.0j
        result = calculate_permittivity(n)
        expected = (1.0 + 1.0j) ** 2  # = 0 + 2j
        assert result == pytest.approx(expected, abs=1e-12)

    def test_n_1_5(self):
        n = 1.5 + 0j
        result = calculate_permittivity(n)
        assert result.real == pytest.approx(2.25, rel=1e-12)
        assert abs(result.imag) < 1e-12

    def test_zero_n(self):
        result = calculate_permittivity(0.0 + 0j)
        assert result == pytest.approx(0.0 + 0j)


# ---------------------------------------------------------------------------
# D2: calculate_refractive_index  (takes 3x3 tensor, returns scalar)
# ---------------------------------------------------------------------------

class TestD2CalculateRefractiveIndex:
    """D2: calculate_refractive_index traces the 3x3 tensor and returns sqrt."""

    def _make_isotropic_eps(self, n):
        """Build an isotropic 3x3 permittivity tensor for refractive index n."""
        eps = (n * n) * np.eye(3, dtype=complex)
        return eps

    def test_vacuum_eps_gives_n1(self):
        eps = np.eye(3, dtype=complex)
        n = calculate_refractive_index(eps)
        assert abs(n - 1.0) < 1e-12

    def test_n_2_from_eps_4(self):
        eps = 4.0 * np.eye(3, dtype=complex)
        n = calculate_refractive_index(eps)
        assert n.real == pytest.approx(2.0, rel=1e-12)
        assert abs(n.imag) < 1e-12

    def test_round_trip_real(self):
        """calculate_permittivity(calculate_refractive_index(eps)) ≈ eps_trace/3."""
        eps = 2.25 * np.eye(3, dtype=complex)
        n = calculate_refractive_index(eps)
        eps_back = calculate_permittivity(n)
        np.testing.assert_allclose(eps_back.real, 2.25, rtol=1e-12)

    def test_imaginary_part_non_negative(self):
        """Im(n) must be >= 0 (causal sign convention)."""
        # absorbing medium
        eps = (2.0 + 1.0j) * np.eye(3, dtype=complex)
        n = calculate_refractive_index(eps)
        assert n.imag >= -1e-12

    def test_returns_scalar_not_array(self):
        eps = np.eye(3, dtype=complex)
        n = calculate_refractive_index(eps)
        # should be a complex scalar (0-d or Python complex)
        assert np.ndim(n) == 0 or isinstance(n, (complex, np.complexfloating))

    def test_anisotropic_uses_trace(self):
        """For anisotropic eps, result uses trace/3."""
        eps = np.diag([1.0, 4.0, 9.0]).astype(complex)
        n = calculate_refractive_index(eps)
        expected_eps_trace = np.trace(eps) / 3.0  # = (1+4+9)/3 = 14/3
        expected_n = np.sqrt(expected_eps_trace)
        assert abs(n - expected_n) < 1e-10


# ---------------------------------------------------------------------------
# D3: calculate_refractive_index_scalar
# ---------------------------------------------------------------------------

class TestD3CalculateRefractiveIndexScalar:
    """D3: calculate_refractive_index_scalar takes and returns a scalar."""

    def test_one_gives_one(self):
        n = calculate_refractive_index_scalar(1.0 + 0j)
        assert abs(n - 1.0) < 1e-12

    def test_four_gives_two(self):
        n = calculate_refractive_index_scalar(4.0 + 0j)
        assert n.real == pytest.approx(2.0, rel=1e-12)

    def test_complex_eps(self):
        eps = 2.0 + 2.0j
        n = calculate_refractive_index_scalar(eps)
        # n^2 should equal eps
        np.testing.assert_allclose(n * n, eps, atol=1e-12)

    def test_imaginary_part_non_negative(self):
        """Im(n) >= 0 for causal media."""
        eps = 1.0 + 3.0j
        n = calculate_refractive_index_scalar(eps)
        assert n.imag >= -1e-12

    def test_round_trip(self):
        """n = sqrt(eps), eps = n^2, round-trip consistent."""
        eps = complex(3.5, 0.5)
        n = calculate_refractive_index_scalar(eps)
        eps_back = n * n
        np.testing.assert_allclose(eps_back, eps, atol=1e-12)

    def test_real_eps_gives_real_n(self):
        """For real positive eps, n should be real (Im ≈ 0)."""
        n = calculate_refractive_index_scalar(2.25 + 0j)
        assert abs(n.imag) < 1e-12
        assert n.real == pytest.approx(1.5, rel=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
