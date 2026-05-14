"""Tests B1-B4: TabulateScalar, Tabulate1, Tabulate3, Tabulate6.

Plan sections B1, B2, B3, B4.

Each tabulated DielectricFunction wraps a cubic-spline interpolator.
Tests verify:

- Values at tabulation points are reproduced (within spline tolerance).
- Correct scalar/tensor flags.
- epsilon_infinity is (or is not) added correctly.
- Symmetry constraints on the tensor classes.
- Frequency range helpers return the correct bounds.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec import DielectricFunction as DF


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _linspace_vs():
    """Return a simple 10-point frequency axis."""
    return np.linspace(100.0, 1000.0, 10)


def _real_eps(vs):
    """Monotonically increasing real permittivity for testing."""
    return 2.0 + 0.001 * vs


def _complex_eps(vs):
    """Complex permittivity: real part grows, imaginary part is small positive."""
    return (2.0 + 0.001 * vs) + 1j * (0.01 + 0.0001 * vs)


# ---------------------------------------------------------------------------
# B1: TabulateScalar
# ---------------------------------------------------------------------------

class TestB1TabulateScalar:
    """B1: TabulateScalar — scalar spline interpolation."""

    def setup_method(self):
        self.vs = _linspace_vs()
        self.eps = _complex_eps(self.vs)
        self.ts = DF.TabulateScalar(self.vs, self.eps)

    def test_is_scalar_true(self):
        assert self.ts.is_scalar() is True

    def test_is_tensor_false(self):
        assert self.ts.is_tensor() is False

    def test_at_tabulation_points_real_part(self):
        for v, e in zip(self.vs, self.eps):
            result = self.ts.calculate(v)
            assert abs(result.real - e.real) < 1e-6, f"real mismatch at v={v}"

    def test_at_tabulation_points_imag_part(self):
        for v, e in zip(self.vs, self.eps):
            result = self.ts.calculate(v)
            assert abs(result.imag - e.imag) < 1e-6, f"imag mismatch at v={v}"

    def test_interpolated_midpoint_bounded(self):
        v_mid = 0.5 * (self.vs[0] + self.vs[1])
        result = self.ts.calculate(v_mid)
        lo = min(self.eps[0].real, self.eps[1].real)
        hi = max(self.eps[0].real, self.eps[1].real)
        # Allow a small spline overshoot margin
        assert lo - 0.1 < result.real < hi + 0.1

    def test_get_lowest_frequency(self):
        assert abs(self.ts.get_lowest_frequency() - self.vs[0]) < 1e-12

    def test_get_highest_frequency(self):
        assert abs(self.ts.get_highest_frequency() - self.vs[-1]) < 1e-12

    def test_return_type_is_complex(self):
        result = self.ts.calculate(500.0)
        assert isinstance(result, complex)


# ---------------------------------------------------------------------------
# B2: Tabulate1 (isotropic tensor)
# ---------------------------------------------------------------------------

class TestB2Tabulate1:
    """B2: Tabulate1 — isotropic 3x3 tensor with epsilon_infinity."""

    def setup_method(self):
        self.vs = _linspace_vs()
        self.eps = _complex_eps(self.vs)
        self.t1 = DF.Tabulate1(self.vs, self.eps)
        self.eps_inf = np.diag([0.5, 0.5, 0.5])
        self.t1.set_epsilon_infinity(self.eps_inf)

    def test_is_tensor_true(self):
        assert self.t1.is_tensor() is True

    def test_is_scalar_false(self):
        assert self.t1.is_scalar() is False

    def test_shape_is_3x3(self):
        result = self.t1.calculate(500.0)
        assert result.shape == (3, 3)

    def test_diagonal_elements_equal(self):
        result = self.t1.calculate(500.0)
        assert abs(result[0, 0] - result[1, 1]) < 1e-10
        assert abs(result[0, 0] - result[2, 2]) < 1e-10

    def test_off_diagonal_zero(self):
        result = self.t1.calculate(500.0)
        assert abs(result[0, 1]) < 1e-12
        assert abs(result[0, 2]) < 1e-12
        assert abs(result[1, 2]) < 1e-12

    def test_epsilon_infinity_added(self):
        """Without eps_inf a midpoint value, then with: difference == eps_inf."""
        t1_bare = DF.Tabulate1(self.vs, self.eps)
        v = 500.0
        bare = t1_bare.calculate(v)
        with_inf = self.t1.calculate(v)
        diff = with_inf - bare
        np.testing.assert_allclose(diff.real, self.eps_inf, atol=1e-10)

    def test_at_tabulation_point_diagonal_matches_input(self):
        # Without eps_inf, diagonal should match tabulated value
        t1_bare = DF.Tabulate1(self.vs, self.eps)
        for v, e in zip(self.vs, self.eps):
            result = t1_bare.calculate(v)
            assert abs(result[0, 0] - e) < 1e-6


# ---------------------------------------------------------------------------
# B3: Tabulate3 (diagonal anisotropic)
# ---------------------------------------------------------------------------

class TestB3Tabulate3:
    """B3: Tabulate3 — diagonal anisotropic 3x3 tensor."""

    def setup_method(self):
        self.vs = _linspace_vs()
        self.epsxx = _complex_eps(self.vs) * 1.0
        self.epsyy = _complex_eps(self.vs) * 1.5
        self.epszz = _complex_eps(self.vs) * 2.0
        self.t3 = DF.Tabulate3(self.vs, self.epsxx, self.epsyy, self.epszz)

    def test_is_tensor_true(self):
        assert self.t3.is_tensor() is True

    def test_diagonal_components_distinct(self):
        result = self.t3.calculate(500.0)
        # xx != yy != zz because the scaling factors differ
        assert abs(result[0, 0] - result[1, 1]) > 0.1
        assert abs(result[0, 0] - result[2, 2]) > 0.1

    def test_off_diagonal_zero(self):
        result = self.t3.calculate(500.0)
        assert abs(result[0, 1]) < 1e-12
        assert abs(result[0, 2]) < 1e-12
        assert abs(result[1, 2]) < 1e-12

    def test_at_tabulation_point_xx(self):
        for v, e in zip(self.vs, self.epsxx):
            result = self.t3.calculate(v)
            assert abs(result[0, 0] - e) < 1e-6, f"xx mismatch at v={v}"

    def test_at_tabulation_point_yy(self):
        for v, e in zip(self.vs, self.epsyy):
            result = self.t3.calculate(v)
            assert abs(result[1, 1] - e) < 1e-6, f"yy mismatch at v={v}"

    def test_at_tabulation_point_zz(self):
        for v, e in zip(self.vs, self.epszz):
            result = self.t3.calculate(v)
            assert abs(result[2, 2] - e) < 1e-6, f"zz mismatch at v={v}"

    def test_epsilon_infinity_shifts_result(self):
        eps_inf = np.diag([1.0, 2.0, 3.0])
        self.t3.set_epsilon_infinity(eps_inf)
        t3_bare = DF.Tabulate3(self.vs, self.epsxx, self.epsyy, self.epszz)
        v = 500.0
        diff = self.t3.calculate(v) - t3_bare.calculate(v)
        np.testing.assert_allclose(diff.real, eps_inf, atol=1e-10)


# ---------------------------------------------------------------------------
# B4: Tabulate6 (full symmetric tensor)
# ---------------------------------------------------------------------------

class TestB4Tabulate6:
    """B4: Tabulate6 — full symmetric 3x3 tensor."""

    def setup_method(self):
        self.vs = _linspace_vs()
        n = len(self.vs)
        # Create 6 distinct components; off-diagonals are small
        self.epsxx = _complex_eps(self.vs) * 1.0
        self.epsyy = _complex_eps(self.vs) * 1.2
        self.epszz = _complex_eps(self.vs) * 1.4
        self.epsxy = 0.1 * np.ones(n, dtype=complex)
        self.epsxz = 0.2 * np.ones(n, dtype=complex)
        self.epsyz = 0.3 * np.ones(n, dtype=complex)
        self.t6 = DF.Tabulate6(
            self.vs,
            self.epsxx, self.epsyy, self.epszz,
            self.epsxy, self.epsxz, self.epsyz,
        )

    def test_is_tensor_true(self):
        assert self.t6.is_tensor() is True

    def test_tensor_is_symmetric(self):
        result = self.t6.calculate(500.0)
        assert abs(result[0, 1] - result[1, 0]) < 1e-12
        assert abs(result[0, 2] - result[2, 0]) < 1e-12
        assert abs(result[1, 2] - result[2, 1]) < 1e-12

    def test_diagonal_components_at_tabulation_point(self):
        for v, exx, eyy, ezz in zip(self.vs, self.epsxx, self.epsyy, self.epszz):
            result = self.t6.calculate(v)
            assert abs(result[0, 0] - exx) < 1e-6, f"xx at v={v}"
            assert abs(result[1, 1] - eyy) < 1e-6, f"yy at v={v}"
            assert abs(result[2, 2] - ezz) < 1e-6, f"zz at v={v}"

    def test_off_diagonal_at_tabulation_point(self):
        for v, exy, exz, eyz in zip(self.vs, self.epsxy, self.epsxz, self.epsyz):
            result = self.t6.calculate(v)
            assert abs(result[0, 1] - exy) < 1e-6, f"xy at v={v}"
            assert abs(result[0, 2] - exz) < 1e-6, f"xz at v={v}"
            assert abs(result[1, 2] - eyz) < 1e-6, f"yz at v={v}"

    def test_epsilon_infinity_shifts_result(self):
        eps_inf = np.diag([1.0, 2.0, 3.0])
        self.t6.set_epsilon_infinity(eps_inf)
        t6_bare = DF.Tabulate6(
            self.vs,
            self.epsxx, self.epsyy, self.epszz,
            self.epsxy, self.epsxz, self.epsyz,
        )
        v = 500.0
        diff = self.t6.calculate(v) - t6_bare.calculate(v)
        np.testing.assert_allclose(diff.real, eps_inf, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
