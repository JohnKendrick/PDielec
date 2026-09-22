"""Tests A1-A3: ConstantScalar, ConstantTensor and DielectricFunction base helpers.

Plan sections A1, A2, A3.

Pure algebra — no file I/O or DFT data required.
Tests verify:

- ConstantScalar returns the same value at all frequencies and has the
  correct scalar/tensor flags.
- ConstantTensor returns value*I + eps_inf and has the correct flags.
- Base class helpers (set_volume, set_epsilon_infinity,
  dielectric_contributions_from_modes/drude) behave correctly.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec import DielectricFunction as DF
from PDielec.Constants import angstrom


# ---------------------------------------------------------------------------
# A1: ConstantScalar
# ---------------------------------------------------------------------------

class TestA1ConstantScalar:
    """A1: ConstantScalar returns a constant value regardless of frequency."""

    def test_calculate_real(self):
        cs = DF.ConstantScalar(2.5)
        for v in [0.0, 100.0, 500.0, 1000.0, 5000.0]:
            assert cs.calculate(v) == 2.5

    def test_calculate_complex(self):
        val = complex(2.0, 0.3)
        cs = DF.ConstantScalar(val)
        for v in [0.0, 200.0, 1500.0]:
            assert cs.calculate(v) == val

    def test_is_scalar_true(self):
        cs = DF.ConstantScalar(1.5)
        assert cs.is_scalar() is True

    def test_is_tensor_false(self):
        cs = DF.ConstantScalar(1.5)
        assert cs.is_tensor() is False

    def test_optical_permittivity_real_value(self):
        cs = DF.ConstantScalar(3.0)
        op = cs.optical_permittivity()
        assert isinstance(op, float)
        assert abs(op - 3.0) < 1e-12

    def test_optical_permittivity_complex_takes_real_part(self):
        cs = DF.ConstantScalar(complex(4.0, 1.0))
        op = cs.optical_permittivity()
        assert isinstance(op, float)
        assert abs(op - 4.0) < 1e-12

    def test_function_callable_matches_calculate(self):
        cs = DF.ConstantScalar(2.25)
        fn = cs.function()
        for v in [0.0, 100.0, 500.0]:
            assert fn(v) == cs.calculate(v)


# ---------------------------------------------------------------------------
# A2: ConstantTensor
# ---------------------------------------------------------------------------

class TestA2ConstantTensor:
    """A2: ConstantTensor returns value*I + eps_inf at every frequency."""

    _value = 3.0
    _eps_inf = np.diag([1.0, 2.0, 3.0])

    def _make(self):
        ct = DF.ConstantTensor(self._value)
        ct.set_epsilon_infinity(self._eps_inf)
        return ct

    def test_calculate_shape(self):
        ct = self._make()
        result = ct.calculate(0.0)
        assert result.shape == (3, 3)

    def test_calculate_diagonal_values(self):
        ct = self._make()
        result = ct.calculate(500.0)
        expected_diag = np.array([1.0 + 3.0, 2.0 + 3.0, 3.0 + 3.0])
        np.testing.assert_allclose(np.diag(result).real, expected_diag, atol=1e-12)

    def test_calculate_off_diagonal_zero(self):
        ct = DF.ConstantTensor(2.0)
        result = ct.calculate(0.0)
        off = result - np.diag(np.diag(result))
        np.testing.assert_allclose(np.abs(off), 0.0, atol=1e-12)

    def test_calculate_independent_of_frequency(self):
        ct = self._make()
        r1 = ct.calculate(0.0)
        r2 = ct.calculate(1000.0)
        np.testing.assert_array_equal(r1, r2)

    def test_is_scalar_false(self):
        ct = DF.ConstantTensor(1.0)
        assert ct.is_scalar() is False

    def test_is_tensor_true(self):
        ct = DF.ConstantTensor(1.0)
        assert ct.is_tensor() is True

    def test_optical_permittivity_returns_eps_inf(self):
        ct = self._make()
        op = ct.optical_permittivity()
        np.testing.assert_allclose(op, self._eps_inf, atol=1e-12)

    def test_optical_permittivity_zero_eps_inf(self):
        ct = DF.ConstantTensor(5.0)
        # eps_inf defaults to zeros
        op = ct.optical_permittivity()
        np.testing.assert_allclose(op, np.zeros((3, 3)), atol=1e-12)

    def test_set_epsilon_infinity_updates_calculate(self):
        ct = DF.ConstantTensor(1.0)
        new_eps = 2.0 * np.eye(3)
        ct.set_epsilon_infinity(new_eps)
        result = ct.calculate(0.0)
        np.testing.assert_allclose(np.diag(result).real, [3.0, 3.0, 3.0], atol=1e-12)


# ---------------------------------------------------------------------------
# A3: DielectricFunction base helpers
# ---------------------------------------------------------------------------

class TestA3BaseHelpers:
    """A3: set_volume, set_epsilon_infinity, contribution helpers."""

    def test_set_volume_stores_angstrom_value(self):
        cs = DF.ConstantScalar(1.0)
        cs.set_volume(100.0)
        assert abs(cs.volume_angs - 100.0) < 1e-12

    def test_set_volume_converts_to_au(self):
        cs = DF.ConstantScalar(1.0)
        volume_angs = 50.0
        cs.set_volume(volume_angs)
        expected_au = volume_angs * angstrom**3
        assert abs(cs.volume_au - expected_au) / expected_au < 1e-12

    def test_set_epsilon_infinity_stores_array(self):
        cs = DF.ConstantScalar(1.0)
        eps = np.diag([2.0, 3.0, 4.0])
        cs.set_epsilon_infinity(eps)
        np.testing.assert_array_equal(cs.epsilon_infinity, eps)

    def test_set_epsilon_infinity_from_list(self):
        cs = DF.ConstantScalar(1.0)
        eps_list = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        cs.set_epsilon_infinity(eps_list)
        assert cs.epsilon_infinity.shape == (3, 3)

    def test_contributions_from_modes_zero_modes(self):
        """With an empty mode list, contribution should be zero."""
        base = DF.ConstantScalar(1.0)
        base.set_volume(100.0)
        result = base.dielectric_contributions_from_modes(
            100.0, [], [], [], [], base.volume_au
        )
        np.testing.assert_allclose(np.abs(result), 0.0, atol=1e-30)

    def test_contributions_from_drude_far_from_resonance(self):
        """Far above the Drude frequency, the contribution should be small."""
        base = DF.ConstantScalar(1.0)
        base.set_volume(1000.0)
        from PDielec.Constants import wavenumber
        drude_freq_au = 100.0 * wavenumber
        sigma_au = 5.0 * wavenumber
        vol_au = base.volume_angs * angstrom**3 if hasattr(base, "volume_angs") else 1.0
        # evaluate at a very high frequency
        result = base.dielectric_contributions_from_drude(
            1e8 * wavenumber, drude_freq_au, sigma_au, 1.0
        )
        assert np.max(np.abs(result)) < 1e-10


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
