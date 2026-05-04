"""Tests C1-C7: DrudeLorentz, FPSQ, Sellmeier and DFT DielectricFunctions.

Plan sections C1, C2, C3, C4, C5, C6, C7.

Tests verify:
  - DrudeLorentz formula at zero frequency (static limit) and far above resonance.
  - FPSQ Lyddane-Sachs-Teller (LST) static limit and tensor behaviour.
  - Sellmeier formula correctness and optical_permittivity() convention.
  - DFT unit-conversion round-trips for frequencies, sigmas, strengths.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec import DielectricFunction as DF
from PDielec.Constants import wavenumber


# ---------------------------------------------------------------------------
# C1-C3: DrudeLorentz
# ---------------------------------------------------------------------------

class TestC1DrudeLorentzIsotropic:
    """C1: DrudeLorentz isotropic — formula and causal behaviour."""

    # Single isotropic mode: v0=400 cm-1, strength=50, sigma=5
    _v0 = 400.0
    _strength = 50.0
    _sigma = 5.0
    _eps_inf_val = 4.0

    def _make(self):
        eps_inf = self._eps_inf_val * np.eye(3)
        dl = DF.DrudeLorentz([self._v0], [self._strength], [self._sigma])
        dl.set_epsilon_infinity(eps_inf)
        return dl

    def test_is_tensor_true(self):
        assert self._make().is_tensor() is True

    def test_is_scalar_false(self):
        assert self._make().is_scalar() is False

    def test_static_limit_analytic(self):
        """At f=0: eps_ii = eps_inf + strength²/v0²."""
        dl = self._make()
        result = dl.calculate(0.0)
        expected = self._eps_inf_val + self._strength**2 / self._v0**2
        assert abs(result[0, 0].real - expected) < 1e-8

    def test_far_above_resonance_approaches_eps_inf(self):
        """Far above resonance the oscillator contributes almost nothing."""
        dl = self._make()
        result = dl.calculate(1e6)
        assert abs(result[0, 0].real - self._eps_inf_val) < 1e-3

    def test_result_is_complex_3x3(self):
        result = self._make().calculate(500.0)
        assert result.shape == (3, 3)
        assert result.dtype == complex

    def test_imaginary_part_positive_near_resonance(self):
        """Lorentzian oscillator must be dissipative: Im(eps) > 0 near v0."""
        dl = self._make()
        result = dl.calculate(self._v0)
        assert result[0, 0].imag > 0

    def test_diagonal_elements_equal_for_isotropic(self):
        """Isotropic constructor copies parameters to all three axes."""
        dl = self._make()
        result = dl.calculate(300.0)
        assert abs(result[0, 0] - result[1, 1]) < 1e-12
        assert abs(result[0, 0] - result[2, 2]) < 1e-12

    def test_off_diagonal_zero(self):
        result = self._make().calculate(300.0)
        assert abs(result[0, 1]) < 1e-12
        assert abs(result[0, 2]) < 1e-12
        assert abs(result[1, 2]) < 1e-12


class TestC2DrudeLorentzRhombic:
    """C2: DrudeLorentz rhombic — each axis has different parameters."""

    def _make(self):
        # Three independent axes: different frequencies and strengths
        omegas    = [[300.0], [500.0], [700.0]]
        strengths = [[40.0],  [60.0],  [80.0]]
        sigmas    = [[5.0],   [5.0],   [5.0]]
        eps_inf = np.diag([2.0, 3.0, 4.0])
        dl = DF.DrudeLorentz(omegas, strengths, sigmas)
        dl.set_epsilon_infinity(eps_inf)
        return dl

    def test_diagonal_values_differ(self):
        result = self._make().calculate(200.0)
        assert abs(result[0, 0] - result[1, 1]) > 0.01
        assert abs(result[0, 0] - result[2, 2]) > 0.01

    def test_off_diagonal_zero(self):
        result = self._make().calculate(200.0)
        assert abs(result[0, 1]) < 1e-12
        assert abs(result[0, 2]) < 1e-12
        assert abs(result[1, 2]) < 1e-12

    def test_static_limit_analytic_per_axis(self):
        """Each axis: eps_ii = eps_inf_ii + strength_i²/omega_i²."""
        dl = self._make()
        result = dl.calculate(0.0)
        expected_xx = 2.0 + 40.0**2 / 300.0**2
        expected_yy = 3.0 + 60.0**2 / 500.0**2
        expected_zz = 4.0 + 80.0**2 / 700.0**2
        assert abs(result[0, 0].real - expected_xx) < 1e-8
        assert abs(result[1, 1].real - expected_yy) < 1e-8
        assert abs(result[2, 2].real - expected_zz) < 1e-8


class TestC3DrudeLorentzRoundTrips:
    """C3: get/set frequencies and sigmas round-trips in cm⁻¹."""

    def _make(self):
        return DF.DrudeLorentz([300.0, 500.0], [40.0, 60.0], [5.0, 8.0])

    def test_set_get_frequencies_roundtrip(self):
        dl = self._make()
        new_freqs = np.array([350.0, 450.0])
        dl.set_frequencies(new_freqs)
        got = dl.get_frequencies()
        np.testing.assert_allclose(got, new_freqs, rtol=1e-10)

    def test_set_get_sigmas_roundtrip(self):
        dl = self._make()
        new_sigmas = [6.0, 9.0]
        dl.set_sigmas(new_sigmas)
        got = dl.get_sigmas()
        np.testing.assert_allclose(got, new_sigmas, rtol=1e-10)


# ---------------------------------------------------------------------------
# C4-C5: FPSQ
# ---------------------------------------------------------------------------

class TestC4FPSQIsotropic:
    """C4: FPSQ isotropic — LST relation and causal behaviour."""

    _eps_inf_val = 4.0
    _omega_to = 200.0   # cm-1
    _gamma_to = 5.0
    _omega_lo = 250.0   # cm-1
    _gamma_lo = 5.0

    def _make(self):
        eps_inf = self._eps_inf_val * np.eye(3)
        fpsq = DF.FPSQ([self._omega_to], [self._gamma_to],
                       [self._omega_lo], [self._gamma_lo])
        fpsq.set_epsilon_infinity(eps_inf)
        return fpsq

    def test_is_tensor_true(self):
        assert self._make().is_tensor() is True

    def test_is_scalar_false(self):
        assert self._make().is_scalar() is False

    def test_lst_static_limit(self):
        """At v→0: eps_ii ≈ eps_inf * (omega_lo/omega_to)²."""
        fpsq = self._make()
        result = fpsq.calculate(1e-3)   # near-zero frequency
        expected = self._eps_inf_val * (self._omega_lo / self._omega_to)**2
        assert abs(result[0, 0].real - expected) < 1e-3

    def test_high_frequency_approaches_eps_inf(self):
        """Far above resonance: eps → eps_inf (contribution ≈ 1)."""
        fpsq = self._make()
        result = fpsq.calculate(1e6)
        assert abs(result[0, 0].real - self._eps_inf_val) < 0.01

    def test_diagonal_elements_equal_for_isotropic(self):
        fpsq = self._make()
        result = fpsq.calculate(100.0)
        assert abs(result[0, 0] - result[1, 1]) < 1e-10
        assert abs(result[0, 0] - result[2, 2]) < 1e-10

    def test_imaginary_part_positive_near_to(self):
        """FPSQ must be causal: Im(eps) >= 0 near omega_TO."""
        fpsq = self._make()
        result = fpsq.calculate(self._omega_to)
        assert result[0, 0].imag >= 0.0

    def test_result_shape(self):
        result = self._make().calculate(100.0)
        assert result.shape == (3, 3)


class TestC5FPSQRhombic:
    """C5: FPSQ rhombic — different parameters per axis."""

    def _make(self):
        # omega_tos and omega_los as lists-of-lists for each axis
        omega_tos = [[200.0], [300.0], [400.0]]
        gamma_tos = [[5.0],   [5.0],   [5.0]]
        omega_los = [[250.0], [380.0], [500.0]]
        gamma_los = [[5.0],   [5.0],   [5.0]]
        eps_inf = np.diag([2.0, 3.0, 4.0])
        fpsq = DF.FPSQ(omega_tos, gamma_tos, omega_los, gamma_los)
        fpsq.set_epsilon_infinity(eps_inf)
        return fpsq

    def test_diagonal_values_differ(self):
        result = self._make().calculate(100.0)
        assert abs(result[0, 0] - result[1, 1]) > 0.01

    def test_off_diagonal_zero(self):
        result = self._make().calculate(100.0)
        assert abs(result[0, 1]) < 1e-12
        assert abs(result[0, 2]) < 1e-12
        assert abs(result[1, 2]) < 1e-12

    def test_lst_per_axis(self):
        """Each axis satisfies its own LST relation."""
        fpsq = self._make()
        result = fpsq.calculate(1e-3)
        eps_infs = [2.0, 3.0, 4.0]
        ratios   = [(250.0/200.0)**2, (380.0/300.0)**2, (500.0/400.0)**2]
        for i, (ei, r) in enumerate(zip(eps_infs, ratios)):
            expected = ei * r
            assert abs(result[i, i].real - expected) < 1e-3, f"axis {i}"


# ---------------------------------------------------------------------------
# C6: Sellmeier
# ---------------------------------------------------------------------------

class TestC6Sellmeier:
    """C6: Sellmeier — formula correctness and optical_permittivity()."""

    # Simple single-term parameters: B=0.5, C=0.1 (micron²)
    _Bs = [0.5]
    _Cs = [0.1]

    def _make(self):
        return DF.Sellmeier(self._Bs, self._Cs)

    def test_is_scalar_true(self):
        assert self._make().is_scalar() is True

    def test_is_tensor_false(self):
        assert self._make().is_tensor() is False

    def test_zero_frequency_limit(self):
        """At v=0: n² = 1 + Σ Bᵢ."""
        sm = self._make()
        expected = 1.0 + sum(self._Bs)
        assert abs(sm.calculate(0.0) - expected) < 1e-12

    @pytest.mark.parametrize("lambda_um", [2.0, 5.0, 10.0])
    def test_formula_at_known_wavelength(self, lambda_um):
        """n² = 1 + B * λ² / (λ² - C) for our single-term Sellmeier."""
        sm = self._make()
        v = 1.0e4 / lambda_um          # cm-1
        lam2 = lambda_um**2
        expected = 1.0 + self._Bs[0] * lam2 / (lam2 - self._Cs[0])
        assert abs(sm.calculate(v) - expected) < 1e-10

    def test_result_is_real(self):
        sm = self._make()
        result = sm.calculate(1000.0)
        assert abs(np.imag(result)) < 1e-15

    def test_optical_permittivity_at_589nm(self):
        """optical_permittivity() must agree with calculate() at 589 nm."""
        sm = self._make()
        v_589nm = 1.0e4 / 0.589
        expected = float(np.real(sm.calculate(v_589nm)))
        assert abs(sm.optical_permittivity() - expected) < 1e-10

    def test_multi_term_zero_frequency(self):
        """Multi-term: n²(0) = 1 + ΣBᵢ."""
        sm = DF.Sellmeier([0.6, 0.3, 0.1], [0.004, 0.01, 100.0])
        expected = 1.0 + 0.6 + 0.3 + 0.1
        assert abs(sm.calculate(0.0) - expected) < 1e-12


# ---------------------------------------------------------------------------
# C7: DFT
# ---------------------------------------------------------------------------

class TestC7DFT:
    """C7: DFT — unit-conversion round-trips and basic calculation."""

    def _make_strengths(self, n_modes):
        """Return isotropic oscillator strengths as n_modes 3x3 arrays."""
        return [np.eye(3) * float(i + 1) * wavenumber for i in range(n_modes)]

    def _make(self, n_modes=3):
        mode_list = list(range(n_modes))
        freqs_au  = np.array([100.0, 200.0, 300.0][:n_modes]) * wavenumber
        sigmas_au = np.array([5.0, 5.0, 5.0][:n_modes]) * wavenumber
        strengths = self._make_strengths(n_modes)
        volume_au = (10.0)**3 * (1e-10)**3   # 10 Å cube in m³ — will be re-used as a.u. placeholder
        dft = DF.DFT(
            mode_list=mode_list,
            mode_frequencies=freqs_au,
            mode_sigmas=sigmas_au,
            mode_oscillator_strengths=strengths,
            crystal_volume=volume_au,
            drude=False,
            drude_plasma=0.0,
            drude_sigma=0.0,
        )
        dft.set_epsilon_infinity(np.diag([2.0, 2.0, 2.0]))
        return dft

    def test_calculate_returns_3x3_complex(self):
        result = self._make().calculate(50.0)
        assert result.shape == (3, 3)
        assert np.issubdtype(result.dtype, complex)

    def test_set_get_frequencies_roundtrip(self):
        dft = self._make()
        new_freqs = np.array([110.0, 210.0, 310.0])
        dft.set_frequencies(new_freqs)
        got = dft.get_frequencies()
        np.testing.assert_allclose(got, new_freqs, rtol=1e-10)

    def test_set_get_sigmas_roundtrip(self):
        dft = self._make()
        new_sigmas = np.array([6.0, 7.0, 8.0])
        dft.set_sigmas(new_sigmas)
        got = dft.get_sigmas()
        np.testing.assert_allclose(got, new_sigmas, rtol=1e-10)

    def test_zero_strengths_gives_eps_inf(self):
        """With zero oscillator strengths the result should equal epsilon_infinity."""
        dft = DF.DFT(
            mode_list=[0],
            mode_frequencies=np.array([200.0 * wavenumber]),
            mode_sigmas=np.array([5.0 * wavenumber]),
            mode_oscillator_strengths=[np.zeros((3, 3))],
            crystal_volume=1.0,
            drude=False,
            drude_plasma=0.0,
            drude_sigma=0.0,
        )
        eps_inf = np.diag([3.0, 3.0, 3.0])
        dft.set_epsilon_infinity(eps_inf)
        result = dft.calculate(50.0)
        np.testing.assert_allclose(result.real, eps_inf, atol=1e-10)
        np.testing.assert_allclose(result.imag, np.zeros((3, 3)), atol=1e-10)

    def test_is_tensor_true(self):
        assert self._make().is_tensor() is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
