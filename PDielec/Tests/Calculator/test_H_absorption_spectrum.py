"""Tests H: Absorption spectrum utilities in Calculator.py.

Covers: absorption_from_mode_intensities, calculate_size_factor.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import absorption_from_mode_intensities, calculate_size_factor


# ---------------------------------------------------------------------------
# H1: absorption_from_mode_intensities
# ---------------------------------------------------------------------------

class TestH1AbsorptionFromModeIntensities:
    """H1: Lorentzian absorption spectrum from mode frequencies and intensities."""

    def _setup_single_mode(self, v0=1000.0, sigma=5.0, intensity=1.0):
        """Return args for a single mode at v0 with width sigma."""
        modes = [0]
        frequencies = [v0]
        sigmas = [sigma]
        intensities = [intensity]
        return modes, frequencies, sigmas, intensities

    def test_zero_intensity_gives_zero(self):
        modes, freqs, sigmas, _ = self._setup_single_mode()
        intensities = [0.0]
        result = absorption_from_mode_intensities(1000.0, modes, freqs, sigmas, intensities)
        assert result == pytest.approx(0.0, abs=1e-20)

    def test_peak_at_mode_frequency(self):
        """Absorption at the mode frequency is the maximum."""
        v0 = 500.0
        modes, freqs, sigmas, intensities = self._setup_single_mode(v0=v0)
        at_peak = absorption_from_mode_intensities(v0, modes, freqs, sigmas, intensities)
        at_off = absorption_from_mode_intensities(v0 + 20.0, modes, freqs, sigmas, intensities)
        assert at_peak > at_off

    def test_non_negative_across_spectrum(self):
        """Absorption must be non-negative everywhere."""
        modes, freqs, sigmas, intensities = self._setup_single_mode()
        for v in np.linspace(900.0, 1100.0, 50):
            result = absorption_from_mode_intensities(v, modes, freqs, sigmas, intensities)
            assert result >= -1e-15

    def test_returns_scalar(self):
        modes, freqs, sigmas, intensities = self._setup_single_mode()
        result = absorption_from_mode_intensities(1000.0, modes, freqs, sigmas, intensities)
        assert np.ndim(result) == 0

    def test_empty_mode_list_gives_zero(self):
        result = absorption_from_mode_intensities(1000.0, [], [1000.0], [5.0], [1.0])
        assert result == pytest.approx(0.0, abs=1e-20)

    def test_lorentzian_halfwidth_at_half_max(self):
        """At f = v0 ± sigma/2, absorption should be ~50% of peak (Lorentzian HWHM)."""
        v0 = 1000.0
        sigma = 10.0
        modes, freqs, sigmas, intensities = self._setup_single_mode(v0=v0, sigma=sigma)
        at_peak = absorption_from_mode_intensities(v0, modes, freqs, sigmas, intensities)
        at_hwhm = absorption_from_mode_intensities(v0 + sigma / 2.0, modes, freqs, sigmas, intensities)
        # For a Lorentzian, at half-width = sigma/2, intensity should be ~0.5 of peak
        assert abs(at_hwhm / at_peak - 0.5) < 0.01

    def test_two_mode_superposition(self):
        """Two separated modes: total = sum of individual Lorentzians."""
        modes = [0, 1]
        frequencies = [900.0, 1100.0]
        sigmas = [5.0, 5.0]
        intensities = [1.0, 1.0]
        # Far from both peaks, contributions should be small
        result = absorption_from_mode_intensities(1000.0, modes, frequencies, sigmas, intensities)
        single1 = absorption_from_mode_intensities(1000.0, [0], frequencies, sigmas, intensities)
        single2 = absorption_from_mode_intensities(1000.0, [1], frequencies, sigmas, intensities)
        assert result == pytest.approx(single1 + single2, rel=1e-10)


# ---------------------------------------------------------------------------
# H2: calculate_size_factor
# ---------------------------------------------------------------------------

class TestH2CalculateSizeFactor:
    """H2: Size correction factor (Sihvola Eqs 10.38/10.39)."""

    def test_zero_x_gives_one(self):
        """At x=0 (no size correction), factor should be exactly 1."""
        result = calculate_size_factor(0.0)
        assert result == pytest.approx(1.0, abs=1e-15)

    def test_very_small_x_near_one(self):
        """For very small x, result should be very close to 1."""
        result = calculate_size_factor(1e-10)
        assert abs(result - 1.0) < 1e-6

    def test_returns_complex_or_real(self):
        """Result should be a number (complex or real)."""
        result = calculate_size_factor(0.1)
        assert np.isscalar(result) or np.ndim(result) == 0

    def test_below_threshold_returns_one(self):
        """x < 1e-12: returns exactly 1."""
        result = calculate_size_factor(1e-13)
        assert result == 1.0

    def test_positive_x_not_absurdly_large(self):
        """For reasonable x (0 < x < 10), result should not diverge."""
        for x in [0.01, 0.1, 0.5, 1.0, 2.0]:
            result = calculate_size_factor(x)
            assert np.isfinite(result) or np.isfinite(abs(result))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
