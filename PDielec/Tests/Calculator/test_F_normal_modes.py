"""Tests F: Normal modes and oscillator/Raman/IR intensities in Calculator.py.

Covers: oscillator_strengths, normal_modes, infrared_intensities,
        raman_intensities, calculate_normal_modes_and_frequencies.

Notes
-----
- calculate_normal_modes_and_frequencies returns (mass_weighted_normal_modes, frequencies).
- oscillator_strengths returns array of shape (nmodes, 3, 3).
- infrared_intensities returns 1-D array of length nmodes.
- raman_intensities returns array of shape (nmodes, 3).
- normal_modes transforms mass-weighted to un-mass-weighted coordinates.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import math
import numpy as np
import pytest
from PDielec.Calculator import (
    oscillator_strengths,
    normal_modes,
    infrared_intensities,
    raman_intensities,
    calculate_normal_modes_and_frequencies,
)
from PDielec.Constants import wavenumber


# ---------------------------------------------------------------------------
# F1: oscillator_strengths
# ---------------------------------------------------------------------------

class TestF1OscillatorStrengths:
    """F1: oscillator_strengths(normal_modes, born_charges) -> (nmodes, 3, 3)."""

    def _make_simple_system(self):
        """One atom, one mode, diagonal Born charge Z, mode displacement along x."""
        # 1 atom, 1 mode: normal mode is [1, 0, 0] (x displacement)
        # Born charge is identity (1 e along each direction)
        nmodes = 3  # 3 degrees of freedom for 1 atom
        natoms = 1
        mass_wt_modes = np.zeros((nmodes, natoms, 3))
        mass_wt_modes[0, 0, :] = [1.0, 0.0, 0.0]  # mode 0: x displacement
        mass_wt_modes[1, 0, :] = [0.0, 1.0, 0.0]  # mode 1: y displacement
        mass_wt_modes[2, 0, :] = [0.0, 0.0, 1.0]  # mode 2: z displacement
        born = [np.eye(3)]  # Born charge for atom 0: identity
        return mass_wt_modes, born

    def test_shape(self):
        modes, born = self._make_simple_system()
        strengths = oscillator_strengths(modes, born)
        assert strengths.shape == (3, 3, 3)

    def test_zero_born_charges_zero_strengths(self):
        modes, _ = self._make_simple_system()
        born_zero = [np.zeros((3, 3))]
        strengths = oscillator_strengths(modes, born_zero)
        np.testing.assert_allclose(np.abs(strengths), 0.0, atol=1e-15)

    def test_identity_born_diagonal_mode_x(self):
        """Mode 0 (x displacement) with identity Born charge → outer([1,0,0],[1,0,0])."""
        modes, born = self._make_simple_system()
        strengths = oscillator_strengths(modes, born)
        expected = np.outer([1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        np.testing.assert_allclose(strengths[0], expected, atol=1e-12)

    def test_symmetry(self):
        """Oscillator strength tensor for each mode should be symmetric."""
        modes, born = self._make_simple_system()
        strengths = oscillator_strengths(modes, born)
        for s in strengths:
            np.testing.assert_allclose(s, s.T, atol=1e-12)


# ---------------------------------------------------------------------------
# F2: normal_modes
# ---------------------------------------------------------------------------

class TestF2NormalModes:
    """F2: normal_modes divides mass-weighted modes by sqrt(mass)."""

    def test_shape(self):
        natoms = 2
        nmodes = 6
        masses = np.array([1.0, 4.0])  # in atomic units
        mwmodes = np.random.default_rng(0).random((nmodes, natoms, 3))
        result = normal_modes(masses, mwmodes)
        assert result.shape == (nmodes, natoms, 3)

    def test_single_atom_unit_mass(self):
        """For mass=1.0, normal_modes == mass_weighted_normal_modes."""
        masses = np.array([1.0])
        mwmodes = np.ones((3, 1, 3))
        result = normal_modes(masses, mwmodes)
        np.testing.assert_allclose(result, mwmodes, atol=1e-12)

    def test_mass_scaling(self):
        """Atom with mass m: result = mw_mode / sqrt(m)."""
        masses = np.array([4.0])
        mwmodes = np.ones((1, 1, 3))
        result = normal_modes(masses, mwmodes)
        np.testing.assert_allclose(result, mwmodes / 2.0, atol=1e-12)


# ---------------------------------------------------------------------------
# F3: infrared_intensities
# ---------------------------------------------------------------------------

class TestF3InfraredIntensities:
    """F3: infrared_intensities returns the trace-based IR intensities."""

    def test_shape(self):
        strengths = np.zeros((5, 3, 3))
        result = infrared_intensities(strengths)
        assert result.shape == (5,)

    def test_zero_strengths_zero_intensities(self):
        strengths = np.zeros((4, 3, 3))
        result = infrared_intensities(strengths)
        np.testing.assert_allclose(result, np.zeros(4), atol=1e-30)

    def test_positive_diagonal_strengths(self):
        """Each mode has a positive diagonal → positive intensity."""
        strengths = np.zeros((3, 3, 3))
        for i in range(3):
            strengths[i, i, i] = 1.0
        result = infrared_intensities(strengths)
        assert np.all(result >= 0.0)


# ---------------------------------------------------------------------------
# F4: raman_intensities
# ---------------------------------------------------------------------------

class TestF4RamanIntensities:
    """F4: raman_intensities returns (nmodes, 3) array of Raman activities."""

    def test_shape(self):
        tensors = [np.eye(3) for _ in range(4)]
        result = raman_intensities(tensors, volume=100.0)
        assert result.shape == (4, 3)

    def test_zero_tensors_zero_intensities(self):
        tensors = [np.zeros((3, 3)) for _ in range(3)]
        result = raman_intensities(tensors, volume=100.0)
        np.testing.assert_allclose(result, np.zeros((3, 3)), atol=1e-15)

    def test_all_non_negative(self):
        rng = np.random.default_rng(7)
        tensors = [rng.random((3, 3)) for _ in range(5)]
        result = raman_intensities(tensors, volume=50.0)
        assert np.all(result >= 0.0)

    def test_total_equals_parallel_plus_perp(self):
        """Total = parallel + perpendicular for symmetric tensors (kappa=0)."""
        R = np.diag([1.0, 2.0, 3.0])  # symmetric, traceless antisymmetric part = 0
        result = raman_intensities([R], volume=1.0)
        # total == parallel + perpendicular only when kappa^2 = 0 (symmetric tensor)
        assert result[0, 0] == pytest.approx(result[0, 1] + result[0, 2], rel=1e-10)

    def test_isotropic_tensor_no_depolarisation(self):
        """Isotropic tensor: parallel >> perpendicular."""
        R = np.eye(3) * 5.0
        result = raman_intensities([R], volume=1.0)
        assert result[0, 1] > result[0, 2]  # parallel > perp for isotropic


# ---------------------------------------------------------------------------
# F5: calculate_normal_modes_and_frequencies
# ---------------------------------------------------------------------------

class TestF5NormalModesAndFrequencies:
    """F5: calculate_normal_modes_and_frequencies diagonalises a mass-weighted Hessian."""

    def _make_1d_harmonic(self, k_au, mass_au):
        """One-atom 1-D Hessian: 1x1 = k / mass (mass-weighted).
        Frequency = sqrt(k/mass) / wavenumber."""
        # 3 modes, but only one non-zero (the x direction with spring k)
        n = 3
        H = np.zeros((n, n))
        H[0, 0] = k_au / mass_au  # mass-weighted spring in x
        return H

    def test_returns_two_items(self):
        H = np.diag([1.0, 1.0, 1.0])
        result = calculate_normal_modes_and_frequencies(H)
        assert len(result) == 2

    def test_modes_come_first(self):
        """Return order is (mass_weighted_normal_modes, frequencies)."""
        H = np.diag([1.0, 2.0, 3.0])
        modes, freqs = calculate_normal_modes_and_frequencies(H)
        # modes should be (3, 1, 3) for 3 modes, 1 atom
        assert hasattr(modes, '__len__')
        assert len(freqs) == 3

    def test_positive_definite_gives_real_positive_freqs(self):
        """For a positive-definite Hessian, all frequencies must be positive."""
        H = np.diag([1.0, 4.0, 9.0])  # eigenvalues 1, 4, 9
        modes, freqs = calculate_normal_modes_and_frequencies(H)
        assert np.all(freqs > 0.0)

    def test_known_harmonic_frequency(self):
        """H = k/m (1D): expected freq = sqrt(k/m) / wavenumber."""
        k_force_au = 1.0e-4
        mass_au = 1.0
        k_mw = k_force_au / mass_au
        H = np.array([[k_mw, 0, 0],
                      [0, 0,   0],
                      [0, 0,   0]])
        modes, freqs = calculate_normal_modes_and_frequencies(H)
        expected_cm1 = math.sqrt(k_force_au / mass_au) / wavenumber
        # The non-zero mode should match
        non_zero_freqs = sorted(freqs)
        np.testing.assert_allclose(non_zero_freqs[-1], expected_cm1, rtol=1e-8)

    def test_frequencies_length_matches_hessian(self):
        n = 6  # 2 atoms × 3 DOF
        H = np.eye(n) * 0.01
        modes, freqs = calculate_normal_modes_and_frequencies(H)
        assert len(freqs) == n

    def test_acoustic_modes_near_zero(self):
        """A Hessian with zero eigenvalues gives near-zero (acoustic) frequencies."""
        H = np.zeros((3, 3))  # all zero → zero eigenvalues
        modes, freqs = calculate_normal_modes_and_frequencies(H)
        np.testing.assert_allclose(np.abs(freqs), 0.0, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
