"""Tests C1-C4: Unit-conversion self-consistency.

Plan sections C1-C4.

These tests verify that conversion factors are mutually consistent without
relying on knowing the exact CODATA values.  They are complementary to the
absolute-value checks in test_A and test_B.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
import PDielec.Constants as C


# ---------------------------------------------------------------------------
# C1: Energy unit round-trips
# ---------------------------------------------------------------------------

class TestC1EnergyRoundTrips:
    """C1: hartree and hartree2ev round-trip consistency."""

    def test_hartree_times_avogadro_positive(self):
        # 1 hartree/mol is a sensible molar energy (should be ~2600 kJ/mol)
        hartree_J = 1.0 / C.joule
        molar_energy = hartree_J * C.avogadro_si
        assert molar_energy > 0
        assert molar_energy > 1e5  # definitely > 100 kJ/mol

    def test_hartree2ev_times_charge_matches_hartree(self):
        # hartree (in J) = hartree2ev * elementary_charge_si
        hartree_J = 1.0 / C.joule
        from_ev = C.hartree2ev * C.elementary_charge_si
        np.testing.assert_allclose(hartree_J, from_ev, rtol=1e-9)

    def test_joule_and_hartree_reciprocal(self):
        # hartree = 1/joule  →  hartree * joule = 1
        np.testing.assert_allclose(C.hartree * C.joule, 1.0, rtol=1e-9)


# ---------------------------------------------------------------------------
# C2: Length unit consistency
# ---------------------------------------------------------------------------

class TestC2LengthUnits:
    """C2: angstrom in atomic units and angs2bohr round-trip."""

    def test_angstrom_in_si_is_1e_minus_10_metres(self):
        # angstrom (in atomic units) / metre (in atomic units) = 1e-10
        ratio = C.angstrom / C.metre
        np.testing.assert_allclose(ratio, 1.0e-10, rtol=1e-9)

    def test_angs2bohr_times_bohr_si_equals_1e_minus10(self):
        # 1 Å = angs2bohr bohr radii; 1 bohr_si metres → angs2bohr * bohr_si = 1e-10
        np.testing.assert_allclose(C.angs2bohr * C.bohr_si, 1.0e-10, rtol=1e-9)

    def test_angstrom_consistent_with_angs2bohr(self):
        # angstrom (AU length) = angs2bohr^{-1} in Bohr units, but since
        # the Bohr radius = 1 in atomic length units, angstrom in AU = 1/angs2bohr
        # Actually: angstrom = 1e-10 * metre (AU length)
        # angs2bohr = 1e-10 / bohr_si (number of bohr per Angstrom)
        # Check that C.angstrom = C.angs2bohr (since 1 AU length = 1 bohr)
        # angstrom in AU = angs2bohr (number of bohr in 1 Angstrom)
        np.testing.assert_allclose(C.angstrom, C.angs2bohr, rtol=1e-8)


# ---------------------------------------------------------------------------
# C3: Pressure unit
# ---------------------------------------------------------------------------

class TestC3PressureUnit:
    """C3: au2GPa converts pressure from atomic units to GPa."""

    def test_au2GPa_positive(self):
        assert C.au2GPa > 0

    def test_au2GPa_approx_29421(self):
        # 1 Hartree/Bohr³ ≈ 29421 GPa
        np.testing.assert_allclose(C.au2GPa, 29421.0, rtol=5e-3)

    def test_au2GPa_from_hartree_and_bohr(self):
        # au2GPa = (1 hartree / bohr_si^3) in GPa
        hartree_J = 1.0 / C.joule
        expected_Pa = hartree_J / C.bohr_si**3
        expected_GPa = expected_Pa * 1.0e-9
        np.testing.assert_allclose(C.au2GPa, expected_GPa, rtol=1e-8)


# ---------------------------------------------------------------------------
# C4: Speed of light and wavenumber
# ---------------------------------------------------------------------------

class TestC4SpeedOfLight:
    """C4: Consistency of wavenumber with speed_light_si."""

    def test_speed_of_light_in_cm_per_s(self):
        # c ≈ 3e10 cm/s
        c_cms = C.speed_light_si * 100.0
        np.testing.assert_allclose(c_cms, 2.99792458e10, rtol=1e-9)

    def test_thz2cm1_from_speed_of_light(self):
        # thz2cm1 = 1e10 / speed_light_si  (1 THz = ? cm⁻¹)
        # 1 THz = 1e12 Hz; ν = f/c where c in cm/s = 1e12/(3e10) ≈ 33.36 cm⁻¹
        expected = 1.0e10 / C.speed_light_si
        np.testing.assert_allclose(C.thz2cm1, expected, rtol=1e-9)

    def test_thz2cm1_approx_33(self):
        # 1 THz ≈ 33.36 cm⁻¹
        np.testing.assert_allclose(C.thz2cm1, 33.36, rtol=1e-3)

    def test_wavenumber_positive_and_small(self):
        # 1 cm⁻¹ in Hartree is very small (~4.6e-6 Ha)
        assert 0 < C.wavenumber < 1e-4


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
