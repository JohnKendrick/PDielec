"""Tests B1-B8: Derived constants are consistent with their fundamental definitions.

Plan sections B1-B8.

Each derived constant is re-derived from the fundamental constants and compared
to the value stored in Constants.py.  All comparisons use rtol=1e-9 unless noted.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
import PDielec.Constants as C


RTOL = 1e-9


# ---------------------------------------------------------------------------
# B1: hbar_si
# ---------------------------------------------------------------------------

class TestB1Hbar:
    """B1: hbar_si = planck_si / (2*pi)."""

    def test_hbar_from_planck(self):
        expected = C.planck_si / (2.0 * np.pi)
        np.testing.assert_allclose(C.hbar_si, expected, rtol=RTOL)

    def test_hbar_positive(self):
        assert C.hbar_si > 0


# ---------------------------------------------------------------------------
# B2: molar_gas_si
# ---------------------------------------------------------------------------

class TestB2MolarGas:
    """B2: molar_gas_si = boltzmann_si * avogadro_si."""

    def test_molar_gas_from_fundamentals(self):
        expected = C.boltzmann_si * C.avogadro_si
        np.testing.assert_allclose(C.molar_gas_si, expected, rtol=RTOL)

    def test_molar_gas_approx_8314(self):
        # R ≈ 8.314 J/(mol·K) — sanity check on magnitude
        assert 8.0 < C.molar_gas_si < 8.5


# ---------------------------------------------------------------------------
# B3: fine_structure
# ---------------------------------------------------------------------------

class TestB3FineStructure:
    """B3: fine_structure ≈ 1/137."""

    def test_fine_structure_approx_1_over_137(self):
        np.testing.assert_allclose(C.fine_structure, 1.0 / 137.0, rtol=1e-3)

    def test_fine_structure_from_fundamentals(self):
        # alpha = e^2 / (4*pi*eps0*hbar*c)
        # Using stored eps0 (epsilon_0_si)
        expected = (
            C.elementary_charge_si ** 2
            / (4.0 * np.pi * C.epsilon_0_si * C.hbar_si * C.speed_light_si)
        )
        np.testing.assert_allclose(C.fine_structure, expected, rtol=1e-8)

    def test_fine_structure_positive(self):
        assert C.fine_structure > 0


# ---------------------------------------------------------------------------
# B4: amu
# ---------------------------------------------------------------------------

class TestB4Amu:
    """B4: amu = 1e-3 / avogadro_si / electron_mass_si (in atomic mass units)."""

    def test_amu_from_definition(self):
        # The module stores amu in atomic units (electron masses)
        expected = 1.0e-3 / C.avogadro_si / C.electron_mass_si
        np.testing.assert_allclose(C.amu, expected, rtol=RTOL)

    def test_amu_positive(self):
        assert C.amu > 0

    def test_amu_roughly_1823(self):
        # 1 AMU ≈ 1822.888 electron masses
        assert 1800 < C.amu < 1850


# ---------------------------------------------------------------------------
# B5: wavenumber
# ---------------------------------------------------------------------------

class TestB5Wavenumber:
    """B5: wavenumber is the energy of 1 cm⁻¹ in Hartree atomic units."""

    def test_wavenumber_positive(self):
        assert C.wavenumber > 0

    def test_wavenumber_from_hertz(self):
        # wavenumber = hertz * speed_light_si * 100 (as coded)
        expected = C.hertz * C.speed_light_si * 100.0
        np.testing.assert_allclose(C.wavenumber, expected, rtol=RTOL)

    def test_wavenumber_consistent_with_planck(self):
        # 1 cm⁻¹ photon energy in SI: h * c * 100
        energy_si = C.planck_si * C.speed_light_si * 100.0  # Joules
        # wavenumber is in Hartree; hartree (in Joules) = 1/joule
        hartree_J = 1.0 / C.joule
        energy_hartree = energy_si * C.joule
        np.testing.assert_allclose(C.wavenumber, energy_hartree, rtol=1e-8)


# ---------------------------------------------------------------------------
# B6: angs2bohr
# ---------------------------------------------------------------------------

class TestB6Angs2Bohr:
    """B6: angs2bohr = 1e-10 / bohr_si."""

    def test_angs2bohr_from_bohr_si(self):
        expected = 1.0e-10 / C.bohr_si
        np.testing.assert_allclose(C.angs2bohr, expected, rtol=RTOL)

    def test_angs2bohr_approx_1_889(self):
        # 1 Å = 1/0.529177 Bohr ≈ 1.889 Bohr
        np.testing.assert_allclose(C.angs2bohr, 1.0 / 0.529177, rtol=1e-4)

    def test_round_trip_angstrom_bohr(self):
        # angs2bohr * bohr_si = 1e-10 (1 Angstrom in metres)
        np.testing.assert_allclose(C.angs2bohr * C.bohr_si, 1.0e-10, rtol=RTOL)


# ---------------------------------------------------------------------------
# B7: hartree (in SI joules)
# ---------------------------------------------------------------------------

class TestB7Hartree:
    """B7: hartree ≈ hartree2ev * elementary_charge_si (in Joules)."""

    def test_hartree_from_ev(self):
        expected = C.hartree2ev * C.elementary_charge_si
        # hartree = 1/joule in the code; compare to expected J value
        hartree_J = 1.0 / C.joule
        np.testing.assert_allclose(hartree_J, expected, rtol=1e-8)

    def test_hartree_magnitude(self):
        # 1 Hartree ≈ 4.3597e-18 J
        hartree_J = 1.0 / C.joule
        np.testing.assert_allclose(hartree_J, 4.3597e-18, rtol=1e-3)


# ---------------------------------------------------------------------------
# B8: d2byamuang2
# ---------------------------------------------------------------------------

class TestB8D2ByAmuAng2:
    """B8: d2byamuang2 = debye² / (angstrom² * amu)."""

    def test_d2byamuang2_positive(self):
        assert C.d2byamuang2 > 0

    def test_d2byamuang2_from_components(self):
        expected = C.debye * C.debye / (C.angstrom * C.angstrom * C.amu)
        np.testing.assert_allclose(C.d2byamuang2, expected, rtol=RTOL)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
