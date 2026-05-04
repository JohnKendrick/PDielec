"""Tests D1-D4: Element data tables spot-checked against reference values.

Plan sections D1-D4.

The tables atomic_number_to_element, element_to_atomic_number, average_masses,
isotope_masses (most-abundant isotope), and covalent_radii are verified against
IUPAC / NIST reference values.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
import PDielec.Constants as C


# ---------------------------------------------------------------------------
# D1: atomic_number_to_element and element_to_atomic_number
# ---------------------------------------------------------------------------

class TestD1ElementMappings:
    """D1: Forward and inverse element/atomic-number mappings."""

    def test_element_1_is_H(self):
        assert C.atomic_number_to_element[1] == "H"

    def test_element_6_is_C(self):
        assert C.atomic_number_to_element[6] == "C"

    def test_element_8_is_O(self):
        assert C.atomic_number_to_element[8] == "O"

    def test_element_26_is_Fe(self):
        assert C.atomic_number_to_element[26] == "Fe"

    def test_inverse_H_is_1(self):
        assert C.element_to_atomic_number["H"] == 1

    def test_inverse_C_is_6(self):
        assert C.element_to_atomic_number["C"] == 6

    def test_inverse_O_is_8(self):
        assert C.element_to_atomic_number["O"] == 8

    def test_forward_inverse_round_trip(self):
        for z, sym in C.atomic_number_to_element.items():
            if sym in C.element_to_atomic_number:
                assert C.element_to_atomic_number[sym] == z, (
                    f"Round-trip failed for Z={z}, sym={sym}"
                )

    def test_no_zero_or_negative_atomic_numbers_in_values(self):
        # All keys in element_to_atomic_number should be >= 0
        for sym, z in C.element_to_atomic_number.items():
            assert z >= 0, f"{sym} has atomic number {z} < 0"

    def test_common_elements_present(self):
        for sym in ("H", "C", "N", "O", "F", "Na", "Mg", "Si", "P", "S", "Cl",
                    "Ca", "Fe", "Cu", "Zn", "Br", "Ag", "Au"):
            assert sym in C.element_to_atomic_number, f"{sym} missing from element_to_atomic_number"


# ---------------------------------------------------------------------------
# D2: average_masses
# ---------------------------------------------------------------------------

class TestD2AverageMasses:
    """D2: Spot-check average atomic masses (in u)."""

    def test_H_mass(self):
        np.testing.assert_allclose(C.average_masses["H"], 1.008, atol=0.001)

    def test_C_mass(self):
        np.testing.assert_allclose(C.average_masses["C"], 12.011, atol=0.001)

    def test_O_mass(self):
        np.testing.assert_allclose(C.average_masses["O"], 15.999, atol=0.001)

    def test_Fe_mass(self):
        np.testing.assert_allclose(C.average_masses["Fe"], 55.845, atol=0.01)

    def test_Na_mass(self):
        np.testing.assert_allclose(C.average_masses["Na"], 22.990, atol=0.001)

    def test_all_defined_masses_positive(self):
        for sym, mass in C.average_masses.items():
            if mass is not None:
                assert mass > 0, f"{sym} has non-positive mass: {mass}"

    def test_H_mass_not_less_than_1(self):
        # Hydrogen is the lightest element, ~1 u
        assert C.average_masses["H"] >= 1.0


# ---------------------------------------------------------------------------
# D3: isotope_masses (most-abundant isotope, stored in isotope_masses dict)
# ---------------------------------------------------------------------------

class TestD3IsotopeMasses:
    """D3: Spot-check dominant isotope masses (isotope_masses dictionary)."""

    def test_H_most_abundant_isotope(self):
        # H-1: 1.00782503207 u
        np.testing.assert_allclose(C.isotope_masses["H"], 1.007825, atol=0.0001)

    def test_C_most_abundant_isotope(self):
        # C-12: exactly 12.000 u
        np.testing.assert_allclose(C.isotope_masses["C"], 12.0, atol=0.001)

    def test_O_most_abundant_isotope(self):
        # O-16: 15.994915 u
        np.testing.assert_allclose(C.isotope_masses["O"], 15.9949, atol=0.001)

    def test_all_defined_isotope_masses_positive(self):
        for sym, mass in C.isotope_masses.items():
            if mass is not None:
                assert mass > 0, f"{sym} has non-positive isotope mass: {mass}"


# ---------------------------------------------------------------------------
# D4: covalent_radii
# ---------------------------------------------------------------------------

class TestD4CovalentRadii:
    """D4: Covalent radii spot-check (in Å, from Cordero et al. 2008)."""

    def test_H_covalent_radius(self):
        # H: 0.31 Å
        np.testing.assert_allclose(C.covalent_radii["H"], 0.31, atol=0.05)

    def test_C_covalent_radius(self):
        # C: 0.76 Å
        np.testing.assert_allclose(C.covalent_radii["C"], 0.76, atol=0.05)

    def test_O_covalent_radius(self):
        # O: 0.66 Å
        np.testing.assert_allclose(C.covalent_radii["O"], 0.66, atol=0.05)

    def test_all_radii_non_negative(self):
        # "X" is a placeholder element with radius 0.0; all real elements should be > 0
        for sym, r in C.covalent_radii.items():
            if r is not None:
                assert r >= 0, f"{sym} has negative covalent radius: {r}"

    def test_real_elements_positive_radius(self):
        # Spot-check: real elements must have a strictly positive radius
        for sym in ("H", "C", "N", "O", "F", "Na", "Si", "Fe", "Cu"):
            r = C.covalent_radii[sym]
            assert r > 0, f"{sym} has non-positive covalent radius: {r}"

    def test_radii_increase_down_group_1(self):
        # Li < Na < K (alkali metals, group 1)
        assert C.covalent_radii["Li"] < C.covalent_radii["Na"]
        assert C.covalent_radii["Na"] < C.covalent_radii["K"]

    def test_radii_in_angstrom_range(self):
        # All covalent radii should be between 0.2 Å and 3.0 Å
        for sym, r in C.covalent_radii.items():
            if r is not None and r > 0:
                assert 0.2 <= r <= 3.0, f"{sym} covalent radius {r} Å out of expected range"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
