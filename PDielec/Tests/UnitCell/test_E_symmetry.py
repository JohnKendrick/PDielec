"""Tests E: Symmetry finding in UnitCell.py.

Covers: find_symmetry, get_atomic_numbers, set_element_names.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import pytest
from PDielec.UnitCell import UnitCell


def _make_fcc_al():
    """Build a conventional FCC Al cell (Fm-3m, #225)."""
    a = 4.046  # Angstrom — experimental Al lattice constant
    cell = UnitCell(a, a, a, alpha=90.0, beta=90.0, gamma=90.0)
    cell.set_element_names(["Al", "Al", "Al", "Al"])
    cell.set_atomic_masses([26.98, 26.98, 26.98, 26.98])
    frac = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ]
    cell.set_fractional_coordinates(frac)
    return cell


def _make_sc_si():
    """Build a simple-cubic Si cell (P m -3 m, #221) — fictitious but symmetry is clear."""
    a = 5.0
    cell = UnitCell(a, a, a, alpha=90.0, beta=90.0, gamma=90.0)
    cell.set_element_names(["Si"])
    cell.set_atomic_masses([28.09])
    cell.set_fractional_coordinates([[0.0, 0.0, 0.0]])
    return cell


# ---------------------------------------------------------------------------
# E1: Empty cell symmetry
# ---------------------------------------------------------------------------

class TestE1EmptyCellSymmetry:
    """E1: An empty UnitCell (no atoms) should return P1 (space group 1)."""

    def test_empty_returns_P1(self):
        cell = UnitCell(5.0, 5.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)
        symbol, number = cell.find_symmetry()
        assert symbol == "P 1"
        assert number == 1


# ---------------------------------------------------------------------------
# E2: Simple cubic single-atom cell
# ---------------------------------------------------------------------------

class TestE2SingleAtomCubic:
    """E2: One atom at origin in a cubic cell → P m -3 m (#221)."""

    def test_sc_symmetry_symbol(self):
        cell = _make_sc_si()
        symbol, number = cell.find_symmetry()
        assert number == 221

    def test_returns_tuple(self):
        cell = _make_sc_si()
        result = cell.find_symmetry()
        assert len(result) == 2


# ---------------------------------------------------------------------------
# E3: FCC cell
# ---------------------------------------------------------------------------

class TestE3FccCell:
    """E3: Conventional FCC Al → Fm-3m (#225)."""

    def test_fcc_number(self):
        cell = _make_fcc_al()
        symbol, number = cell.find_symmetry()
        assert number == 225

    def test_fcc_symbol_contains_F(self):
        cell = _make_fcc_al()
        symbol, number = cell.find_symmetry()
        assert symbol.startswith("F")


# ---------------------------------------------------------------------------
# E4: get_atomic_numbers
# ---------------------------------------------------------------------------

class TestE4GetAtomicNumbers:
    """E4: get_atomic_numbers returns correct Z values for element names."""

    def test_hydrogen_z_is_1(self):
        cell = UnitCell(3.0, 3.0, 3.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["H"])
        cell.set_fractional_coordinates([[0.0, 0.0, 0.0]])
        nos = cell.get_atomic_numbers()
        assert nos == [1]

    def test_carbon_z_is_6(self):
        cell = UnitCell(3.0, 3.0, 3.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["C"])
        cell.set_fractional_coordinates([[0.0, 0.0, 0.0]])
        nos = cell.get_atomic_numbers()
        assert nos == [6]

    def test_mixed_elements(self):
        cell = UnitCell(4.0, 4.0, 4.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["Na", "Cl"])
        cell.set_fractional_coordinates([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        nos = cell.get_atomic_numbers()
        assert nos[0] == 11  # Na
        assert nos[1] == 17  # Cl

    def test_length_matches_nions(self):
        cell = _make_fcc_al()
        nos = cell.get_atomic_numbers()
        assert len(nos) == cell.nions


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
