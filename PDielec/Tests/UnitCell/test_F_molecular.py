"""Tests F: Molecular contents and centre-of-mass in UnitCell.py.

Covers: calculate_molecular_contents, calculate_centre_of_mass.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.UnitCell import UnitCell


def _make_h2_cell():
    """A 10 Å cubic cell with one H2 molecule (H-H bond ~0.74 Å)."""
    cell = UnitCell(10.0, 10.0, 10.0, alpha=90.0, beta=90.0, gamma=90.0)
    cell.set_element_names(["H", "H"])
    cell.set_atomic_masses([1.00794, 1.00794])
    # Place both H atoms near the centre, separated by ~0.74 Å in x
    cell.set_fractional_coordinates([
        [0.50, 0.50, 0.50],
        [0.574, 0.50, 0.50],  # 0.74 Å apart in a 10 Å cell
    ])
    return cell


def _make_two_molecule_cell():
    """A 10 Å cubic cell with two isolated H atoms (each its own 'molecule')."""
    cell = UnitCell(10.0, 10.0, 10.0, alpha=90.0, beta=90.0, gamma=90.0)
    cell.set_element_names(["H", "H"])
    cell.set_atomic_masses([1.00794, 1.00794])
    # Atoms separated by 5 Å — much larger than any H-H bond
    cell.set_fractional_coordinates([
        [0.10, 0.10, 0.10],
        [0.60, 0.60, 0.60],
    ])
    return cell


# ---------------------------------------------------------------------------
# F1: calculate_molecular_contents
# ---------------------------------------------------------------------------

class TestF1MolecularContents:
    """F1: calculate_molecular_contents finds bonded molecules."""

    def test_h2_gives_one_molecule(self):
        """Two bonded H atoms → 1 molecule."""
        cell = _make_h2_cell()
        n_mol = cell.calculate_molecular_contents()
        assert n_mol == 1

    def test_two_isolated_atoms_give_two_molecules(self):
        """Two far-apart H atoms → 2 molecules."""
        cell = _make_two_molecule_cell()
        n_mol = cell.calculate_molecular_contents()
        assert n_mol == 2

    def test_returns_integer(self):
        cell = _make_h2_cell()
        result = cell.calculate_molecular_contents()
        assert isinstance(result, int)

    def test_nions_unchanged(self):
        """calculate_molecular_contents should not change the number of ions."""
        cell = _make_h2_cell()
        nions_before = cell.nions
        cell.calculate_molecular_contents()
        assert cell.nions == nions_before


# ---------------------------------------------------------------------------
# F2: calculate_centre_of_mass
# ---------------------------------------------------------------------------

class TestF2CentreOfMass:
    """F2: calculate_centre_of_mass returns correct xyz and fractional coordinates."""

    def setup_method(self):
        """Two equal-mass atoms at opposite corners → CoM at (0.5,0.5,0.5) frac."""
        self.cell = UnitCell(10.0, 10.0, 10.0, alpha=90.0, beta=90.0, gamma=90.0)
        self.cell.set_element_names(["H", "H"])
        self.cell.set_atomic_masses([1.0, 1.0])
        self.cell.set_fractional_coordinates([
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        ])

    def test_com_xyz_midpoint(self):
        cm = self.cell.calculate_centre_of_mass(output="xyz")
        np.testing.assert_allclose(cm, [5.0, 5.0, 5.0], atol=1e-10)

    def test_com_fractional_midpoint(self):
        cm = self.cell.calculate_centre_of_mass(output="abc")
        np.testing.assert_allclose(cm, [0.5, 0.5, 0.5], atol=1e-10)

    def test_com_mass(self):
        mass = self.cell.calculate_centre_of_mass(output="mass")
        assert mass == pytest.approx(2.0, rel=1e-12)

    def test_com_xyz_returns_array(self):
        cm = self.cell.calculate_centre_of_mass(output="xyz")
        assert hasattr(cm, "__len__")
        assert len(cm) == 3

    def test_com_subset_of_atoms(self):
        """CoM of just atom 0 is at the position of atom 0."""
        cm = self.cell.calculate_centre_of_mass(atom_list=[0], output="xyz")
        np.testing.assert_allclose(cm, [0.0, 0.0, 0.0], atol=1e-12)

    def test_com_unequal_masses(self):
        """Unequal masses shift the centre of mass toward the heavier atom."""
        cell = UnitCell(10.0, 10.0, 10.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["H", "O"])
        cell.set_atomic_masses([1.0, 16.0])
        cell.set_fractional_coordinates([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        cm = cell.calculate_centre_of_mass(output="abc")
        # CoM_x = (1*0 + 16*1)/(1+16) = 16/17
        assert cm[0] == pytest.approx(16.0 / 17.0, rel=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
