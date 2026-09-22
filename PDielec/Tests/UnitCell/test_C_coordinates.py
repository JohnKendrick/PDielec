"""Tests C: Coordinate conversion in UnitCell.py.

Covers: convert_abc_to_xyz, convert_xyz_to_abc, convert_hkl_to_xyz,
        set_fractional_coordinates, set_xyz_coordinates.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.UnitCell import UnitCell


# ---------------------------------------------------------------------------
# C1: convert_abc_to_xyz for a cubic cell
# ---------------------------------------------------------------------------

class TestC1AbcToXyzCubic:
    """C1: Fractional → Cartesian in a cubic 5 Å cell."""

    def setup_method(self):
        self.cell = UnitCell(5.0, 5.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)

    def test_origin_stays_origin(self):
        xyz = self.cell.convert_abc_to_xyz([0.0, 0.0, 0.0])
        np.testing.assert_allclose(xyz, [0.0, 0.0, 0.0], atol=1e-12)

    def test_corner_at_one(self):
        xyz = self.cell.convert_abc_to_xyz([1.0, 1.0, 1.0])
        np.testing.assert_allclose(xyz, [5.0, 5.0, 5.0], atol=1e-12)

    def test_face_centre(self):
        xyz = self.cell.convert_abc_to_xyz([0.5, 0.5, 0.0])
        np.testing.assert_allclose(xyz, [2.5, 2.5, 0.0], atol=1e-12)

    def test_a_axis_unit(self):
        """abc=(1,0,0) → xyz = first lattice vector."""
        xyz = self.cell.convert_abc_to_xyz([1.0, 0.0, 0.0])
        np.testing.assert_allclose(xyz, self.cell.lattice[0], atol=1e-12)

    def test_b_axis_unit(self):
        xyz = self.cell.convert_abc_to_xyz([0.0, 1.0, 0.0])
        np.testing.assert_allclose(xyz, self.cell.lattice[1], atol=1e-12)

    def test_c_axis_unit(self):
        xyz = self.cell.convert_abc_to_xyz([0.0, 0.0, 1.0])
        np.testing.assert_allclose(xyz, self.cell.lattice[2], atol=1e-12)


# ---------------------------------------------------------------------------
# C2: Round-trip abc → xyz → abc
# ---------------------------------------------------------------------------

class TestC2RoundTrip:
    """C2: convert_abc_to_xyz followed by convert_xyz_to_abc recovers original."""

    def setup_method(self):
        # Use a monoclinic cell to make it non-trivial
        self.cell = UnitCell(4.0, 5.0, 6.0, alpha=90.0, beta=105.0, gamma=90.0)

    def test_roundtrip_origin(self):
        abc_in = [0.0, 0.0, 0.0]
        abc_out = self.cell.convert_xyz_to_abc(self.cell.convert_abc_to_xyz(abc_in))
        np.testing.assert_allclose(abc_out, abc_in, atol=1e-12)

    def test_roundtrip_general_point(self):
        abc_in = [0.3, 0.7, 0.2]
        xyz = self.cell.convert_abc_to_xyz(abc_in)
        abc_out = self.cell.convert_xyz_to_abc(xyz)
        np.testing.assert_allclose(abc_out, abc_in, atol=1e-10)

    def test_roundtrip_corner(self):
        abc_in = [1.0, 1.0, 1.0]
        xyz = self.cell.convert_abc_to_xyz(abc_in)
        abc_out = self.cell.convert_xyz_to_abc(xyz)
        np.testing.assert_allclose(abc_out, abc_in, atol=1e-10)


# ---------------------------------------------------------------------------
# C3: convert_hkl_to_xyz
# ---------------------------------------------------------------------------

class TestC3HklToXyz:
    """C3: convert_hkl_to_xyz maps Miller indices to reciprocal-space vectors."""

    def setup_method(self):
        self.cell = UnitCell(3.0, 4.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)

    def test_hkl_100(self):
        """(1,0,0) reciprocal vector for orthorhombic cell."""
        xyz = self.cell.convert_hkl_to_xyz([1, 0, 0])
        # For orthorhombic: reciprocal lattice = inv(diag([3,4,5])) = diag([1/3,1/4,1/5])
        expected = np.dot([1, 0, 0], np.linalg.inv(self.cell.lattice).T)
        np.testing.assert_allclose(xyz, expected, atol=1e-12)

    def test_hkl_000_is_origin(self):
        xyz = self.cell.convert_hkl_to_xyz([0, 0, 0])
        np.testing.assert_allclose(xyz, [0.0, 0.0, 0.0], atol=1e-12)

    def test_hkl_110(self):
        xyz = self.cell.convert_hkl_to_xyz([1, 1, 0])
        expected = np.dot([1, 1, 0], np.linalg.inv(self.cell.lattice).T)
        np.testing.assert_allclose(xyz, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# C4: set_fractional_coordinates and set_xyz_coordinates
# ---------------------------------------------------------------------------

class TestC4SetCoordinates:
    """C4: setting coordinates updates both fractional and Cartesian attributes."""

    def setup_method(self):
        self.cell = UnitCell(4.0, 4.0, 4.0, alpha=90.0, beta=90.0, gamma=90.0)

    def test_set_fractional_updates_xyz(self):
        frac = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        self.cell.set_fractional_coordinates(frac)
        expected_xyz = np.array([[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]])
        np.testing.assert_allclose(self.cell.xyz_coordinates, expected_xyz, atol=1e-12)

    def test_set_fractional_updates_nions(self):
        frac = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        self.cell.set_fractional_coordinates(frac)
        assert self.cell.nions == 3

    def test_set_xyz_updates_fractional(self):
        xyz = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
        self.cell.set_xyz_coordinates(xyz)
        expected_frac = np.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0]])
        np.testing.assert_allclose(self.cell.fractional_coordinates, expected_frac, atol=1e-12)

    def test_set_xyz_updates_nions(self):
        self.cell.set_xyz_coordinates([[1.0, 0.0, 0.0]])
        assert self.cell.nions == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
