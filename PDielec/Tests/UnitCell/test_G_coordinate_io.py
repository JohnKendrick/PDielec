"""Tests G: Coordinate I/O helpers in UnitCell.py.

Covers: get_bounding_box, set_element_names, set_atomic_masses, get_atomic_masses.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.UnitCell import UnitCell


# ---------------------------------------------------------------------------
# G1: set_element_names / get_atomic_masses
# ---------------------------------------------------------------------------

class TestG1ElementsAndMasses:
    """G1: set_element_names and set/get_atomic_masses."""

    def test_set_element_names_stores_names(self):
        cell = UnitCell(4.0, 4.0, 4.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["Na", "Cl"])
        assert cell.element_names == ["Na", "Cl"]

    def test_set_element_names_cleans_symbols(self):
        """set_element_names should strip and capitalise symbols."""
        cell = UnitCell(4.0, 4.0, 4.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["na", "cl"])
        # After cleanup_symbol the names should match known element names
        for name in cell.element_names:
            assert name[0].isupper()

    def test_set_get_atomic_masses_roundtrip(self):
        cell = UnitCell(4.0, 4.0, 4.0, alpha=90.0, beta=90.0, gamma=90.0)
        masses_in = [22.99, 35.45]
        cell.set_element_names(["Na", "Cl"])
        cell.set_atomic_masses(masses_in)
        masses_out = cell.get_atomic_masses()
        np.testing.assert_allclose(masses_out, masses_in, rtol=1e-12)

    def test_masses_length_matches_elements(self):
        cell = UnitCell(4.0, 4.0, 4.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["C", "O", "O"])
        cell.set_atomic_masses([12.011, 15.999, 15.999])
        assert len(cell.get_atomic_masses()) == 3


# ---------------------------------------------------------------------------
# G2: get_bounding_box
# ---------------------------------------------------------------------------

class TestG2BoundingBox:
    """G2: get_bounding_box returns 8 corners and 12 edges."""

    def setup_method(self):
        self.cell = UnitCell(4.0, 5.0, 6.0, alpha=90.0, beta=90.0, gamma=90.0)

    def test_returns_three_items(self):
        result = self.cell.get_bounding_box()
        assert len(result) == 3

    def test_corners_count(self):
        corners, edges, labels = self.cell.get_bounding_box()
        assert len(corners) == 8

    def test_edges_count(self):
        corners, edges, labels = self.cell.get_bounding_box()
        assert len(edges) == 12

    def test_labels_count(self):
        corners, edges, labels = self.cell.get_bounding_box()
        assert len(labels) == 8

    def test_corners_are_3d(self):
        corners, edges, labels = self.cell.get_bounding_box()
        for corner in corners:
            assert len(corner) == 3

    def test_origin_option(self):
        """Specifying originABC should shift the box centre."""
        corners_default, _, _ = self.cell.get_bounding_box()
        corners_shifted, _, _ = self.cell.get_bounding_box(originABC=[0.5, 0.5, 0.5])
        # Centres should differ (shifted vs default)
        centre_default = np.mean(corners_default, axis=0)
        centre_shifted = np.mean(corners_shifted, axis=0)
        # They should not be identical
        assert not np.allclose(centre_default, centre_shifted, atol=1e-3)


# ---------------------------------------------------------------------------
# G3: convert_unitcell_to_abc / convert_abc_to_unitcell
# ---------------------------------------------------------------------------

class TestG3AbcRoundTrip:
    """G3: Triclinic lattice round-trip through convert_abc_to_unitcell."""

    def test_triclinic_abc_round_trip(self):
        a0, b0, c0 = 5.0, 6.0, 7.0
        al0, be0, ga0 = 80.0, 85.0, 95.0
        cell = UnitCell(a0, b0, c0, alpha=al0, beta=be0, gamma=ga0)
        a, b, c, al, be, ga = cell.convert_unitcell_to_abc()
        assert a  == pytest.approx(a0, rel=1e-8)
        assert b  == pytest.approx(b0, rel=1e-8)
        assert c  == pytest.approx(c0, rel=1e-8)
        assert al == pytest.approx(al0, abs=1e-7)
        assert be == pytest.approx(be0, abs=1e-7)
        assert ga == pytest.approx(ga0, abs=1e-7)

    def test_convert_abc_to_unitcell_returns_lattice(self):
        cell = UnitCell()
        lat = cell.convert_abc_to_unitcell(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
        np.testing.assert_allclose(lat, 5.0 * np.eye(3), atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
