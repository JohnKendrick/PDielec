"""Tests H: CIF file I/O in UnitCell.py.

Covers: write_cif (to file and to string buffer).
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import io
import tempfile
import pytest
from PDielec.UnitCell import UnitCell


def _make_nacl_cell():
    """NaCl-like rock-salt cell (simplified 2-atom)."""
    a = 5.64  # Angstrom
    cell = UnitCell(a, a, a, alpha=90.0, beta=90.0, gamma=90.0)
    cell.set_element_names(["Na", "Cl"])
    cell.set_atomic_masses([22.99, 35.45])
    cell.set_fractional_coordinates([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ])
    return cell


# ---------------------------------------------------------------------------
# H1: write_cif to a StringIO buffer
# ---------------------------------------------------------------------------

class TestH1WriteCif:
    """H1: write_cif produces valid CIF text."""

    def setup_method(self):
        self.cell = _make_nacl_cell()
        self.buf = io.StringIO()
        self.cell.write_cif(filedescriptor=self.buf)
        self.buf.seek(0)
        self.content = self.buf.read()

    def test_cif_contains_data_line(self):
        assert "data_" in self.content

    def test_cif_contains_cell_length_a(self):
        assert "_cell_length_a" in self.content

    def test_cif_contains_cell_length_b(self):
        assert "_cell_length_b" in self.content

    def test_cif_contains_cell_length_c(self):
        assert "_cell_length_c" in self.content

    def test_cif_contains_cell_angle_alpha(self):
        assert "_cell_angle_alpha" in self.content

    def test_cif_contains_atom_site_fract_x(self):
        assert "_atom_site_fract_x" in self.content

    def test_cif_contains_end_marker(self):
        assert "#END" in self.content

    def test_cif_contains_element_na(self):
        assert "Na" in self.content

    def test_cif_contains_element_cl(self):
        assert "Cl" in self.content

    def test_cif_contains_correct_a_value(self):
        """The cell length a (5.64) must appear in the CIF."""
        assert "5.640000" in self.content

    def test_cif_loop_block(self):
        assert "loop_" in self.content


# ---------------------------------------------------------------------------
# H2: write_cif to a file
# ---------------------------------------------------------------------------

class TestH2WriteCifToFile:
    """H2: write_cif can write to a named file."""

    def test_file_is_created(self, tmp_path):
        cell = _make_nacl_cell()
        cif_file = str(tmp_path / "test.cif")
        cell.write_cif(filename=cif_file)
        assert os.path.exists(cif_file)

    def test_file_content_contains_end(self, tmp_path):
        cell = _make_nacl_cell()
        cif_file = str(tmp_path / "test.cif")
        cell.write_cif(filename=cif_file)
        with open(cif_file) as f:
            content = f.read()
        assert "#END" in content

    def test_file_with_description(self, tmp_path):
        cell = _make_nacl_cell()
        cif_file = str(tmp_path / "described.cif")
        cell.write_cif(filename=cif_file, description="NaCl_test")
        with open(cif_file) as f:
            content = f.read()
        assert "NaCl_test" in content


# ---------------------------------------------------------------------------
# H3: CIF content correctness
# ---------------------------------------------------------------------------

class TestH3CifContent:
    """H3: Check that CIF numeric values are correct."""

    def test_cell_volume_line(self):
        cell = _make_nacl_cell()
        buf = io.StringIO()
        cell.write_cif(filedescriptor=buf)
        content = buf.getvalue()
        # Volume of 5.64³ ≈ 179.406...
        assert "_cell_volume" in content
        # Check the value is present (rounded to 6 dp)
        lines = content.splitlines()
        vol_line = next(l for l in lines if "_cell_volume" in l)
        vol_value = float(vol_line.split()[-1])
        assert vol_value == pytest.approx(5.64 ** 3, rel=1e-4)

    def test_fractional_coordinates_present(self):
        cell = _make_nacl_cell()
        buf = io.StringIO()
        cell.write_cif(filedescriptor=buf)
        content = buf.getvalue()
        # Origin atom should have 0.000000 coordinates
        assert "0.000000" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
