"""Tests A: convert_length_units in UnitCell.py.

Tests the conversion of lengths between bohr, angstrom, nm, um, mm, cm, m.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.UnitCell import convert_length_units


# ---------------------------------------------------------------------------
# A1: convert_length_units
# ---------------------------------------------------------------------------

class TestA1ConvertLengthUnits:
    """A1: convert_length_units — unit round-trips and known conversions."""

    def test_same_units_unchanged(self):
        assert convert_length_units(1.0, "angstrom", "angstrom") == pytest.approx(1.0)

    def test_same_units_unchanged_nm(self):
        assert convert_length_units(5.0, "nm", "nm") == pytest.approx(5.0)

    def test_bohr_to_angstrom_value(self):
        """1 bohr ≈ 0.5291772 Å."""
        result = convert_length_units(1.0, "bohr", "angstrom")
        assert result == pytest.approx(0.5291772, rel=1e-5)

    def test_angstrom_to_nm(self):
        """10 Å = 1 nm."""
        result = convert_length_units(10.0, "angstrom", "nm")
        assert result == pytest.approx(1.0, rel=1e-12)

    def test_nm_to_angstrom(self):
        """1 nm = 10 Å."""
        result = convert_length_units(1.0, "nm", "angstrom")
        assert result == pytest.approx(10.0, rel=1e-12)

    def test_angstrom_to_m(self):
        """1 Å = 1e-10 m."""
        result = convert_length_units(1.0, "angstrom", "m")
        assert result == pytest.approx(1.0e-10, rel=1e-10)

    def test_round_trip_angstrom_bohr(self):
        """Å → bohr → Å gives original value."""
        original = 5.0
        in_bohr = convert_length_units(original, "angstrom", "bohr")
        back = convert_length_units(in_bohr, "bohr", "angstrom")
        assert back == pytest.approx(original, rel=1e-10)

    def test_case_insensitivity(self):
        """'ANG' and 'ang' and 'Angstrom' are equivalent."""
        r1 = convert_length_units(1.0, "ANG", "NM")
        r2 = convert_length_units(1.0, "ang", "nm")
        r3 = convert_length_units(1.0, "Angstrom", "nm")
        assert r1 == pytest.approx(r2, rel=1e-12)
        assert r1 == pytest.approx(r3, rel=1e-12)

    def test_au_bohr_synonyms(self):
        """'a.u.', 'au', 'bohr' give the same conversion factor."""
        r1 = convert_length_units(1.0, "a.u.", "angstrom")
        r2 = convert_length_units(1.0, "au", "angstrom")
        r3 = convert_length_units(1.0, "bohr", "angstrom")
        assert r1 == pytest.approx(r2, rel=1e-12)
        assert r1 == pytest.approx(r3, rel=1e-12)

    def test_numpy_array_input(self):
        """Array input: each element is converted independently."""
        arr = np.array([1.0, 2.0, 3.0])
        result = convert_length_units(arr, "angstrom", "nm")
        expected = arr * 0.1
        np.testing.assert_allclose(result, expected, rtol=1e-12)

    def test_list_input_returns_array(self):
        """List input should also work and return an array."""
        result = convert_length_units([1.0, 2.0], "angstrom", "nm")
        assert len(result) == 2
        assert result[0] == pytest.approx(0.1, rel=1e-12)

    def test_cm_to_angstrom(self):
        """1 cm = 1e8 Å."""
        result = convert_length_units(1.0, "cm", "angstrom")
        assert result == pytest.approx(1.0e8, rel=1e-10)

    def test_mm_to_nm(self):
        """1 mm = 1e6 nm."""
        result = convert_length_units(1.0, "mm", "nm")
        assert result == pytest.approx(1.0e6, rel=1e-10)

    def test_um_to_angstrom(self):
        """1 µm = 1e4 Å."""
        result = convert_length_units(1.0, "um", "angstrom")
        assert result == pytest.approx(1.0e4, rel=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
