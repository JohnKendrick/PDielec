"""Tests D: Volume and density calculations in UnitCell.py.

Covers: get_volume, get_density.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.UnitCell import UnitCell


# ---------------------------------------------------------------------------
# D1: get_volume
# ---------------------------------------------------------------------------

class TestD1GetVolume:
    """D1: get_volume returns the cell volume in the requested units."""

    def test_cubic_volume_angstrom(self):
        """5 Å cube → 125 Å³."""
        cell = UnitCell(5.0, 5.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)
        vol = cell.get_volume(units="Angstrom")
        assert vol == pytest.approx(125.0, rel=1e-10)

    def test_orthorhombic_volume_angstrom(self):
        """3 × 4 × 5 cell → 60 Å³."""
        cell = UnitCell(3.0, 4.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)
        vol = cell.get_volume(units="Angstrom")
        assert vol == pytest.approx(60.0, rel=1e-10)

    def test_volume_scales_with_units(self):
        """Volume in nm³ should equal volume in Å³ × (0.1)³ = 1e-3."""
        cell = UnitCell(5.0, 5.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)
        vol_ang = cell.get_volume(units="Angstrom")
        vol_nm  = cell.get_volume(units="nm")
        assert vol_nm == pytest.approx(vol_ang * (0.1 ** 3), rel=1e-9)

    def test_volume_positive(self):
        cell = UnitCell(3.0, 4.0, 5.0, alpha=90.0, beta=80.0, gamma=90.0)
        assert cell.get_volume(units="Angstrom") > 0.0

    def test_monoclinic_volume(self):
        """For monoclinic cell (alpha=gamma=90, beta), volume = a*b*c*sin(beta)."""
        a, b, c, beta = 3.0, 4.0, 5.0, 100.0
        cell = UnitCell(a, b, c, alpha=90.0, beta=beta, gamma=90.0)
        expected = a * b * c * np.sin(np.radians(beta))
        vol = cell.get_volume(units="Angstrom")
        assert vol == pytest.approx(expected, rel=1e-8)


# ---------------------------------------------------------------------------
# D2: get_density
# ---------------------------------------------------------------------------

class TestD2GetDensity:
    """D2: get_density returns 1.0 when no atoms are set; otherwise mass/volume."""

    def test_no_atoms_returns_one(self):
        """Without atoms, density defaults to 1.0."""
        cell = UnitCell(5.0, 5.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)
        assert cell.get_density() == pytest.approx(1.0)

    def test_density_positive_with_atoms(self):
        """With atoms set, density should be positive."""
        cell = UnitCell(4.0, 4.0, 4.0, alpha=90.0, beta=90.0, gamma=90.0)
        cell.set_element_names(["Na", "Cl"])
        cell.set_fractional_coordinates([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        cell.set_atomic_masses([22.99, 35.45])
        density = cell.get_density(units="cm")
        assert density > 0.0

    def test_density_units_consistency(self):
        """Density in g/cm³ from a cell with known masses and volume."""
        from PDielec.Constants import avogadro_si
        a = 4.0  # Angstrom
        cell = UnitCell(a, a, a, alpha=90.0, beta=90.0, gamma=90.0)
        # One atom of mass 20 amu
        mass_amu = 20.0
        cell.set_element_names(["Ne"])
        cell.set_fractional_coordinates([[0.0, 0.0, 0.0]])
        cell.set_atomic_masses([mass_amu])
        density = cell.get_density(units="cm")
        # volume in cm³: (4e-8)³ = 64e-24 cm³
        vol_cm3 = (a * 1e-8) ** 3
        mass_g = mass_amu * 1.66053906660e-24  # 1 amu in grams
        expected = mass_g / vol_cm3
        assert density == pytest.approx(expected, rel=1e-5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
