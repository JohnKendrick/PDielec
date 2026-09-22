"""Tests B: Lattice construction and parameter round-trips in UnitCell.py.

Covers: UnitCell.__init__ with lattice vectors and angles,
        set_lattice, convert_unitcell_to_abc, convert_abc_to_unitcell.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import math
import numpy as np
import pytest
from PDielec.UnitCell import UnitCell


# ---------------------------------------------------------------------------
# B1: Cubic lattice construction
# ---------------------------------------------------------------------------

class TestB1CubicLattice:
    """B1: Cubic cell constructed with (a,b,c,α,β,γ) = (5,5,5,90,90,90)."""

    def setup_method(self):
        self.cell = UnitCell(5.0, 5.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)

    def test_lattice_is_5_times_identity(self):
        np.testing.assert_allclose(self.cell.lattice, 5.0 * np.eye(3), atol=1e-10)

    def test_a_b_c_correct(self):
        assert self.cell.a == pytest.approx(5.0, rel=1e-10)
        assert self.cell.b == pytest.approx(5.0, rel=1e-10)
        assert self.cell.c == pytest.approx(5.0, rel=1e-10)

    def test_angles_correct(self):
        assert self.cell.alpha == pytest.approx(90.0, abs=1e-9)
        assert self.cell.beta  == pytest.approx(90.0, abs=1e-9)
        assert self.cell.gamma == pytest.approx(90.0, abs=1e-9)

    def test_convert_unitcell_to_abc_roundtrip(self):
        a, b, c, alpha, beta, gamma = self.cell.convert_unitcell_to_abc()
        assert a     == pytest.approx(5.0, rel=1e-10)
        assert b     == pytest.approx(5.0, rel=1e-10)
        assert c     == pytest.approx(5.0, rel=1e-10)
        assert alpha == pytest.approx(90.0, abs=1e-9)
        assert beta  == pytest.approx(90.0, abs=1e-9)
        assert gamma == pytest.approx(90.0, abs=1e-9)


# ---------------------------------------------------------------------------
# B2: Orthorhombic lattice
# ---------------------------------------------------------------------------

class TestB2OrthorhombicLattice:
    """B2: Orthorhombic cell (a=3, b=4, c=5, all 90°)."""

    def setup_method(self):
        self.cell = UnitCell(3.0, 4.0, 5.0, alpha=90.0, beta=90.0, gamma=90.0)

    def test_lattice_is_diagonal(self):
        np.testing.assert_allclose(self.cell.lattice, np.diag([3.0, 4.0, 5.0]), atol=1e-10)

    def test_a_b_c_correct(self):
        assert self.cell.a == pytest.approx(3.0, rel=1e-10)
        assert self.cell.b == pytest.approx(4.0, rel=1e-10)
        assert self.cell.c == pytest.approx(5.0, rel=1e-10)

    def test_round_trip(self):
        a, b, c, alpha, beta, gamma = self.cell.convert_unitcell_to_abc()
        assert a     == pytest.approx(3.0, rel=1e-10)
        assert b     == pytest.approx(4.0, rel=1e-10)
        assert c     == pytest.approx(5.0, rel=1e-10)
        assert alpha == pytest.approx(90.0, abs=1e-9)
        assert beta  == pytest.approx(90.0, abs=1e-9)
        assert gamma == pytest.approx(90.0, abs=1e-9)


# ---------------------------------------------------------------------------
# B3: Monoclinic lattice round-trip
# ---------------------------------------------------------------------------

class TestB3MonoclinicLattice:
    """B3: Monoclinic cell (alpha=gamma=90°, beta=110°)."""

    def setup_method(self):
        self.a0, self.b0, self.c0 = 5.0, 4.0, 6.0
        self.alpha0, self.beta0, self.gamma0 = 90.0, 110.0, 90.0
        self.cell = UnitCell(self.a0, self.b0, self.c0,
                             alpha=self.alpha0, beta=self.beta0, gamma=self.gamma0)

    def test_round_trip_abc(self):
        a, b, c, alpha, beta, gamma = self.cell.convert_unitcell_to_abc()
        assert a     == pytest.approx(self.a0, rel=1e-8)
        assert b     == pytest.approx(self.b0, rel=1e-8)
        assert c     == pytest.approx(self.c0, rel=1e-8)
        assert alpha == pytest.approx(self.alpha0, abs=1e-7)
        assert beta  == pytest.approx(self.beta0, abs=1e-7)
        assert gamma == pytest.approx(self.gamma0, abs=1e-7)

    def test_convert_abc_to_unitcell_round_trip(self):
        """abc → lattice → abc stays consistent."""
        cell2 = UnitCell()
        cell2.convert_abc_to_unitcell(self.a0, self.b0, self.c0,
                                       self.alpha0, self.beta0, self.gamma0)
        a, b, c, alpha, beta, gamma = cell2.convert_unitcell_to_abc()
        assert a     == pytest.approx(self.a0, rel=1e-8)
        assert beta  == pytest.approx(self.beta0, abs=1e-7)


# ---------------------------------------------------------------------------
# B4: set_lattice
# ---------------------------------------------------------------------------

class TestB4SetLattice:
    """B4: set_lattice stores the matrix and recomputes reciprocal lattice."""

    def test_set_cubic_lattice(self):
        cell = UnitCell()
        cell.set_lattice([[4, 0, 0], [0, 4, 0], [0, 0, 4]])
        np.testing.assert_allclose(cell.lattice, 4.0 * np.eye(3), atol=1e-12)

    def test_reciprocal_lattice_updated(self):
        """After set_lattice, reciprocal_lattice == inv(lattice)."""
        cell = UnitCell()
        lattice = np.diag([3.0, 4.0, 5.0])
        cell.set_lattice(lattice.tolist())
        expected_recip = np.linalg.inv(lattice)
        np.testing.assert_allclose(cell.reciprocal_lattice, expected_recip, atol=1e-12)

    def test_set_lattice_in_bohr(self):
        """set_lattice with units='bohr' converts to Angstrom internally."""
        cell = UnitCell()
        # 1 bohr ≈ 0.5291772 Å
        cell.set_lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]], units="bohr")
        assert cell.lattice[0, 0] == pytest.approx(0.5291772, rel=1e-4)


# ---------------------------------------------------------------------------
# B5: Default UnitCell
# ---------------------------------------------------------------------------

class TestB5DefaultUnitCell:
    """B5: Default UnitCell is a unit cube (1 Å sides)."""

    def test_default_is_identity(self):
        cell = UnitCell()
        np.testing.assert_allclose(cell.lattice, np.eye(3), atol=1e-12)

    def test_default_a_equals_one(self):
        cell = UnitCell()
        assert cell.a == pytest.approx(1.0, abs=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
