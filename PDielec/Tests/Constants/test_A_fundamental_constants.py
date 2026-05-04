"""Tests A1-A2: Fundamental physical constants from NIST CODATA 2018.

Plan sections A1, A2.

Verifies that the numeric literals in Constants.py match the CODATA 2018
recommended values and that all fundamental constants have the expected type
and sign.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
import PDielec.Constants as C


# ---------------------------------------------------------------------------
# A1: Fundamental constants from NIST CODATA 2018
# ---------------------------------------------------------------------------

class TestA1FundamentalValues:
    """A1: Each constant must match its CODATA 2018 value within rtol=1e-9."""

    rtol = 1e-9

    def test_speed_of_light(self):
        # Exact by definition since 2019 SI
        assert abs(C.speed_light_si - 299792458.0) < 1e-3

    def test_planck_constant(self):
        # Exact by definition since 2019 SI: 6.62607015e-34 J·Hz⁻¹
        np.testing.assert_allclose(C.planck_si, 6.62607015e-34, rtol=self.rtol)

    def test_elementary_charge(self):
        # Exact by definition since 2019 SI: 1.602176634e-19 C
        np.testing.assert_allclose(C.elementary_charge_si, 1.602176634e-19, rtol=self.rtol)

    def test_avogadro_number(self):
        # Exact by definition since 2019 SI: 6.02214076e23 mol⁻¹
        np.testing.assert_allclose(C.avogadro_si, 6.02214076e23, rtol=self.rtol)

    def test_boltzmann_constant(self):
        # Exact by definition since 2019 SI: 1.380649e-23 J/K
        np.testing.assert_allclose(C.boltzmann_si, 1.380649e-23, rtol=self.rtol)

    def test_vacuum_permeability(self):
        # CODATA 2018: 1.25663706212e-6 N/A²
        np.testing.assert_allclose(C.mu_0_si, 1.25663706212e-6, rtol=1e-8)

    def test_electron_mass(self):
        # CODATA 2018: 9.1093837015e-31 kg
        np.testing.assert_allclose(C.electron_mass_si, 9.1093837015e-31, rtol=self.rtol)

    def test_hartree_to_ev(self):
        # CODATA 2018: 27.211386245988 eV
        np.testing.assert_allclose(C.hartree2ev, 27.211386245988, rtol=self.rtol)

    def test_bohr_radius_si(self):
        # CODATA 2018: 5.29177210903e-11 m
        np.testing.assert_allclose(C.bohr_si, 5.29177210903e-11, rtol=self.rtol)


# ---------------------------------------------------------------------------
# A2: Types and signs
# ---------------------------------------------------------------------------

class TestA2TypesAndSigns:
    """A2: All fundamental constants must be numeric scalars and positive."""

    _constants = [
        "speed_light_si",
        "planck_si",
        "elementary_charge_si",
        "avogadro_si",
        "boltzmann_si",
        "mu_0_si",
        "electron_mass_si",
        "hartree2ev",
        "bohr_si",
    ]

    def test_all_are_numeric(self):
        for name in self._constants:
            val = getattr(C, name)
            assert isinstance(val, (int, float, np.floating)), (
                f"{name} is not numeric: {type(val)}"
            )

    def test_all_are_positive(self):
        for name in self._constants:
            val = getattr(C, name)
            assert val > 0, f"{name} = {val} is not positive"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
