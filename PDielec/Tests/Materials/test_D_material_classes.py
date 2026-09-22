"""Tests D1-D8: Material subclasses (Constant, NoMatrix, External,
DrudeLorentz, FPSQ, Sellmeier, Tabulated) and Material base helpers.

Plan sections D1–D8.

Tests verify:
  - Correct name, density, type string for each subclass.
  - Correct scalar/tensor classification delegated to the underlying
    DielectricFunction.
  - get_permittivity_object() returns the expected DielectricFunction subclass.
  - get_optical_permittivity() returns physically meaningful values.
  - Tabulated dispatches to TabulateScalar / Tabulate3 / Tabulate6 depending on
    the shape of the permittivity array.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec import DielectricFunction as DF
from PDielec.Materials import (
    Constant, NoMatrix, External, DrudeLorentz, FPSQ, Sellmeier, Tabulated, Material
)


# ---------------------------------------------------------------------------
# D1: Constant
# ---------------------------------------------------------------------------

class TestD1Constant:
    """D1: Constant material — scalar permittivity."""

    def test_get_name(self):
        m = Constant("water", permittivity=1.78, density=1.0)
        assert m.get_name() == "water"

    def test_get_density(self):
        m = Constant("water", permittivity=1.78, density=1.0)
        assert abs(m.get_density() - 1.0) < 1e-12

    def test_is_scalar_true(self):
        m = Constant("test", permittivity=2.0, density=1.0)
        assert m.is_scalar() is True

    def test_is_tensor_false(self):
        m = Constant("test", permittivity=2.0, density=1.0)
        assert m.is_tensor() is False

    def test_permittivity_function_returns_constant(self):
        val = complex(2.25, 0.1)
        m = Constant("test", permittivity=val, density=1.0)
        fn = m.get_permittivity_function()
        for v in [0.0, 100.0, 500.0, 1000.0]:
            assert fn(v) == val

    def test_get_optical_permittivity_returns_float(self):
        m = Constant("test", permittivity=complex(3.0, 0.5), density=1.0)
        op = m.get_optical_permittivity()
        assert isinstance(op, float)
        assert abs(op - 3.0) < 1e-12

    def test_type_string(self):
        m = Constant("test", permittivity=1.5, density=1.0)
        assert m.type == "Constant permittivity"

    def test_complex_permittivity_stored(self):
        val = complex(2.0, 0.3)
        m = Constant("test", permittivity=val, density=1.0)
        fn = m.get_permittivity_function()
        assert fn(0.0) == val


# ---------------------------------------------------------------------------
# D2: NoMatrix
# ---------------------------------------------------------------------------

class TestD2NoMatrix:
    """D2: NoMatrix — sentinel for bypassing EMT."""

    def test_get_name(self):
        assert NoMatrix().get_name() == "none"

    def test_get_density_zero(self):
        assert NoMatrix().get_density() == 0.0

    def test_is_scalar_true(self):
        assert NoMatrix().is_scalar() is True

    def test_permittivity_is_one(self):
        fn = NoMatrix().get_permittivity_function()
        assert fn(0.0) == 1.0
        assert fn(500.0) == 1.0

    def test_type_string(self):
        assert NoMatrix().type == "No matrix (DFT only)"


# ---------------------------------------------------------------------------
# D3: External
# ---------------------------------------------------------------------------

class TestD3External:
    """D3: External material wraps a caller-supplied DielectricFunction."""

    def test_is_scalar_delegates(self):
        obj = DF.ConstantScalar(2.0)
        m = External("ext", permittivity_object=obj, density=1.5)
        assert m.is_scalar() is True

    def test_is_tensor_delegates(self):
        obj = DF.ConstantTensor(2.0)
        m = External("ext", permittivity_object=obj, density=1.5)
        assert m.is_tensor() is True

    def test_get_permittivity_object_identity(self):
        obj = DF.ConstantScalar(3.0)
        m = External("ext", permittivity_object=obj, density=1.5)
        assert m.get_permittivity_object() is obj

    def test_type_string(self):
        obj = DF.ConstantScalar(1.0)
        m = External("ext", permittivity_object=obj, density=1.0)
        assert m.type == "External permittivity"

    def test_get_name(self):
        obj = DF.ConstantScalar(1.0)
        m = External("mymat", permittivity_object=obj, density=1.0)
        assert m.get_name() == "mymat"


# ---------------------------------------------------------------------------
# D4: DrudeLorentz (Material)
# ---------------------------------------------------------------------------

class TestD4DrudeLorentzMaterial:
    """D4: DrudeLorentz Material — tensor, eps_inf, correct DF type."""

    _eps_inf = np.diag([2.0, 2.0, 2.0])

    def _make(self):
        return DrudeLorentz(
            "AlN",
            epsinf=self._eps_inf.tolist(),
            omegas=[[400.0], [400.0], [400.0]],
            strengths=[[50.0], [50.0], [50.0]],
            gammas=[[5.0],   [5.0],   [5.0]],
            density=3.26,
        )

    def test_is_tensor_true(self):
        assert self._make().is_tensor() is True

    def test_is_scalar_false(self):
        assert self._make().is_scalar() is False

    def test_permittivity_object_type(self):
        m = self._make()
        assert isinstance(m.get_permittivity_object(), DF.DrudeLorentz)

    def test_epsilon_infinity_set(self):
        m = self._make()
        obj = m.get_permittivity_object()
        np.testing.assert_allclose(obj.epsilon_infinity, self._eps_inf, atol=1e-12)

    def test_type_string(self):
        assert self._make().type == "Drude-Lorentz"

    def test_get_optical_permittivity_equals_eps_inf(self):
        m = self._make()
        op = m.get_optical_permittivity()
        np.testing.assert_allclose(op, self._eps_inf, atol=1e-12)

    def test_density_set(self):
        assert abs(self._make().get_density() - 3.26) < 1e-10


# ---------------------------------------------------------------------------
# D5: FPSQ (Material)
# ---------------------------------------------------------------------------

class TestD5FPSQMaterial:
    """D5: FPSQ Material — tensor, eps_inf, correct DF type."""

    _eps_inf = np.diag([3.0, 3.0, 3.0])

    def _make(self):
        return FPSQ(
            "MgO",
            epsinf=self._eps_inf.tolist(),
            omega_tos=[[401.0], [401.0], [401.0]],
            gamma_tos=[[7.0],   [7.0],   [7.0]],
            omega_los=[[718.0], [718.0], [718.0]],
            gamma_los=[[7.0],   [7.0],   [7.0]],
            density=3.58,
        )

    def test_is_tensor_true(self):
        assert self._make().is_tensor() is True

    def test_permittivity_object_type(self):
        m = self._make()
        assert isinstance(m.get_permittivity_object(), DF.FPSQ)

    def test_epsilon_infinity_set(self):
        m = self._make()
        obj = m.get_permittivity_object()
        np.testing.assert_allclose(obj.epsilon_infinity, self._eps_inf, atol=1e-12)

    def test_type_string(self):
        assert self._make().type == "FPSQ"

    def test_density_set(self):
        assert abs(self._make().get_density() - 3.58) < 1e-10


# ---------------------------------------------------------------------------
# D6: Sellmeier (Material)
# ---------------------------------------------------------------------------

class TestD6SellmeierMaterial:
    """D6: Sellmeier Material — scalar, correct DF type, optical_permittivity."""

    def _make(self):
        return Sellmeier("SiO2", Bs=[0.696, 0.408, 0.897], Cs=[0.00468, 0.01351, 97.934], density=2.20)

    def test_is_scalar_true(self):
        assert self._make().is_scalar() is True

    def test_is_tensor_false(self):
        assert self._make().is_tensor() is False

    def test_permittivity_object_type(self):
        m = self._make()
        assert isinstance(m.get_permittivity_object(), DF.Sellmeier)

    def test_type_string(self):
        assert self._make().type == "Sellmeier"

    def test_optical_permittivity_reasonable(self):
        """SiO2 n ≈ 1.46 at visible, so n² ≈ 2.13."""
        m = self._make()
        op = m.get_optical_permittivity()
        assert 2.0 < op < 2.5

    def test_density_set(self):
        assert abs(self._make().get_density() - 2.20) < 1e-10


# ---------------------------------------------------------------------------
# D7: Tabulated (Material)
# ---------------------------------------------------------------------------

class TestD7TabulatedMaterial:
    """D7: Tabulated dispatches to the right DielectricFunction subclass."""

    _vs = np.linspace(100.0, 1000.0, 8)

    def test_1d_gives_tabulate_scalar(self):
        eps = np.linspace(2.0, 2.5, 8)
        m = Tabulated("mat1", vs_cm1=self._vs, permittivities=eps, density=1.0)
        assert isinstance(m.get_permittivity_object(), DF.TabulateScalar)
        assert m.is_scalar() is True

    def test_3row_gives_tabulate3(self):
        eps = np.array([
            np.linspace(2.0, 2.5, 8),
            np.linspace(2.5, 3.0, 8),
            np.linspace(3.0, 3.5, 8),
        ])
        m = Tabulated("mat3", vs_cm1=self._vs, permittivities=eps, density=1.0)
        assert isinstance(m.get_permittivity_object(), DF.Tabulate3)
        assert m.is_tensor() is True

    def test_6row_gives_tabulate6(self):
        eps = np.ones((6, 8)) * 2.0
        m = Tabulated("mat6", vs_cm1=self._vs, permittivities=eps, density=1.0)
        assert isinstance(m.get_permittivity_object(), DF.Tabulate6)
        assert m.is_tensor() is True

    def test_type_string(self):
        eps = np.linspace(2.0, 2.5, 8)
        m = Tabulated("mat1", vs_cm1=self._vs, permittivities=eps, density=1.0)
        assert m.type == "Tabulated permittivity"

    def test_get_information_contains_freq_range(self):
        eps = np.linspace(2.0, 2.5, 8)
        m = Tabulated("mat1", vs_cm1=self._vs, permittivities=eps, density=1.0)
        info = m.get_information()
        assert "freq range" in info.lower() or "Tabulate" in info


# ---------------------------------------------------------------------------
# D8: Material base helpers
# ---------------------------------------------------------------------------

class TestD8MaterialBaseHelpers:
    """D8: set_density, get_permittivity_function, get_information."""

    def test_set_get_density_roundtrip(self):
        m = Constant("test", permittivity=1.5, density=1.0)
        m.set_density(2.5)
        assert abs(m.get_density() - 2.5) < 1e-12

    def test_get_permittivity_function_callable(self):
        m = Constant("test", permittivity=2.0, density=1.0)
        fn = m.get_permittivity_function()
        assert callable(fn)

    def test_get_information_returns_string(self):
        m = Constant("test", permittivity=2.0, density=1.0)
        info = m.get_information()
        assert isinstance(info, str)
        assert len(info) > 0

    def test_get_permittivity_object_not_none(self):
        m = Constant("test", permittivity=2.0, density=1.0)
        assert m.get_permittivity_object() is not None

    def test_set_permittivity_object(self):
        m = Constant("test", permittivity=2.0, density=1.0)
        new_obj = DF.ConstantScalar(5.0)
        m.set_permittivity_object(new_obj)
        assert m.get_permittivity_object() is new_obj


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
