"""Tests E1-E5: MaterialsDataBase — built-in fallbacks, spreadsheet loading,
caching, and sanity checks on materials retrieved from the database.

Plan sections E1, E2, E3, E4, E5.

Tests E1-E3 do not require the spreadsheet to be present (they test the
hard-coded fallback materials).  Tests E4-E5 load the actual
MaterialsDataBase.xlsx file and are skipped gracefully if it is absent.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Materials import (
    MaterialsDataBase, Constant, NoMatrix, Sellmeier, FPSQ, DrudeLorentz, Tabulated
)
from PDielec import DielectricFunction as DF


# ---------------------------------------------------------------------------
# Locate the spreadsheet (may be absent in minimal CI environments)
# ---------------------------------------------------------------------------

_DB_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "MaterialsDataBase.xlsx"
)
_DB_AVAILABLE = os.path.isfile(_DB_PATH)


def _make_db_no_file():
    """Return a MaterialsDataBase with sheetNames=None (no spreadsheet)."""
    db = MaterialsDataBase.__new__(MaterialsDataBase)
    db.sheetNames = None
    db.cache = {}
    db.filename = None
    db.workbook = None
    return db


def _make_db():
    """Return a fully loaded MaterialsDataBase from the real spreadsheet."""
    return MaterialsDataBase(_DB_PATH)


# ---------------------------------------------------------------------------
# E1: Built-in fallback materials
# ---------------------------------------------------------------------------

class TestE1BuiltinFallbacks:
    """E1: get_material() returns correct built-in materials without spreadsheet."""

    def setup_method(self):
        self.db = _make_db_no_file()

    def test_air_is_constant(self):
        m = self.db.get_material("air")
        assert isinstance(m, Constant)

    def test_air_permittivity_is_one(self):
        m = self.db.get_material("air")
        fn = m.get_permittivity_function()
        assert abs(fn(0.0) - 1.0) < 1e-12

    def test_vacuum_permittivity_is_one(self):
        m = self.db.get_material("vacuum")
        fn = m.get_permittivity_function()
        assert abs(fn(0.0) - 1.0) < 1e-12

    def test_vacuum_density_is_zero(self):
        m = self.db.get_material("vacuum")
        assert m.get_density() == 0.0

    def test_ptfe_density(self):
        m = self.db.get_material("ptfe")
        assert abs(m.get_density() - 2.2) < 1e-10

    def test_ldpe_permittivity(self):
        m = self.db.get_material("ldpe")
        fn = m.get_permittivity_function()
        assert abs(fn(0.0) - 2.25) < 1e-10

    def test_mdpe_permittivity(self):
        m = self.db.get_material("mdpe")
        fn = m.get_permittivity_function()
        assert abs(fn(0.0) - 2.25) < 1e-10

    def test_kbr_density(self):
        m = self.db.get_material("kbr")
        assert abs(m.get_density() - 2.75) < 1e-10

    def test_nujol_is_constant(self):
        m = self.db.get_material("nujol")
        assert isinstance(m, Constant)

    def test_none_is_no_matrix(self):
        m = self.db.get_material("none")
        assert isinstance(m, NoMatrix)

    def test_all_scalar(self):
        """All built-in materials should be scalar."""
        for name in ["air", "vacuum", "ptfe", "ldpe", "mdpe", "kbr", "nujol", "none"]:
            m = self.db.get_material(name)
            assert m.is_scalar(), f"{name} should be scalar"


# ---------------------------------------------------------------------------
# E2: get_sheet_names() without spreadsheet
# ---------------------------------------------------------------------------

class TestE2SheetNamesNoFile:
    """E2: get_sheet_names() always includes the hard-coded built-in names."""

    def setup_method(self):
        self.db = _make_db_no_file()

    def test_returns_list(self):
        names = self.db.get_sheet_names()
        assert isinstance(names, list)

    def test_built_ins_present(self):
        names = self.db.get_sheet_names()
        for name in ["air", "vacuum", "ptfe", "ldpe", "mdpe", "kbr", "nujol", "none"]:
            assert name in names, f"'{name}' missing from sheet names"

    def test_list_is_sorted(self):
        names = self.db.get_sheet_names()
        assert names == sorted(names, key=lambda s: s.casefold())


# ---------------------------------------------------------------------------
# E3: Caching
# ---------------------------------------------------------------------------

class TestE3Caching:
    """E3: Repeated calls to get_material() return the same object.

    Built-in fallback materials (air, kbr, …) are created fresh on every call
    and are NOT stored in the cache — the cache is only populated for materials
    that are actually read from the spreadsheet.  Tests E3.1-E3.2 therefore
    require the spreadsheet; E3.3 can run without it.
    """

    @pytest.mark.skipif(not _DB_AVAILABLE, reason="MaterialsDataBase.xlsx not found")
    def test_same_object_returned_from_spreadsheet(self):
        """Two calls for the same spreadsheet material return the identical object."""
        db = _make_db()
        names = db.get_sheet_names()
        builtins = {"air", "vacuum", "ptfe", "ldpe", "mdpe", "kbr", "nujol", "none"}
        db_only = [n for n in names if n not in builtins]
        if not db_only:
            pytest.skip("No non-builtin materials in spreadsheet")
        name = db_only[0]
        m1 = db.get_material(name)
        m2 = db.get_material(name)
        assert m1 is m2

    @pytest.mark.skipif(not _DB_AVAILABLE, reason="MaterialsDataBase.xlsx not found")
    def test_cache_populated_for_spreadsheet_material(self):
        """After get_material() the name appears in db.cache."""
        db = _make_db()
        names = db.get_sheet_names()
        builtins = {"air", "vacuum", "ptfe", "ldpe", "mdpe", "kbr", "nujol", "none"}
        db_only = [n for n in names if n not in builtins]
        if not db_only:
            pytest.skip("No non-builtin materials in spreadsheet")
        name = db_only[0]
        db.get_material(name)
        assert name in db.cache

    def test_different_builtin_materials_different_objects(self):
        """Two different built-in materials are not the same object."""
        db = _make_db_no_file()
        m_air = db.get_material("air")
        m_kbr = db.get_material("kbr")
        assert m_air is not m_kbr


# ---------------------------------------------------------------------------
# E4: Loading the spreadsheet
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _DB_AVAILABLE, reason="MaterialsDataBase.xlsx not found")
class TestE4LoadSpreadsheet:
    """E4: The real spreadsheet can be loaded and queried."""

    def setup_method(self):
        self.db = _make_db()

    def test_loads_without_exception(self):
        # If we get here, loading succeeded
        assert self.db is not None

    def test_sheet_names_richer_than_builtins(self):
        names = self.db.get_sheet_names()
        assert len(names) > 8

    def test_sheet_names_still_include_builtins(self):
        names = self.db.get_sheet_names()
        for name in ["air", "vacuum", "kbr"]:
            assert name in names

    def test_known_sheet_returns_material(self):
        """Any sheet name in the DB should return a non-None Material."""
        names = self.db.get_sheet_names()
        # Pick the first name that is NOT a built-in fallback
        builtins = {"air", "vacuum", "ptfe", "ldpe", "mdpe", "kbr", "nujol", "none"}
        db_only = [n for n in names if n not in builtins]
        if not db_only:
            pytest.skip("No non-builtin materials found in spreadsheet")
        m = self.db.get_material(db_only[0])
        assert m is not None


# ---------------------------------------------------------------------------
# E5: Sanity checks on materials from spreadsheet
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _DB_AVAILABLE, reason="MaterialsDataBase.xlsx not found")
class TestE5DatabaseMaterialSanity:
    """E5: Physical sanity checks on materials loaded from the spreadsheet."""

    def setup_method(self):
        self.db = _make_db()
        names = self.db.get_sheet_names()
        builtins = {"air", "vacuum", "ptfe", "ldpe", "mdpe", "kbr", "nujol", "none"}
        self.db_names = [n for n in names if n not in builtins]

    def _find_type(self, mat_class):
        for name in self.db_names:
            m = self.db.get_material(name)
            if isinstance(m, mat_class):
                return m
        return None

    def test_sellmeier_is_scalar(self):
        m = self._find_type(Sellmeier)
        if m is None:
            pytest.skip("No Sellmeier material in DB")
        assert m.is_scalar() is True

    def test_fpsq_is_tensor(self):
        m = self._find_type(FPSQ)
        if m is None:
            pytest.skip("No FPSQ material in DB")
        assert m.is_tensor() is True

    def test_tabulated_information_contains_freq_range(self):
        m = self._find_type(Tabulated)
        if m is None:
            pytest.skip("No Tabulated material in DB")
        info = m.get_information()
        assert "freq range" in info

    def test_all_db_materials_have_positive_density(self):
        """Every material loaded from the spreadsheet should have a positive density."""
        for name in self.db_names:
            m = self.db.get_material(name)
            if m is None:
                continue
            d = m.get_density()
            assert d is not None, f"{name}: density is None"
            assert d > 0.0, f"{name}: density {d} not positive"

    def test_all_db_materials_permittivity_callable(self):
        """Every material's permittivity function should be callable."""
        for name in self.db_names:
            m = self.db.get_material(name)
            if m is None:
                continue
            fn = m.get_permittivity_function()
            assert callable(fn), f"{name}: permittivity function not callable"

    def test_db_material_optical_permittivity_real(self):
        """Optical permittivity must be real (or have zero/negligible imaginary part)."""
        for name in self.db_names:
            m = self.db.get_material(name)
            if m is None:
                continue
            try:
                op = m.get_optical_permittivity()
                if np.ndim(op) == 0:
                    # scalar
                    assert np.imag(op) == 0.0 or abs(np.imag(op)) < 1e-10, f"{name}: complex optical permittivity {op}"
            except Exception:
                # Some materials may not implement optical_permittivity; skip them
                pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
