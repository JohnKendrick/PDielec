"""Finite-field Raman JSON reader checks."""

import os

import numpy as np
import pytest

from PDielec import Calculator
from PDielec.Constants import average_masses
from PDielec.FiniteFieldOutputReader import FiniteFieldOutputReader
from PDielec.Utilities import find_program_from_name, get_reader

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
_RAMAN_JSON = os.path.join(
    _REPO,
    "Examples",
    "Crystal_Raman",
    "Finite_field",
    "raman_results.json",
)
_FINITE_FIELD_JSON = os.path.join(
    _REPO,
    "Examples",
    "Crystal_Raman",
    "Finite_field",
    "finite_field_dataset.json",
)


@pytest.mark.skipif(not os.path.exists(_RAMAN_JSON), reason="BaTiO3 finite-field Raman JSON example not present")
def test_finite_field_reader_detects_json_schema():
    """Finite-field JSON files are detected as finite_field datasets."""
    assert find_program_from_name(_RAMAN_JSON) == "finite_field"
    assert isinstance(get_reader(_RAMAN_JSON, "finite_field"), FiniteFieldOutputReader)


@pytest.mark.skipif(not os.path.exists(_RAMAN_JSON), reason="BaTiO3 finite-field Raman JSON example not present")
def test_finite_field_reader_loads_core_tensors_and_average_masses():
    """The reader loads structure/tensors and defaults to PDielec average masses."""
    reader = FiniteFieldOutputReader([_RAMAN_JSON])
    reader.read_output()

    assert reader.type == "Finite field Raman JSON"
    assert reader.nions == 5
    assert reader.volume == pytest.approx(64.953592751027)
    assert np.asarray(reader.zerof_optical_dielectric).shape == (3, 3)
    assert np.asarray(reader.born_charges).shape == (5, 3, 3)
    assert np.asarray(reader.nonlinear_optical_susceptibility).shape == (3, 3, 3)
    assert len(reader.frequencies) == 15
    assert len(reader.mass_weighted_normal_modes) == 15
    assert len(reader.get_raman_tensors()) == 15
    assert reader.masses_per_type == pytest.approx([average_masses["Ba"], average_masses["O"], average_masses["Ti"]])
    assert reader.program_mass_dictionary["Ba"] == pytest.approx(average_masses["Ba"])


@pytest.mark.skipif(not os.path.exists(_RAMAN_JSON), reason="BaTiO3 finite-field Raman JSON example not present")
def test_finite_field_reader_recalculates_modes_from_force_constants():
    """The active modes are regenerated from force constants and remain Raman-aligned."""
    reader = FiniteFieldOutputReader([_RAMAN_JSON])
    reader.read_output()

    assert reader.nomass_hessian_has_been_set
    assert np.asarray(reader.hessian).shape == (15, 15)
    assert len(reader.get_raman_tensors()) == len(reader.frequencies)
    assert np.all(np.diff(reader.frequencies) >= -1.0e-8)

    activities = Calculator.raman_intensities(reader.get_raman_tensors(), reader.volume)
    assert activities.shape[0] == 15
    assert np.isfinite(activities).all()


@pytest.mark.skipif(not os.path.exists(_FINITE_FIELD_JSON), reason="BaTiO3 lower-level finite-field JSON not present")
def test_finite_field_reader_derives_raw_displacement_dataset_tensors():
    """The lower-level finite-displacement dataset can be read directly."""
    reader = FiniteFieldOutputReader([_FINITE_FIELD_JSON])
    reader.read_output()

    assert reader.nomass_hessian_has_been_set
    assert reader._deps_dr is not None
    assert len(reader.frequencies) == 15
    assert len(reader.get_raman_tensors()) == 15
    assert np.asarray(reader.born_charges).shape == (5, 3, 3)
