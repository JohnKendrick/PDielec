"""Finite-field Raman JSON reader checks."""

import json
import os

import numpy as np
import pytest

from PDielec import Calculator
from PDielec.Constants import angs2bohr, average_masses, hartree2ev
from PDielec.FiniteFieldOutputReader import FiniteFieldOutputReader
from PDielec.Utilities import find_program_from_name, get_reader

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
_RAMAN_JSON = os.path.join(
    _REPO,
    "Examples",
    "Powder_Raman",
    "Finite_difference",
    "finite_difference.json",
)


@pytest.mark.skipif(not os.path.exists(_RAMAN_JSON), reason="ZnO finite-difference JSON example not present")
def test_finite_field_reader_detects_json_schema():
    """Finite-field JSON files are detected as finite_field datasets."""
    assert find_program_from_name(_RAMAN_JSON) == "finite_field"
    assert isinstance(get_reader(_RAMAN_JSON, "finite_field"), FiniteFieldOutputReader)


@pytest.mark.skipif(not os.path.exists(_RAMAN_JSON), reason="ZnO finite-difference JSON example not present")
def test_finite_field_reader_loads_core_tensors_and_average_masses():
    """The reader loads structure/tensors and defaults to PDielec average masses."""
    reader = FiniteFieldOutputReader([_RAMAN_JSON])
    reader.read_output()

    assert reader.type == "Finite field Raman JSON"
    assert reader.nions == 4
    assert reader.volume == pytest.approx(45.62538772357219)
    assert np.asarray(reader.zerof_optical_dielectric).shape == (3, 3)
    assert np.asarray(reader.born_charges).shape == (4, 3, 3)
    assert np.asarray(reader.nonlinear_optical_susceptibility).shape == (3, 3, 3)
    assert len(reader.frequencies) == 12
    assert len(reader.mass_weighted_normal_modes) == 12
    assert len(reader.get_raman_tensors()) == 12
    assert reader.masses_per_type == pytest.approx([average_masses["Zn"], average_masses["O"]])
    assert reader.program_mass_dictionary["Zn"] == pytest.approx(65.38)
    assert reader.program_mass_dictionary["O"] == pytest.approx(15.999)


@pytest.mark.skipif(not os.path.exists(_RAMAN_JSON), reason="ZnO finite-difference JSON example not present")
def test_finite_field_reader_recalculates_modes_from_force_constants():
    """The active modes are regenerated from force constants and remain Raman-aligned."""
    reader = FiniteFieldOutputReader([_RAMAN_JSON])
    reader.read_output()

    assert reader.nomass_hessian_has_been_set
    assert np.asarray(reader.hessian).shape == (12, 12)
    assert len(reader.get_raman_tensors()) == len(reader.frequencies)
    assert np.all(np.diff(reader.frequencies) >= -1.0e-8)

    activities = Calculator.raman_intensities(reader.get_raman_tensors(), reader.volume)
    assert activities.shape[0] == 12
    assert np.isfinite(activities).all()


@pytest.mark.skipif(not os.path.exists(_RAMAN_JSON), reason="ZnO finite-difference JSON example not present")
def test_finite_field_reader_derives_raw_displacement_dataset_tensors(tmp_path):
    """The lower-level finite-displacement dataset can be read directly."""
    with open(_RAMAN_JSON) as fd:
        data = json.load(fd)
    # Keep only raw input records so precomputed modes, derivatives and force
    # constants cannot bypass the finite-difference reconstruction paths.
    raw_data = {
        "schema": data["schema"],
        "structure": data["structure"],
        "reference": data["reference"],
        "cartesian_displacements": data["cartesian_displacements"],
    }
    raw_file = tmp_path / "finite_field_dataset.json"
    raw_file.write_text(json.dumps(raw_data))
    reader = FiniteFieldOutputReader([str(raw_file)])
    reader.read_output()

    assert reader.nomass_hessian_has_been_set
    assert reader._deps_dr is not None
    assert len(reader.frequencies) == 12
    assert len(reader.get_raman_tensors()) == 12
    assert np.asarray(reader.born_charges).shape == (4, 3, 3)

    # Compare with the saved raw derivatives, before the production dataset's
    # sum-rule corrections to the top-level deps_dr tensor.
    for item in data["cartesian_displacements"]:
        np.testing.assert_allclose(reader._deps_dr[item["atom"], item["cart"]], item["deps_dr"], atol=1.0e-12)
    np.testing.assert_allclose(reader.born_charges, data["reference"]["born_effective_charges"])
    # This raw grid uses 0.010 Angstrom displacements; the production force
    # constants use a separate 0.015 Angstrom grid, so compare with the saved
    # force-constant columns belonging to these actual raw records.
    force_constants = np.column_stack([
        item["force_constant_column_eV_per_angstrom2"]
        for item in data["cartesian_displacements"]
    ])
    expected_hessian = (force_constants + force_constants.T) / (2.0 * hartree2ev * angs2bohr**2)
    np.testing.assert_allclose(reader.nomass_hessian, expected_hessian, atol=1.0e-12)
    assert np.isfinite(reader.get_raman_tensors()).all()
