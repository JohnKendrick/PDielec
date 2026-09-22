"""Tests for CASTEP atom-resolved DFPT Raman derivatives."""

import io
import math
import os

import numpy as np
import pytest

from PDielec import Calculator
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.Constants import angs2bohr
from PDielec.GUI.SettingsTab import SettingsTab
from PDielec.UnitCell import UnitCell

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
_CASTEP_25_FILE = os.path.join(
    _REPO,
    "Examples",
    "Crystal_Raman",
    "Castep-25.12",
    "raman.castep",
)


def _new_reader(tmp_path):
    """Return a fully initialised reader for an empty temporary CASTEP seed."""
    castep_file = tmp_path / "synthetic.castep"
    castep_file.write_text("")
    return CastepOutputReader([str(castep_file)])


def _reader_with_polar_deriv_block(tmp_path, nions=1):
    """Return a reader after parsing a complete synthetic POLAR_DERIV block."""
    reader = _new_reader(tmp_path)
    reader.nions = nions
    rows = []
    value = 1.0
    for atom in range(1, nions + 1):
        for displacement in range(1, 4):
            for tensor_row in range(1, 4):
                rows.append(
                    f"{atom:6d}{displacement:6d}{tensor_row:6d}"
                    f" {value:25.15E} {value + 1.0:25.15E} {value + 2.0:25.15E}"
                )
                value += 3.0
    block = (
        " polar_deriv(atom,displacement,tensor_row,tensor_column) in internal atomic units\n"
        " atom displacement tensor_row tensor_column_1 tensor_column_2 tensor_column_3\n"
        + "\n".join(rows)
        + "\n END DFPT RAMAN POLAR_DERIV\n"
    )
    reader.file_descriptor = io.StringIO(block)
    reader._read_raman_polar_derivatives(" BEGIN DFPT RAMAN POLAR_DERIV")
    return reader


def test_castep_atomic_raman_derivative_block_is_read_in_cartesian_layout(tmp_path):
    """The machine-readable block should become [atom, displacement, row, column]."""
    reader = _reader_with_polar_deriv_block(tmp_path, nions=2)

    derivatives = reader._castep_raman_polar_deriv
    assert derivatives.shape == (2, 3, 3, 3)
    np.testing.assert_allclose(derivatives[0, 0, 0], [1.0, 2.0, 3.0])
    np.testing.assert_allclose(derivatives[0, 2, 2], [25.0, 26.0, 27.0])
    np.testing.assert_allclose(derivatives[1, 0, 0], [28.0, 29.0, 30.0])
    np.testing.assert_allclose(derivatives[1, 2, 2], [52.0, 53.0, 54.0])


def test_castep_atomic_derivatives_are_projected_with_the_phonon_mode_phase(tmp_path):
    """Reversing a mode must reverse its reconstructed Raman tensor."""
    reader = _new_reader(tmp_path)
    reader.nions = 1
    reader.masses = [4.0]
    reader.volume = 4.0
    derivatives = np.zeros((1, 3, 3, 3))
    # One Angstrom^2 in CASTEP's internal Bohr^2 representation.
    derivatives[0, 0, 0, 0] = angs2bohr**2
    reader._castep_raman_polar_deriv = derivatives

    reader.mass_weighted_normal_modes = np.array([[[1.0, 0.0, 0.0]]])
    reader._calculate_raman_tensors_from_polar_derivatives()
    positive_tensor = np.asarray(reader.raman_tensors[0])
    np.testing.assert_allclose(positive_tensor[0, 0], math.pi)

    reader.mass_weighted_normal_modes = np.array([[[-1.0, 0.0, 0.0]]])
    reader._update_raman_tensors_after_mode_recalculation(None, None)
    negative_tensor = np.asarray(reader.raman_tensors[0])
    np.testing.assert_allclose(negative_tensor, -positive_tensor)


def test_castep_atomic_derivatives_replace_preprojected_fallback_tensors(tmp_path):
    """Raw derivatives should take precedence over CASTEP's projected mode block."""
    reader = _new_reader(tmp_path)
    reader.nions = 1
    reader.masses = [1.0]
    reader.volume = 1.0
    reader.mass_weighted_normal_modes = np.array([[[1.0, 0.0, 0.0]]])
    reader.raman_tensors = [99.0 * np.eye(3)]
    reader._castep_raman_polar_deriv = np.zeros((1, 3, 3, 3))
    reader._castep_raman_polar_deriv[0, 0] = angs2bohr**2 * np.eye(3)

    reader._calculate_raman_tensors_from_polar_derivatives()

    np.testing.assert_allclose(reader.raman_tensors[0], 4.0 * math.pi * np.eye(3))


def test_official_atomic_polar_tensor_table_is_preferred_and_mapped_by_species(tmp_path):
    """The official Angstrom-squared table should map species-local atom indices."""
    reader = _new_reader(tmp_path)
    reader.nions = 2
    cell = UnitCell([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
    cell.set_element_names(["Zn", "O"])
    reader.unit_cells = [cell]

    rows = [" +----------------------------------------------------------------------------+"]
    for species, base in (("O", 10.0), ("Zn", 20.0)):
        for displacement, offset in zip("XYZ", (0.0, 1.0, 2.0)):
            value = base + offset
            rows.extend(
                [
                    f" + {species:2s}      1   {displacement}  {value:12.6f} 0.000000 0.000000 +",
                    " +                              0.000000 0.000000 0.000000 +",
                    " +                              0.000000 0.000000 0.000000 +",
                    " +                                                                            +",
                ]
            )
    rows.append(" +============================================================================+")
    reader.file_descriptor = io.StringIO("\n".join(rows) + "\n")
    reader._read_raman_atomic_polar_tensors(
        " + Raman atomic polar tensor V*dChi(1)/dR (A**2) +"
    )

    assert reader._castep_raman_polar_deriv_source == "official"
    derivatives_a2 = reader._castep_raman_polar_deriv / angs2bohr**2
    np.testing.assert_allclose(derivatives_a2[0, :, 0, 0], [20.0, 21.0, 22.0])
    np.testing.assert_allclose(derivatives_a2[1, :, 0, 0], [10.0, 11.0, 12.0])

    official_derivatives = reader._castep_raman_polar_deriv.copy()
    reader._set_raman_polar_derivatives(np.full_like(official_derivatives, 99.0), "developer")
    np.testing.assert_allclose(reader._castep_raman_polar_deriv, official_derivatives)
    assert reader._castep_raman_polar_deriv_source == "official"


@pytest.mark.skipif(not os.path.exists(_CASTEP_25_FILE), reason="CASTEP 25.12 example not present")
def test_castep_25_official_atomic_polar_tensors_are_read_and_translationally_invariant():
    """The official CASTEP 25.12 table should satisfy its imposed Raman sum rule."""
    reader = CastepOutputReader([_CASTEP_25_FILE])
    reader.read_output()

    assert reader._castep_raman_polar_deriv_source == "official"
    assert reader._castep_raman_polar_deriv.shape == (4, 3, 3, 3)
    translation_residual_a2 = np.sum(reader._castep_raman_polar_deriv, axis=0) / angs2bohr**2
    np.testing.assert_allclose(translation_residual_a2, 0.0, atol=1.1e-6)

    # Reprojection from the official atom-resolved table recovers the CASTEP
    # mode-7 susceptibility tensor to its six-decimal printed precision.
    expected_mode_7 = 4.0 * math.pi * np.diag([0.051242, 0.051259, 0.110664])
    np.testing.assert_allclose(np.diag(reader.raman_tensors[6]), np.diag(expected_mode_7), atol=7.0e-6)


@pytest.mark.skipif(not os.path.exists(_CASTEP_25_FILE), reason="CASTEP 25.12 example not present")
def test_castep_25_raman_tensors_follow_the_phonon_mode_order():
    """The official tensors should remain aligned with CASTEP's ordered phonon modes."""
    reader = CastepOutputReader([_CASTEP_25_FILE])
    reader.read_output()

    assert len(reader.frequencies) == 12
    assert len(reader.raman_tensors) == len(reader.frequencies)
    assert reader.frequencies[5:8] == pytest.approx(
        [259.437429, 310.551568, 340.074893],
        abs=1.0e-6,
    )
    expected_mode_7 = 4.0 * math.pi * np.array(
        [
            [0.051242, -0.000003, -0.000079],
            [-0.000003, 0.051259, -0.000019],
            [-0.000079, -0.000019, 0.110664],
        ]
    )
    np.testing.assert_allclose(reader.raman_tensors[6], expected_mode_7, atol=7.0e-6)


@pytest.mark.skipif(not os.path.exists(_CASTEP_25_FILE), reason="CASTEP 25.12 example not present")
def test_castep_25_mode_7_activity_uses_pdgui_display_conversion():
    """PDGui should display CASTEP's mode-7 activity in polarizability-volume units."""
    reader = CastepOutputReader([_CASTEP_25_FILE])
    reader.read_output()
    activities = Calculator.raman_intensities(reader.get_raman_tensors(), reader.volume)

    settings_tab = SettingsTab.__new__(SettingsTab)
    settings_tab.reader = reader
    settings_tab.settings = {"Raman activity units": "polarizability"}
    display = settings_tab._raman_activity_display_metadata()

    assert display["unit_label"] == "Å⁴/amu"
    assert display["factor"] == pytest.approx(reader.volume / (16.0 * math.pi * math.pi))
    assert activities[6, 0] * display["factor"] == pytest.approx(12.7504798, abs=2.0e-6)

    settings_tab.settings["Raman activity units"] = "epsilon"
    epsilon_display = settings_tab._raman_activity_display_metadata()
    assert epsilon_display["unit_label"] == "Å/amu"
    assert epsilon_display["factor"] == 1.0
