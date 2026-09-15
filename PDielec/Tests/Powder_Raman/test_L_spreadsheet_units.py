"""Check Raman activity units in exported Settings worksheets."""

from types import SimpleNamespace

import numpy as np
import pytest
from openpyxl import load_workbook

from PDielec.GUI.SettingsTab import SettingsTab
from PDielec.GUI.SpreadSheetManager import SpreadSheetManager


@pytest.mark.parametrize("units,label", [("polarizability", "Å⁴/amu"), ("epsilon", "Å/amu")])
@pytest.mark.parametrize("volume_factor", [2.0, 3.0])
def test_exported_raman_activities_use_settings_table_units(tmp_path, units, label, volume_factor):
    """Saved numbers and headers follow the selected units without changing stored activities."""
    filename = tmp_path / "settings.xlsx"
    spreadsheet = SpreadSheetManager(str(filename))
    tab = SettingsTab.__new__(SettingsTab)
    tab.notebook = SimpleNamespace(spreadsheet=spreadsheet)
    tab.reader = SimpleNamespace(volume=volume_factor * 16.0 * np.pi**2)
    tab.settings = {"Raman activity units": units}
    tab.frequencies_cm1 = [100.0, 200.0]
    tab.intensities = [1.0, 2.0]
    tab.sigmas_cm1 = [5.0, 6.0]
    tab.modes_selected = [True, False]
    tab.raman_intensities = [10.0, 20.0]
    tab.raman_intensities_par = [7.0, 14.0]
    tab.raman_intensities_perp = [3.0, 6.0]
    original = np.array([
        tab.raman_intensities, tab.raman_intensities_par, tab.raman_intensities_perp,
    ])

    try:
        tab.write_spreadsheet()
    finally:
        spreadsheet.close()

    workbook = load_workbook(filename, read_only=True, data_only=True)
    try:
        rows = list(workbook["Settings"].values)
        header_index = next(i for i, row in enumerate(rows) if row[1] == "Mode")
        assert rows[header_index][8:11] == tuple(
            f"Raman R_{component} ({label})" for component in ("total", "parallel", "perp")
        )
        factor = volume_factor if units == "polarizability" else 1.0
        for row, expected in zip(rows[header_index + 1:], original.T):
            np.testing.assert_allclose(row[8:11], factor * expected)
        assert len(rows[header_index + 1:]) == 2
        assert rows[header_index + 1][2] == "Yes"
        assert rows[header_index + 2][2] == "No"
        assert rows[header_index + 1][4:6] == (100.0, 1.0)
    finally:
        workbook.close()

    np.testing.assert_array_equal(
        [tab.raman_intensities, tab.raman_intensities_par, tab.raman_intensities_perp], original,
    )
