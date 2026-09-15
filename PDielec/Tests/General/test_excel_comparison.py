"""Regression coverage for Settings comparisons and missing result sheets."""

import sys

import pytest
from openpyxl import Workbook

from PDielec import checkexcel, pdmake


def _workbook(path, sheets):
    """Write a small workbook with explicitly supplied sheet data."""
    workbook = Workbook()
    workbook.remove(workbook.active)
    for name, rows in sheets.items():
        sheet = workbook.create_sheet(name)
        for row in rows:
            sheet.append(row)
    workbook.save(path)
    workbook.close()


def _compare(tmp_path, monkeypatch, left, right, *flags):
    """Compare temporary workbooks using the command-line interface."""
    reference = tmp_path / "reference.xlsx"
    result = tmp_path / "result.xlsx"
    _workbook(reference, left)
    _workbook(result, right)
    monkeypatch.setattr(sys, "argv", ["checkexcel", str(reference), str(result), *flags])
    return checkexcel.main()


@pytest.mark.parametrize("changed", [[["Raman (Å⁴/amu)", 10.0]], [["Raman (Å/amu)", 1.0]]])
def test_settings_values_and_headers_are_checked(tmp_path, monkeypatch, changed):
    """Settings-only comparison detects both numeric and unit-label differences."""
    reference = {"Settings": [["Raman (Å⁴/amu)", 1.0]], "Main": [["old path"]]}
    result = {"Settings": changed, "Main": [["new path"]]}
    assert _compare(tmp_path, monkeypatch, reference, result)[0] == 0
    assert _compare(tmp_path, monkeypatch, reference, result, "-settings")[0] == 1
    assert _compare(tmp_path, monkeypatch, reference, result, "-f")[0] == 2


@pytest.mark.parametrize("sheet", ["Powder Raman", "Crystal Raman", "Settings"])
@pytest.mark.parametrize("reverse", [False, True])
def test_missing_checked_sheet_is_an_error(tmp_path, monkeypatch, capsys, sheet, reverse):
    """A sheet present in only one workbook fails, in either comparison direction."""
    left = {"Main": [["same"]], sheet: [[1.0]]}
    right = {"Main": [["same"]]}
    if reverse:
        left, right = right, left
    result = _compare(tmp_path, monkeypatch, left, right, "-settings")
    assert result[0] == 1
    assert result[3] == sheet
    output = capsys.readouterr().out
    assert f"missing expected sheet '{sheet}'" in output
    assert ("reference.xlsx" if reverse else "result.xlsx") in output


def test_optional_sheets_absent_from_both_are_not_errors(tmp_path, monkeypatch):
    """Infrared workbooks do not require Raman sheets, or vice versa."""
    data = {"Settings": [["Units", "Å⁴/amu"]], "Powder Raman": [[1.0]]}
    assert _compare(tmp_path, monkeypatch, data, data, "-settings")[0] == 0


def test_missing_sheet_does_not_hide_spectrum_errors(tmp_path, monkeypatch, capsys):
    """Continue numeric comparison after finding a missing sheet."""
    left = {"Powder Raman": [[1.0]], "Crystal Raman": [[1.0]]}
    right = {"Crystal Raman": [[10.0]]}
    assert _compare(tmp_path, monkeypatch, left, right)[0] == 2
    assert "missing expected sheet 'Powder Raman'" in capsys.readouterr().out


@pytest.mark.parametrize("settings_value,omit_raman", [(10.0, False), (1.0, True)])
def test_pdmake_checks_settings_and_reports_missing_sheets(
    tmp_path, monkeypatch, capsys, settings_value, omit_raman,
):
    """Exercise pdmake's real comparator without launching the GUI or updating references."""
    monkeypatch.chdir(tmp_path)
    reference = tmp_path / "results.ref.xlsx"
    _workbook(reference, {"Settings": [[1.0]], "Crystal Raman": [[1.0]]})
    original_reference = reference.read_bytes()
    monkeypatch.setattr(pdmake, "useLocal", False)
    monkeypatch.setattr(pdmake, "viewing", False)
    monkeypatch.setattr(pdmake, "debug", False)
    monkeypatch.setitem(pdmake.settings, "PDGui tokens", [])
    monkeypatch.setattr(sys, "argv", ["pdmake"])

    def generate_result(*args, **kwargs):
        sheets = {"Settings": [[settings_value]]}
        if not omit_raman:
            sheets["Crystal Raman"] = [[1.0]]
        _workbook(tmp_path / "results.xlsx", sheets)

    monkeypatch.setattr(pdmake.subprocess, "run", generate_result)
    pdmake.run_pdgui_test("Example", [], regenerate=False)
    output = capsys.readouterr().out
    assert "1 ERRORS:" in output
    if omit_raman:
        assert "missing expected sheet 'Crystal Raman'" in output
        assert "results.xlsx" in output
    else:
        assert "OF SHEET Settings" in output
    assert reference.read_bytes() == original_reference
