"""Regression checks for command failures found during clean-install validation."""

import sys

import numpy as np
import pytest
from openpyxl import Workbook

from PDielec import pdcompare, pickled_reader


@pytest.mark.parametrize("command", [pdcompare, pickled_reader])
@pytest.mark.parametrize("flag", ["-h", "-help", "--help"])
def test_help_exits_without_reading_input(monkeypatch, capsys, command, flag):
    """Help flags print usage and exit successfully without treating flags as files."""
    name = command.__name__.rsplit(".", 1)[-1]
    monkeypatch.setattr(sys, "argv", [name, flag])
    with pytest.raises(SystemExit) as error:
        command.main()
    assert error.value.code in (None, 0)
    captured = capsys.readouterr()
    assert name in captured.out + captured.err


def test_pdcompare_without_excel_output(tmp_path, monkeypatch, capsys):
    """Compare real spectra without requiring the optional output workbook."""
    filename = tmp_path / "spectrum.xlsx"
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "Molar Absorption"
    worksheet.append(["", "Frequency", "", "Intensity"])
    for frequency, intensity in [(100, 1), (200, 3), (300, 2)]:
        worksheet.append(["", frequency, "", intensity])
    workbook.save(filename)
    workbook.close()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["pdcompare", str(filename), str(filename)])
    pdcompare.main()
    output = capsys.readouterr().out
    assert "correlations" in output
    correlations = output.rsplit("correlations", 1)[1].replace("[", "").replace("]", "")
    np.testing.assert_allclose(np.fromstring(correlations, sep=" "), np.ones(4), atol=1e-12)
    assert list(tmp_path.iterdir()) == [filename]
