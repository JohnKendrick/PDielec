"""Focused regression tests for the Irmer finite-q verification."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from analyse_verification import load_pdgui_spectra, summarise_pdgui
from polariton_response import (
    PAPER_LO_E1,
    PAPER_TO_E1,
    Scan,
    aperture_angles,
    paper_permittivity,
    raman_wavevector,
    response_point,
)


ROOT = Path(__file__).resolve().parent


def read_peaks():
    with (ROOT / "Data" / "polariton_peaks.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_pdgui_workbook_has_all_paper_aperture_scenarios():
    frequencies, spectra = load_pdgui_spectra(ROOT / "Data" / "irmer_pdgui.xlsx")
    assert len(spectra) == 26
    assert frequencies[0] == 150.0
    assert frequencies[-1] == 780.0
    assert np.allclose(np.diff(frequencies), 0.5)


def test_y_equals_two_mm_refracts_to_near_forward_angle():
    external, internal = aperture_angles(2.0, "y", "paper")
    assert np.isclose(external, 1.432096, atol=1.0e-5)
    assert 0.5 < internal < 0.7


def test_paper_dielectric_has_to_pole_and_lo_zero_in_lossless_limit():
    near_to = paper_permittivity(PAPER_TO_E1 + 1.0e-3, linewidth_cm1=0.0)[0, 0]
    at_lo = paper_permittivity(PAPER_LO_E1, linewidth_cm1=0.0)[0, 0]
    assert abs(near_to) > 1.0e5
    assert abs(at_lo) < 1.0e-10


def test_raman_q_is_frequency_dependent_and_finite():
    scan = Scan("9", 6.0, "xy", "z", "y")
    q_450 = raman_wavevector(scan, 450.0, "paper")
    q_500 = raman_wavevector(scan, 500.0, "paper")
    assert np.linalg.norm(q_450) > 0.0
    assert not np.allclose(q_450, q_500)


def test_response_peak_is_also_a_maxwell_pole():
    scan = Scan("9", 6.0, "xy", "z", "y")
    response_at_peak, locator_at_peak = response_point(scan, 477.25, "paper")
    response_off_peak, locator_off_peak = response_point(scan, 450.0, "paper")
    assert response_at_peak > response_off_peak
    assert locator_at_peak > locator_off_peak


def test_finite_q_peaks_move_monotonically_with_aperture():
    rows = read_peaks()
    for model in ("paper", "dft"):
        for figure, configuration in (
            ("9", "x(zy)x"),
            ("10", "x(yz)x"),
            ("11a", "x(zz)x"),
            ("11b", "x(yy)x"),
        ):
            selected = [
                row for row in rows
                if row["model"] == model and row["figure"] == figure
                and row["configuration"] == configuration
            ]
            selected.sort(key=lambda row: float(row["offset_mm"]))
            peaks = np.asarray([float(row["peak_frequency_cm-1"]) for row in selected])
            assert np.all(np.diff(peaks) > 0.0)


def test_directional_nac_control_does_not_follow_finite_q_branch():
    frequencies, spectra = load_pdgui_spectra(ROOT / "Data" / "irmer_pdgui.xlsx")
    rows = summarise_pdgui(frequencies, spectra)
    for figure in ("9", "10"):
        peaks = [row["nac_peak_frequency_cm-1"] for row in rows if row["figure"] == figure]
        assert max(peaks) - min(peaks) <= 0.5
        assert min(peaks) > 700.0


def test_dft_response_visits_485_region_in_figure_9_geometry():
    rows = read_peaks()
    y6 = next(
        row for row in rows
        if row["model"] == "dft" and row["figure"] == "9"
        and float(row["offset_mm"]) == 6.0
    )
    assert abs(float(y6["peak_frequency_cm-1"]) - 485.0) < 5.0
