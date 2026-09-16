"""Regression checks for the two incoherent modal-pairs models."""

import numpy as np

from analyse_incoherent_models import load_spectra, summarise


def scenario(angle, model, tolerance, channel):
    """Return the workbook legend for one scenario."""
    return (
        f"theta_ext={angle:g} model={model} tol={tolerance:g} "
        f"channel={channel}"
    )


def test_complete_scenario_matrix():
    frequencies, spectra = load_spectra()
    assert len(spectra) == 32
    assert frequencies[0] == 400.0
    assert frequencies[-1] == 800.0
    assert np.allclose(np.diff(frequencies), 0.5)


def test_local_incoherent_is_tolerance_invariant():
    _, spectra = load_spectra()
    for angle in (0.0, 2.0):
        for channel in ("e-e", "o-o", "e-o", "o-e"):
            tol_0 = spectra[scenario(angle, "local", 0.0, channel)]
            tol_90 = spectra[scenario(angle, "local", 90.0, channel)]
            assert np.array_equal(tol_0, tol_90)


def test_bulk_tolerance_selects_internal_pairs_at_two_degrees():
    frequencies, spectra = load_spectra()
    rows = summarise(frequencies, spectra)
    index = {
        (row["model"], row["q_tolerance_deg"], row["channel"]): row
        for row in rows
        if row["angle_external_deg"] == 2.0
    }
    for channel in ("e-e", "e-o", "o-e"):
        tol_0 = index[("bulk", 0.0, channel)]["integrated_intensity_450_780"]
        tol_90 = index[("bulk", 90.0, channel)]["integrated_intensity_450_780"]
        assert tol_0 < 1.0e-8 * tol_90
    assert (
        index[("bulk", 0.0, "o-o")]["integrated_intensity_450_780"]
        < 0.6 * index[("bulk", 90.0, "o-o")]["integrated_intensity_450_780"]
    )


def test_models_change_intensity_not_accepted_line_frequency():
    frequencies, spectra = load_spectra()
    rows = summarise(frequencies, spectra)
    index = {
        (row["model"], row["channel"]): row
        for row in rows
        if row["angle_external_deg"] == 2.0 and row["q_tolerance_deg"] == 90.0
    }
    for channel in ("e-e", "o-o", "e-o", "o-e"):
        assert (
            index[("bulk", channel)]["peak_frequency_cm-1"]
            == index[("local", channel)]["peak_frequency_cm-1"]
        )
        assert (
            index[("bulk", channel)]["integrated_intensity_450_780"]
            < index[("local", channel)]["integrated_intensity_450_780"]
        )

