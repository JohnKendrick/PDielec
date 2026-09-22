"""Regression checks for the generated GaN modal-pairs workbook."""

from pathlib import Path

import numpy as np

from analyse_results import load_spectra, summarise


WORKBOOK = Path(__file__).with_name("Data") / "gan_modal_pairs.xlsx"


def _indexed_summary():
    frequencies, spectra = load_spectra(WORKBOOK)
    rows = summarise(frequencies, spectra)
    index = {
        (row["angle_external_deg"], row["nac_mode"], row["channel"]): row
        for row in rows
    }
    return frequencies, spectra, rows, index


def test_workbook_has_complete_scenario_matrix():
    frequencies, spectra, rows, _ = _indexed_summary()
    assert len(spectra) == 24
    assert len(rows) == 24
    assert np.allclose(np.diff(frequencies), 0.5)
    assert frequencies[0] == 400.0
    assert frequencies[-1] == 800.0


def test_single_q_methods_agree_in_this_homogeneous_geometry():
    _, spectra, _, _ = _indexed_summary()
    for angle in (0, 2):
        for channel in ("e-e", "o-o", "e-o", "o-e"):
            geometry = spectra[f"theta_ext={angle:g} geometry {channel}"]
            dominant = spectra[f"theta_ext={angle:g} dominant {channel}"]
            peak_scale = np.max(np.abs(geometry))
            assert np.max(np.abs(geometry - dominant)) < 1.0e-5 * peak_scale


def test_modal_pairs_resolves_perfect_forward_crossed_channels():
    _, _, _, index = _indexed_summary()
    for channel in ("e-o", "o-e"):
        geometry = index[(0.0, "geometry", channel)]
        modal = index[(0.0, "modal_pairs", channel)]
        # Pair-resolved amplitudes nearly cancel in perfect-forward crossed
        # polarization; the single-q approximation cannot represent this.
        assert modal["integrated_intensity_450_780"] < (
            1.0e-8 * geometry["integrated_intensity_450_780"]
        )


def test_modal_pairs_detects_crossed_channel_swap_at_small_angle():
    _, _, _, index = _indexed_summary()
    e_to_o = index[(2.0, "modal_pairs", "e-o")]["integrated_intensity_450_780"]
    o_to_e = index[(2.0, "modal_pairs", "o-e")]["integrated_intensity_450_780"]
    assert e_to_o > 1.4 * o_to_e
