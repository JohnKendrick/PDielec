"""Crystal Raman electro-optic toggle tests."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np

import PDielec.GUI.CrystalScenarioTab as crystal_tab
from PDielec.GUI.CrystalScenarioTab import CrystalScenarioTab, _selection_for_nac_modes
from PDielec.LayeredRamanCalculator import MODAL_PAIR_GROUP_Q, LayeredRamanCalculator, RamanLayer

from .conftest import build_system


class _Reader:
    nions = 1
    volume = 100.0
    masses = [1.0]
    born_charges = [np.eye(3)]
    hessian = np.eye(3)
    mass_weighted_normal_modes = [
        [[1.0, 0.0, 0.0]],
        [[0.0, 1.0, 0.0]],
        [[0.0, 0.0, 1.0]],
    ]
    nonlinear_optical_susceptibility = np.ones((3, 3, 3))


def _scenario(include_eo):
    scenario = CrystalScenarioTab.__new__(CrystalScenarioTab)
    scenario.reader = _Reader()
    scenario.settings = {
        "Raman electro-optic term": include_eo,
    }
    return scenario


def test_crystal_eo_setting_disabled_passes_no_chi2(monkeypatch):
    """Crystal Raman disables the separable χ² tensor before NAC calculation."""
    captured = {}

    def fake_standalone(*args, **kwargs):
        captured["chi2"] = kwargs["chi2_repsilon"]
        return np.array([1.0, 2.0, 3.0]), [np.eye(3)] * 3, np.ones(3)

    monkeypatch.setattr(crystal_tab, "_compute_nac_dynamical_matrix_standalone", fake_standalone)

    scenario = _scenario(include_eo=False)
    scenario._compute_nac_dynamical_matrix(
        raman_tensors_physical=[np.eye(3)] * 3,
        frequencies_cm1=np.array([100.0, 200.0, 300.0]),
        sigmas_cm1=np.ones(3),
        q_hat_crystal=np.array([0.0, 0.0, 1.0]),
        eps_inf=np.eye(3),
    )

    assert captured["chi2"] is None


def test_crystal_eo_setting_enabled_passes_reader_chi2(monkeypatch):
    """Crystal Raman includes reader χ² when the setting is enabled."""
    captured = {}

    def fake_standalone(*args, **kwargs):
        captured["chi2"] = kwargs["chi2_repsilon"]
        return np.array([1.0, 2.0, 3.0]), [np.eye(3)] * 3, np.ones(3)

    monkeypatch.setattr(crystal_tab, "_compute_nac_dynamical_matrix_standalone", fake_standalone)

    scenario = _scenario(include_eo=True)
    scenario._compute_nac_dynamical_matrix(
        raman_tensors_physical=[np.eye(3)] * 3,
        frequencies_cm1=np.array([100.0, 200.0, 300.0]),
        sigmas_cm1=np.ones(3),
        q_hat_crystal=np.array([0.0, 0.0, 1.0]),
        eps_inf=np.eye(3),
    )

    assert captured["chi2"] is scenario.reader.nonlinear_optical_susceptibility


def test_crystal_raman_respects_mode_selection_mask():
    """Crystal Raman should skip Settings-tab modes that are deselected."""
    frequencies_cm1 = np.array([100.0, 200.0])
    raman_tensors = [np.eye(3), 2.0 * np.eye(3)]
    system = build_system([(100e-9, 1.5)])
    raman_layer = RamanLayer(
        layer_index=0,
        phonon_frequencies_cm1=frequencies_cm1,
        raman_tensors=raman_tensors,
        rotation_matrix=np.eye(3),
    )

    calc = LayeredRamanCalculator(
        system=system,
        raman_layers=[raman_layer],
        laser_frequency_cm1=20000.0,
        incident_angle_rad=0.0,
        incident_pol="p",
        detected_pol="p",
        temperature_K=0.0,
        linewidths_cm1=np.array([5.0, 5.0]),
        n_gauss=10,
        approximate_es=True,
        modes_selected=[True, False],
    )

    active_freqs, active_intensities, active_sigmas = calc.calculate_mode_intensities()

    np.testing.assert_allclose(active_freqs, [100.0])
    assert len(active_intensities) == 1
    np.testing.assert_allclose(active_sigmas, [5.0])


def test_crystal_raman_modal_pairs_respects_mode_selection_mask():
    """The q-channel modal-pair path should use the same selected-mode mask."""
    frequencies_cm1 = np.array([100.0, 200.0])
    raman_tensors = [np.eye(3), 2.0 * np.eye(3)]
    system = build_system([(100e-9, 1.5)])
    raman_layer = RamanLayer(
        layer_index=0,
        phonon_frequencies_cm1=frequencies_cm1,
        raman_tensors=raman_tensors,
        rotation_matrix=np.eye(3),
    )

    calc = LayeredRamanCalculator(
        system=system,
        raman_layers=[raman_layer],
        laser_frequency_cm1=20000.0,
        incident_angle_rad=0.0,
        incident_pol="p",
        detected_pol="p",
        temperature_K=0.0,
        linewidths_cm1=np.array([5.0, 5.0]),
        n_gauss=10,
        approximate_es=True,
        modal_pairs=True,
        modal_pair_combination=MODAL_PAIR_GROUP_Q,
        modes_selected=[False, True],
    )

    active_freqs, active_intensities, active_sigmas = calc.calculate_mode_intensities()

    np.testing.assert_allclose(active_freqs, [200.0])
    assert len(active_intensities) == 1
    np.testing.assert_allclose(active_sigmas, [5.0])


def test_nac_selection_uses_dominant_to_mode_map():
    """A TO-mode GUI selection should follow its sorted NAC branch."""
    selected_to_modes = [True, False, False]
    dominant_to_by_nac = np.array([2, 1, 0])

    nac_selected = _selection_for_nac_modes(selected_to_modes, dominant_to_by_nac)

    np.testing.assert_array_equal(nac_selected, [False, False, True])


def test_crystal_raman_layer_selection_can_follow_reordered_nac_modes():
    """Layer-specific NAC masks select reordered branches rather than raw TO indices."""
    frequencies_cm1 = np.array([400.0, 560.0])
    raman_tensors = [np.eye(3), 2.0 * np.eye(3)]
    system = build_system([(100e-9, 1.5)])
    raman_layer = RamanLayer(
        layer_index=0,
        phonon_frequencies_cm1=frequencies_cm1,
        raman_tensors=raman_tensors,
        rotation_matrix=np.eye(3),
        modes_selected=[False, True],
    )

    calc = LayeredRamanCalculator(
        system=system,
        raman_layers=[raman_layer],
        laser_frequency_cm1=20000.0,
        incident_angle_rad=0.0,
        incident_pol="p",
        detected_pol="p",
        temperature_K=0.0,
        linewidths_cm1=np.array([5.0, 5.0]),
        n_gauss=10,
        approximate_es=True,
        modes_selected=[True, False],
    )

    active_freqs, active_intensities, active_sigmas = calc.calculate_mode_intensities()

    np.testing.assert_allclose(active_freqs, [560.0])
    assert len(active_intensities) == 1
    np.testing.assert_allclose(active_sigmas, [5.0])


def test_crystal_raman_modal_pairs_selection_can_follow_reordered_nac_modes():
    """Modal-pair NAC masks must be applied after each q-resolved branch mapping."""
    frequencies_cm1 = np.array([400.0, 500.0])
    fallback_tensors = [np.eye(3), 2.0 * np.eye(3)]
    nac_tensors = [3.0 * np.eye(3), 4.0 * np.eye(3)]

    def nac_function(_q_hat_lab):
        return (
            np.array([390.0, 560.0]),
            nac_tensors,
            np.array([5.0, 5.0]),
            np.array([False, True]),
        )

    system = build_system([(100e-9, 1.5)])
    raman_layer = RamanLayer(
        layer_index=0,
        phonon_frequencies_cm1=frequencies_cm1,
        raman_tensors=fallback_tensors,
        rotation_matrix=np.eye(3),
        nac_function=nac_function,
        modes_selected=[True, False],
    )

    calc = LayeredRamanCalculator(
        system=system,
        raman_layers=[raman_layer],
        laser_frequency_cm1=20000.0,
        incident_angle_rad=0.0,
        incident_pol="p",
        detected_pol="p",
        temperature_K=0.0,
        linewidths_cm1=np.array([5.0, 5.0]),
        n_gauss=10,
        approximate_es=True,
        modal_pairs=True,
        modal_pair_combination=MODAL_PAIR_GROUP_Q,
        modes_selected=[True, False],
    )

    active_freqs, active_intensities, active_sigmas = calc.calculate_mode_intensities()

    assert np.any(np.isclose(active_freqs, 560.0))
    assert not np.any(np.isclose(active_freqs, 390.0))
    np.testing.assert_allclose(active_sigmas, 5.0 * np.ones(len(active_sigmas)))
