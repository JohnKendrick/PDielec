"""Tests for Crystal Raman GUI layer coalescing."""

import pytest

from PDielec.GUI.CrystalScenarioTab import CrystalScenarioTab
from PDielec.GUI.SingleCrystalLayer import SingleCrystalLayer
from PDielec.Materials import Constant


def _layer(thickness_nm, dielectric_flag=True, azimuthal=0.0):
    material = Constant("Dielectric layer", permittivity=2.25)
    return SingleCrystalLayer(
        material,
        hkl=[0, 0, 0],
        azimuthal=azimuthal,
        thickness=thickness_nm,
        thickness_unit="nm",
        incoherentOption="Coherent",
        dielectricFlag=dielectric_flag,
    )


def _stack():
    return [
        _layer(1000.0, dielectric_flag=False),
        _layer(50.0),
        _layer(50.0),
        _layer(1000.0, dielectric_flag=False),
    ]


def test_equivalent_adjacent_dft_layers_are_coalesced():
    layers = [_layer(50.0), _layer(50.0)]

    coalesced = CrystalScenarioTab._coalesce_equivalent_layers(layers)

    assert len(coalesced) == 1
    assert coalesced[0].get_thickness() == 100.0
    assert coalesced[0].get_thickness_unit() == "nm"
    assert coalesced[0].is_dielectric()


def test_select_finite_layers_respects_disabled_coalescing():
    selected = CrystalScenarioTab._select_finite_raman_layers(_stack(), coalesce_equivalent=False)

    assert len(selected) == 2
    assert [layer.get_thickness() for layer in selected] == [50.0, 50.0]


def test_select_finite_layers_coalesces_when_enabled():
    selected = CrystalScenarioTab._select_finite_raman_layers(_stack(), coalesce_equivalent=True)

    assert len(selected) == 1
    assert selected[0].get_thickness() == 100.0


def test_layers_with_different_raman_activity_are_not_coalesced():
    layers = [_layer(50.0, dielectric_flag=True), _layer(50.0, dielectric_flag=False)]

    coalesced = CrystalScenarioTab._coalesce_equivalent_layers(layers)

    assert len(coalesced) == 2


def test_rotated_layers_are_not_coalesced():
    layers = [_layer(50.0, azimuthal=0.0), _layer(50.0, azimuthal=30.0)]

    coalesced = CrystalScenarioTab._coalesce_equivalent_layers(layers)

    assert len(coalesced) == 2


def test_normalised_depth_merges_equivalent_neighbours_even_if_coalescing_disabled():
    """Normalisation regions combine every mergeable run of GUI layers."""
    from PDielec.LayeredRamanCalculator import DEPTH_INTEGRATION_NORMALISED_INCOHERENT

    selected = CrystalScenarioTab._select_finite_raman_layers(
        _stack(), coalesce_equivalent=False,
        depth_integration=DEPTH_INTEGRATION_NORMALISED_INCOHERENT,
    )
    assert len(selected) == 1
    assert selected[0].get_thickness_in_metres() == pytest.approx(100e-9)


def test_normalised_depth_preserves_non_equivalent_layers():
    """Only layers allowed by the existing equivalence rules can merge."""
    from PDielec.LayeredRamanCalculator import DEPTH_INTEGRATION_NORMALISED_INCOHERENT

    layers = [_layer(1000), _layer(25), _layer(75, azimuthal=45), _layer(1000)]
    selected = CrystalScenarioTab._select_finite_raman_layers(
        layers, depth_integration=DEPTH_INTEGRATION_NORMALISED_INCOHERENT,
    )
    assert len(selected) == 2
    assert [layer.get_thickness() for layer in selected] == [25, 75]


def test_normalised_depth_gui_disables_layer_coherence_and_enables_final_state_controls(monkeypatch):
    """The new option retains incoherent-depth control behaviour when selected."""
    from types import SimpleNamespace

    from qtpy.QtWidgets import QApplication, QComboBox, QLabel, QSpinBox

    from PDielec.LayeredRamanCalculator import (
        DEPTH_INTEGRATION_COHERENT,
        DEPTH_INTEGRATION_NORMALISED_INCOHERENT,
        DEPTH_INTEGRATION_OPTIONS,
        FINAL_STATE_BULK_PHASE_MATCHED,
    )

    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance() or QApplication([])
    scenario = SimpleNamespace(
        settings={"Layer NAC mode": "modal_pairs", "Modal pair final-state model": FINAL_STATE_BULK_PHASE_MATCHED},
        depth_integration_cb=QComboBox(), layer_combination_cb=QComboBox(), layer_combination_label=QLabel(),
        modal_pair_final_state_cb=QComboBox(), modal_pair_final_state_label=QLabel(),
        modal_pair_q_tol_sb=QSpinBox(), modal_pair_q_tol_label=QLabel(),
    )
    scenario._update_modal_final_state_controls = lambda: CrystalScenarioTab._update_modal_final_state_controls(scenario)
    scenario.depth_integration_cb.addItems(DEPTH_INTEGRATION_OPTIONS)
    scenario.depth_integration_cb.setCurrentText(DEPTH_INTEGRATION_NORMALISED_INCOHERENT)
    CrystalScenarioTab.on_depth_integration_cb_activated(scenario, scenario.depth_integration_cb.currentIndex())
    assert scenario.settings["Depth coherence"] == DEPTH_INTEGRATION_NORMALISED_INCOHERENT
    assert not scenario.layer_combination_cb.isEnabled()
    assert scenario.modal_pair_final_state_cb.isEnabled()
    assert scenario.modal_pair_q_tol_sb.isEnabled()
    assert scenario.calculation_required and scenario.refresh_required
    scenario.depth_integration_cb.setCurrentText(DEPTH_INTEGRATION_COHERENT)
    CrystalScenarioTab.on_depth_integration_cb_activated(scenario, scenario.depth_integration_cb.currentIndex())
    assert scenario.layer_combination_cb.isEnabled()
    assert not scenario.modal_pair_final_state_cb.isEnabled()
    app.processEvents()
