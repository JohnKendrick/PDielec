"""Tests for Crystal Raman GUI layer coalescing."""

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
