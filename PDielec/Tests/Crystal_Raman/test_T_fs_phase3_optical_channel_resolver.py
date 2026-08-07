#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Phase-3 synthetic tests for OpticalChannelResolver diagnostics."""

import numpy as np
import numpy.testing as npt
import pytest

from PDielec.LayeredRamanCalculator import MODAL_PAIR_COHERENT_ALL, MODAL_PAIR_INCOHERENT, LayeredRamanCalculator
from PDielec.OpticalChannelResolver import OpticalChannelResolver
from PDielec.RamanGeometry import resolve_collection_angle
from PDielec.Tests.Crystal_Raman.conftest import build_system, make_raman_layer


def _calculator(
    collection_side="superstrate",
    incident_angle_rad=0.0,
    layer_specs=None,
    layer_index=0,
    modal_pairs=False,
    modal_pair_combination=MODAL_PAIR_COHERENT_ALL,
    approximate_es=True,
    collection_angle_rad=None,
):
    if layer_specs is None:
        layer_specs = [(1.0e-6, 1.5)]
    system = build_system(layer_specs, n_sup=1.0, n_sub=1.0)
    raman_layer = make_raman_layer(layer_index, np.eye(3), 100.0)
    return LayeredRamanCalculator(
        system=system,
        raman_layers=[raman_layer],
        laser_frequency_cm1=20000.0,
        incident_angle_rad=incident_angle_rad,
        incident_pol="p",
        detected_pol="p",
        temperature_K=0.0,
        linewidths_cm1=np.array([5.0]),
        n_gauss=8,
        approximate_es=approximate_es,
        collection_side=collection_side,
        collection_angle_rad=collection_angle_rad,
        modal_pairs=modal_pairs,
        modal_pair_combination=modal_pair_combination,
        coherent_layers=True,
    )


@pytest.mark.parametrize(
    ("incident_angle", "collection_angle", "collection_side", "expected"),
    [
        (0.0, None, "superstrate", 0.0),
        (0.0, None, "substrate", 0.0),
        (0.4, None, "superstrate", -0.4),
        (0.4, 0.4, "superstrate", 0.4),
        (0.4, None, "substrate", 0.4),
        (0.4, -0.2, "superstrate", -0.2),
    ],
)
def test_collection_angle_resolver(incident_angle, collection_angle, collection_side, expected):
    """Automatic and explicit collection geometries retain their signed conventions."""
    assert resolve_collection_angle(incident_angle, collection_angle, collection_side) == pytest.approx(expected)


def test_optical_resolver_modal_channels_reconstruct_total_field():
    """Summing q-channel fields reproduces the full GTM field for each polarisation."""
    calc = _calculator()
    resolver = OpticalChannelResolver()
    layer_indices = [rl.layer_index for rl in calc.raman_layers]

    solution = resolver.solve(
        calc.system,
        calc.laser_frequency_cm1,
        calc.incident_angle_rad,
        calc._gl_z,
        layer_indices,
        calc._gl_layer_slices,
        pol_indices=(0, 1),
    )
    layer_slice = calc._gl_layer_slices[0]

    for pol_idx in (0, 1):
        channel_sum = sum(channel.field for channel in solution.channels_by_pol[pol_idx][0])
        expected = solution.total_fields[3 * pol_idx:3 * pol_idx + 3, :]
        npt.assert_allclose(channel_sum[:, layer_slice], expected[:, layer_slice], atol=1e-10, rtol=1e-10)


def test_optical_resolver_external_q_ext_backscatter_and_forward_limits():
    """External momentum transfer distinguishes backscatter from forward scattering."""
    resolver = OpticalChannelResolver()
    backscatter = _calculator(collection_side="superstrate")
    forward = _calculator(collection_side="substrate")
    reversed_system = forward.system.reversed_system()

    q_back = resolver.external_momentum_transfer(
        backscatter.system,
        20000.0,
        0.0,
        backscatter.system,
        19900.0,
        0.0,
        "superstrate",
    )
    q_forward = resolver.external_momentum_transfer(
        forward.system,
        20000.0,
        0.0,
        reversed_system,
        19900.0,
        0.0,
        "substrate",
    )

    npt.assert_allclose(q_back, [0.0, 0.0, 2.0], atol=1e-12)
    npt.assert_allclose(q_forward, [0.0, 0.0, 0.0], atol=1e-12)


def test_oblique_automatic_retro_specular_and_forward_q_ext():
    """Oblique automatic collection is retro or forward while explicit reflection stays specular."""
    angle = np.radians(30.0)
    automatic_retro = _calculator(collection_side="superstrate", incident_angle_rad=angle)
    explicit_specular = _calculator(
        collection_side="superstrate",
        incident_angle_rad=angle,
        collection_angle_rad=angle,
    )
    automatic_forward = _calculator(collection_side="substrate", incident_angle_rad=angle)

    q_retro = automatic_retro.resolve_optical_channels(19900.0)["q_ext"]
    q_specular = explicit_specular.resolve_optical_channels(19900.0)["q_ext"]
    q_forward = automatic_forward.resolve_optical_channels(19900.0)["q_ext"]

    npt.assert_allclose(q_retro, [1.0, 0.0, np.sqrt(3.0)], atol=1e-12)
    npt.assert_allclose(q_specular, [0.0, 0.0, np.sqrt(3.0)], atol=1e-12)
    npt.assert_allclose(q_forward, [0.0, 0.0, 0.0], atol=1e-12)


def test_layered_calculator_exposes_phase3_optical_diagnostics():
    """LayeredRamanCalculator exposes incident/detected optical channel diagnostics."""
    calc = _calculator()

    diagnostics = calc.resolve_optical_channels(scattered_frequency_cm1=19900.0)

    assert set(diagnostics) == {"incident", "detected", "q_ext", "z_detected", "collection_system"}
    assert diagnostics["incident"].total_fields.shape == (6, len(calc._gl_z))
    assert diagnostics["detected"].total_fields.shape == (6, len(calc._gl_z))
    assert 0 in diagnostics["incident"].channels_by_pol[0]
    npt.assert_allclose(diagnostics["q_ext"], [0.0, 0.0, 2.0], atol=1e-12)


def test_forward_scatter_diagnostics_map_physical_layer_to_reversed_system_layer():
    """A Raman layer with nonzero physical index is retained in the reversed stack."""
    calc = _calculator(collection_side="substrate", layer_specs=[(0.4e-6, 1.3), (0.6e-6, 1.7)], layer_index=1)

    diagnostics = calc.resolve_optical_channels(scattered_frequency_cm1=19900.0)

    assert 1 in diagnostics["detected"].channels_by_pol[0]
    assert diagnostics["detected"].channels_by_pol[0][1]
    qz_values = [np.real(channel.qz) for channel in diagnostics["detected"].channels_by_pol[0][1]]
    assert any(abs(abs(qz) - 1.7) < 1e-2 for qz in qz_values)
    assert not any(abs(abs(qz) - 1.3) < 1e-2 for qz in qz_values)


def test_phase3_diagnostic_path_does_not_change_existing_intensity_calculation():
    """Calling optical diagnostics leaves the existing Raman calculation unchanged."""
    calc = _calculator()

    before = calc.calculate_mode_intensities()
    _diagnostics = calc.resolve_optical_channels(scattered_frequency_cm1=19900.0)
    after = calc.calculate_mode_intensities()

    for before_arr, after_arr in zip(before, after, strict=True):
        npt.assert_allclose(after_arr, before_arr, atol=0.0, rtol=0.0)


def test_phase3_modal_diagnostic_matches_existing_coherent_all_modal_pairs():
    """Phase-3 resolver/accumulator diagnostic agrees with simple modal-pairs case."""
    calc = _calculator(modal_pairs=True)

    existing = calc.calculate_mode_intensities()
    diagnostic = calc.calculate_phase3_modal_pair_diagnostic_intensities()

    for existing_arr, diagnostic_arr in zip(existing, diagnostic, strict=True):
        npt.assert_allclose(diagnostic_arr, existing_arr, atol=1e-8, rtol=1e-10)


def test_phase3_modal_diagnostic_rejects_incoherent_degenerate_channel_case():
    """Subspace channels cannot reproduce raw per-mode incoherent modal-pair sums."""
    calc = _calculator(modal_pairs=True, modal_pair_combination=MODAL_PAIR_INCOHERENT)

    with pytest.raises(NotImplementedError, match="per-mode incoherent"):
        calc.calculate_phase3_modal_pair_diagnostic_intensities()


def test_phase3_modal_diagnostic_rejects_non_approximate_es_for_now():
    """Diagnostic does not silently use laser-frequency detector fields when nu_S is needed."""
    calc = _calculator(modal_pairs=True, approximate_es=False)

    with pytest.raises(NotImplementedError, match="approximate_es=True"):
        calc.calculate_phase3_modal_pair_diagnostic_intensities()
