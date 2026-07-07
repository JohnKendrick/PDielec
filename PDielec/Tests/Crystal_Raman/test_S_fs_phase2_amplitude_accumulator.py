#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Phase-2 synthetic tests for RamanAmplitudeAccumulator."""

import numpy as np
import pytest

from PDielec.PhononFinalStateResolver import (
    COHERENT_FILM,
    CoherentExternalChannelContribution,
    DiscardedInternalComponent,
    FinalState,
    ResolvedFinalState,
)
from PDielec.RamanAmplitudeAccumulator import RamanAmplitudeAccumulator, RamanAmplitudeRecord


def test_same_final_state_and_group_sums_amplitudes_coherently():
    """Amplitudes assigned to one key/group are summed before squaring."""
    accumulator = RamanAmplitudeAccumulator()

    accumulator.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="layer-group")
    accumulator.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="layer-group")

    assert accumulator.amplitude("mode-a", "layer-group") == 2.0 + 0.0j
    assert accumulator.intensities()["mode-a"] == pytest.approx(4.0)
    assert accumulator.total_intensity() == pytest.approx(4.0)


def test_different_final_states_sum_incoherently():
    """Different final-state keys are squared separately and summed."""
    accumulator = RamanAmplitudeAccumulator()

    accumulator.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="same-layer")
    accumulator.add_amplitude("mode-b", 1.0 + 0.0j, coherence_group="same-layer")

    assert accumulator.intensities() == {"mode-a": 1.0, "mode-b": 1.0}
    assert accumulator.total_intensity() == pytest.approx(2.0)


def test_fs11_layer_coherence_group_controls_cross_layer_terms():
    """Same final state interferes only across layers in the same coherence group."""
    coherent = RamanAmplitudeAccumulator()
    coherent.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="coherent-stack")
    coherent.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="coherent-stack")

    incoherent = RamanAmplitudeAccumulator()
    incoherent.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="layer-0")
    incoherent.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="layer-1")

    assert coherent.intensities()["mode-a"] == pytest.approx(4.0)
    assert incoherent.intensities()["mode-a"] == pytest.approx(2.0)


def test_incoherent_depth_local_intensity_is_added_after_squaring():
    """Local intensities are added without reintroducing phase interference."""
    accumulator = RamanAmplitudeAccumulator()

    accumulator.add_incoherent_intensity("mode-a", 2.0)
    accumulator.add_incoherent_intensity("mode-a", 3.0)
    accumulator.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="thin-slice")

    assert accumulator.intensities()["mode-a"] == pytest.approx(6.0)


def test_records_can_build_accumulator():
    """Batch records use the same coherent/incoherent boundaries."""
    records = [
        RamanAmplitudeRecord("mode-a", amplitude=1.0 + 0.0j, coherence_group="g0"),
        RamanAmplitudeRecord("mode-a", amplitude=2.0 + 0.0j, coherence_group="g0"),
        RamanAmplitudeRecord("mode-a", amplitude=3.0 + 0.0j, coherence_group="g1"),
        RamanAmplitudeRecord("mode-b", local_intensity=5.0),
    ]

    accumulator = RamanAmplitudeAccumulator.from_records(records)

    assert accumulator.intensities()["mode-a"] == pytest.approx(18.0)
    assert accumulator.intensities()["mode-b"] == pytest.approx(5.0)
    assert accumulator.total_intensity() == pytest.approx(23.0)


def test_legacy_amplitudes_raises_for_multi_group_final_state():
    """Legacy amplitude view must not silently drop coherence groups."""
    accumulator = RamanAmplitudeAccumulator()
    accumulator.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="g0")
    accumulator.add_amplitude("mode-a", 2.0 + 0.0j, coherence_group="g1")

    with pytest.raises(ValueError, match="multi-group"):
        accumulator.amplitudes()


def test_final_state_amplitudes_raises_for_multi_group_final_state():
    """Final-state amplitude view must not coherently collapse incoherent groups."""
    accumulator = RamanAmplitudeAccumulator()
    accumulator.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="g0")
    accumulator.add_amplitude("mode-a", 2.0 + 0.0j, coherence_group="g1")

    with pytest.raises(ValueError, match="multi-group"):
        accumulator.final_state_amplitudes()


def test_coherent_final_state_amplitudes_is_explicit_diagnostic_collapse():
    """The dangerous coherent group collapse is available only under an explicit name."""
    accumulator = RamanAmplitudeAccumulator()
    accumulator.add_amplitude("mode-a", 1.0 + 0.0j, coherence_group="g0")
    accumulator.add_amplitude("mode-a", 2.0 + 0.0j, coherence_group="g1")

    assert accumulator.coherent_final_state_amplitudes()["mode-a"] == 3.0 + 0.0j


def test_resolver_classifications_can_be_added_directly():
    """Resolved and coherent external-channel classifications map to accumulator keys."""
    final_state = FinalState(("layer", "det", 0, 100.0, ("q_dir", 0, 0)))
    resolved = ResolvedFinalState(final_state, q_ph=np.array([0.0, 0.0, 1.0]), phase_weight=2.0 + 0.0j)
    coherent = CoherentExternalChannelContribution(
        ("layer", "det", 0, COHERENT_FILM),
        q_ph=np.array([0.0, 0.0, 0.0]),
        phase_weight=0.5 + 0.0j,
    )
    discarded = DiscardedInternalComponent("not retained")
    accumulator = RamanAmplitudeAccumulator()

    accumulator.add_classified(resolved, 1.0 + 0.0j, coherence_group="g0")
    accumulator.add_classified(coherent, 2.0 + 0.0j, coherence_group="g1")
    accumulator.add_classified(discarded, 100.0 + 0.0j, coherence_group="g2")

    assert accumulator.intensities()[final_state.key] == pytest.approx(4.0)
    assert accumulator.intensities()[coherent.key] == pytest.approx(1.0)
    assert accumulator.total_intensity() == pytest.approx(5.0)


def test_legacy_add_api_still_sums_coherently():
    """The original add/key API remains compatible with FS3 tests."""
    accumulator = RamanAmplitudeAccumulator()

    accumulator.add("mode-a", 0.5 + 0.0j)
    accumulator.add("mode-a", 0.5 + 0.0j)

    assert accumulator.amplitude("mode-a") == 1.0 + 0.0j
    assert accumulator.total_intensity() == pytest.approx(1.0)
