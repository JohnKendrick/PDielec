#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Phase-1 synthetic tests for PhononFinalStateResolver."""

import numpy as np
import numpy.testing as npt
import pytest

from PDielec.PhononFinalStateResolver import (
    BULK_PHASE_MATCHED,
    COHERENT_EXTERNAL_CHANNEL,
    COHERENT_FILM,
    DISCARDED_INTERNAL_COMPONENT,
    INCOHERENT_DEPTH,
    RESOLVED_FINAL_STATE,
    PhononFinalStateResolver,
)


def test_fs4_nac_branch_mapping_uses_overlap_not_sorted_index():
    """FS4: swapped NAC eigenvalue order must still map to the correct branch."""
    resolver = PhononFinalStateResolver()
    to_modes = np.eye(3, dtype=complex)
    nac_modes = np.column_stack([to_modes[:, 1], to_modes[:, 0], to_modes[:, 2]])
    frequencies = np.array([300.0, 550.0, 700.0])
    tensors = [np.eye(3) * value for value in [3.0, 5.5, 7.0]]
    linewidths = [3.0, 5.0, 7.0]
    selected = [False, True, True]

    branch = resolver.map_nac_branch(
        0,
        to_modes,
        nac_modes,
        frequencies,
        tensors,
        nac_linewidths_cm1=linewidths,
        nac_selected=selected,
    )

    assert branch.branch_index == 1
    assert branch.frequency_cm1 == 550.0
    assert branch.linewidth_cm1 == 5.0
    assert branch.selected is True
    npt.assert_allclose(branch.raman_tensor, 5.5 * np.eye(3))
    npt.assert_allclose(branch.overlap, 1.0)


def test_fs9_selection_follows_dominant_to_by_nac_map():
    """FS9: selection is read from the mapped NAC branch, not raw sorted index."""
    resolver = PhononFinalStateResolver()
    to_modes = np.eye(2, dtype=complex)
    nac_modes = np.column_stack([to_modes[:, 1], to_modes[:, 0]])
    dominant_to_by_nac = [1, 0]
    selected = [False, True]

    branch = resolver.map_nac_branch(
        0,
        to_modes,
        nac_modes,
        [300.0, 550.0],
        [np.eye(3), 2.0 * np.eye(3)],
        nac_selected=selected,
        dominant_to_by_nac=dominant_to_by_nac,
    )

    assert branch.branch_index == 1
    assert branch.selected is True
    npt.assert_allclose(branch.overlap, 1.0)


def test_dominant_to_by_nac_path_reports_overlap():
    """Precomputed NAC maps still expose the mapped branch overlap for diagnostics."""
    resolver = PhononFinalStateResolver()
    to_modes = np.eye(2, dtype=complex)
    nac_modes = np.column_stack([to_modes[:, 1], to_modes[:, 0]])

    branch_index, overlap = resolver.nac_branch_index(
        0,
        to_modes,
        nac_modes,
        dominant_to_by_nac=[1, 0],
    )

    assert branch_index == 1
    npt.assert_allclose(overlap, 1.0)


def test_low_overlap_nac_mapping_warns_and_falls_back_to_frequency():
    """Low-overlap NAC matching falls back to closest frequency when available."""
    resolver = PhononFinalStateResolver(min_nac_overlap=0.5)
    to_modes = np.eye(3, dtype=complex)
    nac_modes = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=complex,
    )
    tensors = [np.eye(3), 2.0 * np.eye(3), 3.0 * np.eye(3)]

    with pytest.warns(UserWarning, match="below threshold"):
        branch = resolver.map_nac_branch(
            0,
            to_modes,
            nac_modes,
            [100.0, 210.0, 300.0],
            tensors,
            to_frequencies_cm1=[205.0, 400.0, 500.0],
        )

    assert branch.branch_index == 1
    assert branch.frequency_cm1 == 210.0
    npt.assert_allclose(branch.overlap, 0.0)


def test_fs6_nonpolar_key_ignores_q_direction_and_magnitude():
    """FS6: non-polar modes are q independent final states."""
    resolver = PhononFinalStateResolver(coherence_regime=BULK_PHASE_MATCHED, q_ext=[0.0, 0.0, 2.0])

    state_a = resolver.resolve_pair(2, [0.0, 0.0, 2.0], is_polar=False, frequency_cm1=120.0)
    state_b = resolver.resolve_pair(2, [1.0, 0.0, 0.0], is_polar=False, frequency_cm1=120.0)

    assert state_a.kind == RESOLVED_FINAL_STATE
    assert state_b.kind == RESOLVED_FINAL_STATE
    assert state_a.final_state.key == state_b.final_state.key
    assert state_a.final_state.key[-1] == ("q_independent",)


def test_nonpolar_classification_uses_relative_oscillator_strength():
    """Weak oscillator strength below threshold is classified as non-polar."""
    resolver = PhononFinalStateResolver(nonpolar_strength_threshold=1e-4)

    assert resolver.is_polar_mode(1.0, 10.0)
    assert not resolver.is_polar_mode(1.0e-6, 10.0)


def test_fs7_q_class_key_uses_direction_not_magnitude():
    """FS7: q class is stable for same direction and different magnitude."""
    resolver = PhononFinalStateResolver(angular_tolerance_deg=5.0)

    assert resolver.q_class_key([0.0, 0.0, 1.0]) == resolver.q_class_key([0.0, 0.0, 10.0])
    assert resolver.q_class_key([1.0, 0.0, 0.0]) != resolver.q_class_key([-1.0, 0.0, 0.0])


def test_fs14_near_zero_q_rules_depend_on_coherence_regime():
    """FS14: internal zero-q backscatter cross-pairs are not independent TO states."""
    q_ext = [0.0, 0.0, 2.0]
    kwargs = {"mode_idx": 0, "q_ph": [0.0, 0.0, 0.0], "is_polar": True, "reference_k": 1.0}

    coherent = PhononFinalStateResolver(coherence_regime=COHERENT_FILM, q_ext=q_ext).resolve_pair(**kwargs)
    incoherent = PhononFinalStateResolver(coherence_regime=INCOHERENT_DEPTH, q_ext=q_ext).resolve_pair(**kwargs)
    bulk = PhononFinalStateResolver(coherence_regime=BULK_PHASE_MATCHED, q_ext=q_ext).resolve_pair(**kwargs)

    assert coherent.kind == COHERENT_EXTERNAL_CHANNEL
    assert incoherent.kind == DISCARDED_INTERNAL_COMPONENT
    assert bulk.kind == DISCARDED_INTERNAL_COMPONENT


def test_forward_zero_q_polar_mode_is_resolved_final_state_in_bulk():
    """Near-zero q is a genuine TO-like final state when q_ext is also near zero."""
    resolver = PhononFinalStateResolver(coherence_regime=BULK_PHASE_MATCHED, q_ext=[0.0, 0.0, 0.0])

    state = resolver.resolve_pair(0, [0.0, 0.0, 0.0], is_polar=True, reference_k=1.0)

    assert state.kind == RESOLVED_FINAL_STATE
    assert state.final_state.key[-1] == ("q_zero",)


def test_forward_nonzero_internal_q_is_not_resolved_when_external_q_is_zero():
    """A non-zero internal polar q is not a separate forward-scattering final state."""
    resolver = PhononFinalStateResolver(coherence_regime=BULK_PHASE_MATCHED, q_ext=[0.0, 0.0, 0.0])

    state = resolver.resolve_pair(0, [0.0, 0.0, 2.0], is_polar=True, reference_k=1.0)

    assert state.kind == DISCARDED_INTERNAL_COMPONENT


def test_bulk_phase_matched_keeps_only_q_matching_external_direction():
    """Bulk phase matching keeps q parallel to q_ext and discards other directions."""
    resolver = PhononFinalStateResolver(coherence_regime=BULK_PHASE_MATCHED, q_ext=[0.0, 0.0, 2.0])

    matched = resolver.resolve_pair(0, [0.0, 0.0, 4.0], is_polar=True)
    mismatched = resolver.resolve_pair(0, [2.0, 0.0, 0.0], is_polar=True)

    assert matched.kind == RESOLVED_FINAL_STATE
    assert mismatched.kind == DISCARDED_INTERNAL_COMPONENT


def test_bulk_phase_matched_keys_accepted_pairs_by_external_final_state():
    """Accepted internal pairs share the externally selected final-state key."""
    resolver = PhononFinalStateResolver(
        coherence_regime=BULK_PHASE_MATCHED,
        q_ext=[0.0, 0.0, 2.0],
        angular_tolerance_deg=90.0,
    )

    parallel = resolver.resolve_pair(0, [0.0, 0.0, 4.0], is_polar=True)
    oblique = resolver.resolve_pair(0, [1.0, 0.0, 1.0], is_polar=True)

    assert parallel.kind == RESOLVED_FINAL_STATE
    assert oblique.kind == RESOLVED_FINAL_STATE
    assert parallel.final_state.key == oblique.final_state.key
    npt.assert_allclose(parallel.q_ph, [0.0, 0.0, 4.0])
    npt.assert_allclose(oblique.q_ph, [1.0, 0.0, 1.0])


def test_bulk_matching_angle_is_separate_from_q_grouping_angle():
    """The 90 degree bulk acceptance must not coarsen q-direction bins."""
    resolver = PhononFinalStateResolver(
        coherence_regime=BULK_PHASE_MATCHED,
        q_ext=[0.0, 0.0, 1.0],
        angular_tolerance_deg=1.0,
        matching_tolerance_deg=90.0,
    )

    perpendicular = resolver.resolve_pair(0, [1.0, 0.0, 0.0], is_polar=True)

    assert perpendicular.kind == RESOLVED_FINAL_STATE
    assert resolver.q_class_key([0.0, 0.0, 1.0]) != resolver.q_class_key([1.0, 0.0, 0.0])


def test_phonon_subspace_sums_tensors_for_nearly_degenerate_branches():
    """Degenerate phonon branches are represented by a tensor-summed subspace."""
    resolver = PhononFinalStateResolver(eps_phonon=1e-3)
    tensors = [np.eye(3), 2.0 * np.eye(3), 10.0 * np.eye(3)]

    subspaces = resolver.classify_phonon_subspaces([100.0, 100.05, 120.0], tensors)

    assert [space.indices for space in subspaces] == [(0, 1), (2,)]
    npt.assert_allclose(subspaces[0].raman_tensor, 3.0 * np.eye(3))


def test_phonon_subspace_grouping_uses_connected_components_not_seed_only():
    """Transitive near-degenerate phonon branches form one branch subspace."""
    resolver = PhononFinalStateResolver(eps_phonon=1.0e-3)
    tensors = [value * np.eye(3) for value in [1.0, 2.0, 3.0, 10.0]]

    subspaces = resolver.classify_phonon_subspaces([100.0, 100.075, 100.15, 120.0], tensors)

    assert [space.indices for space in subspaces] == [(0, 1, 2), (3,)]
    npt.assert_allclose(subspaces[0].raman_tensor, 6.0 * np.eye(3))


def test_optical_subspace_grouping_uses_connected_components_not_seed_only():
    """Three chained near-degenerate kz values must form one optical subspace."""
    resolver = PhononFinalStateResolver(eps_optical=1.0e-4)
    kz_values = [1.0, 1.0 + 0.75e-4, 1.0 + 1.5e-4, 1.2]
    modes = np.eye(4, dtype=complex)

    subspaces = resolver.classify_optical_subspaces(kz_values, modes=modes)

    assert [space.indices for space in subspaces] == [(0, 1, 2), (3,)]
    npt.assert_allclose(subspaces[0].projector, np.diag([1.0, 1.0, 1.0, 0.0]), atol=1e-14)


def test_subspace_amplitude_projects_detector_degenerate_subspace():
    """Detector-side basis rotations must not change the projected amplitude."""
    resolver = PhononFinalStateResolver()
    raman_tensor = np.array([[0.1, 0.2, 0.0], [0.2, 0.5, 0.0], [0.0, 0.0, 2.0]], dtype=complex)
    incident_modes = np.eye(3, 2, dtype=complex)
    incident_field = np.array([1.0 + 0.5j, -0.3 + 0.2j, 3.0 + 4.0j])
    detector_field = np.array([0.7 - 0.1j, 0.2 + 0.9j, 5.0 - 2.0j])
    detector_modes = np.eye(3, 2, dtype=complex)
    theta = np.radians(37.0)
    rotation = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    detector_modes_rot = detector_modes @ rotation.T

    amplitude = resolver.subspace_amplitude(
        raman_tensor,
        incident_modes,
        incident_field,
        detector_field,
        detector_modes=detector_modes,
    )
    amplitude_rot = resolver.subspace_amplitude(
        raman_tensor,
        incident_modes,
        incident_field,
        detector_field,
        detector_modes=detector_modes_rot,
    )

    npt.assert_allclose(amplitude_rot, amplitude, atol=1e-14)
