"""Stage 2 tests for complete NAC branch diagnostics."""

import numpy as np

from PDielec.Constants import amu, wavenumber
from PDielec.GUI.CrystalScenarioTab import _compute_nac_dynamical_matrix_standalone
from PDielec.NACDiagnostics import (
    compute_nac_direction_diagnostics,
    group_degenerate_modes,
    track_nac_branches,
)


def _synthetic_inputs():
    """Return a three-mode polar system with an exactly degenerate TO pair."""
    frequencies = np.array([100.0, 100.0, 300.0])
    hessian = np.diag((frequencies * wavenumber) ** 2)
    born_charges = np.array([np.diag([2.0, 1.0, 0.5])])
    epsilon_inf = np.diag([4.0, 5.0, 6.0])
    volume_au = 100.0
    masses_au = np.array([20.0 * amu])
    u_to = np.eye(3)
    raman_tensors = [
        np.diag([1.0, 0.0, 0.0]),
        np.diag([0.0, 2.0, 0.0]),
        np.diag([0.0, 0.0, 3.0]),
    ]
    linewidths = np.array([4.0, 4.0, 6.0])
    chi2 = np.zeros((3, 3, 3))
    chi2[0, 0, 0] = 0.02
    chi2[1, 1, 0] = -0.01
    return (
        hessian,
        born_charges,
        epsilon_inf,
        volume_au,
        masses_au,
        u_to,
        raman_tensors,
        linewidths,
        chi2,
    )


def _diagnostics(q=(1.0, 0.0, 0.0), inputs=None):
    """Calculate diagnostics for the synthetic system."""
    if inputs is None:
        inputs = _synthetic_inputs()
    return compute_nac_direction_diagnostics(np.asarray(q), *inputs[:-1], chi2_repsilon=inputs[-1])


def test_every_nac_mode_has_complete_to_composition():
    """Every sorted NAC branch appears once and has unit total squared overlap."""
    result = _diagnostics()

    assert result["overlap_matrix"].shape == (3, 3)
    assert result["squared_overlap_matrix"].shape == (3, 3)
    np.testing.assert_allclose(result["squared_overlap_sums"], np.ones(3), atol=1.0e-14)
    assert sorted(np.concatenate(result["degenerate_groups"]).tolist()) == [0, 1, 2]
    assert len(result["raman_tensors_without_eo"]) == 3
    assert len(result["raman_tensors_with_eo"]) == 3
    assert result["raman_activities_without_eo"].shape == (3, 3)
    assert result["mode_effective_charge_vectors_e_sqrt_amu"].shape == (3, 3)


def test_legacy_nac_result_is_preserved():
    """Refactoring the NAC algebra into diagnostics preserves the legacy API."""
    inputs = _synthetic_inputs()
    q = np.array([1.0, 0.0, 0.0])
    result = _diagnostics(q, inputs)
    frequencies, tensors, linewidths, dominant = _compute_nac_dynamical_matrix_standalone(
        q, *inputs[:-1], chi2_repsilon=inputs[-1], return_mode_map=True
    )

    np.testing.assert_allclose(frequencies, result["frequencies_cm1"])
    np.testing.assert_allclose(tensors, result["raman_tensors_with_eo"])
    np.testing.assert_allclose(linewidths, result["linewidths_cm1"])
    np.testing.assert_array_equal(dominant, result["dominant_to_by_nac"])


def test_to_reordering_does_not_change_observable_nac_results():
    """Reordering the TO basis only reorders rows of the composition matrix."""
    inputs = list(_synthetic_inputs())
    reference = _diagnostics(inputs=inputs)
    permutation = np.array([2, 0, 1])
    inputs[5] = inputs[5][permutation]
    inputs[6] = [inputs[6][index] for index in permutation]
    inputs[7] = inputs[7][permutation]
    reordered = _diagnostics(inputs=inputs)

    np.testing.assert_allclose(reordered["frequencies_cm1"], reference["frequencies_cm1"])
    np.testing.assert_allclose(reordered["raman_tensors_without_eo"], reference["raman_tensors_without_eo"])
    np.testing.assert_allclose(reordered["raman_tensors_with_eo"], reference["raman_tensors_with_eo"])
    np.testing.assert_allclose(
        reordered["squared_overlap_matrix"], reference["squared_overlap_matrix"][permutation]
    )


def test_degenerate_to_rotation_does_not_change_observable_nac_results():
    """Rotating the exactly degenerate TO pair leaves tensors and activities invariant."""
    inputs = list(_synthetic_inputs())
    reference = _diagnostics(q=(0.3, 0.4, 0.5), inputs=inputs)
    angle = np.radians(37.0)
    rotation = np.array([
        [np.cos(angle), np.sin(angle), 0.0],
        [-np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    inputs[5] = rotation @ inputs[5]
    inputs[6] = list(np.einsum("ij,jab->iab", rotation, np.asarray(inputs[6])))
    rotated = _diagnostics(q=(0.3, 0.4, 0.5), inputs=inputs)

    np.testing.assert_allclose(rotated["frequencies_cm1"], reference["frequencies_cm1"])
    np.testing.assert_allclose(rotated["raman_tensors_without_eo"], reference["raman_tensors_without_eo"])
    np.testing.assert_allclose(rotated["raman_tensors_with_eo"], reference["raman_tensors_with_eo"])
    np.testing.assert_allclose(
        rotated["raman_activities_with_eo"], reference["raman_activities_with_eo"]
    )


def test_subspace_tracking_handles_a_rotated_degenerate_endpoint():
    """Projector continuity remains one when a degenerate pair rotates then splits."""
    angle = np.radians(43.0)
    rotated = np.array([
        [np.cos(angle), np.sin(angle), 0.0],
        [-np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    sequence = [
        {"frequencies_cm1": np.array([100.0, 100.0, 300.0]), "eigenvectors": np.eye(3)},
        {"frequencies_cm1": np.array([99.0, 101.0, 300.0]), "eigenvectors": rotated},
    ]

    tracking = track_nac_branches(sequence)

    assert sorted(tracking[1]["branch_ids_by_sorted_mode"].tolist()) == [0, 1, 2]
    assert all(match["projector_overlap"] > 1.0 - 1.0e-14 for match in tracking[1]["subspace_matches"])
    assert group_degenerate_modes(sequence[0]["frequencies_cm1"]) == [[0, 1], [2]]
