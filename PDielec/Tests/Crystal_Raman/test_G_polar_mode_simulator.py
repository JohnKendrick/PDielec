"""Tests for the tensor-level Raman polar mode simulator."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np

from PDielec.RamanPolarCalculator import (
    calculate_polar_curves,
    intensity_for_modes,
    raman_active_mode_indices,
)


def _rz(phi):
    """Rotation matrix for angle phi (rad) about the z axis."""
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def test_rotating_tensor_matches_rotating_polarisations():
    """Rotating the tensor and rotating the polarisations are equivalent."""
    tensor = np.array([[2.0, 0.3, 0.0], [0.3, 1.0, 0.0], [0.0, 0.0, 0.2]])
    e_l = np.array([1.0, 0.0, 0.0])
    e_s = np.array([0.0, 1.0, 0.0])
    for phi in np.linspace(0.0, 2.0 * np.pi, 17):
        q = _rz(phi)
        rotated_tensor = q @ tensor @ q.T
        intensity_tensor = abs(e_s @ rotated_tensor @ e_l) ** 2
        e_l_rot = q.T @ e_l
        e_s_rot = q.T @ e_s
        intensity_pol = abs(e_s_rot @ tensor @ e_l_rot) ** 2
        assert np.isclose(intensity_tensor, intensity_pol, atol=1.0e-12)


def test_summed_curve_is_sum_of_individual_curves():
    """Multiple selected modes are summed incoherently."""
    tensors = [
        np.diag([2.0, 1.0, 0.0]),
        np.array([[0.0, 1.5, 0.0], [1.5, 0.0, 0.0], [0.0, 0.0, 0.0]]),
    ]
    angles = np.linspace(0.0, 2.0 * np.pi, 49)
    combined = intensity_for_modes(tensors, [0, 1], [0, 0, 1], [1, 0, 0], angles, "VV")
    first = intensity_for_modes(tensors, [0], [0, 0, 1], [1, 0, 0], angles, "VV")
    second = intensity_for_modes(tensors, [1], [0, 0, 1], [1, 0, 0], angles, "VV")
    assert np.allclose(combined, first + second, atol=1.0e-12)


def test_active_mode_indices_preserve_one_based_table_numbering_basis():
    """Active filtering preserves original zero-based indices for one-based labels."""
    frequencies = [0.0, 8.0, 50.0, 60.0, 70.0]
    activities = [0.0, 10.0, 0.0, 1.0e-8, 2.0]
    tensors = [np.eye(3) for _i in frequencies]
    active = raman_active_mode_indices(frequencies, activities, tensors)
    assert active == [4]
    assert [index + 1 for index in active] == [5]


def test_normalised_curves_share_raw_maximum():
    """VV/HV/theta curves are normalised with one shared raw maximum."""
    tensors = [np.diag([2.0, 1.0, 0.0])]
    result = calculate_polar_curves(tensors, [0], angles_rad=np.linspace(0.0, 2.0 * np.pi, 37))
    assert set(result["curves"]) == {"VV", "HV", "theta"}
    assert np.isclose(result["raw_max_intensity"], 4.0, atol=1.0e-12)
    assert max(np.max(curve) for curve in result["curves"].values()) == 1.0
