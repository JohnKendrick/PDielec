#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""Tensor-level Raman polar plot calculations."""

import numpy as np

RAMAN_ACTIVE_THRESHOLD = 1.0e-6
ACOUSTIC_THRESHOLD_CM1 = 10.0
POLAR_GEOMETRIES = ("VV", "HV", "theta")


def normalise_vector(vector, name="vector"):
    """Return a normalised 3-vector."""
    vector = np.asarray(vector, dtype=float)
    if vector.shape != (3,):
        raise ValueError(f"{name} must be a 3-vector")
    norm = np.linalg.norm(vector)
    if norm <= 0.0:
        raise ValueError(f"{name} must be non-zero")
    return vector / norm


def build_polarisation_basis(rotation_axis, zero_angle_polarisation):
    """Build the polar plot basis from the light direction and zero-angle polarisation.

    Parameters
    ----------
    rotation_axis : array_like, shape (3,)
        Light propagation direction.
    zero_angle_polarisation : array_like, shape (3,)
        Incoming polarisation direction at zero polar angle.

    Returns
    -------
    tuple
        ``(rz, pza, pp, orthogonality)`` where all vectors are normalised and
        ``orthogonality = dot(rz, pza_input_normalised)``.

    """
    rz = normalise_vector(rotation_axis, "rotation_axis")
    pza_input = normalise_vector(zero_angle_polarisation, "zero_angle_polarisation")
    orthogonality = float(np.dot(rz, pza_input))

    # Project Pza onto the plane perpendicular to Rz so the plotted basis is
    # physically meaningful even if the input vector is slightly non-orthogonal.
    pza = pza_input - orthogonality * rz
    pza_norm = np.linalg.norm(pza)
    if pza_norm <= 1.0e-12:
        raise ValueError("zero_angle_polarisation must not be parallel to rotation_axis")
    pza = pza / pza_norm
    pp = np.cross(rz, pza)
    pp = pp / np.linalg.norm(pp)
    return rz, pza, pp, orthogonality


def rotate_around_axis(vector, axis, angle_rad):
    """Rotate a vector around an axis using Rodrigues' rotation formula."""
    vector = np.asarray(vector, dtype=float)
    axis = normalise_vector(axis, "axis")
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return vector * c + np.cross(axis, vector) * s + axis * np.dot(axis, vector) * (1.0 - c)


def polarisation_vectors(rotation_axis, zero_angle_polarisation, angle_rad, geometry):
    """Return incident and scattered polarisation vectors for one polar angle."""
    if geometry not in POLAR_GEOMETRIES:
        raise ValueError(f"Unknown polar geometry {geometry}")
    rz, pza, pp, _orthogonality = build_polarisation_basis(rotation_axis, zero_angle_polarisation)
    v = rotate_around_axis(pza, rz, angle_rad)
    h = rotate_around_axis(pp, rz, angle_rad)
    if geometry == "VV":
        return v, v
    if geometry == "HV":
        return v, h
    return pza, v


def intensity_for_modes(raman_tensors, mode_indices, rotation_axis, zero_angle_polarisation, angles_rad, geometry):
    """Calculate a polar intensity curve for selected Raman modes.

    Contributions from multiple modes are summed incoherently.
    """
    if geometry not in POLAR_GEOMETRIES:
        raise ValueError(f"Unknown polar geometry {geometry}")
    tensors = [np.asarray(raman_tensors[index], dtype=complex) for index in mode_indices]
    intensities = np.zeros(len(angles_rad), dtype=float)
    for i, angle in enumerate(angles_rad):
        e_l, e_s = polarisation_vectors(rotation_axis, zero_angle_polarisation, angle, geometry)
        value = 0.0
        for tensor in tensors:
            amplitude = e_s @ tensor @ e_l
            value += float(np.real(amplitude * np.conj(amplitude)))
        intensities[i] = value
    return intensities


def calculate_polar_curves(
        raman_tensors,
        mode_indices,
        rotation_axis=(0.0, 0.0, 1.0),
        zero_angle_polarisation=(1.0, 0.0, 0.0),
        angles_rad=None,
        geometries=POLAR_GEOMETRIES,
        normalise=True):
    """Calculate Raman polar curves for selected modes and geometries.

    Returns
    -------
    dict
        Keys are ``angles_rad``, ``raw_curves``, ``curves``,
        ``raw_max_intensity`` and ``orthogonality``.

    """
    if angles_rad is None:
        angles_rad = np.linspace(0.0, 2.0 * np.pi, 361)
    angles_rad = np.asarray(angles_rad, dtype=float)
    mode_indices = list(mode_indices)
    _rz, _pza, _pp, orthogonality = build_polarisation_basis(rotation_axis, zero_angle_polarisation)

    raw_curves = {}
    for geometry in geometries:
        raw_curves[geometry] = intensity_for_modes(
            raman_tensors, mode_indices, rotation_axis, zero_angle_polarisation, angles_rad, geometry)

    raw_max = max((float(np.max(values)) for values in raw_curves.values()), default=0.0)
    if normalise and raw_max > 0.0:
        curves = {geometry: values / raw_max for geometry, values in raw_curves.items()}
    else:
        curves = {geometry: values.copy() for geometry, values in raw_curves.items()}

    return {
        "angles_rad": angles_rad,
        "raw_curves": raw_curves,
        "curves": curves,
        "raw_max_intensity": raw_max,
        "orthogonality": orthogonality,
    }


def raman_active_mode_indices(frequencies_cm1, raman_activities, raman_tensors, threshold=RAMAN_ACTIVE_THRESHOLD):
    """Return zero-based indices for Raman-active modes, preserving table order."""
    active = []
    n_tensors = len(raman_tensors) if raman_tensors is not None else 0
    for index, (frequency, activity) in enumerate(zip(frequencies_cm1, raman_activities)):
        if index >= n_tensors:
            continue
        if frequency > ACOUSTIC_THRESHOLD_CM1 and activity > threshold:
            active.append(index)
    return active
