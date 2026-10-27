#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Optical channel decomposition helpers for layered Raman calculations."""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence
from copy import deepcopy

import numpy as np

from PDielec.Constants import speed_light_si
from PDielec.PhononFinalStateResolver import PhononFinalStateResolver


def reciprocal_system(system, collection_side):
    """Return the detector stack in coordinates (x, y, L-z).

    Reflect permittivity tensors as P epsilon P, P=diag(1,1,-1), on
    independent layer copies. Electric vectors are mapped back by P; the
    optical solver handles magnetic axial vectors in its own coordinate frame.
    """
    if collection_side == "superstrate":
        return system
    if collection_side != "substrate":
        raise ValueError("unknown collection side")
    reversed_stack = deepcopy(system.reversed_system())
    reflection = np.diag([1.0, 1.0, -1.0])
    for layer in [reversed_stack.superstrate, *reversed_stack.layers, reversed_stack.substrate]:
        layer.euler = reflection @ layer.euler
        layer.euler_inverse = layer.euler_inverse @ reflection
    reversed_stack._raman_depth_reversed = True
    return reversed_stack


def electric_fields_to_lab(fields, system, component_axis=0):
    """Map electric components from a depth-reflected solver to the lab frame."""
    if not getattr(system, "_raman_depth_reversed", False):
        return fields
    result = fields.copy()
    indices = [slice(None)] * result.ndim
    for component in range(2, result.shape[component_axis], 3):
        indices[component_axis] = component
        result[tuple(indices)] *= -1
    return result


def modal_pair_momentum(channel_l, channel_s, stokes_ratio, collection_side):
    """Return real overlap momentum in common laser-vacuum wavevector units.

    Channel zeta/qz belong to the incident and reciprocal solver frames.
    This real-vector convention is for NAC/legacy angular classification;
    attenuation remains in the complex fields, not in a phonon-state label.
    """
    sign = 1.0 if collection_side == "superstrate" else -1.0
    return np.array([
        np.real(channel_l.get("zeta", 0.0) + stokes_ratio * channel_s.get("zeta", 0.0)),
        0.0,
        np.real(channel_l["qz"] + sign * stokes_ratio * channel_s["qz"]),
    ])


@dataclasses.dataclass(frozen=True)
class OpticalChannel:
    """One qz-resolved optical field component inside a finite layer."""

    layer_index: int
    qz: complex
    field: np.ndarray
    mode_indices: tuple[int, ...]


@dataclasses.dataclass(frozen=True)
class OpticalChannelSolution:
    """Optical fields and modal channels for one external optical problem.

    ``channels_by_pol`` is nested as
    ``channels_by_pol[pol_idx][physical_layer_index] -> list[OpticalChannel]``.
    The physical layer index is the index in the original, unreversed stack.
    Electric components are in the laboratory frame. Coordinates, zeta and qz
    remain in the solver frame; substrate collection uses z_prime = L - z.
    """

    frequency_cm1: float
    zeta: complex
    z_points: np.ndarray
    total_fields: np.ndarray
    modal_fields: np.ndarray
    qs_by_layer: dict[int, np.ndarray]
    channels_by_pol: dict[int, dict[int, list[OpticalChannel]]]


class OpticalChannelResolver:
    """Resolve full fields and Berreman q-channel components using GTM."""

    def __init__(
        self,
        optical_subspace_tolerance: float = 1.0e-4,
    ) -> None:
        """Initialise optical subspace grouping tolerance."""
        self._subspace_resolver = PhononFinalStateResolver(eps_optical=optical_subspace_tolerance)

    def solve(
        self,
        system,
        frequency_cm1: float,
        angle_rad: float,
        z_points: Sequence[float],
        layer_indices: Sequence[int],
        layer_slices: Sequence[slice],
        pol_indices: Sequence[int] = (0, 1),
        system_layer_indices: Sequence[int] | None = None,
    ) -> OpticalChannelSolution:
        """Compute total fields, modal fields, and q-channel groups.

        ``layer_indices`` are physical layer ids used in the returned channel
        dictionaries.  ``system_layer_indices`` can be supplied when solving on
        a transformed optical system, such as the reversed stack used for
        forward-scattering reciprocal fields.
        """
        z_array = np.asarray(z_points, dtype=float)
        freq_hz = float(frequency_cm1) * speed_light_si * 1.0e2
        system.initialize_sys(freq_hz)
        zeta = np.sin(float(angle_rad)) * np.sqrt(system.superstrate.epsilon[0, 0])
        _, total_fields, _ = system.calculate_Efield(freq_hz, zeta, z_vect=z_array)

        total_fields = electric_fields_to_lab(total_fields, system)
        modal_fields, qs_by_layer = compute_modal_fields_at_points(
            system,
            freq_hz,
            zeta,
            z_array,
            layer_indices,
            layer_slices,
            system_layer_indices=system_layer_indices,
        )
        channels_by_pol = {
            pol_idx: self._channels_for_pol(modal_fields, qs_by_layer, layer_indices, layer_slices, pol_idx)
            for pol_idx in pol_indices
        }

        return OpticalChannelSolution(
            frequency_cm1=float(frequency_cm1),
            zeta=zeta,
            z_points=z_array,
            total_fields=total_fields,
            modal_fields=modal_fields,
            qs_by_layer=qs_by_layer,
            channels_by_pol=channels_by_pol,
        )

    def external_momentum_transfer(
        self,
        incident_system,
        incident_frequency_cm1: float,
        incident_angle_rad: float,
        detected_system,
        detected_frequency_cm1: float,
        detected_angle_rad: float,
        collection_side: str,
    ) -> np.ndarray:
        """Return the macroscopic external momentum transfer ``q_ext = k_inc - k_det``.

        The returned vector is normalized to the incident vacuum wavevector.
        Detector angles specify the outgoing ray in laboratory coordinates.
        The Stokes/laser ratio is retained, including in forward scattering.

        This helper follows the current GTM boundary convention and assumes the
        external half-space wavevector is determined by ``epsilon[0, 0]``.  It is
        therefore exact for isotropic superstrate/substrate media and a
        documented approximation for anisotropic external media.
        """
        k_inc = self._external_wavevector(incident_system, incident_frequency_cm1, incident_angle_rad, z_sign=1.0)
        if collection_side == "superstrate":
            k_det = self._external_wavevector(detected_system, detected_frequency_cm1, detected_angle_rad, z_sign=-1.0)
        elif collection_side == "substrate":
            k_det = self._external_wavevector(detected_system, detected_frequency_cm1, detected_angle_rad, z_sign=1.0)
        else:
            msg = f"collection_side must be 'superstrate' or 'substrate', got {collection_side!r}"
            raise ValueError(msg)
        return k_inc - (detected_frequency_cm1 / incident_frequency_cm1) * k_det

    def _channels_for_pol(self, modal_fields, qs_by_layer, layer_indices, layer_slices, pol_idx):
        channels_by_layer = {}
        for layer_index, layer_slice in zip(layer_indices, layer_slices, strict=True):
            if layer_index not in qs_by_layer:
                continue
            subspaces = self._subspace_resolver.classify_optical_subspaces(qs_by_layer[layer_index])
            layer_channels = []
            for subspace in subspaces:
                field = np.sum(modal_fields[list(subspace.indices), pol_idx, :, :], axis=0)
                channel_field = np.zeros((3, modal_fields.shape[-1]), dtype=np.complex128)
                channel_field[:, layer_slice] = field[:, layer_slice]
                layer_channels.append(
                    OpticalChannel(
                        layer_index=layer_index,
                        qz=subspace.kz_reference,
                        field=channel_field,
                        mode_indices=subspace.indices,
                    )
                )
            channels_by_layer[layer_index] = sorted(layer_channels, key=lambda item: np.real(item.qz), reverse=True)
        return channels_by_layer

    @staticmethod
    def _external_wavevector(system, frequency_cm1, angle_rad, z_sign):
        freq_hz = float(frequency_cm1) * speed_light_si * 1.0e2
        system.initialize_sys(freq_hz)
        eps_xx = system.superstrate.epsilon[0, 0]
        zeta = np.sin(float(angle_rad)) * np.sqrt(eps_xx)
        qz_sq = eps_xx - zeta**2
        qz = np.sqrt(qz_sq)
        return np.array([np.real(zeta), 0.0, float(z_sign) * np.real(qz)], dtype=float)


def compute_modal_fields_at_points(
    system,
    freq_hz: float,
    zeta: complex,
    z_array: np.ndarray,
    layer_indices: Sequence[int],
    layer_slices: Sequence[slice],
    system_layer_indices: Sequence[int] | None = None,
) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """Compute per-Berreman-mode electric fields at selected z-points.

    ``layer_indices`` are the physical layer ids used in returned dictionaries.
    ``system_layer_indices`` are the corresponding indices in ``system.layers``.
    Passing distinct values is required for reversed systems so physical layer
    ids remain stable while the optical stack order changes.
    """
    system_indices = layer_indices if system_layer_indices is None else system_layer_indices
    system.calculate_GammaStar(freq_hz, zeta)
    modal_amps, _zn = system.calculate_modal_amplitudes(freq_hz, zeta)
    boundaries = system.get_layers_boundaries()
    modal_fields = np.zeros((4, 2, 3, len(z_array)), dtype=np.complex128)
    qs_by_layer = {}

    for physical_layer_index, system_layer_index, layer_slice in zip(
        layer_indices,
        system_indices,
        layer_slices,
        strict=True,
    ):
        if system_layer_index not in modal_amps or system_layer_index >= len(system.layers):
            continue
        layer = system.layers[system_layer_index]
        amps = modal_amps[system_layer_index]
        qs_by_layer[physical_layer_index] = layer.qs.copy()
        z_reference = boundaries[system_layer_index + 2]

        for offset, z_value in enumerate(z_array[layer_slice]):
            z_pos = layer_slice.start + offset
            propagation = np.array(
                [
                    np.exp(-layer.propagation_exponents[mode_idx] * (z_value - z_reference) / layer.thick)
                    for mode_idx in range(4)
                ],
                dtype=np.complex128,
            )
            for mode_idx in range(4):
                modal_fields[mode_idx, 0, :, z_pos] = propagation[mode_idx] * amps[mode_idx] * layer.gamma[mode_idx]
                modal_fields[mode_idx, 1, :, z_pos] = (
                    propagation[mode_idx] * amps[mode_idx + 4] * layer.gamma[mode_idx]
                )

    return electric_fields_to_lab(modal_fields, system, component_axis=2), qs_by_layer
