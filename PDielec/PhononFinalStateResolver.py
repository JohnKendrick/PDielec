#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Final-state classification helpers for modal-pair Raman calculations.

This module is intentionally independent of GTM, the GUI, and DFT readers.  It
contains only the phonon-side decisions needed by the modal-pairs final-state
refactor: optical subspace grouping, near-zero-q classification, q-direction
keys, and finite-thickness phase-matching weights.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Hashable, Iterable, Sequence

import numpy as np

COHERENT_FILM = "coherent_film"
INCOHERENT_DEPTH = "incoherent_depth"
BULK_PHASE_MATCHED = "bulk_phase_matched"
COHERENCE_REGIMES = (COHERENT_FILM, INCOHERENT_DEPTH, BULK_PHASE_MATCHED)

RESOLVED_FINAL_STATE = "resolved_final_state"
COHERENT_EXTERNAL_CHANNEL = "coherent_external_channel"
DISCARDED_INTERNAL_COMPONENT = "discarded_internal_component"


@dataclasses.dataclass(frozen=True)
class FinalState:
    """Resolved phonon final-state data."""

    key: Hashable
    frequency_cm1: float | None = None
    raman_tensor: np.ndarray | None = None
    linewidth_cm1: float | None = None
    selected: bool = True


@dataclasses.dataclass(frozen=True)
class ResolvedFinalState:
    """Classification for a physically resolved phonon final state."""

    final_state: FinalState
    q_ph: np.ndarray
    phase_weight: complex = 1.0 + 0.0j
    kind: str = RESOLVED_FINAL_STATE


@dataclasses.dataclass(frozen=True)
class CoherentExternalChannelContribution:
    """Classification for a coherent component of one external-channel amplitude."""

    key: Hashable
    q_ph: np.ndarray
    raman_tensor: np.ndarray | None = None
    phase_weight: complex = 1.0 + 0.0j
    reason: str = ""
    kind: str = COHERENT_EXTERNAL_CHANNEL


@dataclasses.dataclass(frozen=True)
class DiscardedInternalComponent:
    """Classification for an internal component that is not retained."""

    reason: str
    q_ph: np.ndarray | None = None
    phase_weight: complex = 0.0 + 0.0j
    kind: str = DISCARDED_INTERNAL_COMPONENT


@dataclasses.dataclass(frozen=True)
class OpticalSubspace:
    """A numerically degenerate optical propagation subspace."""

    indices: tuple[int, ...]
    kz_reference: complex
    projector: np.ndarray | None = None


class PhononFinalStateResolver:
    """Classify modal-pair phonon final states without optical-solver dependencies."""

    def __init__(
        self,
        coherence_regime: str = COHERENT_FILM,
        q_ext: Sequence[float] | None = None,
        angular_tolerance_deg: float = 5.0,
        q_zero_tol: float = 1.0e-3,
        eps_optical: float = 1.0e-4,
        layer_coherence_group: Hashable = "layer",
        detected_channel_key: Hashable = "detected",
    ) -> None:
        """Initialise resolver thresholds and external-channel identity."""
        if coherence_regime not in COHERENCE_REGIMES:
            msg = f"Unknown coherence_regime {coherence_regime!r}; expected one of {COHERENCE_REGIMES}"
            raise ValueError(msg)
        self.coherence_regime = coherence_regime
        self.q_ext = self._as_vector(q_ext) if q_ext is not None else np.zeros(3, dtype=float)
        self.angular_tolerance_deg = float(angular_tolerance_deg)
        self.q_zero_tol = float(q_zero_tol)
        self.eps_optical = float(eps_optical)
        self.layer_coherence_group = layer_coherence_group
        self.detected_channel_key = detected_channel_key

    @staticmethod
    def _as_vector(vector: Sequence[float] | np.ndarray) -> np.ndarray:
        arr = np.asarray(vector, dtype=float)
        if arr.shape != (3,):
            msg = f"Expected a 3-vector, got shape {arr.shape}"
            raise ValueError(msg)
        return arr

    @staticmethod
    def sinc(x: float | complex) -> complex:
        """Return unnormalised sinc, ``sin(x) / x`` with a continuous zero limit."""
        if abs(x) < 1.0e-14:
            return 1.0 + 0.0j
        return np.sin(x) / x

    @classmethod
    def phase_matching_factor(cls, delta_k: float, thickness: float) -> complex:
        """Return ``integral_0^L exp(i delta_k z) dz`` for one pair."""
        arg = delta_k * thickness / 2.0
        return thickness * np.exp(1j * arg) * cls.sinc(arg)

    @classmethod
    def normalised_phase_matching_factor(cls, delta_k: float, thickness: float) -> complex:
        """Return the finite-thickness factor divided by layer thickness."""
        return cls.phase_matching_factor(delta_k, thickness) / thickness

    @staticmethod
    def phonon_wavevector(k_l_physical: Sequence[float], k_s_physical: Sequence[float]) -> np.ndarray:
        """Return the physical phonon wavevector ``q_ph = k_L - k_S``."""
        return np.asarray(k_l_physical, dtype=float) - np.asarray(k_s_physical, dtype=float)

    def is_near_zero_q(self, q_ph: Sequence[float], reference_k: float | None = None) -> bool:
        """Return whether ``q_ph`` is below the near-zero threshold."""
        scale = 1.0 if reference_k is None else max(abs(reference_k), 1.0e-30)
        return np.linalg.norm(self._as_vector(q_ph)) <= self.q_zero_tol * scale

    def q_class_key(self, q_ph: Sequence[float]) -> tuple[str, int, int] | tuple[str]:
        """Return an angular q-direction key.

        The key stores rounded spherical angles using bins no finer than
        ``angular_tolerance_deg``.  The vector magnitude is deliberately absent.
        """
        q_vec = self._as_vector(q_ph)
        q_norm = np.linalg.norm(q_vec)
        if q_norm <= self.q_zero_tol:
            return ("q_zero",)
        q_hat = q_vec / q_norm
        theta = np.degrees(np.arccos(np.clip(q_hat[2], -1.0, 1.0)))
        phi = np.degrees(np.arctan2(q_hat[1], q_hat[0]))
        if phi < 0.0:
            phi += 360.0
        bin_size = max(self.angular_tolerance_deg, 1.0e-12)
        return ("q_dir", int(np.round(theta / bin_size)), int(np.round(phi / bin_size)))

    def final_state_key(
        self,
        mode_idx: int,
        q_ph: Sequence[float],
        is_polar: bool = True,
        frequency_key: Hashable | None = None,
    ) -> tuple[Hashable, Hashable, int, Hashable, Hashable | tuple[str]]:
        """Return the physical final-state key for a mode and q-vector."""
        q_key: Hashable | tuple[str] = self.q_class_key(q_ph) if is_polar else ("q_independent",)
        return (self.layer_coherence_group, self.detected_channel_key, int(mode_idx), frequency_key, q_key)

    def resolve_pair(  # noqa: PLR0911
        self,
        mode_idx: int,
        q_ph: Sequence[float],
        *,
        is_polar: bool = True,
        frequency_cm1: float | None = None,
        raman_tensor: np.ndarray | None = None,
        linewidth_cm1: float | None = None,
        selected: bool = True,
        q_ext: Sequence[float] | None = None,
        reference_k: float | None = None,
        delta_k: float | None = None,
        thickness: float | None = None,
    ) -> ResolvedFinalState | CoherentExternalChannelContribution | DiscardedInternalComponent:
        """Classify one optical pair's phonon contribution."""
        q_vec = self._as_vector(q_ph)
        q_external = self._as_vector(q_ext) if q_ext is not None else self.q_ext
        phase_weight = 1.0 + 0.0j
        if delta_k is not None and thickness is not None:
            phase_weight = self.phase_matching_factor(float(delta_k), float(thickness))

        if not selected:
            return DiscardedInternalComponent("mode is not selected", q_vec)

        if not is_polar:
            return self._resolved(mode_idx, q_vec, False, frequency_cm1, raman_tensor, linewidth_cm1, selected, phase_weight)

        q_near_zero = self.is_near_zero_q(q_vec, reference_k=reference_k)
        q_ext_near_zero = self.is_near_zero_q(q_external, reference_k=reference_k)

        if q_near_zero and not q_ext_near_zero:
            if self.coherence_regime == COHERENT_FILM:
                return CoherentExternalChannelContribution(
                    self._coherent_key(mode_idx),
                    q_vec,
                    raman_tensor=raman_tensor,
                    phase_weight=phase_weight,
                    reason="internal near-zero-q component of coherent external-channel amplitude",
                )
            return DiscardedInternalComponent("cross-pair phase-averages to zero", q_vec)

        if self.coherence_regime == COHERENT_FILM:
            return CoherentExternalChannelContribution(
                self._coherent_key(mode_idx),
                q_vec,
                raman_tensor=raman_tensor,
                phase_weight=phase_weight,
                reason="pair-resolved coherent finite-film contribution",
            )

        if self.coherence_regime == BULK_PHASE_MATCHED and not self._matches_external_q(q_vec, q_external, reference_k):
            return DiscardedInternalComponent("phonon q does not match external momentum transfer", q_vec)

        return self._resolved(mode_idx, q_vec, True, frequency_cm1, raman_tensor, linewidth_cm1, selected, phase_weight)

    def _resolved(
        self,
        mode_idx: int,
        q_vec: np.ndarray,
        is_polar: bool,
        frequency_cm1: float | None,
        raman_tensor: np.ndarray | None,
        linewidth_cm1: float | None,
        selected: bool,
        phase_weight: complex,
    ) -> ResolvedFinalState:
        frequency_key = None if frequency_cm1 is None else round(float(frequency_cm1), 8)
        key = self.final_state_key(mode_idx, q_vec, is_polar=is_polar, frequency_key=frequency_key)
        final_state = FinalState(key, frequency_cm1, raman_tensor, linewidth_cm1, selected)
        return ResolvedFinalState(final_state, q_vec, phase_weight=phase_weight)

    def _coherent_key(self, mode_idx: int) -> tuple[Hashable, Hashable, int, str]:
        return (self.layer_coherence_group, self.detected_channel_key, int(mode_idx), COHERENT_FILM)

    def _matches_external_q(
        self,
        q_ph: np.ndarray,
        q_ext: np.ndarray,
        reference_k: float | None = None,
    ) -> bool:
        if self.is_near_zero_q(q_ph, reference_k=reference_k) and self.is_near_zero_q(q_ext, reference_k=reference_k):
            return True
        if self.is_near_zero_q(q_ph, reference_k=reference_k) != self.is_near_zero_q(q_ext, reference_k=reference_k):
            return False
        q_norm = np.linalg.norm(q_ph)
        q_ext_norm = np.linalg.norm(q_ext)
        if q_norm <= 0.0 or q_ext_norm <= 0.0:
            return False
        cos_angle = np.dot(q_ph, q_ext) / (q_norm * q_ext_norm)
        angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
        return angle <= self.angular_tolerance_deg

    def classify_optical_subspaces(
        self,
        kz_values: Sequence[complex],
        modes: np.ndarray | None = None,
    ) -> list[OpticalSubspace]:
        """Group nearly degenerate optical kz values and optionally attach projectors."""
        kz_array = np.asarray(kz_values, dtype=complex)
        unused = set(range(len(kz_array)))
        groups: list[OpticalSubspace] = []
        while unused:
            seed = min(unused)
            group = [seed]
            unused.remove(seed)
            for idx in list(unused):
                if self._kz_degenerate(kz_array[seed], kz_array[idx]):
                    group.append(idx)
                    unused.remove(idx)
            projector = None if modes is None else self.subspace_projector(np.asarray(modes)[:, group])
            groups.append(OpticalSubspace(tuple(group), kz_array[group].mean(), projector))
        return groups

    def _kz_degenerate(self, kz_a: complex, kz_b: complex) -> bool:
        mean = max(0.5 * (abs(kz_a) + abs(kz_b)), 1.0e-30)
        return abs(kz_a - kz_b) / mean < self.eps_optical

    @staticmethod
    def subspace_projector(modes: np.ndarray) -> np.ndarray:
        """Return the Euclidean projector onto the column span of ``modes``."""
        mode_array = np.asarray(modes, dtype=complex)
        q_basis, _ = np.linalg.qr(mode_array)
        rank = np.linalg.matrix_rank(mode_array)
        q_basis = q_basis[:, :rank]
        return q_basis @ np.conj(q_basis).T

    @classmethod
    def project_field_onto_subspace(cls, modes: np.ndarray, field: Sequence[complex]) -> np.ndarray:
        """Project an electric field onto the span of an optical subspace."""
        return cls.subspace_projector(modes) @ np.asarray(field, dtype=complex)

    @classmethod
    def subspace_amplitude(
        cls,
        raman_tensor: np.ndarray,
        incident_modes: np.ndarray,
        incident_field: Sequence[complex],
        detector_field: Sequence[complex],
        detector_amplitude: complex = 1.0 + 0.0j,
    ) -> complex:
        """Return ``d * E_S^T R P E_L`` for a degenerate incident subspace."""
        projected_incident = cls.project_field_onto_subspace(incident_modes, incident_field)
        return detector_amplitude * (np.asarray(detector_field, dtype=complex) @ raman_tensor @ projected_incident)

    @staticmethod
    def coherent_amplitude(
        raman_tensor: np.ndarray,
        incident_amplitudes: Sequence[complex],
        incident_modes: np.ndarray,
        detector_amplitudes: Sequence[complex],
        detector_modes: np.ndarray,
    ) -> complex:
        """Return the full external-channel bilinear amplitude from modal components."""
        e_l = np.asarray(incident_modes, dtype=complex) @ np.asarray(incident_amplitudes, dtype=complex)
        e_s = np.asarray(detector_modes, dtype=complex) @ np.asarray(detector_amplitudes, dtype=complex)
        return e_s @ raman_tensor @ e_l

    @staticmethod
    def coherent_pair_amplitude(pairs: Iterable[tuple[float, complex]], thickness: float) -> complex:
        """Return ``sum a_ij F_ij(L)`` for synthetic finite-thickness pair tests."""
        return sum(a_ij * PhononFinalStateResolver.phase_matching_factor(delta_k, thickness) for delta_k, a_ij in pairs)
