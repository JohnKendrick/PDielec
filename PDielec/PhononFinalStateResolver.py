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
import logging
import warnings
from collections.abc import Hashable, Iterable, Sequence

import numpy as np

logger = logging.getLogger(__name__)

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


@dataclasses.dataclass(frozen=True)
class NacBranchData:
    """NAC branch data mapped back to one original TO mode."""

    branch_index: int
    frequency_cm1: float
    raman_tensor: np.ndarray
    linewidth_cm1: float | None = None
    selected: bool = True
    overlap: float | None = None


@dataclasses.dataclass(frozen=True)
class PhononSubspace:
    """A nearly degenerate phonon branch subspace."""

    indices: tuple[int, ...]
    frequency_cm1: float
    raman_tensor: np.ndarray


class PhononFinalStateResolver:
    """Classify modal-pair phonon final states without optical-solver dependencies."""

    def __init__(
        self,
        coherence_regime: str = COHERENT_FILM,
        q_ext: Sequence[float] | None = None,
        angular_tolerance_deg: float = 5.0,
        matching_tolerance_deg: float | None = None,
        q_zero_tol: float = 1.0e-3,
        eps_optical: float = 1.0e-4,
        eps_phonon: float = 1.0e-3,
        min_nac_overlap: float = 0.5,
        nonpolar_strength_threshold: float = 1.0e-4,
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
        self.matching_tolerance_deg = (
            self.angular_tolerance_deg
            if matching_tolerance_deg is None
            else float(matching_tolerance_deg)
        )
        self.q_zero_tol = float(q_zero_tol)
        self.eps_optical = float(eps_optical)
        self.eps_phonon = float(eps_phonon)
        self.min_nac_overlap = float(min_nac_overlap)
        self.nonpolar_strength_threshold = float(nonpolar_strength_threshold)
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

    def is_polar_mode(self, oscillator_strength: float, max_oscillator_strength: float) -> bool:
        """Return whether a mode is polar according to the configured strength threshold."""
        max_strength = max(float(max_oscillator_strength), 1.0e-30)
        return float(oscillator_strength) / max_strength >= self.nonpolar_strength_threshold

    def map_nac_branch(
        self,
        mode_idx: int,
        to_eigenvectors: np.ndarray,
        nac_eigenvectors: np.ndarray,
        nac_frequencies_cm1: Sequence[float],
        nac_raman_tensors: Sequence[np.ndarray],
        to_frequencies_cm1: Sequence[float] | None = None,
        nac_linewidths_cm1: Sequence[float] | None = None,
        nac_selected: Sequence[bool] | None = None,
        dominant_to_by_nac: Sequence[int] | None = None,
    ) -> NacBranchData:
        """Map an original TO mode to the matching NAC branch by subspace overlap."""
        branch_index, overlap = self.nac_branch_index(
            mode_idx,
            to_eigenvectors,
            nac_eigenvectors,
            nac_frequencies_cm1=nac_frequencies_cm1,
            fallback_frequency_cm1=None if to_frequencies_cm1 is None else float(to_frequencies_cm1[mode_idx]),
            dominant_to_by_nac=dominant_to_by_nac,
        )
        linewidth = None if nac_linewidths_cm1 is None else float(nac_linewidths_cm1[branch_index])
        selected = True if nac_selected is None else bool(nac_selected[branch_index])
        return NacBranchData(
            branch_index=branch_index,
            frequency_cm1=float(nac_frequencies_cm1[branch_index]),
            raman_tensor=np.asarray(nac_raman_tensors[branch_index], dtype=complex),
            linewidth_cm1=linewidth,
            selected=selected,
            overlap=overlap,
        )

    def nac_branch_index(
        self,
        mode_idx: int,
        to_eigenvectors: np.ndarray,
        nac_eigenvectors: np.ndarray,
        nac_frequencies_cm1: Sequence[float] | None = None,
        fallback_frequency_cm1: float | None = None,
        dominant_to_by_nac: Sequence[int] | None = None,
    ) -> tuple[int, float]:
        """Return the NAC branch index associated with an original TO mode."""
        to_modes = self._modes_as_columns(to_eigenvectors)
        nac_modes = self._modes_as_columns(nac_eigenvectors)
        to_vec = to_modes[:, int(mode_idx)]
        to_norm = np.linalg.norm(to_vec)
        nac_norms = np.linalg.norm(nac_modes, axis=0)
        overlaps = np.abs(np.conj(nac_modes).T @ to_vec) ** 2
        overlaps /= np.maximum(nac_norms**2 * to_norm**2, 1.0e-30)
        branch_index = self._dominant_map_branch(mode_idx, dominant_to_by_nac)
        if branch_index is None:
            branch_index = int(np.argmax(overlaps))

        if overlaps[branch_index] < self.min_nac_overlap:
            branch_index = self._fallback_nac_branch(
                mode_idx,
                branch_index,
                float(overlaps[branch_index]),
                nac_frequencies_cm1,
                fallback_frequency_cm1,
            )
        return branch_index, float(overlaps[branch_index])

    @staticmethod
    def _dominant_map_branch(mode_idx: int, dominant_to_by_nac: Sequence[int] | None) -> int | None:
        if dominant_to_by_nac is None:
            return None
        candidates = np.where(np.asarray(dominant_to_by_nac, dtype=int) == int(mode_idx))[0]
        if len(candidates) == 0:
            return None
        return int(candidates[0])

    def _fallback_nac_branch(
        self,
        mode_idx: int,
        branch_index: int,
        overlap: float,
        nac_frequencies_cm1: Sequence[float] | None,
        fallback_frequency_cm1: float | None,
    ) -> int:
        msg = (
            f"NAC branch overlap {overlap:.3f} for TO mode {mode_idx} is below "
            f"threshold {self.min_nac_overlap:.3f}"
        )
        if nac_frequencies_cm1 is None or fallback_frequency_cm1 is None:
            warnings.warn(f"{msg}; retaining branch {branch_index} because no frequency fallback was supplied", stacklevel=2)
            logger.warning("%s; retaining branch %s because no frequency fallback was supplied", msg, branch_index)
            return branch_index

        frequencies = np.asarray(nac_frequencies_cm1, dtype=float)
        fallback_index = int(np.argmin(np.abs(frequencies - float(fallback_frequency_cm1))))
        warnings.warn(
            f"{msg}; falling back from NAC branch {branch_index} to closest-frequency branch {fallback_index}",
            stacklevel=2,
        )
        logger.warning(
            "%s; falling back from NAC branch %s to closest-frequency branch %s",
            msg,
            branch_index,
            fallback_index,
        )
        return fallback_index

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

        if self.coherence_regime == BULK_PHASE_MATCHED:
            if not self._matches_external_q(q_vec, q_external, reference_k):
                return DiscardedInternalComponent("phonon q does not match external momentum transfer", q_vec)
            return self._resolved(
                mode_idx,
                q_vec,
                True,
                frequency_cm1,
                raman_tensor,
                linewidth_cm1,
                selected,
                phase_weight,
                final_state_q=q_external,
            )

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
        final_state_q: np.ndarray | None = None,
    ) -> ResolvedFinalState:
        frequency_key = None if frequency_cm1 is None else round(float(frequency_cm1), 8)
        key_q = q_vec if final_state_q is None else final_state_q
        key = self.final_state_key(mode_idx, key_q, is_polar=is_polar, frequency_key=frequency_key)
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
        return angle <= self.matching_tolerance_deg

    def classify_optical_subspaces(
        self,
        kz_values: Sequence[complex],
        modes: np.ndarray | None = None,
    ) -> list[OpticalSubspace]:
        """Group nearly degenerate optical kz values and optionally attach projectors."""
        kz_array = np.asarray(kz_values, dtype=complex)
        neighbours = self._connected_neighbours(len(kz_array))
        for idx_a in range(len(kz_array)):
            for idx_b in range(idx_a + 1, len(kz_array)):
                if self._kz_degenerate(kz_array[idx_a], kz_array[idx_b]):
                    neighbours[idx_a].add(idx_b)
                    neighbours[idx_b].add(idx_a)

        groups: list[OpticalSubspace] = []
        for group in self._connected_components(neighbours):
            projector = None if modes is None else self.subspace_projector(np.asarray(modes)[:, group])
            groups.append(OpticalSubspace(tuple(group), kz_array[group].mean(), projector))
        return groups

    def classify_phonon_subspaces(
        self,
        frequencies_cm1: Sequence[float],
        raman_tensors: Sequence[np.ndarray],
    ) -> list[PhononSubspace]:
        """Group nearly degenerate phonon branches and sum their Raman tensors."""
        frequencies = np.asarray(frequencies_cm1, dtype=float)
        neighbours = self._connected_neighbours(len(frequencies))
        for idx_a in range(len(frequencies)):
            for idx_b in range(idx_a + 1, len(frequencies)):
                if self._phonon_degenerate(frequencies[idx_a], frequencies[idx_b]):
                    neighbours[idx_a].add(idx_b)
                    neighbours[idx_b].add(idx_a)

        subspaces: list[PhononSubspace] = []
        for group in self._connected_components(neighbours):
            tensor = sum(np.asarray(raman_tensors[idx], dtype=complex) for idx in group)
            subspaces.append(PhononSubspace(tuple(group), float(np.mean(frequencies[group])), tensor))
        return subspaces

    def _kz_degenerate(self, kz_a: complex, kz_b: complex) -> bool:
        mean = max(0.5 * (abs(kz_a) + abs(kz_b)), 1.0e-30)
        return abs(kz_a - kz_b) / mean < self.eps_optical

    def _phonon_degenerate(self, frequency_a: float, frequency_b: float) -> bool:
        mean = max(0.5 * (abs(frequency_a) + abs(frequency_b)), 1.0e-30)
        return abs(frequency_a - frequency_b) / mean < self.eps_phonon

    @staticmethod
    def _connected_neighbours(size: int) -> dict[int, set[int]]:
        return {idx: set() for idx in range(size)}

    @staticmethod
    def _connected_components(neighbours: dict[int, set[int]]) -> list[list[int]]:
        unused = set(neighbours)
        components = []
        while unused:
            stack = [min(unused)]
            group = []
            while stack:
                idx = stack.pop()
                if idx not in unused:
                    continue
                unused.remove(idx)
                group.append(idx)
                stack.extend(sorted(neighbours[idx] & unused, reverse=True))
            components.append(sorted(group))
        return components

    @staticmethod
    def subspace_projector(modes: np.ndarray) -> np.ndarray:
        """Return the Euclidean projector onto the column span of ``modes``."""
        mode_array = np.asarray(modes, dtype=complex)
        q_basis, _ = np.linalg.qr(mode_array)
        rank = np.linalg.matrix_rank(mode_array)
        q_basis = q_basis[:, :rank]
        return q_basis @ np.conj(q_basis).T

    @staticmethod
    def _modes_as_columns(modes: np.ndarray) -> np.ndarray:
        mode_array = np.asarray(modes, dtype=complex)
        if mode_array.ndim != 2:
            msg = f"Expected a 2D mode array, got shape {mode_array.shape}"
            raise ValueError(msg)
        return mode_array

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
        detector_modes: np.ndarray | None = None,
    ) -> complex:
        """Return ``d * (P_S E_S)^T R P_L E_L`` for optical subspaces."""
        projected_incident = cls.project_field_onto_subspace(incident_modes, incident_field)
        projected_detector = np.asarray(detector_field, dtype=complex)
        if detector_modes is not None:
            projected_detector = cls.project_field_onto_subspace(detector_modes, projected_detector)
        return detector_amplitude * (projected_detector @ raman_tensor @ projected_incident)

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
