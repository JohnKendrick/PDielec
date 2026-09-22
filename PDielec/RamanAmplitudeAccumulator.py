#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Coherent Raman amplitude accumulator for final-state grouped records."""

from __future__ import annotations

import dataclasses
from collections import defaultdict
from collections.abc import Hashable, Mapping

from PDielec.PhononFinalStateResolver import (
    CoherentExternalChannelContribution,
    DiscardedInternalComponent,
    ResolvedFinalState,
)


@dataclasses.dataclass(frozen=True)
class RamanAmplitudeRecord:
    """One classified Raman amplitude or local-intensity contribution."""

    final_state_key: Hashable
    amplitude: complex = 0.0 + 0.0j
    coherence_group: Hashable | None = None
    local_intensity: float = 0.0


class RamanAmplitudeAccumulator:
    """Accumulate Raman signal with coherent and incoherent boundaries explicit.

    Contributions with the same ``final_state_key`` and ``coherence_group`` are
    summed as complex amplitudes before squaring.  Different final-state keys, or
    different coherence groups for the same final state, are summed incoherently.
    """

    def __init__(self) -> None:
        """Initialise an empty accumulator."""
        self._coherent_amplitudes = defaultdict(complex)
        self._local_intensities = defaultdict(float)

    def add(self, key: Hashable, amplitude: complex) -> None:
        """Add one coherent contribution using the legacy key-as-group API."""
        self.add_amplitude(key, amplitude, coherence_group=key)

    def add_amplitude(
        self,
        final_state_key: Hashable,
        amplitude: complex,
        coherence_group: Hashable | None = None,
    ) -> None:
        """Add an amplitude to a final state and coherent accumulation group."""
        group = final_state_key if coherence_group is None else coherence_group
        self._coherent_amplitudes[(final_state_key, group)] += complex(amplitude)

    def add_incoherent_intensity(self, final_state_key: Hashable, intensity: float) -> None:
        """Add a local intensity that has already crossed an incoherent boundary."""
        self._local_intensities[final_state_key] += float(intensity)

    def add_record(self, record: RamanAmplitudeRecord) -> None:
        """Add a pre-classified amplitude record."""
        if record.local_intensity:
            self.add_incoherent_intensity(record.final_state_key, record.local_intensity)
        if record.amplitude:
            self.add_amplitude(record.final_state_key, record.amplitude, record.coherence_group)

    def add_classified(
        self,
        classification: ResolvedFinalState | CoherentExternalChannelContribution | DiscardedInternalComponent,
        amplitude: complex,
        coherence_group: Hashable | None = None,
    ) -> None:
        """Add an amplitude using a resolver classification object."""
        if isinstance(classification, DiscardedInternalComponent):
            return
        if isinstance(classification, ResolvedFinalState):
            key = classification.final_state.key
        elif isinstance(classification, CoherentExternalChannelContribution):
            key = classification.key
        else:
            msg = f"Unsupported classification type: {type(classification)!r}"
            raise TypeError(msg)
        self.add_amplitude(key, amplitude * classification.phase_weight, coherence_group)

    def amplitude(self, key: Hashable, coherence_group: Hashable | None = None) -> complex:
        """Return the accumulated amplitude for a final state/group pair."""
        group = key if coherence_group is None else coherence_group
        return self._coherent_amplitudes[(key, group)]

    def amplitudes(self) -> dict[Hashable, complex]:
        """Return amplitudes for legacy callers only when each final state has one group."""
        self._raise_if_multiple_groups("amplitudes")
        return {key: amplitude for (key, _group), amplitude in self._coherent_amplitudes.items()}

    def grouped_amplitudes(self) -> dict[tuple[Hashable, Hashable], complex]:
        """Return amplitudes keyed by ``(final_state_key, coherence_group)``."""
        return dict(self._coherent_amplitudes)

    def coherent_final_state_amplitudes(self) -> dict[Hashable, complex]:
        """Return diagnostic amplitudes after coherently collapsing all groups.

        This intentionally ignores coherence-group boundaries and is not a
        physical intensity path when a final state has multiple groups.
        """
        amplitudes = defaultdict(complex)
        for (key, _group), amplitude in self._coherent_amplitudes.items():
            amplitudes[key] += amplitude
        return dict(amplitudes)

    def final_state_amplitudes(self) -> dict[Hashable, complex]:
        """Return amplitudes only when each final state has one coherent group."""
        self._raise_if_multiple_groups("final_state_amplitudes")
        return self.coherent_final_state_amplitudes()

    def intensities(self) -> dict[Hashable, float]:
        """Return final-state intensities after coherent-group squaring."""
        intensities = defaultdict(float, self._local_intensities)
        for (key, _group), amplitude in self._coherent_amplitudes.items():
            intensities[key] += abs(amplitude) ** 2
        return dict(intensities)

    def total_intensity(self) -> float:
        """Return the incoherent sum over all resolved final-state intensities."""
        return sum(self.intensities().values())

    @classmethod
    def from_records(cls, records: list[RamanAmplitudeRecord]) -> RamanAmplitudeAccumulator:
        """Build an accumulator from amplitude records."""
        accumulator = cls()
        for record in records:
            accumulator.add_record(record)
        return accumulator

    def update(self, records: Mapping[Hashable, complex]) -> None:
        """Add legacy ``key -> amplitude`` contributions."""
        for key, amplitude in records.items():
            self.add(key, amplitude)

    def _raise_if_multiple_groups(self, caller: str) -> None:
        groups_by_key = defaultdict(set)
        for key, group in self._coherent_amplitudes:
            groups_by_key[key].add(group)
        multi_group_keys = [key for key, groups in groups_by_key.items() if len(groups) > 1]
        if multi_group_keys:
            msg = (
                f"{caller}() is only valid when each final state has one coherence group; "
                f"use grouped_amplitudes() or intensities() for multi-group keys {multi_group_keys!r}"
            )
            raise ValueError(msg)
