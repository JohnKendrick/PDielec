#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Small coherent-amplitude accumulator for final-state resolver tests."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Hashable


class RamanAmplitudeAccumulator:
    """Accumulate complex Raman amplitudes by physical final-state key."""

    def __init__(self) -> None:
        """Initialise an empty accumulator."""
        self._amplitudes = defaultdict(complex)

    def add(self, key: Hashable, amplitude: complex) -> None:
        """Add one coherent contribution to ``key``."""
        self._amplitudes[key] += complex(amplitude)

    def amplitude(self, key: Hashable) -> complex:
        """Return the accumulated amplitude for one key."""
        return self._amplitudes[key]

    def amplitudes(self) -> dict[Hashable, complex]:
        """Return a plain dictionary of accumulated amplitudes."""
        return dict(self._amplitudes)

    def intensities(self) -> dict[Hashable, float]:
        """Return ``|A_key|^2`` for each resolved final-state key."""
        return {key: abs(amplitude) ** 2 for key, amplitude in self._amplitudes.items()}

    def total_intensity(self) -> float:
        """Return the incoherent sum over resolved final-state intensities."""
        return sum(self.intensities().values())
