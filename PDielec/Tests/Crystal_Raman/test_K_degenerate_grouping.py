#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Tests for group_degenerate_modes in AzimuthalSweepWindow.

The key physical property under test: for a degenerate subspace the
individual per-eigenmode intensities are basis-dependent (they change when
the eigenvectors are rotated within the subspace), but the sum of
intensities over the subspace must be invariant to that rotation.
"""
import numpy as np
import pytest

from PDielec.GUI.AzimuthalSweepWindow import group_degenerate_modes


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_sweep(n_angles, frequencies, raw_ints):
    """Return intensities_2d of shape (n_angles, M) by tiling raw_ints."""
    raw = np.asarray(raw_ints, dtype=float)
    return np.tile(raw, (n_angles, 1))


# ---------------------------------------------------------------------------
# Basic grouping behaviour
# ---------------------------------------------------------------------------

class TestGroupDegenerateModes:
    """Unit tests for group_degenerate_modes."""

    def test_no_grouping_all_distinct(self):
        """Modes separated > tol_cm1 must each remain in their own group."""
        freqs = np.array([100.0, 200.0, 300.0])
        ints = np.eye(3, dtype=float)  # 3 angles × 3 modes
        gf, gi, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert len(gf) == 3
        assert gi.shape == (3, 3)
        assert len(labels) == 3

    def test_all_merged_when_within_tol(self):
        """All modes within tol_cm1 collapse to a single group."""
        freqs = np.array([100.0, 100.5, 101.0])
        ints = _make_sweep(5, freqs, [1.0, 2.0, 3.0])
        gf, gi, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert len(gf) == 1
        np.testing.assert_allclose(gi[:, 0], 6.0)
        assert "3 modes" in labels[0]

    def test_partial_grouping(self):
        """Only neighbouring modes within tol_cm1 are merged."""
        # Modes 0+1 are close; mode 2 is far away
        freqs = np.array([100.0, 100.8, 200.0])
        ints = _make_sweep(4, freqs, [1.0, 1.0, 5.0])
        gf, gi, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert len(gf) == 2
        np.testing.assert_allclose(gi[:, 0], 2.0)   # merged pair
        np.testing.assert_allclose(gi[:, 1], 5.0)   # standalone

    def test_empty_input(self):
        """Empty frequency array must return empty outputs without error."""
        freqs = np.array([])
        ints = np.zeros((10, 0))
        gf, gi, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert len(gf) == 0
        assert gi.shape == (10, 0)
        assert labels == []

    def test_single_mode(self):
        """Single mode must pass through unchanged."""
        freqs = np.array([123.4])
        ints = np.ones((5, 1))
        gf, gi, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        np.testing.assert_allclose(gf, [123.4])
        np.testing.assert_allclose(gi, np.ones((5, 1)))
        assert "123.4" in labels[0]

    def test_label_format_single(self):
        """Single-mode group label should contain just the frequency."""
        freqs = np.array([456.7])
        ints = np.ones((2, 1))
        _, _, labels = group_degenerate_modes(freqs, ints)
        assert "456.7" in labels[0]
        assert "modes" not in labels[0]

    def test_label_format_group(self):
        """Multi-mode group label must include both boundary frequencies."""
        freqs = np.array([100.0, 101.0])
        ints = np.ones((2, 2))
        _, _, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert "100.0" in labels[0]
        assert "101.0" in labels[0]
        assert "2 modes" in labels[0]

    def test_output_shape_is_angles_by_groups(self):
        """grouped_ints must have shape (N_angles, N_groups)."""
        n_angles = 37
        freqs = np.array([100.0, 100.5, 200.0, 201.0, 300.0])
        ints = np.random.default_rng(0).random((n_angles, len(freqs)))
        gf, gi, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        # Groups: {100.0, 100.5}, {200.0, 201.0}, {300.0} → 3 groups
        assert gi.shape == (n_angles, 3)

    def test_unsorted_input_is_handled(self):
        """Input frequencies need not be sorted; output must still be correct."""
        # Give frequencies in reverse order
        freqs = np.array([200.0, 100.0, 100.5])
        # ints columns correspond to: 200.0, 100.0, 100.5
        ints = np.array([[10.0, 1.0, 2.0],
                         [20.0, 3.0, 4.0]])
        gf, gi, labels = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert len(gf) == 2
        # Group at ~100.25 should sum columns for 100.0 and 100.5
        g_low = gi[:, 0]   # lower mean freq group
        g_high = gi[:, 1]  # higher mean freq group
        np.testing.assert_allclose(g_low, [3.0, 7.0])   # 1+2 and 3+4
        np.testing.assert_allclose(g_high, [10.0, 20.0])


# ---------------------------------------------------------------------------
# Physical invariance: degenerate-subspace rotation
# ---------------------------------------------------------------------------

class TestDegenerateSubspaceInvariance:
    """The grouped intensity must be invariant to in-subspace basis rotation.

    Physical motivation: two eigenvectors that span a doubly-degenerate
    phonon subspace are not uniquely defined — any unitary rotation within
    the subspace gives an equally valid pair.  Individual per-mode
    intensities therefore change with this choice, but their sum must not.
    """

    @pytest.fixture
    def degenerate_pair(self):
        """Return (freqs, base_ints) for a doubly-degenerate pair.

        base_ints has shape (N_angles, 2) with deliberately asymmetric
        column intensities so that rotation of the basis visibly changes
        individual values.
        """
        rng = np.random.default_rng(42)
        n_angles = 36
        # Two columns with different intensities at each angle
        col0 = rng.random(n_angles) * 2.0
        col1 = rng.random(n_angles) * 3.0
        ints = np.column_stack([col0, col1])
        freqs = np.array([380.0, 380.5])  # within default 2 cm⁻¹ tol
        return freqs, ints

    def test_grouped_sum_equals_raw_sum(self, degenerate_pair):
        """Grouped intensity must equal the direct column sum."""
        freqs, ints = degenerate_pair
        gf, gi, _ = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert gi.shape[1] == 1, "Should merge into a single group"
        expected = ints.sum(axis=1)
        np.testing.assert_allclose(gi[:, 0], expected)

    def test_invariant_under_basis_rotation(self, degenerate_pair):
        """Summed intensity must be invariant to in-subspace amplitude rotation.

        Intensities are |amplitude|².  When the pair of eigenvector amplitudes
        is rotated within the degenerate subspace, the individual per-mode
        intensities change but their sum is preserved by Parseval's theorem.
        """
        freqs, ints = degenerate_pair
        # Recover unsigned amplitudes from the intensity data
        amps = np.sqrt(ints)            # shape (N_angles, 2), non-negative

        theta = np.radians(37.0)
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta),  np.cos(theta)]])
        # Rotate the amplitude pair; re-compute intensities as squares
        amps_rotated = amps @ R.T
        ints_rotated = amps_rotated ** 2

        gf1, gi1, _ = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        gf2, gi2, _ = group_degenerate_modes(freqs, ints_rotated, tol_cm1=2.0)

        # Individual columns must differ after rotation (sanity check)
        assert not np.allclose(ints[:, 0], ints_rotated[:, 0]), (
            "Sanity check: rotation should change individual per-mode intensities"
        )
        # The grouped sums must be equal (sum of squares preserved by rotation)
        np.testing.assert_allclose(gi1, gi2, rtol=1e-12)

    def test_non_degenerate_mode_unaffected(self, degenerate_pair):
        """A well-separated third mode must not be merged with the pair."""
        freqs_pair, ints_pair = degenerate_pair
        freqs = np.append(freqs_pair, 500.0)
        rng = np.random.default_rng(7)
        col2 = rng.random(ints_pair.shape[0]) * 4.0
        ints = np.column_stack([ints_pair, col2])

        gf, gi, _ = group_degenerate_modes(freqs, ints, tol_cm1=2.0)
        assert gi.shape[1] == 2
        np.testing.assert_allclose(gi[:, 1], col2)

    def test_zero_tol_never_merges(self):
        """With tol_cm1=0 no two distinct-frequency modes should merge."""
        freqs = np.array([100.0, 100.0001, 200.0])
        ints = np.ones((3, 3))
        gf, gi, _ = group_degenerate_modes(freqs, ints, tol_cm1=0.0)
        # 100.0 and 100.0001 differ by 0.0001 > 0, so should be separate
        assert len(gf) == 3

    def test_zno_style_e_symmetry(self):
        """Mimic ZnO E-symmetry doublet: pair invariant, A1 unchanged.

        ZnO has E modes (doubly degenerate) and A1 modes (non-degenerate).
        The E-mode pair should collapse to a single group; the A1 modes
        should remain as separate groups.
        """
        n_angles = 72
        rng = np.random.default_rng(99)

        # Simulate ZnO-like frequencies: two E-mode pairs + two A1 modes
        freqs = np.array([380.0, 380.3,   # E(TO) doublet
                          574.0, 574.2,   # E(LO) doublet
                          410.0,          # A1(TO)
                          579.0])         # A1(LO)

        raw = rng.random((n_angles, 6))
        gf, gi, labels = group_degenerate_modes(freqs, raw, tol_cm1=2.0)

        # After sorting by freq: 380.0, 380.3, 410.0, 574.0, 574.2, 579.0
        # Groups: {380.0, 380.3}, {410.0}, {574.0, 574.2}, {579.0} → 4 groups
        assert gi.shape[1] == 4, f"Expected 4 groups, got {gi.shape[1]}"

        # Confirm the E-doublet groups sum correctly
        sort_idx = np.argsort(freqs)
        s = raw[:, sort_idx]
        # Group 0: cols 0+1 of sorted → E(TO) pair
        np.testing.assert_allclose(gi[:, 0], s[:, 0] + s[:, 1])
        # Group 1: col 2 of sorted → A1(TO)
        np.testing.assert_allclose(gi[:, 1], s[:, 2])
