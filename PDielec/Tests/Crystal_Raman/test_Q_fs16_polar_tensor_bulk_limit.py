#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""FS16 — Coherent pair-resolved amplitude recovers correct bulk LO/TO limits.

Design reference: modal_pairs_final_state_design.md, test FS16.

Physical requirement
--------------------
For a polar mode in a coherent finite film the Raman amplitude is

    A_m(L) = sum_{i,j} F_ij(L) * (E_S_rec,j^T @ R_m(q_ph(i,j)) @ E_L,i)

where the finite-thickness phase-matching factor is

    F_ij(L) = integral_0^L exp(i Delta_k_ij z) dz
            = L * exp(i Delta_k_ij L / 2) * sinc(Delta_k_ij L / 2)

and sinc here is the UNNORMALIZED sinc: sinc(x) = sin(x) / x, sinc(0) = 1.

The phase mismatch for pair (i, j) is

    Delta_k_ij = q_ph(i,j) - q_ext

where q_ph(i,j) = k_L,z(i) - k_S,z(j) is the phonon wavevector carried by
the pair and q_ext is the externally imposed momentum transfer.

Bulk limits
-----------
- If Delta_k_ij = 0: F_ij(L) / L = 1 for all L (pair survives in bulk).
- If Delta_k_ij != 0: F_ij(L) / L = sinc(Delta_k_ij L/2) * phase -> 0
  as L -> infinity, at rate ~1/(|Delta_k_ij| L) (algebraic, not exponential).

Therefore the normalized amplitude A_m(L) / L converges to the sum over
only the phase-matched pairs, selecting the correct tensor R_m(q_ph) for
each geometry automatically, without any geometry switch.

Pair model for the tests
------------------------
The tests use a simplified two-pair model (one matched, one mismatched)
and a four-pair backscattering model that follows the design document
exactly.  All tests are pure Python/numpy with no GTM or DFT data.

Two-pair model:
    matched pair:     Delta_k = 0,      raman_contribution = a_matched
    mismatched pair:  Delta_k = dk ≠ 0, raman_contribution = a_mismatch

A(L) / L = a_matched * 1 + a_mismatch * sinc(dk L/2) * exp(i dk L/2)

As L -> infty: A(L)/L -> a_matched.

For backscattering, the matched pair carries R_LO; for forward scattering
the matched pair carries R_TO.  The formula is identical; only the pair
assignments differ.

Four-pair backscattering model (follows design document step 1):
    LO co-prop fwd-fwd: Delta_k = 0     (phase-matched), R = R_LO
    LO co-prop bwd-bwd: Delta_k = -4k   (mismatched),    R = R_LO
    TO cross-pair fwd-bwd: Delta_k = -2k (mismatched),   R = R_TO
    TO cross-pair bwd-fwd: Delta_k = -2k (mismatched),   R = R_TO

Only the fwd-fwd pair survives in the bulk backscattering limit.
"""
import numpy as np
import numpy.testing as npt
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sinc(x):
    """Unnormalized sinc: sin(x)/x, with sinc(0) = 1.

    Note: numpy.sinc is normalised (sinc(x) = sin(pi x)/(pi x)).
    This function uses the unnormalized form needed for the phase-matching
    factor derivation.
    """
    if np.abs(x) < 1e-14:
        return 1.0 + 0j
    return np.sin(x) / x


def _F(delta_k, L):
    """Finite-thickness phase-matching factor for a single pair.

    F_ij(L) = L * exp(i Delta_k L/2) * sinc(Delta_k L/2)

    This equals integral_0^L exp(i Delta_k z) dz exactly.

    Parameters
    ----------
    delta_k : float
        Phase mismatch Delta_k_ij = q_ph(i,j) - q_ext.
    L : float
        Layer thickness in the same length units as 1/delta_k.

    Returns
    -------
    complex
        F_ij(L).
    """
    arg = delta_k * L / 2.0
    return L * np.exp(1j * arg) * _sinc(arg)


def _coherent_amplitude(pairs, L):
    """Total coherent amplitude A(L) = sum_{i,j} a_ij * F_ij(L).

    Parameters
    ----------
    pairs : list of (delta_k, raman_contribution)
        Each pair's phase mismatch and complex amplitude contribution
        (field overlap × tensor contraction, constant in z).
    L : float
        Layer thickness.

    Returns
    -------
    complex
        Total amplitude A(L).
    """
    return sum(a_ij * _F(dk, L) for dk, a_ij in pairs)


def _upper_bound_residual(pairs_mismatched, L):
    """Upper bound on |A_mismatch(L)/L|.

    For each mismatched pair (Delta_k != 0):
        |F_ij(L)/L| = |sinc(Delta_k L/2)| <= 2 / |Delta_k L|

    Returns the sum of per-pair bounds.
    """
    bound = 0.0
    for dk, a_ij in pairs_mismatched:
        if abs(dk) < 1e-14:
            continue
        bound += abs(a_ij) * 2.0 / (abs(dk) * L)
    return bound


# ---------------------------------------------------------------------------
# FS16a: Phase-matching factor properties
# ---------------------------------------------------------------------------

class TestFS16PhasMatchingFactor:
    """Direct tests of F_ij(L) and F_ij(L)/L.

    These tests establish the mathematical properties of the factor without
    any pair summation.
    """

    def test_F_equals_integral_exactly(self):
        """F_ij(L) equals the numerical integral of exp(i Delta_k z) over [0, L]."""
        rng = np.random.default_rng(7)
        for _ in range(20):
            delta_k = rng.uniform(-5.0, 5.0)
            L = rng.uniform(0.1, 10.0)
            F_analytic = _F(delta_k, L)

            # Numerical integration
            n = 10000
            zz = np.linspace(0, L, n)
            F_numeric = np.trapezoid(np.exp(1j * delta_k * zz), zz)

            assert abs(F_numeric - F_analytic) < 1e-5, (
                f"delta_k={delta_k:.3f}, L={L:.3f}: "
                f"|F_numeric - F_analytic| = {abs(F_numeric - F_analytic):.2e}"
            )

    def test_F_normalized_is_unity_for_zero_delta_k(self):
        """F(0, L) / L = 1 for all positive L."""
        for L in [0.01, 1.0, 100.0, 1e6]:
            F_over_L = _F(0.0, L) / L
            npt.assert_allclose(
                F_over_L.real, 1.0, atol=1e-14,
                err_msg=f"F(0, L={L})/L must be exactly 1"
            )
            npt.assert_allclose(
                F_over_L.imag, 0.0, atol=1e-14,
                err_msg=f"F(0, L={L})/L must be real"
            )

    @pytest.mark.parametrize("delta_k", [0.5, 1.0, 2.0, 5.0, -1.0, -3.0])
    def test_F_normalized_decays_for_mismatched_pairs(self, delta_k):
        """F(delta_k, L)/L -> 0 as L -> infinity for delta_k != 0."""
        L_values = [1e2, 1e4, 1e6]
        magnitudes = [abs(_F(delta_k, L)) / L for L in L_values]

        # Must be monotonically decreasing (up to oscillation envelope)
        # Tighter check: must be below the 2/(|dk| L) bound
        for L, mag in zip(L_values, magnitudes):
            bound = 2.0 / (abs(delta_k) * L)
            assert mag <= bound + 1e-14, (
                f"delta_k={delta_k}, L={L}: |F|/L = {mag:.2e} exceeds bound {bound:.2e}"
            )

        # At L=1e6, the factor should be very small
        assert magnitudes[-1] < 1e-5, (
            f"delta_k={delta_k}: |F|/L at L=1e6 = {magnitudes[-1]:.2e}, expected < 1e-5"
        )

    def test_F_continuity_at_zero_delta_k(self):
        """F(delta_k, L)/L is continuous at delta_k=0 (no branch cut)."""
        L = 1.0
        eps_values = [1e-3, 1e-6, 1e-9, 1e-12]
        for eps in eps_values:
            F_pos = _F(eps, L) / L
            F_neg = _F(-eps, L) / L
            F_zero = _F(0.0, L) / L
            npt.assert_allclose(
                F_pos.real, F_zero.real, atol=1e-6,
                err_msg=f"F not continuous from above at delta_k=+{eps}"
            )
            npt.assert_allclose(
                F_neg.real, F_zero.real, atol=1e-6,
                err_msg=f"F not continuous from below at delta_k=-{eps}"
            )

    def test_suppression_is_algebraic_not_exponential(self):
        """Suppression of mismatched pairs scales as 1/L not exp(-L).

        |F(delta_k, L)/L| = |sinc(delta_k L/2)| has local maxima near
        delta_k L/2 = (n+0.5)pi, i.e. L = (2n+1)pi/delta_k.  At those
        peaks, |sinc| = 1/((n+0.5)pi) which scales exactly as 1/L because
        L_peak ~ (2n+1) and 1/((n+0.5)pi) ~ 1/n ~ 1/L_peak.

        By sampling at sinc peaks the oscillation is eliminated and the
        log-log slope must be -1 (algebraic) rather than << -1 (exponential).
        """
        delta_k = 1.0
        # Peak L values: L = (2n+1)*pi for n = 10, 100, 1000, 10000
        ns = np.array([10, 100, 1000, 10000])
        L_peaks = (2 * ns + 1) * np.pi / delta_k
        magnitudes = np.array([abs(_F(delta_k, L)) / L for L in L_peaks])

        # Log-log fit: slope should be -1 (algebraic 1/L decay)
        slope = np.polyfit(np.log(L_peaks), np.log(magnitudes), 1)[0]

        npt.assert_allclose(
            slope, -1.0, atol=0.01,
            err_msg=f"Expected 1/L algebraic decay (slope=-1), got slope={slope:.4f}"
        )


# ---------------------------------------------------------------------------
# FS16b: Bulk limits of the coherent pair sum
# ---------------------------------------------------------------------------

class TestFS16BulkLimits:
    """Tests that A(L)/L -> correct tensor in the bulk limit.

    These tests use the two-pair and four-pair models described in the
    module docstring.  The pair assignments for each geometry are made
    explicitly; no geometry switch is used in the formula itself.
    """

    # Canonical pair amplitudes (field overlap × tensor contraction)
    a_LO = 1.0 + 0.3j   # LO-like pair amplitude (R_LO contribution)
    a_TO = 0.5 - 0.7j   # TO-like pair amplitude (R_TO contribution)
    dk = 2.0             # Phase mismatch for mismatched pairs (dimensionless)
    L_bulk = 1e6         # Large L where bulk limit is well established

    def test_backscatter_selects_LO_at_large_L(self):
        """In backscattering, the phase-matched pair carries R_LO.

        Two-pair model:
            LO pair: Delta_k = 0 (phase-matched),   a = a_LO
            TO pair: Delta_k = dk (mismatched),      a = a_TO
        """
        pairs = [(0.0, self.a_LO), (self.dk, self.a_TO)]
        A_large = _coherent_amplitude(pairs, self.L_bulk)
        A_normalized = A_large / self.L_bulk

        # Upper bound on residual from the mismatched pair
        residual_bound = _upper_bound_residual([(self.dk, self.a_TO)], self.L_bulk)

        npt.assert_allclose(
            A_normalized, self.a_LO, atol=residual_bound * 2,
            err_msg=(
                f"Backscatter bulk limit: A/L should approach a_LO={self.a_LO}, "
                f"got {A_normalized:.6g}"
            )
        )

    def test_forward_scatter_selects_TO_at_large_L(self):
        """In forward scattering, the phase-matched pair carries R_TO.

        Two-pair model:
            TO pair: Delta_k = 0 (phase-matched),   a = a_TO
            LO pair: Delta_k = dk (mismatched),      a = a_LO
        """
        pairs = [(0.0, self.a_TO), (self.dk, self.a_LO)]
        A_large = _coherent_amplitude(pairs, self.L_bulk)
        A_normalized = A_large / self.L_bulk

        residual_bound = _upper_bound_residual([(self.dk, self.a_LO)], self.L_bulk)

        npt.assert_allclose(
            A_normalized, self.a_TO, atol=residual_bound * 2,
            err_msg=(
                f"Forward scatter bulk limit: A/L should approach a_TO={self.a_TO}, "
                f"got {A_normalized:.6g}"
            )
        )

    def test_same_formula_for_both_geometries(self):
        """Backscatter and forward scatter use an identical formula.

        Only the pair assignments (which pair has Delta_k = 0) differ.
        This test makes that explicit by constructing both with the same
        function and different pair lists.
        """
        L = self.L_bulk

        pairs_bs = [(0.0, self.a_LO), (self.dk, self.a_TO)]
        pairs_fs = [(0.0, self.a_TO), (self.dk, self.a_LO)]

        A_bs = _coherent_amplitude(pairs_bs, L) / L
        A_fs = _coherent_amplitude(pairs_fs, L) / L

        # They must be different (different geometries, different limits)
        assert abs(A_bs - A_fs) > 0.1, (
            "Backscatter and forward scatter gave the same bulk limit — "
            "they should differ because their phase-matched tensors differ."
        )

        # Each must be close to its respective limit
        assert abs(A_bs - self.a_LO) < 1e-3, f"|A_bs - a_LO| = {abs(A_bs - self.a_LO):.2e}"
        assert abs(A_fs - self.a_TO) < 1e-3, f"|A_fs - a_TO| = {abs(A_fs - self.a_TO):.2e}"

    def test_four_pair_backscattering_model(self):
        """Four-pair model matching the design document specification.

        Design doc step 1: 'two LO-like co-propagating pairs and two
        TO-like cross-pairs'.  For backscattering with q_ext = 2k (k=1):

            fwd-fwd LO pair:  Delta_k = 0    (phase-matched)
            bwd-bwd LO pair:  Delta_k = -4k  (mismatched)
            fwd-bwd TO pair:  Delta_k = -2k  (mismatched)
            bwd-fwd TO pair:  Delta_k = -2k  (mismatched)

        Only the fwd-fwd LO pair survives in the bulk limit.
        The amplitude contributions use distinct complex values to ensure
        the test is not trivially satisfied.
        """
        k = 1.0    # wavevector magnitude (dimensionless)
        a_LO_ff = 1.0 + 0.3j   # fwd-fwd LO pair
        a_LO_bb = 0.6 - 0.2j   # bwd-bwd LO pair (different amplitude)
        a_TO_fb = 0.5 - 0.7j   # fwd-bwd TO cross-pair
        a_TO_bf = 0.4 + 0.1j   # bwd-fwd TO cross-pair

        pairs = [
            (0.0,  a_LO_ff),   # fwd-fwd, phase-matched
            (-4*k, a_LO_bb),   # bwd-bwd, mismatched by 4k
            (-2*k, a_TO_fb),   # fwd-bwd cross-pair, mismatched by 2k
            (-2*k, a_TO_bf),   # bwd-fwd cross-pair, mismatched by 2k
        ]

        L = self.L_bulk
        A_normalized = _coherent_amplitude(pairs, L) / L

        # Bulk limit must be only the fwd-fwd LO pair contribution
        mismatched = [(-4*k, a_LO_bb), (-2*k, a_TO_fb), (-2*k, a_TO_bf)]
        bound = _upper_bound_residual(mismatched, L)

        npt.assert_allclose(
            A_normalized, a_LO_ff, atol=bound * 2,
            err_msg=(
                f"Four-pair backscatter: bulk limit must equal a_LO_ff={a_LO_ff}, "
                f"got {A_normalized:.6g}, bound={bound:.2e}"
            )
        )

    def test_no_zero_q_suppression_switch_needed(self):
        """The cross-pair suppression is automatic from F_ij(L); no explicit switch required.

        This test verifies that:
        1.  A small-L (coherent film) computation retains all four pairs.
        2.  A large-L (bulk) computation suppresses the mismatched pairs.
        3.  Both use the identical formula _coherent_amplitude; only L differs.

        The 'hard zero-q suppression' in the current code is the limit of
        setting F_ij = 0 by hand for cross-pairs.  This test shows that
        the same effect emerges automatically from the phase factor.
        """
        k = 1.0
        a_LO = 1.0 + 0.0j
        a_TO = 0.0 + 1.0j    # purely imaginary so LO and TO are easily distinguished

        pairs = [
            (0.0,  a_LO),    # phase-matched LO
            (-2*k, a_TO),    # mismatched TO cross-pair
        ]

        # Small L: cross-pair contributes significantly
        L_small = 0.5
        A_small = _coherent_amplitude(pairs, L_small) / L_small
        cross_contribution_small = abs(a_TO * _F(-2*k, L_small) / L_small)
        assert cross_contribution_small > 0.1, (
            f"At small L={L_small}, cross-pair should contribute noticeably; "
            f"|cross/L| = {cross_contribution_small:.3f}"
        )

        # Large L: cross-pair is suppressed automatically
        L_large = 1e5
        cross_contribution_large = abs(a_TO * _F(-2*k, L_large) / L_large)
        assert cross_contribution_large < 1e-4, (
            f"At large L={L_large}, cross-pair should be suppressed; "
            f"|cross/L| = {cross_contribution_large:.2e}"
        )

        # Both use the exact same _coherent_amplitude — no switch
        A_large = _coherent_amplitude(pairs, L_large) / L_large
        npt.assert_allclose(
            A_large, a_LO, atol=1e-3,
            err_msg="Large-L limit should equal a_LO without any explicit switch"
        )

    @pytest.mark.parametrize("L_factor", [1, 10, 100, 1000])
    def test_convergence_to_bulk_limit_is_monotone_in_envelope(self, L_factor):
        """The residual |A(L)/L - a_matched| is bounded by 1/(dk L) at each L.

        This checks that the bound is satisfied for several L values,
        confirming the algebraic (1/L) convergence rate.
        """
        dk = 2.0
        a_matched = 1.0 + 0.3j
        a_mismatch = 0.5 - 0.7j

        pairs = [(0.0, a_matched), (dk, a_mismatch)]
        L = 100.0 * L_factor

        A_normalized = _coherent_amplitude(pairs, L) / L
        residual = abs(A_normalized - a_matched)
        bound = _upper_bound_residual([(dk, a_mismatch)], L)

        assert residual <= bound + 1e-14, (
            f"L={L}: residual {residual:.2e} exceeds bound {bound:.2e}"
        )


