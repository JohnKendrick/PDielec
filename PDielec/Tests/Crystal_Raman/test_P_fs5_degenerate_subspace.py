#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""FS5 — Raman intensity invariant to rotation within a degenerate Berreman subspace.

Design reference: modal_pairs_final_state_design.md, test FS5.

Physical requirement
--------------------
At normal incidence on a uniaxial or isotropic medium, the two forward
Berreman modes are exactly degenerate (same kz eigenvalue).  Any orthogonal
rotation within this 2D subspace leaves the dispersion unchanged but
permutes the eigenvectors.  The Raman amplitude, and therefore the
intensity, must be independent of this basis choice.

The coherent amplitude

    A = d * [c_a * (v^T @ R @ u_a) + c_b * (v^T @ R @ u_b)]

equals d * v^T @ R @ E_L_degen where E_L_degen = c_a u_a + c_b u_b is the
total field in the degenerate subspace.  This is invariant to rotations
that keep E_L_degen fixed.

The incoherent pair-sum

    I = |d * c_a * (v^T @ R @ u_a)|^2 + |d * c_b * (v^T @ R @ u_b)|^2

is generally NOT invariant.

Rotation convention
-------------------
For a 2D real rotation by angle theta:

    u'_a =  cos(theta) u_a + sin(theta) u_b
    u'_b = -sin(theta) u_a + cos(theta) u_b

To keep E_L_degen = c_a u_a + c_b u_b = c'_a u'_a + c'_b u'_b fixed:

    c'_a =  cos(theta) c_a + sin(theta) c_b
    c'_b = -sin(theta) c_a + cos(theta) c_b

(For real orthogonal U, conj(U) = U, so the amplitude transform c' = U c
uses the same rotation matrix as the mode transform [u'_a, u'_b] = U [u_a, u_b].)

The projector approach
-----------------------
The design document specifies that the refactored implementation must
detect degenerate optical subspaces and compute their contribution via
the subspace projector

    P = u_a u_a^H + u_b u_b^H

rather than via individual modes.  The amplitude contribution from P is
then sum-of-basis-amplitudes, which equals v^T @ R @ P E_L and is
basis-invariant by construction.

These tests verify the mathematical property directly.  The implementation
stubs at the end of this file will exercise the PhononFinalStateResolver
degenerate-subspace code once it is written.
"""
import numpy as np
import numpy.testing as npt
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rotate_degenerate_pair(u_a, u_b, c_a, c_b, theta_deg):
    """Rotate the degenerate forward-mode pair by theta_deg.

    Returns (u_a', u_b', c_a', c_b') with E_L_degen preserved.
    """
    theta = np.radians(theta_deg)
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    u_a_rot = cos_t * u_a + sin_t * u_b
    u_b_rot = -sin_t * u_a + cos_t * u_b
    c_a_rot = cos_t * c_a + sin_t * c_b
    c_b_rot = -sin_t * c_a + cos_t * c_b

    return u_a_rot, u_b_rot, c_a_rot, c_b_rot


def _coherent_amplitude_degen(R, u_a, u_b, c_a, c_b, v, d):
    """Coherent amplitude from a degenerate forward-mode pair.

        A = d * [c_a * (v^T @ R @ u_a) + c_b * (v^T @ R @ u_b)]
    """
    return d * (c_a * (v @ R @ u_a) + c_b * (v @ R @ u_b))


def _incoherent_intensity_degen(R, u_a, u_b, c_a, c_b, v, d):
    """Incoherent pair intensity from a degenerate forward-mode pair.

        I = |d * c_a * (v^T @ R @ u_a)|^2 + |d * c_b * (v^T @ R @ u_b)|^2
    """
    amp_a = d * c_a * (v @ R @ u_a)
    amp_b = d * c_b * (v @ R @ u_b)
    return abs(amp_a) ** 2 + abs(amp_b) ** 2


def _projector_amplitude(R, u_a, u_b, E_L_degen, v, d):
    """Amplitude computed via the subspace projector P = u_a u_a^H + u_b u_b^H.

    A = d * v^T @ R @ P @ E_L_degen

    This is the implementation path the refactor must use.
    """
    P_E_L = u_a * np.dot(np.conj(u_a), E_L_degen) + u_b * np.dot(np.conj(u_b), E_L_degen)
    return d * (v @ R @ P_E_L)


# ---------------------------------------------------------------------------
# FS5: Degenerate subspace invariance — pure numpy, no GTM, no DFT data
# ---------------------------------------------------------------------------

class TestFS5DegenerateSubspaceRotation:
    """FS5 — Coherent amplitude invariant to in-subspace basis rotation.

    Uses two orthonormal degenerate modes in the xy plane to mimic the
    normal-incidence uniaxial geometry.  All tests pass without GTM or DFT.
    """

    @pytest.fixture
    def rng(self):
        return np.random.default_rng(123)

    @pytest.fixture
    def degenerate_modes(self, rng):
        """Two orthonormal 3-vectors spanning the xy plane (degenerate at normal incidence).

        Amplitudes and detector field are random complex.
        """
        u_a = np.array([1.0, 0.0, 0.0], dtype=complex)
        u_b = np.array([0.0, 1.0, 0.0], dtype=complex)

        c_a = rng.standard_normal() + 1j * rng.standard_normal()
        c_b = rng.standard_normal() + 1j * rng.standard_normal()

        v = np.array([1.0, 0.5, 0.0], dtype=complex)   # detector mode
        d = rng.standard_normal() + 1j * rng.standard_normal()

        return u_a, u_b, c_a, c_b, v, d

    @pytest.fixture
    def off_diagonal_R(self):
        """Off-diagonal Raman tensor: reveals basis-dependence in incoherent sum."""
        return np.array([
            [0.0,  1.0,  0.0],
            [1.0,  0.0,  0.0],
            [0.0,  0.0,  0.0],
        ], dtype=complex)

    @pytest.fixture
    def general_R(self):
        """Generic Raman tensor for completeness tests."""
        return np.array([
            [1.0,  0.5,  0.2],
            [0.5,  0.3, -0.1],
            [0.2, -0.1,  0.0],
        ], dtype=complex)

    # --- Sanity check ---

    def test_rotation_preserves_degenerate_field(self, degenerate_modes):
        """Rotated modes and amplitudes must give the same E_L_degen."""
        u_a, u_b, c_a, c_b, v, d = degenerate_modes

        E_orig = c_a * u_a + c_b * u_b

        for theta_deg in [0.0, 30.0, 45.0, 90.0]:
            u_a_r, u_b_r, c_a_r, c_b_r = _rotate_degenerate_pair(
                u_a, u_b, c_a, c_b, theta_deg
            )
            E_rot = c_a_r * u_a_r + c_b_r * u_b_r
            npt.assert_allclose(
                E_rot, E_orig, atol=1e-14,
                err_msg=f"theta={theta_deg}: rotation does not preserve E_L_degen"
            )

    # --- Core coherent invariance ---

    @pytest.mark.parametrize("theta_deg", [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0])
    def test_coherent_amplitude_invariant_off_diagonal_R(
            self, degenerate_modes, off_diagonal_R, theta_deg):
        """Coherent amplitude unchanged for off-diagonal R (most sensitive case)."""
        u_a, u_b, c_a, c_b, v, d = degenerate_modes
        R = off_diagonal_R

        A_ref = _coherent_amplitude_degen(R, u_a, u_b, c_a, c_b, v, d)

        u_a_r, u_b_r, c_a_r, c_b_r = _rotate_degenerate_pair(
            u_a, u_b, c_a, c_b, theta_deg
        )
        A_rot = _coherent_amplitude_degen(R, u_a_r, u_b_r, c_a_r, c_b_r, v, d)

        assert abs(A_rot - A_ref) < 1e-12, (
            f"theta={theta_deg} deg, off-diagonal R: "
            f"A_ref={A_ref:.6g}, A_rot={A_rot:.6g}, "
            f"|delta|={abs(A_rot - A_ref):.2e}"
        )

    @pytest.mark.parametrize("theta_deg", [0.0, 22.5, 45.0, 67.5, 90.0])
    def test_coherent_amplitude_invariant_general_R(
            self, degenerate_modes, general_R, theta_deg):
        """Coherent amplitude unchanged for a general Raman tensor."""
        u_a, u_b, c_a, c_b, v, d = degenerate_modes
        R = general_R

        A_ref = _coherent_amplitude_degen(R, u_a, u_b, c_a, c_b, v, d)

        u_a_r, u_b_r, c_a_r, c_b_r = _rotate_degenerate_pair(
            u_a, u_b, c_a, c_b, theta_deg
        )
        A_rot = _coherent_amplitude_degen(R, u_a_r, u_b_r, c_a_r, c_b_r, v, d)

        assert abs(A_rot - A_ref) < 1e-12, (
            f"theta={theta_deg} deg, general R: |delta|={abs(A_rot - A_ref):.2e}"
        )

    # --- Projector approach gives same result as coherent sum ---

    @pytest.mark.parametrize("theta_deg", [0.0, 30.0, 60.0, 90.0])
    def test_projector_equals_coherent_amplitude(
            self, degenerate_modes, general_R, theta_deg):
        """Subspace projector P gives the same amplitude as explicit coherent sum.

        This tests the key property of the degenerate-subspace implementation
        path: P @ E_L = sum_i u_i (u_i^H E_L) for a complete orthonormal set,
        so v^T @ R @ P @ E_L = v^T @ R @ E_L_degen.
        """
        u_a, u_b, c_a, c_b, v, d = degenerate_modes
        R = general_R

        u_a_r, u_b_r, c_a_r, c_b_r = _rotate_degenerate_pair(
            u_a, u_b, c_a, c_b, theta_deg
        )
        E_L_degen = c_a_r * u_a_r + c_b_r * u_b_r

        A_coherent = _coherent_amplitude_degen(
            R, u_a_r, u_b_r, c_a_r, c_b_r, v, d
        )
        A_projector = _projector_amplitude(R, u_a_r, u_b_r, E_L_degen, v, d)

        assert abs(A_projector - A_coherent) < 1e-12, (
            f"theta={theta_deg} deg: projector and coherent sum disagree, "
            f"|delta|={abs(A_projector - A_coherent):.2e}"
        )

    # --- Demonstration that incoherent sum FAILS the invariance ---

    def test_incoherent_sum_changes_at_45_degrees(
            self, degenerate_modes, off_diagonal_R):
        """Incoherent sum changes under 45° rotation for off-diagonal R.

        Setup: c_b = 0, u_a = [1,0,0], R = [[0,1,0],[1,0,0],[0,0,0]].
        At theta=0, only mode a contributes and the overlap v^T R u_a
        picks out the off-diagonal element.  At theta=45, both modes carry
        equal amplitude and their overlaps differ, giving a different sum of
        squared pair amplitudes.

        This is the failure mode that motivates the projector approach.
        """
        u_a, u_b, c_a, _, v, d = degenerate_modes
        R = off_diagonal_R
        c_b_zero = 0.0 + 0.0j

        I_orig = _incoherent_intensity_degen(
            R, u_a, u_b, c_a, c_b_zero, v, d
        )

        u_a_r, u_b_r, c_a_r, c_b_r = _rotate_degenerate_pair(
            u_a, u_b, c_a, c_b_zero, 45.0
        )
        I_rot = _incoherent_intensity_degen(
            R, u_a_r, u_b_r, c_a_r, c_b_r, v, d
        )

        # The incoherent sum must differ (by ~factor 2 for this setup)
        assert abs(I_rot - I_orig) > 1e-8, (
            "Incoherent sum was unexpectedly invariant at 45°.  "
            f"I_orig={I_orig:.6g}, I_rot={I_rot:.6g}"
        )

    def test_incoherent_failure_is_large_not_numerical_noise(
            self, off_diagonal_R):
        """The incoherent failure for the canonical example is a factor of 2.

        Setup: u_a = [1,0,0], u_b = [0,1,0], R = [[0,1,0],[1,0,0],[0,0,0]],
        v = [0,1,0], c_a = 1, c_b = 0, d = 1.

        At theta = 0:
            R @ u_a = [0,1,0]; v^T R u_a = [0,1,0].[0,1,0] = 1
            R @ u_b = [1,0,0]; v^T R u_b = [0,1,0].[1,0,0] = 0
            I_0 = |1*1*1|^2 + |1*0*0|^2 = 1

        At theta = 45°:
            u_a' = [1/√2, 1/√2, 0],  c_a' =  1/√2
            u_b' = [-1/√2, 1/√2, 0], c_b' = -1/√2
            R @ u_a' = [1/√2, 1/√2, 0]; v^T R u_a' = 1/√2
            R @ u_b' = [1/√2, -1/√2, 0]; v^T R u_b' = -1/√2
            amp_a = (1/√2)(1/√2) = 1/2
            amp_b = (-1/√2)(-1/√2) = 1/2
            I_45 = |1/2|^2 + |1/2|^2 = 0.5  (factor 2 error vs I_0 = 1)

        The coherent amplitude is constant = 1 for all theta (invariant).
        """
        R = off_diagonal_R
        u_a = np.array([1.0, 0.0, 0.0], dtype=complex)
        u_b = np.array([0.0, 1.0, 0.0], dtype=complex)
        v = np.array([0.0, 1.0, 0.0], dtype=complex)
        c_a = 1.0 + 0j
        c_b = 0.0 + 0j
        d = 1.0 + 0j

        I_0 = _incoherent_intensity_degen(R, u_a, u_b, c_a, c_b, v, d)
        npt.assert_allclose(I_0, 1.0, rtol=1e-14, err_msg="I at theta=0 must be 1")

        u_a_r, u_b_r, c_a_r, c_b_r = _rotate_degenerate_pair(
            u_a, u_b, c_a, c_b, 45.0
        )
        I_45 = _incoherent_intensity_degen(R, u_a_r, u_b_r, c_a_r, c_b_r, v, d)
        npt.assert_allclose(I_45, 0.5, rtol=1e-14,
                            err_msg="I at theta=45 must be 0.5 (factor 2 error)")

        assert abs(I_45 - I_0) / I_0 >= 0.4, (
            "Incoherent failure is smaller than expected — check the tensor setup."
        )


# ---------------------------------------------------------------------------
# Stub: FS5 using the future PhononFinalStateResolver degenerate-subspace code
# ---------------------------------------------------------------------------

class TestFS5WithImplementation:
    """Stubs for FS5 that will exercise the degenerate-subspace projector
    inside PhononFinalStateResolver once Phase 1 is complete.

    Replace the skip marker and fill in the body at that point.
    """

    @pytest.mark.skip(reason="PhononFinalStateResolver not yet implemented — Phase 1")
    def test_resolver_detects_degenerate_optical_subspace(self):
        """PhononFinalStateResolver must detect kz degeneracy and use the projector."""
        from PDielec.PhononFinalStateResolver import PhononFinalStateResolver  # noqa: F401
        raise NotImplementedError(
            "Construct two Berreman modes with |kz_a - kz_b|/mean < eps_optical "
            "and verify the resolver returns a subspace-projected contribution "
            "rather than two separate pair amplitudes."
        )

    @pytest.mark.skip(reason="PhononFinalStateResolver not yet implemented — Phase 1")
    def test_resolver_gives_basis_independent_intensity_at_normal_incidence(self):
        """At normal incidence, intensity must not change when Berreman basis is rotated."""
        from PDielec.PhononFinalStateResolver import PhononFinalStateResolver  # noqa: F401
        raise NotImplementedError(
            "Run PhononFinalStateResolver at normal incidence with a uniaxial "
            "permittivity tensor and two different orthonormal bases for the "
            "degenerate forward-mode subspace.  Assert intensities agree to 1e-8."
        )
