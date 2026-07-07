#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""FS3 — External-channel Raman amplitude invariant to Berreman basis rotation.

Design reference: modal_pairs_final_state_design.md, test FS3.

Physical requirement
--------------------
For a finite coherent layer, the Raman amplitude is

    A = E_S_rec(z)^T @ R @ E_L(z)

where E_L and E_S_rec are the total internal fields.  When E_L is expanded
in Berreman eigenmodes as E_L = sum_i c_i u_i, the coherent amplitude
sum

    A_coherent = sum_{i,j} d_j * c_i * (v_j^T @ R @ u_i)

equals E_S_rec^T @ R @ E_L exactly, so it is unchanged by any unitary
rotation of the mode basis that leaves the total field unchanged.

The current incoherent modal-pairs approach,

    I_incoherent = sum_{i,j} |d_j * c_i * (v_j^T @ R @ u_i)|^2,

does NOT satisfy this invariance.  The tests here document both properties
so that the refactored implementation must pass the coherent invariant.

Basis-rotation convention
-------------------------
If forward laser modes (columns) are transformed as

    [u'_0, u'_1] = [u_0, u_1] @ U^T      (U is 2x2 unitary)

then to keep E_L = sum_i c_i u_i unchanged, amplitudes must transform as

    [c'_0, c'_1] = conj(U) @ [c_0, c_1]

Proof: E_L fixed requires U^T @ conj(U) = I, which holds for any unitary U
because U^T @ conj(U) = U^T @ (U^H)^T = (U^H @ U)^T = I^T = I.

All tests work at a single z = 0 (no depth integration needed to establish
the invariance; the integral factorises for fixed z).
"""
import numpy as np
import numpy.testing as npt
import pytest
from numpy.linalg import qr


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _random_unitary_2x2(rng):
    """Return a random 2x2 complex unitary matrix via QR decomposition."""
    A = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
    Q, _ = qr(A)
    return Q


def _coherent_amplitude(R, c_L, modes_L, d_S, modes_S):
    """Coherent external-channel Raman amplitude at z = 0.

    Parameters
    ----------
    R : (3, 3) complex ndarray
        Raman tensor.
    c_L : (n_L,) complex ndarray
        Scalar amplitudes for each laser Berreman mode.
    modes_L : (3, n_L) complex ndarray
        Laser Berreman E-field 3-vectors as columns.
    d_S : (n_S,) complex ndarray
        Scalar amplitudes for each reciprocal-detector Berreman mode.
    modes_S : (3, n_S) complex ndarray
        Reciprocal-detector Berreman E-field 3-vectors as columns.

    Returns
    -------
    complex
        A = sum_{i,j} d_j * c_i * (v_j^T @ R @ u_i)
    """
    # Reconstruct total fields and use the direct bilinear product.
    # This is exactly the modal sum by linearity.
    E_L = modes_L @ c_L          # (3,)
    E_S = modes_S @ d_S          # (3,)
    return E_S @ R @ E_L


def _modal_coherent_amplitude(R, c_L, modes_L, d_S, modes_S):
    """Same as _coherent_amplitude but computed by explicit double loop.

    Used to verify that the loop form equals the direct bilinear product.
    """
    A = 0.0 + 0.0j
    for i in range(len(c_L)):
        for j in range(len(d_S)):
            A += d_S[j] * c_L[i] * (modes_S[:, j] @ R @ modes_L[:, i])
    return A


def _incoherent_intensity(R, c_L, modes_L, d_S, modes_S):
    """Incoherent modal-pairs intensity (current approach).

    I = sum_{i,j} |d_j * c_i * (v_j^T @ R @ u_i)|^2
    """
    I = 0.0
    for i in range(len(c_L)):
        for j in range(len(d_S)):
            pair = d_S[j] * c_L[i] * (modes_S[:, j] @ R @ modes_L[:, i])
            I += abs(pair) ** 2
    return I


def _rotate_forward_modes(modes_L, c_L, U):
    """Apply 2x2 unitary U to forward laser modes (columns 0 and 1).

    Returns
    -------
    modes_L_rot : (3, n_L) ndarray
        Modes with forward pair rotated: [u'_0, u'_1] = [u_0, u_1] @ U^T
    c_L_rot : (n_L,) ndarray
        Amplitudes with forward pair transformed: c'_fwd = conj(U) @ c_fwd
        (keeps total field E_L unchanged)
    """
    modes_L_rot = modes_L.copy()
    c_L_rot = c_L.copy()
    modes_L_rot[:, :2] = modes_L[:, :2] @ U.T
    c_L_rot[:2] = np.conj(U) @ c_L[:2]
    return modes_L_rot, c_L_rot


# ---------------------------------------------------------------------------
# FS3: Basis invariance — pure numpy, no GTM, no DFT data
# ---------------------------------------------------------------------------

class TestFS3BasisInvariance:
    """FS3 — Coherent amplitude is invariant to Berreman basis rotation.

    These tests use synthetic random inputs and pass without any GTM or DFT
    data.  They define the mathematical invariant that the new implementation
    must satisfy.
    """

    @pytest.fixture
    def rng(self):
        return np.random.default_rng(42)

    @pytest.fixture
    def raman_tensor(self):
        """Non-trivial off-diagonal Raman tensor."""
        R = np.array([
            [1.0,  0.5,  0.0],
            [0.5,  0.3, -0.2],
            [0.0, -0.2,  0.1],
        ], dtype=complex)
        return R

    @pytest.fixture
    def field_setup(self, rng):
        """Synthetic Berreman field: 2 forward + 2 backward laser modes,
        2 forward + 2 backward reciprocal-detector modes.

        modes are (3, 4) arrays; amplitudes are (4,) arrays.
        Columns 0,1 are 'forward'; columns 2,3 are 'backward'.
        """
        modes_L = (rng.standard_normal((3, 4)) +
                   1j * rng.standard_normal((3, 4)))
        c_L = (rng.standard_normal(4) +
               1j * rng.standard_normal(4))

        modes_S = (rng.standard_normal((3, 4)) +
                   1j * rng.standard_normal((3, 4)))
        d_S = (rng.standard_normal(4) +
               1j * rng.standard_normal(4))

        return modes_L, c_L, modes_S, d_S

    # --- Sanity checks for helper functions ---

    def test_modal_loop_equals_direct_product(self, raman_tensor, field_setup):
        """Loop-form modal sum and direct bilinear product must agree."""
        modes_L, c_L, modes_S, d_S = field_setup
        R = raman_tensor

        A_direct = _coherent_amplitude(R, c_L, modes_L, d_S, modes_S)
        A_loop = _modal_coherent_amplitude(R, c_L, modes_L, d_S, modes_S)

        assert abs(A_loop - A_direct) < 1e-12, (
            f"Loop and direct amplitudes disagree: |delta| = {abs(A_loop - A_direct):.2e}"
        )

    def test_rotation_preserves_total_field(self, rng, field_setup):
        """After rotation, E_L_rotated must equal E_L_original."""
        modes_L, c_L, modes_S, d_S = field_setup
        U = _random_unitary_2x2(rng)

        E_L_original = modes_L @ c_L
        modes_L_rot, c_L_rot = _rotate_forward_modes(modes_L, c_L, U)
        E_L_rotated = modes_L_rot @ c_L_rot

        npt.assert_allclose(
            E_L_rotated, E_L_original, atol=1e-12,
            err_msg="Rotation convention does not preserve E_L total field"
        )

    # --- Core invariance tests ---

    def test_coherent_amplitude_invariant_under_real_rotation(
            self, raman_tensor, field_setup):
        """Coherent amplitude is unchanged by a real rotation of forward modes."""
        modes_L, c_L, modes_S, d_S = field_setup
        R = raman_tensor

        theta = np.radians(37.0)
        U_real = np.array([[ np.cos(theta), np.sin(theta)],
                           [-np.sin(theta), np.cos(theta)]])

        A_orig = _coherent_amplitude(R, c_L, modes_L, d_S, modes_S)
        modes_rot, c_rot = _rotate_forward_modes(modes_L, c_L, U_real)
        A_rot = _coherent_amplitude(R, c_rot, modes_rot, d_S, modes_S)

        assert abs(A_rot - A_orig) < 1e-12, (
            f"Coherent amplitude changed under real rotation: |delta| = {abs(A_rot - A_orig):.2e}"
        )

    def test_coherent_amplitude_invariant_under_complex_unitary_rotation(
            self, rng, raman_tensor, field_setup):
        """Coherent amplitude is unchanged by a complex unitary rotation of forward modes."""
        modes_L, c_L, modes_S, d_S = field_setup
        R = raman_tensor
        U = _random_unitary_2x2(rng)

        A_orig = _coherent_amplitude(R, c_L, modes_L, d_S, modes_S)
        modes_rot, c_rot = _rotate_forward_modes(modes_L, c_L, U)
        A_rot = _coherent_amplitude(R, c_rot, modes_rot, d_S, modes_S)

        assert abs(A_rot - A_orig) < 1e-12, (
            f"Coherent amplitude changed under complex unitary rotation: "
            f"|delta| = {abs(A_rot - A_orig):.2e}"
        )

    @pytest.mark.parametrize("theta_deg", [15.0, 37.0, 63.0, 90.0, 135.0])
    def test_coherent_amplitude_invariant_for_multiple_rotations(
            self, theta_deg, raman_tensor, field_setup):
        """Coherent amplitude is unchanged for several rotation angles."""
        modes_L, c_L, modes_S, d_S = field_setup
        R = raman_tensor

        theta = np.radians(theta_deg)
        U = np.array([[ np.cos(theta), np.sin(theta)],
                      [-np.sin(theta), np.cos(theta)]])

        A_orig = _coherent_amplitude(R, c_L, modes_L, d_S, modes_S)
        modes_rot, c_rot = _rotate_forward_modes(modes_L, c_L, U)
        A_rot = _coherent_amplitude(R, c_rot, modes_rot, d_S, modes_S)

        assert abs(A_rot - A_orig) < 1e-12, (
            f"theta={theta_deg} deg: |delta| = {abs(A_rot - A_orig):.2e}"
        )

    # --- Demonstration that incoherent sum FAILS the invariance ---

    def test_incoherent_sum_changes_under_forward_rotation(
            self, rng, field_setup):
        """Incoherent modal-pair sum changes under forward-mode rotation.

        This is the failure mode that motivates the coherent-amplitude
        approach.  The test uses a setup designed to give a large relative
        change: c_L[1] = 0 so all laser intensity is in mode 0 at theta=0,
        but is shared after a 45-degree rotation.
        """
        modes_L, c_L, modes_S, d_S = field_setup

        # Force c_L[1] = 0 so mode 1 carries no amplitude at theta=0
        c_L_skewed = c_L.copy()
        c_L_skewed[1] = 0.0

        R = np.diag([1.0 + 0j, 0.0, 0.0])   # only Ex contributes

        I_orig = _incoherent_intensity(R, c_L_skewed, modes_L, d_S, modes_S)

        theta = np.radians(45.0)
        U = np.array([[ np.cos(theta), np.sin(theta)],
                      [-np.sin(theta), np.cos(theta)]])
        modes_rot, c_rot = _rotate_forward_modes(modes_L, c_L_skewed, U)
        I_rot = _incoherent_intensity(R, c_rot, modes_rot, d_S, modes_S)

        assert abs(I_rot - I_orig) > 1e-8, (
            "Incoherent sum was unexpectedly invariant — the failure mode "
            "used to motivate the refactor was not observed.  "
            f"I_orig={I_orig:.6g}, I_rot={I_rot:.6g}"
        )


# ---------------------------------------------------------------------------
# Stub: FS3 using the future PhononFinalStateResolver + RamanAmplitudeAccumulator
# ---------------------------------------------------------------------------

class TestFS3WithImplementation:
    """Stub tests for FS3 that will exercise the refactored implementation.

    These are marked skip until Phase 1 (PhononFinalStateResolver) and
    Phase 2 (RamanAmplitudeAccumulator) are complete.  Replace the skip
    marker and fill in the body at that point.

    The tests in TestFS3BasisInvariance above establish the mathematical
    property; these tests will verify that the implementation enforces it.
    """

    @pytest.mark.skip(reason="PhononFinalStateResolver not yet implemented — Phase 1")
    def test_accumulator_sums_coherently_across_berreman_pairs(self):
        """RamanAmplitudeAccumulator must sum amplitudes, not intensities."""
        # Imports will succeed once Phase 1/2 are complete.
        from PDielec.PhononFinalStateResolver import PhononFinalStateResolver   # noqa: F401
        from PDielec.RamanAmplitudeAccumulator import RamanAmplitudeAccumulator  # noqa: F401
        raise NotImplementedError(
            "Fill in using the API defined in modal_pairs_final_state_design.md "
            "once PhononFinalStateResolver and RamanAmplitudeAccumulator exist."
        )

    @pytest.mark.skip(reason="PhononFinalStateResolver not yet implemented — Phase 1")
    def test_coherent_film_regime_gives_basis_independent_result(self):
        """In coherent_film regime, two basis-rotated setups give the same intensity."""
        from PDielec.PhononFinalStateResolver import PhononFinalStateResolver   # noqa: F401
        raise NotImplementedError("Fill in after Phase 1")
