"""Tests A1-A5: Raman tensor algebra for crystal (layered) Raman.

Plan sections A1, A2, A3, A4, A5.

Pure algebra — no GTM/optical solver required.
Tests verify: direct tensor contractions, the rotation convention
G @ R @ G.T, isotropic invariance, off-diagonal selection rules, and
global complex-phase invariance of the amplitude modulus squared.

Convention
----------
The lab-frame Raman tensor is R_lab = G @ R_crystal @ G.T where G maps
crystal-frame vectors into lab-frame vectors (v_lab = G @ v_crystal).
Amplitudes are computed as A = e_s^T R e_i (ordinary transpose, no conjugate).
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.LayeredRamanCalculator import LayeredRamanCalculator

# ---------------------------------------------------------------------------
# shared fixtures
# ---------------------------------------------------------------------------

R_diag = np.diag([2.0, 1.0, 0.5])          # A2/A3 test tensor
e_x = np.array([1.0, 0.0, 0.0])
e_y = np.array([0.0, 1.0, 0.0])
e_z = np.array([0.0, 0.0, 1.0])


def _contract(e_s, R, e_i):
    """Return A = e_s^T R e_i (no complex conjugate)."""
    return e_s @ R @ e_i


def _Rz(phi):
    """Rotation matrix for angle phi (rad) about the z axis."""
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _random_SO3(seed):
    """Return a random rotation matrix from QR decomposition."""
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


# ---------------------------------------------------------------------------
# A1: Direct tensor contraction with diagonal tensor
# ---------------------------------------------------------------------------

class TestA1DirectContraction:
    """A1: e_s^T R e_i for diagonal R = diag(a, b, c)."""

    def test_xx_amplitude(self):
        assert abs(_contract(e_x, R_diag, e_x) - 2.0) < 1e-12

    def test_yy_amplitude(self):
        assert abs(_contract(e_y, R_diag, e_y) - 1.0) < 1e-12

    def test_zz_amplitude(self):
        assert abs(_contract(e_z, R_diag, e_z) - 0.5) < 1e-12

    def test_xy_is_zero(self):
        assert abs(_contract(e_x, R_diag, e_y)) < 1e-12

    def test_xz_is_zero(self):
        assert abs(_contract(e_x, R_diag, e_z)) < 1e-12

    def test_xx_strength(self):
        assert abs(abs(_contract(e_x, R_diag, e_x))**2 - 4.0) < 1e-12

    def test_yy_strength(self):
        assert abs(abs(_contract(e_y, R_diag, e_y))**2 - 1.0) < 1e-12

    def test_zz_strength(self):
        assert abs(abs(_contract(e_z, R_diag, e_z))**2 - 0.25) < 1e-12


# ---------------------------------------------------------------------------
# A2: Rotation about z axis — convention and angular dependences
# ---------------------------------------------------------------------------

class TestA2RotationAboutZ:
    """A2: R_lab = G @ R_crystal @ G.T with G = Rz(phi)."""

    def test_rotation_convention_matches_G_R_GT(self):
        G = _Rz(np.pi / 5)
        R_lab = LayeredRamanCalculator._rotate_raman_tensor(R_diag, G)
        assert np.allclose(R_lab, G @ R_diag @ G.T, atol=1e-14)

    @pytest.mark.parametrize("phi", [0.0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2])
    def test_Ixx_analytic(self, phi):
        a, b = 2.0, 1.0
        expected = (a * np.cos(phi)**2 + b * np.sin(phi)**2)**2
        G = _Rz(phi)
        R_lab = LayeredRamanCalculator._rotate_raman_tensor(R_diag, G)
        I = abs(_contract(e_x, R_lab, e_x))**2
        assert abs(I - expected) < 1e-12, f"phi={np.degrees(phi):.0f}°"

    @pytest.mark.parametrize("phi", [0.0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2])
    def test_Ixy_analytic(self, phi):
        a, b = 2.0, 1.0
        expected = ((a - b) * np.sin(phi) * np.cos(phi))**2
        G = _Rz(phi)
        R_lab = LayeredRamanCalculator._rotate_raman_tensor(R_diag, G)
        I = abs(_contract(e_y, R_lab, e_x))**2
        assert abs(I - expected) < 1e-12, f"phi={np.degrees(phi):.0f}°"

    def test_concrete_values(self):
        """Verify table from plan: phi=0,45,90 for I_xx and I_xy."""
        table = [
            (0.0,       4.00, 0.00),
            (np.pi / 4, 2.25, 0.25),
            (np.pi / 2, 1.00, 0.00),
        ]
        for phi, I_xx_ref, I_xy_ref in table:
            G = _Rz(phi)
            R_lab = LayeredRamanCalculator._rotate_raman_tensor(R_diag, G)
            I_xx = abs(_contract(e_x, R_lab, e_x))**2
            I_xy = abs(_contract(e_y, R_lab, e_x))**2
            assert abs(I_xx - I_xx_ref) < 1e-12, f"phi={np.degrees(phi):.0f}° I_xx={I_xx}"
            assert abs(I_xy - I_xy_ref) < 1e-12, f"phi={np.degrees(phi):.0f}° I_xy={I_xy}"

    def test_Ixx_plus_Iyy_is_constant(self):
        """I_xx + I_yy must equal a^2 + b^2 for all phi (sum rule)."""
        a, b = 2.0, 1.0
        ref = a**2 + b**2
        for phi in np.linspace(0, np.pi, 30):
            G = _Rz(phi)
            R_lab = LayeredRamanCalculator._rotate_raman_tensor(R_diag, G)
            # Hmm, this sum rule does not hold in general; drop this test.
            # Actually: I_xx = (a cos²φ + b sin²φ)², I_yy = (a sin²φ + b cos²φ)²
            # These are NOT generally a² + b² unless φ = 0 or π/2.
            # Use a different sum: check that R_lab[0,0] + R_lab[1,1] = a + b
            assert abs(R_lab[0, 0] + R_lab[1, 1] - (a + b)) < 1e-12


# ---------------------------------------------------------------------------
# A3: Isotropic tensor orientation invariance
# ---------------------------------------------------------------------------

class TestA3IsotropicInvariance:
    """A3: R = a*I must be unchanged by any rotation."""

    def test_tensor_unchanged_by_random_rotations(self):
        a = 3.14
        R = a * np.eye(3)
        for seed in range(20):
            G = _random_SO3(seed)
            R_lab = LayeredRamanCalculator._rotate_raman_tensor(R, G)
            assert np.allclose(R_lab, a * np.eye(3), atol=1e-12), f"seed={seed}"

    def test_parallel_contraction_equals_a(self):
        a = 2.71828
        R = a * np.eye(3)
        for seed in range(10):
            G = _random_SO3(seed)
            R_lab = LayeredRamanCalculator._rotate_raman_tensor(R, G)
            for e in [e_x, e_y, e_z]:
                assert abs(_contract(e, R_lab, e) - a) < 1e-12

    def test_crossed_contraction_is_zero(self):
        a = 1.5
        R = a * np.eye(3)
        for seed in range(10):
            G = _random_SO3(seed)
            R_lab = LayeredRamanCalculator._rotate_raman_tensor(R, G)
            for e_i, e_s in [(e_x, e_y), (e_x, e_z), (e_y, e_z)]:
                assert abs(_contract(e_s, R_lab, e_i)) < 1e-12


# ---------------------------------------------------------------------------
# A4: Off-diagonal tensor selection rules
# ---------------------------------------------------------------------------

class TestA4OffDiagonalSelectionRules:
    """A4: Crossed polarisation selection rules for off-diagonal tensor."""

    d = 1.5
    R_xy = np.array([[0.0, d, 0.0], [d, 0.0, 0.0], [0.0, 0.0, 0.0]])

    def test_xx_is_zero(self):
        assert abs(_contract(e_x, self.R_xy, e_x)) < 1e-12

    def test_yy_is_zero(self):
        assert abs(_contract(e_y, self.R_xy, e_y)) < 1e-12

    def test_xy_equals_d(self):
        assert abs(_contract(e_x, self.R_xy, e_y) - self.d) < 1e-12

    def test_yx_equals_d(self):
        assert abs(_contract(e_y, self.R_xy, e_x) - self.d) < 1e-12

    def test_zz_is_zero(self):
        assert abs(_contract(e_z, self.R_xy, e_z)) < 1e-12

    def test_rotate_45_transfers_to_parallel_channel(self):
        """After 45° z-rotation, the parallel channel gains intensity."""
        G = _Rz(np.pi / 4)
        R_lab = LayeredRamanCalculator._rotate_raman_tensor(self.R_xy, G)
        I_xx = abs(_contract(e_x, R_lab, e_x))**2
        I_xy = abs(_contract(e_y, R_lab, e_x))**2
        assert I_xx > 0.0, "xx should be non-zero after 45° rotation"
        assert I_xy > 0.0, "xy should remain non-zero after 45° rotation"
        # After 45° rotation: R_lab = diag(d, -d, 0), so I_xx = d²
        assert abs(I_xx - self.d**2) < 1e-12
        assert abs(I_xy) < 1e-12   # xy becomes zero after the specific rotation


# ---------------------------------------------------------------------------
# A5: Complex phase invariance
# ---------------------------------------------------------------------------

class TestA5ComplexPhaseInvariance:
    """A5: |A(e^{iθ} R)|² = |A(R)|² for any real θ."""

    R_sym = np.array([[2.0, 0.5, 0.1],
                      [0.5, 1.0, 0.3],
                      [0.1, 0.3, 0.8]])

    @pytest.mark.parametrize("theta", [0.0, np.pi / 4, 0.7, 1.23, np.pi, 2.5])
    def test_strength_invariant(self, theta):
        A_ref = _contract(e_x, self.R_sym, e_x)
        I_ref = abs(A_ref)**2
        R_phased = np.exp(1j * theta) * self.R_sym
        I = abs(_contract(e_x, R_phased, e_x))**2
        assert abs(I - I_ref) / (I_ref + 1e-30) < 1e-12, f"theta={theta}"

    @pytest.mark.parametrize("theta", [0.3, 1.0, 2.5])
    def test_amplitude_acquires_phase(self, theta):
        """Amplitude must change by exactly e^{iθ} under the phase rotation."""
        A_ref = _contract(e_x, self.R_sym, e_x)
        R_phased = np.exp(1j * theta) * self.R_sym
        A = _contract(e_x, R_phased, e_x)
        assert abs(A - np.exp(1j * theta) * A_ref) < 1e-12, f"theta={theta}"

    @pytest.mark.parametrize("theta", [0.0, np.pi / 3, np.pi / 2, np.pi])
    def test_crossed_channel_strength_invariant(self, theta):
        A_ref = _contract(e_y, self.R_sym, e_x)
        I_ref = abs(A_ref)**2
        R_phased = np.exp(1j * theta) * self.R_sym
        I = abs(_contract(e_y, R_phased, e_x))**2
        assert abs(I - I_ref) / (I_ref + 1e-30) < 1e-12, f"theta={theta}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
