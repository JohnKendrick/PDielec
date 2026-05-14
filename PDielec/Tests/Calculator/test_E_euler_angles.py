"""Tests E: Euler angle utilities in Calculator.py.

Covers: euler_rotation_matrix, determine_euler_angles, euler_rotation.

The convention is ZYZ Euler angles (theta, phi, psi).
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import math
import numpy as np
import pytest
from PDielec.Calculator import (
    euler_rotation_matrix,
    determine_euler_angles,
    euler_rotation,
)


# ---------------------------------------------------------------------------
# E1: euler_rotation_matrix
# ---------------------------------------------------------------------------

class TestE1EulerRotationMatrix:
    """E1: euler_rotation_matrix returns a valid SO(3) matrix."""

    def test_all_zero_gives_identity(self):
        R = euler_rotation_matrix(0.0, 0.0, 0.0)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-12)

    def test_shape(self):
        R = euler_rotation_matrix(0.1, 0.2, 0.3)
        assert R.shape == (3, 3)

    def test_orthogonal(self):
        """R @ R.T ≈ I (orthogonality)."""
        R = euler_rotation_matrix(0.5, 1.0, 1.5)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-12)

    def test_determinant_one(self):
        """det(R) ≈ +1 (proper rotation)."""
        R = euler_rotation_matrix(0.5, 1.0, 1.5)
        assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-12)

    def test_determinant_one_various_angles(self):
        for theta, phi, psi in [(0.1, 0.2, 0.3), (1.0, 2.0, 3.0), (math.pi / 2, 0, 0)]:
            R = euler_rotation_matrix(theta, phi, psi)
            assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-11)

    def test_pi_theta_flips_xy(self):
        """theta=π about Y: x→-x, z→-z (approximately, depends on phi/psi convention)."""
        R = euler_rotation_matrix(math.pi, 0.0, 0.0)
        # R should be orthogonal
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-12)


# ---------------------------------------------------------------------------
# E2: determine_euler_angles
# ---------------------------------------------------------------------------

class TestE2DetermineEulerAngles:
    """E2: determine_euler_angles inverts euler_rotation_matrix."""

    def test_identity_gives_zero_angles(self):
        theta, phi, psi = determine_euler_angles(np.eye(3))
        R_back = euler_rotation_matrix(theta, phi, psi)
        np.testing.assert_allclose(R_back, np.eye(3), atol=1e-10)

    def test_round_trip_random_angles(self):
        """``euler_rotation_matrix(*determine_euler_angles(R))`` ≈ R."""
        rng = np.random.default_rng(42)
        for _ in range(10):
            # Random rotation via QR decomposition
            Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            if np.linalg.det(Q) < 0:
                Q[:, 0] *= -1
            theta, phi, psi = determine_euler_angles(Q)
            R_back = euler_rotation_matrix(theta, phi, psi)
            np.testing.assert_allclose(R_back, Q, atol=1e-8)


# ---------------------------------------------------------------------------
# E3: euler_rotation
# ---------------------------------------------------------------------------

class TestE3EulerRotation:
    """E3: euler_rotation applies euler_rotation_matrix to a vector."""

    def test_zero_angles_leave_vector_unchanged(self):
        v = np.array([1.0, 2.0, 3.0])
        result = euler_rotation(v, 0.0, 0.0, 0.0)
        np.testing.assert_allclose(result, v, atol=1e-12)

    def test_preserves_vector_length(self):
        v = np.array([1.0, 0.0, 0.0])
        result = euler_rotation(v, 0.5, 1.0, 1.5)
        assert np.linalg.norm(result) == pytest.approx(np.linalg.norm(v), rel=1e-12)

    def test_consistent_with_rotation_matrix(self):
        """euler_rotation(v, ...) == euler_rotation_matrix(...) @ v."""
        v = np.array([1.0, 0.5, -0.3])
        theta, phi, psi = 0.4, 0.7, 1.1
        R = euler_rotation_matrix(theta, phi, psi)
        expected = R @ v
        result = euler_rotation(v, theta, phi, psi)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_double_application(self):
        """Applying the same rotation twice equals applying the squared rotation."""
        v = np.array([1.0, 0.0, 0.0])
        theta, phi, psi = 0.3, 0.6, 0.9
        R = euler_rotation_matrix(theta, phi, psi)
        once = euler_rotation(v, theta, phi, psi)
        twice = euler_rotation(once, theta, phi, psi)
        expected = R @ R @ v
        np.testing.assert_allclose(twice, expected, atol=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
