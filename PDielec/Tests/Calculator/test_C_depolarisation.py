"""Tests C: Depolarisation matrices in Calculator.py.

Covers: initialise_sphere_depolarisation_matrix,
        initialise_plate_depolarisation_matrix,
        initialise_needle_depolarisation_matrix,
        initialise_ellipsoid_depolarisation_matrix.

All depolarisation matrices L must satisfy:
  - trace(L) == 1
  - L is symmetric
  - All eigenvalues in [0, 1]
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import (
    initialise_sphere_depolarisation_matrix,
    initialise_plate_depolarisation_matrix,
    initialise_needle_depolarisation_matrix,
    initialise_ellipsoid_depolarisation_matrix,
)


def _check_common_properties(L):
    """Assert that L is a valid depolarisation matrix."""
    assert L.shape == (3, 3), "Shape must be (3,3)"
    np.testing.assert_allclose(np.trace(L), 1.0, atol=1e-12, err_msg="trace must be 1")
    np.testing.assert_allclose(L, L.T, atol=1e-12, err_msg="L must be symmetric")
    eigs = np.linalg.eigvalsh(L)
    assert np.all(eigs >= -1e-12), f"Eigenvalues must be >= 0, got {eigs}"
    assert np.all(eigs <= 1.0 + 1e-12), f"Eigenvalues must be <= 1, got {eigs}"


# ---------------------------------------------------------------------------
# C1: initialise_sphere_depolarisation_matrix
# ---------------------------------------------------------------------------

class TestC1SpherePolarisation:
    """C1: Sphere depolarisation matrix is 1/3 * I."""

    def test_shape(self):
        L = initialise_sphere_depolarisation_matrix()
        assert L.shape == (3, 3)

    def test_diagonal_one_third(self):
        L = initialise_sphere_depolarisation_matrix()
        np.testing.assert_allclose(np.diag(L), [1.0/3, 1.0/3, 1.0/3], atol=1e-12)

    def test_off_diagonal_zero(self):
        L = initialise_sphere_depolarisation_matrix()
        off = L - np.diag(np.diag(L))
        np.testing.assert_allclose(off, 0.0, atol=1e-12)

    def test_trace_is_one(self):
        L = initialise_sphere_depolarisation_matrix()
        np.testing.assert_allclose(np.trace(L), 1.0, atol=1e-12)

    def test_common_properties(self):
        _check_common_properties(initialise_sphere_depolarisation_matrix())


# ---------------------------------------------------------------------------
# C2: initialise_plate_depolarisation_matrix
# ---------------------------------------------------------------------------

class TestC2PlatePolarisation:
    """C2: Plate depolarisation matrix concentrates all weight along the normal."""

    def test_normal_z_diagonal(self):
        """Normal along z: L[2,2] = 1, L[0,0] = L[1,1] = 0."""
        L = initialise_plate_depolarisation_matrix(np.array([0, 0, 1]))
        np.testing.assert_allclose(L[2, 2], 1.0, atol=1e-12)
        np.testing.assert_allclose(L[0, 0], 0.0, atol=1e-12)
        np.testing.assert_allclose(L[1, 1], 0.0, atol=1e-12)

    def test_normal_x_diagonal(self):
        """Normal along x: L[0,0] = 1, L[1,1] = L[2,2] = 0."""
        L = initialise_plate_depolarisation_matrix(np.array([1, 0, 0]))
        np.testing.assert_allclose(L[0, 0], 1.0, atol=1e-12)
        np.testing.assert_allclose(L[1, 1], 0.0, atol=1e-12)
        np.testing.assert_allclose(L[2, 2], 0.0, atol=1e-12)

    def test_trace_is_one(self):
        L = initialise_plate_depolarisation_matrix(np.array([0, 0, 1]))
        np.testing.assert_allclose(np.trace(L), 1.0, atol=1e-12)

    def test_symmetric(self):
        L = initialise_plate_depolarisation_matrix(np.array([1, 1, 1]))
        np.testing.assert_allclose(L, L.T, atol=1e-12)

    def test_unnormalised_normal_same_result(self):
        """An unnormalised normal gives the same L as the normalised one."""
        L1 = initialise_plate_depolarisation_matrix(np.array([0, 0, 1]))
        L2 = initialise_plate_depolarisation_matrix(np.array([0, 0, 5]))
        np.testing.assert_allclose(L1, L2, atol=1e-12)

    def test_common_properties(self):
        _check_common_properties(initialise_plate_depolarisation_matrix(np.array([0, 0, 1])))


# ---------------------------------------------------------------------------
# C3: initialise_needle_depolarisation_matrix
# ---------------------------------------------------------------------------

class TestC3NeedlePolarisation:
    """C3: Needle depolarisation matrix splits weight equally across transverse directions."""

    def test_unique_z_diagonal(self):
        """Unique direction along z: L[0,0]=L[1,1]=0.5, L[2,2]=0."""
        L = initialise_needle_depolarisation_matrix(np.array([0, 0, 1]))
        np.testing.assert_allclose(L[0, 0], 0.5, atol=1e-12)
        np.testing.assert_allclose(L[1, 1], 0.5, atol=1e-12)
        np.testing.assert_allclose(L[2, 2], 0.0, atol=1e-12)

    def test_unique_x_diagonal(self):
        """Unique direction along x: L[1,1]=L[2,2]=0.5, L[0,0]=0."""
        L = initialise_needle_depolarisation_matrix(np.array([1, 0, 0]))
        np.testing.assert_allclose(L[0, 0], 0.0, atol=1e-12)
        np.testing.assert_allclose(L[1, 1], 0.5, atol=1e-12)
        np.testing.assert_allclose(L[2, 2], 0.5, atol=1e-12)

    def test_trace_is_one(self):
        L = initialise_needle_depolarisation_matrix(np.array([0, 0, 1]))
        np.testing.assert_allclose(np.trace(L), 1.0, atol=1e-12)

    def test_symmetric(self):
        L = initialise_needle_depolarisation_matrix(np.array([1, 1, 0]))
        np.testing.assert_allclose(L, L.T, atol=1e-12)

    def test_common_properties(self):
        _check_common_properties(initialise_needle_depolarisation_matrix(np.array([0, 0, 1])))


# ---------------------------------------------------------------------------
# C4: initialise_ellipsoid_depolarisation_matrix
# ---------------------------------------------------------------------------

class TestC4EllipsoidPolarisation:
    """C4: Ellipsoid reduces to sphere (aoverb=1), and satisfies trace=1."""

    def test_aoverb_1_sphere(self):
        """aoverb == 1 should reproduce the sphere result (1/3, 1/3, 1/3)."""
        L = initialise_ellipsoid_depolarisation_matrix(np.array([0.0, 0.0, 1.0]), 1.0)
        np.testing.assert_allclose(np.trace(L), 1.0, atol=1e-12)
        # All eigenvalues should be 1/3
        eigs = sorted(np.linalg.eigvalsh(L))
        np.testing.assert_allclose(eigs, [1/3, 1/3, 1/3], atol=1e-8)

    def test_aoverb_very_large_needle_limit(self):
        """Very large aoverb (prolate) → unique axis approaches 0 (needle)."""
        L = initialise_ellipsoid_depolarisation_matrix(np.array([0.0, 0.0, 1.0]), 1000.0)
        # L[2,2] should approach 0
        assert L[2, 2] < 0.05

    def test_aoverb_very_small_plate_limit(self):
        """Very small aoverb (oblate) → unique axis approaches 1 (plate)."""
        L = initialise_ellipsoid_depolarisation_matrix(np.array([0.0, 0.0, 1.0]), 0.001)
        # L[2,2] should approach 1
        assert L[2, 2] > 0.9

    def test_trace_is_one_for_various_ratios(self):
        unique = np.array([0.0, 0.0, 1.0])
        for aoverb in [0.1, 0.5, 1.0, 2.0, 10.0]:
            L = initialise_ellipsoid_depolarisation_matrix(unique, aoverb)
            np.testing.assert_allclose(np.trace(L), 1.0, atol=1e-10,
                                       err_msg=f"trace != 1 for aoverb={aoverb}")

    def test_symmetric(self):
        L = initialise_ellipsoid_depolarisation_matrix(np.array([0.0, 0.0, 1.0]), 2.0)
        np.testing.assert_allclose(L, L.T, atol=1e-12)

    def test_common_properties_prolate(self):
        _check_common_properties(
            initialise_ellipsoid_depolarisation_matrix(np.array([0.0, 0.0, 1.0]), 2.0)
        )

    def test_common_properties_oblate(self):
        _check_common_properties(
            initialise_ellipsoid_depolarisation_matrix(np.array([0.0, 0.0, 1.0]), 0.5)
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
