"""Tests A: Tensor utilities in Calculator.py.

Covers: initialise_diagonal_tensor, initialise_complex_diagonal_tensor,
        initialise_unit_tensor, average_tensor, similarity_transform.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import (
    initialise_diagonal_tensor,
    initialise_complex_diagonal_tensor,
    initialise_unit_tensor,
    average_tensor,
    similarity_transform,
)


# ---------------------------------------------------------------------------
# A1: initialise_diagonal_tensor
# ---------------------------------------------------------------------------

class TestA1DiagonalTensor:
    """A1: initialise_diagonal_tensor returns a real 3x3 diagonal array."""

    def test_shape(self):
        t = initialise_diagonal_tensor([1.0, 2.0, 3.0])
        assert t.shape == (3, 3)

    def test_dtype_real(self):
        t = initialise_diagonal_tensor([1.0, 2.0, 3.0])
        assert np.issubdtype(t.dtype, np.floating)

    def test_diagonal_values(self):
        t = initialise_diagonal_tensor([4.0, 5.0, 6.0])
        np.testing.assert_array_equal(np.diag(t), [4.0, 5.0, 6.0])

    def test_off_diagonal_zero(self):
        t = initialise_diagonal_tensor([1.0, 2.0, 3.0])
        off = t - np.diag(np.diag(t))
        np.testing.assert_array_equal(off, np.zeros((3, 3)))

    def test_specific_values(self):
        t = initialise_diagonal_tensor([1.0, 2.0, 3.0])
        assert t[0, 0] == 1.0
        assert t[1, 1] == 2.0
        assert t[2, 2] == 3.0


# ---------------------------------------------------------------------------
# A2: initialise_complex_diagonal_tensor
# ---------------------------------------------------------------------------

class TestA2ComplexDiagonalTensor:
    """A2: initialise_complex_diagonal_tensor returns complex 3x3 diagonal."""

    def test_shape(self):
        t = initialise_complex_diagonal_tensor([1.0, 2.0, 3.0])
        assert t.shape == (3, 3)

    def test_dtype_complex(self):
        t = initialise_complex_diagonal_tensor([1.0, 2.0, 3.0])
        assert np.issubdtype(t.dtype, np.complexfloating)

    def test_diagonal_values(self):
        t = initialise_complex_diagonal_tensor([1.5, 2.5, 3.5])
        np.testing.assert_array_equal(np.diag(t), [1.5 + 0j, 2.5 + 0j, 3.5 + 0j])

    def test_off_diagonal_zero(self):
        t = initialise_complex_diagonal_tensor([1.0, 2.0, 3.0])
        off = t - np.diag(np.diag(t))
        np.testing.assert_array_equal(off, np.zeros((3, 3), dtype=complex))

    def test_zero_input(self):
        t = initialise_complex_diagonal_tensor([0.0, 0.0, 0.0])
        np.testing.assert_array_equal(t, np.zeros((3, 3), dtype=complex))


# ---------------------------------------------------------------------------
# A3: initialise_unit_tensor
# ---------------------------------------------------------------------------

class TestA3UnitTensor:
    """A3: initialise_unit_tensor returns 3x3 identity matrix."""

    def test_equals_identity(self):
        t = initialise_unit_tensor()
        np.testing.assert_array_equal(t, np.eye(3))

    def test_shape(self):
        t = initialise_unit_tensor()
        assert t.shape == (3, 3)

    def test_dtype_real(self):
        t = initialise_unit_tensor()
        assert np.issubdtype(t.dtype, np.floating)

    def test_diagonal_ones(self):
        t = initialise_unit_tensor()
        np.testing.assert_array_equal(np.diag(t), [1.0, 1.0, 1.0])


# ---------------------------------------------------------------------------
# A4: average_tensor
# ---------------------------------------------------------------------------

class TestA4AverageTensor:
    """A4: average_tensor returns a diagonal 3x3 with trace/3 on diagonal."""

    def test_identity_tensor_average_is_one(self):
        t = np.eye(3)
        result = average_tensor(t)
        np.testing.assert_allclose(np.diag(result), [1.0, 1.0, 1.0], atol=1e-12)

    def test_diagonal_abc_average(self):
        a, b, c = 3.0, 6.0, 9.0
        t = np.diag([a, b, c])
        result = average_tensor(t)
        expected = (a + b + c) / 3.0
        np.testing.assert_allclose(np.diag(result), [expected, expected, expected], atol=1e-12)

    def test_result_shape(self):
        result = average_tensor(np.eye(3))
        assert result.shape == (3, 3)

    def test_off_diagonal_zero(self):
        t = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        result = average_tensor(t)
        off = result - np.diag(np.diag(result))
        np.testing.assert_allclose(np.abs(off), 0.0, atol=1e-12)

    def test_zero_tensor_average_is_zero(self):
        result = average_tensor(np.zeros((3, 3)))
        np.testing.assert_allclose(result, np.zeros((3, 3)), atol=1e-12)

    def test_diagonal_elements_all_equal(self):
        """Diagonal elements of the averaged tensor must all be equal."""
        t = np.diag([2.0, 4.0, 6.0])
        result = average_tensor(t)
        assert result[0, 0] == result[1, 1] == result[2, 2]


# ---------------------------------------------------------------------------
# A5: similarity_transform
# ---------------------------------------------------------------------------

class TestA5SimilarityTransform:
    """A5: similarity_transform(rot, mat) == rot @ mat @ inv(rot)."""

    def test_identity_rotation_leaves_matrix_unchanged(self):
        M = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        result = similarity_transform(np.eye(3), M)
        np.testing.assert_allclose(result, M, atol=1e-12)

    def test_orthogonal_rotation_preserves_trace(self):
        """For orthogonal R, trace(R M R^-1) == trace(M)."""
        R = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], dtype=float)  # 90° about z
        M = np.diag([1.0, 2.0, 3.0])
        result = similarity_transform(R, M)
        np.testing.assert_allclose(np.trace(result), np.trace(M), atol=1e-12)

    def test_permutation_matrix_permutes_diagonal(self):
        """Cyclic permutation matrix rotates the diagonal elements cyclically."""
        # P maps (x,y,z) -> (y,z,x)
        P = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)
        M = np.diag([1.0, 2.0, 3.0])
        result = similarity_transform(P, M)
        np.testing.assert_allclose(np.diag(result), [2.0, 3.0, 1.0], atol=1e-12)

    def test_result_shape(self):
        result = similarity_transform(np.eye(3), np.eye(3))
        assert result.shape == (3, 3)

    def test_double_rotation_composes(self):
        """Applying same rotation twice: R(R M R^-1)R^-1 = R^2 M R^-2."""
        R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)  # 90° about z
        M = np.diag([1.0, 2.0, 3.0])
        step1 = similarity_transform(R, M)
        step2 = similarity_transform(R, step1)
        expected = similarity_transform(R @ R, M)
        np.testing.assert_allclose(step2, expected, atol=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
