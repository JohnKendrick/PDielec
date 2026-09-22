"""Tests B: Geometry functions in Calculator.py.

Covers: calculate_distance, calculate_angle, calculate_torsion.

Note: calculate_angle and calculate_torsion return values in DEGREES.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import math
import numpy as np
import pytest
from PDielec.Calculator import calculate_distance, calculate_angle, calculate_torsion


# ---------------------------------------------------------------------------
# B1: calculate_distance
# ---------------------------------------------------------------------------

class TestB1Distance:
    """B1: calculate_distance returns Euclidean distance between two points."""

    def test_identical_points_distance_zero(self):
        assert calculate_distance([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == pytest.approx(0.0)

    def test_along_x_axis(self):
        assert calculate_distance([0, 0, 0], [3, 0, 0]) == pytest.approx(3.0)

    def test_pythagorean_triple_3_4_5(self):
        assert calculate_distance([0, 0, 0], [3, 4, 0]) == pytest.approx(5.0)

    def test_3d_unit_diagonal(self):
        assert calculate_distance([0, 0, 0], [1, 1, 1]) == pytest.approx(math.sqrt(3.0))

    def test_symmetric(self):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        assert calculate_distance(a, b) == pytest.approx(calculate_distance(b, a))

    def test_negative_coords(self):
        assert calculate_distance([0, 0, 0], [-3, -4, 0]) == pytest.approx(5.0)

    def test_numpy_arrays(self):
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([1.0, 0.0, 0.0])
        assert calculate_distance(a, b) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# B2: calculate_angle
# ---------------------------------------------------------------------------

class TestB2Angle:
    """B2: calculate_angle returns the a-b-c bond angle in DEGREES."""

    def test_right_angle(self):
        # b at origin, a along x, c along y → 90°
        a = [1, 0, 0]
        b = [0, 0, 0]
        c = [0, 1, 0]
        assert calculate_angle(a, b, c) == pytest.approx(90.0, abs=1e-10)

    def test_collinear_180(self):
        a = [-1, 0, 0]
        b = [0, 0, 0]
        c = [1, 0, 0]
        assert calculate_angle(a, b, c) == pytest.approx(180.0, abs=1e-8)

    def test_equilateral_60(self):
        # Equilateral triangle: vertices at (0,0,0), (1,0,0), (0.5, sqrt(3)/2, 0)
        b = [0.0, 0.0, 0.0]
        a = [1.0, 0.0, 0.0]
        c = [0.5, math.sqrt(3) / 2, 0.0]
        assert calculate_angle(a, b, c) == pytest.approx(60.0, abs=1e-8)

    def test_result_in_range(self):
        a = [1, 0, 0]
        b = [0, 0, 0]
        c = [0.5, 0.5, 0.0]
        angle = calculate_angle(a, b, c)
        assert 0.0 <= angle <= 180.0

    def test_returns_degrees_not_radians(self):
        # A right angle is 90° not π/2 ≈ 1.57
        a = [1, 0, 0]
        b = [0, 0, 0]
        c = [0, 1, 0]
        angle = calculate_angle(a, b, c)
        assert angle > 1.0  # definitely degrees, not radians


# ---------------------------------------------------------------------------
# B3: calculate_torsion
# ---------------------------------------------------------------------------

class TestB3Torsion:
    """B3: calculate_torsion returns the dihedral angle a-b-c-d in DEGREES."""

    def test_trans_configuration_180(self):
        # trans: a and d on opposite sides of the b-c bond
        a = [1.0,  1.0, 0.0]
        b = [0.0,  0.0, 0.0]
        c = [1.0,  0.0, 0.0]
        d = [2.0, -1.0, 0.0]
        angle = calculate_torsion(a, b, c, d)
        assert abs(angle) == pytest.approx(180.0, abs=1e-8)

    def test_cis_configuration_0(self):
        # cis: a and d on same side
        a = [0.0, 1.0, 0.0]
        b = [0.0, 0.0, 0.0]
        c = [1.0, 0.0, 0.0]
        d = [1.0, 1.0, 0.0]
        angle = calculate_torsion(a, b, c, d)
        assert abs(angle) == pytest.approx(0.0, abs=1e-8)

    def test_result_in_range(self):
        a = [0.0, 1.0, 0.0]
        b = [0.0, 0.0, 0.0]
        c = [1.0, 0.0, 0.0]
        d = [2.0, 0.0, 1.0]
        angle = calculate_torsion(a, b, c, d)
        assert -180.0 <= angle <= 180.0

    def test_returns_degrees_not_radians(self):
        # A 90° dihedral should return ~90, not ~1.57
        a = [0.0, 1.0, 0.0]
        b = [0.0, 0.0, 0.0]
        c = [1.0, 0.0, 0.0]
        d = [1.0, 0.0, 1.0]  # perpendicular plane
        angle = calculate_torsion(a, b, c, d)
        # 90° in radians would be ~1.57; the result should be around ±90
        assert abs(angle) > 1.0

    def test_numpy_input(self):
        a = np.array([0.0, 1.0, 0.0])
        b = np.array([0.0, 0.0, 0.0])
        c = np.array([1.0, 0.0, 0.0])
        d = np.array([2.0, 0.0, 0.0])
        # All four in same plane → 180° or 0° dihedral
        angle = calculate_torsion(a, b, c, d)
        assert abs(angle) == pytest.approx(180.0, abs=1e-6) or abs(angle) == pytest.approx(0.0, abs=1e-6)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
