"""Tests B1-B3: Layer class construction and basic properties.

Plan sections B1 (instantiation), B2 (calculate_epsilon), B3 (calculate_q
and calculate_matrices).

Uses simple frequency-independent permittivity lambdas to avoid needing
any QM/MM output reader data.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import Layer


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _const_scalar_eps(n_complex):
    """Return a frequency-independent scalar permittivity function (ε = n²)."""
    eps = complex(n_complex) ** 2
    return lambda freq_cm1: eps


def _const_tensor_eps(n_complex):
    """Return a frequency-independent diagonal tensor permittivity (ε = n²·I)."""
    eps = complex(n_complex) ** 2
    return lambda freq_cm1: eps * np.eye(3, dtype=complex)


def _make_layer_with_flags(thickness_m, eps_func):
    """Create a Layer with the coherent mode flags that System.calculate_GammaStar expects."""
    layer = Layer(thickness=thickness_m, epsilon=eps_func)
    layer.inCoherentIntensity = False
    layer.inCoherentPhase = False
    layer.inCoherentAveragePhase = False
    layer.inCoherentThick = False
    layer.SMatrix = None
    return layer


# ---------------------------------------------------------------------------
# B1: Layer instantiation
# ---------------------------------------------------------------------------

class TestB1Instantiation:
    """B1: Default and explicit Layer construction."""

    def test_default_thickness(self):
        layer = Layer(epsilon=_const_scalar_eps(1.0))
        assert layer.thick == 1.0e-6

    def test_custom_thickness(self):
        layer = Layer(thickness=5.0e-7, epsilon=_const_scalar_eps(1.0))
        assert layer.thick == 5.0e-7

    def test_set_thickness(self):
        layer = Layer(epsilon=_const_scalar_eps(1.5))
        layer.set_thickness(2.0e-6)
        assert layer.thick == 2.0e-6

    def test_set_euler_stores_angles(self):
        layer = Layer(epsilon=_const_scalar_eps(1.0))
        layer.set_euler(0.1, 0.2, 0.3)
        np.testing.assert_allclose(layer.theta, 0.1)
        np.testing.assert_allclose(layer.phi, 0.2)
        np.testing.assert_allclose(layer.psi, 0.3)

    def test_is_coherent_true_for_base_layer(self):
        layer = Layer(epsilon=_const_scalar_eps(1.0))
        assert layer.is_coherent() is True

    def test_set_epsilon_stores_function(self):
        eps_fn = _const_scalar_eps(2.0)
        layer = Layer(epsilon=_const_scalar_eps(1.0))
        layer.set_epsilon(eps_fn)
        assert layer.epsilon_tensor_function is eps_fn


# ---------------------------------------------------------------------------
# B2: calculate_epsilon
# ---------------------------------------------------------------------------

class TestB2CalculateEpsilon:
    """B2: calculate_epsilon returns the correct 3×3 tensor in the lab frame."""

    def test_scalar_eps_gives_isotropic_tensor(self):
        n = 1.5
        eps_val = n**2
        layer = Layer(thickness=1e-6, epsilon=_const_scalar_eps(n))
        # frequency in Hz (internally converted to cm⁻¹)
        result = layer.calculate_epsilon(3e12)
        assert result.shape == (3, 3)
        np.testing.assert_allclose(np.diag(result).real, [eps_val, eps_val, eps_val], atol=1e-12)

    def test_tensor_eps_preserved(self):
        n = 1.5
        eps_val = n**2
        layer = Layer(thickness=1e-6, epsilon=_const_tensor_eps(n))
        result = layer.calculate_epsilon(3e12)
        np.testing.assert_allclose(np.abs(result - eps_val * np.eye(3)), 0.0, atol=1e-12)

    def test_result_is_3x3_complex(self):
        layer = Layer(thickness=1e-6, epsilon=_const_tensor_eps(1.5))
        result = layer.calculate_epsilon(1e13)
        assert result.shape == (3, 3)
        assert np.iscomplexobj(result)

    def test_constant_eps_independent_of_frequency(self):
        layer = Layer(thickness=1e-6, epsilon=_const_tensor_eps(2.0))
        r1 = layer.calculate_epsilon(1e12)
        r2 = layer.calculate_epsilon(1e15)
        np.testing.assert_array_equal(r1, r2)

    def test_zero_euler_rotation_leaves_diagonal_tensor_diagonal(self):
        # With zero Euler angles, diagonal tensor stays diagonal
        layer = Layer(thickness=1e-6, epsilon=_const_tensor_eps(1.5))
        layer.set_euler(0.0, 0.0, 0.0)
        result = layer.calculate_epsilon(1e13)
        off_diag = result - np.diag(np.diag(result))
        np.testing.assert_allclose(np.abs(off_diag), 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# B3: calculate_q and calculate_matrices
# ---------------------------------------------------------------------------

class TestB3QAndMatrices:
    """B3: calculate_q / calculate_matrices give correct shapes and symmetry."""

    def _setup_vacuum_layer(self, zeta=0.0):
        """Return a vacuum layer ready for calculate_matrices."""
        layer = Layer(thickness=1e-6, epsilon=_const_scalar_eps(1.0))
        layer.calculate_epsilon(3e14)   # set epsilon before calling matrices
        layer.calculate_matrices(zeta)
        return layer

    def test_calculate_matrices_sets_delta_shape(self):
        layer = self._setup_vacuum_layer()
        assert layer.Delta.shape == (4, 4)

    def test_calculate_q_vacuum_normal_incidence(self):
        layer = self._setup_vacuum_layer(zeta=0.0)
        layer.calculate_q()
        # For vacuum at normal incidence: qs should contain ±1 (complex)
        qs = layer.qs
        assert qs.shape == (4,)
        # The four eigenvalues should come in ±1 pairs
        magnitudes = np.abs(qs)
        np.testing.assert_allclose(magnitudes, 1.0, atol=1e-10)

    def test_calculate_q_glass_normal_incidence(self):
        n = 1.5
        layer = Layer(thickness=1e-6, epsilon=_const_scalar_eps(n))
        layer.calculate_epsilon(3e14)
        layer.calculate_matrices(0.0)
        layer.calculate_q()
        qs = layer.qs
        # Eigenvalues for isotropic medium at normal incidence are ±n
        magnitudes = np.abs(qs)
        np.testing.assert_allclose(magnitudes, n, atol=1e-9)

    def test_ai_matrix_shape(self):
        layer = Layer(thickness=1e-6, epsilon=_const_scalar_eps(1.5))
        layer.calculate_epsilon(3e14)
        layer.calculate_matrices(0.0)
        layer.calculate_q()
        layer.calculate_gamma(0.0)
        layer.calculate_ai(0.0)
        assert layer.Ai.shape == (4, 4)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
