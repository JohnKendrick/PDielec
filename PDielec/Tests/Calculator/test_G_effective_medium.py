"""Tests G: Effective medium theories and internal field tensor in Calculator.py.

Covers: averaged_permittivity, maxwell, maxwell_sihvola, balan,
        compute_internal_field_tensor, compute_particle_raman_tensor.

Common invariants for all EMT functions:
  - At vf == 0, result approaches the medium permittivity.
  - Result is a complex 3x3 tensor.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.Calculator import (
    initialise_sphere_depolarisation_matrix,
    averaged_permittivity,
    maxwell,
    maxwell_sihvola,
    balan,
    compute_internal_field_tensor,
    compute_particle_raman_tensor,
)


def _isotropic_eps(n_real):
    """Return isotropic 3x3 permittivity tensor for refractive index n."""
    return (n_real ** 2) * np.eye(3, dtype=complex)


_MEDIUM = _isotropic_eps(1.5)    # eps of the matrix
_CRYSTAL = _isotropic_eps(2.0)   # eps of the inclusion
_L_SPHERE = initialise_sphere_depolarisation_matrix().astype(complex)
_SHAPE = "sphere"
_SIZE = 0.0  # zero size → no size correction


# ---------------------------------------------------------------------------
# G1: averaged_permittivity
# ---------------------------------------------------------------------------

class TestG1AveragedPermittivity:
    """G1: averaged_permittivity = vf*eps_crystal + (1-vf)*eps_medium, isotropic avg."""

    def test_vf_zero_returns_medium(self):
        result = averaged_permittivity(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.0, _SIZE)
        trace_med = np.trace(_MEDIUM) / 3.0
        expected = np.array([[trace_med, 0, 0], [0, trace_med, 0], [0, 0, trace_med]])
        np.testing.assert_allclose(result.real, expected.real, atol=1e-12)

    def test_vf_one_returns_crystal_avg(self):
        result = averaged_permittivity(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 1.0, _SIZE)
        trace_crystal = np.trace(_CRYSTAL) / 3.0
        np.testing.assert_allclose(result[0, 0].real, trace_crystal.real, atol=1e-12)

    def test_shape_is_3x3(self):
        result = averaged_permittivity(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.5, _SIZE)
        assert result.shape == (3, 3)

    def test_isotropic_result_diagonal(self):
        """Isotropic inputs → diagonal output."""
        result = averaged_permittivity(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.3, _SIZE)
        off = result - np.diag(np.diag(result))
        np.testing.assert_allclose(np.abs(off), 0.0, atol=1e-12)

    def test_identical_media_trivial(self):
        """If medium == crystal, result == medium regardless of vf."""
        result = averaged_permittivity(_MEDIUM, _MEDIUM, _SHAPE, _L_SPHERE, 0.5, _SIZE)
        trace_med = np.trace(_MEDIUM) / 3.0
        np.testing.assert_allclose(result[0, 0].real, trace_med.real, atol=1e-12)


# ---------------------------------------------------------------------------
# G2: maxwell
# ---------------------------------------------------------------------------

class TestG2Maxwell:
    """G2: Maxwell-Garnett effective medium theory."""

    def test_vf_zero_returns_medium(self):
        result = maxwell(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.0, _SIZE)
        trace_med = np.trace(_MEDIUM) / 3.0
        np.testing.assert_allclose(result[0, 0].real, trace_med.real, atol=1e-8)

    def test_shape_is_3x3(self):
        result = maxwell(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.3, _SIZE)
        assert result.shape == (3, 3)

    def test_isotropic_result(self):
        result = maxwell(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.3, _SIZE)
        off = result - np.diag(np.diag(result))
        np.testing.assert_allclose(np.abs(off), 0.0, atol=1e-12)

    def test_identical_media_trivial(self):
        """When eps_crystal == eps_medium, effective eps == eps_medium."""
        result = maxwell(_MEDIUM, _MEDIUM, _SHAPE, _L_SPHERE, 0.5, _SIZE)
        trace_med = np.trace(_MEDIUM) / 3.0
        np.testing.assert_allclose(result[0, 0].real, trace_med.real, atol=1e-8)

    def test_result_between_media(self):
        """Effective eps should lie between the two component eps values."""
        vf = 0.3
        result = maxwell(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, vf, _SIZE)
        eps_eff = result[0, 0].real
        eps_med = (_MEDIUM[0, 0].real)
        eps_crys = (_CRYSTAL[0, 0].real)
        assert min(eps_med, eps_crys) <= eps_eff <= max(eps_med, eps_crys) + 1e-6


# ---------------------------------------------------------------------------
# G3: maxwell_sihvola
# ---------------------------------------------------------------------------

class TestG3MaxwellSihvola:
    """G3: Maxwell-Sihvola effective medium theory."""

    def test_vf_zero_returns_medium(self):
        result = maxwell_sihvola(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.0, _SIZE)
        trace_med = np.trace(_MEDIUM) / 3.0
        np.testing.assert_allclose(result[0, 0].real, trace_med.real, atol=1e-8)

    def test_shape_is_3x3(self):
        result = maxwell_sihvola(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.3, _SIZE)
        assert result.shape == (3, 3)

    def test_identical_media_trivial(self):
        result = maxwell_sihvola(_MEDIUM, _MEDIUM, _SHAPE, _L_SPHERE, 0.5, _SIZE)
        trace_med = np.trace(_MEDIUM) / 3.0
        np.testing.assert_allclose(result[0, 0].real, trace_med.real, atol=1e-8)


# ---------------------------------------------------------------------------
# G4: balan
# ---------------------------------------------------------------------------

class TestG4Balan:
    """G4: Balan effective medium theory."""

    def test_vf_zero_returns_near_zero(self):
        """At vf=0, Balan gives a small contribution (not the full medium eps)."""
        result = balan(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.0, _SIZE)
        # At vf=0, the balan formula gives vf * something = 0
        np.testing.assert_allclose(result[0, 0].real, 0.0, atol=1e-12)

    def test_shape_is_3x3(self):
        result = balan(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.3, _SIZE)
        assert result.shape == (3, 3)

    def test_isotropic_result(self):
        result = balan(_MEDIUM, _CRYSTAL, _SHAPE, _L_SPHERE, 0.3, _SIZE)
        off = result - np.diag(np.diag(result))
        np.testing.assert_allclose(np.abs(off), 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# G5: compute_internal_field_tensor / compute_particle_raman_tensor
# ---------------------------------------------------------------------------

class TestG5InternalFieldTensor:
    """G5: compute_internal_field_tensor and compute_particle_raman_tensor."""

    def test_no_contrast_gives_identity(self):
        """When eps_i == eps_e (scalar), the internal field tensor N = I."""
        eps_e = 2.25  # scalar
        eps_i = eps_e * np.eye(3, dtype=complex)
        L = _L_SPHERE
        N = compute_internal_field_tensor(L, eps_i, eps_e)
        np.testing.assert_allclose(N, np.eye(3, dtype=complex), atol=1e-12)

    def test_identity_internal_field_gives_crystal_tensor(self):
        """With N=I (no contrast), particle Raman tensor == crystal Raman tensor."""
        eps_e = 2.25
        eps_i = eps_e * np.eye(3, dtype=complex)
        L = _L_SPHERE
        N = compute_internal_field_tensor(L, eps_i, eps_e)
        R_crystal = np.diag([1.0, 2.0, 3.0]).astype(complex)
        R_particle = compute_particle_raman_tensor(R_crystal, N, L, eps_i, eps_e)
        np.testing.assert_allclose(R_particle, R_crystal, atol=1e-12)

    def test_internal_field_tensor_shape(self):
        N = compute_internal_field_tensor(_L_SPHERE, _CRYSTAL, 1.5 ** 2)
        assert N.shape == (3, 3)

    def test_particle_raman_tensor_shape(self):
        N = compute_internal_field_tensor(_L_SPHERE, _CRYSTAL, 1.5 ** 2)
        R = np.eye(3, dtype=complex)
        result = compute_particle_raman_tensor(R, N, _L_SPHERE, _CRYSTAL, 1.5 ** 2)
        assert result.shape == (3, 3)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
