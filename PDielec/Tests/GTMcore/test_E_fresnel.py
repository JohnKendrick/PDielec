"""Tests E1-E3: Full optical stack verified against Fresnel analytic formulae.

Plan sections E1 (single interface), E2 (thin transparent slab Fabry-Pérot
fringes), E3 (absorbing medium behaviour).

Uses TransferMatrixSystem with two semi-infinite-like layers (very thick) and
no finite layers, calling calculate_GammaStar + calculate_r_t at normal
incidence (zeta=0).

The Fresnel formulae for normal incidence between media n1 (superstrate) and
n2 (substrate) are:
    r = (n1 - n2) / (n1 + n2)
    R = |r|²
    T = 1 - R   (energy conservation in lossless media)
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import Layer, TransferMatrixSystem
from PDielec.Constants import speed_light_si


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _iso_eps(n):
    eps = complex(n) ** 2
    return lambda freq_cm1: eps * np.eye(3, dtype=complex)


def _make_layer(thickness_m, n=1.0):
    layer = Layer(thickness=thickness_m, epsilon=_iso_eps(n))
    layer.inCoherentIntensity = False
    layer.inCoherentPhase = False
    layer.inCoherentAveragePhase = False
    layer.inCoherentThick = False
    layer.SMatrix = None
    return layer


def _build_system(n_sup, n_sub, layer_specs=None):
    """Build a TransferMatrixSystem.

    Parameters
    ----------
    n_sup : float or complex
        Superstrate refractive index.
    n_sub : float or complex
        Substrate refractive index.
    layer_specs : list of (float, float) or None
        Each entry is (thickness_m, n) for finite layers.
    """
    sup = _make_layer(1e-3, n=n_sup)
    sub = _make_layer(1e-3, n=n_sub)
    layers = [] if layer_specs is None else [_make_layer(d, n) for d, n in layer_specs]
    return TransferMatrixSystem(substrate=sub, superstrate=sup, layers=layers)


def _rt_at_normal_incidence(system, f_hz=3e14):
    """Run the GTM calculation and return (r_out, R_out, t_out, T_out)."""
    zeta = 0.0
    system.calculate_GammaStar(f_hz, zeta)
    return system.calculate_r_t(zeta)


# ---------------------------------------------------------------------------
# E1: Normal incidence, single air/glass interface
# ---------------------------------------------------------------------------

class TestE1SingleInterface:
    """E1: Air (n=1) / glass (n=1.5) at normal incidence."""

    # Analytical values
    n1, n2 = 1.0, 1.5
    r_analytic = (n1 - n2) / (n1 + n2)   # = -0.2
    R_analytic = r_analytic ** 2           # = 0.04

    def _make_sys(self):
        return _build_system(self.n1, self.n2)

    def test_rpp_matches_fresnel(self):
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(self._make_sys())
        # r_out = [rpp, rps, rss, rsp]; at normal incidence rpp = rss = Fresnel r
        np.testing.assert_allclose(np.abs(r_out[0]), abs(self.r_analytic), atol=1e-4)

    def test_rss_matches_fresnel(self):
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(self._make_sys())
        np.testing.assert_allclose(np.abs(r_out[2]), abs(self.r_analytic), atol=1e-4)

    def test_Rpp_matches_fresnel(self):
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(self._make_sys())
        # R_out = [Rpp, Rss, Rsp, Rps]
        np.testing.assert_allclose(R_out[0], self.R_analytic, atol=1e-4)

    def test_Rss_matches_fresnel(self):
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(self._make_sys())
        np.testing.assert_allclose(R_out[1], self.R_analytic, atol=1e-4)

    def test_matched_media_zero_reflectance(self):
        # n1 == n2: no reflection
        sys_ = _build_system(n_sup=1.5, n_sub=1.5)
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(sys_)
        np.testing.assert_allclose(R_out[0], 0.0, atol=1e-6)
        np.testing.assert_allclose(R_out[1], 0.0, atol=1e-6)

    def test_energy_conservation_lossless(self):
        # R_pp + T_pp ≈ 1 for lossless media
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(self._make_sys())
        np.testing.assert_allclose(R_out[0] + T_out[0], 1.0, atol=1e-4)

    def test_cross_polarisation_near_zero(self):
        # Normal incidence, isotropic media: rps = rsp ≈ 0
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(self._make_sys())
        # R_out = [Rpp, Rss, Rsp, Rps]
        np.testing.assert_allclose(R_out[2], 0.0, atol=1e-6)
        np.testing.assert_allclose(R_out[3], 0.0, atol=1e-6)


# ---------------------------------------------------------------------------
# E2: Single transparent layer — Fabry-Pérot fringes
# ---------------------------------------------------------------------------

class TestE2FabryPerotLayer:
    """E2: A thin glass slab surrounded by air shows Fabry-Pérot oscillations."""

    n_slab = 1.5
    n_air  = 1.0
    f_base = 3e14  # Hz

    def _R_at_thickness(self, thickness_m):
        sys_ = _build_system(self.n_air, self.n_air, [(thickness_m, self.n_slab)])
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(sys_, f_hz=self.f_base)
        return R_out[0]

    def test_reflectance_varies_with_thickness(self):
        # Fabry-Pérot: R oscillates as thickness changes
        thicknesses = np.linspace(50e-9, 500e-9, 20)
        Rs = [self._R_at_thickness(d) for d in thicknesses]
        # There should be a non-trivial spread in values (oscillations present)
        assert max(Rs) - min(Rs) > 0.01

    def test_R_T_sum_to_one_for_any_thickness(self):
        for d in [100e-9, 200e-9, 300e-9]:
            sys_ = _build_system(self.n_air, self.n_air, [(d, self.n_slab)])
            r_out, R_out, t_out, T_out = _rt_at_normal_incidence(sys_, f_hz=self.f_base)
            np.testing.assert_allclose(R_out[0] + T_out[0], 1.0, atol=1e-4,
                                       err_msg=f"Energy conservation failed at d={d}")

    def test_zero_thickness_slab_gives_zero_reflectance(self):
        # A slab with zero thickness is pure air → R = 0 (no interface)
        R_thin = self._R_at_thickness(1e-12)
        np.testing.assert_allclose(R_thin, 0.0, atol=1e-6)


# ---------------------------------------------------------------------------
# E3: Absorbing medium
# ---------------------------------------------------------------------------

class TestE3AbsorbingMedium:
    """E3: Increasing absorption increases reflectance at normal incidence."""

    n_real = 1.5

    def _R_at_kappa(self, kappa):
        """Reflectance for superstrate=air, substrate=n_real+i*kappa."""
        n_complex = complex(self.n_real, kappa)
        sys_ = _build_system(n_sup=1.0, n_sub=n_complex)
        r_out, R_out, t_out, T_out = _rt_at_normal_incidence(sys_)
        return R_out[0]

    def test_increasing_kappa_increases_R(self):
        kappas = [0.0, 0.5, 1.0, 3.0, 10.0]
        Rs = [self._R_at_kappa(k) for k in kappas]
        for i in range(len(Rs) - 1):
            assert Rs[i] <= Rs[i + 1] + 1e-6, (
                f"R did not increase: R({kappas[i]}) = {Rs[i]:.4f}, "
                f"R({kappas[i+1]}) = {Rs[i+1]:.4f}"
            )

    def test_lossless_matches_analytic_fresnel(self):
        n1, n2 = 1.0, self.n_real
        R_analytic = ((n1 - n2) / (n1 + n2)) ** 2
        R_gtm = self._R_at_kappa(0.0)
        np.testing.assert_allclose(R_gtm, R_analytic, atol=1e-4)

    def test_high_absorption_R_approaches_1(self):
        R_large_k = self._R_at_kappa(1000.0)
        assert R_large_k > 0.95, f"R = {R_large_k:.4f} for large kappa, expected > 0.95"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
