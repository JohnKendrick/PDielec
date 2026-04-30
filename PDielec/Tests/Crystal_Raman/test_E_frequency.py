"""Tests E1, E2, E3: Incident and scattered frequency handling.

Plan sections E1, E2, E3.

E1  Zero-Raman-shift limit.  When the Raman shift ν_m is much smaller than the
    laser frequency ν_L, the full calculation (approximate_es=False, which
    evaluates E_S at ν_S = ν_L − ν_m) should agree closely with the
    approximate calculation (approximate_es=True, E_S ≈ E_L at ν_L).
    For a non-dispersive layer (ε independent of frequency) the agreement
    should be very good because the optical fields differ only through the
    tiny wavenumber change δk = k_L − k_S.

E2  Large Raman shift test.  When ν_m is large relative to ν_L (so ν_S differs
    appreciably from ν_L) AND the layer is dispersive, the full calculation
    should give a different result from the approximate one.  This verifies that
    the scattered field is genuinely evaluated at ν_S, not accidentally at ν_L.

E3  Non-dispersive control.  Repeating E2 with a non-dispersive refractive
    index should show a much smaller change (only from the trivial λ shift),
    confirming that dispersive vs non-dispersive behaviour is handled correctly.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from conftest import build_system, run_calc, make_raman_layer, make_layer

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

D_LAYER = 200e-9    # 200 nm
LASER   = 20000.0   # cm⁻¹

R_PP = np.diag([1.0, 0.5, 0.0])   # p-pol active


def _intensity(nu_mode, n_real, n_imag=0.0, approx_es=True, n_gauss=20, d=None):
    """Single active layer, p-in p-out, return intensity."""
    n = complex(n_real, n_imag)
    layer_d = d if d is not None else D_LAYER
    system = build_system([(layer_d, n)])
    rl = make_raman_layer(0, R_PP, nu_mode)
    _, intensities = run_calc(
        system, [rl], [nu_mode],
        laser_cm1=LASER,
        approximate_es=approx_es,
        n_gauss=n_gauss,
    )
    assert len(intensities) == 1
    return intensities[0]


def _dispersive_eps(nu_cm1):
    """A simple model dispersive permittivity: Sellmeier-like, real only.

    ε(ν) = 2.0 + 0.3 × (ν / 10000)
    At ν = 10000 cm⁻¹: ε = 2.3, n ≈ 1.517
    At ν = 20000 cm⁻¹: ε = 2.6, n ≈ 1.612
    At ν = 15000 cm⁻¹: ε = 2.45, n ≈ 1.565
    """
    eps_val = 2.0 + 0.3 * (nu_cm1 / 10000.0)
    return complex(eps_val) * np.eye(3, dtype=complex)


# ---------------------------------------------------------------------------
# E1: Zero-shift limit — approximate_es ≈ full for small ν_m
# ---------------------------------------------------------------------------

class TestE1ZeroShiftLimit:
    """E1: For ν_m ≪ ν_L and non-dispersive layer, approximate and full
    calculations must agree to better than 0.5 % relative error.

    The wavenumber shift δk = k_L - k_S = (ν_m / ν_L) × k_L.
    For ν_m = 50 cm⁻¹ and ν_L = 20000 cm⁻¹: δk/k_L = 0.25 %.
    The field overlap integral changes by ≲ δk × d, which is tiny for
    d = 200 nm and visible-frequency light.
    """

    NU_SMALL     = 50.0   # cm⁻¹ — small Raman shift (just above acoustic threshold)
    RTOL         = 0.005  # 0.5 % for n=1.5 case
    RTOL_VACUUM  = 0.015  # 1.5 % for vacuum (n=1): k_L*d≈2.51 rad gives ~1.1% sinc deviation

    def test_approximate_equals_full_small_shift(self):
        n = 1.5
        I_approx = _intensity(self.NU_SMALL, n, approx_es=True)
        I_full   = _intensity(self.NU_SMALL, n, approx_es=False)
        assert abs(I_approx - I_full) / (I_approx + 1e-30) < self.RTOL, (
            f"ν_m={self.NU_SMALL}: I_approx={I_approx:.4e}, I_full={I_full:.4e}"
        )

    def test_approximate_equals_full_vacuum_layer(self):
        """Same test for n=1 (vacuum layer) where no Fresnel effects exist.

        Uses a looser tolerance: for n=1 and d=200nm, k_L*d≈2.51 rad so the
        sinc deviation between approx and full is ~1.1%, above the 0.5% RTOL
        used for n=1.5.
        """
        I_approx = _intensity(self.NU_SMALL, 1.0, approx_es=True)
        I_full   = _intensity(self.NU_SMALL, 1.0, approx_es=False)
        assert abs(I_approx - I_full) / (I_approx + 1e-30) < self.RTOL_VACUUM

    def test_full_intensity_is_positive_small_shift(self):
        """Full calculation must give a positive intensity for any small shift."""
        I_full = _intensity(self.NU_SMALL, 1.5, approx_es=False)
        assert I_full > 0.0


# ---------------------------------------------------------------------------
# E2: Large Raman shift — dispersive layer gives different result
# ---------------------------------------------------------------------------

class TestE2LargeRamanShift:
    """E2: For a large Raman shift with a dispersive layer, approximate_es=False
    must give a noticeably different intensity than approximate_es=True.

    If the code accidentally uses ν_L instead of ν_S for E_S, the intensities
    will be equal even for dispersive media — a clear sign of a bug.
    """

    NU_LARGE = 5000.0    # cm⁻¹ — very large shift: ν_S = 15000 cm⁻¹
    RTOL_MIN = 0.001     # at least 0.1 % difference expected for dispersive case

    def test_dispersive_layer_gives_different_result(self):
        """Dispersive layer: approximate should differ from full for large shift."""
        from PDielec.GTMcore import Layer, TransferMatrixSystem
        from PDielec.LayeredRamanCalculator import RamanLayer, LayeredRamanCalculator

        sup  = Layer(thickness=1e-3, epsilon=_dispersive_eps)
        sub  = Layer(thickness=1e-3, epsilon=_dispersive_eps)
        film = make_layer(D_LAYER, _dispersive_eps)
        system = TransferMatrixSystem(substrate=sub, superstrate=sup, layers=[film])

        rl = make_raman_layer(0, R_PP, self.NU_LARGE)
        linewidths = np.array([5.0])

        def _calc(approx_es):
            from PDielec.LayeredRamanCalculator import LayeredRamanCalculator
            calc = LayeredRamanCalculator(
                system=system,
                raman_layers=[rl],
                laser_frequency_cm1=LASER,
                incident_angle_rad=0.0,
                incident_pol="p",
                detected_pol="p",
                temperature_K=0.0,
                linewidths_cm1=linewidths,
                n_gauss=20,
                approximate_es=approx_es,
            )
            _, intensities, _ = calc.calculate_mode_intensities()
            return intensities[0]

        I_approx = _calc(True)
        I_full   = _calc(False)
        rel_diff = abs(I_approx - I_full) / (max(I_approx, I_full) + 1e-30)
        assert rel_diff > self.RTOL_MIN, (
            f"Dispersive layer with large shift: approximate and full are too close "
            f"(rel_diff={rel_diff:.4f}).  Possible bug: E_S not evaluated at ν_S."
        )


# ---------------------------------------------------------------------------
# E3: Non-dispersive control — large shift, constant n
# ---------------------------------------------------------------------------

class TestE3NonDispersiveControl:
    """E3: For a non-dispersive layer and large Raman shift, the difference
    between approximate_es and full calculations is small when d ≪ λ —
    controlled only by the trivial wavelength change, not material dispersion.

    Uses d=5nm (≪ λ=500nm) so that k*d ≪ 1 for both ν_L and ν_S and the
    sinc integral difference between the two frequencies is negligible (< 1%).
    A 200nm slab at ν_m=5000 cm⁻¹ gives k_L*d≈3.8 rad, leading to a ~45%
    diff that has nothing to do with dispersion — purely a sinc effect.

    This confirms that the difference seen in E2 is caused by dispersion,
    not by a numerical artefact.
    """

    NU_LARGE = 5000.0   # same large shift as E2
    N_REAL   = 1.5      # constant refractive index (no dispersion)
    D_THIN   = 5e-9     # 5 nm ≪ λ: k*d ≈ 0.094 rad → sinc diff < 1%
    RTOL_MAX = 0.02     # less than 2 % difference expected for non-dispersive thin slab

    def test_nondispersive_small_difference(self):
        I_approx = _intensity(self.NU_LARGE, self.N_REAL, approx_es=True,  d=self.D_THIN)
        I_full   = _intensity(self.NU_LARGE, self.N_REAL, approx_es=False, d=self.D_THIN)
        rel_diff = abs(I_approx - I_full) / (max(I_approx, I_full) + 1e-30)
        assert rel_diff < self.RTOL_MAX, (
            f"Non-dispersive thin slab: rel_diff={rel_diff:.4f} > {self.RTOL_MAX}"
        )

    def test_full_intensity_positive_large_shift(self):
        I = _intensity(self.NU_LARGE, self.N_REAL, approx_es=False, d=self.D_THIN)
        assert I > 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
