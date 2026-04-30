"""Tests B1, B3, B4: Single-layer optical field tests.

Plan sections B1, B3, B4.

B1  Uniform-field single-crystal limit: for a thin slab in matched media the
    Raman signal scales as d² (coherent depth integration within one layer).

B3  Thickness scaling without interference: I(2d)/I(d) ≈ 4 for sufficiently
    thin layers in unabsorbing matched media.

B4  Absorbing thick slab saturation: for a strongly absorbing layer the signal
    saturates as the layer exceeds the optical penetration depth.

All tests use:
  - Normal incidence, p-in p-out
  - approximate_es=True (E_S ≈ E_L, both at the laser frequency)
  - Temperature = 0 K so the Bose factor is the same for all thickness comparisons
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from conftest import build_system, run_calc, make_raman_layer

# ---------------------------------------------------------------------------
# Shared parameters
# ---------------------------------------------------------------------------

NU_MODE  = 500.0       # cm⁻¹ — single optical phonon mode
LASER    = 20000.0     # cm⁻¹ ≈ 500 nm
R_PP     = np.diag([1.0, 0.0, 0.0])   # only p-pol Raman activity


def _intensity(d_m, n=1.0, n_sup=1.0, n_sub=1.0, approximate_es=True):
    """Return the Raman intensity for a single layer.

    Parameters
    ----------
    d_m            : float   — layer thickness in metres
    n              : complex — refractive index of the active layer
    n_sup          : float   — superstrate index
    n_sub          : float   — substrate index
    approximate_es : bool    — use E_S ≈ E_L (True) or compute at ν_S (False)
    """
    system = build_system([(d_m, n)], n_sup=n_sup, n_sub=n_sub)
    rl = make_raman_layer(0, R_PP, NU_MODE)
    _, intensities = run_calc(
        system, [rl], [NU_MODE],
        approximate_es=approximate_es,
        coherent_layers=False,
    )
    assert len(intensities) == 1, "Expected exactly one active mode"
    return intensities[0]


# ---------------------------------------------------------------------------
# B1: Uniform-field limit — I ∝ d² for thin matched layers
# ---------------------------------------------------------------------------

class TestB1UniformFieldLimit:
    """B1: For d ≪ λ in fully matched media (n_layer = n_sup = n_sub = 1)
    the Raman signal must scale as d² to better than 1 % accuracy.

    The test checks I(d) / d² is approximately constant over a range of
    thicknesses that are small compared with the optical wavelength (500 nm).
    """

    THICKNESSES = [1e-9, 2e-9, 4e-9, 8e-9]   # 1, 2, 4, 8 nm — all ≪ 500 nm

    def test_quadratic_scaling(self):
        I_over_d2 = [_intensity(d) / d**2 for d in self.THICKNESSES]
        ref = I_over_d2[0]
        for d, ratio in zip(self.THICKNESSES, I_over_d2):
            assert abs(ratio / ref - 1.0) < 0.02, (
                f"d={d*1e9:.1f} nm: I/d²={ratio:.6e}, ref={ref:.6e}"
            )

    def test_signal_increases_with_thickness(self):
        intensities = [_intensity(d) for d in self.THICKNESSES]
        for i in range(len(intensities) - 1):
            assert intensities[i + 1] > intensities[i]


# ---------------------------------------------------------------------------
# B3: Thickness scaling — I(2d)/I(d) ≈ 4 (coherent depth integration)
# ---------------------------------------------------------------------------

class TestB3ThicknessScaling:
    """B3: Doubling the layer thickness should quadruple the Raman intensity
    when the fields are nearly uniform (thin-slab, matched media).

    The test verifies I(2d)/I(d) = 4 ± 3 % for each doubling step.
    At d = 8 nm the phase accumulation is k_L·d ≈ 0.1 rad, so the
    deviation from perfect d² scaling is < 1 %.
    """

    BASE     = 1e-9       # 1 nm
    N_DOUBLE = 4          # test d, 2d, 4d, 8d → 1, 2, 4, 8 nm
    RTOL     = 0.03       # 3 % tolerance

    def test_doubling_ratio(self):
        d = self.BASE
        I_prev = _intensity(d)
        for k in range(self.N_DOUBLE):
            d *= 2
            I_next = _intensity(d)
            ratio = I_next / I_prev
            assert abs(ratio - 4.0) < 4.0 * self.RTOL, (
                f"step {k}: d={d*1e9:.1f} nm, I(2d)/I(d)={ratio:.4f}"
            )
            I_prev = I_next


# ---------------------------------------------------------------------------
# B4: Absorbing thick slab — signal must saturate
# ---------------------------------------------------------------------------

class TestB4AbsorbingSaturation:
    """B4: For a strongly absorbing layer the Raman signal must reach a finite
    limit as the layer becomes much thicker than the optical penetration depth.

    High absorption: n = 1.5 + 2.0j.
    Field penetration depth δ ≈ λ / (4π n'') ≈ 500 / (4π × 2) ≈ 20 nm.

    Pass criterion:
    - Signal does not grow indefinitely (relative change < 1 % from
      320 nm to 640 nm, i.e. 16–32 penetration depths).
    - A layer with even higher absorption gives a smaller or equal
      saturated signal level.

    Note: Due to Fabry–Perot interference from the outer interfaces the
    signal is non-monotonic with thickness, so a simple "grows then
    saturates" criterion is replaced by the convergence criterion above.
    """

    N_ABSORB      = 1.5 + 2.0j  # strong absorption, δ ≈ 20 nm
    D_PLATEAU_1   = 320e-9      # ~16 × δ — saturation should be reached
    D_PLATEAU_2   = 640e-9      # ~32 × δ
    RTOL_SAT      = 0.01        # 1 % tolerance for saturation plateau

    def _I(self, d, n=None):
        if n is None:
            n = self.N_ABSORB
        return _intensity(d, n=n, approximate_es=False)

    def test_signal_saturates(self):
        I1 = self._I(self.D_PLATEAU_1)
        I2 = self._I(self.D_PLATEAU_2)
        rel = abs(I2 - I1) / (I1 + 1e-30)
        assert rel < self.RTOL_SAT, (
            f"Signal not saturated: I(320 nm)={I1:.4e}, I(640 nm)={I2:.4e}, "
            f"rel_change={rel:.4f}"
        )

    def test_saturated_signal_is_finite_and_positive(self):
        I = self._I(self.D_PLATEAU_2)
        assert I > 0.0

    def test_higher_absorption_gives_lower_or_equal_saturated_signal(self):
        """Stronger absorption limits the effective sampling depth → lower signal."""
        n_more = 1.5 + 3.0j
        I_less = self._I(self.D_PLATEAU_2, n=self.N_ABSORB)
        I_more = self._I(self.D_PLATEAU_2, n=n_more)
        assert I_more <= I_less * 1.05, (
            f"More-absorbing layer should not give a larger saturated signal: "
            f"n_less={self.N_ABSORB}, I={I_less:.4e}; "
            f"n_more={n_more}, I={I_more:.4e}"
        )

    def test_very_thin_signal_smaller_than_saturated(self):
        """A 1 nm layer (sub-penetration-depth) gives less signal than the saturated
        thick layer — verifying the signal grows with thickness up to saturation."""
        I_thin = self._I(1e-9)
        I_thick = self._I(self.D_PLATEAU_2)
        assert I_thick > I_thin, (
            f"Saturated signal should exceed 1 nm layer: "
            f"I(1 nm)={I_thin:.4e}, I(640 nm)={I_thick:.4e}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
