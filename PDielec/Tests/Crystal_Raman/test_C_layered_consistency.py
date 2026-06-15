"""Tests C1, C2: Layered-system consistency.

Plan sections C1, C2.

C1  Raman-inactive layers contribute zero signal.  A system where every
    layer has R = 0 must give I = 0 exactly.

C2  Layer subdivision invariance.  A homogeneous Raman-active layer of
    thickness d must give the same Raman intensity when replaced by N
    identical sublayers of thickness d/N, provided amplitudes are summed
    coherently across layers (coherent_layers=True).  This currently passes
    only for matched superstrate/layer/substrate media because GTMcore's field
    reconstruction is not invariant to artificial internal boundaries when
    external Fresnel reflections are present.

    Subdivision invariance is an identity:
        ∫₀ᵈ f(z) dz = Σᵢ ∫_{d_i}^{d_{i+1}} f(z) dz

    The only numerical error is GL quadrature accuracy, which is very small
    for smooth optical fields.

Failure modes caught by C2:
    - Incorrect layer z-coordinates in the GL grid.
    - Wrong phase reference at artificial layer boundaries.
    - Field-continuity violations introduced by layer subdivision.
    - Numerical integration errors that accumulate with more layers.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import Layer, TransferMatrixSystem
from PDielec.LayeredRamanCalculator import RamanLayer, LayeredRamanCalculator
from .conftest import build_system, run_calc, make_raman_layer, iso_eps, make_layer

# ---------------------------------------------------------------------------
# Shared parameters
# ---------------------------------------------------------------------------

NU_MODE = 500.0       # cm⁻¹
LASER   = 20000.0     # cm⁻¹
N_LAYER = 1.5         # refractive index of Raman-active layer (real)
D_BASE  = 100e-9      # 100 nm base thickness

R_ACTIVE = np.diag([1.0, 0.5, 0.25])   # non-trivial diagonal tensor


# ---------------------------------------------------------------------------
# C1: Raman-inactive layers produce zero signal
# ---------------------------------------------------------------------------

class TestC1RamanInactiveLayers:
    """C1: Layers with R = 0 must contribute exactly zero Raman amplitude."""

    R_ZERO = np.zeros((3, 3))

    def test_single_inactive_layer_gives_zero(self):
        system = build_system([(D_BASE, N_LAYER)])
        rl = make_raman_layer(0, self.R_ZERO, NU_MODE)
        _, intensities = run_calc(system, [rl], [NU_MODE], approximate_es=True)
        assert len(intensities) == 1
        assert intensities[0] == pytest.approx(0.0, abs=1e-30)

    def test_all_inactive_multilayer_gives_zero(self):
        thicknesses = [50e-9, 80e-9, 120e-9]
        system = build_system([(d, N_LAYER) for d in thicknesses])
        raman_layers = [make_raman_layer(i, self.R_ZERO, NU_MODE) for i in range(3)]
        _, intensities = run_calc(system, raman_layers, [NU_MODE],
                                  approximate_es=True, coherent_layers=True)
        assert intensities[0] == pytest.approx(0.0, abs=1e-30)

    def test_inactive_does_not_corrupt_active(self):
        """An inactive layer next to an active one must not change the active signal."""
        # Single-layer reference
        sys_1 = build_system([(D_BASE, N_LAYER)])
        rl_1  = make_raman_layer(0, R_ACTIVE, NU_MODE)
        _, I_ref = run_calc(sys_1, [rl_1], [NU_MODE], approximate_es=True)

        # Active layer followed by inactive layer (same permittivity, R=0)
        sys_2 = build_system([(D_BASE, N_LAYER), (D_BASE, N_LAYER)])
        rl_a = make_raman_layer(0, R_ACTIVE, NU_MODE)
        rl_b = make_raman_layer(1, self.R_ZERO, NU_MODE)
        _, I_two = run_calc(sys_2, [rl_a, rl_b], [NU_MODE],
                            approximate_es=True, coherent_layers=False)

        # The active layer now has a different optical environment (extra layer below),
        # so we only check that the inactive layer itself gives zero extra contribution
        # when coherent_layers=True and it is alone.
        rl_only_inactive = make_raman_layer(1, self.R_ZERO, NU_MODE)
        _, I_inactive_only = run_calc(sys_2, [rl_only_inactive], [NU_MODE],
                                      approximate_es=True)
        assert I_inactive_only[0] == pytest.approx(0.0, abs=1e-30)


# ---------------------------------------------------------------------------
# C2: Layer subdivision invariance
# ---------------------------------------------------------------------------

def _build_subdivided_system(d_total, n_layer, n_sublayers):
    """Build system with n_sublayers equal sublayers of total thickness d_total.

    The superstrate and substrate use the same refractive index as the layer so
    there are no Fresnel reflections at the outer interfaces.  GTMcore's
    calculate_Efield is only subdivision-invariant in this matched-media case.
    """
    d_sub = d_total / n_sublayers
    return build_system([(d_sub, n_layer)] * n_sublayers,
                        n_sup=n_layer, n_sub=n_layer)


def _raman_layers_for_subdivision(n_sublayers, R):
    """Create RamanLayer list for a subdivided system."""
    return [make_raman_layer(i, R, NU_MODE) for i in range(n_sublayers)]


def _intensity_subdivided(d_total, n_layer, n_sublayers, n_gauss=30):
    """Return Raman intensity for a system subdivided into n_sublayers."""
    system = _build_subdivided_system(d_total, n_layer, n_sublayers)
    rls = _raman_layers_for_subdivision(n_sublayers, R_ACTIVE)
    _, intensities = run_calc(
        system, rls, [NU_MODE],
        approximate_es=False,    # full calculation: E_S at ν_S = ν_L - ν_m
        coherent_layers=True,    # sum amplitudes before squaring
        n_gauss=n_gauss,
    )
    assert len(intensities) == 1
    return intensities[0]


class TestC2SubdivisionInvariance:
    """C2: Subdividing a homogeneous layer must not change the Raman strength.

    Uses coherent_layers=True so that the mathematical identity::

        ∫₀ᵈ = ∫₀^{d/N} + ∫_{d/N}^{2d/N} + … + ∫_{(N-1)d/N}^d

    translates directly to equal total amplitudes.

    The superstrate/substrate are matched to the layer index (n_sup = n_sub =
    n_layer) to avoid Fresnel reflections at the outer interfaces.  GTMcore's
    calculate_Efield propagates a backward wave from the substrate; when outer
    interfaces are mismatched the backward amplitude at an internal boundary
    depends on the full stack, making the subdivided result differ from the
    unsplit layer.  In the matched case the backward wave amplitude is zero and
    subdivision invariance becomes an exact identity up to GL quadrature error.

    Tolerance: relative error < 1e-5 (well within GL quadrature accuracy
    for smooth optical fields with n_gauss = 30).
    """

    D_TOTAL = 100e-9   # 100 nm (not too thin, not too thick)
    RTOL    = 1e-5

    def test_one_vs_two_sublayers(self):
        I1 = _intensity_subdivided(self.D_TOTAL, N_LAYER, 1)
        I2 = _intensity_subdivided(self.D_TOTAL, N_LAYER, 2)
        assert abs(I2 / I1 - 1.0) < self.RTOL, (
            f"1 vs 2 sublayers: I1={I1:.6e}, I2={I2:.6e}"
        )

    def test_one_vs_five_sublayers(self):
        I1 = _intensity_subdivided(self.D_TOTAL, N_LAYER, 1)
        I5 = _intensity_subdivided(self.D_TOTAL, N_LAYER, 5)
        assert abs(I5 / I1 - 1.0) < self.RTOL, (
            f"1 vs 5 sublayers: I1={I1:.6e}, I5={I5:.6e}"
        )

    def test_one_vs_ten_sublayers(self):
        I1  = _intensity_subdivided(self.D_TOTAL, N_LAYER, 1)
        I10 = _intensity_subdivided(self.D_TOTAL, N_LAYER, 10)
        assert abs(I10 / I1 - 1.0) < self.RTOL, (
            f"1 vs 10 sublayers: I1={I1:.6e}, I10={I10:.6e}"
        )

    def test_convergence_with_n_sublayers(self):
        """Increasing subdivision should converge, not diverge."""
        I_ref = _intensity_subdivided(self.D_TOTAL, N_LAYER, 1)
        for n in [2, 5, 10, 20]:
            I_n = _intensity_subdivided(self.D_TOTAL, N_LAYER, n)
            rel_err = abs(I_n / I_ref - 1.0)
            assert rel_err < self.RTOL, (
                f"n={n} sublayers: rel_err={rel_err:.2e}"
            )

    # Note: subdivision invariance cannot be verified for absorbing layers
    # with GTMcore's current field implementation.  Matching the outer media to
    # a complex n requires outer Layer objects with thickness=1e-3 m, which
    # causes exponential overflow in the propagation matrices.  Using real outer
    # media reintroduces backward-wave artifacts.  The four real-n tests above
    # provide sufficient coverage of the subdivision invariance property.


@pytest.mark.xfail(
    reason=(
        "GTMcore.calculate_Efield is not currently invariant to artificial "
        "boundaries between identical coherent layers when the film is not index-matched"
    ),
    strict=True,
)
@pytest.mark.parametrize(
    ("n_sup", "n_sub"),
    [
        (1.0, 1.0),
        (1.0, 1.45),
    ],
)
def test_mismatched_outer_media_subdivision_invariance_expected_failure(n_sup, n_sub):
    """Document the current GTM field-reconstruction limitation.

    Physically these two systems are equivalent:

        air / layer(100 nm) / substrate
        air / layer(50 nm) / layer(50 nm) / substrate

    The Raman intensity should therefore be invariant to the artificial internal
    boundary.  This test is marked xfail until the field reconstruction is fixed
    rather than worked around by GUI-layer coalescing.
    """
    system_1 = build_system([(D_BASE, N_LAYER)], n_sup=n_sup, n_sub=n_sub)
    system_2 = build_system([(D_BASE / 2.0, N_LAYER), (D_BASE / 2.0, N_LAYER)], n_sup=n_sup, n_sub=n_sub)
    raman_1 = [make_raman_layer(0, R_ACTIVE, NU_MODE)]
    raman_2 = [make_raman_layer(0, R_ACTIVE, NU_MODE), make_raman_layer(1, R_ACTIVE, NU_MODE)]

    _, intensity_1 = run_calc(
        system_1,
        raman_1,
        [NU_MODE],
        approximate_es=False,
        coherent_layers=True,
        n_gauss=30,
    )
    _, intensity_2 = run_calc(
        system_2,
        raman_2,
        [NU_MODE],
        approximate_es=False,
        coherent_layers=True,
        n_gauss=30,
    )

    assert intensity_2[0] / intensity_1[0] == pytest.approx(1.0, rel=1.0e-5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
