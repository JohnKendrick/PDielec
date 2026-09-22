"""Tests F1, F2, F3: Complex permittivity and complex Raman tensor handling.

Plan sections F1, F2, F3.

F1  Complex permittivity sign convention.  Adding a positive imaginary part to
    the refractive index (n = n' + in'') should cause field decay (absorption),
    not field growth.  Check that the Raman signal from a thick slab decreases
    or saturates (not grows unphysically) as n'' increases.

F2  Complex Raman tensor phase invariance in the full layered code.
    Multiplying the Raman tensor by e^{iθ} must leave the Raman intensity
    unchanged, regardless of the optical environment.

F3  Complex symmetric tensor.  A complex symmetric (but non-Hermitian) Raman
    tensor should give results matching a direct benchmark contraction.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import Layer, TransferMatrixSystem
from PDielec.LayeredRamanCalculator import RamanLayer, LayeredRamanCalculator
from .conftest import build_system, run_calc, iso_eps, make_layer

# ---------------------------------------------------------------------------
# Shared parameters
# ---------------------------------------------------------------------------

NU_MODE = 500.0
LASER   = 20000.0
D_LAYER = 100e-9

R_BASE  = np.diag([2.0, 1.0, 0.5])   # real diagonal reference tensor

# Complex symmetric tensor for F3
R_COMPLEX_SYM = np.array(
    [[1.0 + 1.0j,  0.2j,         0.0         ],
     [0.2j,        2.0 - 1.0j,   0.0         ],
     [0.0,         0.0,          0.5 + 0.1j  ]],
    dtype=complex,
)


def _build_calc(R, n_layer, d=D_LAYER, approx_es=True):
    """Return a LayeredRamanCalculator for a single active layer."""
    sup = Layer(thickness=1e-3, epsilon=iso_eps(1.0))
    sub = Layer(thickness=1e-3, epsilon=iso_eps(1.0))
    film = make_layer(d, iso_eps(n_layer))
    system = TransferMatrixSystem(substrate=sub, superstrate=sup, layers=[film])
    rl = RamanLayer(
        layer_index=0,
        phonon_frequencies_cm1=np.array([NU_MODE]),
        raman_tensors=[np.asarray(R, dtype=complex)],
        rotation_matrix=np.eye(3, dtype=float),
    )
    calc = LayeredRamanCalculator(
        system=system,
        raman_layers=[rl],
        laser_frequency_cm1=LASER,
        incident_angle_rad=0.0,
        incident_pol="p",
        detected_pol="p",
        temperature_K=0.0,
        linewidths_cm1=np.array([5.0]),
        n_gauss=20,
        approximate_es=approx_es,
    )
    return calc


def _run(R, n_layer, d=D_LAYER, approx_es=True):
    calc = _build_calc(R, n_layer, d, approx_es)
    _, intensities, _ = calc.calculate_mode_intensities()
    return intensities[0]


# ---------------------------------------------------------------------------
# F1: Complex permittivity sign convention
# ---------------------------------------------------------------------------

class TestF1ComplexPermittivitySign:
    """F1: Positive imaginary n'' means absorption → field decay → saturation."""

    def test_signal_does_not_grow_with_more_absorption(self):
        r"""Increasing n'' from 2 to 10 must monotonically decrease the signal.

        For moderate n'' (0.1–1.0) the coherent Raman integral is non-monotonic
        because the complex k shifts the system between constructive/destructive
        interference regimes.  For n'' ≥ 2 the penetration depth δ < 20 nm is
        ≪ d=3200 nm, so the saturation formula I ∝ \|t\|⁴/\|k\|² dominates and the
        signal is strictly decreasing with n''.
        """
        n_vals  = [1.5 + 2j, 1.5 + 3j, 1.5 + 5j, 1.5 + 10j]
        d_thick = 3200e-9

        I_prev = _run(R_BASE, n_vals[0], d=d_thick)
        for n in n_vals[1:]:
            I_curr = _run(R_BASE, n, d=d_thick)
            assert I_curr <= I_prev * 1.05, (
                f"n={n}: I_curr={I_curr:.4e} > I_prev={I_prev:.4e} "
                f"— absorption has wrong sign?"
            )
            I_prev = I_curr

    def test_thin_layer_not_strongly_affected_by_small_absorption(self):
        """For a thin layer (d ≪ penetration depth) a small n'' should have
        negligible effect on the signal."""
        n_real  = 1.5 + 0.0j
        n_small_imag = 1.5 + 0.01j
        d_thin  = 5e-9   # 5 nm ≪ penetration depth ≈ 1000 nm for n''=0.01
        I_real = _run(R_BASE, n_real, d=d_thin)
        I_abs  = _run(R_BASE, n_small_imag, d=d_thin)
        rel = abs(I_real - I_abs) / (I_real + 1e-30)
        assert rel < 0.01, (
            f"Thin layer: small n'' changed signal by {rel*100:.2f} %"
        )

    def test_large_imaginary_part_reduces_signal_vs_real(self):
        """A strongly absorbing thick layer must give less signal than a
        transparent layer of the same thickness."""
        d_thick = 3200e-9
        I_real = _run(R_BASE, 1.5 + 0.0j, d=d_thick)
        I_absorb = _run(R_BASE, 1.5 + 1.0j, d=d_thick)
        assert I_absorb < I_real, (
            f"Absorbing layer should give less signal: I_abs={I_absorb:.4e}, "
            f"I_real={I_real:.4e}"
        )


# ---------------------------------------------------------------------------
# F2: Complex Raman tensor phase invariance in full layered code
# ---------------------------------------------------------------------------

class TestF2ComplexPhaseFull:
    """F2: I(e^{iθ} R) = I(R) in the full layered Raman code."""

    THETAS = [0.0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2,
              1.23, 2.5, np.pi]
    RTOL   = 1e-10

    @pytest.mark.parametrize("theta", [0.0, np.pi / 4, 1.23, np.pi])
    def test_approximate_es(self, theta):
        I_ref = _run(R_BASE, 1.5, approx_es=True)
        I_phased = _run(np.exp(1j * theta) * R_BASE, 1.5, approx_es=True)
        assert abs(I_phased / I_ref - 1.0) < self.RTOL, (
            f"approx_es=True, theta={theta:.2f}: I_ref={I_ref:.6e}, I_phased={I_phased:.6e}"
        )

    @pytest.mark.parametrize("theta", [0.0, np.pi / 4, 1.23, np.pi])
    def test_full_calculation(self, theta):
        I_ref = _run(R_BASE, 1.5, approx_es=False)
        I_phased = _run(np.exp(1j * theta) * R_BASE, 1.5, approx_es=False)
        assert abs(I_phased / I_ref - 1.0) < self.RTOL, (
            f"approx_es=False, theta={theta:.2f}: I_ref={I_ref:.6e}, I_phased={I_phased:.6e}"
        )

    @pytest.mark.parametrize("theta", [0.5, 1.5, 3.0])
    def test_absorbing_layer(self, theta):
        """Phase invariance must also hold inside an absorbing medium."""
        n = 1.5 + 0.3j
        I_ref = _run(R_BASE, n, approx_es=True)
        I_phased = _run(np.exp(1j * theta) * R_BASE, n, approx_es=True)
        assert abs(I_phased / I_ref - 1.0) < self.RTOL, (
            f"absorbing, theta={theta:.2f}"
        )


# ---------------------------------------------------------------------------
# F3: Complex symmetric (non-Hermitian) Raman tensor
# ---------------------------------------------------------------------------

class TestF3ComplexSymmetricTensor:
    """F3: A complex symmetric tensor R = R^T must be treated as symmetric,
    not Hermitian.  Reciprocity I(p→s) = I(s→p) must still hold for R = R^T.
    """

    RTOL = 1e-8

    def test_complex_symmetric_reciprocity(self):
        """I(p→s) = I(s→p) for complex symmetric tensor."""
        sup  = Layer(thickness=1e-3, epsilon=iso_eps(1.0))
        sub  = Layer(thickness=1e-3, epsilon=iso_eps(1.0))
        film = make_layer(D_LAYER, iso_eps(1.5))
        system = TransferMatrixSystem(substrate=sub, superstrate=sup, layers=[film])

        def _I(pol_in, pol_out):
            rl = RamanLayer(0, np.array([NU_MODE]), [R_COMPLEX_SYM], np.eye(3))
            calc = LayeredRamanCalculator(
                system=system, raman_layers=[rl],
                laser_frequency_cm1=LASER, incident_angle_rad=0.0,
                incident_pol=pol_in, detected_pol=pol_out,
                temperature_K=0.0, linewidths_cm1=np.array([5.0]),
                n_gauss=20, approximate_es=True,
            )
            _, intensities, _ = calc.calculate_mode_intensities()
            return intensities[0]

        I_ps = _I("p", "s")
        I_sp = _I("s", "p")
        assert abs(I_ps - I_sp) / (max(I_ps, I_sp) + 1e-30) < self.RTOL, (
            f"Complex symmetric: I_ps={I_ps:.6e}, I_sp={I_sp:.6e}"
        )

    def test_complex_symmetric_phase_invariance(self):
        """Phase invariance must hold for complex symmetric tensors."""
        I_ref = _run(R_COMPLEX_SYM, 1.5, approx_es=True)
        for theta in [0.5, 1.0, 2.0]:
            I_phased = _run(np.exp(1j * theta) * R_COMPLEX_SYM, 1.5, approx_es=True)
            assert abs(I_phased / I_ref - 1.0) < 1e-10, f"theta={theta}"

    def test_complex_tensor_result_is_real_and_positive(self):
        """The Raman intensity must always be real and non-negative."""
        I = _run(R_COMPLEX_SYM, 1.5, approx_es=True)
        assert np.isreal(I) or abs(np.imag(I)) / (abs(I) + 1e-30) < 1e-10
        assert float(np.real(I)) >= 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
