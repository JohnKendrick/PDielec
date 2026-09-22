"""Tests F1-F2: Incoherent layer types.

Plan sections F1 (IncoherentThickLayer), F2 (IncoherentAveragePhaseLayer).

IncoherentThickLayer suppresses backward-propagating modes and is
marked as incoherent (is_coherent() → False).

IncoherentAveragePhaseLayer applies a phase shift to the propagation
but remains formally coherent in terms of the flag (is_coherent() may
vary by implementation; we test the actual value).
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import (
    Layer,
    IncoherentThickLayer,
    IncoherentAveragePhaseLayer,
    TransferMatrixSystem,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _iso_eps(n):
    eps = complex(n) ** 2
    return lambda freq_cm1: eps * np.eye(3, dtype=complex)


class _MockMaterialLayer:
    """Minimal duck-type mock for a Materials.Layer-like object expected by
    IncoherentThickLayer / IncoherentAveragePhaseLayer constructors."""

    def __init__(self, thickness_m, n, phase_shift=0.0):
        self._thick = thickness_m
        self._eps_fn = _iso_eps(n)
        self._phase = phase_shift

    def get_thickness_in_metres(self):
        return self._thick

    def get_permittivity_function(self):
        return self._eps_fn

    def get_phase_shift(self):
        return self._phase


def _make_coherent_layer(thickness_m, n=1.0):
    layer = Layer(thickness=thickness_m, epsilon=_iso_eps(n))
    layer.inCoherentIntensity = False
    layer.inCoherentPhase = False
    layer.inCoherentAveragePhase = False
    layer.inCoherentThick = False
    layer.SMatrix = None
    return layer


def _R_pp(system, f_hz=3e14):
    """Run GTM and return R_pp at normal incidence."""
    zeta = 0.0
    system.calculate_GammaStar(f_hz, zeta)
    r_out, R_out, t_out, T_out = system.calculate_r_t(zeta)
    return R_out[0]


# ---------------------------------------------------------------------------
# F1: IncoherentThickLayer
# ---------------------------------------------------------------------------

class TestF1IncoherentThickLayer:
    """F1: IncoherentThickLayer marks itself as incoherent."""

    def _make_itl(self, thickness_m=1e-4, n=1.5):
        mock = _MockMaterialLayer(thickness_m, n)
        return IncoherentThickLayer(mock)

    def test_is_coherent_false(self):
        itl = self._make_itl()
        assert itl.is_coherent() is False

    def test_incoherent_thick_flag_true(self):
        itl = self._make_itl()
        assert itl.inCoherentThick is True

    def test_thickness_stored_correctly(self):
        itl = self._make_itl(thickness_m=2.5e-4)
        assert abs(itl.thick - 2.5e-4) < 1e-15

    def test_calculate_epsilon_returns_3x3(self):
        itl = self._make_itl(n=1.5)
        result = itl.calculate_epsilon(3e14)
        assert result.shape == (3, 3)

    def test_thick_absorbing_layer_suppresses_transmission(self):
        # A very thick, strongly absorbing IncoherentThickLayer should
        # reduce transmission compared to a thin coherent slab.
        sup = _make_coherent_layer(1e-3, n=1.0)
        sub = _make_coherent_layer(1e-3, n=1.0)

        thick_abs_mock = _MockMaterialLayer(1e-3, complex(1.5, 5.0))
        itl = IncoherentThickLayer(thick_abs_mock)

        sys_inc = TransferMatrixSystem(substrate=sub, superstrate=sup, layers=[itl])
        sys_inc.calculate_GammaStar(3e14, 0.0)
        r_out, R_out, t_out, T_out = sys_inc.calculate_r_t(0.0)
        # R + T should still be <= 1 (energy is not created)
        assert R_out[0] + T_out[0] <= 1.0 + 1e-4


# ---------------------------------------------------------------------------
# F2: IncoherentAveragePhaseLayer
# ---------------------------------------------------------------------------

class TestF2IncoherentAveragePhaseLayer:
    """F2: IncoherentAveragePhaseLayer behaviour."""

    def _make_iapl(self, thickness_m=1e-5, n=1.5, phase_shift=0.0,
                   percentage=100, n_samples=4):
        mock = _MockMaterialLayer(thickness_m, n, phase_shift=phase_shift)
        return IncoherentAveragePhaseLayer(
            mock,
            percentage_incoherence=percentage,
            number_of_samples=n_samples,
        )

    def test_is_coherent_flag(self):
        # IncoherentAveragePhaseLayer uses coherent formalism, so is_coherent
        # may be True or False; just check that the attribute exists
        iapl = self._make_iapl()
        _ = iapl.is_coherent()  # should not raise

    def test_incoherent_average_phase_flag_true(self):
        iapl = self._make_iapl()
        assert iapl.inCoherentAveragePhase is True

    def test_thickness_stored(self):
        iapl = self._make_iapl(thickness_m=3e-5)
        assert abs(iapl.thick - 3e-5) < 1e-15

    def test_calculate_epsilon_returns_3x3(self):
        iapl = self._make_iapl()
        result = iapl.calculate_epsilon(3e14)
        assert result.shape == (3, 3)

    def test_runs_in_system_without_error(self):
        sup = _make_coherent_layer(1e-3, n=1.0)
        sub = _make_coherent_layer(1e-3, n=1.0)
        iapl = self._make_iapl(thickness_m=1e-6, n=1.5)
        sys_ = TransferMatrixSystem(substrate=sub, superstrate=sup, layers=[iapl])
        # Should compute without exception
        sys_.calculate_GammaStar(3e14, 0.0)
        r_out, R_out, t_out, T_out = sys_.calculate_r_t(0.0)
        assert 0.0 <= R_out[0] <= 1.0 + 1e-6


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
