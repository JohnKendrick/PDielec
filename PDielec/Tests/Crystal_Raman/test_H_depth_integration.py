"""Tests for Crystal Raman depth-integration modes."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
from .conftest import build_system, make_raman_layer, run_calc

from PDielec.LayeredRamanCalculator import (
    DEPTH_INTEGRATION_COHERENT,
    DEPTH_INTEGRATION_INCOHERENT,
)

NU_MODE = 500.0
N_LAYER = 1.0
R_PP = np.diag([1.0, 0.0, 0.0])


def _intensity(thickness_m, depth_integration):
    """Return a single-layer Raman intensity for the given depth mode."""
    system = build_system([(thickness_m, N_LAYER)], n_sup=N_LAYER, n_sub=N_LAYER)
    layer = make_raman_layer(0, R_PP, NU_MODE)
    _freqs, intensities = run_calc(
        system,
        [layer],
        [NU_MODE],
        n_gauss=31,
        approximate_es=True,
        coherent_layers=False,
        depth_integration=depth_integration,
    )
    return intensities[0]


def test_incoherent_depth_scales_linearly_for_uniform_field():
    """Incoherent depth integration sums local intensity, so I ∝ thickness."""
    i1 = _intensity(1.0e-9, DEPTH_INTEGRATION_INCOHERENT)
    i2 = _intensity(2.0e-9, DEPTH_INTEGRATION_INCOHERENT)
    assert abs(i2 / i1 - 2.0) < 1.0e-10


def test_coherent_depth_scales_quadratically_for_uniform_field():
    """Coherent amplitude depth integration gives I ∝ thickness²."""
    i1 = _intensity(1.0e-9, DEPTH_INTEGRATION_COHERENT)
    i2 = _intensity(2.0e-9, DEPTH_INTEGRATION_COHERENT)
    assert abs(i2 / i1 - 4.0) < 1.0e-2


def test_incoherent_depth_is_subdivision_stable():
    """Splitting a homogeneous layer should not change incoherent depth intensity."""
    total = 100.0e-9
    one = build_system([(total, N_LAYER)], n_sup=N_LAYER, n_sub=N_LAYER)
    two = build_system([(0.5 * total, N_LAYER), (0.5 * total, N_LAYER)], n_sup=N_LAYER, n_sub=N_LAYER)
    one_layers = [make_raman_layer(0, R_PP, NU_MODE)]
    two_layers = [make_raman_layer(0, R_PP, NU_MODE), make_raman_layer(1, R_PP, NU_MODE)]
    _freqs, i_one = run_calc(
        one,
        one_layers,
        [NU_MODE],
        n_gauss=31,
        approximate_es=True,
        coherent_layers=False,
        depth_integration=DEPTH_INTEGRATION_INCOHERENT,
    )
    _freqs, i_two = run_calc(
        two,
        two_layers,
        [NU_MODE],
        n_gauss=31,
        approximate_es=True,
        coherent_layers=False,
        depth_integration=DEPTH_INTEGRATION_INCOHERENT,
    )
    assert abs(i_two[0] / i_one[0] - 1.0) < 1.0e-10
