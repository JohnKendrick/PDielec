"""Tests for Crystal Raman depth-integration modes."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest

from PDielec.LayeredRamanCalculator import (
    DEPTH_INTEGRATION_COHERENT,
    DEPTH_INTEGRATION_INCOHERENT,
    DEPTH_INTEGRATION_NORMALISED_INCOHERENT,
    LayeredRamanCalculator,
)

from .conftest import build_system, make_raman_layer, run_calc

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


@pytest.mark.parametrize("modal", [False, True])
@pytest.mark.parametrize("worker", [False, True])
@pytest.mark.parametrize("detected", ["p", "unpolarised"])
@pytest.mark.parametrize("thickness", [0.0, 1e-9, 2e-7])
def test_normalised_depth_uniform_field(thickness, detected, worker, modal):
    """A uniform single layer has unit overlap at any positive thickness."""
    from .test_W_review_regressions import InlinePool

    calc = LayeredRamanCalculator(
        build_system([(thickness, 1.0)]),
        [make_raman_layer(0, R_PP, NU_MODE)],
        laser_frequency_cm1=20000.0, incident_angle_rad=0.0,
        incident_pol="p", detected_pol=detected, temperature_K=0.0,
        linewidths_cm1=[5.0], n_gauss=21, approximate_es=True,
        depth_integration=DEPTH_INTEGRATION_NORMALISED_INCOHERENT,
        modal_pairs=modal, modal_pair_use_nac=[False],
    )
    _, intensity, _ = calc.calculate_mode_intensities(pool=InlinePool() if worker else None)
    expected = (20000.0 - NU_MODE)**4 / NU_MODE if thickness > 0 else 0.0
    np.testing.assert_allclose(np.sum(intensity), expected, rtol=1e-10, atol=0.0)
    _, detailed, _, contributions = calc.calculate_mode_intensities(return_contributions=True)
    np.testing.assert_allclose(sum(c.intensity for c in contributions), np.sum(detailed), rtol=1e-10)


@pytest.mark.parametrize("modal", [False, True])
@pytest.mark.parametrize("worker", [False, True])
def test_normalised_depth_averages_each_layer_independently(modal, worker):
    """Unequal, unmerged layers contribute averages, rather than a stack average."""
    from .test_W_review_regressions import InlinePool

    calc = LayeredRamanCalculator(
        build_system([(1e-9, 1.0), (3e-9, 1.0)]),
        [make_raman_layer(0, R_PP, NU_MODE), make_raman_layer(1, 2 * R_PP, NU_MODE)],
        laser_frequency_cm1=20000.0, incident_angle_rad=0.0,
        incident_pol="p", detected_pol="p", temperature_K=0.0,
        linewidths_cm1=[5.0], n_gauss=21, approximate_es=True,
        depth_integration=DEPTH_INTEGRATION_NORMALISED_INCOHERENT,
        modal_pairs=modal, modal_pair_use_nac=[False],
    )
    _, intensity, _ = calc.calculate_mode_intensities(pool=InlinePool() if worker else None)
    np.testing.assert_allclose(intensity, [5 * (20000.0 - NU_MODE)**4 / NU_MODE], rtol=1e-10)


def test_normalised_absorbing_layer_is_raw_integral_divided_by_thickness():
    """Geometric thickness, rather than optical penetration, sets the denominator."""
    calc = LayeredRamanCalculator(
        build_system([(2e-7, 2.0 + 0.5j)]), [make_raman_layer(0, R_PP, NU_MODE)],
        laser_frequency_cm1=20000.0, incident_angle_rad=0.0,
        incident_pol="p", detected_pol="p", temperature_K=0.0,
        linewidths_cm1=[5.0], n_gauss=61,
        depth_integration=DEPTH_INTEGRATION_INCOHERENT,
    )
    _, raw, _ = calc.calculate_mode_intensities()
    calc.depth_integration = DEPTH_INTEGRATION_NORMALISED_INCOHERENT
    _, averaged, _ = calc.calculate_mode_intensities()
    np.testing.assert_allclose(averaged, raw / 2000.0, rtol=1e-12)
