"""Independent discrete-source and analytic checks of normalized depth weights."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from PDielec.DepthCorrelation import depth_weight
from PDielec.LayeredRamanCalculator import LayeredRamanCalculator

from .conftest import build_system, make_raman_layer


def grid(edges, order=24):
    """Integrate each physical domain separately."""
    x, w = np.polynomial.legendre.leggauss(order)
    z, weights, ids = [], [], []
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        half = (right - left) / 2
        z.extend(left + half + half * x)
        weights.extend(half * w)
        ids.extend([index] * order)
    return np.array(z), np.array(weights), np.array(ids)


@pytest.mark.parametrize("domains", [1, 3, 17])
def test_constant_source_collective_normalization(domains):
    """A normalized cell collective coordinate scales as sqrt(N), not N."""
    edges = np.linspace(0, 12, domains + 1)
    z, w, ids = grid(edges)
    f = np.full((1, len(z)), 2 + 3j)
    # Discrete dipole derivative sqrt(Vcell) R, N=A*L/Vcell.
    area, cell_volume, cells = 5.0, 0.25, 240
    u = np.ones(cells) / np.sqrt(cells)
    discrete = cell_volume / area * abs(np.sum((2 + 3j) * u)) ** 2
    assert_allclose(depth_weight(f, w, ids, np.diff(edges)), discrete)
    assert_allclose(depth_weight(f, w), discrete)


@pytest.mark.parametrize("q", [0.0, 0.3, 4.2])
@pytest.mark.parametrize("domains", [1, 7])
def test_plane_wave_analytic(q, domains):
    """Normalized domain response is L sinc(q L/2)^2, summed over domains."""
    edges = np.linspace(0, 3, domains + 1)
    z, w, ids = grid(edges, 40)
    f = np.exp(1j * q * z)[None, :]
    lengths = np.diff(edges)
    expected = np.sum(lengths * np.sinc(q * lengths / (2 * np.pi)) ** 2)
    assert_allclose(depth_weight(f, w, ids, lengths), expected, atol=1e-14)


def test_local_limit_and_quadrature_convergence():
    """Domain refinement and numerical integration are independent limits."""
    values = []
    for count in (1, 8, 64, 512):
        z, w, ids = grid(np.linspace(0, 2, count + 1), 8)
        values.append(depth_weight(np.exp(2j * z)[None, :], w, ids, np.full(count, 2 / count)))
    assert np.all(np.diff(values) > 0)
    assert abs(values[-1] - 2) < 1.1e-5
    results = []
    for order in (2, 4, 12, 24):
        z, w, ids = grid([0, 2], order)
        results.append(depth_weight(np.exp((0.3 + 2j) * z)[None, :], w, ids, [2]))
    exact = abs(np.expm1((0.3 + 2j) * 2) / (0.3 + 2j)) ** 2 / 2
    assert abs(results[-1] - exact) < 1e-14
    assert abs(results[1] - exact) < abs(results[0] - exact)
    assert_allclose(results[-2:], exact, atol=1e-13)


def test_positivity_channels_and_length_units():
    """The PSD projector is bounded by the local norm; units scale linearly."""
    rng = np.random.default_rng(314)
    z, w, ids = grid([0, 1, 4])
    f = rng.normal(size=(2, len(z))) + 1j * rng.normal(size=(2, len(z)))
    value = depth_weight(f, w, ids, [1, 3])
    assert 0 <= value <= depth_weight(f, w)
    assert_allclose(value, sum(depth_weight(c[None, :], w, ids, [1, 3]) for c in f))
    assert_allclose(depth_weight(f, w * 1e-10, ids, np.array([1, 3]) * 1e-10), value * 1e-10)
    assert depth_weight(np.empty((1, 0)), [], np.array([], dtype=int), [0]) == 0


@pytest.mark.parametrize("weights,ids,lengths", [([-1], [0], [1]), ([1], [1], [1]),
                                                ([1], [0], [0]), ([1], [0], [2]),
                                                ([1], [0.0], [1])])
def test_invalid_domains(weights, ids, lengths):
    """Reject negative measures, unnormalised domains and invalid labels."""
    with pytest.raises(ValueError):
        depth_weight([[1]], weights, ids, lengths)


def calculator(thicknesses, **kwargs):
    """Transparent matched film with one fixed tensor and no NAC."""
    return LayeredRamanCalculator(
        system=build_system([(d, 1.0) for d in thicknesses]),
        raman_layers=[make_raman_layer(i, np.diag([1, 1, 0]), 500) for i in range(len(thicknesses))],
        laser_frequency_cm1=20000, incident_angle_rad=0, incident_pol="p",
        detected_pol="p", temperature_K=0, linewidths_cm1=[5], n_gauss=24, **kwargs,
    )


@pytest.mark.parametrize("side", ["superstrate", "substrate"])
@pytest.mark.parametrize("domains", [None, [0, 2.3e-7, 6e-7], [0, 6e-7]])
def test_optical_subdivision_invariance(side, domains):
    """Physical domains survive an optical interface inserted inside a domain."""
    one = calculator([6e-7], collection_side=side)
    split = calculator([1.7e-7, 4.3e-7], collection_side=side)
    a = one.calculate_depth_correlated_weights(domains)
    b = split.calculate_depth_correlated_weights(domains)
    assert_allclose(a, b, rtol=1e-11, atol=1e-10)
    # Exact vacuum propagation phases, with independent Stokes wavevector.
    k_l = 2 * np.pi * 20000 * 100
    k_s = 2 * np.pi * 19500 * 100
    q = k_l + k_s if side == "superstrate" else k_l - k_s
    lengths = np.diff(domains) if domains is not None else np.array([6e-7])
    expected = 6000 if domains is None else np.sum(lengths * 1e10 * np.sinc(q * lengths / (2*np.pi))**2)
    assert_allclose(a[1], [expected], rtol=1e-11)


def test_physical_independent_domains_change_interference():
    """Changing physical domains differs from merely subdividing optical layers."""
    calc = calculator([6e-7])
    independent = calc.calculate_depth_correlated_weights([0, 3e-7, 6e-7])[1]
    collective = calc.calculate_depth_correlated_weights([0, 6e-7])[1]
    assert not np.allclose(independent, collective)


def test_zero_thickness_and_unsupported_model():
    """Zero volume is zero weight; modal grouping requires separate physics."""
    calc = calculator([0])
    assert_allclose(calc.calculate_depth_correlated_weights([0, 0])[1], [0])
    calc = calculator([1e-8], modal_pairs=True)
    with pytest.raises(ValueError, match="total fields"):
        calc.calculate_depth_correlated_weights()


def test_legacy_conversion_and_detector_channels():
    """Only the derived normalization changes; Stokes factors are applied once."""
    from PDielec.RamanSpectrum import stokes_prefactor

    calc = calculator([1e-7], approximate_es=True)
    prefactor = stokes_prefactor(500, 20000, 0)
    legacy = calc.calculate_mode_intensities()[1]
    normalized = calc.calculate_depth_correlated_weights([0, 1e-7])[1]
    assert_allclose(normalized * 1000 * prefactor, legacy, rtol=1e-12)
    calc.depth_integration = "Incoherent intensity"
    legacy_local = calc.calculate_mode_intensities()[1]
    local = calc.calculate_depth_correlated_weights()[1]
    assert_allclose(local * prefactor, legacy_local, rtol=1e-12)
    calc._detected_jones = None
    assert_allclose(calc.calculate_depth_correlated_weights()[1], local, rtol=1e-12)


def test_invalid_optical_domain_configuration():
    """Unsupported optical averages and incomplete physical domains are explicit."""
    calc = calculator([1e-8])
    with pytest.raises(ValueError, match="domain edges"):
        calc.calculate_depth_correlated_weights([0, 0.5e-8])
    calc.system.layers[0].inCoherentIntensity = True
    with pytest.raises(ValueError, match="coherent optical"):
        calc.calculate_depth_correlated_weights()
