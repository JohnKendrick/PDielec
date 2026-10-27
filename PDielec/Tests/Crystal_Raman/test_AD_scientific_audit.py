"""Independent detector-completeness and matched-slab Raman checks."""

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np
import pytest

from PDielec.LayeredRamanCalculator import (
    DEPTH_INTEGRATION_COHERENT,
    DEPTH_INTEGRATION_INCOHERENT,
    MODAL_PAIR_COHERENT_ALL,
    MODAL_PAIR_GROUP_Q,
    MODAL_PAIR_INCOHERENT,
    LayeredRamanCalculator,
    RamanLayer,
)

from .conftest import build_system


@pytest.mark.parametrize("combination", [MODAL_PAIR_GROUP_Q, MODAL_PAIR_INCOHERENT, MODAL_PAIR_COHERENT_ALL])
@pytest.mark.parametrize("polar", [False, True])
@pytest.mark.parametrize("depth", [DEPTH_INTEGRATION_COHERENT, DEPTH_INTEGRATION_INCOHERENT])
@pytest.mark.parametrize("parallel", [False, True])
@pytest.mark.parametrize("cross_sign", [-1.0, 1.0])
def test_unpolarised_matched_slab_analytic(combination, polar, depth, parallel, cross_sign):
    """Orthogonal detector channels add powers, with no relative-phase cross term."""
    laser, phonon, thickness = 20000.0, 500.0, 300e-9
    tensor = np.array([[1.0, cross_sign * 0.6, 0.0], [cross_sign * 0.6, 0.7, 0.0], [0.0, 0.0, 0.3]])

    def intensity(detector):
        system = build_system([(thickness, 1.0)], n_sup=1.0, n_sub=1.0)
        layer = RamanLayer(0, [phonon], [tensor], np.eye(3))
        calculator = LayeredRamanCalculator(
            system, [layer], laser_frequency_cm1=laser, incident_angle_rad=0.0,
            incident_pol="p", detected_pol=detector, temperature_K=0.0,
            linewidths_cm1=[5.0], n_gauss=40, approximate_es=False,
            modal_pairs=True, modal_pair_combination=combination,
            modal_pair_use_nac=[polar], depth_integration=depth,
        )
        if parallel:
            with ThreadPool(2) as pool:
                result = calculator.calculate_mode_intensities(pool=pool)
        else:
            result = calculator.calculate_mode_intensities()
        return np.sum(result[1])

    # Unit incident/detector fields in vacuum; retro-backscatter phase is
    # exp(i (k_L + k_S) z). Integrate exactly, independently of the GTM solver.
    length_angstrom = thickness * 1e10
    phase = 2 * np.pi * 100 * (2 * laser - phonon) * thickness
    if depth == DEPTH_INTEGRATION_COHERENT:
        weight = length_angstrom**2 * np.sinc(phase / (2 * np.pi))**2
    else:
        weight = length_angstrom
    prefactor = (laser - phonon)**4 / phonon
    expected = prefactor * weight * (1.0 + 0.6**2)
    p_value, s_value, unpolarised = [intensity(pol) for pol in ("p", "s", "unpolarised")]
    assert p_value + s_value == pytest.approx(expected, rel=2e-11)
    assert unpolarised == pytest.approx(expected, rel=2e-11)


@pytest.mark.parametrize("phase", [0.0, 0.7, 2 * np.pi, 13.0])
def test_uniform_domain_is_a_projection_not_a_complete_state_sum(phase):
    """A uniform collective state loses the analytically known phase variance."""
    from PDielec.DepthCorrelation import depth_weight

    nodes, weights = np.polynomial.legendre.leggauss(48)
    length = 1200.0
    z, weights = (nodes + 1) * length / 2, weights * length / 2
    field = np.array([1.0, 0.4j])[:, None] * np.exp(1j * phase * z / length)
    local = depth_weight(field, weights)
    projected = depth_weight(field, weights, np.zeros(48, dtype=int), [length])
    expected_local = 1.16 * length
    expected_projection = expected_local * np.sinc(phase / (2 * np.pi))**2
    assert local == pytest.approx(expected_local, rel=1e-12)
    assert projected == pytest.approx(expected_projection, rel=1e-12, abs=1e-12)
    mean = (field @ weights)[:, None] / length
    omitted = np.sum(np.abs(field - mean)**2 * weights)
    assert local - projected == pytest.approx(omitted, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_nac_diatomic_si_splitting_and_tensor_completeness(axis):
    """Compare a neutral two-ion crystal with its SI longitudinal oscillator."""
    from scipy.constants import c, e, epsilon_0, physical_constants

    from PDielec.Constants import amu, angs2bohr, wavenumber
    from PDielec.NACDiagnostics import compute_nac_direction_diagnostics

    mass_amu, volume_angstrom3, charge = 30.0, 80.0, 2.0
    frequencies = np.array([0, 0, 0, 300, 400, 500], dtype=float)
    eye = np.eye(3)
    modes = np.block([[eye, eye], [eye, -eye]]) / np.sqrt(2)
    hessian = modes.T @ np.diag((frequencies * wavenumber)**2) @ modes
    epsilon = np.diag([4.0, 5.0, 6.0])
    rng = np.random.default_rng(17)
    tensors = rng.normal(size=(6, 3, 3)) + 1j * rng.normal(size=(6, 3, 3))
    data = compute_nac_direction_diagnostics(
        eye[axis], hessian, [charge * eye, -charge * eye], epsilon,
        volume_angstrom3 * angs2bohr**3, np.full(2, mass_amu * amu),
        modes, tensors, np.full(6, 5.0),
    )
    mass_kg = mass_amu * physical_constants["atomic mass constant"][0]
    delta_omega_squared = 2 * (charge * e)**2 / (
        mass_kg * epsilon_0 * volume_angstrom3 * 1e-30 * epsilon[axis, axis])
    expected = frequencies[3:].copy()
    expected[axis] = np.sqrt(expected[axis]**2 + delta_omega_squared / (2 * np.pi * c * 100)**2)
    # Different CODATA releases in scipy and package constants permit small differences.
    np.testing.assert_allclose(data["frequencies_cm1"][3:], np.sort(expected), rtol=2e-7)
    assert np.sum(np.abs(data["raman_tensors_without_eo"])**2) == pytest.approx(
        np.sum(np.abs(tensors)**2), rel=2e-14)
