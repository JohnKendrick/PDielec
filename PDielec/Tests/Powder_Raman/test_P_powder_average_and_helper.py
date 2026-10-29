"""Independent rotational averages and harmonic Stokes checks for scripting."""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import Boltzmann, Planck, speed_of_light

from PDielec.Calculator import compute_powder_raman_intensities
from PDielec.HelperRoutines import calculate_powder_raman_spectrum


def _reader(frequencies):
    """Return a minimal reader with a fixed diagonal Raman tensor."""
    tensors = [np.diag([1., 2., 3.]) for _ in frequencies]
    return SimpleNamespace(frequencies=frequencies, volume=100., zerof_optical_dielectric=np.eye(3),
                           get_raman_tensors=lambda: tensors)


def _z_rotation(angle):
    """Return a Cartesian rotation about z."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])


@pytest.mark.parametrize("tensor", [
    np.eye(3), np.diag([1., -2., 1.]),
    np.array([[0., 2., -1.], [-2., 0., 3.], [1., -3., 0.]]),
    np.array([[1.+2j, 3.-1j, 0.2j], [-2., 0.7j, 1.], [0.5, -1j, 4.]]),
])
def test_invariants_match_deterministic_rotation_integral(tensor):
    """Check eq-invariants1/eq-Intensities by integrating amplitudes over Euler angles."""
    # Haar measure d(alpha) d(cos beta) d(gamma)/(8 pi^2).
    # Eight periodic nodes resolve harmonics through order four exactly.
    cosines, weights = np.polynomial.legendre.leggauss(6)
    angles = 2. * np.pi * np.arange(8) / 8.
    average = np.zeros(2)
    for cosine, weight in zip(cosines, weights):
        sine = np.sqrt(1. - cosine**2)
        y_rotation = np.array([[cosine, 0., sine], [0., 1., 0.], [-sine, 0., cosine]])
        for alpha in angles:
            for gamma in angles:
                rotation = _z_rotation(alpha) @ y_rotation @ _z_rotation(gamma)
                lab = rotation @ tensor @ rotation.T
                average += weight * np.abs([lab[1, 1], lab[0, 1]])**2 / (2. * len(angles)**2)
    np.testing.assert_allclose(compute_powder_raman_intensities(tensor), 45. * average,
                               rtol=3.e-14, atol=2.e-13)


@pytest.mark.parametrize("invalid", [-400., np.nan, np.inf, -np.inf, 0., 0.5, 20000., 21000.])
def test_explicit_selection_cannot_admit_invalid_stokes_modes(invalid):
    """Invalid selected modes neither add a line nor contaminate a stable line."""
    axis = np.array([395., 400., 405.])
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        reference = calculate_powder_raman_spectrum(axis, _reader([400.]),
                                                    modes_selected=[True], polarisation="VV",
                                                    laser_wavelength_nm=500.)
        actual = calculate_powder_raman_spectrum(axis, _reader([invalid, 400.]),
                                                 modes_selected=[True, True], polarisation="VV",
                                                 laser_wavelength_nm=500.)
    np.testing.assert_array_equal(actual, reference)


@pytest.mark.parametrize("temperature", [0., 0.01, 300.])
@pytest.mark.parametrize("polarisation,activity", [("VV", 192.), ("VH", 9.), ("HV", 9.), ("Unpolarised", 201.)])
def test_helper_line_height_matches_thermal_stokes_oracle(temperature, polarisation, activity):
    """Diagonal tensor invariants and a Boltzmann population give an absolute reduced line height."""
    frequency, laser, width = 400., 20000., 5.
    population = 1. if temperature == 0 else 1. / (
        1. - np.exp(-Planck * speed_of_light * 100. * frequency / (Boltzmann * temperature)))
    expected = activity * (laser-frequency)**4 * population / frequency / width
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        actual = calculate_powder_raman_spectrum([frequency], _reader([frequency]),
                                                 modes_selected=[True], polarisation=polarisation,
                                                 laser_wavelength_nm=500., temperature=temperature, sigma=width)
    assert actual[0] == pytest.approx(expected, rel=2.e-8)


@pytest.mark.parametrize("wavelength", [0., -500., np.nan, np.inf])
def test_helper_rejects_invalid_laser_wavelength(wavelength):
    """The laser must define a finite positive incoming photon energy."""
    with pytest.raises(ValueError, match="finite.*positive"):
        calculate_powder_raman_spectrum([400.], _reader([400.]), modes_selected=[True],
                                       laser_wavelength_nm=wavelength)
