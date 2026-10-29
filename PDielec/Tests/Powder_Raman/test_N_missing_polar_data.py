"""Powder scattering remains available without polar phonon-response data."""

from types import SimpleNamespace

import numpy as np
import pytest

from PDielec.GUI.PowderScenarioTab import PowderScenarioTab
from PDielec.RamanSpectrum import stokes_prefactor


def powder_tab(shape, frequency=100.0, matrix="vacuum", missing="born"):
    """Create the non-widget calculation state with one Raman-active TO mode."""
    tab = PowderScenarioTab.__new__(PowderScenarioTab)
    tab.settings = {
        "Legend": "fallback test",
        "Matrix": matrix,
        "Volume fraction": 1.0,
        "Raman laser wavelength": 500.0,
        "Raman laser polarisation": "VV",
        "Raman temperature": 0.0,
        "Raman orientation samples": 32,
    }
    tab.reader = SimpleNamespace(
        get_raman_tensors=lambda: [np.eye(3)],
        born_charges=[] if missing == "born" else [np.eye(3)],
        mass_weighted_normal_modes=[[[1, 0, 0]]],
    )
    if missing == "born":
        tab.reader.hessian = np.eye(3)
    tab.notebook = SimpleNamespace(
        plottingTab=object(),
        settingsTab=SimpleNamespace(
            frequencies_cm1=[frequency],
            sigmas_cm1=[5.0],
            modes_selected=[True],
            settings={"Optical permittivity": 2 * np.eye(3)},
        ),
        progressbars_update=lambda **kw: None,
    )
    tab.matrixMaterial = SimpleNamespace(get_optical_permittivity=lambda: 1.0)
    tab.calculate_depolarisation_tensor = lambda: np.diag(shape)
    tab.calculation_required = True
    return tab


@pytest.mark.parametrize("shape", [[1 / 3] * 3, [0.2, 0.3, 0.5]])
@pytest.mark.parametrize("missing", ["born", "hessian"])
@pytest.mark.parametrize("matrix", ["vacuum", "none"])
def test_fallback_matches_diagonal_ellipsoid_optics(shape, missing, matrix):
    """Analytic local fields and rotational invariants supply a nonzero physical oracle."""
    tab = powder_tab(shape, matrix=matrix, missing=missing)
    tab._calculate_raman(np.array([100.0]))
    diagonal = np.ones(3) if matrix == "none" else 1 / (1 + np.array(shape)) ** 2
    alpha = diagonal.mean()
    gamma2 = 1.5 * np.sum((diagonal - alpha) ** 2)
    expected = (45 * alpha**2 + 4 * gamma2) * stokes_prefactor(100.0, 20000.0, 0.0) / 5
    assert tab.raman_spectrum[0] == pytest.approx(expected, rel=1e-12)
    assert expected > 0


@pytest.mark.parametrize("frequency", [-100.0, 0.0, 0.999, np.nan, np.inf, 20000.0, 21000.0])
def test_powder_fallback_excludes_invalid_frequencies(frequency):
    """Unstable and impossible Stokes modes cannot generate a spectrum."""
    tab = powder_tab([0.2, 0.3, 0.5], frequency)
    tab._calculate_raman(np.array([100.0]))
    assert tab.raman_spectrum == [0.0]


@pytest.mark.parametrize("wavelength_nm", [500.0, 785.0, 1064.0])
def test_powder_laser_wavelength_controls_stokes_weight(wavelength_nm):
    """The wavelength setting is converted from nm to excitation wavenumber."""
    tab = powder_tab([1 / 3] * 3, matrix="none")
    tab.settings["Raman laser wavelength"] = wavelength_nm
    tab._calculate_raman(np.array([100.0]))
    expected = 45.0 * stokes_prefactor(100.0, 1.0e7 / wavelength_nm, 0.0) / 5.0
    assert tab.raman_spectrum[0] == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("wavelength_nm", [0.0, -1.0, np.nan, np.inf])
def test_powder_invalid_laser_wavelength_is_rejected(wavelength_nm):
    """Reject invalid wavelengths before dividing to obtain the wavenumber."""
    tab = powder_tab([1 / 3] * 3, matrix="none")
    tab.settings["Raman laser wavelength"] = wavelength_nm
    with pytest.raises(ValueError, match="laser wavelength"):
        tab._calculate_raman(np.array([100.0]))


def test_missing_born_fallback_matches_zero_charge_orientation_average():
    """With zero mode charges, numerical orientation averaging reduces to optical invariants."""
    from PDielec.Constants import wavenumber

    tab = powder_tab([0.2, 0.3, 0.5])
    tab._calculate_raman(np.array([100.0]))
    expected = tab.raman_spectrum[0]
    tab.reader.nions = 1
    tab.reader.volume = 100.0
    tab.reader.masses = [10.0]
    tab.reader.hessian = np.diag((np.array([100.0, 200.0, 300.0]) * wavenumber) ** 2)
    tab.reader.mass_weighted_normal_modes = np.eye(3).reshape(3, 1, 3)
    tab.reader.born_charges = np.zeros((1, 3, 3))
    tensors = [np.eye(3), np.zeros((3, 3)), np.zeros((3, 3))]
    values = tab._compute_orientation_sampled_spectrum(
        np.diag([0.2, 0.3, 0.5]),
        1.0,
        2 * np.eye(3),
        np.eye(3),
        tensors,
        [100.0, 200.0, 300.0],
        [5.0] * 3,
        [True, False, False],
        "VV",
        20000.0,
        0.0,
        4096,
        np.array([100.0]),
        include_eo=False,
    )
    assert values[0] == pytest.approx(expected, rel=2e-3)


@pytest.mark.parametrize("frequency", [-100.0, 0.0, 21000.0])
def test_orientation_path_excludes_invalid_stokes_modes(frequency):
    """The sampled path follows the same validity policy even without particle corrections."""
    tab = powder_tab([0.2, 0.3, 0.5])
    values = tab._compute_orientation_sampled_spectrum(
        np.diag([0.2, 0.3, 0.5]),
        1.0,
        2 * np.eye(3),
        np.eye(3),
        [np.eye(3)],
        [frequency],
        [5.0],
        [True],
        "VV",
        20000.0,
        0.0,
        8,
        np.array([100.0]),
        no_matrix=True,
    )
    np.testing.assert_array_equal(values, [0.0])


@pytest.mark.parametrize("polarisation", ["VV", "VH", "Unpolarised"])
def test_particle_ellipsoid_uses_analytic_average(polarisation):
    """Corrected particle spectra agree with quadrature without a sample setting."""
    from PDielec.Constants import wavenumber

    tab = powder_tab([0.2, 0.3, 0.5])
    del tab.settings["Raman orientation samples"]
    tab.settings["Raman laser polarisation"] = polarisation
    tab.settings["Raman electro-optic term"] = True
    tab.reader.nions = 1
    tab.reader.volume = 100.0
    tab.reader.masses = [10.0]
    frequencies = np.array([100.0, 200.0, 300.0])
    tab.reader.hessian = np.diag((frequencies * wavenumber) ** 2)
    tab.reader.mass_weighted_normal_modes = np.eye(3).reshape(3, 1, 3)
    tab.reader.born_charges = [np.diag([0.1, 0.2, 0.3])]
    tab.reader.nonlinear_optical_susceptibility = np.ones((3, 3, 3)) * 0.1
    tensors = [np.diag([1.0, 0.2, 0.3]), np.ones((3, 3)), np.diag([0.1, 0.4, 1.0])]
    tab.reader.get_raman_tensors = lambda: tensors
    settings = tab.notebook.settingsTab
    settings.frequencies_cm1 = frequencies
    settings.sigmas_cm1 = [5.0] * 3
    settings.modes_selected = [True] * 3
    axis = np.linspace(80, 400, 161)
    tab._calculate_raman(axis)
    analytic = np.array(tab.raman_spectrum)
    sampled = tab._compute_orientation_sampled_spectrum(
        np.diag([0.2, 0.3, 0.5]), 1.0, 2 * np.eye(3), np.eye(3),
        tensors, frequencies, [5.0] * 3, [True] * 3,
        polarisation, 20000.0, 0.0, 4096, axis,
        include_eo=True, chi2_repsilon=tab.reader.nonlinear_optical_susceptibility,
    )
    np.testing.assert_allclose(analytic, sampled, rtol=0.01)
    # Legacy scripts may still supply this key, but it cannot alter the result.
    tab.settings["Raman orientation samples"] = 0
    tab.calculation_required = True
    tab._calculate_raman(axis)
    np.testing.assert_array_equal(tab.raman_spectrum, analytic)
