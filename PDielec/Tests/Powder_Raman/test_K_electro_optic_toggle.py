"""Powder Raman electro-optic toggle tests."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np

import PDielec.GUI.PowderScenarioTab as powder_tab
from PDielec.GUI.PowderScenarioTab import PowderScenarioTab


class _Reader:
    nions = 1
    volume = 100.0
    masses = [1.0]
    born_charges = [np.diag([1.0, 0.5, 2.0])]
    mass_weighted_normal_modes = [
        [[1.0, 0.0, 0.0]],
        [[0.0, 1.0, 0.0]],
        [[0.0, 0.0, 1.0]],
    ]
    hessian = np.diag((np.array([100.0, 150.0, 220.0]) * powder_tab.wavenumber) ** 2)


def _tab_with_reader(reader=None):
    tab = PowderScenarioTab.__new__(PowderScenarioTab)
    tab.reader = _Reader() if reader is None else reader
    return tab


def _spectrum(tab, include_eo, chi2, born_charges=None, no_matrix=False):
    if born_charges is not None:
        tab.reader.born_charges = born_charges
    raman_tensors = [
        np.diag([1.0, 0.2, 0.1]),
        np.array([[0.0, 0.7, 0.0], [0.7, 0.0, 0.1], [0.0, 0.1, 0.0]]),
        np.diag([0.3, 0.4, 1.2]),
    ]
    return tab._compute_orientation_sampled_spectrum(
        L=np.eye(3) / 3.0,
        epsilon_e=None if no_matrix else 1.5,
        epsilon_inf_i=np.diag([2.0, 2.5, 3.0]),
        I3=np.eye(3, dtype=complex),
        raman_tensors=raman_tensors,
        frequencies_cm1=np.array([100.0, 150.0, 220.0]),
        sigmas_cm1=np.array([5.0, 5.0, 5.0]),
        modes_selected=[True, True, True],
        polarisation="VV",
        nu_L=1.0e7 / 532.0,
        temperature=298.0,
        n_samples=32,
        vs_cm1=np.linspace(80.0, 240.0, 41),
        no_matrix=no_matrix,
        include_eo=include_eo,
        chi2_repsilon=chi2,
    )


def test_powder_eo_disabled_matches_zero_chi2():
    """Disabling EO must remove the separable χ² correction."""
    tab = _tab_with_reader()
    chi2 = np.zeros((3, 3, 3))
    chi2[0, 0, 2] = 4.0
    chi2[1, 1, 2] = -2.0

    disabled = _spectrum(tab, include_eo=False, chi2=chi2)
    zero = _spectrum(tab, include_eo=True, chi2=np.zeros((3, 3, 3)))

    np.testing.assert_allclose(disabled, zero, rtol=0.0, atol=1.0e-12)


def test_powder_eo_enabled_changes_nonzero_chi2_spectrum():
    """A non-zero χ² tensor changes the deterministic backscatter powder average."""
    tab = _tab_with_reader()
    chi2 = np.zeros((3, 3, 3))
    chi2[0, 0, 2] = 4.0
    chi2[1, 1, 2] = -2.0

    disabled = _spectrum(tab, include_eo=False, chi2=chi2)
    enabled = _spectrum(tab, include_eo=True, chi2=chi2)

    assert np.max(np.abs(enabled - disabled)) > 1.0e-6


def test_powder_eo_zero_born_charges_no_change():
    """The EO correction vanishes when particle-mode effective charges vanish."""
    tab = _tab_with_reader()
    chi2 = np.zeros((3, 3, 3))
    chi2[0, 0, 2] = 4.0

    zero_born = [np.zeros((3, 3))]
    disabled = _spectrum(
        tab, include_eo=False, chi2=chi2, born_charges=zero_born)
    enabled = _spectrum(
        tab,
        include_eo=True,
        chi2=chi2,
        born_charges=zero_born,
    )

    np.testing.assert_allclose(enabled, disabled, rtol=0.0, atol=1.0e-12)


def test_powder_eo_is_bypassed_when_matrix_is_none():
    """Matrix=none retains the documented no-particle-correction behavior."""
    tab = _tab_with_reader()
    chi2 = np.zeros((3, 3, 3))
    chi2[0, 0, 2] = 4.0

    disabled = _spectrum(tab, include_eo=False, chi2=chi2, no_matrix=True)
    enabled = _spectrum(tab, include_eo=True, chi2=chi2, no_matrix=True)

    np.testing.assert_allclose(enabled, disabled, rtol=0.0, atol=1.0e-12)


def test_powder_eo_sampling_is_deterministic():
    """The Sobol particle-EO model gives reproducible EO-enabled spectra."""
    tab = _tab_with_reader()
    chi2 = np.zeros((3, 3, 3))
    chi2[0, 0, 2] = 4.0

    first = _spectrum(tab, include_eo=True, chi2=chi2)
    second = _spectrum(tab, include_eo=True, chi2=chi2)

    np.testing.assert_allclose(first, second, rtol=0.0, atol=0.0)


def test_powder_particle_modes_selection_follows_dominant_to_mode(monkeypatch):
    """Particle-mode selections should inherit the dominant Settings-tab TO mode."""
    tab = _tab_with_reader()
    eig_vec = np.array([
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    eig_val = (np.array([400.0, 560.0, 700.0]) * powder_tab.wavenumber) ** 2

    monkeypatch.setattr(powder_tab.np.linalg, "eigh", lambda _matrix: (eig_val, eig_vec))

    particle_freqs, _particle_raman, particle_sigmas, particle_selected = tab._compute_particle_modes(
        N_bg=np.eye(3),
        L=np.eye(3) / 3.0,
        epsilon_e=2.0,
        epsilon_inf_i=np.diag([2.0, 2.5, 3.0]),
        I3=np.eye(3, dtype=complex),
        raman_tensors=[np.eye(3), 2.0 * np.eye(3), 3.0 * np.eye(3)],
        frequencies_cm1=np.array([100.0, 150.0, 220.0]),
        sigmas_cm1=np.array([3.0, 5.0, 7.0]),
        modes_selected=[True, False, False],
    )

    np.testing.assert_allclose(particle_freqs, [400.0, 560.0, 700.0])
    assert particle_selected == [False, True, False]
    np.testing.assert_allclose(particle_sigmas, [5.0, 3.0, 7.0])


def test_powder_eo_changes_tensors_but_not_particle_frequencies():
    """EO uses particle eigenvectors but does not alter the particle dynamical matrix."""
    tab = _tab_with_reader()
    epsilon_e = 1.5
    epsilon_inf = np.diag([2.0, 2.5, 3.0])
    depolarisation = np.eye(3) / 3.0
    n_background = powder_tab.Calculator.compute_internal_field_tensor(
        depolarisation, epsilon_inf, epsilon_e)
    tensors = [
        np.diag([1.0, 0.2, 0.1]),
        np.array([[0.0, 0.7, 0.0], [0.7, 0.0, 0.1], [0.0, 0.1, 0.0]]),
        np.diag([0.3, 0.4, 1.2]),
    ]
    chi2 = np.zeros((3, 3, 3))
    chi2[0, 0, 2] = 4.0
    chi2[1, 1, 2] = -2.0
    arguments = (
        n_background,
        depolarisation,
        epsilon_e,
        epsilon_inf,
        np.eye(3, dtype=complex),
        tensors,
        np.array([100.0, 150.0, 220.0]),
        np.array([5.0, 5.0, 5.0]),
        [True, True, True],
    )

    frequencies_off, tensors_off, _, _ = tab._compute_particle_modes(*arguments)
    frequencies_on, tensors_on, _, _ = tab._compute_particle_modes(
        *arguments, include_eo=True, chi2_repsilon=chi2)

    np.testing.assert_allclose(frequencies_on, frequencies_off, rtol=0.0, atol=0.0)
    assert max(np.max(np.abs(on - off)) for on, off in zip(tensors_on, tensors_off)) > 0.0


def test_powder_eo_orientation_selection_follows_dominant_to_mode(monkeypatch):
    """The particle EO powder path should not apply the GUI mask by sorted mode index."""
    tab = _tab_with_reader()
    eig_vec = np.array([
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    eig_val = (np.array([400.0, 560.0, 700.0]) * powder_tab.wavenumber) ** 2

    monkeypatch.setattr(powder_tab.np.linalg, "eigh", lambda _matrix: (eig_val, eig_vec))
    monkeypatch.setattr(PowderScenarioTab, "_get_sobol_rotations", staticmethod(lambda _n: [np.eye(3)]))

    spectrum = tab._compute_orientation_sampled_spectrum(
        L=np.eye(3) / 3.0,
        epsilon_e=2.0,
        epsilon_inf_i=np.diag([2.0, 2.5, 3.0]),
        I3=np.eye(3, dtype=complex),
        raman_tensors=[np.eye(3), np.eye(3), np.eye(3)],
        frequencies_cm1=np.array([100.0, 150.0, 220.0]),
        sigmas_cm1=np.array([0.1, 0.1, 0.1]),
        modes_selected=[True, False, False],
        polarisation="VV",
        nu_L=1.0e7 / 532.0,
        temperature=298.0,
        n_samples=1,
        vs_cm1=np.array([400.0, 560.0]),
        no_matrix=False,
        include_eo=True,
        chi2_repsilon=np.zeros((3, 3, 3)),
    )

    assert spectrum[1] > 1.0e4 * spectrum[0]
