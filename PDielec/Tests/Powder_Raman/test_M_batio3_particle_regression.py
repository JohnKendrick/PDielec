"""BaTiO3 finite-particle powder Raman regression tests."""

import json
from pathlib import Path

import numpy as np
import pytest

from PDielec import Calculator
from PDielec.Constants import amu, angs2bohr
from PDielec.GUI.PowderScenarioTab import PowderScenarioTab

BATIO3_FIXTURE = (
    Path(__file__).parents[1] / "Crystal_Raman" / "data" / "batio3_stage8.json"
)


@pytest.fixture(scope="module")
def particle_case():
    """Build the air-embedded spherical BaTiO3 case from the compact PEAD fixture."""
    fixture = json.loads(BATIO3_FIXTURE.read_text(encoding="utf-8"))
    inputs = fixture["inputs"]

    reader = type("BaTiO3Reader", (), {})()
    reader.nions = 5
    reader.volume = inputs["volume_bohr3"] / angs2bohr**3
    reader.masses = (np.asarray(inputs["masses_electron_mass"]) / amu).tolist()
    reader.born_charges = inputs["born_charges"]
    reader.mass_weighted_normal_modes = np.asarray(
        inputs["to_eigenvectors_mass_weighted"]
    ).reshape(15, 5, 3).tolist()
    reader.hessian = np.asarray(inputs["hessian"])

    tab = PowderScenarioTab.__new__(PowderScenarioTab)
    tab.reader = reader
    depolarisation = np.eye(3) / 3.0
    epsilon_infinity = np.asarray(inputs["epsilon_infinity"])
    epsilon_external = 1.0
    internal_field = Calculator.compute_internal_field_tensor(
        depolarisation, epsilon_infinity, epsilon_external)
    arguments = (
        internal_field,
        depolarisation,
        epsilon_external,
        epsilon_infinity,
        np.eye(3, dtype=complex),
        np.asarray(inputs["to_raman_tensors_R_epsilon"]),
        np.asarray(fixture["anaddb"]["to"]["frequencies_cm1"]),
        np.asarray(inputs["linewidths_cm1"]),
        [False, False, False, True, True, True, True, True, True, False, True, True, True, True, True],
    )
    return tab, arguments, internal_field, np.asarray(inputs["chi2_R_epsilon"])


def test_batio3_particle_eo_preserves_frequencies_and_high_surface_mode(particle_case):
    """EO changes tensors only, retaining the 627 cm-1 spherical particle mode."""
    tab, arguments, _, chi2 = particle_case

    without_eo = tab._compute_particle_modes(*arguments)
    with_eo = tab._compute_particle_modes(
        *arguments, include_eo=True, chi2_repsilon=chi2)

    np.testing.assert_allclose(with_eo[0], without_eo[0], rtol=0.0, atol=0.0)
    high_mode = int(np.argmax(without_eo[0]))
    assert without_eo[0][high_mode] == pytest.approx(627.00501147, abs=2.0e-6)


def test_batio3_corrected_high_mode_modal_intensities(particle_case):
    """Regress the two-factor local-field and finite-particle EO high-mode intensities."""
    tab, arguments, internal_field, chi2 = particle_case
    without_eo = tab._compute_particle_modes(*arguments)
    with_eo = tab._compute_particle_modes(
        *arguments, include_eo=True, chi2_repsilon=chi2)
    high_mode = int(np.argmax(without_eo[0]))

    effective_without_eo = Calculator.compute_effective_raman_tensor(
        without_eo[1][high_mode], internal_field, internal_field)
    effective_with_eo = Calculator.compute_effective_raman_tensor(
        with_eo[1][high_mode], internal_field, internal_field)
    intensity_without_eo = Calculator.compute_powder_raman_intensities(effective_without_eo)
    intensity_with_eo = Calculator.compute_powder_raman_intensities(effective_with_eo)

    np.testing.assert_allclose(
        intensity_without_eo, [0.5285216048721, 0.01249086629222], rtol=2.0e-10)
    np.testing.assert_allclose(
        intensity_with_eo, [0.2985764829650, 0.004745719012318], rtol=2.0e-10)
    assert intensity_with_eo[0] != pytest.approx(intensity_without_eo[0])


def test_batio3_particle_modes_are_repeatable(particle_case):
    """The compact BaTiO3 particle calculation is exactly repeatable."""
    tab, arguments, _, chi2 = particle_case

    first = tab._compute_particle_modes(
        *arguments, include_eo=True, chi2_repsilon=chi2)
    second = tab._compute_particle_modes(
        *arguments, include_eo=True, chi2_repsilon=chi2)

    np.testing.assert_allclose(first[0], second[0], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(first[1], second[1], rtol=0.0, atol=0.0)
