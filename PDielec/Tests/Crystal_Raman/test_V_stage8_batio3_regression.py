"""Stage 8 numerical regressions for the BaTiO3 NAC/EO calculation."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from PDielec.NACDiagnostics import compute_nac_direction_diagnostics

FIXTURE = Path(__file__).with_name("data") / "batio3_stage8.json"
SPECTRUM_GRID_CM1 = np.linspace(100.0, 750.0, 651)
LINEWIDTH_CM1 = 4.0


@pytest.fixture(scope="module")
def batio3():
    """Load the compact PDielec inputs and independent ANADDB references."""
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def inputs(batio3):
    """Return production NAC diagnostic arguments from the fixture."""
    values = batio3["inputs"]
    return (
        np.asarray(values["hessian"]),
        np.asarray(values["born_charges"]),
        np.asarray(values["epsilon_infinity"]),
        values["volume_bohr3"],
        np.asarray(values["masses_electron_mass"]),
        np.asarray(values["to_eigenvectors_mass_weighted"]),
        np.asarray(values["to_raman_tensors_R_epsilon"]),
        np.asarray(values["linewidths_cm1"]),
    )


@pytest.fixture(scope="module")
def chi2(batio3):
    """Return the independently unit-audited internal chi(2) tensor."""
    return np.asarray(batio3["inputs"]["chi2_R_epsilon"])


def _endpoint_results(batio3, inputs, chi2):
    """Calculate both ANADDB endpoint directions using the production routine."""
    return [
        compute_nac_direction_diagnostics(endpoint["q_cartesian"], *inputs, chi2_repsilon=chi2)
        for endpoint in batio3["anaddb"]["endpoints"]
    ]


def _fibonacci_directions(number):
    """Return deterministic equal-area directions on the unit sphere."""
    indices = np.arange(number)
    z = 1.0 - 2.0 * (indices + 0.5) / number
    phi = math.pi * (3.0 - math.sqrt(5.0)) * indices
    radius = np.sqrt(1.0 - z * z)
    return np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))


def _powder_spectrum(inputs, chi2, number, eo_enabled):
    """Return a normalised deterministic q-spherical Raman reference spectrum."""
    spectrum = np.zeros_like(SPECTRUM_GRID_CM1)
    for q_hat in _fibonacci_directions(number):
        result = compute_nac_direction_diagnostics(
            q_hat, *inputs, chi2_repsilon=chi2 if eo_enabled else None
        )
        frequencies = np.asarray(result["frequencies_cm1"])
        activity_key = "raman_activities_with_eo" if eo_enabled else "raman_activities_without_eo"
        activities = np.asarray(result[activity_key])[:, 0]
        optical = frequencies > 10.0
        spectrum += np.sum(
            activities[optical, np.newaxis] * LINEWIDTH_CM1
            / (
                (SPECTRUM_GRID_CM1[np.newaxis, :] - frequencies[optical, np.newaxis]) ** 2
                + LINEWIDTH_CM1 ** 2
            ),
            axis=0,
        )
    spectrum /= number
    return spectrum / np.max(spectrum)


def test_batio3_endpoint_frequencies_and_soft_mode_ancestry(batio3, inputs, chi2):
    """Regress the E/A1 endpoints and the mixed 700 cm-1 soft-E branch."""
    perpendicular, parallel = _endpoint_results(batio3, inputs, chi2)
    frequency_tolerance = batio3["metadata"]["frequency_tolerance_cm1"]

    for target in (183.60, 298.80, 453.90, 700.34):
        assert np.min(np.abs(perpendicular["frequencies_cm1"] - target)) < frequency_tolerance
    for target in (186.63, 466.58, 693.28):
        assert np.min(np.abs(parallel["frequencies_cm1"] - target)) < frequency_tolerance

    high_branch = int(np.argmax(perpendicular["frequencies_cm1"]))
    composition = perpendicular["squared_overlap_matrix"][:, high_branch]
    assert np.sum(composition[3:5]) == pytest.approx(0.74175, abs=2.0e-4)
    assert composition[12] == pytest.approx(0.2224, abs=2.0e-4)


def test_batio3_every_endpoint_mode_is_accounted_for_once(batio3, inputs, chi2):
    """Every sorted endpoint branch has a complete unit-normalised TO composition."""
    for result in _endpoint_results(batio3, inputs, chi2):
        assert result["overlap_matrix"].shape == (15, 15)
        np.testing.assert_allclose(result["squared_overlap_sums"], np.ones(15), atol=5.0e-14)
        assert sorted(np.concatenate(result["degenerate_groups"]).tolist()) == list(range(15))


def test_batio3_to_and_nac_tensors_agree_with_anaddb(batio3, inputs, chi2):
    """Reproduce native ANADDB TO and EO-corrected NAC tensors.

    Both the mechanical and EO terms use
    ``R_epsilon = 4*pi*sqrt(Vcell)*dchi/dQ``.  The q||x and q||z endpoint
    tensors come from a fresh ANADDB 10.6.5 analysis of the merged DDB and are
    compared directly, without a fixture-only scale correction.
    """
    tolerance = batio3["metadata"]["tensor_tolerance_R_epsilon"]
    anaddb_to = batio3["anaddb"]["to"]
    anaddb_to_vectors = np.asarray(anaddb_to["eigenvectors_mass_weighted"])
    anaddb_to_tensors = np.asarray(anaddb_to["raman_tensors_R_epsilon"])

    to_overlap = anaddb_to_vectors @ inputs[5].T
    expected_to = np.einsum("ij,iab->jab", to_overlap, anaddb_to_tensors)
    np.testing.assert_allclose(expected_to, inputs[6], atol=tolerance, rtol=0.0)

    for endpoint, result in zip(batio3["anaddb"]["endpoints"], _endpoint_results(batio3, inputs, chi2)):
        pdielec_nac_vectors = np.asarray(result["eigenvectors"])
        expected_without_eo = np.einsum(
            "ij,iab->jab", anaddb_to_vectors @ pdielec_nac_vectors.T, anaddb_to_tensors
        )
        np.testing.assert_allclose(
            expected_without_eo, result["raman_tensors_without_eo"], atol=tolerance, rtol=0.0
        )

        anaddb_nac_vectors = np.asarray(endpoint["eigenvectors_mass_weighted"])
        anaddb_nac_tensors = np.asarray(endpoint["raman_tensors_R_epsilon"])
        expected_with_eo = np.einsum(
            "ij,iab->jab", anaddb_nac_vectors @ pdielec_nac_vectors.T, anaddb_nac_tensors
        )
        np.testing.assert_allclose(
            expected_with_eo,
            result["raman_tensors_with_eo"],
            atol=tolerance,
            rtol=0.0,
        )


def test_batio3_degenerate_soft_e_rotation_preserves_observables(inputs, chi2):
    """Rotate the 160.65 cm-1 E(TO) pair without changing physical NAC results."""
    q_hat = np.array([0.3, 0.4, 0.5])
    reference = compute_nac_direction_diagnostics(q_hat, *inputs, chi2_repsilon=chi2)
    angle = math.radians(37.0)
    rotation = np.eye(15)
    rotation[3:5, 3:5] = [
        [math.cos(angle), math.sin(angle)],
        [-math.sin(angle), math.cos(angle)],
    ]
    rotated = list(inputs)
    rotated[5] = rotation @ inputs[5]
    rotated[6] = np.einsum("ij,jab->iab", rotation, inputs[6])
    candidate = compute_nac_direction_diagnostics(q_hat, *rotated, chi2_repsilon=chi2)

    np.testing.assert_allclose(candidate["frequencies_cm1"], reference["frequencies_cm1"], atol=1.0e-12)
    np.testing.assert_allclose(candidate["raman_tensors_with_eo"], reference["raman_tensors_with_eo"], atol=1.0e-13)
    np.testing.assert_allclose(candidate["raman_activities_with_eo"], reference["raman_activities_with_eo"], atol=2.0e-14)


def test_batio3_consistent_to_reordering_preserves_observables(inputs, chi2):
    """Consistently reversing the TO basis must not alter the NAC spectrum or tensors."""
    q_hat = np.array([0.2, 0.7, 0.4])
    reference = compute_nac_direction_diagnostics(q_hat, *inputs, chi2_repsilon=chi2)
    permutation = np.arange(15)[::-1]
    reordered = list(inputs)
    reordered[5] = inputs[5][permutation]
    reordered[6] = inputs[6][permutation]
    reordered[7] = inputs[7][permutation]
    candidate = compute_nac_direction_diagnostics(q_hat, *reordered, chi2_repsilon=chi2)

    np.testing.assert_allclose(candidate["frequencies_cm1"], reference["frequencies_cm1"], atol=1.0e-12)
    np.testing.assert_allclose(candidate["raman_tensors_with_eo"], reference["raman_tensors_with_eo"], atol=1.0e-13)
    np.testing.assert_allclose(candidate["raman_activities_with_eo"], reference["raman_activities_with_eo"], atol=2.0e-14)


def test_batio3_q_spherical_spectrum_is_grid_converged(inputs, chi2):
    """Doubling a deterministic spherical grid changes the normalised spectrum negligibly."""
    coarse = _powder_spectrum(inputs, chi2, 256, eo_enabled=True)
    fine = _powder_spectrum(inputs, chi2, 512, eo_enabled=True)
    assert np.max(np.abs(fine - coarse)) < 3.0e-5


def test_batio3_eo_changes_powder_intensity_and_retains_700_cm1_feature(inputs, chi2):
    """EO produces a benchmarked nonzero change while retaining the high longitudinal feature."""
    without_eo = _powder_spectrum(inputs, chi2, 256, eo_enabled=False)
    with_eo = _powder_spectrum(inputs, chi2, 256, eo_enabled=True)
    assert np.max(np.abs(with_eo - without_eo)) == pytest.approx(0.4260, abs=2.0e-3)
    high_frequency = SPECTRUM_GRID_CM1 > 650.0
    high_peak = SPECTRUM_GRID_CM1[high_frequency][np.argmax(with_eo[high_frequency])]
    assert high_peak == pytest.approx(696.0, abs=2.0)
