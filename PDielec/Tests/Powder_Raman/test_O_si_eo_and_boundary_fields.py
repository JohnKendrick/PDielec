"""Independent SI and electrostatic boundary checks for powder Raman tensors."""

import numpy as np
import pytest
from scipy.constants import atomic_mass, elementary_charge, epsilon_0

from PDielec.Calculator import compute_effective_raman_tensor, compute_internal_field_tensor
from PDielec.Constants import amu
from PDielec.GenericOutputReader import chi2_pm_per_v_to_repsilon
from PDielec.RamanPolarCalculator import apply_eo_correction, apply_particle_eo_correction, build_Z_mat


@pytest.mark.parametrize("volume", [25.0, 125.0, 500.0])
@pytest.mark.parametrize("host", [1.0, 2.25])
@pytest.mark.parametrize("bulk", [False, True])
def test_eo_conversion_matches_si_phonon_field(volume, host, bulk):
    """Compare the reader-to-EO chain with SI polarization and induced field."""
    masses = np.array([16.0, 65.0])
    charge = np.array([[1.2, 0.2, -0.1], [0.2, 1.7, 0.3], [-0.1, 0.3, 2.1]])
    born = np.array([charge, -charge])
    mode = np.array([1.0, 2.0, -1.0, -2.0, 0.5, 1.5])
    mode /= np.linalg.norm(mode)
    epsilon = np.array([[4.0, 0.4, 0.2], [0.4, 5.0, -0.3], [0.2, -0.3, 6.0]])
    shape = np.diag([0.15, 0.3, 0.55])
    direction = np.array([1.0, -2.0, 3.0])
    direction /= np.linalg.norm(direction)
    chi_pm = np.arange(27.0).reshape(3, 3, 3) - 10.0
    chi_pm = (chi_pm + chi_pm.swapaxes(0, 1)) / 2.0

    # SI unit-cell coordinate: displacements = u Q / sqrt(M in kg).
    coordinate = 1.0e-25
    displacements = mode.reshape(2, 3) * coordinate / np.sqrt(masses[:, None] * atomic_mass)
    polarization = sum(elementary_charge * born[a] @ displacements[a] for a in range(2))
    polarization /= volume * 1.0e-30
    if bulk:
        field = -direction * (direction @ polarization) / (epsilon_0 * (direction @ epsilon @ direction))
    else:
        # Boundary equation: (host I + L (epsilon-host I)) E = -L P/epsilon_0.
        field = np.linalg.solve(host * np.eye(3) + shape @ (epsilon - host * np.eye(3)),
                                -shape @ polarization / epsilon_0)
    delta_epsilon = 2.0 * np.einsum("ijl,l->ij", chi_pm * 1.0e-12, field)
    expected_si = np.sqrt(volume * 1.0e-30) * delta_epsilon / coordinate

    charges_au = build_Z_mat(born, masses * amu)
    chi_internal = chi_pm * chi2_pm_per_v_to_repsilon(volume)
    if bulk:
        actual = apply_eo_correction([np.zeros((3, 3))], chi_internal, direction,
                                    charges_au, mode[:, None], epsilon)[0]
    else:
        kernel = compute_internal_field_tensor(shape, epsilon, host) @ shape / host
        actual = apply_particle_eo_correction([np.zeros((3, 3))], chi_internal, kernel,
                                             charges_au, mode[:, None])[0]
    actual_si = actual * np.sqrt(1.0e-10 / atomic_mass)
    # Independent constants tables differ slightly in their atomic-mass values.
    np.testing.assert_allclose(actual_si, expected_si, rtol=2.0e-8, atol=1.0e-12)


@pytest.mark.parametrize("loss", [0.0, 0.7])
@pytest.mark.parametrize("host", [1.0, 2.25])
def test_local_fields_match_polarizability_finite_difference(loss, host):
    """Differentiate the ellipsoid boundary solution, without using N in the oracle."""
    shape = np.diag([0.1, 0.3, 0.6])
    epsilon = np.array([[4.0, 0.4, 0.2], [0.4, 5.0, -0.3], [0.2, -0.3, 6.0]], dtype=complex)
    epsilon += 1j * loss * np.diag([1.0, 2.0, 3.0])
    raman = np.array([[1.0, 0.2, -0.4], [0.2, 0.5, 0.3], [-0.4, 0.3, -0.7]])

    def dipole_response(permittivity):
        contrast = permittivity - host * np.eye(3)
        # Solve for excess polarization under three independent external fields.
        return np.linalg.solve(np.eye(3) + contrast @ shape / host, contrast)

    step = 1.0e-5
    expected = (dipole_response(epsilon + step * raman)
                - dipole_response(epsilon - step * raman)) / (2.0 * step)
    internal = compute_internal_field_tensor(shape, epsilon, host)
    actual = compute_effective_raman_tensor(raman, internal)
    np.testing.assert_allclose(actual, expected, rtol=2.0e-8, atol=2.0e-10)


@pytest.mark.parametrize("epsilon", [4.0, 4.0 + 0.7j])
@pytest.mark.parametrize("host", [1.0, 2.25])
def test_scalar_sphere_field_matches_analytic_solution(epsilon, host):
    """A scalar inclusion has the analytic field 3 epsilon_e/(epsilon_i+2 epsilon_e)."""
    actual = compute_internal_field_tensor(np.eye(3) / 3.0, epsilon, host)
    expected = np.eye(3) * (3.0 * host / (epsilon + 2.0 * host))
    np.testing.assert_allclose(actual, expected, rtol=1.0e-14, atol=1.0e-14)
