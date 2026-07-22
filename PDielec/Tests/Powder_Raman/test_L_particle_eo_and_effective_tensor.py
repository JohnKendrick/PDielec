"""Finite-particle EO and two-frequency powder Raman tensor tests."""

import inspect

import numpy as np

from PDielec.Calculator import (
    compute_effective_raman_tensor,
    compute_internal_field_tensor,
    compute_powder_raman_intensities,
    sobol_rotations,
)
from PDielec.RamanPolarCalculator import apply_eo_correction, apply_particle_eo_correction


def _rotation():
    """Return a fixed proper orthogonal matrix."""
    axis = np.array([1.0, 2.0, -0.5])
    axis /= np.linalg.norm(axis)
    angle = 0.71
    cross = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) + np.sin(angle) * cross + (1.0 - np.cos(angle)) * (cross @ cross)


def test_effective_tensor_has_one_input_and_one_output_factor():
    """The effective tensor is N_S.T R N_L, with no third local-field factor."""
    raman = np.array([[1.0, 0.2, -0.1], [0.4, 2.0, 0.3], [0.0, -0.2, 0.7]])
    n_laser = np.array([[0.8, 0.1, 0.0], [0.0, 0.7, 0.2], [0.1, 0.0, 0.9]])
    n_scattered = np.array([[0.6, -0.1, 0.2], [0.0, 0.75, 0.0], [0.05, 0.1, 0.85]])

    result = compute_effective_raman_tensor(raman, n_laser, n_scattered)

    np.testing.assert_allclose(result, n_scattered.T @ raman @ n_laser)


def test_effective_tensor_uses_transpose_not_hermitian_conjugate():
    """Reciprocal source screening requires an ordinary transpose for complex N."""
    raman = np.diag([1.0, 2.0, 3.0]).astype(complex)
    n_complex = np.array([
        [0.8 + 0.1j, 0.2 - 0.05j, 0.0],
        [0.0, 0.7 + 0.2j, 0.1j],
        [0.05, 0.0, 0.9 - 0.1j],
    ])

    result = compute_effective_raman_tensor(raman, n_complex, n_complex)

    np.testing.assert_allclose(result, n_complex.T @ raman @ n_complex)
    assert not np.allclose(result, n_complex.conj().T @ raman @ n_complex)


def test_effective_tensor_is_rotation_covariant_with_distinct_optical_factors():
    """Co-rotating R, N_L, and N_S co-rotates the two-frequency effective tensor."""
    rotation = _rotation()
    raman = np.arange(9.0).reshape(3, 3) / 7.0
    n_laser = np.array([[0.8, 0.1, 0.0], [0.0, 0.7, 0.2], [0.1, 0.0, 0.9]])
    n_scattered = np.array([[0.6, -0.1, 0.2], [0.0, 0.75, 0.0], [0.05, 0.1, 0.85]])
    reference = compute_effective_raman_tensor(raman, n_laser, n_scattered)

    candidate = compute_effective_raman_tensor(
        rotation @ raman @ rotation.T,
        rotation @ n_laser @ rotation.T,
        rotation @ n_scattered @ rotation.T,
    )

    np.testing.assert_allclose(candidate, rotation @ reference @ rotation.T, atol=1.0e-12)


def test_particle_eo_matches_direct_contraction_and_does_not_modify_input():
    """Particle EO uses -2 chi2 contracted with K_particle Z_mw u_m."""
    tensors = [np.diag([1.0, 2.0, 3.0]), np.eye(3) * 0.5]
    original = [tensor.copy() for tensor in tensors]
    chi2 = np.arange(27.0).reshape(3, 3, 3) / 10.0
    kernel = np.array([[0.4, 0.1, 0.0], [0.1, 0.3, 0.0], [0.0, 0.0, 0.2]])
    z_mat = np.array([[1.0, 0.0], [0.2, 0.5], [0.0, -0.4]])
    eigvecs = np.array([[1.0, 0.0], [0.0, 1.0]])

    corrected = apply_particle_eo_correction(tensors, chi2, kernel, z_mat, eigvecs)

    for mode in range(2):
        field_vector = kernel @ (z_mat @ eigvecs[:, mode])
        expected = original[mode] - 2.0 * np.einsum("ijl,l->ij", chi2, field_vector)
        np.testing.assert_allclose(corrected[mode], expected)
        np.testing.assert_array_equal(tensors[mode], original[mode])


def test_particle_eo_zero_limits_and_linear_scaling():
    """The correction vanishes with zero chi2 or charges and scales bilinearly."""
    tensors = [np.eye(3)]
    chi2 = np.ones((3, 3, 3))
    kernel = np.diag([0.4, 0.3, 0.2])
    z_mat = np.array([[1.0], [0.5], [-0.2]])
    eigvecs = np.ones((1, 1))

    zero_chi = apply_particle_eo_correction(tensors, np.zeros_like(chi2), kernel, z_mat, eigvecs)
    zero_z = apply_particle_eo_correction(tensors, chi2, kernel, np.zeros_like(z_mat), eigvecs)
    once = apply_particle_eo_correction(tensors, chi2, kernel, z_mat, eigvecs)[0] - tensors[0]
    twice_chi = apply_particle_eo_correction(tensors, 2.0 * chi2, kernel, z_mat, eigvecs)[0] - tensors[0]
    twice_z = apply_particle_eo_correction(tensors, chi2, kernel, 2.0 * z_mat, eigvecs)[0] - tensors[0]

    np.testing.assert_allclose(zero_chi[0], tensors[0])
    np.testing.assert_allclose(zero_z[0], tensors[0])
    np.testing.assert_allclose(twice_chi, 2.0 * once)
    np.testing.assert_allclose(twice_z, 2.0 * once)


def test_particle_eo_is_rotation_covariant():
    """Co-rotating chi2, K, charges, and Raman tensors rotates the result."""
    rotation = _rotation()
    tensor = np.array([[1.0, 0.2, 0.0], [0.2, 0.5, 0.1], [0.0, 0.1, 0.8]])
    chi2 = np.arange(27.0).reshape(3, 3, 3) / 30.0
    kernel = np.array([[0.4, 0.05, 0.0], [0.05, 0.25, 0.02], [0.0, 0.02, 0.3]])
    z_mat = np.array([[1.0], [0.3], [-0.2]])
    eigvecs = np.ones((1, 1))
    reference = apply_particle_eo_correction([tensor], chi2, kernel, z_mat, eigvecs)[0]

    tensor_rot = rotation @ tensor @ rotation.T
    chi2_rot = np.einsum("ia,jb,kc,abc->ijk", rotation, rotation, rotation, chi2)
    kernel_rot = rotation @ kernel @ rotation.T
    z_mat_rot = rotation @ z_mat
    candidate = apply_particle_eo_correction(
        [tensor_rot], chi2_rot, kernel_rot, z_mat_rot, eigvecs)[0]

    np.testing.assert_allclose(candidate, rotation @ reference @ rotation.T, atol=1.0e-12)


def test_particle_eo_respects_complete_mode_phase_reversal():
    """Reversing both a mode and its zero-field tensor only reverses the final tensor."""
    tensor = np.array([[1.0, 0.2, 0.0], [0.2, 0.5, 0.1], [0.0, 0.1, 0.8]])
    chi2 = np.arange(27.0).reshape(3, 3, 3) / 30.0
    kernel = np.diag([0.4, 0.25, 0.3])
    z_mat = np.array([[1.0], [0.3], [-0.2]])
    eigvecs = np.ones((1, 1))

    positive = apply_particle_eo_correction([tensor], chi2, kernel, z_mat, eigvecs)[0]
    negative = apply_particle_eo_correction([-tensor], chi2, kernel, z_mat, -eigvecs)[0]

    np.testing.assert_allclose(negative, -positive)
    np.testing.assert_allclose(
        compute_powder_raman_intensities(negative),
        compute_powder_raman_intensities(positive),
    )


def test_particle_eo_api_has_no_propagation_direction():
    """The finite-particle EO API contains no bulk propagation-vector argument."""
    parameters = inspect.signature(apply_particle_eo_correction).parameters
    assert not {"q", "q_hat", "q_hat_crystal"}.intersection(parameters)


def test_eo_corrected_sphere_analytic_and_sobol_averages_agree():
    """An EO-corrected sphere tensor gives the same analytic and sampled VV/VH averages."""
    tensor = np.array([[1.0, 0.2, 0.0], [0.2, 0.5, 0.1], [0.0, 0.1, 0.8]])
    chi2 = np.arange(27.0).reshape(3, 3, 3) / 50.0
    corrected = apply_particle_eo_correction(
        [tensor], chi2, np.eye(3) / 9.0, np.array([[1.0], [0.3], [-0.2]]), np.ones((1, 1))
    )[0]
    vv_analytic, vh_analytic = compute_powder_raman_intensities(corrected)
    e_laser = np.array([0.0, 1.0, 0.0])
    e_vv = e_laser
    e_vh = np.array([1.0, 0.0, 0.0])
    vv_sampled = 0.0
    vh_sampled = 0.0
    rotations = sobol_rotations(8192, seed=42)
    for rotation in rotations:
        rotated = rotation @ corrected @ rotation.T
        vv_sampled += 45.0 * abs(e_vv @ rotated @ e_laser) ** 2
        vh_sampled += 45.0 * abs(e_vh @ rotated @ e_laser) ** 2

    np.testing.assert_allclose(
        [vv_sampled / len(rotations), vh_sampled / len(rotations)],
        [vv_analytic, vh_analytic],
        rtol=2.0e-3,
    )


def test_particle_plate_kernel_matches_bulk_eo_kernel():
    """The particle helper reproduces bulk EO when K_particle equals the plate kernel."""
    q_hat = np.array([0.0, 0.0, 1.0])
    eps_inf = np.diag([2.0, 3.0, 5.0])
    kernel = np.outer(q_hat, q_hat) / (q_hat @ eps_inf @ q_hat)
    tensors = [np.diag([1.0, 0.5, 2.0])]
    chi2 = np.arange(27.0).reshape(3, 3, 3) / 50.0
    z_mat = np.array([[0.2], [-0.1], [0.7]])
    eigvecs = np.ones((1, 1))

    particle = apply_particle_eo_correction(tensors, chi2, kernel, z_mat, eigvecs)
    bulk = apply_eo_correction(tensors, chi2, q_hat, z_mat, eigvecs, eps_inf)

    np.testing.assert_allclose(particle, bulk)


def test_isotropic_sphere_particle_kernel_has_surface_optical_denominator():
    """A scalar sphere gives K_particle = I/(epsilon_i_b + 2 epsilon_e_b)."""
    epsilon_i = 5.0
    epsilon_e = 2.0
    depolarisation = np.eye(3) / 3.0
    n_background = compute_internal_field_tensor(
        depolarisation, epsilon_i * np.eye(3), epsilon_e)

    kernel = n_background @ depolarisation / epsilon_e

    np.testing.assert_allclose(
        kernel, np.eye(3) / (epsilon_i + 2.0 * epsilon_e), atol=1.0e-12)
