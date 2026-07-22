#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""Tensor-level Raman polar plot calculations."""

import numpy as np

RAMAN_ACTIVE_THRESHOLD = 1.0e-6
ACOUSTIC_THRESHOLD_CM1 = 10.0
POLAR_GEOMETRIES = ("VV", "HV", "Theta")


def normalise_vector(vector, name="vector"):
    """Return a normalised 3-vector."""
    vector = np.asarray(vector, dtype=float)
    if vector.shape != (3,):
        raise ValueError(f"{name} must be a 3-vector")
    norm = np.linalg.norm(vector)
    if norm <= 0.0:
        raise ValueError(f"{name} must be non-zero")
    return vector / norm


def build_polarisation_basis(rotation_axis, zero_angle_polarisation):
    """Build the polar plot basis from the light direction and zero-angle polarisation.

    Parameters
    ----------
    rotation_axis : array_like, shape (3,)
        Light propagation direction.
    zero_angle_polarisation : array_like, shape (3,)
        Incoming polarisation direction at zero polar angle.

    Returns
    -------
    tuple
        ``(rz, pza, pp, orthogonality)`` where all vectors are normalised and
        ``orthogonality = dot(rz, pza_input_normalised)``.

    """
    rz = normalise_vector(rotation_axis, "rotation_axis")
    pza_input = normalise_vector(zero_angle_polarisation, "zero_angle_polarisation")
    orthogonality = float(np.dot(rz, pza_input))

    # Project Pza onto the plane perpendicular to Rz so the plotted basis is
    # physically meaningful even if the input vector is slightly non-orthogonal.
    pza = pza_input - orthogonality * rz
    pza_norm = np.linalg.norm(pza)
    if pza_norm <= 1.0e-12:
        raise ValueError("zero_angle_polarisation must not be parallel to rotation_axis")
    pza = pza / pza_norm
    pp = np.cross(rz, pza)
    pp = pp / np.linalg.norm(pp)
    return rz, pza, pp, orthogonality


def rotate_around_axis(vector, axis, angle_rad):
    """Rotate a vector around an axis using Rodrigues' rotation formula."""
    vector = np.asarray(vector, dtype=float)
    axis = normalise_vector(axis, "axis")
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return vector * c + np.cross(axis, vector) * s + axis * np.dot(axis, vector) * (1.0 - c)


def polarisation_vectors(rotation_axis, zero_angle_polarisation, angle_rad, geometry):
    """Return incident and scattered polarisation vectors for one polar angle."""
    if geometry not in POLAR_GEOMETRIES:
        raise ValueError(f"Unknown polar geometry {geometry}")
    rz, pza, pp, _orthogonality = build_polarisation_basis(rotation_axis, zero_angle_polarisation)
    v = rotate_around_axis(pza, rz, angle_rad)
    h = rotate_around_axis(pp, rz, angle_rad)
    if geometry == "VV":
        return v, v
    if geometry == "HV":
        return v, h
    return pza, v


def intensity_for_modes(raman_tensors, mode_indices, rotation_axis, zero_angle_polarisation, angles_rad, geometry):
    """Calculate a polar intensity curve for selected Raman modes.

    Contributions from multiple modes are summed incoherently.
    """
    if geometry not in POLAR_GEOMETRIES:
        raise ValueError(f"Unknown polar geometry {geometry}")
    tensors = [np.asarray(raman_tensors[index], dtype=complex) for index in mode_indices]
    intensities = np.zeros(len(angles_rad), dtype=float)
    for i, angle in enumerate(angles_rad):
        e_l, e_s = polarisation_vectors(rotation_axis, zero_angle_polarisation, angle, geometry)
        value = 0.0
        for tensor in tensors:
            amplitude = e_s @ tensor @ e_l
            value += float(np.real(amplitude * np.conj(amplitude)))
        intensities[i] = value
    return intensities


def calculate_polar_curves(
        raman_tensors,
        mode_indices,
        rotation_axis=(0.0, 0.0, 1.0),
        zero_angle_polarisation=(1.0, 0.0, 0.0),
        angles_rad=None,
        geometries=POLAR_GEOMETRIES,
        normalise=True):
    """Calculate Raman polar curves for selected modes and geometries.

    Returns
    -------
    dict
        Keys are ``angles_rad``, ``raw_curves``, ``curves``,
        ``raw_max_intensity`` and ``orthogonality``.

    """
    if angles_rad is None:
        angles_rad = np.linspace(0.0, 2.0 * np.pi, 361)
    angles_rad = np.asarray(angles_rad, dtype=float)
    mode_indices = list(mode_indices)
    _rz, _pza, _pp, orthogonality = build_polarisation_basis(rotation_axis, zero_angle_polarisation)

    raw_curves = {}
    for geometry in geometries:
        raw_curves[geometry] = intensity_for_modes(
            raman_tensors, mode_indices, rotation_axis, zero_angle_polarisation, angles_rad, geometry)

    raw_max = max((float(np.max(values)) for values in raw_curves.values()), default=0.0)
    if normalise and raw_max > 0.0:
        curves = {geometry: values / raw_max for geometry, values in raw_curves.items()}
    else:
        curves = {geometry: values.copy() for geometry, values in raw_curves.items()}

    return {
        "angles_rad": angles_rad,
        "raw_curves": raw_curves,
        "curves": curves,
        "raw_max_intensity": raw_max,
        "orthogonality": orthogonality,
    }


def raman_active_mode_indices(frequencies_cm1, raman_activities, raman_tensors, threshold=RAMAN_ACTIVE_THRESHOLD):
    """Return zero-based indices for Raman-active modes, preserving table order."""
    active = []
    n_tensors = len(raman_tensors) if raman_tensors is not None else 0
    for index, (frequency, activity) in enumerate(zip(frequencies_cm1, raman_activities)):
        if index >= n_tensors:
            continue
        if frequency > ACOUSTIC_THRESHOLD_CM1 and activity > threshold:
            active.append(index)
    return active


def compute_lo_frequencies(q_hat_crystal, hessian, U_TO, Z_mat, eps_inf, volume_au):
    """Return the LO frequency (cm⁻¹) for each TO mode at phonon wavevector q̂.

    Diagonalises ``D_NAC(q̂) = D_TO + ΔD(q̂)`` and maps each NAC eigenmode back
    to the TO mode with the highest overlap ``C[n,m] = ⟨u_n^TO | u_m^NAC⟩``.

    Parameters
    ----------
    q_hat_crystal : ndarray, shape (3,)
        Unit phonon wavevector in the crystal frame (must be normalised).
    hessian : ndarray, shape (3N, 3N)
        Mass-weighted TO dynamical (Hessian) matrix.
    U_TO : ndarray, shape (n_to_modes, 3N)
        TO eigenvectors as rows (e.g. ``build_eigvecs_from_normal_modes(...).T``).
    Z_mat : ndarray, shape (3, 3N)
        Mass-weighted Born charge matrix from :func:`build_Z_mat`.
    eps_inf : ndarray, shape (3, 3)
        High-frequency optical dielectric tensor in the crystal frame.
    volume_au : float
        Unit-cell volume in Bohr³.

    Returns
    -------
    lo_freqs : ndarray, shape (n_to_modes,)
        LO frequency in cm⁻¹ for each TO mode (matched by dominant overlap
        with the NAC eigenmodes).  Modes unaffected by NAC (zero Born-charge
        projection onto q̂) return a frequency very close to their TO value.

    """
    import math

    from PDielec.Constants import wavenumber as _wn

    # NAC screening factor and correction matrix
    eps_b_q = float(q_hat_crystal @ eps_inf @ q_hat_crystal)
    if abs(eps_b_q) < 1e-12:
        eps_b_q = 1.0
    S_nac = np.outer(q_hat_crystal, q_hat_crystal) / eps_b_q
    delta_D = (4.0 * np.pi / volume_au) * (Z_mat.T @ S_nac @ Z_mat)

    # Diagonalise D_NAC
    eig_val, eig_vec = np.linalg.eigh(np.array(hessian, dtype=float) + delta_D)
    nac_freqs = np.array([
        (math.sqrt(abs(ev)) / _wn) * (1.0 if ev >= 0.0 else -1.0)
        for ev in eig_val
    ])

    # Overlap C[n_to, m_nac] = <u_n^TO | u_m^NAC>; pick best-matching NAC mode
    C = U_TO @ eig_vec   # (n_to_modes, 3N)
    n_to_modes = U_TO.shape[0]
    lo_freqs = np.zeros(n_to_modes)
    for n in range(n_to_modes):
        m_best = int(np.argmax(np.abs(C[n, :])))
        lo_freqs[n] = nac_freqs[m_best]

    return lo_freqs


def compute_lo_modes(q_hat_crystal, hessian, U_TO, Z_mat, eps_inf, volume_au):
    """Like :func:`compute_lo_frequencies` but also returns matched LO eigenvectors.

    Parameters
    ----------
    q_hat_crystal, hessian, U_TO, Z_mat, eps_inf, volume_au
        Same as :func:`compute_lo_frequencies`.

    Returns
    -------
    lo_freqs : ndarray, shape (n_to_modes,)
        LO frequency in cm⁻¹ for each TO mode.
    lo_eigvecs : ndarray, shape (3N, n_to_modes)
        Columns are the mass-weighted LO eigenvectors best-matched to each TO mode.

    """
    import math

    from PDielec.Constants import wavenumber as _wn

    eps_b_q = float(q_hat_crystal @ eps_inf @ q_hat_crystal)
    if abs(eps_b_q) < 1e-12:
        eps_b_q = 1.0
    S_nac = np.outer(q_hat_crystal, q_hat_crystal) / eps_b_q
    delta_D = (4.0 * np.pi / volume_au) * (Z_mat.T @ S_nac @ Z_mat)

    eig_val, eig_vec = np.linalg.eigh(np.array(hessian, dtype=float) + delta_D)
    nac_freqs = np.array([
        (math.sqrt(abs(ev)) / _wn) * (1.0 if ev >= 0.0 else -1.0)
        for ev in eig_val
    ])

    C = U_TO @ eig_vec   # (n_to_modes, 3N)
    n_to_modes = U_TO.shape[0]
    lo_freqs = np.zeros(n_to_modes)
    lo_eigvecs = np.zeros((eig_vec.shape[0], n_to_modes))
    for n in range(n_to_modes):
        m_best = int(np.argmax(np.abs(C[n, :])))
        lo_freqs[n] = nac_freqs[m_best]
        lo_eigvecs[:, n] = eig_vec[:, m_best]

    return lo_freqs, lo_eigvecs


def build_Z_mat(born_charges, masses_au):
    """Build the mass-weighted Born charge matrix.

    Parameters
    ----------
    born_charges : ndarray, shape (N, 3, 3)
        Born effective charge tensors.
    masses_au : ndarray, shape (N,)
        Atomic masses in atomic units (a.u.).

    Returns
    -------
    Z_mat : ndarray, shape (3, 3N)
        Mass-weighted Born charges: ``Z_mat[α, κβ] = Z*[κ, α, β] / √M_κ``.

    """
    import math
    nAtoms = len(masses_au)
    n_modes = 3 * nAtoms
    Z_mat = np.zeros((3, n_modes))
    for kappa in range(nAtoms):
        inv_sqrtM = 1.0 / math.sqrt(masses_au[kappa])
        for beta in range(3):
            Z_mat[:, kappa * 3 + beta] = born_charges[kappa, :, beta] * inv_sqrtM
    return Z_mat


def build_eigvecs_from_normal_modes(mass_weighted_normal_modes):
    """Convert ``reader.mass_weighted_normal_modes`` to an eigenvector matrix.

    Parameters
    ----------
    mass_weighted_normal_modes : list of array_like, each shape (N, 3)
        Mass-weighted TO normal modes as stored by the output readers.

    Returns
    -------
    eigvecs : ndarray, shape (3N, n_to_modes)
        Columns are the TO eigenvectors (convention consistent with
        ``numpy.linalg.eigh`` output), ready for use in
        :func:`apply_eo_correction`.

    """
    n_to_modes = len(mass_weighted_normal_modes)
    if n_to_modes == 0:
        return np.zeros((0, 0))
    n_atoms = len(mass_weighted_normal_modes[0])
    n_dof = 3 * n_atoms
    U_TO = np.zeros((n_to_modes, n_dof))
    for imode, mode in enumerate(mass_weighted_normal_modes):
        col = 0
        for atom in mode:
            U_TO[imode, col:col + 3] = atom
            col += 3
    return U_TO.T   # (3N, n_to_modes) — columns are TO eigenvectors


def apply_eo_correction(tensors, chi2_repsilon, q_hat_crystal, Z_mat, eigvecs, eps_inf):
    """Apply the electro-optic (EO) q-dependent correction to Raman tensors.

    Implements the electro-optic term in ``eq-nonanalytic`` and ``eq-nac_ramantensor``.
    For each mode *p* the correction is::

        ΔR_p = −2 · f(q̂) · [Z^mw(q̂) · u_p] / (q̂ᵀ ε_∞ q̂)

    where ``f_ij(q̂) = Σ_l χ^(2)_ijl q̂_l`` and
    ``Z^mw(q̂)[n] = (Z_mat^T q̂)[n]``.

    Parameters
    ----------
    tensors : list of ndarray, each (3, 3)
        Input physical bulk Raman tensors ``R_epsilon`` (not modified; corrected
        copies are returned).
    chi2_repsilon : ndarray, shape (3, 3, 3)
        Second-order NLO susceptibility χ^(2), already converted by the
        readers to the same Angstrom-based internal convention as
        ``R_epsilon`` Raman tensors.
    q_hat_crystal : ndarray, shape (3,)
        Unit phonon wavevector in the crystal frame.
    Z_mat : ndarray, shape (3, 3N)
        Mass-weighted Born charge matrix from :func:`build_Z_mat`.
    eigvecs : ndarray, shape (3N, n_modes)
        Eigenvector matrix — columns are the mode eigenvectors.  For TO modes
        pass ``U_TO.T`` (i.e. the output of
        :func:`build_eigvecs_from_normal_modes`); for NAC modes pass the
        ``eig_vec`` returned by ``numpy.linalg.eigh``.
    eps_inf : ndarray, shape (3, 3)
        High-frequency optical dielectric tensor in the crystal frame.

    Returns
    -------
    list of ndarray, each (3, 3)
        EO-corrected copies of *tensors*.

    Notes
    -----
    Reader Raman tensors and reader χ^(2) tensors use the same Angstrom-based
    ``R_epsilon`` convention used by
    :func:`PDielec.Calculator.raman_intensities`, so no additional unit
    conversion is applied here.

    """
    chi2_repsilon = np.asarray(chi2_repsilon, dtype=float)

    # f_ij[i,j] = Σ_l χ^(2)_ijl * q̂_l   (3×3 EO tensor contracted with q̂)
    f_ij = np.einsum("ijl,l->ij", chi2_repsilon, q_hat_crystal)

    # Dielectric screening factor
    eps_b_q = float(q_hat_crystal @ eps_inf @ q_hat_crystal)
    if abs(eps_b_q) < 1e-12:
        eps_b_q = 1.0

    # Z_q[n] = (Z_mat^T q̂)[n]  — project mass-weighted Born charges onto q̂
    Z_q = Z_mat.T @ q_hat_crystal   # shape (3N,)

    corrected = []
    n_modes = len(tensors)
    for p_idx in range(n_modes):
        scalar_p = float(np.dot(Z_q, eigvecs[:, p_idx]))
        corrected.append(tensors[p_idx] + (-2.0 * f_ij * scalar_p / eps_b_q))
    return corrected


def apply_particle_eo_correction(tensors, chi2_repsilon, K_particle, Z_mat, eigvecs):
    """Apply the finite-particle macroscopic EO correction to Raman tensors.

    Implements ``eq-particle-eo-raman-correction`` using the same electrostatic
    kernel as the particle dynamical matrix.  For particle mode *m*::

        z_m = Z_mat @ eigvecs[:, m]
        Delta R_m = -2 chi2_repsilon : (K_particle @ z_m)

    ``chi2_repsilon`` is already converted by the output readers to the
    internal ``R_epsilon`` convention, so no further volume or unit factor is
    applied here.  Unlike :func:`apply_eo_correction`, this function has no
    bulk propagation direction.

    Parameters
    ----------
    tensors : sequence of array_like, each shape (3, 3)
        Zero-field Raman tensors in the particle-mode basis.
    chi2_repsilon : array_like, shape (3, 3, 3)
        Second-order susceptibility in the reader ``R_epsilon`` convention.
    K_particle : array_like, shape (3, 3)
        Finite-particle kernel ``N_bg @ L / epsilon_e_b``.
    Z_mat : array_like, shape (3, n_dof)
        Mass-weighted Born-charge matrix.
    eigvecs : array_like, shape (n_dof, n_modes)
        Particle eigenvectors as columns.

    Returns
    -------
    list of ndarray
        EO-corrected tensor copies, one per particle mode.

    """
    chi2_repsilon = np.asarray(chi2_repsilon)
    K_particle = np.asarray(K_particle)
    Z_mat = np.asarray(Z_mat)
    eigvecs = np.asarray(eigvecs)

    if chi2_repsilon.shape != (3, 3, 3):
        raise ValueError(
            f"chi2_repsilon must have shape (3, 3, 3), got {chi2_repsilon.shape}")
    if K_particle.shape != (3, 3):
        raise ValueError(f"K_particle must have shape (3, 3), got {K_particle.shape}")
    if Z_mat.ndim != 2 or Z_mat.shape[0] != 3:
        raise ValueError(f"Z_mat must have shape (3, n_dof), got {Z_mat.shape}")
    if eigvecs.ndim != 2 or eigvecs.shape[0] != Z_mat.shape[1]:
        raise ValueError(
            "eigvecs must have shape (n_dof, n_modes) consistent with Z_mat; "
            f"got Z_mat {Z_mat.shape} and eigvecs {eigvecs.shape}")
    if len(tensors) != eigvecs.shape[1]:
        raise ValueError(
            "the number of tensors must equal the number of eigenvectors; "
            f"got {len(tensors)} tensors and {eigvecs.shape[1]} eigenvectors")

    corrected = []
    for mode, tensor in enumerate(tensors):
        tensor = np.asarray(tensor, dtype=complex)
        if tensor.shape != (3, 3):
            raise ValueError(f"tensor {mode} must have shape (3, 3), got {tensor.shape}")
        mode_charge = Z_mat @ eigvecs[:, mode]
        field_vector = K_particle @ mode_charge
        delta_raman = -2.0 * np.einsum("ijl,l->ij", chi2_repsilon, field_vector)
        corrected.append(tensor + delta_raman)
    return corrected
