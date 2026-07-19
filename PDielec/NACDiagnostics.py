#!/usr/bin/env python
#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec.
#
"""Direction-dependent non-analytic phonon diagnostics and branch tracking."""

import math

import numpy as np
from scipy.optimize import linear_sum_assignment

from PDielec import Calculator
from PDielec.Constants import amu, wavenumber
from PDielec.RamanPolarCalculator import apply_eo_correction


def group_degenerate_modes(frequencies_cm1, tolerance_cm1=0.05):
    """Return contiguous mode-index groups degenerate within ``tolerance_cm1``."""
    frequencies = np.asarray(frequencies_cm1, dtype=float)
    if frequencies.size == 0:
        return []
    groups = [[0]]
    for index in range(1, len(frequencies)):
        if abs(frequencies[index] - frequencies[index - 1]) <= tolerance_cm1:
            groups[-1].append(index)
        else:
            groups.append([index])
    return groups


def _normalised_q(q_hat_crystal):
    """Return a finite unit Cartesian direction."""
    q_hat = np.asarray(q_hat_crystal, dtype=float)
    norm = np.linalg.norm(q_hat)
    if norm <= 1.0e-14:
        raise ValueError("The NAC direction must be non-zero")
    return q_hat / norm


def compute_nac_direction_diagnostics(
    q_hat_crystal,
    hessian,
    born_charges,
    eps_inf,
    volume_au,
    masses_au,
    u_to,
    raman_tensors,
    to_sigmas,
    chi2_repsilon=None,
    degeneracy_tolerance_cm1=0.05,
):
    """Calculate a complete NAC eigensystem and its TO/Raman diagnostics.

    The overlap convention is ``C[n, m] = <u_n(TO)|u_m(NAC)>``. Raman tensors
    without the electro-optic term are obtained by transforming the complete TO
    tensor basis with ``C``. When ``chi2_repsilon`` is supplied, a second tensor
    set and activity table containing the EO correction are returned.

    Parameters are equivalent to
    :func:`PDielec.GUI.CrystalScenarioTab._compute_nac_dynamical_matrix_standalone`.

    Returns
    -------
    dict
        NumPy-backed diagnostic record containing frequencies, eigenvectors,
        full overlaps, effective charges, tensors, activities, linewidths, and
        degenerate-subspace groups.

    """
    q_hat = _normalised_q(q_hat_crystal)
    hessian = np.asarray(hessian, dtype=float)
    born_charges = np.asarray(born_charges, dtype=float)
    eps_inf = np.asarray(eps_inf, dtype=float)
    masses_au = np.asarray(masses_au, dtype=float)
    u_to = np.asarray(u_to, dtype=float)
    to_sigmas = np.asarray(to_sigmas, dtype=float)
    n_atoms = len(masses_au)
    n_modes = 3 * n_atoms

    if hessian.shape != (n_modes, n_modes):
        raise ValueError(f"Expected a {(n_modes, n_modes)} Hessian, got {hessian.shape}")
    if u_to.shape != (n_modes, n_modes):
        raise ValueError(f"Expected a {(n_modes, n_modes)} TO eigenvector matrix, got {u_to.shape}")

    # Z'[alpha, kappa beta] = Z*[kappa, alpha, beta] / sqrt(M_kappa).
    z_matrix = np.zeros((3, n_modes))
    for atom in range(n_atoms):
        start = 3 * atom
        z_matrix[:, start:start + 3] = born_charges[atom] / math.sqrt(masses_au[atom])

    epsilon_q = float(q_hat @ eps_inf @ q_hat)
    if abs(epsilon_q) <= 1.0e-12:
        raise ValueError("The optical dielectric tensor gives zero screening along the NAC direction")
    screening = np.outer(q_hat, q_hat) / epsilon_q
    delta_dynamical_matrix = (4.0 * np.pi / volume_au) * (z_matrix.T @ screening @ z_matrix)
    eigenvalues, eigenvectors_columns = np.linalg.eigh(hessian + delta_dynamical_matrix)
    eigenvectors = eigenvectors_columns.T
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) / wavenumber

    overlaps = u_to @ eigenvectors_columns
    squared_overlaps = np.abs(overlaps) ** 2
    overlap_sums = np.sum(squared_overlaps, axis=0)
    dominant_to = np.argmax(squared_overlaps, axis=0)

    tensors_without_eo = []
    linewidths = np.zeros(n_modes)
    tensors = [np.asarray(tensor, dtype=float) for tensor in raman_tensors]
    for nac_index in range(n_modes):
        tensor = np.zeros((3, 3), dtype=float)
        for to_index, weight in enumerate(overlaps[:, nac_index]):
            if to_index < len(tensors):
                tensor += weight * tensors[to_index]
        tensors_without_eo.append(tensor)
        to_index = dominant_to[nac_index]
        linewidths[nac_index] = to_sigmas[to_index] if to_index < len(to_sigmas) else 5.0

    if chi2_repsilon is None:
        tensors_with_eo = [tensor.copy() for tensor in tensors_without_eo]
    else:
        tensors_with_eo = apply_eo_correction(
            tensors_without_eo,
            chi2_repsilon,
            q_hat,
            z_matrix,
            eigenvectors_columns,
            eps_inf,
        )

    # Multiply by sqrt(amu_au) to report conventional e/sqrt(amu) units.
    effective_charge_vectors = (z_matrix @ eigenvectors_columns) * math.sqrt(amu)
    return {
        "q_hat_crystal": q_hat,
        "frequencies_cm1": frequencies,
        "eigenvectors": eigenvectors,
        "overlap_matrix": overlaps,
        "squared_overlap_matrix": squared_overlaps,
        "squared_overlap_sums": overlap_sums,
        "dominant_to_by_nac": dominant_to,
        "mode_effective_charge_vectors_e_sqrt_amu": effective_charge_vectors.T,
        "mode_effective_charge_magnitudes_e_sqrt_amu": np.linalg.norm(effective_charge_vectors, axis=0),
        "raman_tensors_without_eo": tensors_without_eo,
        "raman_tensors_with_eo": tensors_with_eo,
        "raman_activities_without_eo": Calculator.raman_intensities(tensors_without_eo, 1.0),
        "raman_activities_with_eo": Calculator.raman_intensities(tensors_with_eo, 1.0),
        "linewidths_cm1": linewidths,
        "degenerate_groups": group_degenerate_modes(frequencies, degeneracy_tolerance_cm1),
    }


def track_nac_branches(direction_diagnostics, degeneracy_tolerance_cm1=0.05):
    """Track NAC branches by mode overlaps and degenerate-subspace projectors."""
    if not direction_diagnostics:
        return []
    n_modes = len(direction_diagnostics[0]["frequencies_cm1"])
    tracked = [{
        "branch_ids_by_sorted_mode": np.arange(n_modes, dtype=int),
        "predecessor_sorted_mode": np.full(n_modes, -1, dtype=int),
        "adjacent_squared_overlap": np.ones(n_modes),
        "subspace_matches": [],
    }]

    for previous, current in zip(direction_diagnostics[:-1], direction_diagnostics[1:]):
        previous_tracking = tracked[-1]
        previous_vectors = np.asarray(previous["eigenvectors"])
        current_vectors = np.asarray(current["eigenvectors"])
        continuity = np.abs(previous_vectors @ current_vectors.T) ** 2
        previous_indices, current_indices = linear_sum_assignment(-continuity)
        predecessor = np.empty(n_modes, dtype=int)
        branch_ids = np.empty(n_modes, dtype=int)
        assigned_overlap = np.empty(n_modes)
        for previous_index, current_index in zip(previous_indices, current_indices):
            predecessor[current_index] = previous_index
            branch_ids[current_index] = previous_tracking["branch_ids_by_sorted_mode"][previous_index]
            assigned_overlap[current_index] = continuity[previous_index, current_index]

        previous_groups = group_degenerate_modes(
            previous["frequencies_cm1"], degeneracy_tolerance_cm1
        )
        current_groups = group_degenerate_modes(current["frequencies_cm1"], degeneracy_tolerance_cm1)
        subspace_matches = []
        for current_group in current_groups:
            scores = []
            current_subspace = current_vectors[current_group]
            for previous_group in previous_groups:
                previous_subspace = previous_vectors[previous_group]
                raw_score = np.linalg.norm(previous_subspace @ current_subspace.T, ord="fro") ** 2
                scores.append(raw_score / min(len(previous_group), len(current_group)))
            best_group = int(np.argmax(scores))
            subspace_matches.append({
                "current_modes": list(current_group),
                "previous_modes": list(previous_groups[best_group]),
                "projector_overlap": float(scores[best_group]),
            })

        tracked.append({
            "branch_ids_by_sorted_mode": branch_ids,
            "predecessor_sorted_mode": predecessor,
            "adjacent_squared_overlap": assigned_overlap,
            "subspace_matches": subspace_matches,
        })
    return tracked
