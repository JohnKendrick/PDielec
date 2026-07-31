#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Phase-4 tests for the modal-pair final-state replacement path."""

from multiprocessing.pool import ThreadPool

import numpy as np
import numpy.testing as npt

from PDielec.LayeredRamanCalculator import (
    DEPTH_INTEGRATION_INCOHERENT,
    MODAL_PAIR_COHERENT_ALL,
    MODAL_PAIR_GROUP_Q,
    MODAL_PAIR_INCOHERENT,
    LayeredRamanCalculator,
    RamanLayer,
    bose_factor,
)
from PDielec.Tests.Crystal_Raman.conftest import build_system

NU_MODE = 500.0
LASER_CM1 = 20000.0
THICKNESS_M = 1.0e-6


def _calculator(
    *,
    detected_pol="p",
    collection_side="superstrate",
    modal_pair_combination=MODAL_PAIR_GROUP_Q,
    modal_pair_use_nac=None,
    depth_integration=None,
    n_gauss=1,
):
    system = build_system([(THICKNESS_M, 2.0)], n_sup=1.0, n_sub=2.0)
    tensor = np.eye(3, dtype=complex)

    def nac_function(_q_hat_lab):
        return np.array([NU_MODE]), [tensor], np.array([5.0])

    layer = RamanLayer(
        layer_index=0,
        phonon_frequencies_cm1=np.array([NU_MODE]),
        raman_tensors=[tensor],
        rotation_matrix=np.eye(3),
        nac_function=nac_function,
    )
    kwargs = {}
    if depth_integration is not None:
        kwargs["depth_integration"] = depth_integration
    return LayeredRamanCalculator(
        system=system,
        raman_layers=[layer],
        laser_frequency_cm1=LASER_CM1,
        incident_angle_rad=0.0,
        collection_side=collection_side,
        incident_pol="p",
        detected_pol=detected_pol,
        temperature_K=0.0,
        linewidths_cm1=np.array([5.0]),
        n_gauss=n_gauss,
        approximate_es=True,
        modal_pairs=True,
        modal_pair_combination=modal_pair_combination,
        modal_pair_use_nac=modal_pair_use_nac,
        **kwargs,
    )


def _install_fake_modal_fields(calc, modal_fields, qs_by_layer, zeta=0.0):
    def fake_modal_fields(_freq_cm1, _system, _angle_rad, z_arr, system_layer_indices=None):
        del system_layer_indices
        assert len(z_arr) == modal_fields.shape[-1]
        return modal_fields.copy(), {0: np.array(qs_by_layer, dtype=complex)}, zeta

    calc._get_modal_fields_at_gl_points = fake_modal_fields


def test_phase4_group_q_sums_same_final_state_coherently():
    """GROUP_Q coherently combines optical-channel pairs with the same q direction."""
    calc = _calculator(detected_pol="s", modal_pair_combination=MODAL_PAIR_GROUP_Q, n_gauss=1)
    modal_fields = np.zeros((4, 2, 3, 1), dtype=complex)
    modal_fields[0, 0, 0, :] = 1.0
    modal_fields[1, 0, 0, :] = 1.0
    modal_fields[2, 1, 0, :] = 1.0
    _install_fake_modal_fields(calc, modal_fields, [2.0, 1.0, 0.0, -10.0])

    _freqs, intensities, _sigmas = calc.calculate_mode_intensities()

    weight = calc._gl_phys_weights[0]
    expected = 4.0 * weight**2 * bose_factor(NU_MODE, 0.0)
    npt.assert_allclose(intensities[0], expected, rtol=1.0e-12, atol=1.0e-12)


def test_phase4_incoherent_nonpolar_mode_groups_coherently():
    """Non-polar modal-pair modes are grouped coherently for INCOHERENT.

    When ``modal_pair_use_nac`` is False (non-polar mode), all Berreman pairs
    are q-independent and must be summed coherently regardless of the combination
    option.  The MODAL_PAIR_INCOHERENT path therefore uses the same coherent
    grouping as MODAL_PAIR_COHERENT_ALL, giving the same result as the total-field
    integral.  This is the corrected behaviour after the non-polar bypass fix.
    """
    calc = _calculator(
        modal_pair_combination=MODAL_PAIR_INCOHERENT,
        modal_pair_use_nac=[False],
        n_gauss=1,
    )
    modal_fields = np.zeros((4, 2, 3, 1), dtype=complex)
    modal_fields[0, 0, 0, :] = 1.0
    modal_fields[2, 0, 0, :] = 1.0
    _install_fake_modal_fields(calc, modal_fields, [1.0, 0.0, -1.0, 0.0])

    _freqs, intensities, _sigmas = calc.calculate_mode_intensities()

    weight = calc._gl_phys_weights[0]
    # Coherent sum of both pairs: amplitude = weight*(1+1)*(1+1) → squared = 16*weight²
    expected_coherent_grouping = 16.0 * weight**2 * bose_factor(NU_MODE, 0.0)
    npt.assert_allclose(intensities[0], expected_coherent_grouping, rtol=1.0e-12, atol=1.0e-12)


def test_phase4_coherent_all_nonpolar_mode_matches_total_field_integral():
    """COHERENT_ALL non-polar channel sum equals the total-field Raman integral."""
    calc = _calculator(
        modal_pair_combination=MODAL_PAIR_COHERENT_ALL,
        modal_pair_use_nac=[False],
        n_gauss=1,
    )
    modal_fields = np.zeros((4, 2, 3, 1), dtype=complex)
    modal_fields[0, 0, 0, :] = 1.0
    modal_fields[2, 0, 0, :] = 2.0
    _install_fake_modal_fields(calc, modal_fields, [1.0, 0.0, -1.0, 0.0])

    _freqs, intensities, _sigmas = calc.calculate_mode_intensities()

    weight = calc._gl_phys_weights[0]
    total_field_amplitude = weight * (1.0 + 2.0) * (1.0 + 2.0)
    expected = abs(total_field_amplitude) ** 2 * bose_factor(NU_MODE, 0.0)
    npt.assert_allclose(intensities[0], expected, rtol=1.0e-12, atol=1.0e-12)


def test_phase4_incoherent_depth_squares_local_final_state_integrand():
    """Incoherent depth integration uses the new local-integrand accumulation path."""
    calc = _calculator(
        detected_pol="s",
        modal_pair_combination=MODAL_PAIR_GROUP_Q,
        depth_integration=DEPTH_INTEGRATION_INCOHERENT,
        n_gauss=3,
    )
    modal_fields = np.zeros((4, 2, 3, 3), dtype=complex)
    modal_fields[0, 0, 0, :] = np.array([1.0, 2.0, 3.0])
    modal_fields[1, 0, 0, :] = np.array([4.0, 5.0, 6.0])
    modal_fields[2, 1, 0, :] = 1.0
    _install_fake_modal_fields(calc, modal_fields, [2.0, 1.0, 0.0, -10.0])

    _freqs, intensities, _sigmas = calc.calculate_mode_intensities()
    with ThreadPool(2) as pool:
        _parallel_freqs, parallel_intensities, _parallel_sigmas = calc.calculate_mode_intensities(pool=pool)

    expected_local_integrand = np.array([5.0, 7.0, 9.0])
    expected = np.dot(calc._gl_phys_weights, expected_local_integrand**2) * bose_factor(NU_MODE, 0.0)
    npt.assert_allclose(intensities[0], expected, rtol=1.0e-12, atol=1.0e-12)
    npt.assert_allclose(parallel_intensities[0], expected, rtol=1.0e-12, atol=1.0e-12)


def test_phase4_incoherent_forward_cross_pairs_share_external_final_state():
    """Internal p-s q differences remain local-field terms in forward scattering."""
    calc = _calculator(
        detected_pol="s",
        collection_side="substrate",
        modal_pair_combination=MODAL_PAIR_GROUP_Q,
        depth_integration=DEPTH_INTEGRATION_INCOHERENT,
        n_gauss=3,
    )
    modal_fields = np.zeros((4, 2, 3, 3), dtype=complex)
    modal_fields[0, 0, 0, :] = np.array([1.0, 2.0, 3.0])
    modal_fields[1, 0, 0, :] = np.array([4.0, 5.0, 6.0])
    modal_fields[2, 1, 0, :] = 1.0
    _install_fake_modal_fields(calc, modal_fields, [2.0, 1.0, 0.5, -10.0])

    _freqs, intensities, _sigmas = calc.calculate_mode_intensities()
    with ThreadPool(2) as pool:
        _parallel_freqs, parallel_intensities, _parallel_sigmas = calc.calculate_mode_intensities(pool=pool)

    expected_local_integrand = np.array([5.0, 7.0, 9.0])
    expected = np.dot(calc._gl_phys_weights, expected_local_integrand**2) * bose_factor(NU_MODE, 0.0)
    npt.assert_allclose(intensities[0], expected, rtol=1.0e-12, atol=1.0e-12)
    npt.assert_allclose(parallel_intensities[0], expected, rtol=1.0e-12, atol=1.0e-12)
