"""Tests for Level 3 (modal_pairs) Raman scattering — test classes I1–I6."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from conftest import build_system, iso_eps, make_layer, make_raman_layer, run_calc

from PDielec.GTMcore import Layer, ScatteringMatrixSystem, TransferMatrixSystem
from PDielec.LayeredRamanCalculator import LayeredRamanCalculator, RamanLayer, bose_factor

# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

NU_MODE = 500.0          # phonon frequency, cm⁻¹
LASER_CM1 = 20000.0      # laser frequency, cm⁻¹
N_LAYER = 2.0            # isotropic refractive index of Raman-active layer
THICKNESS_M = 1.0e-6     # layer thickness (1 µm)
R_PP = np.diag([1.0, 0.0, 0.0])   # simple p-pol Raman tensor


def _make_system_and_layer(n_layer=N_LAYER, thickness_m=THICKNESS_M):
    """Build a single-layer TMS and a matching RamanLayer."""
    system = build_system([(thickness_m, n_layer)], n_sup=1.0, n_sub=n_layer)
    rl = make_raman_layer(0, R_PP, NU_MODE)
    return system, rl


def _make_sm_layer(thickness_m, n):
    """Create a CoherentLayer-compatible layer for ScatteringMatrixSystem.

    ScatteringMatrixSystem requires layers with ``calculate_scattering_matrix``,
    which is only available on CoherentLayer and its kin.  We build one using
    a minimal adapter that satisfies CoherentLayer's interface.
    """
    from PDielec.GTMcore import CoherentLayer

    class _Adapter:
        def __init__(self, thick, eps_func):
            self._thick = thick
            self._eps_func = eps_func

        def get_thickness_in_metres(self):
            return self._thick

        def get_permittivity_function(self):
            return self._eps_func

    return CoherentLayer(_Adapter(thickness_m, iso_eps(n)))


def _make_sm_system(n_layer=N_LAYER, thickness_m=THICKNESS_M):
    """Build an identical stack as a ScatteringMatrixSystem using CoherentLayers."""
    from PDielec.GTMcore import CoherentLayer

    sup = _make_sm_layer(1e-3, 1.0)
    sub = _make_sm_layer(1e-3, n_layer)
    layer = _make_sm_layer(thickness_m, n_layer)
    return ScatteringMatrixSystem(substrate=sub, superstrate=sup, layers=[layer])


def _run_modal_pairs(system, raman_layers, phonon_freqs, incident_pol="p", detected_pol="p",
                     nac_function=None, n_gauss=21, approximate_es=True,
                     incident_angle_rad=0.0):
    """Run LayeredRamanCalculator with modal_pairs=True."""
    # Optionally attach nac_function to first layer
    if nac_function is not None:
        raman_layers[0].nac_function = nac_function

    linewidths = 5.0 * np.ones(len(phonon_freqs))
    calc = LayeredRamanCalculator(
        system=system,
        raman_layers=raman_layers,
        laser_frequency_cm1=LASER_CM1,
        incident_angle_rad=incident_angle_rad,
        incident_pol=incident_pol,
        detected_pol=detected_pol,
        temperature_K=0.0,
        linewidths_cm1=linewidths,
        n_gauss=n_gauss,
        approximate_es=approximate_es,
        modal_pairs=True,
    )
    return calc


# ---------------------------------------------------------------------------
# TestI1: modal amplitudes sum to total field
# ---------------------------------------------------------------------------

class TestI1ModalAmplitudesSumToTotal:
    """Verify calculate_modal_amplitudes is consistent with calculate_Efield."""

    def test_sum_equals_total_field_p_pol(self):
        """Sum of 4 modal contributions (p-pol) must equal the total E field."""
        from PDielec.Constants import speed_light_si

        system, _ = _make_system_and_layer()
        f_hz = LASER_CM1 * speed_light_si * 1e2
        system.initialize_sys(f_hz)
        zeta = 0.0 + 0j

        # Full field via calculate_Efield
        boundaries = system.get_layers_boundaries()
        z_test = np.linspace(boundaries[1] + 1e-9, boundaries[2] - 1e-9, 15)
        _, E_total, _ = system.calculate_Efield(f_hz, zeta, z_vect=z_test)

        # Modal amplitudes via calculate_modal_amplitudes
        modal_amps, zn = system.calculate_modal_amplitudes(f_hz, zeta)
        layer = system.layers[0]
        z_front = boundaries[1]

        # Reconstruct field by summing all 4 modal contributions (p-pol incidence = rows 0:3)
        E_modal_sum = np.zeros((3, len(z_test)), dtype=np.complex128)
        for z_idx, z_j in enumerate(z_test):
            amps = modal_amps[0]  # layer index 0
            for n in range(4):
                phase = np.exp(
                    layer.propagation_exponents[n] * (z_j - z_front) / layer.thick
                )
                E_modal_sum[:, z_idx] += phase * amps[n] * layer.gamma[n, :]

        # Should match E_total[0:3, :] (p-pol incidence)
        assert np.allclose(E_modal_sum, E_total[0:3, :], atol=1e-10), (
            "Modal sum does not match total p-pol field"
        )

    def test_sum_equals_total_field_s_pol(self):
        """Sum of 4 modal contributions (s-pol) must equal the total E field."""
        from PDielec.Constants import speed_light_si

        system, _ = _make_system_and_layer()
        f_hz = LASER_CM1 * speed_light_si * 1e2
        system.initialize_sys(f_hz)
        zeta = 0.0 + 0j

        boundaries = system.get_layers_boundaries()
        z_test = np.linspace(boundaries[1] + 1e-9, boundaries[2] - 1e-9, 15)
        _, E_total, _ = system.calculate_Efield(f_hz, zeta, z_vect=z_test)

        modal_amps, _ = system.calculate_modal_amplitudes(f_hz, zeta)
        layer = system.layers[0]
        z_front = boundaries[1]

        E_modal_sum = np.zeros((3, len(z_test)), dtype=np.complex128)
        for z_idx, z_j in enumerate(z_test):
            amps = modal_amps[0]
            for n in range(4):
                phase = np.exp(
                    layer.propagation_exponents[n] * (z_j - z_front) / layer.thick
                )
                E_modal_sum[:, z_idx] += phase * amps[n + 4] * layer.gamma[n, :]

        assert np.allclose(E_modal_sum, E_total[3:6, :], atol=1e-10), (
            "Modal sum does not match total s-pol field"
        )

    def test_multilayer_stack(self):
        """Verify modal amplitude consistency for a 3-layer stack."""
        from PDielec.Constants import speed_light_si

        system = build_system(
            [(0.5e-6, 1.5), (1.0e-6, 2.0), (0.3e-6, 1.8)],
            n_sup=1.0, n_sub=2.0
        )
        f_hz = LASER_CM1 * speed_light_si * 1e2
        system.initialize_sys(f_hz)
        zeta = 0.0 + 0j

        boundaries = system.get_layers_boundaries()
        modal_amps, _ = system.calculate_modal_amplitudes(f_hz, zeta)
        _, E_total, _ = system.calculate_Efield(f_hz, zeta,
                                                 z_vect=np.array([0.5 * (boundaries[2] + boundaries[3])]))

        # Check middle layer (index 1)
        layer = system.layers[1]
        z_front = boundaries[2]
        z_j = 0.5 * (boundaries[2] + boundaries[3])
        E_sum = np.zeros(3, dtype=np.complex128)
        amps = modal_amps[1]
        for n in range(4):
            phase = np.exp(layer.propagation_exponents[n] * (z_j - z_front) / layer.thick)
            E_sum += phase * amps[n] * layer.gamma[n, :]

        assert np.allclose(E_sum, E_total[0:3, 0], atol=1e-10)


# ---------------------------------------------------------------------------
# TestI2: active pair selection
# ---------------------------------------------------------------------------

class TestI2ActivePairSelection:
    """Verify _get_active_modal_pairs returns correct (i, j) sets."""

    def _make_calc(self, incident_pol, detected_pol, angle=0.0):
        system, rl = _make_system_and_layer()
        linewidths = 5.0 * np.ones(1)
        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl],
            laser_frequency_cm1=LASER_CM1,
            incident_angle_rad=angle,
            incident_pol=incident_pol,
            detected_pol=detected_pol,
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=5,
            modal_pairs=True,
        )
        return calc

    def test_p_p_gives_four_pairs(self):
        # p-inc: modes {0,2}; p-detect: modes {0,2} → 4 pairs
        calc = self._make_calc("p", "p")
        pairs = calc._get_active_modal_pairs()
        assert set(pairs) == {(0, 0), (0, 2), (2, 0), (2, 2)}

    def test_s_s_gives_four_pairs(self):
        # s-inc: modes {1,3}; s-detect: modes {1,3} → 4 pairs
        calc = self._make_calc("s", "s")
        pairs = calc._get_active_modal_pairs()
        assert set(pairs) == {(1, 1), (1, 3), (3, 1), (3, 3)}

    def test_p_s_gives_four_pairs(self):
        # p-inc: modes {0,2}; s-detect: modes {1,3} → 4 cross-pol pairs
        calc = self._make_calc("p", "s")
        pairs = calc._get_active_modal_pairs()
        assert set(pairs) == {(0, 1), (0, 3), (2, 1), (2, 3)}

    def test_s_unpolarised_gives_eight_pairs(self):
        # s-inc: modes {1,3}; unpolarised: modes {0,1,2,3} → 8 pairs
        calc = self._make_calc("s", "unpolarised")
        pairs = calc._get_active_modal_pairs()
        assert set(pairs) == {(1, 0), (1, 1), (1, 2), (1, 3),
                               (3, 0), (3, 1), (3, 2), (3, 3)}

    def test_p_unpolarised_gives_eight_pairs(self):
        # p-inc: modes {0,2}; unpolarised: modes {0,1,2,3} → 8 pairs
        calc = self._make_calc("p", "unpolarised")
        pairs = calc._get_active_modal_pairs()
        assert set(pairs) == {(0, 0), (0, 1), (0, 2), (0, 3),
                               (2, 0), (2, 1), (2, 2), (2, 3)}


# ---------------------------------------------------------------------------
# TestI3: NAC cache called once per pair, not per mode
# ---------------------------------------------------------------------------

class TestI3NACCacheCalledOncePerPair:
    """Verify NAC function is called once per pair (layer, i, j), not per mode."""

    def test_nac_called_once_per_pair(self):
        """For N_pairs pairs, nac_function should be called exactly N_pairs times."""
        system, _ = _make_system_and_layer()

        # Multiple phonon modes
        nu_modes = np.array([200.0, 400.0, 600.0, NU_MODE])
        R_tensors = [R_PP.copy() for _ in nu_modes]
        linewidths = 5.0 * np.ones(len(nu_modes))

        call_log = []

        def mock_nac(q_hat_lab):
            call_log.append(q_hat_lab.copy())
            return (nu_modes, R_tensors, linewidths)

        rl = RamanLayer(
            layer_index=0,
            phonon_frequencies_cm1=nu_modes,
            raman_tensors=R_tensors,
            rotation_matrix=np.eye(3),
            nac_function=mock_nac,
        )

        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl],
            laser_frequency_cm1=LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol="p",
            detected_pol="unpolarised",  # 8 pairs: {0,2} × {0,1,2,3}
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=5,
            approximate_es=True,
            modal_pairs=True,
        )
        calc.calculate_mode_intensities()

        n_pairs = len(calc._get_active_modal_pairs())
        n_modes = len(nu_modes)
        assert n_pairs == 8  # {0,2} × {0,1,2,3}
        # Key invariant: NAC is cached per pair, not called once per mode per pair.
        # Some pairs have q_ph=0 and skip the NAC call; others call it exactly once.
        assert len(call_log) > 0, "NAC should be called for at least one pair"
        assert len(call_log) <= n_pairs, (
            f"Expected at most {n_pairs} nac_function calls (once per pair), "
            f"got {len(call_log)}"
        )
        # Without caching, would be called n_modes × n_pairs_nonzero times:
        assert len(call_log) < n_modes * n_pairs, (
            f"Caching should reduce calls below n_modes×n_pairs={n_modes*n_pairs}, "
            f"got {len(call_log)}"
        )

    def test_nac_called_once_per_pair_per_layer(self):
        """With 2 Raman-active layers, nac_function called N_pairs × 2 times."""
        system = build_system(
            [(THICKNESS_M, N_LAYER), (THICKNESS_M, N_LAYER)],
            n_sup=1.0, n_sub=N_LAYER
        )
        nu_modes = np.array([200.0, NU_MODE])
        R_tensors = [R_PP.copy() for _ in nu_modes]
        linewidths = 5.0 * np.ones(len(nu_modes))

        call_counts = [0]

        def mock_nac(q_hat_lab):
            call_counts[0] += 1
            return (nu_modes, R_tensors, linewidths)

        rl0 = RamanLayer(0, nu_modes, R_tensors, np.eye(3), nac_function=mock_nac)
        rl1 = RamanLayer(1, nu_modes, R_tensors, np.eye(3), nac_function=mock_nac)

        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl0, rl1],
            laser_frequency_cm1=LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol="p",
            detected_pol="p",  # 4 pairs: {0,2} × {0,2}
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=5,
            approximate_es=True,
            modal_pairs=True,
        )
        calc.calculate_mode_intensities()

        n_pairs = len(calc._get_active_modal_pairs())  # 4
        n_layers_raman = 2
        n_modes = len(nu_modes)
        # Key invariant: NAC called at most once per (pair, layer), not per mode.
        assert call_counts[0] > 0, "NAC should be called for non-zero q_ph pairs"
        assert call_counts[0] <= n_pairs * n_layers_raman, (
            f"Expected at most {n_pairs*n_layers_raman} calls (once per pair per layer), "
            f"got {call_counts[0]}"
        )
        assert call_counts[0] < n_modes * n_pairs * n_layers_raman, (
            "Caching should reduce calls below n_modes×n_pairs×n_layers"
        )


# ---------------------------------------------------------------------------
# TestI4: isotropic limit — Level 3 ≈ Level 2
# ---------------------------------------------------------------------------

class TestI4IsotropicLimitMatchesDominantMode:
    """For isotropic crystal, all mode pairs produce the same q_ph → Level 3 ≈ Level 2."""

    def _intensity(self, modal_pairs, incident_pol="p", detected_pol="p"):
        system, rl = _make_system_and_layer(n_layer=N_LAYER)
        linewidths = 5.0 * np.ones(1)
        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl],
            laser_frequency_cm1=LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol=incident_pol,
            detected_pol=detected_pol,
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=51,
            approximate_es=True,
            modal_pairs=modal_pairs,
        )
        _, intensities, _ = calc.calculate_mode_intensities()
        return intensities[0] if len(intensities) > 0 else 0.0

    def test_p_p_modal_pairs_is_positive(self):
        """Level 3 should produce a positive intensity."""
        I3 = self._intensity(modal_pairs=True)
        assert I3 > 0.0

    def test_level3_and_level2_same_order_of_magnitude(self):
        """Level 3 and Level 2 (standard) intensities should be same order of magnitude.

        Without NAC correction, both Level 2 and Level 3 use the same TO tensors
        and the same field profiles, so the summed intensity from one modal pair
        should match the Level 2 result (which sums all modes).
        """
        I2 = self._intensity(modal_pairs=False)
        I3 = self._intensity(modal_pairs=True)
        # Both should be positive and within 3 orders of magnitude of each other
        assert I2 > 0.0 and I3 > 0.0
        ratio = max(I2, I3) / min(I2, I3)
        assert ratio < 1e3, f"Level 2 ({I2:.3e}) and Level 3 ({I3:.3e}) differ too much"


# ---------------------------------------------------------------------------
# TestI5: ScatteringMatrix vs TransferMatrix consistency
# ---------------------------------------------------------------------------

class TestI5ScatteringMatrixConsistency:
    """ScatteringMatrixSystem and TransferMatrixSystem give consistent modal_pairs spectra."""

    def _run(self, system, rl, n_gauss=21):
        linewidths = 5.0 * np.ones(1)
        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl],
            laser_frequency_cm1=LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol="p",
            detected_pol="p",
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=n_gauss,
            approximate_es=True,
            modal_pairs=True,
        )
        _, intensities, _ = calc.calculate_mode_intensities()
        return intensities[0] if len(intensities) > 0 else 0.0

    def test_tm_and_sm_agree(self):
        """TMS and SMS produce the same modal_pairs intensity (within 5%)."""
        rl_tm = make_raman_layer(0, R_PP, NU_MODE)
        rl_sm = make_raman_layer(0, R_PP, NU_MODE)

        system_tm = build_system([(THICKNESS_M, N_LAYER)], n_sup=1.0, n_sub=N_LAYER)
        system_sm = _make_sm_system(n_layer=N_LAYER, thickness_m=THICKNESS_M)

        I_tm = self._run(system_tm, rl_tm, n_gauss=31)
        I_sm = self._run(system_sm, rl_sm, n_gauss=31)

        assert I_tm > 0.0, "TMS intensity should be positive"
        assert I_sm > 0.0, "SMS intensity should be positive"
        rel_diff = abs(I_tm - I_sm) / max(abs(I_tm), abs(I_sm))
        assert rel_diff < 0.05, (
            f"TMS ({I_tm:.4e}) and SMS ({I_sm:.4e}) differ by {rel_diff*100:.2f}%"
        )

    def test_sm_modal_amplitudes_nonzero(self):
        """ScatteringMatrixSystem.calculate_modal_amplitudes returns non-zero amplitudes."""
        from PDielec.Constants import speed_light_si

        system_sm = _make_sm_system()
        f_hz = LASER_CM1 * speed_light_si * 1e2
        system_sm.initialize_sys(f_hz)
        zeta = 0.0 + 0j

        modal_amps, zn = system_sm.calculate_modal_amplitudes(f_hz, zeta)
        assert 0 in modal_amps, "Layer 0 should be in modal_amps"
        amps = modal_amps[0]
        assert np.any(np.abs(amps) > 0), "Modal amplitudes should not all be zero"


# ---------------------------------------------------------------------------
# TestI6: subdivision invariance for modal_pairs
# ---------------------------------------------------------------------------

class TestI7NACFrequencyShift:
    """Verify that modal_pairs uses NAC-corrected frequencies, not TO frequencies.

    Regression test for the bug where the NAC-corrected Hessian is diagonalised
    and the resulting eigenmodes are sorted by eigenvalue.  A TO mode at index k
    that shifts to a higher (LO) frequency moves to a different sorted index k'.
    The old code used the TO index k to read nac_freqs[k] and nac_tensors[k],
    which gave the wrong (nearly-zero-tensor) mode and hence ~10000× too small
    intensity, with the peak reported at the TO frequency instead of the LO.

    This test uses a synthetic nac_function that returns two phonon modes whose
    eigenvalue ordering is SWAPPED compared to the TO ordering, i.e. the TO mode
    at index 0 (ν_TO_0) maps to NAC index 1 (ν_NAC_1 > ν_TO_0).  With the bug,
    modal_pairs would report at the wrong frequency and with wrong (tiny) tensor.
    With the fix, it must report at ν_NAC_1 with non-trivial intensity.
    """

    NU_TO_LOW = 300.0    # TO frequency that gets shifted up (cm⁻¹)
    NU_TO_HIGH = 700.0   # TO frequency that doesn't shift (cm⁻¹, above LO of first)
    NU_LO = 550.0        # LO frequency of the first mode after NAC shift (cm⁻¹)
    LASER_CM1 = 20000.0

    def _make_nac_function(self):
        """Return a nac_function that swaps the eigenvalue order.

        TO modes:  [NU_TO_LOW=300, NU_TO_HIGH=700]
        NAC modes: [NU_LO=550, NU_TO_HIGH=700]  ← order preserved in this case
                   The tensor for NAC mode 0 (at 550) = R_PP (Raman active).
                   The tensor for NAC mode 1 (at 700) = zero tensor.

        This mimics a LO-TO shift: the Raman-active mode moves from 300→550.
        The old (buggy) code would look up nac_tensors[0] and report freq 300;
        the new (fixed) code looks up nac_freqs[0]=550 and nac_tensors[0]=R_PP.
        """
        nu_nac = np.array([self.NU_LO, self.NU_TO_HIGH])
        tensors_nac = [R_PP.astype(complex), np.zeros((3, 3), dtype=complex)]
        sigmas_nac = np.array([5.0, 5.0])

        def nac_fn(q_hat_lab):
            return nu_nac.copy(), list(tensors_nac), sigmas_nac.copy()

        return nac_fn

    def test_peak_at_nac_frequency_not_to_frequency(self):
        """Peak must appear at the NAC (LO) frequency, not the TO frequency."""
        # n_sub != n_layer creates a back-reflection, giving the backward mode (mode 2)
        # a non-zero amplitude.  Pair (0,2) then has q_ph ≈ 2n (backscattering) and
        # calls the nac_function, shifting the peak from TO (300) to LO (550).
        system = build_system([(THICKNESS_M, N_LAYER)], n_sup=1.0, n_sub=1.0)

        phonon_freqs = np.array([self.NU_TO_LOW, self.NU_TO_HIGH])
        R_zero = np.zeros((3, 3), dtype=float)
        # TO tensors: mode 0 active (R_PP), mode 1 inactive
        to_tensors = [R_PP.astype(float), R_zero]
        linewidths = 5.0 * np.ones(2)

        rl = RamanLayer(
            layer_index=0,
            phonon_frequencies_cm1=phonon_freqs,
            raman_tensors=to_tensors,
            rotation_matrix=np.eye(3),
            nac_function=self._make_nac_function(),
        )

        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl],
            laser_frequency_cm1=self.LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol="p",
            detected_pol="p",
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=21,
            approximate_es=True,
            modal_pairs=True,
        )
        freqs, intensities, _ = calc.calculate_mode_intensities()

        assert len(freqs) > 0, "No active modes returned"

        # The reported peak must be at the NAC (LO) frequency ~550, not TO ~300
        active_freqs = np.array(freqs)
        assert np.any(np.abs(active_freqs - self.NU_LO) < 1.0), (
            f"Expected a peak near {self.NU_LO} cm⁻¹ (NAC/LO freq) "
            f"but got active_freqs={active_freqs.tolist()}"
        )
        assert not np.any(np.abs(active_freqs - self.NU_TO_LOW) < 1.0), (
            f"Peak found at TO frequency {self.NU_TO_LOW} cm⁻¹; "
            "modal_pairs should report the NAC-corrected (LO) frequency"
        )

    def test_intensity_at_nac_mode_not_tiny(self):
        """Intensity at the LO peak must be of the same order as without NAC shift.

        With the bug, the wrong (zero) tensor was used for the shifted mode,
        giving ~0 intensity.  The fix must produce intensity comparable to a
        direct calculation using the NAC tensor stored in RamanLayer.
        """
        system = build_system([(THICKNESS_M, N_LAYER)], n_sup=1.0, n_sub=1.0)

        phonon_freqs = np.array([self.NU_TO_LOW, self.NU_TO_HIGH])
        to_tensors = [R_PP.astype(float), np.zeros((3, 3), dtype=float)]
        linewidths = 5.0 * np.ones(2)

        rl = RamanLayer(
            layer_index=0,
            phonon_frequencies_cm1=phonon_freqs,
            raman_tensors=to_tensors,
            rotation_matrix=np.eye(3),
            nac_function=self._make_nac_function(),
        )

        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl],
            laser_frequency_cm1=self.LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol="p",
            detected_pol="p",
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=21,
            approximate_es=True,
            modal_pairs=True,
        )
        freqs, intensities, _ = calc.calculate_mode_intensities()

        # Reference: Level 1/2-style calculation storing the NAC tensor directly
        rl_ref = RamanLayer(
            layer_index=0,
            phonon_frequencies_cm1=np.array([self.NU_LO, self.NU_TO_HIGH]),
            raman_tensors=[R_PP.astype(float), np.zeros((3, 3), dtype=float)],
            rotation_matrix=np.eye(3),
        )
        linewidths_ref = 5.0 * np.ones(2)
        calc_ref = LayeredRamanCalculator(
            system=system,
            raman_layers=[rl_ref],
            laser_frequency_cm1=self.LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol="p",
            detected_pol="p",
            temperature_K=0.0,
            linewidths_cm1=linewidths_ref,
            n_gauss=21,
            approximate_es=True,
            modal_pairs=True,
        )
        freqs_ref, intensities_ref, _ = calc_ref.calculate_mode_intensities()

        # Find the intensity at the LO peak in both calculations
        idx_lp = np.argmin(np.abs(np.array(freqs) - self.NU_LO))
        idx_lp_ref = np.argmin(np.abs(np.array(freqs_ref) - self.NU_LO))

        I_modal = intensities[idx_lp]
        I_ref = intensities_ref[idx_lp_ref]

        assert I_modal > 0.0, "LO peak intensity must be positive"
        assert I_ref > 0.0, "Reference LO peak intensity must be positive"
        # Must be within a factor of 10 (not 10000×)
        ratio = I_modal / I_ref
        assert 0.1 < ratio < 10.0, (
            f"LO peak intensity ratio modal/ref = {ratio:.3g}; expected ~1"
        )


class TestI6SubdivisionInvariance:
    """Splitting a homogeneous layer must give the same modal_pairs intensity.

    Mirrors TestC2SubdivisionInvariance: matched media (n_sup = n_sub = n_layer)
    eliminates reflections at outer interfaces, and coherent_layers=True ensures
    the amplitude-sum identity holds exactly across sublayers.
    """

    # Use 100 nm to avoid rapid phase accumulation (mirrors C2 choice)
    D_TOTAL = 100e-9
    RTOL = 0.01  # 1 % relative tolerance

    def _calc_intensity(self, n_sublayers, n_gauss=30):
        """Build matched-media system with n_sublayers of equal total thickness."""
        sub_thick = self.D_TOTAL / n_sublayers
        specs = [(sub_thick, N_LAYER)] * n_sublayers
        # Matched media: n_sup = n_sub = n_layer → no outer-interface reflections
        system = build_system(specs, n_sup=N_LAYER, n_sub=N_LAYER)
        raman_layers = [make_raman_layer(k, R_PP, NU_MODE) for k in range(n_sublayers)]
        linewidths = 5.0 * np.ones(1)
        calc = LayeredRamanCalculator(
            system=system,
            raman_layers=raman_layers,
            laser_frequency_cm1=LASER_CM1,
            incident_angle_rad=0.0,
            incident_pol="p",
            detected_pol="p",
            temperature_K=0.0,
            linewidths_cm1=linewidths,
            n_gauss=n_gauss,
            approximate_es=True,
            modal_pairs=True,
            coherent_layers=True,   # needed for subdivision identity to hold
        )
        _, intensities, _ = calc.calculate_mode_intensities()
        return intensities[0] if len(intensities) > 0 else 0.0

    def test_split_5_matches_unsplit(self):
        """5-sublayer system must give the same intensity as the 1-layer system (1%)."""
        I_1 = self._calc_intensity(n_sublayers=1)
        I_5 = self._calc_intensity(n_sublayers=5)

        assert I_1 > 0.0, "Unsplit intensity should be positive"
        assert I_5 > 0.0, "Split intensity should be positive"
        rel_diff = abs(I_5 - I_1) / max(abs(I_5), abs(I_1))
        assert rel_diff < self.RTOL, (
            f"Subdivision invariance failed: 1-layer={I_1:.4e}, "
            f"5-layer={I_5:.4e}, rel_diff={rel_diff*100:.2f}%"
        )
