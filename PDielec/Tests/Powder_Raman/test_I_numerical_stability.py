"""Tests I1-I5: Numerical stability, mode ordering, and scaling invariance.

Plan sections: I1, I2, I3, I4, I5, J10, J11, J12, M7, M8, M9, M10.

Tests verify:
- stability near local-field resonances (I1);
- convergence with number of orientation samples (I2);
- frequency grid convergence (I3);
- invariance to mode ordering (I4 / J12);
- tensor scaling laws (I5 / J10);
- zero Born charge gives no polar frequency shift (J10);
- surface-mode resonance position for scalar sphere (J11).
"""
import os
import sys
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest

from PDielec.Calculator import (
    compute_internal_field_tensor,
    compute_particle_raman_tensor,
    compute_powder_raman_intensities,
    initialise_sphere_depolarisation_matrix,
    raman_intensities,
    sobol_rotations,
)


def _powder_average_numerical(raman_tensor, n_samples=2048, seed=42):
    R_cryst = np.asarray(raman_tensor, dtype=complex)
    e_L = np.array([0.0, 1.0, 0.0])
    e_VV = np.array([0.0, 1.0, 0.0])
    e_VH = np.array([1.0, 0.0, 0.0])
    vv_sum = 0.0
    vh_sum = 0.0
    for G in sobol_rotations(n_samples, seed=seed):
        R_lab = G @ R_cryst @ G.T
        vv_sum += 45.0 * abs(e_VV @ R_lab @ e_L) ** 2
        vh_sum += 45.0 * abs(e_VH @ R_lab @ e_L) ** 2
    return vv_sum / n_samples, vh_sum / n_samples

L_sphere = initialise_sphere_depolarisation_matrix()


# ---------------------------------------------------------------------------
# I1 — Matrix inversion stability near resonance
# ---------------------------------------------------------------------------

def test_i1_no_nan_near_resonance():
    """I1: No NaN or Inf from compute_N when ε_i is near the surface-mode condition."""
    eps_e = 1.0
    # Near ε_i = -2 ε_e = -2, with a small imaginary part
    for delta in np.linspace(-0.5, 0.5, 21):
        eps_i = ((-2.0 + delta) + 0.05j) * np.eye(3, dtype=complex)
        N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
        assert np.all(np.isfinite(N)), f"N has non-finite elements at Δ={delta}"


def test_i1_N_continuous_across_resonance():
    """I1: N is continuous (no sudden jumps) as ε_i crosses the resonance."""
    eps_e = 1.0
    gamma = 0.1
    eps_real = np.linspace(-4.0, 0.0, 51)
    N_prev = None
    for eps_r in eps_real:
        eps_i = (eps_r + 1j * gamma) * np.eye(3, dtype=complex)
        N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
        if N_prev is not None:
            diff = np.max(np.abs(N - N_prev))
            assert diff < 25.0, f"N jumped discontinuously at ε_r={eps_r}: diff={diff}"
        N_prev = N.copy()


def test_i1_R_particle_finite_near_resonance():
    """I1: R_particle remains finite near the surface-mode resonance."""
    eps_e = 1.0
    R = np.diag([1.0, 1.0, 1.0])
    for eps_r in np.linspace(-4.0, 0.0, 21):
        eps_i = (eps_r + 0.05j) * np.eye(3, dtype=complex)
        N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
        R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i, eps_e)
        vv, vh = compute_powder_raman_intensities(R_p)
        assert np.isfinite(vv), f"VV not finite at ε_r={eps_r}"
        assert np.isfinite(vh), f"VH not finite at ε_r={eps_r}"


# ---------------------------------------------------------------------------
# I2 — Orientation sample convergence
# ---------------------------------------------------------------------------

def test_i2_VV_converges_with_more_samples():
    """I2: VV converges to a stable value as orientation sample count increases."""
    R = np.diag([2.0, 1.0, 0.5])
    vv_analytic, _ = compute_powder_raman_intensities(R)
    errors = []
    for n in [64, 256, 1024, 4096]:
        vv_num, _ = _powder_average_numerical(R, n_samples=n)
        errors.append(abs(vv_num - vv_analytic))
    # The error should generally decrease (allow one non-monotonic step)
    n_increases = sum(1 for i in range(len(errors)-1) if errors[i+1] > errors[i])
    assert n_increases <= 1, f"Error not converging: {errors}"


def test_i2_VH_converges_with_more_samples():
    """I2: VH converges to a stable value as orientation sample count increases."""
    R = np.diag([2.0, 1.0, 0.5])
    _, vh_analytic = compute_powder_raman_intensities(R)
    for n in [256, 2048]:
        _, vh_num = _powder_average_numerical(R, n_samples=n)
        rel_err = abs(vh_num - vh_analytic) / (vh_analytic + 1e-10)
        assert rel_err < 0.1, f"n={n}: VH relative error {rel_err:.4f}"


def test_i2_depolarisation_ratio_stable():
    """I2: Depolarisation ratio converges and is stable for large n."""
    R = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    for n in [1024, 4096]:
        vv_num, vh_num = _powder_average_numerical(R, n_samples=n)
        rho = vh_num / vv_num
        assert abs(rho - 0.75) < 0.1, f"n={n}: rho={rho:.4f}"


# ---------------------------------------------------------------------------
# I3 — Frequency grid convergence (Lorentzian integration)
# ---------------------------------------------------------------------------

def _lorentzian_spectrum(vs_cm1, freq, sigma, amplitude):
    """Lorentzian lineshape: amplitude * sigma / ((v - freq)^2 + sigma^2)."""
    return amplitude * sigma / ((vs_cm1 - freq) ** 2 + sigma ** 2)


def test_i3_peak_position_independent_of_grid_spacing():
    """I3: Lorentzian peak position does not depend on frequency grid spacing."""
    freq = 300.0
    sigma = 5.0
    amplitude = 100.0
    peak_positions = []
    for n_pts in [100, 200, 500, 1000]:
        vs = np.linspace(200.0, 400.0, n_pts)
        spec = _lorentzian_spectrum(vs, freq, sigma, amplitude)
        peak_positions.append(vs[np.argmax(spec)])
    for pos in peak_positions:
        assert abs(pos - freq) < 2.0 * (400.0 - 200.0) / 100, \
            f"Peak position {pos} too far from true {freq}"


def test_i3_peak_height_scales_inversely_with_sigma():
    """I3: Peak height of a Lorentzian scales as amplitude/sigma."""
    freq = 300.0
    amplitude = 100.0
    vs = np.linspace(200.0, 400.0, 1000)
    for sigma in [2.0, 5.0, 10.0]:
        spec = _lorentzian_spectrum(vs, freq, sigma, amplitude)
        peak = np.max(spec)
        expected = amplitude / sigma
        assert abs(peak - expected) / expected < 0.01, \
            f"sigma={sigma}: peak={peak:.2f}, expected={expected:.2f}"


# ---------------------------------------------------------------------------
# I4 / J12 — Invariance to mode ordering
# ---------------------------------------------------------------------------

def test_i4_mode_ordering_invariance_analytical():
    """I4/J12: Shuffling mode list does not change the summed spectrum (analytical path)."""
    # Simulate a spectrum as sum of Lorentzians for several modes
    raman_tensors = [
        np.diag([1.0, 1.0, 1.0]),
        np.diag([2.0, 1.0, 0.5]),
        np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
    ]
    frequencies = [100.0, 200.0, 300.0]
    sigmas = [5.0, 5.0, 5.0]

    vs = np.linspace(50.0, 400.0, 400)

    def compute_spectrum(tensors, freqs, sigs):
        spec = np.zeros(len(vs))
        for R, freq, sigma in zip(tensors, freqs, sigs):
            vv, _ = compute_powder_raman_intensities(R)
            spec += vv * sigma / ((vs - freq) ** 2 + sigma ** 2)
        return spec

    spec_orig = compute_spectrum(raman_tensors, frequencies, sigmas)
    # Shuffle
    order = [2, 0, 1]
    spec_shuffled = compute_spectrum(
        [raman_tensors[i] for i in order],
        [frequencies[i] for i in order],
        [sigmas[i] for i in order],
    )
    assert np.allclose(spec_orig, spec_shuffled, atol=1e-12), \
        "Spectrum changed after mode reordering"


def test_i4_mode_reorder_reverse():
    """I4/J12: Reversing the mode order gives the same spectrum."""
    raman_tensors = [
        np.diag([1.0, 1.0, 1.0]),
        np.diag([2.0, 1.0, 0.5]),
        np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
    ]
    frequencies = [100.0, 200.0, 300.0]
    sigmas = [5.0, 5.0, 5.0]
    vs = np.linspace(50.0, 400.0, 400)

    def compute_spectrum(tensors, freqs, sigs):
        spec = np.zeros(len(vs))
        for R, freq, sigma in zip(tensors, freqs, sigs):
            vv, _ = compute_powder_raman_intensities(R)
            spec += vv * sigma / ((vs - freq) ** 2 + sigma ** 2)
        return spec

    spec_fwd = compute_spectrum(raman_tensors, frequencies, sigmas)
    spec_rev = compute_spectrum(raman_tensors[::-1], frequencies[::-1], sigmas[::-1])
    assert np.allclose(spec_fwd, spec_rev, atol=1e-12)


# ---------------------------------------------------------------------------
# I5 / J10 — Tensor scaling laws
# ---------------------------------------------------------------------------

def test_i5_raman_tensor_scaling_by_2_gives_4x_intensity():
    """I5: Scaling R by 2 gives 4× VV and VH intensities."""
    R = np.diag([1.0, 2.0, 0.5])
    vv1, vh1 = compute_powder_raman_intensities(R)
    vv2, vh2 = compute_powder_raman_intensities(2.0 * R)
    assert abs(vv2 / vv1 - 4.0) < 1e-10
    assert abs(vh2 / vh1 - 4.0) < 1e-10


def test_i5_complex_phase_does_not_change_intensity():
    """I5: Multiplying R by a complex phase leaves intensities unchanged."""
    R = np.diag([1.0, 2.0, 0.5])
    vv_ref, vh_ref = compute_powder_raman_intensities(R)
    for phase_angle in [0.3, 1.0, 1.57, 2.5, 3.14]:
        R_phased = np.exp(1j * phase_angle) * R
        vv, vh = compute_powder_raman_intensities(R_phased)
        assert abs(vv - vv_ref) < 1e-10, f"phase={phase_angle}: VV changed"
        assert abs(vh - vh_ref) < 1e-10, f"phase={phase_angle}: VH changed"


def test_i5_raman_intensities_parallel_column():
    """I5: PDielec.Calculator.raman_intensities parallel column scales correctly."""
    R = np.diag([1.0, 1.0, 1.0])
    volume = 1.0
    act1 = raman_intensities([R], volume)
    act2 = raman_intensities([2.0 * R], volume)
    assert abs(act2[0, 1] / act1[0, 1] - 4.0) < 1e-10, \
        "Parallel activity should scale as 4× when R is doubled"


def test_i5_raman_intensities_perpendicular_column():
    """I5: PDielec.Calculator.raman_intensities perpendicular column scales correctly."""
    R = np.diag([2.0, 1.0, 0.5])
    volume = 1.0
    act1 = raman_intensities([R], volume)
    act2 = raman_intensities([3.0 * R], volume)
    assert abs(act2[0, 2] / act1[0, 2] - 9.0) < 1e-10, \
        "Perpendicular activity should scale as 9× when R is tripled"


def test_i5_raman_intensities_accepts_complex_tensors_without_cast_warning():
    """I5: Complex Raman tensors are handled with magnitude-squared invariants."""
    R = np.array([
        [1.0 + 0.2j, 0.3 - 0.1j, 0.0],
        [0.3 - 0.1j, 2.0 + 0.4j, 0.2j],
        [0.0, 0.2j, 0.5 - 0.3j],
    ])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        activities = raman_intensities([R], volume=2.0)

    assert not any("Casting complex values to real" in str(w.message) for w in caught)
    assert np.all(np.isfinite(activities))
    assert np.all(activities >= 0.0)


def test_i5_raman_intensities_complex_phase_does_not_change_activity():
    """I5: Global complex phase does not change Raman activities."""
    R = np.array([
        [1.0, 0.2, 0.1],
        [0.2, 2.0, -0.3],
        [0.1, -0.3, 0.5],
    ], dtype=complex)
    reference = raman_intensities([R], volume=1.5)
    phased = raman_intensities([np.exp(0.7j) * R], volume=1.5)
    np.testing.assert_allclose(phased, reference, rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# J11 — Surface-mode resonance position for isotropic sphere
# ---------------------------------------------------------------------------

def test_j11_surface_mode_resonance_near_minus_2():
    """J11: Raman enhancement peaks near ε_i = -2 ε_e for a sphere in a scalar matrix."""
    eps_e = 1.0
    gamma = 0.05
    R = np.eye(3)
    vv_ref, _ = compute_powder_raman_intensities(R)

    eps_sweep = np.linspace(-4.0, 0.0, 101)
    enhancements = []
    for eps_r in eps_sweep:
        eps_i = (eps_r + 1j * gamma) * np.eye(3, dtype=complex)
        N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
        R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i, eps_e)
        vv, _ = compute_powder_raman_intensities(R_p)
        enhancements.append(vv / (vv_ref + 1e-30))

    peak_idx = np.argmax(enhancements)
    eps_peak = eps_sweep[peak_idx]
    assert abs(eps_peak - (-2.0 * eps_e)) < 0.15, \
        f"Resonance peak at ε_r={eps_peak:.3f}, expected near {-2*eps_e}"


def test_j11_resonance_shifts_with_matrix_permittivity():
    """J11: Sphere resonance shifts to ε_i = -2 ε_e as ε_e changes."""
    R = np.eye(3)
    gamma = 0.05
    for eps_e in [1.0, 2.0, 4.0]:
        eps_sweep = np.linspace(-4.0 * eps_e - 1.0, 0.0, 201)
        enhancements = []
        for eps_r in eps_sweep:
            eps_i = (eps_r + 1j * gamma) * np.eye(3, dtype=complex)
            N = compute_internal_field_tensor(L_sphere, eps_i, eps_e)
            R_p = compute_particle_raman_tensor(R, N, L_sphere, eps_i, eps_e)
            vv, _ = compute_powder_raman_intensities(R_p)
            enhancements.append(vv)
        peak_idx = np.argmax(enhancements)
        eps_peak = eps_sweep[peak_idx]
        assert abs(eps_peak - (-2.0 * eps_e)) < 0.3, \
            f"ε_e={eps_e}: peak at {eps_peak:.3f}, expected near {-2*eps_e}"


# ---------------------------------------------------------------------------
# Additional regression: raman_intensities function from Calculator
# ---------------------------------------------------------------------------

def test_raman_intensities_isotropic_zero_perp():
    """Calculator.raman_intensities: Isotropic R = a*I gives zero perpendicular activity."""
    a = 2.0
    R = a * np.eye(3)
    volume = 10.0
    act = raman_intensities([R], volume)
    assert act[0, 1] == pytest.approx(act[0, 0])
    assert act[0, 2] == pytest.approx(0.0, abs=1e-12)


def test_raman_intensities_sum_parallel_perp():
    """Calculator.raman_intensities: total equals VV plus VH for a general tensor."""
    R = np.array([[2.0, 1.0, 0.0], [-0.25, 1.0, 0.0], [0.0, 0.0, 0.5]])
    volume = 1.0
    act = raman_intensities([R], volume)
    assert abs(act[0, 0] - act[0, 1] - act[0, 2]) < 1e-10, \
        "total != parallel + perpendicular"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
