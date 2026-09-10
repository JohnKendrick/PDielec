"""Independent physics regressions for the September 2026 Raman review."""

import numpy as np
import pytest

from PDielec.Constants import amu, wavenumber
from PDielec.NACDiagnostics import compute_nac_direction_diagnostics


def nac(tensors, born=None, chi2=None):
    """Evaluate a small nondegenerate system with an explicit TO basis."""
    return compute_nac_direction_diagnostics(
        np.array([1.0, 2.0, 3.0]) / np.sqrt(14),
        np.diag((np.array([100.0, 200.0, 300.0]) * wavenumber) ** 2),
        np.zeros((1, 3, 3)) if born is None else born,
        np.eye(3) * 2,
        100.0,
        np.array([amu]),
        np.eye(3),
        tensors,
        np.ones(3),
        chi2_repsilon=chi2,
    )


def test_complex_nac_linearity_and_phase():
    """Zero NAC preserves complex phases, and real mode mixing is complex-linear."""
    rng = np.random.default_rng(11)
    a, b = rng.normal(size=(2, 3, 3, 3))
    np.testing.assert_allclose(nac(a + 1j * b)["raman_tensors_without_eo"], a + 1j * b)
    np.testing.assert_allclose(nac(a)["raman_activities_without_eo"], nac(1j * a)["raman_activities_without_eo"])
    born = np.array([[[1.0, 0.2, 0.3], [0.4, 2.0, 0.1], [0.2, 0.5, 3.0]]])
    mixed = nac(a + 1j * b, born)
    expected = np.asarray(nac(a, born)["raman_tensors_without_eo"]) + 1j * np.asarray(
        nac(b, born)["raman_tensors_without_eo"]
    )
    np.testing.assert_allclose(mixed["raman_tensors_without_eo"], expected, atol=1e-12)
    chi2 = np.arange(27.0).reshape(3, 3, 3) / 100
    with_eo = nac(a + 1j * b, born, chi2)
    eo_only = nac(np.zeros_like(a), born, chi2)["raman_tensors_with_eo"]
    np.testing.assert_allclose(with_eo["raman_tensors_with_eo"], expected + eo_only, atol=1e-12)


class InlinePool:
    """Exercise worker dispatch deterministically without shared optical state races."""

    def imap(self, function, values, chunksize=1):
        """Map tasks using the same worker API as a process pool."""
        return map(function, values)


def slab(frequencies, laser=20000.0, approximate=True, modal=False):
    """Construct a vacuum slab whose optical overlap has no spectral dispersion."""
    from PDielec.LayeredRamanCalculator import LayeredRamanCalculator, RamanLayer
    from PDielec.Tests.Crystal_Raman.conftest import build_system

    return LayeredRamanCalculator(
        system=build_system([(1e-7, 1.0)]),
        raman_layers=[RamanLayer(0, frequencies, [np.eye(3)] * len(frequencies), np.eye(3))],
        laser_frequency_cm1=laser,
        incident_angle_rad=0.0,
        incident_pol="s",
        detected_pol="s",
        temperature_K=0.0,
        linewidths_cm1=[5.0] * len(frequencies),
        n_gauss=10,
        approximate_es=approximate,
        depth_integration="Incoherent intensity",
        modal_pairs=modal,
        modal_pair_use_nac=[False] * len(frequencies),
    )


@pytest.mark.parametrize("modal", [False, True])
@pytest.mark.parametrize("approximate", [False, True])
@pytest.mark.parametrize("worker", [False, True])
def test_invalid_stokes_modes_are_excluded(modal, approximate, worker):
    """Field approximations and workers cannot admit unstable or impossible Stokes lines."""
    calc = slab(
        [-100.0, 0.0, 9.99, np.nan, np.inf, 20000.0, 21000.0, 10.0, 100.0], approximate=approximate, modal=modal
    )
    freq, intensity, _ = calc.calculate_mode_intensities(pool=InlinePool() if worker else None)
    np.testing.assert_allclose(sorted(freq), [10.0, 100.0])
    assert np.all(np.isfinite(intensity)) and np.all(intensity >= 0)


@pytest.mark.parametrize("laser", [0.0, -1.0, np.nan, np.inf])
def test_invalid_laser_is_rejected(laser):
    """Laser frequency is validated before solving optics."""
    with pytest.raises(ValueError, match="laser|Laser"):
        slab([100.0], laser=laser)


@pytest.mark.parametrize("modal", [False, True])
@pytest.mark.parametrize("approximate", [False, True])
@pytest.mark.parametrize("worker", [False, True])
@pytest.mark.parametrize("laser", [10000.0, 20000.0])
def test_crystal_stokes_radiation_prefactor(modal, approximate, worker, laser):
    """An analytic slab checks absolute weighting, line ratios and contribution totals."""
    calc = slab([100.0, 1000.0], laser, approximate, modal)
    freq, intensity, _ = calc.calculate_mode_intensities(pool=InlinePool() if worker else None)
    # The depth integration uses Angstrom weights: this slab is 1000 Angstrom thick.
    # Vacuum unit-illumination fields give a raw incoherent-depth overlap of 1000.
    np.testing.assert_allclose(intensity, 1000 * (laser - freq)**4 / freq, rtol=1e-10)
    expected_ratio = 0.1 * ((laser - 1000.0) / (laser - 100.0)) ** 4
    assert intensity[1] / intensity[0] == pytest.approx(expected_ratio, rel=1e-10)
    other = slab([100.0, 1000.0], laser * 2, approximate, modal)
    _, intensity_other, _ = other.calculate_mode_intensities(pool=InlinePool() if worker else None)
    np.testing.assert_allclose(intensity_other / intensity, ((2 * laser - freq) / (laser - freq)) ** 4, rtol=1e-10)
    f, i, _, contributions = calc.calculate_mode_intensities(return_contributions=True)
    for frequency, total in zip(f, i):
        assert sum(c.intensity for c in contributions if c.frequency == frequency) == pytest.approx(total)


@pytest.mark.parametrize("worker", [False, True])
@pytest.mark.parametrize("corrected", [-100.0, 800.0, 21000.0])
def test_modal_nac_uses_final_frequency_for_guard_and_weight(worker, corrected):
    """Final NAC states, rather than their TO seed, determine allowed lines and weights."""
    calc = slab([100.0], modal=True)
    calc.modal_pair_use_nac = [True]
    calc.raman_layers[0].nac_function = lambda q: (np.array([corrected]), [np.eye(3)], np.array([5.0]))
    freq, intensity, _ = calc.calculate_mode_intensities(pool=InlinePool() if worker else None)
    if corrected <= 0 or corrected >= 20000:
        assert len(freq) == 0
    else:
        np.testing.assert_allclose(freq, [corrected])
        _, _, _, contributions = calc.calculate_mode_intensities(return_contributions=True)
        raw = sum(c.local_intensity for c in contributions)
        assert raw > 0
        np.testing.assert_allclose(intensity, [raw * (20000 - corrected) ** 4 / corrected], rtol=1e-10)


def test_bose_factor_limits_and_domain():
    """The harmonic Bose factor has a stable low-frequency and zero-temperature limit."""
    from PDielec.Constants import boltzmann_si, planck_si, speed_light_si
    from PDielec.RamanSpectrum import bose_factor

    assert bose_factor(100.0, 0.0) == pytest.approx(0.01)
    assert bose_factor(1e6, 1.0) == pytest.approx(1e-6)
    nu = 1e-5
    classical = boltzmann_si * 300 / (planck_si * speed_light_si * 100 * nu**2)
    assert bose_factor(nu, 300.0) == pytest.approx(classical, rel=1e-7)
    for invalid in [-1.0, 0.0, np.nan, np.inf]:
        with pytest.raises(ValueError):
            bose_factor(invalid, 300.0)
