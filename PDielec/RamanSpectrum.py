"""Shared harmonic Stokes frequency policy and spectral weighting."""

import logging

import numpy as np

from PDielec.Constants import boltzmann_si, planck_si, speed_light_si

logger = logging.getLogger(__name__)


def validate_laser_frequency(laser_cm1):
    """Require a finite positive excitation wavenumber before evaluating optics."""
    if not np.isfinite(laser_cm1) or laser_cm1 <= 0:
        raise ValueError("Raman laser frequency must be finite and positive")


def valid_stokes_mode(nu_cm1, laser_cm1, acoustic_cutoff=0.0):
    """Accept stable modes with a positive outgoing photon, regardless of field approximation."""
    validate_laser_frequency(laser_cm1)
    valid = np.isfinite(nu_cm1) and 0 < nu_cm1 < laser_cm1 and nu_cm1 >= acoustic_cutoff
    if not valid:
        logger.debug("Skipping Raman mode %s cm-1 outside the harmonic Stokes range", nu_cm1)
    return bool(valid)


def bose_factor(nu_cm1, temperature_K):
    """Return (n+1)/nu for a stable mode, including the zero-temperature limit (eq-bose)."""
    if not np.isfinite(nu_cm1) or nu_cm1 <= 0:
        raise ValueError("Bose factor requires a finite positive mode frequency")
    if not np.isfinite(temperature_K) or temperature_K < 0:
        raise ValueError("Raman temperature must be finite and nonnegative")
    if temperature_K == 0:
        return 1.0 / nu_cm1
    x = planck_si * speed_light_si * 100 * nu_cm1 / (boltzmann_si * temperature_K)
    # -expm1(-x) is accurate at small x and never overflows at large x.
    return 1.0 / (-np.expm1(-x) * nu_cm1)


def stokes_prefactor(nu_cm1, laser_cm1, temperature_K):
    """Return nu_scattered**4 * (n+1)/nu (eq-ramanefficiency_depolarised).

    Wavenumbers are in cm^-1, consistently for powder and crystal spectra.
    Absolute collection/cross-section constants remain outside this factor.
    """
    if not valid_stokes_mode(nu_cm1, laser_cm1):
        return 0.0
    return (laser_cm1 - nu_cm1) ** 4 * bose_factor(nu_cm1, temperature_K)
