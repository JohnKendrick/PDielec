"""Shared helpers for Crystal Raman tests (B–F).

Provides lightweight system-building utilities so each test module
does not need to repeat GTMcore boilerplate.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
from PDielec.GTMcore import Layer, TransferMatrixSystem
from PDielec.LayeredRamanCalculator import RamanLayer, LayeredRamanCalculator


def iso_eps(n):
    """Return a frequency-independent permittivity function for refractive index n.

    Parameters
    ----------
    n : complex or float
        Refractive index (ε = n²).

    Returns
    -------
    callable
        Function f(freq_cm1) → 3×3 complex ndarray.
    """
    eps_val = complex(n) ** 2
    return lambda freq_cm1: eps_val * np.eye(3, dtype=complex)


def make_layer(thickness_m, eps_func):
    """Create a GTMcore Layer with the coherent-mode flags required by
    ``calculate_GammaStar``.

    Plain ``Layer`` objects lack the ``inCoherentIntensity`` etc. flags that
    ``CoherentLayer`` sets; this helper adds them directly so the layer can be
    used as a finite slab in a ``TransferMatrixSystem`` without needing the GUI's
    ``SingleCrystalLayer`` adapter.

    Parameters
    ----------
    thickness_m : float
        Layer thickness in metres.
    eps_func : callable
        Function(freq_cm1) → 3×3 complex ndarray.

    Returns
    -------
    Layer
    """
    layer = Layer(thickness=thickness_m, epsilon=eps_func)
    layer.inCoherentIntensity = False
    layer.inCoherentPhase = False
    layer.inCoherentAveragePhase = False
    layer.inCoherentThick = False
    layer.SMatrix = None
    return layer


def build_system(layer_specs, n_sup=1.0, n_sub=1.0):
    """Build a TransferMatrixSystem from a list of (thickness_m, n) pairs.

    Parameters
    ----------
    layer_specs : list of (float, complex)
        Each entry is (thickness in metres, refractive index).
    n_sup, n_sub : float or complex
        Refractive indices for the semi-infinite superstrate and substrate.

    Returns
    -------
    TransferMatrixSystem
    """
    sup = Layer(thickness=1e-3, epsilon=iso_eps(n_sup))
    sub = Layer(thickness=1e-3, epsilon=iso_eps(n_sub))
    layers = [make_layer(d, iso_eps(n)) for d, n in layer_specs]
    return TransferMatrixSystem(substrate=sub, superstrate=sup, layers=layers)


def run_calc(
    system,
    raman_layers,
    phonon_freqs_cm1,
    incident_pol="p",
    detected_pol="p",
    laser_cm1=20000.0,
    temperature_K=0.0,
    n_gauss=20,
    approximate_es=True,
    coherent_layers=False,
    incident_angle_rad=0.0,
    collection_side="superstrate",
):
    """Run LayeredRamanCalculator and return (active_freqs, intensities).

    Parameters
    ----------
    system : TransferMatrixSystem
    raman_layers : list of RamanLayer
    phonon_freqs_cm1 : array_like
        All phonon frequencies including acoustic (> 10 cm⁻¹ are kept).
    """
    linewidths = 5.0 * np.ones(len(phonon_freqs_cm1))
    calc = LayeredRamanCalculator(
        system=system,
        raman_layers=raman_layers,
        laser_frequency_cm1=laser_cm1,
        incident_angle_rad=incident_angle_rad,
        incident_pol=incident_pol,
        detected_pol=detected_pol,
        temperature_K=temperature_K,
        linewidths_cm1=linewidths,
        n_gauss=n_gauss,
        approximate_es=approximate_es,
        coherent_layers=coherent_layers,
        collection_side=collection_side,
    )
    freqs, intensities, _ = calc.calculate_mode_intensities()
    return freqs, intensities


def make_raman_layer(layer_index, R_crystal, nu_cm1, G=None):
    """Convenience wrapper for RamanLayer with a single phonon mode.

    Parameters
    ----------
    layer_index : int
    R_crystal : ndarray, shape (3, 3)
        Raman tensor in the crystal frame.
    nu_cm1 : float
        Phonon frequency in cm⁻¹ (must be > 10 cm⁻¹ to be kept as active).
    G : ndarray or None
        Rotation matrix (identity if None).
    """
    if G is None:
        G = np.eye(3, dtype=float)
    return RamanLayer(
        layer_index=layer_index,
        phonon_frequencies_cm1=np.array([nu_cm1]),
        raman_tensors=[np.asarray(R_crystal, dtype=float)],
        rotation_matrix=G,
    )
