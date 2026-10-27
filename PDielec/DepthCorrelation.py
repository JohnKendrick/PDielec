"""Reduced depth weights for explicitly normalized phonon correlation domains.

Optical fields are unit-launch fields. Overlaps contain the physical R_epsilon
in sqrt(angstrom/amu); weights and domain lengths use angstrom. Both models
return angstrom**2/amu, before occupation/radiation factors. No absolute power
or cross section is implied. See eq-depth-domain-kernel in raman-development.tex.
"""

import numpy as np


def depth_weight(overlap, weights, domain_indices=None, domain_lengths=None):
    """Integrate local sources or one normalized collective mode per domain.

    Parameters
    ----------
    overlap : array_like, shape (channels, points)
        Bilinear reciprocal overlaps; detector channels add incoherently.
    weights : array_like, shape (points,)
        Nonnegative quadrature weights, in the same length unit as domains.
    domain_indices : array_like of int, optional
        Physical domain of each point. Omit for the local delta kernel.
    domain_lengths : array_like, optional
        Physical normalization lengths, independent of quadrature or optical
        subdivision. Domains must partition the source support. A zero length
        domain contributes zero and cannot contain positive-weight samples.

    Returns
    -------
    float
        Reduced weight with units overlap squared times length.

    """
    f = np.asarray(overlap, dtype=complex)
    w = np.asarray(weights, dtype=float)
    if f.ndim != 2 or w.ndim != 1 or f.shape[1] != w.size:
        raise ValueError("overlap must be (channels, points), weights must be (points,)")
    if not np.all(np.isfinite(f)) or not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("overlaps must be finite and weights finite and nonnegative")
    if domain_indices is None:
        if domain_lengths is not None:
            raise ValueError("domain_lengths requires domain_indices")
        return float(np.sum(np.abs(f) ** 2 * w))
    ids = np.asarray(domain_indices)
    lengths = np.asarray(domain_lengths, dtype=float)
    if (ids.shape != w.shape or not np.issubdtype(ids.dtype, np.integer)
            or lengths.ndim != 1 or not np.all(np.isfinite(lengths))
            or np.any(lengths < 0) or np.any(ids < 0) or np.any(ids >= lengths.size)):
        raise ValueError("invalid domain indices or lengths")
    result = 0.0
    for index, length in enumerate(lengths):
        selected = ids == index
        measure = np.sum(w[selected])
        if not np.isclose(measure, length, rtol=1e-10, atol=0.0):
            raise ValueError("quadrature must cover each physical domain exactly")
        if length > 0:
            amplitude = f[:, selected] @ w[selected]
            result += np.sum(np.abs(amplitude) ** 2) / length
    return float(result)
