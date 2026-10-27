"""Complete modal-pair NAC enumeration with dominant-TO GUI ownership."""

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np
import pytest

from PDielec.LayeredRamanCalculator import LayeredRamanCalculator, RamanLayer

from .conftest import build_system


def _spectrum(frequencies, tensors, nac_data, selected, depth, combination, parallel):
    """Run a matched slab with independently supplied NAC eigendata."""
    layer = RamanLayer(0, frequencies, tensors, np.eye(3),
                       nac_function=lambda _q: nac_data, modes_selected=selected)
    calc = LayeredRamanCalculator(
        build_system([(300e-9, 1.0)], n_sup=1.0, n_sub=1.0), [layer],
        laser_frequency_cm1=20000.0, incident_angle_rad=0.0,
        incident_pol="p", detected_pol="p", temperature_K=0.0,
        linewidths_cm1=[3.0, 7.0, 11.0], n_gauss=40, approximate_es=False,
        modal_pairs=True, modal_pair_use_nac=[True] * len(frequencies),
        depth_integration=depth, modal_pair_combination=combination,
    )
    if parallel:
        with ThreadPool(2) as pool:
            result = calc.calculate_mode_intensities(pool=pool)
    else:
        result = calc.calculate_mode_intensities(return_contributions=True)
        assert sum(c.intensity for c in result[3]) == pytest.approx(sum(result[1]), rel=1e-12)
    for seed, vectors in calc._modal_pair_q_vectors_by_branch.items():
        for branch, _pair in vectors:
            assert nac_data[4][branch] == seed
    return result[:3]


def _analytic_line(frequency, rxx, depth):
    """Exact vacuum overlap with finite-Stokes backscattering phase."""
    length = 3000.0  # angstrom
    weight = length
    if depth == "Coherent amplitude":
        weight = length**2 * np.sinc(100 * (40000 - frequency) * 300e-9)**2
    return abs(rxx)**2 * weight * (20000 - frequency)**4 / frequency


@pytest.mark.parametrize("selected", [[True, True, True], [True, False, False],
                                     [False, True, False], [False, False, True]])
@pytest.mark.parametrize("depth", ["Coherent amplitude", "Incoherent intensity"])
@pytest.mark.parametrize("parallel", [False, True])
def test_rank_one_nac_retains_every_selected_branch(selected, depth, parallel):
    """A physical positive rank-one correction must not lose its highest branch."""
    to_frequencies = np.array([100.0, 200.0, 300.0])
    vector = np.array([341.3150353965, -265.6599463666, 205.3665076454])
    eigenvalues, overlaps = np.linalg.eigh(np.diag(to_frequencies**2) + np.outer(vector, vector))
    frequencies = np.sqrt(eigenvalues)
    owners = np.argmax(abs(overlaps), axis=0)
    np.testing.assert_array_equal(owners, [0, 2, 0])
    tensors = np.array([np.diag(row) for row in np.eye(3)])
    nac_tensors = np.einsum("np,nij->pij", overlaps, tensors)
    widths = np.array([3.0, 7.0, 11.0])[owners]
    nac_data = (frequencies, nac_tensors, widths, None, owners)
    f, intensities, sigmas = _spectrum(to_frequencies, tensors, nac_data, selected,
                                     depth, "Group q channels", parallel)
    keep = np.asarray(selected)[owners]
    order = np.argsort(f)
    np.testing.assert_allclose(f[order], frequencies[keep], rtol=1e-12)
    np.testing.assert_allclose(sigmas[order], widths[keep], rtol=1e-12)
    expected = [_analytic_line(nu, tensor[0, 0], depth)
                for nu, tensor in zip(frequencies[keep], nac_tensors[keep])]
    np.testing.assert_allclose(intensities[order], expected, rtol=2e-11)


@pytest.mark.parametrize("combination", ["Group q channels", "Incoherent pairs", "Coherent all pairs"])
@pytest.mark.parametrize("depth", ["Coherent amplitude", "Incoherent intensity"])
@pytest.mark.parametrize("parallel", [False, True])
@pytest.mark.parametrize("branch_selected", [[True, True], [False, True]])
def test_degenerate_branches_add_powers_not_amplitudes(combination, depth, parallel, branch_selected):
    """Equal line frequencies do not make orthogonal phonon states coherent."""
    tensor = np.diag([1.0, 0.0, 0.0])
    data = ([250.0, 250.0], [tensor, -tensor], [5.0, 5.0], branch_selected, [0, 0])
    f, intensities, widths = _spectrum([100.0, 200.0], [tensor, tensor], data, [True, True],
                                     depth, combination, parallel)
    np.testing.assert_allclose(f, [250.0])
    np.testing.assert_allclose(widths, [5.0])
    assert sum(intensities) == pytest.approx(
        sum(branch_selected) * _analytic_line(250.0, 1.0, depth), rel=2e-11)
