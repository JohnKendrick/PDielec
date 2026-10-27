"""Independent momentum, reciprocal-field and constant-response regressions."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from PDielec.LayeredRamanCalculator import LayeredRamanCalculator
from PDielec.OpticalChannelResolver import OpticalChannelResolver
from PDielec.PhononFinalStateResolver import PhononFinalStateResolver
from PDielec.RamanSpectrum import stokes_prefactor
from PDielec.Tests.Crystal_Raman.conftest import build_system, make_raman_layer

LASER = 20000.0
STOKES = 19500.0
K_L = 2 * np.pi * 100 * LASER
K_S = 2 * np.pi * 100 * STOKES
LENGTH = 100e-9


def calculator(side="superstrate", angle=0.0, modal=False, local=False,
               q_resolved=False, final_state_model="Bulk phase matched"):
    """Build a matched slab with a single constant Raman response."""
    return LayeredRamanCalculator(
        system=build_system([(LENGTH, 1.0)]),
        raman_layers=[make_raman_layer(0, np.diag([1, 1, 0]), LASER-STOKES)],
        laser_frequency_cm1=LASER, incident_angle_rad=angle,
        incident_pol="s", detected_pol="s", temperature_K=0,
        linewidths_cm1=[5], n_gauss=32, approximate_es=False,
        collection_side=side, modal_pairs=modal,
        modal_pair_use_nac=[q_resolved],
        depth_integration="Incoherent intensity" if local else "Coherent amplitude",
        final_state_model=final_state_model,
    )


def phase_slope(calc):
    """Fit the phase of the bilinear reciprocal field overlap."""
    data = calc.resolve_optical_channels(scattered_frequency_cm1=STOKES)
    # s polarization isolates a unit transverse field at normal/oblique incidence.
    overlap = data["incident"].total_fields[4] * data["detected"].total_fields[4]
    slope = np.polyfit(calc._gl_z, np.unwrap(np.angle(overlap)), 1)[0]
    return float(slope), data


@pytest.mark.parametrize("side", ["superstrate", "substrate"])
@pytest.mark.parametrize("angle", [0.0, np.pi/6])
def test_actual_depth_phase_has_reciprocal_sign_and_stokes_frequency(side, angle):
    """Actual depth phase has reciprocal sign and stokes frequency."""
    calc = calculator(side, angle)
    slope, _ = phase_slope(calc)
    sign = 1 if side == "superstrate" else -1
    expected = (K_L + sign*K_S) * np.cos(angle)
    assert_allclose(slope, expected, rtol=1e-10)


@pytest.mark.parametrize("side", ["superstrate", "substrate"])
def test_production_cached_pair_matches_actual_depth_phase(side):
    """Production cached pair matches actual depth phase."""
    calc = calculator(side, modal=True)
    calc.calculate_mode_intensities()
    # Positive-qz travelling channel is first; s launch is polarization index 1.
    cached = calc._modal_pair_q_vectors[(0, 0, 1, 0)].copy()
    slope, _ = phase_slope(calc)
    assert_allclose(cached[2] * K_L, slope, rtol=1e-10, atol=1e-6)


@pytest.mark.parametrize("side", ["superstrate", "substrate"])
def test_external_momentum_in_common_laser_units_keeps_stokes_shift(side):
    """External momentum in common laser units keeps stokes shift."""
    calc = calculator(side)
    detector = calc.system if side == "superstrate" else calc.system.reversed_system()
    actual = OpticalChannelResolver().external_momentum_transfer(
        calc.system, LASER, 0, detector, STOKES, 0, side)
    sign = 1 if side == "superstrate" else -1
    assert_allclose(actual, [0, 0, 1 + sign * STOKES/LASER], atol=1e-12)


@pytest.mark.parametrize("side", ["superstrate", "substrate"])
def test_reciprocal_launch_opposes_documented_outgoing_tangential_component(side):
    """Reciprocal launch opposes documented outgoing tangential component."""
    calc = calculator(side, angle=np.pi/6)
    _, data = phase_slope(calc)
    # Documented collection angle labels outgoing x component; reciprocity reverses it.
    expected_reciprocal_kx = -K_S * np.sin(calc.collection_angle_rad)
    assert_allclose(K_S * data["detected"].zeta.real, expected_reciprocal_kx, atol=1e-8)


@pytest.mark.parametrize("modal", [False, True])
def test_nonpolar_local_path_matches_complete_state_sum(modal):
    """Nonpolar local path matches complete state sum."""
    calc = calculator(modal=modal, local=True)
    raw = calc.calculate_mode_intensities()[1][0] / stokes_prefactor(500, LASER, 0)
    # Equal-weight complete orthonormal state sum of a unit-magnitude overlap = L.
    assert_allclose(raw, LENGTH*1e10, rtol=1e-11)


def test_bulk_filter_preserves_matched_slab_constant_response_limit():
    # Deliberately enable q-resolved grouping with a constant tensor/frequency.
    # This is a bookkeeping limit, not a polar-material prediction.
    """Bulk filter preserves matched slab constant response limit."""
    calc = calculator(modal=True, local=True, q_resolved=True)
    raw = calc.calculate_mode_intensities()[1].sum() / stokes_prefactor(500, LASER, 0)
    assert_allclose(raw, LENGTH*1e10, rtol=1e-11)


def test_local_model_preserves_same_constant_response_limit():
    """Local model preserves same constant response limit."""
    calc = calculator(modal=True, local=True, q_resolved=True,
                      final_state_model="Local incoherent")
    raw = calc.calculate_mode_intensities()[1].sum() / stokes_prefactor(500, LASER, 0)
    assert_allclose(raw, LENGTH*1e10, rtol=1e-11)


def test_distinct_same_direction_periodic_states_are_orthogonal():
    """Distinct same direction periodic states are orthogonal."""
    cells = np.arange(64)
    u1 = np.exp(2j*np.pi*cells/64)/8
    u2 = np.exp(4j*np.pi*cells/64)/8
    assert_allclose(np.vdot(u1, u2), 0, atol=1e-14)
    resolver = PhononFinalStateResolver()
    # Establish the collision without treating angular keys as erroneous NAC keys.
    assert resolver.q_class_key([0, 0, 1]) == resolver.q_class_key([0, 0, 2])


def test_equal_weight_complete_sum_is_invariant_under_basis_rotation():
    """Equal weight complete sum is invariant under basis rotation."""
    rng = np.random.default_rng(271027)
    coefficients = rng.normal(size=8) + 1j*rng.normal(size=8)
    rotation, _ = np.linalg.qr(rng.normal(size=(8, 8)) + 1j*rng.normal(size=(8, 8)))
    assert_allclose(np.linalg.norm(rotation.conj().T @ coefficients)**2,
                    np.linalg.norm(coefficients)**2, rtol=1e-13)
    # A pair of vectors in one degenerate optical subspace reconstructs unchanged fields.
    basis = rng.normal(size=(3, 2)) + 1j*rng.normal(size=(3, 2))
    optical_rotation, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
    amplitudes = coefficients[:2]
    assert_allclose((basis @ optical_rotation) @ (optical_rotation.conj().T @ amplitudes),
                    basis @ amplitudes, atol=1e-13)


def test_production_degenerate_optical_group_survives_basis_rotation():
    """Production degenerate optical group survives basis rotation."""
    from PDielec.LayeredRamanCalculator import _modal_q_channels_from_fields

    rng = np.random.default_rng(271028)
    basis = rng.normal(size=(3, 2)) + 1j*rng.normal(size=(3, 2))
    amps = np.array([1+0.2j, -0.3+0.7j])
    rotation, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
    fields = []
    for vectors, coefficients in ((basis, amps), (basis @ rotation, rotation.conj().T @ amps)):
        modal = np.zeros((4, 2, 3, 5), dtype=complex)
        for i in range(2):
            modal[i, 0] = (vectors[:, i] * coefficients[i])[:, None]
        grouped = _modal_q_channels_from_fields(modal, {0: np.array([1, 1, -1, -1])},
                                               [0], [slice(0, 5)], 0)
        fields.append(grouped[0][0]["field"])
    assert_allclose(fields[0], fields[1], atol=1e-13)


@pytest.mark.parametrize("side", ["superstrate", "substrate"])
def test_reciprocal_p_vector_is_in_laboratory_frame(side):
    """Unit p vectors are transverse to the reversed physical outgoing ray."""
    angle = np.pi/6
    calc = calculator(side, angle)
    data = calc.resolve_optical_channels(STOKES)
    field = data["detected"].total_fields[:3]
    # Automatic retro and forward both have p vector (cos(theta),0,-sin(theta)).
    assert_allclose(field[2]/field[0], -np.tan(angle), atol=1e-12)
    assert_allclose(field[1], 0, atol=1e-12)


def test_reflected_anisotropic_modes_satisfy_original_lab_maxwell_equation():
    """Tensor reflection and electric-vector return transform must agree."""
    calc = calculator("substrate", np.pi/6)
    epsilon = np.array([[2.8, 0.2, 0.4], [0.2, 2.1, 0.3], [0.4, 0.3, 3.6]])
    layer = calc.system.layers[0]
    layer.epsilon_tensor_function = lambda nu: epsilon.copy()
    euler = layer.euler.copy()
    data = calc.resolve_optical_channels(STOKES)
    solution = data["detected"]
    reflected = data["collection_system"].layers[0]
    reflection = np.diag([1, 1, -1])
    assert reflected is not layer
    assert_allclose(layer.euler, euler)
    assert_allclose(reflected.epsilon, reflection @ epsilon @ reflection)
    for mode, qz in enumerate(solution.qs_by_layer[0]):
        k = np.array([solution.zeta, 0, -qz])
        operator = np.dot(k, k)*np.eye(3) - np.outer(k, k) - epsilon
        field = solution.modal_fields[mode].reshape(2, 3, -1)
        residual = np.einsum("ij,pjn->pin", operator, field)
        assert np.linalg.norm(residual) < 1e-10 * max(1, np.linalg.norm(field))


@pytest.mark.parametrize("modal", [False, True])
@pytest.mark.parametrize("approximate", [False, True])
def test_forward_worker_matches_lab_frame_serial_result(modal, approximate):
    """Both worker field reconstruction routes retain the p-vector transform."""
    from multiprocessing.pool import ThreadPool

    def make():
        calc = calculator("substrate", np.pi/6, modal=modal, local=True)
        calc.incident_pol = calc.detected_pol = "p"
        calc._incident_jones = calc._detected_jones = np.array([1, 0], dtype=complex)
        calc.approximate_es = approximate
        calc.raman_layers[0].raman_tensors = [np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])]
        return calc

    serial = make().calculate_mode_intensities()
    with ThreadPool(2) as pool:
        parallel = make().calculate_mode_intensities(pool=pool)
    assert_allclose(parallel, serial, rtol=1e-11)
    expected_raw = 1000 * np.sin(np.pi/3)**2
    assert_allclose(serial[1]/stokes_prefactor(500, LASER, 0), [expected_raw], rtol=1e-11)


def test_each_branch_and_worker_uses_returned_stokes_frequency():
    """NAC-seed policy is explicit; returned-frequency momentum is not cached at laser."""
    from multiprocessing.pool import ThreadPool

    def make():
        calc = calculator("substrate", modal=True, q_resolved=True)
        layer = calc.raman_layers[0]
        layer.phonon_frequencies_cm1 = np.array([500, 1500])
        layer.raman_tensors = [np.eye(3), np.eye(3)]
        layer.nac_function = lambda q: (np.array([700, 1700]), [np.eye(3), np.eye(3)], [5, 5])
        calc.modal_pair_use_nac = np.array([True, True])
        calc.linewidths_cm1 = np.array([5, 5])
        return calc

    serial = make()
    reference = serial.calculate_mode_intensities()
    parallel = make()
    with ThreadPool(2) as pool:
        actual = parallel.calculate_mode_intensities(pool=pool)
    assert_allclose(actual, reference, rtol=1e-11)
    for mode, frequency in enumerate([700, 1700]):
        for calc in (serial, parallel):
            q = calc._modal_pair_q_vectors_by_mode[mode][(0, 0, 1, 0)]
            assert_allclose(q, [0, 0, frequency/LASER], atol=1e-12)


def test_nac_channel_split_is_rejected_instead_of_misassigned():
    """A changed optical subspace count needs branch tracking, not ordinal reuse."""
    from PDielec.LayeredRamanCalculator import _validate_modal_channel_counts

    with pytest.raises(ValueError, match="tracking"):
        _validate_modal_channel_counts({0: {0: [{}, {}]}}, {0: {0: [{}, {}, {}]}})


@pytest.mark.parametrize("retro", [True, False])
def test_oblique_coherent_p_tensor_contraction_distinguishes_collection_rays(retro):
    """Retro and specular rays have different xx/zz interference at fixed frequency."""
    calc = calculator(angle=np.pi/4)
    calc.incident_pol = calc.detected_pol = "p"
    calc._incident_jones = calc._detected_jones = np.array([1, 0], dtype=complex)
    calc.collection_angle_rad = -np.pi/4 if retro else np.pi/4
    calc.raman_layers[0].raman_tensors = [np.diag([1.0, 0.0, -0.3])]
    frequencies, intensities, _ = calc.calculate_mode_intensities()
    # Unit p fields: (cos(theta), 0, -sin(theta)) for the incident wave.
    # Reciprocal p is the same for retro; its z component flips for specular.
    contraction = 0.5 - 0.15 if retro else 0.5 + 0.15
    phase_span = (K_L + K_S) * np.cos(np.pi/4) * LENGTH
    expected = (LENGTH * 1e10 * contraction * np.sinc(phase_span/(2*np.pi)))**2
    assert_allclose(frequencies, [LASER-STOKES])
    assert_allclose(intensities/stokes_prefactor(500, LASER, 0), [expected], rtol=1e-11)
