"""Independent SI Maxwell and modal-flux regression tests (Stage 2)."""
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.linalg import expm

from PDielec.Constants import epsilon_0_si, speed_light_si
from PDielec.GTMcore import CoherentLayer, Layer, ScatteringMatrixSystem, SemiInfiniteLayer, TransferMatrixSystem

FREQUENCY = 5e14
SYSTEMS = [TransferMatrixSystem, ScatteringMatrixSystem]


def layer(epsilon, thickness=180e-9, exterior=False, mu=1):
    """Construct a coherent layer with constant relative material properties."""
    adapter = SimpleNamespace(get_thickness_in_metres=lambda: thickness,
                              get_permittivity_function=lambda: lambda f: np.asarray(epsilon, dtype=complex))
    result = (SemiInfiniteLayer if exterior else CoherentLayer)(adapter)
    result.mu = mu
    return result


def anisotropic(loss=0):
    """Return a rotated biaxial dielectric, optionally passive."""
    c, s = np.cos(.47), np.sin(.47)
    rotation = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return rotation @ np.diag([2.25, 4., 3.]) @ rotation.T + 1j * loss * np.eye(3)


def tilted_absorbing():
    """Return a passive tensor with unequal cross-polarized reflection powers."""
    rotation, _ = np.linalg.qr(np.array([[1., 2., 3.], [3., 1., 2.], [2., 3., 1.]]))
    return rotation @ np.diag([2.25+.1j, 4.+.3j, 3.+.2j]) @ rotation.T


@pytest.mark.parametrize("system_type", SYSTEMS)
@pytest.mark.parametrize("frontend", ["helper", "gui"])
def test_absorbing_interface_incident_polarization(monkeypatch, system_type, frontend):
    """Both spectrum frontends must group reflected modes by incident channel."""
    from PDielec import HelperRoutines
    from PDielec.GUI import CrystalScenarioTab

    system = system_type(superstrate=layer(np.eye(3), exterior=True),
                         substrate=layer(tilted_absorbing(), exterior=True), layers=[])
    zeta = .3
    system.calculate_GammaStar(FREQUENCY, zeta)
    r, reflectance, t, transmittance = system.calculate_r_t(zeta)
    # r is [p->p, p->s, s->s, s->p]. Avoid the intensity-array
    # aggregation being tested when constructing the independent expectation.
    expected = np.sum(abs(r.reshape(2, 2))**2, axis=1)
    assert abs(reflectance[2]-reflectance[3]) > 1e-4
    assert_allclose(expected + transmittance, 1, atol=2e-12)
    result = (1., r, reflectance, t, transmittance, system.substrate.epsilon, 0, 0.)
    if frontend == "helper":
        monkeypatch.setattr(HelperRoutines, "solve_single_crystal_equations", lambda *args: result)
        reflected, transmitted, absorbed = HelperRoutines.calculate_crystal_infrared_spectrum([1.], [], 0., 0.)
    else:
        notebook = SimpleNamespace(pool=SimpleNamespace(imap=lambda *args, **kwargs: [result]),
                                   progressbars_update=lambda: None)
        tab = SimpleNamespace(settings={"Legend": "test", "Mode": "Transfer matrix"},
                              notebook=notebook, vs_cm1=[1.])
        output = CrystalScenarioTab.CrystalScenarioTab.coherent_calculator(
            tab, [], "Transfer matrix", 0., 0., 0., 0., 700)
        reflected, transmitted, absorbed = np.asarray(output[:6]).reshape(3, 2, 1)
    assert_allclose(reflected[:, 0], expected, atol=2e-12)
    assert_allclose(transmitted[:, 0], transmittance, atol=2e-12)
    assert_allclose(absorbed, 0, atol=2e-12)


@pytest.mark.parametrize("system_type", SYSTEMS)
def test_depth_absorption_uses_incident_reflection_channels(system_type):
    """Depth-resolved absorption agrees with directly normalized Poynting loss."""
    system = system_type(superstrate=layer(np.eye(3), exterior=True),
                         substrate=layer(tilted_absorbing(), exterior=True),
                         layers=[layer(tilted_absorbing(), thickness=300e-9)])
    zeta = .3
    z = np.array([-100e-9, 0., 100e-9, 200e-9])
    _, electric, magnetic, _ = system.calculate_Efield(FREQUENCY, zeta, z_vect=z, magnetic=True)
    _, reflectance, _, _ = system.calculate_r_t(zeta)
    flux, absorption = system.calculate_Poynting_Absorption_vs_z(z, electric, magnetic, reflectance)
    incident_flux = .5*np.sqrt(1-zeta**2)
    expected = (flux[[2, 5], :1]-flux[[2, 5], 2:])/incident_flux
    assert np.all(expected > 0)
    assert_allclose(absorption[:, 2:], expected, atol=2e-12)


def si_generator(epsilon, mu, zeta):
    """Eliminate longitudinal fields directly in SI curl equations.

    The state is (Ex, Hy, Ey, -Hx), deliberately with mixed SI units.
    """
    omega = 2 * np.pi * FREQUENCY
    mu0 = 1 / (epsilon_0_si * speed_light_si**2)
    kx = omega / speed_light_si * zeta

    def derivative(state):
        ex, hy, ey, minus_hx = state
        ez = -(epsilon[2, 0]*ex + epsilon[2, 1]*ey + kx*hy/(omega*epsilon_0_si))/epsilon[2, 2]
        hz = kx*ey/(omega*mu0*mu)
        return [1j*kx*ez + 1j*omega*mu0*mu*hy,
                1j*omega*epsilon_0_si*np.dot(epsilon[0], [ex, ey, ez]),
                1j*omega*mu0*mu*minus_hx,
                1j*omega*epsilon_0_si*np.dot(epsilon[1], [ex, ey, ez]) - 1j*kx*hz]

    return np.column_stack([derivative(column) for column in np.eye(4)])


@pytest.mark.parametrize("epsilon,zeta", [(anisotropic(.12), .4), (np.diag([-2., .5, 2.]), 1.)])
def test_complex_modes_are_finite_unit_vectors(epsilon, zeta):
    """A nonzero electric vector (1,0,-i) must not have zero norm."""
    medium = Layer(epsilon=lambda f: epsilon)
    medium.update_tm(FREQUENCY, zeta)
    assert np.isfinite(medium.Ai).all()
    assert_allclose(np.linalg.norm(medium.gamma, axis=1), 1, atol=2e-14)


@pytest.mark.parametrize("system_type", SYSTEMS)
@pytest.mark.parametrize("mu", [1., 1.3])
def test_birefringent_transmission_conserves_flux(system_type, mu):
    """Each transmitted mode has its own wavevector and magnetic field."""
    system = system_type(superstrate=layer(np.eye(3), exterior=True),
                         substrate=layer(anisotropic(), exterior=True, mu=mu), layers=[])
    system.calculate_GammaStar(FREQUENCY, .3)
    _, reflectance, _, transmittance = system.calculate_r_t(.3)
    assert_allclose(transmittance + reflectance[:2] + reflectance[[3, 2]], 1, atol=2e-12)


@pytest.mark.parametrize("mu", [1., 1.3])
def test_scaled_state_matches_si_maxwell(mu):
    """Check SI similarity transform, Ai permutation and propagation signs."""
    epsilon = anisotropic(.15)
    epsilon[0, 2] = epsilon[2, 0] = .2
    medium = layer(epsilon, mu=mu)
    medium.update_tm(FREQUENCY, .37)
    z0 = 1 / (epsilon_0_si * speed_light_si)
    scaling = np.diag([1, z0, 1, z0])
    generator = scaling @ si_generator(epsilon, mu, .37) @ np.linalg.inv(scaling)
    k0 = 2*np.pi*FREQUENCY/speed_light_si
    assert_allclose(generator/(1j*k0), medium.Delta, atol=2e-14)
    basis = medium.Ai[[0, 2, 1, 3], :]
    assert_allclose(generator @ basis, basis @ np.diag(1j*k0*medium.qs), atol=3e-8)
    assert_allclose(expm(generator*medium.thick),
                    basis @ np.diag(np.exp(1j*k0*medium.qs*medium.thick)) @ np.linalg.inv(basis), atol=2e-13)
    assert_allclose(medium.Ki, np.diag(np.exp(-1j*k0*medium.qs*medium.thick)), atol=2e-14)


def isotropic_basis(index, zeta):
    """Analytic tangential basis in (Ex, Z0 Hy, Ey, -Z0 Hx) order."""
    q = np.sqrt(index**2-zeta**2)
    electric = np.array([[q/index, 0, -zeta/index], [0, 1, 0],
                         [-q/index, 0, -zeta/index], [0, 1, 0]])
    wavevectors = np.array([[zeta, 0, q], [zeta, 0, q], [zeta, 0, -q], [zeta, 0, -q]])
    magnetic = np.cross(wavevectors, electric)
    return np.array([electric[:, 0], magnetic[:, 1], electric[:, 1], -magnetic[:, 0]])


@pytest.mark.parametrize("system_type", SYSTEMS)
@pytest.mark.parametrize("loss", [0., .2])
@pytest.mark.parametrize("tilted", [False, True])
def test_fields_against_direct_si_boundary_solve(system_type, loss, tilted):
    """Solve Maxwell propagation without the GTM eigenvectors or transfer matrix."""
    thickness, zeta = 180e-9, .4
    epsilon = anisotropic(loss)
    if tilted:
        rotation, _ = np.linalg.qr(np.array([[1., 2., 3.], [3., 1., 2.], [2., 3., 1.]]))
        epsilon = rotation @ epsilon @ rotation.T
    system = system_type(superstrate=layer(np.eye(3), exterior=True),
                         substrate=layer(1.7**2*np.eye(3), exterior=True),
                         layers=[layer(epsilon, thickness)])
    z = np.linspace(0, thickness, 13)
    _, electric, magnetic, _ = system.calculate_Efield(FREQUENCY, zeta, z_vect=z, magnetic=True)
    z0 = 1/(epsilon_0_si*speed_light_si)
    scaling = np.diag([1, z0, 1, z0])
    generator = scaling @ si_generator(epsilon, 1, zeta) @ np.linalg.inv(scaling)
    forward = expm(generator*thickness)
    incident, substrate = isotropic_basis(1, zeta), isotropic_basis(1.7, zeta)
    unknown = np.linalg.solve(np.column_stack([-forward@incident[:, 2:], substrate[:, :2]]),
                              forward@incident[:, :2])
    front = incident[:, :2] + incident[:, 2:] @ unknown[:2]
    for j, position in enumerate(z):
        expected = expm(generator*position) @ front
        actual = np.array([electric[[0, 3], j], magnetic[[1, 4], j],
                           electric[[1, 4], j], -magnetic[[0, 3], j]])
        assert_allclose(actual, expected, atol=2e-12)
    _, reflectance, transmission, transmittance = system.calculate_r_t(zeta)
    assert_allclose(transmission.reshape(2, 2).T, unknown[2:], atol=2e-12)
    flux = .5*np.real(np.cross(electric.T.reshape(-1, 2, 3), magnetic.T.reshape(-1, 2, 3).conj()))[:, :, 2]
    incident_flux = .5*np.sqrt(1-zeta**2)
    absorption = 1-reflectance[:2]-reflectance[[3, 2]]-transmittance
    assert_allclose((flux[0]-flux[-1])/incident_flux, absorption, atol=3e-12)
    if loss:
        assert np.all(absorption > 0)
        assert np.all(np.diff(flux, axis=0) < 0)
    else:
        assert_allclose(absorption, 0, atol=3e-12)
    amplitudes, _ = system.calculate_modal_amplitudes(FREQUENCY, zeta)
    medium = system.layers[0]
    k0 = 2*np.pi*FREQUENCY/speed_light_si
    for j, position in enumerate(z[1:], 1):
        modes = amplitudes[0].reshape(2, 4)*np.exp(1j*k0*medium.qs*(position-thickness))
        assert_allclose((modes@medium.gamma).ravel(), electric[:, j], atol=2e-12)


@pytest.mark.parametrize("system_type", SYSTEMS)
@pytest.mark.parametrize("gauge", ["rescale", "degenerate_rotation", "subdivide"])
def test_internal_basis_and_subdivision_invariance(system_type, gauge):
    """Internal mode coordinates must not change observable fields."""
    epsilon = 2.2**2*np.eye(3) if gauge == "degenerate_rotation" else anisotropic(.1)
    thickness, zeta = 180e-9, .4

    def make(layers):
        return system_type(superstrate=layer(np.eye(3), exterior=True),
                           substrate=layer(1.7**2*np.eye(3), exterior=True), layers=layers)

    reference = make([layer(epsilon, thickness)])
    candidate = make([layer(epsilon, thickness/3) for _ in range(3)] if gauge == "subdivide"
                     else [layer(epsilon, thickness)])
    if gauge != "subdivide":
        medium = candidate.layers[0]
        original = medium.calculate_gamma

        def changed_basis(zeta):
            original(zeta)
            if gauge == "rescale":
                medium.gamma *= np.array([2j, -.3+.2j, 1.1-.4j, -.7j])[:, None]
            else:
                rotation = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
                medium.gamma[:2] = rotation @ medium.gamma[:2]
                medium.gamma[2:] = rotation @ medium.gamma[2:]

        medium.calculate_gamma = changed_basis
    z = np.linspace(0, thickness, 19)
    expected = reference.calculate_Efield(FREQUENCY, zeta, z_vect=z, magnetic=True)
    actual = candidate.calculate_Efield(FREQUENCY, zeta, z_vect=z, magnetic=True)
    assert_allclose(actual[1], expected[1], atol=4e-12)
    assert_allclose(actual[2], expected[2], atol=4e-12)
    for result, baseline in zip(candidate.calculate_r_t(zeta), reference.calculate_r_t(zeta)):
        assert_allclose(result, baseline, atol=4e-12)


@pytest.mark.parametrize("system_type", SYSTEMS)
@pytest.mark.parametrize("n1,n2,zeta", [(1., 1., .6), (1., 1.5, .6), (1.5, 1., 1.2)])
def test_oblique_fresnel_and_total_internal_reflection(system_type, n1, n2, zeta):
    """Independent Fresnel powers include a purely evanescent substrate."""
    system = system_type(superstrate=layer(n1**2*np.eye(3), exterior=True),
                         substrate=layer(n2**2*np.eye(3), exterior=True), layers=[])
    system.calculate_GammaStar(FREQUENCY, zeta)
    _, reflectance, _, transmittance = system.calculate_r_t(zeta)
    q1, q2 = np.sqrt(complex(n1**2-zeta**2)), np.sqrt(complex(n2**2-zeta**2))
    rp = (n2**2*q1-n1**2*q2)/(n2**2*q1+n1**2*q2)
    rs = (q1-q2)/(q1+q2)
    assert_allclose(reflectance, [abs(rp)**2, abs(rs)**2, 0, 0], atol=2e-12)
    assert_allclose(transmittance, 1-np.array([abs(rp)**2, abs(rs)**2]), atol=2e-12)
