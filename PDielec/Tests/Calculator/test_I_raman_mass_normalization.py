"""Mass-coordinate Raman normalization against Cartesian derivative oracles."""

import numpy as np
import pytest

from PDielec.GenericOutputReader import GenericOutputReader


def make_reader():
    """Create a complete synthetic mode basis and an independent Cartesian response."""
    r = GenericOutputReader([])
    r.nions = 2
    r.species = ["H", "He"]
    r.atom_type_list = [0, 1]
    r.masses = [2.0, 3.0]
    r.masses_per_type = r.masses.copy()
    rng = np.random.default_rng(37)
    u, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    r.mass_weighted_normal_modes = u.reshape(6, 2, 3)
    r.frequencies = [100.0, 100.0, 200.0, 300.0, 400.0, 500.0]
    response = rng.normal(size=(6, 3, 3)) + 1j * rng.normal(size=(6, 3, 3))
    r.raman_tensors = np.einsum("nk,kij->nij", u / np.repeat(np.sqrt(r.masses), 3), response)
    return r, response


@pytest.mark.parametrize("before_first", [True, False])
@pytest.mark.parametrize("edited", [True, False])
def test_mass_edits_and_reset_match_cartesian_projection(before_first, edited):
    """Changing normal coordinates includes inverse square root masses without drift."""
    r, response = make_reader()
    if not before_first:
        r.calculate_mass_weighted_normal_modes()
    for masses in ([8.0, 12.0], [2.0, 7.0], [2.0, 3.0], [2.0, 3.0]):
        if edited:
            r.set_edited_masses(list(masses))
        r.change_masses(dict(zip(r.species, masses)), {})
        r.calculate_mass_weighted_normal_modes()
        u = np.asarray(r.mass_weighted_normal_modes).reshape(6, 6)
        expected = np.einsum("nk,kij->nij", u / np.repeat(np.sqrt(masses), 3), response)
        np.testing.assert_allclose(r.raman_tensors, expected, atol=2e-12)
        if masses == [8.0, 12.0]:
            np.testing.assert_allclose(r.frequencies, [50.0, 50.0, 100.0, 150.0, 200.0, 250.0], atol=1e-9)


def test_incomplete_reference_cannot_support_mass_reprojection():
    """Missing modes must not silently yield a supposedly exact reprojected tensor."""
    r, _ = make_reader()
    with pytest.raises(ValueError, match="complete"):
        r._transform_raman_tensors_between_mode_bases(
            r.mass_weighted_normal_modes[:3],
            r.mass_weighted_normal_modes,
            r.raman_tensors[:3],
            old_masses=[2.0, 3.0],
            new_masses=[8.0, 12.0],
        )


@pytest.mark.parametrize("force_constants", [True, False])
def test_finite_field_native_to_average_mass_projection(tmp_path, force_constants):
    """JSON projected tensors use program masses, even when loading selects average masses."""
    import json

    from PDielec.Constants import amu, angs2bohr, hartree2ev, wavenumber
    from PDielec.FiniteFieldOutputReader import FiniteFieldOutputReader

    r, response = make_reader()
    # JSON is real-valued, so use the real Cartesian response.
    u = np.asarray(r.mass_weighted_normal_modes).reshape(6, 6)
    d = u.T @ np.diag((np.array(r.frequencies) * wavenumber) ** 2) @ u
    sqrt_mass = np.repeat(np.sqrt(np.array(r.masses) * amu), 3)
    phi = d * np.outer(sqrt_mass, sqrt_mass) * hartree2ev * angs2bohr**2
    data = {
        "schema": "pdielect-finite-field-raman-v1",
        "structure": {
            "lattice_angstrom": (np.eye(3) * 3).tolist(),
            "symbols": r.species,
            "positions_fractional": [[0, 0, 0], [0.5, 0.5, 0.5]],
        },
        "program_masses_amu": r.masses,
        "force_constants_matrix_sym_eV_per_angstrom2": phi.tolist(),
        "modes": [
            {"frequency_cm-1": f, "mass_weighted_eigenvector": v.tolist(), "raman_tensor": t.real.tolist()}
            for f, v, t in zip(r.frequencies, r.mass_weighted_normal_modes, r.raman_tensors)
        ],
    }
    if not force_constants:
        del data["force_constants_matrix_sym_eV_per_angstrom2"]
    path = tmp_path / "raman.json"
    path.write_text(json.dumps(data))
    reader = FiniteFieldOutputReader([str(path)])
    reader.read_output()
    for masses in (reader.masses, [8.0, 12.0], [2.0, 3.0]):
        reader.change_masses(dict(zip(reader.species, masses)), {})
        reader.calculate_mass_weighted_normal_modes()
        modes = np.asarray(reader.mass_weighted_normal_modes).reshape(6, 6)
        expected = np.einsum("nk,kij->nij", modes / np.repeat(np.sqrt(masses), 3), response.real)
        np.testing.assert_allclose(reader.raman_tensors, expected, atol=3e-12)

    # A projected-only file with no reference vectors cannot support even the
    # initial program-to-average mass change; do not silently invent phases.
    for mode in data["modes"]:
        del mode["mass_weighted_eigenvector"]
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="reference eigenvectors"):
        FiniteFieldOutputReader([str(path)]).read_output()


@pytest.mark.parametrize("program", ["Vasp", "AbInit"])
def test_native_reader_uniform_mass_scaling_and_reset(program):
    """A per-mode reader and a raw-derivative control both obey the mass-coordinate law."""
    from pathlib import Path

    from PDielec.AbinitOutputReader import AbinitOutputReader
    from PDielec.VaspOutputReader import VaspOutputReader

    cls, filename = (VaspOutputReader, "OUTCAR") if program == "Vasp" else (AbinitOutputReader, "raman.abo")
    path = Path(__file__).resolve().parents[3] / "Examples/Powder_Raman" / program / filename
    reader = cls([str(path)])
    reader.read_output()
    reader.calculate_mass_weighted_normal_modes()
    frequencies = np.array(reader.frequencies)
    strength = np.sum(np.abs(reader.raman_tensors)**2)
    reader.change_masses({k: 4*v for k, v in reader.mass_dictionary().items()}, {})
    reader.calculate_mass_weighted_normal_modes()
    np.testing.assert_allclose(reader.frequencies, frequencies / 2, atol=1e-7)
    assert np.sum(np.abs(reader.raman_tensors)**2) == pytest.approx(strength / 4, rel=1e-10)
    reader.reset_masses()
    reader.calculate_mass_weighted_normal_modes()
    np.testing.assert_allclose(reader.frequencies, frequencies, atol=1e-7)
    assert np.sum(np.abs(reader.raman_tensors)**2) == pytest.approx(strength, rel=1e-10)


@pytest.mark.parametrize("eckart", [False, True])
def test_aims_native_masses_survive_average_mass_change(tmp_path, eckart):
    """FHI-aims native per-atom masses are available before capturing the mode reference."""
    from PDielec.AimsOutputReader import AimsOutputReader
    from PDielec.Constants import amu, angs2bohr, average_masses, hartree2ev

    # Two isotopes of the same species distinguish per-atom provenance from a
    # species-mass dictionary. FHI-aims force constants are in eV/Angstrom^2.
    masses_file = tmp_path / "masses.dat"
    masses_file.write_text("1.0\n2.0\n")
    geometry = tmp_path / "geometry.in"
    geometry.write_text("lattice_vector 3 0 0\nlattice_vector 0 3 0\nlattice_vector 0 0 3\n"
                        "atom_frac 0 0 0 H\natom_frac .5 .5 .5 H\n")
    rng = np.random.default_rng(19)
    a = rng.normal(size=(6, 6))
    force_constants = a.T @ a + np.eye(6)
    hessian_file = tmp_path / "hessian.dat"
    np.savetxt(hessian_file, force_constants)
    reader = AimsOutputReader([str(tmp_path / "ignore.out"), str(masses_file),
                               str(geometry), str(hessian_file), str(tmp_path / "ignore.born")])
    reader.read_output()
    np.testing.assert_array_equal(reader.masses, [1., 2.])
    np.testing.assert_array_equal(reader._mode_reference_masses, [1., 2.])
    reader.change_masses(average_masses, {})
    reader.eckart = eckart
    reader.calculate_mass_weighted_normal_modes()
    expected_force_constants = force_constants / (hartree2ev * angs2bohr**2)
    np.testing.assert_allclose(reader.nomass_hessian, expected_force_constants, atol=1e-14)
    masses = np.repeat(np.asarray(reader.masses) * amu, 3)
    expected = expected_force_constants / np.sqrt(np.outer(masses, masses))
    if eckart:
        translations = np.tile(np.eye(3), (2, 1)) * np.sqrt(masses[:, None])
        translations /= np.linalg.norm(translations, axis=0)
        projector = np.eye(6) - translations @ translations.T
        expected = projector @ expected @ projector
    np.testing.assert_allclose(reader.hessian, expected, atol=1e-14)
    assert np.all(np.isfinite(reader.frequencies))


@pytest.mark.parametrize("storage", ["old_masses", "species_masses", "program_masses"])
def test_generic_reference_supports_native_mass_representations(storage):
    """Native references must not depend on a reader having initialized active masses."""
    reader, response = make_reader()
    reader.masses = []
    if storage == "old_masses":
        reader._old_masses = [2., 3.]
        # Native per-atom masses take precedence over a species-level fallback.
        reader.masses_per_type = [4., 5.]
    elif storage == "program_masses":
        reader.program_mass_dictionary = {"H": 2., "He": 3.}
        reader.masses_per_type = []
    reader.change_masses({"H": 8., "He": 12.}, {})
    reader.calculate_mass_weighted_normal_modes()
    np.testing.assert_allclose(reader.frequencies, [50., 50., 100., 150., 200., 250.], atol=1e-9)
    u = np.asarray(reader.mass_weighted_normal_modes).reshape(6, 6)
    expected = np.einsum("nk,kij->nij", u / np.repeat(np.sqrt([8., 12.]), 3), response)
    np.testing.assert_allclose(reader.raman_tensors, expected, atol=2e-12)


@pytest.mark.parametrize("invalid_masses", [[], [2.], [2., np.nan], [2., 0.]])
def test_generic_invalid_mass_reference_is_not_cached(invalid_masses):
    """Reject incomplete or unphysical provenance before setting Hessian/reference state."""
    reader, _ = make_reader()
    reader.masses = invalid_masses
    reader.masses_per_type = []
    with pytest.raises(ValueError, match="native.*mass|mass.*native"):
        reader.calculate_mass_weighted_normal_modes()
    assert reader._mode_reference_masses is None
    assert reader.nomass_hessian is None
    assert not reader.nomass_hessian_has_been_set


@pytest.mark.parametrize("sample", ["ZnO", "Na2SO42"])
@pytest.mark.parametrize("mass_definition", ["program", "average"])
def test_phonopy_bundled_mass_recalculation(sample, mass_definition):
    """Real Phonopy inputs retain their original dynamical matrix under mass conversion."""
    from pathlib import Path

    from PDielec.Constants import amu, average_masses, wavenumber
    from PDielec.Utilities import get_reader

    path = Path(__file__).resolve().parents[3] / "Examples/Phonopy" / sample / "phonopy.yaml"
    reader = get_reader(str(path), "phonopy")
    reader.read_output()
    native = np.asarray(reader._old_masses)
    np.testing.assert_array_equal(reader.masses, native)
    u = np.asarray(reader.mass_weighted_normal_modes).reshape(3*reader.nions, -1)
    f = np.asarray(reader.frequencies)*wavenumber
    dynamical = u.T @ np.diag(np.sign(f)*f**2) @ u
    m = np.repeat(native*amu, 3)
    force_constants = dynamical*np.sqrt(np.outer(m, m))
    if mass_definition == "average":
        reader.change_masses(average_masses, {})
    reader.eckart = True
    reader.calculate_mass_weighted_normal_modes()
    np.testing.assert_allclose(reader.nomass_hessian, force_constants, atol=1e-11)
    m = np.repeat(np.asarray(reader.masses)*amu, 3)
    translations = np.tile(np.eye(3), (reader.nions, 1))*np.sqrt(m[:, None])
    translations /= np.linalg.norm(translations, axis=0)
    p = np.eye(len(m)) - translations @ translations.T
    expected = p @ (force_constants/np.sqrt(np.outer(m, m))) @ p
    np.testing.assert_allclose(reader.hessian, expected, atol=1e-13)
