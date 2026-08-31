"""Tests J1–J8: Electro-optic (EO) correction to NAC Raman tensors.

Covers:
  J1  GenericOutputReader attribute and accessor.
  J2  Abinit reader parses χ^(2) = 2d from raman.abo correctly.
  J3  CASTEP reader parses χ^(2) = 2d from raman.castep correctly.
  J4  Backward compatibility — chi2_repsilon=None leaves tensors unchanged.
  J5  Zero χ^(2) tensor → EO correction is exactly zero.
  J6  Non-zero χ^(2) → at least one polar mode Raman tensor changes.
  J7  EO correction scales linearly with χ^(2) magnitude.
  J8  Zero Born charges → EO correction is zero for all modes.
  J9  QE reader parses χ^(2) from the electro-optic tensor block.
"""
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import math

import numpy as np
import pytest

from PDielec.Constants import amu, angs2bohr
from PDielec.GenericOutputReader import GenericOutputReader, chi2_pm_per_v_to_repsilon
from PDielec.GUI.CrystalScenarioTab import _compute_nac_dynamical_matrix_standalone

# ---------------------------------------------------------------------------
# Paths to DFT example files
# ---------------------------------------------------------------------------

_REPO = os.path.join(os.path.dirname(__file__), "..", "..", "..")
_ABINIT_FILE = os.path.join(_REPO, "Examples", "Crystal_Raman", "AbInit", "raman.abo")
_POWDER_ABINIT_FILE = os.path.join(_REPO, "Examples", "Powder_Raman", "AbInit", "raman.abo")
_CASTEP_FILE = os.path.join(_REPO, "Examples", "Crystal_Raman", "Castep", "raman.castep")
_CRYSTAL23_FILE = os.path.join(_REPO, "Examples", "Crystal_Raman", "Crystal23", "opt_raman.out")
_QE_TENSORS_FILE = os.path.join(_REPO, "Examples", "Powder_Raman", "QE", "tensors.xml")

_have_abinit = os.path.exists(_ABINIT_FILE)
_have_powder_abinit = os.path.exists(_POWDER_ABINIT_FILE)
_have_castep = os.path.exists(_CASTEP_FILE)
_have_crystal23 = os.path.exists(_CRYSTAL23_FILE)
_have_qe_tensors = os.path.exists(_QE_TENSORS_FILE)

# ---------------------------------------------------------------------------
# Shared helper: load Abinit reader and extract standalone-function inputs
# ---------------------------------------------------------------------------

def _load_abinit():
    """Return an AbinitOutputReader loaded from raman.abo."""
    from PDielec.AbinitOutputReader import AbinitOutputReader
    r = AbinitOutputReader([_ABINIT_FILE])
    r.read_output()
    return r


def _standalone_inputs(r):
    """Extract all inputs for _compute_nac_dynamical_matrix_standalone from a reader."""
    nAtoms = r.nions
    n_modes = 3 * nAtoms
    volume_au = r.volume * angs2bohr ** 3
    masses_au = np.array(r.masses) * amu
    born_charges = np.array(r.born_charges)
    hessian = np.array(r.hessian, dtype=float)
    eps_inf = np.array(r.zerof_optical_dielectric, dtype=float)
    if eps_inf.ndim == 1:
        eps_inf = np.diag(eps_inf)

    n_to_modes = len(r.mass_weighted_normal_modes)
    U_TO = np.zeros((n_to_modes, n_modes))
    for imode, mode in enumerate(r.mass_weighted_normal_modes):
        col = 0
        for atom in mode:
            U_TO[imode, col:col + 3] = atom
            col += 3

    raman_tensors = [np.asarray(R, dtype=float) for R in r.raman_tensors]
    sigmas = np.ones(n_to_modes) * 5.0
    return hessian, born_charges, eps_inf, volume_au, masses_au, U_TO, raman_tensors, sigmas


# ---------------------------------------------------------------------------
# J1: GenericOutputReader attribute initialisation and accessor
# ---------------------------------------------------------------------------

# Note: GenericOutputReader.__init__ returns early when the file doesn't exist,
# so J1 uses the Abinit reader with a real file to exercise a fully initialised
# instance.  Tests that only need the accessor are pure class-level checks.

@pytest.mark.skipif(not _have_abinit, reason="raman.abo not present")
class TestJ1Attribute:
    """J1: nonlinear_optical_susceptibility attribute and accessor."""

    def test_attribute_initialised_to_none(self):
        """A freshly loaded reader has chi2 = None before any NLO block is parsed.

        We exploit the fact that a reader for a file without an NLO block will
        initialise the attribute to None.  We verify via the Abinit reader which
        calls GenericOutputReader.__init__ with a valid file.
        """
        from PDielec.AbinitOutputReader import AbinitOutputReader
        r = AbinitOutputReader([_ABINIT_FILE])
        # Before read_output() the attribute must already be None
        assert r.nonlinear_optical_susceptibility is None

    def test_accessor_returns_none_before_read(self):
        from PDielec.AbinitOutputReader import AbinitOutputReader
        r = AbinitOutputReader([_ABINIT_FILE])
        assert r.get_nonlinear_optical_susceptibility() is None

    def test_accessor_returns_set_value(self):
        from PDielec.AbinitOutputReader import AbinitOutputReader
        r = AbinitOutputReader([_ABINIT_FILE])
        chi2 = np.zeros((3, 3, 3))
        chi2[0, 0, 0] = 5.0
        r.nonlinear_optical_susceptibility = chi2
        result = r.get_nonlinear_optical_susceptibility()
        assert result is chi2


class TestJ1BChi2UnitConversion:
    """J1B: pm/V conversion follows the documented cell-volume convention."""

    def test_batio3_reference_volume(self):
        """Check the SI/atomic-unit conversion used by the Stage 3 reference."""
        assert chi2_pm_per_v_to_repsilon(61.545914) == pytest.approx(0.7836682249873389, rel=1e-13)

    def test_inverse_square_root_volume_scaling(self):
        """The internal chi(2) factor must scale as inverse square root volume."""
        factor = chi2_pm_per_v_to_repsilon(50.0)
        assert chi2_pm_per_v_to_repsilon(200.0) == pytest.approx(0.5 * factor, rel=1e-15)

    def test_nonpositive_volume_rejected(self):
        """Conversion without a physical cell volume is undefined."""
        with pytest.raises(ValueError, match="positive cell volume"):
            chi2_pm_per_v_to_repsilon(0.0)


# ---------------------------------------------------------------------------
# J2: Abinit reader parses χ^(2) from raman.abo
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_abinit, reason="raman.abo not present")
class TestJ2AbinitParser:
    """J2: Abinit _read_nlo_susceptibility reads χ^(2) = 2d correctly."""

    @pytest.fixture(scope="class")
    def chi2(self):
        return _load_abinit().nonlinear_optical_susceptibility

    def test_tensor_not_none(self, chi2):
        assert chi2 is not None

    def test_shape(self, chi2):
        assert chi2.shape == (3, 3, 3)

    def test_zz_component(self, chi2):
        # raman.abo raw d[2,2,2] = -42.296970629; the reader stores
        # χ^(2)=2d after converting from pm/V to R_epsilon units.
        factor = chi2_pm_per_v_to_repsilon(_load_abinit().volume)
        assert chi2[2, 2, 2] == pytest.approx(-84.593941258 * factor, abs=1e-8)

    def test_last_two_indices_symmetric(self, chi2):
        """Abinit χ^(2) stores d[i,j,k]; last two indices must be symmetric."""
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    assert abs(chi2[i, j, k] - chi2[i, k, j]) < 1e-8, \
                        f"Symmetry broken at [{i},{j},{k}] vs [{i},{k},{j}]"

    def test_d_to_chi2_factor(self, chi2):
        """χ^(2) = 2d before the reader applies point-group symmetrisation."""
        from PDielec.AbinitOutputReader import AbinitOutputReader
        r = AbinitOutputReader([_ABINIT_FILE])
        r.read_output()
        factor = chi2_pm_per_v_to_repsilon(r.volume)
        d_recovered = chi2 / (2.0 * factor)
        assert d_recovered[2, 2, 2] == pytest.approx(-42.296970629, abs=1e-6)

    def test_zno_6mm_symmetry_after_reader_symmetrisation(self, chi2):
        """ZnO should obey 6mm symmetry even if Abinit's raw d table does not."""
        assert chi2[0, 0, 2] == pytest.approx(chi2[1, 1, 2], abs=2e-8)
        assert chi2[0, 1, 2] == pytest.approx(0.0, abs=1e-8)
        assert chi2[1, 0, 2] == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------------------
# J9: QE reader parses χ^(2) from electro-optic tensor block
# ---------------------------------------------------------------------------

class TestJ9QEParser:
    """J9: QE _read_nlo_susceptibility reads chi^2 from electro-optic tensors."""

    def test_qe_electro_optic_block(self, tmp_path):
        """Parse QE electro-optic matrices and convert them to internal chi^2 units."""
        from PDielec.QEOutputReader import QEOutputReader

        output = tmp_path / "qe.raman.log"
        output.write_text(
            """
          Electro-optic tensor in cartesian axis:

          (       1.000000000       2.000000000       3.000000000 )
          (       4.000000000       5.000000000       6.000000000 )
          (       7.000000000       8.000000000       9.000000000 )

          (      10.000000000      11.000000000      12.000000000 )
          (      13.000000000      14.000000000      15.000000000 )
          (      16.000000000      17.000000000      18.000000000 )

          (      19.000000000      20.000000000      21.000000000 )
          (      22.000000000      23.000000000      24.000000000 )
          (      25.000000000      26.000000000      27.000000000 )
            """,
            encoding="utf-8",
        )
        r = QEOutputReader([str(output)])
        r.volume = 100.0
        r.read_output()
        chi2 = r.nonlinear_optical_susceptibility

        assert chi2 is not None
        assert chi2.shape == (3, 3, 3)
        factor = 0.5 * 2.7502 * chi2_pm_per_v_to_repsilon(r.volume)
        assert chi2[0, 0, 0] == pytest.approx(1.0 * factor)
        assert chi2[2, 2, 2] == pytest.approx(27.0 * factor)
        assert chi2[1, 2, 1] == pytest.approx(15.0 * factor)

    @pytest.mark.skipif(not _have_qe_tensors, reason="QE tensors.xml not present")
    def test_qe_elop_tns_xml_block(self):
        """Parse QE ELOP_TNS from tensors.xml and convert it to internal chi^2 units."""
        from PDielec.QEOutputReader import QEOutputReader

        r = QEOutputReader([_QE_TENSORS_FILE])
        # tensors.xml contains no cell geometry; supply a test volume before
        # converting the physical pm/V values to the internal convention.
        r.volume = 100.0
        r.read_output()
        chi2 = r.nonlinear_optical_susceptibility

        assert chi2 is not None
        assert chi2.shape == (3, 3, 3)
        factor = 0.5 * 2.7502 * chi2_pm_per_v_to_repsilon(r.volume)
        assert chi2[0, 0, 2] == pytest.approx(133.8996252820044 * factor)
        assert chi2[2, 2, 2] == pytest.approx(-125.0252942735832 * factor)

    @pytest.mark.skipif(not _have_qe_tensors, reason="QE tensors.xml not present")
    def test_qe_elop_tns_xml_preferred_over_log(self, tmp_path):
        """Keep the higher-precision XML electro-optic tensor when a log block is also present."""
        from PDielec.QEOutputReader import QEOutputReader

        output = tmp_path / "qe.raman.log"
        output.write_text(
            """
          Electro-optic tensor in cartesian axis:

          (       1.000000000       1.000000000       1.000000000 )
          (       1.000000000       1.000000000       1.000000000 )
          (       1.000000000       1.000000000       1.000000000 )

          (       1.000000000       1.000000000       1.000000000 )
          (       1.000000000       1.000000000       1.000000000 )
          (       1.000000000       1.000000000       1.000000000 )

          (       1.000000000       1.000000000       1.000000000 )
          (       1.000000000       1.000000000       1.000000000 )
          (       1.000000000       1.000000000       1.000000000 )
            """,
            encoding="utf-8",
        )
        r = QEOutputReader([_QE_TENSORS_FILE, str(output)])
        r.volume = 100.0
        r.read_output()

        factor = 0.5 * 2.7502 * chi2_pm_per_v_to_repsilon(r.volume)
        assert r.nonlinear_optical_susceptibility[0, 0, 2] == pytest.approx(133.8996252820044 * factor)


# ---------------------------------------------------------------------------
# J3: CASTEP reader parses χ^(2) from raman.castep
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_castep, reason="raman.castep not present")
class TestJ3CastepParser:
    """J3: CASTEP _read_nlo_susceptibility reads χ^(2) = 2d (Voigt) correctly."""

    @pytest.fixture(scope="class")
    def chi2(self):
        from PDielec.CastepOutputReader import CastepOutputReader
        r = CastepOutputReader([_CASTEP_FILE])
        r.read_output()
        return r.nonlinear_optical_susceptibility

    def test_tensor_not_none(self, chi2):
        assert chi2 is not None

    def test_shape(self, chi2):
        assert chi2.shape == (3, 3, 3)

    def test_zz_component(self, chi2):
        # raman.castep: d[2,2,2] = -38.3539; χ^(2) = 2d ≈ -76.71 pm/V before internal conversion.
        from PDielec.CastepOutputReader import CastepOutputReader
        reader = CastepOutputReader([_CASTEP_FILE])
        reader.read_output()
        factor = chi2_pm_per_v_to_repsilon(reader.volume)
        assert abs(chi2[2, 2, 2] - (-76.7078 * factor)) < 1e-4

    def test_last_two_indices_symmetric(self, chi2):
        """Voigt reconstruction must produce a tensor symmetric in last two indices."""
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    assert abs(chi2[i, j, k] - chi2[i, k, j]) < 1e-8, \
                        f"Symmetry broken at [{i},{j},{k}] vs [{i},{k},{j}]"

    def test_d_to_chi2_factor(self, chi2):
        """χ^(2) = 2d; row i=0 col 4 is Voigt pair (0,2): d[0,0,2] = 21.655 → chi2 = 43.31."""
        # raman.castep row 0: ... 21.65509 at col 4 → d[0,0,2] = d[0,2,0] = 21.655
        from PDielec.CastepOutputReader import CastepOutputReader
        reader = CastepOutputReader([_CASTEP_FILE])
        reader.read_output()
        factor = chi2_pm_per_v_to_repsilon(reader.volume)
        # The reader subsequently averages symmetry-equivalent components.
        assert chi2[0, 0, 2] / factor == pytest.approx(2.0 * 21.65509, abs=4e-3)


@pytest.mark.skipif(not _have_crystal23, reason="CRYSTAL23 Raman example not present")
class TestJ3BCrystalParser:
    """J3B: CRYSTAL CHI2.DAT reads χ^(2) = 2d correctly."""

    def test_crystal_chi2_dat_internal_units(self):
        """CRYSTAL d(MKS) values are converted from pm/V to internal units."""
        from PDielec.CrystalOutputReader import CrystalOutputReader

        r = CrystalOutputReader([_CRYSTAL23_FILE])
        r.read_output()
        chi2 = r.nonlinear_optical_susceptibility
        factor = chi2_pm_per_v_to_repsilon(r.volume)

        assert chi2 is not None
        assert chi2.shape == (3, 3, 3)
        assert chi2[0, 0, 2] == pytest.approx(2.0 * 11.691 * factor)
        assert chi2[1, 1, 2] == pytest.approx(2.0 * 11.691 * factor)
        assert chi2[2, 2, 2] == pytest.approx(2.0 * -33.166 * factor)

    def test_crystal_chi2_from_output_without_companion(self, tmp_path):
        """The embedded CRYSTAL output block is sufficient when CHI2.DAT is absent."""
        output_file = tmp_path / "opt_raman.out"
        shutil.copyfile(_CRYSTAL23_FILE, output_file)

        from PDielec.CrystalOutputReader import CrystalOutputReader

        r = CrystalOutputReader([str(output_file)])
        r.read_output()
        chi2 = r.nonlinear_optical_susceptibility
        factor = chi2_pm_per_v_to_repsilon(r.volume)

        assert chi2 is not None
        assert chi2.shape == (3, 3, 3)
        assert chi2[0, 0, 2] == pytest.approx(2.0 * 9.7410 * factor)
        assert chi2[1, 1, 2] == pytest.approx(2.0 * 9.7410 * factor)
        assert chi2[2, 2, 2] == pytest.approx(2.0 * -25.982 * factor)


# ---------------------------------------------------------------------------
# J4: Backward compatibility — chi2=None leaves results unchanged
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_abinit, reason="raman.abo not present")
class TestJ4BackwardCompat:
    """J4: chi2_repsilon=None must give identical results to the old call."""

    def test_frequencies_unchanged(self):
        r = _load_abinit()
        args = _standalone_inputs(r)
        hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas = args
        q = np.array([0.0, 0.0, 1.0])

        freqs_base, tensors_base, sigmas_base = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas)
        freqs_none, tensors_none, sigmas_none = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=None)

        np.testing.assert_array_equal(freqs_base, freqs_none)
        for R1, R2 in zip(tensors_base, tensors_none):
            np.testing.assert_array_equal(R1, R2)

    def test_sigmas_unchanged(self):
        r = _load_abinit()
        args = _standalone_inputs(r)
        hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas = args
        q = np.array([1.0, 0.0, 0.0])

        _, _, sigmas_base = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas)
        _, _, sigmas_none = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=None)

        np.testing.assert_array_equal(sigmas_base, sigmas_none)


# ---------------------------------------------------------------------------
# J5: Zero χ^(2) → EO correction is exactly zero
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_abinit, reason="raman.abo not present")
class TestJ5ZeroChi2:
    """J5: A zero χ^(2) tensor must leave all Raman tensors unchanged."""

    @pytest.mark.parametrize("q", [
        np.array([0.0, 0.0, 1.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([1.0, 1.0, 1.0]) / math.sqrt(3),
    ])
    def test_zero_chi2_no_change(self, q):
        r = _load_abinit()
        args = _standalone_inputs(r)
        hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas = args
        chi2_zero = np.zeros((3, 3, 3))

        _, tensors_no, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=None)
        _, tensors_z, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=chi2_zero)

        for R_no, R_z in zip(tensors_no, tensors_z):
            np.testing.assert_allclose(R_z, R_no, atol=1e-14)


# ---------------------------------------------------------------------------
# J6: Non-zero χ^(2) → at least one polar mode tensor changes
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_abinit, reason="raman.abo not present")
class TestJ6NonZeroChi2:
    """J6: Real χ^(2) from Abinit must change at least one NAC Raman tensor."""

    @pytest.mark.parametrize("q", [
        np.array([0.0, 0.0, 1.0]),
        np.array([1.0, 0.0, 0.0]),
    ])
    def test_at_least_one_tensor_changes(self, q):
        r = _load_abinit()
        args = _standalone_inputs(r)
        hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas = args
        chi2 = r.nonlinear_optical_susceptibility

        _, tensors_no, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=None)
        _, tensors_eo, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=chi2)

        diffs = [np.max(np.abs(R_eo - R_no))
                 for R_eo, R_no in zip(tensors_eo, tensors_no)]
        assert max(diffs) > 1e-10, "Expected at least one tensor to change with non-zero χ^(2)"


# ---------------------------------------------------------------------------
# J7: EO correction scales linearly with χ^(2) magnitude
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_abinit, reason="raman.abo not present")
class TestJ7LinearScaling:
    """J7: Doubling χ^(2) must double the EO correction ΔR for every mode."""

    @pytest.mark.parametrize("q", [
        np.array([0.0, 0.0, 1.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([1.0, 1.0, 0.0]) / math.sqrt(2),
    ])
    def test_linear_scaling(self, q):
        r = _load_abinit()
        args = _standalone_inputs(r)
        hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas = args
        chi2 = r.nonlinear_optical_susceptibility

        _, tensors_no, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=None)
        _, tensors_1x, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=chi2)
        _, tensors_2x, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=2.0 * chi2)

        for p, (R_no, R_1x, R_2x) in enumerate(zip(tensors_no, tensors_1x, tensors_2x)):
            delta_1x = R_1x - R_no
            delta_2x = R_2x - R_no
            # delta_2x should equal 2 * delta_1x
            np.testing.assert_allclose(delta_2x, 2.0 * delta_1x, atol=1e-12,
                                       err_msg=f"Linear scaling failed for mode {p}")


# ---------------------------------------------------------------------------
# J8: Zero Born charges → EO correction vanishes
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_abinit, reason="raman.abo not present")
class TestJ8ZeroBornCharges:
    """J8: When Born charges are zero, Z_q = 0 and the EO correction must be zero."""

    def test_zero_born_charges_no_correction(self):
        r = _load_abinit()
        args = _standalone_inputs(r)
        hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas = args
        chi2 = r.nonlinear_optical_susceptibility

        bc_zero = np.zeros_like(bc)
        q = np.array([0.0, 0.0, 1.0])

        _, tensors_no, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc_zero, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=None)
        _, tensors_eo, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc_zero, eps_inf, vol, masses, U_TO, tensors, sigmas, chi2_repsilon=chi2)

        for p, (R_no, R_eo) in enumerate(zip(tensors_no, tensors_eo)):
            np.testing.assert_allclose(R_eo, R_no, atol=1e-14,
                                       err_msg=f"EO correction non-zero for mode {p} with zero Born charges")


# ---------------------------------------------------------------------------
# J10: Reader tensor convention
# ---------------------------------------------------------------------------

class TestJ10ReaderTensorConvention:
    """J10: Reader tensors are physical R_epsilon tensors for EO correction."""

    def test_eo_correction_uses_angstrom_r_epsilon_units(self):
        """The EO tensor increment must be in the same units as reader R_epsilon tensors."""
        from PDielec.RamanPolarCalculator import apply_eo_correction

        tensors = [np.zeros((3, 3))]
        chi2 = np.zeros((3, 3, 3))
        chi2[0, 0, 2] = 4.0
        q_hat = np.array([0.0, 0.0, 1.0])
        Z_mat = np.array([[0.0], [0.0], [2.0]])
        eigvecs = np.array([[1.0]])
        eps_inf = np.diag([3.0, 3.0, 5.0])

        corrected = apply_eo_correction(tensors, chi2, q_hat, Z_mat, eigvecs, eps_inf)

        expected = -8.0 * np.pi * chi2[0, 0, 2] * 2.0 / eps_inf[2, 2]
        assert corrected[0][0, 0] == pytest.approx(expected)

    def test_eo_correction_operates_directly_on_reader_tensors(self):
        from PDielec.RamanPolarCalculator import apply_eo_correction

        tensors = [
            np.diag([1.0, 2.0, 3.0]),
            np.array([[0.1, 0.2, 0.0], [0.2, 0.3, 0.4], [0.0, 0.4, 0.5]]),
        ]
        chi2 = np.zeros((3, 3, 3))
        chi2[0, 0, 2] = 4.0
        chi2[1, 1, 2] = -2.0
        q_hat = np.array([0.0, 0.0, 1.0])
        Z_mat = np.array([
            [1.0, 0.0],
            [0.0, 0.5],
            [0.25, -0.75],
        ])
        eigvecs = np.eye(2)
        eps_inf = np.diag([3.0, 3.0, 4.0])

        corrected = apply_eo_correction(tensors, chi2, q_hat, Z_mat, eigvecs, eps_inf)
        delta = [corrected_tensor - tensor for corrected_tensor, tensor in zip(corrected, tensors)]

        for scale in (2.0, 5.0):
            scaled_corrected = apply_eo_correction(
                [scale * tensor for tensor in tensors], chi2, q_hat, Z_mat, eigvecs, eps_inf)
            for scaled_tensor, tensor, delta_tensor in zip(scaled_corrected, tensors, delta):
                np.testing.assert_allclose(scaled_tensor, scale * tensor + delta_tensor, atol=1e-14)


# ---------------------------------------------------------------------------
# J11: AbInit Raman tensors stay consistent after GUI mass/mode recalculation
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_powder_abinit, reason="Powder Raman AbInit raman.abo not present")
class TestJ11AbinitRamanTensorRecalculation:
    """J11: Recalculated AbInit normal modes must keep Raman tensors phase-consistent."""

    def test_zno_a1_lo_eo_correction_keeps_its_phase_after_mass_recalculation(self):
        """The AbInit ZnO q||c A1 LO correction should keep its phase after Settings setup."""
        from PDielec.AbinitOutputReader import AbinitOutputReader
        from PDielec.HelperRoutines import calculate_dft_permittivity_object

        r = AbinitOutputReader([_POWDER_ABINIT_FILE])
        r.read_output()

        # Mirrors the Settings/DFT-permittivity setup used by PDGui.  This
        # recalculates mass-weighted normal modes for the chosen mass convention.
        calculate_dft_permittivity_object(r, sigma=5.0, eckart=True, mass_definition="average")

        hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas = _standalone_inputs(r)
        q = np.array([0.0, 0.0, 1.0])
        freqs_no, tensors_no, _, dominant_to = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas,
            chi2_repsilon=None,
            return_mode_map=True,
        )
        _freqs_eo, tensors_eo, _ = _compute_nac_dynamical_matrix_standalone(
            q, hessian, bc, eps_inf, vol, masses, U_TO, tensors, sigmas,
            chi2_repsilon=r.nonlinear_optical_susceptibility,
        )

        a1_lo = [
            idx for idx, (freq, to_idx) in enumerate(zip(freqs_no, dominant_to))
            if 550.0 < freq < 570.0 and to_idx == 6
        ]
        assert len(a1_lo) == 1
        mode_idx = a1_lo[0]

        assert tensors_eo[mode_idx][0, 0] - tensors_no[mode_idx][0, 0] < 0.0
        assert tensors_eo[mode_idx][1, 1] - tensors_no[mode_idx][1, 1] < 0.0


# ---------------------------------------------------------------------------
# J12: Generic per-mode Raman tensor basis transform
# ---------------------------------------------------------------------------

class TestJ12GenericRamanTensorBasisTransform:
    """J12: Generic readers keep per-mode Raman tensors aligned after mode recalculation."""

    def test_generic_transform_corrects_single_mode_phase_flip(self):
        old_modes = np.array([
            [[1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        ])
        new_modes = np.array([
            [[-1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        ])
        old_tensors = [np.eye(3), 2.0 * np.eye(3)]

        transformed = GenericOutputReader._transform_raman_tensors_between_mode_bases(
            old_modes, new_modes, old_tensors)

        np.testing.assert_allclose(transformed[0], -np.eye(3))
        np.testing.assert_allclose(transformed[1], 2.0 * np.eye(3))

    def test_generic_transform_rotates_degenerate_subspace(self):
        angle = np.pi / 4.0
        old_modes = np.array([
            [[1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        ])
        new_modes = np.array([
            [[np.cos(angle), np.sin(angle), 0.0]],
            [[-np.sin(angle), np.cos(angle), 0.0]],
        ])
        R_x = np.diag([1.0, 0.0, 0.0])
        R_y = np.diag([0.0, 2.0, 0.0])

        transformed = GenericOutputReader._transform_raman_tensors_between_mode_bases(
            old_modes, new_modes, [R_x, R_y])

        np.testing.assert_allclose(transformed[0], np.cos(angle) * R_x + np.sin(angle) * R_y)
        np.testing.assert_allclose(transformed[1], -np.sin(angle) * R_x + np.cos(angle) * R_y)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
