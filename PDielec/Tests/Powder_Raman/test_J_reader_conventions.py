"""Reader convention checks for powder Raman examples."""

import math
import os

import numpy as np
import pytest
import yaml

from PDielec import Calculator
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.VaspOutputReader import VaspOutputReader


_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
_CASTEP_FILE = os.path.join(_REPO, "Examples", "Powder_Raman", "Castep", "raman.castep")
_CRYSTAL_FILE = os.path.join(_REPO, "Examples", "Powder_Raman", "Crystal23", "opt_raman.out")
_QE_DIR = os.path.join(_REPO, "Examples", "Powder_Raman", "QE")
_QE_LOG = os.path.join(_QE_DIR, "zno.raman.log")
_QE_DYNG = os.path.join(_QE_DIR, "zno.raman.dynG")
_QE_TENSORS = os.path.join(_QE_DIR, "tensors.xml")
_VASP_OUTCAR = os.path.join(_REPO, "Examples", "Powder_Raman", "Vasp", "OUTCAR")
_VASP_RAMAN_TENSORS = os.path.join(_REPO, "Examples", "Powder_Raman", "Vasp", "Raman-Tensors.yaml")


@pytest.mark.skipif(not os.path.exists(_CASTEP_FILE), reason="CASTEP powder Raman example not present")
def test_castep_raman_tensors_convert_back_to_castep_a4_per_amu_activities():
    """CASTEP printed tensors should reproduce CASTEP A^4/amu activities in display units."""
    reader = CastepOutputReader([_CASTEP_FILE])
    reader.read_output()

    epsilon_activities = Calculator.raman_intensities(reader.get_raman_tensors(), reader.volume)
    polarizability_display_factor = reader.volume / (16.0 * math.pi * math.pi)
    activities_a4_per_amu = epsilon_activities[:, 0] * polarizability_display_factor

    # Reference values from Examples/Powder_Raman/Castep/raman.phonon.
    assert activities_a4_per_amu[3] == pytest.approx(0.4757574, abs=2.0e-2)
    assert activities_a4_per_amu[6] == pytest.approx(13.6090575, abs=2.0e-2)
    assert activities_a4_per_amu[9] == pytest.approx(2.0180457, abs=2.0e-2)


@pytest.mark.skipif(not os.path.exists(_CRYSTAL_FILE), reason="CRYSTAL powder Raman example not present")
def test_crystal_raman_tensors_convert_to_a4_per_amu_display_activities():
    """CRYSTAL TENS_RAMAN.DAT tensors should use the shared R_epsilon convention."""
    reader = CrystalOutputReader([_CRYSTAL_FILE])
    reader.read_output()

    epsilon_activities = Calculator.raman_intensities(reader.get_raman_tensors(), reader.volume)
    polarizability_display_factor = reader.volume / (16.0 * math.pi * math.pi)
    activities_a4_per_amu = epsilon_activities[:, 0] * polarizability_display_factor

    # Reference values are the absolute PDielec polarizability-volume activities
    # obtained from the CRYSTAL23 TENS_RAMAN.DAT projected polarizability tensors.
    assert activities_a4_per_amu[6] == pytest.approx(2.72482, abs=1.0e-4)
    assert activities_a4_per_amu[7] == pytest.approx(0.0428894, abs=1.0e-6)
    assert activities_a4_per_amu[9] == pytest.approx(1.73748, abs=1.0e-4)


@pytest.mark.skipif(
    not (os.path.exists(_QE_LOG) and os.path.exists(_QE_DYNG) and os.path.exists(_QE_TENSORS)),
    reason="QE powder Raman example not present",
)
def test_qe_log_a2_and_tensors_xml_raman_paths_agree_for_matching_dyng():
    """QE log A^2 and tensors.xml dε/du paths should produce the same R_epsilon activities."""
    log_reader = QEOutputReader([_QE_LOG, _QE_DYNG])
    log_reader.read_output()
    xml_reader = QEOutputReader([_QE_DYNG, _QE_TENSORS])
    xml_reader.read_output()

    log_activities = Calculator.raman_intensities(log_reader.get_raman_tensors(), log_reader.volume)
    xml_activities = Calculator.raman_intensities(xml_reader.get_raman_tensors(), xml_reader.volume)
    display_factor = log_reader.volume / (16.0 * math.pi * math.pi)

    assert log_reader.volume == pytest.approx(xml_reader.volume, abs=1.0e-10)
    assert log_activities[8, 0] * display_factor == pytest.approx(xml_activities[8, 0] * display_factor, abs=1.0e-8)
    assert log_activities[9, 0] * display_factor == pytest.approx(xml_activities[9, 0] * display_factor, abs=1.0e-8)
    assert xml_activities[8, 0] * display_factor == pytest.approx(215.944797, abs=1.0e-5)


@pytest.mark.skipif(
    not (os.path.exists(_QE_LOG) and os.path.exists(_QE_DYNG) and os.path.exists(_QE_TENSORS)),
    reason="QE powder Raman example not present",
)
def test_qe_tensors_xml_raman_preferred_over_lower_precision_log_blocks():
    """QE tensors.xml Raman data should remain active even if a log/dynG block is read later."""
    reader = QEOutputReader([_QE_TENSORS, _QE_DYNG, _QE_LOG])
    reader.read_output()

    activities = Calculator.raman_intensities(reader.get_raman_tensors(), reader.volume)
    display_factor = reader.volume / (16.0 * math.pi * math.pi)

    assert reader._qe_raman_suscept_units == "xml_epsilon_per_bohr"
    assert activities[8, 0] * display_factor == pytest.approx(215.944797, abs=1.0e-5)


@pytest.mark.skipif(
    not (os.path.exists(_VASP_OUTCAR) and os.path.exists(_VASP_RAMAN_TENSORS)),
    reason="VASP powder Raman example not present",
)
def test_vasp_skelton_raman_tensors_convert_to_a4_per_amu_display_activities():
    """Skelton VASP tensors are polarizability-volume tensors in Å^2/sqrt(amu)."""
    reader = VaspOutputReader([_VASP_OUTCAR])
    reader.read_output()

    activities = Calculator.raman_intensities(reader.get_raman_tensors(), reader.volume)
    display_factor = reader.volume / (16.0 * math.pi * math.pi)
    activities_a4_per_amu = activities[:, 0] * display_factor

    with open(_VASP_RAMAN_TENSORS) as fd:
        skelton_data = yaml.safe_load(fd)
    raw_entries = sorted(skelton_data["raman_activities"], key=lambda entry: entry["band_index"], reverse=True)
    raw_a4_per_amu = [
        Calculator.raman_intensities([np.array(entry["raman_tensor"], dtype=float)], reader.volume)[0, 0]
        for entry in raw_entries
    ]

    # Reader display values should recover the raw Skelton polarizability-volume
    # invariants in PDielec ordering: acoustic/imaginary modes first, then
    # optical modes in increasing frequency.
    assert reader.frequencies[:3] == pytest.approx([-2.000636, -2.000636, -1.490162], abs=1.0e-6)
    assert activities_a4_per_amu[:3] == pytest.approx([0.0, 0.0, 0.0], abs=1.0e-12)
    assert reader.frequencies[3] == pytest.approx(89.892444, abs=1.0e-6)
    assert activities_a4_per_amu[3] == pytest.approx(raw_a4_per_amu[0], abs=1.0e-10)
    assert activities_a4_per_amu[3] == pytest.approx(0.3389968, abs=1.0e-6)
    assert activities_a4_per_amu[4] == pytest.approx(0.5950713, abs=1.0e-6)
    assert activities_a4_per_amu[6] == pytest.approx(75.238672, abs=1.0e-6)
    assert activities_a4_per_amu[7] == pytest.approx(3.8560467, abs=1.0e-6)
    assert activities_a4_per_amu[8] == pytest.approx(3.8425679, abs=1.0e-6)
