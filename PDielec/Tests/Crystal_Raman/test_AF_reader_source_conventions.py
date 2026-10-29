"""Source-format and constitutive checks for nonlinear-response readers.

QE 7.5: PHonon/PH/el_opt.f90 writes (elop_(ipa,ipb,ipc),ipa=1,3)
inside loops over ipc then ipb. upflib/xmltools.f90 writes rank-three arrays
in Fortran element order. Distinct components expose otherwise hidden swaps.
"""

import io
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.QEOutputReader import QEOutputReader


@pytest.mark.parametrize("source", ["log", "xml"])
def test_qe_tensor_indices_follow_upstream_writer(source):
    """Preserve ipa, ipb, ipc rather than transposing the first two indices."""
    tensor = np.fromfunction(lambda i, j, k: 100.0 * (i + 1) + 10.0 * (j + 1) + k + 1, (3, 3, 3))
    reader = QEOutputReader([])
    if source == "log":
        # Reproduce el_opt's explicit ipc/ipb/ipa loops, not a NumPy reshape.
        lines = ["( " + " ".join(str(tensor[i, j, k]) for i in range(3)) + " )"
                 for k in range(3) for j in range(3)]
        reader.file_descriptor = io.StringIO("\n".join(lines))
        assert reader._read_nlo_susceptibility("")
    else:
        element = ET.Element("ELOP_TNS")
        element.text = " ".join(str(value) for value in tensor.ravel(order="F"))
        assert reader._read_nlo_susceptibility_xml(element)
    physical = reader._nonlinear_optical_susceptibility_pm_per_v
    np.testing.assert_allclose(physical, 0.5 * 2.7502 * tensor, rtol=1.0e-14)


@pytest.mark.parametrize("column", range(6))
def test_castep_contracted_d_reconstructs_polarization(column):
    """Full chi(2) reproduces the contracted d-vector including cross-field multiplicity."""
    d_matrix = np.zeros((3, 6))
    d_matrix[:, column] = [1.0, -2.0, 3.0]
    block = "separator\n" + "\n".join(" ".join(map(str, row)) for row in d_matrix)
    reader = CastepOutputReader(["synthetic.castep"])
    reader.file_descriptor = io.StringIO(block)
    assert reader._read_nlo_susceptibility("")
    field = np.array([0.7, -1.3, 2.1])
    # P/epsilon_0 = 2 d (Ex^2,Ey^2,Ez^2,2EyEz,2ExEz,2ExEy).
    contracted_fields = np.array([field[0]**2, field[1]**2, field[2]**2,
                                  2 * field[1] * field[2], 2 * field[0] * field[2],
                                  2 * field[0] * field[1]])
    expected = 2.0 * d_matrix @ contracted_fields
    actual = np.einsum("ijk,j,k->i", reader._nonlinear_optical_susceptibility_pm_per_v, field, field)
    np.testing.assert_allclose(actual, expected, atol=1.0e-14)


def test_crystal_native_chi2_and_d_columns_agree():
    """Check the native table's independent CHI(2) column and printed unit identity."""
    path = Path(__file__).resolve().parents[3] / "Examples/Powder_Raman/Crystal23/CHI2.DAT"
    reader = CrystalOutputReader([])
    with path.open() as handle:
        reader._read_chi2_block(handle, str(path))
    # CHI2.DAT explicitly defines d(AU)=CHI(2)/2 and d(MKS)=d(AU)/0.514220632.
    # Read CHI(2), not the d(MKS) column used by the parser.
    for line in path.read_text().splitlines():
        if ")" not in line or line.split()[0] not in {"XXZ", "YYZ", "ZZZ"}:
            continue
        label = line.split()[0]
        chi_au = float(line.split(")", 1)[1].split()[0])
        indices = tuple("XYZ".index(axis) for axis in label)
        expected = chi_au / 0.514220632
        assert reader._nonlinear_optical_susceptibility_pm_per_v[indices] == pytest.approx(expected, rel=6.0e-5)
