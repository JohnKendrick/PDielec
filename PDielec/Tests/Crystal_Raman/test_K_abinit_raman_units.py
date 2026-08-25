"""Regression test for ABINIT susceptibility-to-dielectric Raman units."""
from pathlib import Path
import math

import numpy as np
import pytest

from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.Constants import angs2bohr


ABINIT_OUTPUT = (
    Path(__file__).resolve().parents[3]
    / "Examples/Powder_Raman/AbInit/raman.abo"
)


@pytest.mark.skipif(not ABINIT_OUTPUT.is_file(), reason="ABINIT Raman example unavailable")
def test_abinit_raman_tensor_converts_dchi_to_depsilon():
    """Reader tensors use R_epsilon=sqrt(V) d epsilon/dQ, not d chi/dQ."""
    reader = AbinitOutputReader([str(ABINIT_OUTPUT)])
    reader.read_output()

    dchi = reader._susceptibility_derivatives
    assert dchi is not None
    assert reader.mass_weighted_normal_modes

    mode_index = 6  # ZnO A1 transverse mode in the bundled example.
    expected = np.zeros((3, 3))
    for atom in range(reader.nions):
        for direction in range(3):
            expected += (
                dchi[atom, direction]
                * reader.mass_weighted_normal_modes[mode_index][atom][direction]
                / math.sqrt(reader.masses[atom])
            )
    expected *= math.sqrt(reader.volume) * angs2bohr * 4.0 * math.pi

    np.testing.assert_allclose(reader.raman_tensors[mode_index], expected, rtol=1e-12, atol=1e-12)
