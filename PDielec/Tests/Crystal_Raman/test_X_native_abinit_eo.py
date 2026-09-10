"""Independent native ANADDB directional EO reference, with no fitted amplitude."""

import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.Constants import amu, angs2bohr
from PDielec.NACDiagnostics import compute_nac_direction_diagnostics

DATA = Path(__file__).parent / "data/abinit_zno_native"


def native_tensors(text):
    """Extract the native transverse and directional longitudinal printed tensors."""
    start = text.index("Raman susceptibilities of transverse zone-center phonon modes")
    marker = "Raman susceptibility of zone-center phonons, with non-analyticity in the"
    sections = [text[start : text.index("Output of the EO tensor", start)]] + text.split(marker)[1:]
    result = []
    for section in sections:
        modes = {}
        lines = section.splitlines()
        for i, line in enumerate(lines):
            match = re.match(r"\s*Mod(?:e)?\s+(\d+)\s+\(", line)
            if match:
                modes[int(match[1])] = np.array([[float(x) for x in row.split()[1:4]] for row in lines[i + 1 : i + 4]])
        result.append(modes)
    return result


def test_native_anaddb_signed_lo_tensors(tmp_path):
    """Absolute tensor normalization, EO sign and magnitude agree with native ANADDB."""
    provenance = json.loads((DATA / "provenance.json").read_text())
    raw = {}
    for name, metadata in provenance["files"].items():
        raw[name] = gzip.decompress((DATA / (name + ".gz")).read_bytes())
        assert hashlib.sha256(raw[name]).hexdigest() == metadata["sha256_uncompressed"]
    path = tmp_path / "raman.abo"
    path.write_bytes(raw["raman.abo"])
    reader = AbinitOutputReader([str(path)])
    reader.read_output()
    reader.calculate_mass_weighted_normal_modes()
    to, x, z = native_tensors(raw["anaddb_lo.abo"].decode())
    # Native sqrt(V_bohr3) dchi_G/dQ_electron_mass -> reader R_epsilon.
    scale = 4 * np.pi * np.sqrt(amu / angs2bohr)
    for q, native, to_index, tolerance in [([1.0, 0, 0], x, 8, 2e-4), ([0.0, 0, 1], z, 7, 6e-5)]:
        d = compute_nac_direction_diagnostics(
            np.array(q),
            reader.hessian,
            reader.born_charges,
            reader.zerof_optical_dielectric,
            reader.volume * angs2bohr**3,
            np.array(reader.masses) * amu,
            np.array(reader.mass_weighted_normal_modes).reshape(12, 12),
            reader.raman_tensors,
            np.ones(12),
            chi2_repsilon=reader.nonlinear_optical_susceptibility,
        )
        mechanical = np.asarray(d["raman_tensors_without_eo"][10])
        corrected = np.asarray(d["raman_tensors_with_eo"][10])
        sign = np.sign(np.vdot(to[to_index], mechanical).real)
        assert np.linalg.norm(mechanical - sign * scale * to[to_index]) / np.linalg.norm(mechanical) < 1e-6
        expected = sign * scale * native[11]
        assert np.linalg.norm(corrected - expected) / np.linalg.norm(expected) < tolerance
