"""Generate and validate documentation inputs from the existing CASTEP ZnO data."""
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from PDielec.Utilities import get_reader

OUT = ROOT / 'Sphinx/_static/inputs'
OUT.mkdir(parents=True, exist_ok=True)
r = get_reader(str(ROOT / 'Examples/Crystal_Raman/Castep/raman.castep'), 'castep')
r.read_output()
u = np.asarray(r.calculate_mass_weighted_normal_modes())
cell = r.get_unit_cell()
raman = np.asarray(r.raman_tensors).real
masses = np.asarray(r.masses)
derivatives = np.einsum('mak,mij->akij', u, raman) * np.sqrt(masses[:, None, None, None] / r.volume)
data = {
    'schema': 'pdielect-finite-field-raman-v1',
    'provenance': 'Converted from Examples/Crystal_Raman/Castep/raman.castep by Sphinx/tools/make_format_examples.py; not a new finite-difference DFT calculation.',
    'structure': {'lattice_angstrom': np.asarray(cell.lattice).tolist(), 'symbols': cell.element_names, 'positions_fractional': np.asarray(cell.fractional_coordinates).tolist()},
    'program_masses_amu': masses.tolist(),
    'epsilon_infinity': r.zerof_optical_dielectric,
    'born_effective_charges': r.born_charges,
    'chi2': np.asarray(r.nonlinear_optical_susceptibility).tolist(),
    'deps_dr': derivatives.tolist(),
    'modes': [{'frequency_cm-1': float(f), 'mass_weighted_eigenvector': v.tolist(), 'raman_tensor': t.tolist()} for f, v, t in zip(r.frequencies, u, raman)],
}
(OUT / 'zno-raman.json').write_text(json.dumps(data, indent=2) + '\n')
# A complete mode set, with explicit epsilon-convention tensors.
lines = ['# ZnO converted from the existing CASTEP example; no EO block in .exp.', 'lattice', '1.0']
lines += [' '.join(map(str, row)) for row in cell.lattice]
lines += ['species ' + str(r.nspecies)]
lines += [f'{s} {m}' for s, m in zip(r.species, r.masses_per_type)]
lines += ['unitcell ' + str(r.nions)]
lines += [s + ' ' + ' '.join(map(str, xyz)) for s, xyz in zip(cell.element_names, cell.fractional_coordinates)]
lines += ['epsinf'] + [' '.join(map(str, row)) for row in r.zerof_optical_dielectric]
lines += ['frequencies ' + str(len(r.frequencies))] + [str(f) for f in r.frequencies]
lines += ['born_charges ' + str(r.nions)]
lines += [' '.join(map(str, row)) for tensor in r.born_charges for row in tensor]
lines += ['normal_modes ' + str(len(u))]
lines += [' '.join(map(str, row)) for mode in u for row in mode]
lines += ['raman_tensors ' + str(len(raman)) + ' epsilon']
lines += [' '.join(map(str, row)) for tensor in raman for row in tensor]
(OUT / 'zno-raman.exp').write_text('\n'.join(lines) + '\n')
(OUT / 'constant.exp').write_text('''# Illustrative constant dielectric, not measured data.
lattice
1.0
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
species 1
Si 28.0855
unitcell 1
Si 0.0 0.0 0.0
constant
2.0+0.1j 0.0 0.0
0.0 2.0+0.1j 0.0
0.0 0.0 2.0+0.1j
''')
minimal = {'schema': 'pdielect-finite-field-raman-v1', 'structure': {'lattice_angstrom': [[5.,0.,0.],[0.,5.,0.],[0.,0.,5.]], 'symbols': ['Si'], 'positions_fractional': [[0.,0.,0.]]}, 'epsilon_infinity': [[2.,0.,0.],[0.,2.,0.],[0.,0.,2.]]}
(OUT / 'minimal.json').write_text(json.dumps(minimal, indent=2) + '\n')
for filename, program in [('constant.exp', 'experiment'), ('zno-raman.exp', 'experiment'), ('minimal.json', 'finite_field'), ('zno-raman.json', 'finite_field')]:
    reader = get_reader(str(OUT / filename), program)
    reader.read_output()
    if filename.startswith('zno'):
        assert np.asarray(reader.raman_tensors).shape == (12, 3, 3)
        assert np.isfinite(reader.raman_tensors).all()
    print(filename, 'read successfully:', reader.nions, 'atoms;', len(reader.frequencies), 'modes')
