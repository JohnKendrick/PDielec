# UnitCell Module Unit Test Plan

## Overview

`UnitCell.py` handles all crystal geometry: lattice vectors, fractional/Cartesian coordinate
conversions, reciprocal lattice, volume, density, symmetry detection via `spglib`, CIF output,
and molecular content analysis (bonds, angles, torsions).  It also exposes the free function
`convert_length_units`.

Most operations have exact mathematical inverses or known analytic results for simple lattices
(cubic, tetragonal, orthorhombic) making them well-suited to unit testing.

---

## Test A — Length Unit Conversion (`test_A_length_units.py`)

### A1: `convert_length_units`

| ID   | Description | Expected |
|------|-------------|----------|
| A1.1 | Same units in and out → value unchanged | exact |
| A1.2 | Bohr → Angstrom: 1 bohr ≈ 0.5291772 Å | within 1e-6 |
| A1.3 | Angstrom → nm: 10 Å = 1 nm | exact |
| A1.4 | nm → Angstrom: 1 nm = 10 Å | exact |
| A1.5 | Angstrom → m: 1 Å = 1e−10 m | within 1e-20 |
| A1.6 | Round-trip: Å → bohr → Å matches input | within 1e-12 |
| A1.7 | Case insensitivity: 'ANG' and 'ang' give same result | exact |
| A1.8 | Input as numpy array: each element converted | element-wise correct |
| A1.9 | Input as list: output is numpy array with correct values | correct |
| A1.10 | 'a.u.', 'au', 'bohr' are synonymous | identical results |

---

## Test B — Lattice Construction and Parameter Round-Trips (`test_B_lattice.py`)

All tests use a `UnitCell` constructed with known parameters.

### B1: Cubic lattice construction

| ID   | Description | Expected |
|------|-------------|----------|
| B1.1 | `a=b=c=5`, `alpha=beta=gamma=90°` → lattice matrix is `5*I` | correct |
| B1.2 | `convert_unitcell_to_abc()` recovers original a,b,c,α,β,γ | within 1e-10 |
| B1.3 | `set_lattice()` with explicit matrix → same result as constructor | equal |

### B2: Orthorhombic lattice

| ID   | Description | Expected |
|------|-------------|----------|
| B2.1 | `a=3, b=4, c=5`, right angles → lattice matrix is `diag([3,4,5])` | correct |
| B2.2 | `convert_unitcell_to_abc()` recovers parameters | within 1e-10 |

### B3: Monoclinic lattice

| ID   | Description | Expected |
|------|-------------|----------|
| B3.1 | `alpha=gamma=90°, beta≠90°` → correct off-diagonal in lattice | analytic |
| B3.2 | Round-trip: construct → `convert_unitcell_to_abc` → `convert_abc_to_unitcell` → same matrix | within 1e-10 |

### B4: `convert_abc_to_unitcell`

| ID   | Description | Expected |
|------|-------------|----------|
| B4.1 | Cubic case → `a*I` | correct |
| B4.2 | Round-trip with `convert_unitcell_to_abc` | within 1e-10 |

### B5: Units in `set_lattice`

| ID   | Description | Expected |
|------|-------------|----------|
| B5.1 | Setting lattice in Bohr then recovering in Angstrom gives 0.5292 * input | within 1e-4 |

---

## Test C — Coordinate Conversions (`test_C_coordinates.py`)

### C1: `convert_abc_to_xyz`

| ID   | Description | Expected |
|------|-------------|----------|
| C1.1 | Cubic lattice: `[1,0,0]` fractional → `[a,0,0]` Cartesian | correct |
| C1.2 | Cubic lattice: `[0.5,0.5,0.5]` → `[a/2,a/2,a/2]` | correct |
| C1.3 | Orthorhombic: each fractional axis independent | no cross-coupling |
| C1.4 | Round-trip: `convert_abc_to_xyz` → `convert_xyz_to_abc` ≈ identity | within 1e-12 |

### C2: `convert_xyz_to_abc`

| ID   | Description | Expected |
|------|-------------|----------|
| C2.1 | Cubic: `[a,0,0]` Cartesian → `[1,0,0]` fractional | correct |
| C2.2 | Origin is fixed point of both conversions | `convert_xyz_to_abc([0,0,0]) == [0,0,0]` |
| C2.3 | Result is 3-element array | correct shape |

### C3: `convert_hkl_to_xyz`

| ID   | Description | Expected |
|------|-------------|----------|
| C3.1 | Cubic: `[1,0,0]` Miller → parallel to first lattice vector | direction correct |
| C3.2 | Cubic: `[1,0,0]` → perpendicular to (100) plane | dot product with off-plane vector is 0 |
| C3.3 | Normalised output (if applicable) | unit vector |

### C4: `convert_xyz_to_integer_abc`

| ID   | Description | Expected |
|------|-------------|----------|
| C4.1 | Exact rational fractional coords → correct integer Miller-like indices | exact |

---

## Test D — Volume and Density (`test_D_volume_density.py`)

### D1: `get_volume`

| ID   | Description | Expected |
|------|-------------|----------|
| D1.1 | Cubic lattice, `a=5 Å`, units='cm': V = (5e-8)³ cm³ | within 1e-30 |
| D1.2 | `get_volume("Bohr")`: convert Å³ to Bohr³ | within 1e-6 relative |
| D1.3 | Orthorhombic `a=3,b=4,c=5`: V = 60 Å³ | 60.0 |
| D1.4 | Volume is positive for any valid lattice | > 0 |

### D2: `get_density`

| ID   | Description | Expected |
|------|-------------|----------|
| D2.1 | Single-atom cubic cell with known mass: density matches ρ = m/V | within 1e-6 relative |
| D2.2 | Density in g/cm³: result is physically reasonable (> 0) | > 0 |
| D2.3 | Empty cell (no atoms): density returns 0 or raises predictably | 0 or handled |

---

## Test E — Symmetry (`test_E_symmetry.py`)

### E1: `find_symmetry`

| ID   | Description | Expected |
|------|-------------|----------|
| E1.1 | Cubic cell with one atom at origin: space group is Pm-3m (221) | spglib number 221 |
| E1.2 | Two-atom FCC-like structure: correct space group | known number |
| E1.3 | Returns a dictionary or structure containing space group info | not None |
| E1.4 | Distorted cell (broken symmetry): lower space group | lower symmetry |

---

## Test F — Molecular Contents (`test_F_molecular.py`)

### F1: `calculate_molecular_contents`

| ID   | Description | Expected |
|------|-------------|----------|
| F1.1 | Two atoms close together (bonded): detected as one molecule | 1 molecule |
| F1.2 | Atoms far apart: each is its own molecule | N separate molecules |
| F1.3 | Water-like geometry: O bonded to 2 H atoms | molecule of size 3 |

### F2: `get_bonds`

| ID   | Description | Expected |
|------|-------------|----------|
| F2.1 | After `calculate_molecular_contents`, bonds are non-empty for bonded atoms | list not empty |
| F2.2 | Bond list is symmetric: if (i,j) present then (j,i) present or at least one of each pair | consistent |

### F3: `get_bond_angles`

| ID   | Description | Expected |
|------|-------------|----------|
| F3.1 | Linear molecule (3 atoms in line): angle ≈ 180° | within 0.01° |
| F3.2 | Right-angle geometry: angle ≈ 90° | within 0.01° |

### F4: `get_torsions`

| ID   | Description | Expected |
|------|-------------|----------|
| F4.1 | Four atoms in a plane: torsion ≈ 0° or 180° | correct |
| F4.2 | Returns list (may be empty for small cells) | list type |

---

## Test G — Coordinate Setting and Retrieval (`test_G_coordinate_io.py`)

### G1: `set_fractional_coordinates` / `get_fractional_coordinates`

| ID   | Description | Expected |
|------|-------------|----------|
| G1.1 | Set and retrieve fractional coords: round-trip | within 1e-12 |
| G1.2 | `get_number_of_atoms()` reflects number of coordinates set | correct count |

### G2: `set_xyz_coordinates` / `get_xyz_coordinates`

| ID   | Description | Expected |
|------|-------------|----------|
| G2.1 | Set Cartesian, retrieve Cartesian: round-trip | within 1e-12 |
| G2.2 | Set Cartesian, retrieve fractional: consistent with lattice | within 1e-12 |

### G3: `set_element_names` / `get_element_names`

| ID   | Description | Expected |
|------|-------------|----------|
| G3.1 | Round-trip: set and get same list | equal |

### G4: `get_species`

| ID   | Description | Expected |
|------|-------------|----------|
| G4.1 | Unique element list for mixed cell | correct set |
| G4.2 | Single-element cell: one species | list of length 1 |

---

## Test H — CIF Output (`test_H_cif.py`)

### H1: `write_cif`

| ID   | Description | Expected |
|------|-------------|----------|
| H1.1 | Writing to a file descriptor produces non-empty output | len > 0 |
| H1.2 | Output contains lattice parameters `_cell_length_a` etc. | keywords present |
| H1.3 | Cubic cell: `_cell_angle_alpha` line contains 90.0 | 90.0 in output |
| H1.4 | Element names appear in the output | correct elements |

---

## File Layout

```
PDielec/Tests/UnitCell/
    unitcell_verification_plan.md        <- this file
    test_A_length_units.py
    test_B_lattice.py
    test_C_coordinates.py
    test_D_volume_density.py
    test_E_symmetry.py
    test_F_molecular.py
    test_G_coordinate_io.py
    test_H_cif.py
    Makefile
```

Run all tests with `make` or `pytest .` from inside the directory.
