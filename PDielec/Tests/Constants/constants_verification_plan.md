# Constants Module Unit Test Plan

## Overview

`Constants.py` defines the physical constants and element data used throughout PDielec.  All
constants are taken from the NIST CODATA 2018 recommended values.  Although this module contains
no computational logic, testing it protects against accidental typos in numeric literals, wrong
exponents, or incorrect derived-constant formulae that would silently corrupt all downstream
calculations.

The test suite covers four areas:
- Fundamental constants (verified against CODATA 2018 values)
- Derived constants (verified by re-deriving them from fundamentals)
- Unit-conversion self-consistency (round-trips)
- Element data tables (masses, radii — spot-checked against reference values)

---

## Test A — Fundamental Constants (`test_A_fundamental_constants.py`)

All comparisons use relative tolerance `rtol = 1e-9` unless the constant is defined exactly.

### A1: Fundamental constants from NIST CODATA 2018

| ID   | Constant | Expected value | Source |
|------|----------|---------------|--------|
| A1.1 | `speed_light_si` | 299792458.0 m/s (exact) | CODATA 2018 |
| A1.2 | `planck_si` | 6.62607015e-34 J·Hz⁻¹ (exact) | CODATA 2018 |
| A1.3 | `elementary_charge_si` | 1.602176634e-19 C (exact) | CODATA 2018 |
| A1.4 | `avogadro_si` | 6.02214076e23 mol⁻¹ (exact) | CODATA 2018 |
| A1.5 | `boltzmann_si` | 1.380649e-23 J/K (exact) | CODATA 2018 |
| A1.6 | `mu_0_si` | 1.25663706212e-6 N/A² | CODATA 2018 |
| A1.7 | `electron_mass_si` | 9.1093837015e-31 kg | CODATA 2018 |
| A1.8 | `hartree2ev` | 27.211386245988 eV | CODATA 2018 |

### A2: Expected types

| ID   | Description | Expected |
|------|-------------|----------|
| A2.1 | All fundamental constants are Python `float` or `numpy.float64` | numeric |
| A2.2 | All values are positive | > 0 |

---

## Test B — Derived Constants (`test_B_derived_constants.py`)

Derived constants are re-derived from the fundamental constants and compared to the stored values.
Tolerance: `rtol = 1e-9` unless otherwise specified.

### B1: `hbar_si`

| ID   | Description | Expected |
|------|-------------|----------|
| B1.1 | `hbar_si ≈ planck_si / (2 * pi)` | within 1e-9 relative |

### B2: `molar_gas_si`

| ID   | Description | Expected |
|------|-------------|----------|
| B2.1 | `molar_gas_si ≈ boltzmann_si * avogadro_si` | within 1e-9 relative |

### B3: `fine_structure`

| ID   | Description | Expected |
|------|-------------|----------|
| B3.1 | `fine_structure ≈ 1/137` | within 0.01% |
| B3.2 | Derived from `e²/(4π ε₀ ℏ c)` ≈ stored value | within 1e-8 relative |

### B4: `amu`

| ID   | Description | Expected |
|------|-------------|----------|
| B4.1 | `amu = 1e-3 / avogadro_si` (mass of 1 g/mol in kg) | within 1e-9 relative |
| B4.2 | `amu ≈ 1.66053906660e-27` kg | within 1e-9 relative |

### B5: `wavenumber`

| ID   | Description | Expected |
|------|-------------|----------|
| B5.1 | 1 cm⁻¹ → energy in SI: `wavenumber = planck_si * speed_light_si * 100` | within 1e-9 relative |

### B6: `angs2bohr`

| ID   | Description | Expected |
|------|-------------|----------|
| B6.1 | 1 Å = 1e-10 m; 1 bohr ≈ 0.52918 Å; `angs2bohr ≈ 1/0.52918` | within 1e-5 relative |
| B6.2 | `angs2bohr * 0.5291772` ≈ 1.0 (round-trip via bohr→Å definition) | within 1e-5 |

### B7: `hartree` (in SI)

| ID   | Description | Expected |
|------|-------------|----------|
| B7.1 | `hartree ≈ hartree2ev * elementary_charge_si` | within 1e-9 relative |

### B8: `d2byamuang2`

| ID   | Description | Expected |
|------|-------------|----------|
| B8.1 | Value is positive | > 0 |
| B8.2 | Consistent with `elementary_charge_si², avogadro_si, amu` | self-consistent formula |

---

## Test C — Unit Conversion Self-Consistency (`test_C_unit_conversions.py`)

These tests verify that conversions are mutually consistent without needing to know the exact
CODATA values.

### C1: Energy unit round-trips

| ID   | Description | Expected |
|------|-------------|----------|
| C1.1 | `hartree * avogadro_si` has units of J/mol: positive and large | > 0 |
| C1.2 | `hartree2ev * elementary_charge_si` ≈ `hartree` | within 1e-9 relative |

### C2: Length unit consistency

| ID   | Description | Expected |
|------|-------------|----------|
| C2.1 | `angstrom` (in SI) ≈ 1e-10 m | within 1e-15 |
| C2.2 | `angs2bohr` × 0.5291772 Å/bohr ≈ 1 | within 1e-5 |

### C3: Pressure unit

| ID   | Description | Expected |
|------|-------------|----------|
| C3.1 | `au2GPA` > 0 | positive |
| C3.2 | `au2GPA ≈ 29421.0` GPa per Hartree/Bohr³ | within 0.1% |

### C4: Speed of light

| ID   | Description | Expected |
|------|-------------|----------|
| C4.1 | `speed_light_si * 100` (in cm/s) gives `~3e10 cm/s` | within 1e-9 relative |
| C4.2 | `1 wavenumber / (speed_light_si * 100)` = 1 s (period of 1 cm⁻¹ oscillation) | consistent |

---

## Test D — Element Data Tables (`test_D_element_data.py`)

The tables (`atomic_number_to_element`, `element_to_atomic_number`, `average_masses`,
`isotope_masses`, `covalent_radii`) are spot-checked against reference IUPAC/NIST values.

### D1: `atomic_number_to_element` and `element_to_atomic_number`

| ID   | Description | Expected |
|------|-------------|----------|
| D1.1 | Element 1 is 'H' | 'H' |
| D1.2 | Element 6 is 'C' | 'C' |
| D1.3 | Element 8 is 'O' | 'O' |
| D1.4 | Element 26 is 'Fe' | 'Fe' |
| D1.5 | Inverse mapping: `element_to_atomic_number['H'] == 1` | 1 |
| D1.6 | Forward/inverse round-trip for all elements | consistent |
| D1.7 | No element has atomic number ≤ 0 | all > 0 |

### D2: `average_masses`

| ID   | Description | Expected |
|------|-------------|----------|
| D2.1 | H: ≈ 1.008 u | within 0.001 |
| D2.2 | C: ≈ 12.011 u | within 0.001 |
| D2.3 | O: ≈ 15.999 u | within 0.001 |
| D2.4 | Fe: ≈ 55.845 u | within 0.01 |
| D2.5 | All masses are positive | > 0 |
| D2.6 | No element has mass < 1 u | > 1 (except H, which is ~1) |

### D3: `isotope_masses`

| ID   | Description | Expected |
|------|-------------|----------|
| D3.1 | H most abundant isotope ≈ 1.007825 u | within 0.0001 |
| D3.2 | C-12 isotope ≈ 12.000 u | within 0.001 |
| D3.3 | All entries positive | > 0 |

### D4: `covalent_radii`

| ID   | Description | Expected |
|------|-------------|----------|
| D4.1 | H covalent radius ≈ 0.31 Å | within 0.05 |
| D4.2 | C covalent radius ≈ 0.76 Å | within 0.05 |
| D4.3 | All entries positive | > 0 |
| D4.4 | Generally: covalent radius increases down a group (e.g. Li < Na < K) | monotonic |

---

## File Layout

```
PDielec/Tests/Constants/
    constants_verification_plan.md      <- this file
    test_A_fundamental_constants.py
    test_B_derived_constants.py
    test_C_unit_conversions.py
    test_D_element_data.py
    Makefile
```

Run all tests with `make` or `pytest .` from inside the directory.
