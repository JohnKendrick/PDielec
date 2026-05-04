# Materials and DielectricFunction Unit Test Plan

## Overview

Five test modules covering the `DielectricFunction` hierarchy, all `Material`
subclasses, and the `MaterialsDataBase` spreadsheet loader.  Each module is
self-contained (no conftest fixtures) and exercises both the mathematical
correctness of the permittivity models and the integration layer in
`Materials.py`.

---

## Test A — Constant DielectricFunctions (`test_A_dielectric_constant.py`)

### A1: `ConstantScalar`

| ID   | Description | Expected |
|------|-------------|----------|
| A1.1 | `calculate(v)` returns the same complex value at any frequency | value unchanged |
| A1.2 | `is_scalar()` is True | True |
| A1.3 | `is_tensor()` is False | False |
| A1.4 | `optical_permittivity()` returns `float(real(value))` | float |
| A1.5 | `function()` callable returns same as `calculate()` | equal |

### A2: `ConstantTensor`

| ID   | Description | Expected |
|------|-------------|----------|
| A2.1 | `calculate(v)` returns `value * I + eps_inf` as complex 3x3 | correct tensor |
| A2.2 | `is_scalar()` is False | False |
| A2.3 | `is_tensor()` is True | True |
| A2.4 | `optical_permittivity()` returns `eps_inf` (not the value) | equals eps_inf |
| A2.5 | After `set_epsilon_infinity(eps)`, `calculate()` uses new eps_inf | updated result |

### A3: Base class `DielectricFunction`

| ID   | Description | Expected |
|------|-------------|----------|
| A3.1 | `set_volume()` stores volume in Å³ and converts to a.u. | attrs set |
| A3.2 | `set_epsilon_infinity()` stores 3x3 array | stored correctly |
| A3.3 | `dielectric_contributions_from_modes()` — zero modes returns zero | zero tensor |
| A3.4 | `dielectric_contributions_from_drude()` — far from resonance is small | near zero |

---

## Test B — Tabulated DielectricFunctions (`test_B_dielectric_tabulated.py`)

### B1: `TabulateScalar`

| ID   | Description | Expected |
|------|-------------|----------|
| B1.1 | At each tabulation point `calculate(v)` matches input (within spline tolerance) | ~exact |
| B1.2 | `is_scalar()` is True | True |
| B1.3 | `get_lowest_frequency()` / `get_highest_frequency()` | min/max of vs |
| B1.4 | Interpolated midpoint is between neighbour values | bounded |
| B1.5 | Complex permittivity: real and imaginary parts both interpolated | complex result |

### B2: `Tabulate1` (isotropic tensor)

| ID   | Description | Expected |
|------|-------------|----------|
| B2.1 | Diagonal elements all equal at any frequency | eps[0,0]=eps[1,1]=eps[2,2] |
| B2.2 | `is_tensor()` is True | True |
| B2.3 | `epsilon_infinity` is added to the result | shifted by eps_inf |
| B2.4 | Off-diagonal elements are zero | zero |

### B3: `Tabulate3` (diagonal anisotropic)

| ID   | Description | Expected |
|------|-------------|----------|
| B3.1 | Each diagonal component interpolates its own column | distinct values |
| B3.2 | Off-diagonal elements are zero | zero |
| B3.3 | `epsilon_infinity` is added | shifted by eps_inf |

### B4: `Tabulate6` (full symmetric)

| ID   | Description | Expected |
|------|-------------|----------|
| B4.1 | Diagonal components match their input columns | correct |
| B4.2 | Off-diagonal elements are symmetric (eps[i,j] == eps[j,i]) | symmetric |
| B4.3 | `epsilon_infinity` is added | shifted by eps_inf |

---

## Test C — Lorentzian / Sellmeier DielectricFunctions (`test_C_dielectric_lorentzian.py`)

### C1: `DrudeLorentz` — isotropic

| ID   | Description | Expected |
|------|-------------|----------|
| C1.1 | `is_tensor()` is True | True |
| C1.2 | At f=0: eps[i,i] = eps_inf[i,i] + Σ strength²/v_mode² | analytic |
| C1.3 | Far above resonance (f >> v_mode): contribution approaches 0 | near zero |
| C1.4 | Result is complex tensor | complex ndarray 3x3 |
| C1.5 | Imaginary part is positive near resonance (causal) | Im(eps) > 0 |

### C2: `DrudeLorentz` — rhombic

| ID   | Description | Expected |
|------|-------------|----------|
| C2.1 | xx, yy, zz diagonal values differ at same frequency | distinct |
| C2.2 | Off-diagonal elements are zero | zero |

### C3: `DrudeLorentz` — get/set round-trips

| ID   | Description | Expected |
|------|-------------|----------|
| C3.1 | `set_frequencies` / `get_frequencies` round-trip in cm⁻¹ | unchanged |
| C3.2 | `set_sigmas` / `get_sigmas` round-trip in cm⁻¹ | unchanged |

### C4: `FPSQ` — isotropic

| ID   | Description | Expected |
|------|-------------|----------|
| C4.1 | `is_tensor()` is True | True |
| C4.2 | LST limit at v≈0: eps[i,i] ≈ eps_inf[i,i] × (ω_LO/ω_TO)² | LST |
| C4.3 | Far above all resonances: result approaches eps_inf | near eps_inf |
| C4.4 | Imaginary part of denominator gives positive absorption near ω_TO | Im > 0 |

### C5: `FPSQ` — rhombic (different axes)

| ID   | Description | Expected |
|------|-------------|----------|
| C5.1 | Each diagonal component uses its own TO/LO parameters | distinct |
| C5.2 | Off-diagonal elements are zero | zero |

### C6: `Sellmeier`

| ID   | Description | Expected |
|------|-------------|----------|
| C6.1 | `is_scalar()` is True | True |
| C6.2 | At v=0: result = 1 + Σ Bᵢ | analytic |
| C6.3 | Single-term formula at known λ: n² = 1 + B λ²/(λ²−C) | analytic |
| C6.4 | `optical_permittivity()` evaluated at 589 nm (16978 cm⁻¹) | consistent with formula |
| C6.5 | Result is real | Im(result) = 0 |

### C7: `DFT`

| ID   | Description | Expected |
|------|-------------|----------|
| C7.1 | `calculate()` returns complex 3x3 tensor | shape (3,3) complex |
| C7.2 | `set_frequencies` / `get_frequencies` round-trip (cm⁻¹ ↔ a.u.) | within 1e-10 |
| C7.3 | `set_sigmas` / `get_sigmas` round-trip | within 1e-10 |
| C7.4 | With zero oscillator strengths, result equals eps_inf | equals eps_inf |

---

## Test D — Material subclasses (`test_D_material_classes.py`)

### D1: `Constant`

| ID   | Description | Expected |
|------|-------------|----------|
| D1.1 | `get_name()` returns supplied name | correct |
| D1.2 | `get_density()` returns supplied density | correct |
| D1.3 | `is_scalar()` is True | True |
| D1.4 | Permittivity function returns constant value | constant |
| D1.5 | `get_optical_permittivity()` returns float | float |
| D1.6 | `type` == 'Constant permittivity' | correct |

### D2: `NoMatrix`

| ID   | Description | Expected |
|------|-------------|----------|
| D2.1 | `get_name()` == 'none' | 'none' |
| D2.2 | `get_density()` == 0.0 | 0.0 |
| D2.3 | `is_scalar()` is True | True |
| D2.4 | `type` == 'No matrix (DFT only)' | correct |

### D3: `External`

| ID   | Description | Expected |
|------|-------------|----------|
| D3.1 | Wraps a supplied `ConstantScalar` and delegates `is_scalar()` | True |
| D3.2 | `get_permittivity_object()` returns the supplied object | identity |
| D3.3 | `type` == 'External permittivity' | correct |

### D4: `DrudeLorentz` (Material)

| ID   | Description | Expected |
|------|-------------|----------|
| D4.1 | `is_tensor()` is True | True |
| D4.2 | `get_permittivity_object()` is instance of `DielectricFunction.DrudeLorentz` | correct type |
| D4.3 | `epsilon_infinity` is set correctly in the underlying object | matches input |
| D4.4 | `type` == 'Drude-Lorentz' | correct |
| D4.5 | `get_optical_permittivity()` returns the epsilon infinity tensor | equals eps_inf |

### D5: `FPSQ` (Material)

| ID   | Description | Expected |
|------|-------------|----------|
| D5.1 | `is_tensor()` is True | True |
| D5.2 | `get_permittivity_object()` is instance of `DielectricFunction.FPSQ` | correct type |
| D5.3 | `epsilon_infinity` is set in the underlying object | matches input |
| D5.4 | `type` == 'FPSQ' | correct |

### D6: `Sellmeier` (Material)

| ID   | Description | Expected |
|------|-------------|----------|
| D6.1 | `is_scalar()` is True | True |
| D6.2 | `get_permittivity_object()` is instance of `DielectricFunction.Sellmeier` | correct type |
| D6.3 | `type` == 'Sellmeier' | correct |
| D6.4 | `get_optical_permittivity()` matches formula at 589 nm | consistent |

### D7: `Tabulated` (Material)

| ID   | Description | Expected |
|------|-------------|----------|
| D7.1 | 1D permittivities → `TabulateScalar`, `is_scalar()` True | correct |
| D7.2 | (3, n) permittivities → `Tabulate3`, `is_tensor()` True | correct |
| D7.3 | (6, n) permittivities → `Tabulate6`, `is_tensor()` True | correct |
| D7.4 | `type` == 'Tabulated permittivity' | correct |

### D8: `Material` base helpers

| ID   | Description | Expected |
|------|-------------|----------|
| D8.1 | `set_density` / `get_density` round-trip | unchanged |
| D8.2 | `get_permittivity_function()` is callable | callable |
| D8.3 | `get_information()` returns a non-empty string | str |

---

## Test E — `MaterialsDataBase` (`test_E_database.py`)

### E1: Built-in fallback materials (no spreadsheet needed)

| ID   | Description | Expected |
|------|-------------|----------|
| E1.1 | `get_material("air")` returns `Constant`, permittivity = 1.0 | correct |
| E1.2 | `get_material("vacuum")` returns `Constant`, permittivity = 1.0, density = 0.0 | correct |
| E1.3 | `get_material("ptfe")` returns `Constant`, density = 2.2 | correct |
| E1.4 | `get_material("ldpe")` returns `Constant`, permittivity ≈ 2.25 | correct |
| E1.5 | `get_material("kbr")` returns `Constant`, density = 2.75 | correct |
| E1.6 | `get_material("nujol")` returns `Constant` | correct |
| E1.7 | `get_material("none")` returns `NoMatrix` | correct type |

### E2: `get_sheet_names()` without spreadsheet

| ID   | Description | Expected |
|------|-------------|----------|
| E2.1 | Returns a list (not None) | list |
| E2.2 | Built-in names are present: air, vacuum, ptfe, ldpe, mdpe, kbr, nujol, none | all present |
| E2.3 | List is sorted case-insensitively | sorted |

### E3: Caching

| ID   | Description | Expected |
|------|-------------|----------|
| E3.1 | Two calls to `get_material("air")` return the same object (`is`) | same object |

### E4: Loading the spreadsheet

| ID   | Description | Expected |
|------|-------------|----------|
| E4.1 | `MaterialsDataBase(filename)` loads without raising | no exception |
| E4.2 | `get_sheet_names()` returns more entries than the built-in 8 | len > 8 |
| E4.3 | `get_material(name)` for a known sheet returns the correct type | correct Material subclass |

### E5: Materials from spreadsheet (sanity checks)

| ID   | Description | Expected |
|------|-------------|----------|
| E5.1 | A `Sellmeier` material from the DB has `is_scalar()` True | True |
| E5.2 | An FPSQ material from the DB has `is_tensor()` True | True |
| E5.3 | A tabulated material has `get_information()` containing a frequency range | "freq range" |
| E5.4 | Density of loaded material is positive | > 0 |

---

## File layout

```
PDielec/Tests/Materials/
    materials_verification_plan.md      ← this file
    test_A_dielectric_constant.py
    test_B_dielectric_tabulated.py
    test_C_dielectric_lorentzian.py
    test_D_material_classes.py
    test_E_database.py
    Makefile
```

Run all tests with `make` or `pytest .` from inside the directory.
