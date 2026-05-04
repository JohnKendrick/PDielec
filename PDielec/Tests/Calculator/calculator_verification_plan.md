# Calculator Module Unit Test Plan

## Overview

`Calculator.py` is the largest module in PDielec (~2700 lines) and provides the mathematical
core: tensor utilities, geometry, depolarisation matrices, effective-medium theories, normal-mode
calculations, Euler-angle utilities, and miscellaneous helpers.  Most of its public interface
consists of pure functions with well-defined mathematical inputs and outputs, making it ideal for
unit testing.

The plan is organised into eight test modules, each covering a coherent group of functions.

---

## Test A — Tensor Utilities (`test_A_tensor_utilities.py`)

### A1: `initialise_diagonal_tensor`

| ID   | Description | Expected |
|------|-------------|----------|
| A1.1 | Returns 3×3 real ndarray | `dtype` float, `shape` (3,3) |
| A1.2 | Diagonal elements match input list | exact equality |
| A1.3 | Off-diagonal elements are zero | 0.0 |
| A1.4 | Input `[1,2,3]` → `diag([1,2,3])` | correct values |

### A2: `initialise_complex_diagonal_tensor`

| ID   | Description | Expected |
|------|-------------|----------|
| A2.1 | Returns 3×3 complex ndarray | `dtype` complex, `shape` (3,3) |
| A2.2 | Diagonal elements match input | exact equality |
| A2.3 | Off-diagonal elements are zero | 0.0+0j |

### A3: `initialise_unit_tensor`

| ID   | Description | Expected |
|------|-------------|----------|
| A3.1 | Returns 3×3 real identity matrix | `np.eye(3)` |
| A3.2 | `shape` is (3,3) | correct |

### A4: `average_tensor`

| ID   | Description | Expected |
|------|-------------|----------|
| A4.1 | Average of identity tensor is 1/3 | trace/3 = 1.0 |
| A4.2 | Average of diagonal `[a,b,c]` is `(a+b+c)/3` | analytic |
| A4.3 | Average is a scalar (0-d or length-1 result) | correct type |
| A4.4 | Average of off-diagonal tensor with zero diagonal is 0 | 0.0 |

### A5: `similarity_transform`

| ID   | Description | Expected |
|------|-------------|----------|
| A5.1 | `similarity_transform(I, M) == M` | identity rotation leaves M unchanged |
| A5.2 | `similarity_transform(R, M)` where R is orthogonal: det preserved | `det(result) == det(M)` |
| A5.3 | Diagonal tensor rotated by a known permutation matrix permutes diagonal correctly | exact |

---

## Test B — Geometry Functions (`test_B_geometry.py`)

### B1: `calculate_distance`

| ID   | Description | Expected |
|------|-------------|----------|
| B1.1 | Distance between identical points is 0 | 0.0 |
| B1.2 | Distance along x axis: `[0,0,0]`→`[3,0,0]` | 3.0 |
| B1.3 | Pythagorean triple: `[0,0,0]`→`[3,4,0]` | 5.0 |
| B1.4 | 3D case: `[0,0,0]`→`[1,1,1]` | `sqrt(3)` |
| B1.5 | Result is symmetric: `d(a,b) == d(b,a)` | exact |

### B2: `calculate_angle`

| ID   | Description | Expected |
|------|-------------|----------|
| B2.1 | Collinear points (180°) | π radians or 180° |
| B2.2 | Right angle: `[1,0,0]`, `[0,0,0]`, `[0,1,0]` | π/2 radians |
| B2.3 | Equilateral triangle vertex angle | π/3 radians |
| B2.4 | Result is in correct range [0, π] | bounded |

### B3: `calculate_torsion`

| ID   | Description | Expected |
|------|-------------|----------|
| B3.1 | Trans configuration (180°) | π radians |
| B3.2 | Cis configuration (0°) | 0 radians |
| B3.3 | Gauche+ (60°) | π/3 radians |
| B3.4 | Result is in range [−π, π] | bounded |

---

## Test C — Depolarisation Matrices (`test_C_depolarisation.py`)

All depolarisation matrices L must satisfy:
- `trace(L) == 1` (for a single-component matrix in a single direction)
- All eigenvalues in [0, 1]
- L is symmetric

### C1: `initialise_sphere_depolarisation_matrix`

| ID   | Description | Expected |
|------|-------------|----------|
| C1.1 | Returns 3×3 real ndarray | correct shape |
| C1.2 | Diagonal elements all equal 1/3 | 1/3 |
| C1.3 | Off-diagonal elements are zero | 0.0 |
| C1.4 | `trace(L) == 1` | 1.0 |

### C2: `initialise_plate_depolarisation_matrix`

| ID   | Description | Expected |
|------|-------------|----------|
| C2.1 | Normal along z: L[2,2] == 1, L[0,0] == L[1,1] == 0 | correct |
| C2.2 | Normal along x: L[0,0] == 1, L[1,1] == L[2,2] == 0 | correct |
| C2.3 | `trace(L) == 1` for any normal direction | 1.0 |
| C2.4 | L is symmetric | `L == L.T` |
| C2.5 | Unnormalised normal vector gives same result as normalised | equal |

### C3: `initialise_needle_depolarisation_matrix`

| ID   | Description | Expected |
|------|-------------|----------|
| C3.1 | Unique along z: L[0,0] == L[1,1] == 0.5, L[2,2] == 0 | correct |
| C3.2 | `trace(L) == 1` | 1.0 |
| C3.3 | L is symmetric | `L == L.T` |
| C3.4 | Unique along x: L[1,1] == L[2,2] == 0.5, L[0,0] == 0 | correct |

### C4: `initialise_ellipsoid_depolarisation_matrix`

| ID   | Description | Expected |
|------|-------------|----------|
| C4.1 | `aoverb == 1` (sphere) reproduces sphere result | equal to sphere L |
| C4.2 | `aoverb → 0` (plate limit) approaches plate L | near 0 for transverse components |
| C4.3 | `aoverb → ∞` (needle limit) approaches needle L | near 0 for unique axis |
| C4.4 | `trace(L) == 1` for any `aoverb` | 1.0 |
| C4.5 | L is symmetric | `L == L.T` |

---

## Test D — Permittivity / Refractive Index Conversions (`test_D_permittivity_conversions.py`)

### D1: `calculate_permittivity`

| ID   | Description | Expected |
|------|-------------|----------|
| D1.1 | Refractive index 1+0j → permittivity `[[1,0,0],[0,1,0],[0,0,1]]` | identity |
| D1.2 | n = sqrt(2) (real) → eps = 2.0 on diagonal | 2.0 |
| D1.3 | Complex n = 1+1j → eps = (1+1j)² = 2j on diagonal | correct |
| D1.4 | Returns 3×3 complex ndarray | shape (3,3), complex dtype |

### D2: `calculate_refractive_index`

| ID   | Description | Expected |
|------|-------------|----------|
| D2.1 | eps = 1.0 diagonal → n = 1.0+0j on diagonal | correct |
| D2.2 | eps = 4.0 diagonal → n = 2.0+0j on diagonal | correct |
| D2.3 | Round-trip: `calculate_permittivity(calculate_refractive_index(eps)) ≈ eps` | within 1e-12 |
| D2.4 | Returns 3×3 complex ndarray | correct shape |

### D3: `calculate_refractive_index_scalar`

| ID   | Description | Expected |
|------|-------------|----------|
| D3.1 | eps_scalar = 1.0 → n = 1.0+0j | correct |
| D3.2 | eps_scalar = 4.0 → n = 2.0+0j | correct |
| D3.3 | eps = 2+2j → n matches manual sqrt | within 1e-12 |
| D3.4 | Imaginary part of n is non-negative (physical sign convention) | Im(n) >= 0 |

---

## Test E — Euler Angles (`test_E_euler_angles.py`)

### E1: `euler_rotation_matrix`

| ID   | Description | Expected |
|------|-------------|----------|
| E1.1 | All angles zero → identity matrix | `np.eye(3)` |
| E1.2 | Returns 3×3 real ndarray | correct shape/dtype |
| E1.3 | Result is orthogonal: `R @ R.T ≈ I` | within 1e-12 |
| E1.4 | `det(R) ≈ 1` (proper rotation) | 1.0 |
| E1.5 | Pure rotation by π about z: maps x→−x, y→−y, z→z | correct signs |

### E2: `determine_euler_angles`

| ID   | Description | Expected |
|------|-------------|----------|
| E2.1 | Identity matrix → all angles zero (or equivalent canonical form) | zeroes |
| E2.2 | Round-trip: `euler_rotation_matrix(*determine_euler_angles(R)) ≈ R` for random R | within 1e-10 |

### E3: `euler_rotation`

| ID   | Description | Expected |
|------|-------------|----------|
| E3.1 | All zero angles → vector unchanged | input == output |
| E3.2 | 90° about z: `[1,0,0]` → `[0,1,0]` | correct |
| E3.3 | Rotation preserves vector length | `|Rv| == |v|` |

---

## Test F — Normal Modes and Oscillator Strengths (`test_F_normal_modes.py`)

### F1: `oscillator_strengths`

| ID   | Description | Expected |
|------|-------------|----------|
| F1.1 | Returns array of shape `(nmodes, 3, 3)` | correct shape |
| F1.2 | Symmetric Born charges and diagonal normal modes: result is symmetric | `S == S.T` |
| F1.3 | Zero Born charges → zero oscillator strength | 0.0 |
| F1.4 | Single mode with known Born charges: compare to hand-calculated value | within 1e-12 |

### F2: `normal_modes`

| ID   | Description | Expected |
|------|-------------|----------|
| F2.1 | Returns array of shape `(nmodes, natoms, 3)` | correct shape |
| F2.2 | Each eigenvector is mass-unweighted correctly | verify via formula |

### F3: `infrared_intensities`

| ID   | Description | Expected |
|------|-------------|----------|
| F3.1 | Returns 1D array of length `nmodes` | correct shape |
| F3.2 | Zero oscillator strengths → zero intensities | 0.0 |
| F3.3 | All non-negative values | >= 0 |

### F4: `raman_intensities`

| ID   | Description | Expected |
|------|-------------|----------|
| F4.1 | Returns list with correct number of elements | correct length |
| F4.2 | Zero Raman tensors → zero intensities | 0.0 |

### F5: `calculate_normal_modes_and_frequencies`

| ID   | Description | Expected |
|------|-------------|----------|
| F5.1 | Symmetric positive-definite Hessian → all real positive frequencies | Im(f) = 0, f > 0 |
| F5.2 | Mass-spring Hessian with known analytic frequencies | within 1e-10 |
| F5.3 | Returns `(frequencies, normal_modes)` with correct shapes | correct |

---

## Test G — Effective Medium Theories (`test_G_effective_medium.py`)

For all EMTs, physical constraints at any volume fraction `vf`:
- Result is a complex 3×3 tensor
- In the limit `vf → 0`, result approaches the medium permittivity

### G1: `averaged_permittivity`

| ID   | Description | Expected |
|------|-------------|----------|
| G1.1 | `vf = 0` → medium permittivity | equal to eps_medium |
| G1.2 | Isotropic inputs → isotropic (diagonal, equal elements) output | diagonal |
| G1.3 | Result is complex 3×3 | correct shape |

### G2: `maxwell`

| ID   | Description | Expected |
|------|-------------|----------|
| G2.1 | `vf = 0` → medium permittivity | equal to eps_medium |
| G2.2 | Sphere shape with eps_crystal = eps_medium → result = eps_medium (trivial) | equal |
| G2.3 | Result is complex 3×3 | correct shape |

### G3: `maxwell_sihvola`

| ID   | Description | Expected |
|------|-------------|----------|
| G3.1 | `vf = 0` → medium permittivity | equal to eps_medium |
| G3.2 | Sphere shape with identical eps → trivial result | equal |

### G4: `balan`

| ID   | Description | Expected |
|------|-------------|----------|
| G4.1 | `vf = 0` → medium permittivity | equal to eps_medium |
| G4.2 | Returns complex 3×3 tensor | correct shape |

### G5: `compute_internal_field_tensor` / `compute_particle_raman_tensor`

| ID   | Description | Expected |
|------|-------------|----------|
| G5.1 | `eps_i == eps_e` (no contrast) → internal field tensor is identity | `I` |
| G5.2 | `compute_particle_raman_tensor` with identity internal field → equals crystal Raman tensor | equal |

---

## Test H — Absorption Spectrum (`test_H_absorption_spectrum.py`)

### H1: `absorption_from_mode_intensities`

| ID   | Description | Expected |
|------|-------------|----------|
| H1.1 | Zero intensities → zero spectrum everywhere | 0.0 |
| H1.2 | Single mode: peak near mode frequency | max near mode freq |
| H1.3 | Lorentzian peak area proportional to intensity | analytic |
| H1.4 | Spectrum is non-negative | >= 0 |
| H1.5 | Result has same length as input frequency array | correct shape |

### H2: `calculate_size_factor`

| ID   | Description | Expected |
|------|-------------|----------|
| H2.1 | `x = 0` → 1.0 (no size correction) | 1.0 |
| H2.2 | Non-negative for real positive `x` | >= 0 |

---

## File Layout

```
PDielec/Tests/Calculator/
    calculator_verification_plan.md     <- this file
    test_A_tensor_utilities.py
    test_B_geometry.py
    test_C_depolarisation.py
    test_D_permittivity_conversions.py
    test_E_euler_angles.py
    test_F_normal_modes.py
    test_G_effective_medium.py
    test_H_absorption_spectrum.py
    Makefile
```

Run all tests with `make` or `pytest .` from inside the directory.
