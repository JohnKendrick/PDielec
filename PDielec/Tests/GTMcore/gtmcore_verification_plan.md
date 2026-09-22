# GTMcore Module Unit Test Plan

## Overview

`GTMcore.py` implements the generalised 4×4 transfer-matrix (GTM) and scattering-matrix (SM)
methods for multilayer optical stacks.  It is already exercised indirectly by the Crystal Raman
tests, but it has no dedicated unit tests for its mathematical building blocks.

The module contains:
- Pure mathematical utilities (`exact_inv_2x2/3x3/4x4`, `vacuum_eps`)
- `Layer` hierarchy (`Layer`, `CoherentLayer`, `SemiInfiniteLayer`, `IncoherentIntensityLayer`,
  `IncoherentAveragePhaseLayer`, `IncoherentPhaseLayer`, `IncoherentThickLayer`)
- `System` hierarchy (`System`, `TransferMatrixSystem`, `ScatteringMatrixSystem`)
- `SMatrix` class (Redheffer star-product scattering matrix algebra)

Tests are organised from simple mathematical properties up to full optical systems with known
analytic benchmarks.

---

## Test A — Mathematical Utilities (`test_A_math_utilities.py`)

### A1: `vacuum_eps`

| ID   | Description | Expected |
|------|-------------|----------|
| A1.1 | Returns a 3×3 complex ndarray | shape (3,3), complex dtype |
| A1.2 | Off-diagonal elements are zero | 0+0j |
| A1.3 | Diagonal elements are all equal (isotropic vacuum) | all same value |
| A1.4 | Value at any frequency: `eps = 1.0 + 0j` on diagonal | 1.0+0j |

### A2: `exact_inv_2x2`

| ID   | Description | Expected |
|------|-------------|----------|
| A2.1 | Identity 2×2 → identity | `np.eye(2)` |
| A2.2 | `M @ exact_inv_2x2(M) ≈ I` for random invertible M | within 1e-12 |
| A2.3 | Known 2×2 inverse: `[[a,b],[c,d]]` → `[[d,−b],[−c,a]] / det` | analytic |
| A2.4 | Result dtype is complex (if input is complex) | complex |

### A3: `exact_inv_3x3`

| ID   | Description | Expected |
|------|-------------|----------|
| A3.1 | Identity 3×3 → identity | `np.eye(3)` |
| A3.2 | `M @ exact_inv_3x3(M) ≈ I` for random invertible M | within 1e-12 |
| A3.3 | Diagonal matrix: inverse is reciprocal of diagonal | analytic |
| A3.4 | Agrees with `np.linalg.inv` for well-conditioned M | within 1e-10 |

### A4: `exact_inv_4x4`

| ID   | Description | Expected |
|------|-------------|----------|
| A4.1 | Identity 4×4 → identity | `np.eye(4)` |
| A4.2 | `M @ exact_inv_4x4(M) ≈ I` for random invertible M | within 1e-12 |
| A4.3 | Diagonal matrix: inverse is reciprocal of diagonal | analytic |
| A4.4 | Agrees with `np.linalg.inv` for well-conditioned M | within 1e-10 |
| A4.5 | Complex matrix: result matches complex `np.linalg.inv` | within 1e-10 |

---

## Test B — Layer Construction (`test_B_layer.py`)

### B1: `Layer` instantiation

| ID   | Description | Expected |
|------|-------------|----------|
| B1.1 | Default construction: `thickness == 1.0e-6` | 1.0e-6 |
| B1.2 | `set_thickness(d)` → thickness stored | correct |
| B1.3 | `set_euler(theta, phi, psi)` → Euler angles stored | correct |
| B1.4 | `is_coherent()` → True for base Layer | True |
| B1.5 | `set_epsilon` with a constant-function object → `calculate_epsilon` returns correct 3×3 | correct tensor |

### B2: `Layer.calculate_epsilon`

| ID   | Description | Expected |
|------|-------------|----------|
| B2.1 | With `ConstantTensor` DielectricFunction: result is constant at any `f` | constant |
| B2.2 | Result is 3×3 complex ndarray | correct shape/dtype |
| B2.3 | Zero-angle Euler rotation leaves diagonal tensor diagonal | diagonal |

### B3: `Layer.calculate_q` and `calculate_matrices`

| ID   | Description | Expected |
|------|-------------|----------|
| B3.1 | Isotropic layer at normal incidence (`zeta=0`): q has 4 eigenvalues ±n, ±n | symmetric pairs |
| B3.2 | Vacuum layer: eigenvalues are ±1 at normal incidence | ±1+0j |
| B3.3 | `calculate_matrices` returns matrices of correct shape | shape (4,4) |

---

## Test C — SMatrix (Redheffer Product) (`test_C_smatrix.py`)

The Redheffer star product combines two scattering matrices:
- `unit_matrix()` is the identity element: `S ★ I == I ★ S == S`
- For a lossless layer, `|r|² + |t|² == 1` (energy conservation)

### C1: `SMatrix.unit_matrix`

| ID   | Description | Expected |
|------|-------------|----------|
| C1.1 | `unit_matrix()` produces a valid SMatrix | not None |
| C1.2 | Left identity: `S.redheffer(I_mat) == S` | equal |
| C1.3 | Right identity: `I_mat.redheffer(S) == S` | equal |

### C2: `SMatrix.redheffer` (star product)

| ID   | Description | Expected |
|------|-------------|----------|
| C2.1 | Two unit matrices compose to unit matrix | equal to unit |
| C2.2 | Associativity: `(A ★ B) ★ C == A ★ (B ★ C)` for random matrices | within 1e-12 |
| C2.3 | Two identical symmetric layers: transmission should be real and ≤ 1 | bounded |

### C3: `SMatrix.calculate_s`

| ID   | Description | Expected |
|------|-------------|----------|
| C3.1 | Returns a 4×4 complex ndarray | correct shape |

---

## Test D — System Assembly (`test_D_system.py`)

### D1: `System.add_layer` / `del_layer` / `get_layer`

| ID   | Description | Expected |
|------|-------------|----------|
| D1.1 | Empty system: `get_all_layers()` returns empty list | `[]` |
| D1.2 | `add_layer(L)`: `get_all_layers()` has length 1 | 1 |
| D1.3 | `get_layer(0)` returns the first added layer | correct object |
| D1.4 | `del_layer(0)`: list is empty after deletion | `[]` |
| D1.5 | `add_layer` twice → length 2 | 2 |

### D2: `System.get_substrate` / `get_superstrate`

| ID   | Description | Expected |
|------|-------------|----------|
| D2.1 | Set via constructor: `get_substrate()` returns that layer | correct |
| D2.2 | `set_substrate(L)` / `set_superstrate(L)`: round-trip | correct |

### D3: `System.get_layers_boundaries`

| ID   | Description | Expected |
|------|-------------|----------|
| D3.1 | One layer of thickness d: boundaries at [0, d] | correct positions |
| D3.2 | Two layers d1, d2: boundaries at [0, d1, d1+d2] | cumulative |

### D4: `System.reversed_system`

| ID   | Description | Expected |
|------|-------------|----------|
| D4.1 | Reversed system has same number of layers | equal count |
| D4.2 | Layer order is reversed | first becomes last |

---

## Test E — Optical Properties: Fresnel Limits (`test_E_fresnel.py`)

These tests verify the full optical stack against known analytic Fresnel formulae for simple
isotropic systems.

### E1: Normal incidence, single interface (semi-infinite superstrate / semi-infinite substrate)

For two media with refractive indices n₁ and n₂:
```
r = (n₁ - n₂) / (n₁ + n₂)
t = 2n₁ / (n₁ + n₂)
```

| ID   | Description | Expected |
|------|-------------|----------|
| E1.1 | Air/glass (n=1.5) at normal incidence: r ≈ −0.2, |r|² ≈ 0.04 | within 1e-6 |
| E1.2 | n₁ == n₂: r = 0, t = 1 | within 1e-12 |
| E1.3 | Energy conservation: R + T = 1 (where T accounts for refractive index ratio) | within 1e-10 |

### E2: Single transparent layer (Fabry-Pérot fringes)

| ID   | Description | Expected |
|------|-------------|----------|
| E2.1 | Reflectance of a free-standing film with both interfaces air/n=1.5 is periodic in thickness | oscillations present |
| E2.2 | At zero thickness, reflectance matches single interface formula | approaches single-interface value |
| E2.3 | Quarter-wave optical thickness of antireflection film gives lower reflectance | R < single-interface R |

### E3: Absorbing medium

| ID   | Description | Expected |
|------|-------------|----------|
| E3.1 | Increasing imaginary part of n increases reflectance for normal incidence | monotonic |
| E3.2 | Perfect conductor limit (Im(n)→∞): reflectance → 1 | approaches 1 |

---

## Test F — Incoherent Layers (`test_F_incoherent.py`)

### F1: `IncoherentThickLayer`

| ID   | Description | Expected |
|------|-------------|----------|
| F1.1 | `is_coherent()` → False | False |
| F1.2 | For a very thick layer with absorption: transmission → 0 | near 0 |

### F2: `IncoherentAveragePhaseLayer`

| ID   | Description | Expected |
|------|-------------|----------|
| F2.1 | `is_coherent()` → True (averaged phase uses coherent formalism) | True |
| F2.2 | Averaging over many phases should suppress Fabry-Pérot fringes | smoother than coherent |

---

## File Layout

```
PDielec/Tests/GTMcore/
    gtmcore_verification_plan.md        <- this file
    test_A_math_utilities.py
    test_B_layer.py
    test_C_smatrix.py
    test_D_system.py
    test_E_fresnel.py
    test_F_incoherent.py
    Makefile
```

Run all tests with `make` or `pytest .` from inside the directory.
