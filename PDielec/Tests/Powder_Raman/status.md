# Powder Raman Test Status - April 28, 2026

## Overview
- **Total Tests:** 99
- **Passes:** 94
- **Failures:** 5
- **Location:** `PDielec/Tests/Powder_Raman/`

## Failure Analysis

### 1. Invariant Definition Mismatch
- **Issue:** Discrepancy between Frobenius norm usage in code vs. Placzek invariants expected in tests.
- **Impacted Tests:**
    - `test_a4_antisymmetric_kappa_contribution`: Returns 10.0 instead of 5.0 for unit antisymmetric tensors.
    - `test_e2_diagonal_tensor_VH_converges_to_analytic`: ~50% relative error due to missing 1.5x factor in the $\gamma^2$ analytical formula.

### 2. Numerical Convergence ($n=256$)
- **Issue:** Low sample count with Sobol sequence leads to poor statistics for specific orientations.
- **Impacted Tests:**
    - `test_i2_VH_converges_with_more_samples`: Relative error 0.4253 (limit 0.1).
    - `test_i2_depolarisation_ratio_stable`: Ratio 0.9224 (expected 0.75).
- **Note:** Increasing $n \ge 1024$ resolve these failures.

### 3. Resonance Stability
- **Issue:** Sharp gradients near $\epsilon \approx -2$ (surface-mode resonance) exceed arbitrary test thresholds.
- **Impacted Tests:**
    - `test_i1_N_continuous_across_resonance`: Jump of 9.93 exceeds threshold of 5.0.

## Recommendations
- Synchronize $\gamma^2$ and $\kappa^2$ coefficients between `helpers.py`, `PowderScenarioTab.py`, and the verification plan.
- Relax thresholds for $n=256$ or increase the minimum $n$ in stability tests to ensure reproducibility.
- Re-calibrate resonance jump thresholds to accommodate physical peak gradients.
