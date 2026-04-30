# Verification Plan for Crystal Raman Implementation in Layered Systems

## 1. Purpose

This document gives a structured plan for verifying an implementation of crystal Raman scattering from layered systems and single-crystal slabs. The aim is to test the implementation in stages, beginning with cases where the answer is known analytically and progressing towards realistic layered structures.

The implementation is assumed to calculate Raman strength using optical fields obtained from a layered Maxwell solver such as `pyGTM`, with a Raman amplitude of the general form

```markdown
$$
A_m \propto \int_{z_1}^{z_2}
\mathbf{E}_{s}^{\mathrm{rec}}(z)^T
R_m(z)
\mathbf{E}_{i}(z)\,dz,
$$
```

where:

- `m` labels the Raman-active phonon mode,
- `R_m` is the Raman tensor for mode `m`,
- `E_i(z)` is the incident optical field at the laser frequency,
- `E_s^rec(z)` is the reciprocal or time-reversed collection field at the scattered frequency,
- the integral is over the Raman-active material.

The goal is not only to check that the code runs, but to verify that the tensor rotations, optical fields, layer indexing, depth integration, polarisation handling, reciprocity and limiting behaviours are all correct.

---

## 2. Overall Verification Strategy

Use a ladder of tests:

1. **Raman tensor algebra without any layered optics.**
2. **Single homogeneous slab with trivial or nearly trivial optical fields.**
3. **Single isotropic slab with analytically known Fresnel fields.**
4. **Layered systems with internal consistency checks.**
5. **Reciprocity and layer-reversal checks.**
6. **Complex permittivity and complex Raman tensor checks.**
7. **Realistic validation examples.**

Each level should be added to the automated test suite only after the previous level has passed.

---

# Part A: Raman Tensor Algebra Tests

These tests do not require `pyGTM` or any Maxwell solver. They check the algebra of Raman tensors, rotations and polarisation contractions.

---

## A1. Direct Raman Tensor Contraction

### Objective

Check that the implementation computes the basic Raman amplitude correctly:

```math
A = \mathbf{e}_s^T R \mathbf{e}_i.
```

### Test System

Use a diagonal Raman tensor

```math
R =
\begin{pmatrix}
a & 0 & 0 \\
0 & b & 0 \\
0 & 0 & c
\end{pmatrix}.
```

Use simple Cartesian polarisations:

```math
\mathbf{e}_x = (1,0,0)^T,
\quad
\mathbf{e}_y = (0,1,0)^T,
\quad
\mathbf{e}_z = (0,0,1)^T.
```

### Expected Results

| Incident polarisation | Scattered polarisation | Expected amplitude | Expected strength |
|---|---|---:|---:|
| `x` | `x` | `a` | `|a|^2` |
| `y` | `y` | `b` | `|b|^2` |
| `z` | `z` | `c` | `|c|^2` |
| `x` | `y` | `0` | `0` |
| `x` | `z` | `0` | `0` |

### Pass Criterion

The computed amplitude should agree with the analytic value to numerical precision.

Suggested tolerance:

```text
absolute error < 1.0e-12
```

for double precision arithmetic.

---

## A2. Raman Tensor Rotation About the z Axis

### Objective

Check the convention for rotating Raman tensors between crystal and laboratory frames.

If `G` maps crystal-frame vectors into lab-frame vectors, the lab-frame Raman tensor should be

```math
R^{\mathrm{lab}} = G R^{\mathrm{crystal}} G^T.
```

### Test System

Use

```math
R^{\mathrm{crystal}} =
\begin{pmatrix}
a & 0 & 0 \\
0 & b & 0 \\
0 & 0 & c
\end{pmatrix}
```

and rotate the crystal by angle `phi` about the lab `z` axis.

The expected lab-frame tensor components are

```math
R_{xx}^{\mathrm{lab}} = a\cos^2\phi + b\sin^2\phi,
```

```math
R_{yy}^{\mathrm{lab}} = a\sin^2\phi + b\cos^2\phi,
```

```math
R_{xy}^{\mathrm{lab}} = (a-b)\sin\phi\cos\phi.
```

### Expected Raman Strengths

For `x` incident and `x` analysed scattered polarisation:

```math
I_{xx}(\phi) \propto
\left|a\cos^2\phi + b\sin^2\phi\right|^2.
```

For `x` incident and `y` analysed scattered polarisation:

```math
I_{xy}(\phi) \propto
\left|(a-b)\sin\phi\cos\phi\right|^2.
```

### Concrete Numerical Example

Use

```math
R =
\begin{pmatrix}
2 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 0.5
\end{pmatrix}.
```

Then

```math
I_{xx}(\phi) \propto
\left|2\cos^2\phi + \sin^2\phi\right|^2,
```

and

```math
I_{xy}(\phi) \propto
\left|\sin\phi\cos\phi\right|^2.
```

Expected values:

| `phi` | `I_xx` proportional to | `I_xy` proportional to |
|---:|---:|---:|
| `0°` | `4.00` | `0.00` |
| `45°` | `2.25` | `0.25` |
| `90°` | `1.00` | `0.00` |

### Pass Criterion

The numerical result should agree with the analytic expressions within floating-point precision.

---

## A3. Isotropic Raman Tensor Orientation Invariance

### Objective

Check that an isotropic Raman tensor is unaffected by crystal orientation.

### Test System

Use

```math
R = a I.
```

Rotate the crystal by many random rotations.

### Expected Result

The transformed tensor should remain

```math
R^{\mathrm{lab}} = a I.
```

For orthogonal incident and scattered polarisations:

```math
\mathbf{e}_s^T R \mathbf{e}_i = 0.
```

For parallel polarisations:

```math
\mathbf{e}_i^T R \mathbf{e}_i = a.
```

### Pass Criterion

The Raman strength should be independent of orientation.

This is a sensitive test of whether the rotation operation has been applied in the correct order.

---

## A4. Off-Diagonal Raman Tensor Selection Rule

### Objective

Check crossed-polarisation selection rules.

### Test System

Use

```math
R =
\begin{pmatrix}
0 & d & 0 \\
d & 0 & 0 \\
0 & 0 & 0
\end{pmatrix}.
```

### Expected Results at Normal Incidence

| Incident polarisation | Scattered polarisation | Expected amplitude |
|---|---|---:|
| `x` | `x` | `0` |
| `x` | `y` | `d` |
| `y` | `x` | `d` |
| `y` | `y` | `0` |

### Additional Rotation Test

Rotate the crystal by `45°` about `z`. Intensity should transfer between parallel and crossed polarisation channels.

### Pass Criterion

The zero channels should be numerically zero and the non-zero channels should match the analytic values.

---

## A5. Complex Phase Invariance of the Raman Tensor

### Objective

Check that the final Raman strength is invariant under a global complex phase applied to the Raman tensor.

### Test

Replace

```math
R \rightarrow e^{i\theta} R.
```

### Expected Result

The amplitude changes phase:

```math
A \rightarrow e^{i\theta} A,
```

but the strength is unchanged:

```math
|A|^2 \rightarrow |A|^2.
```

### Pass Criterion

For several values of `theta`, the computed Raman strength should remain unchanged.

Suggested tolerance:

```text
relative error < 1.0e-12
```

---

# Part B: Single-Layer Optical Field Tests

These tests use the layered Maxwell solver, but only for the simplest possible optical systems.

---

## B1. Uniform-Field Single-Crystal Limit

### Objective

Check that the full layered Raman code reduces to the ordinary single-crystal Raman expression when the optical fields are uniform.

### Test System

Use a single Raman-active layer with:

- negligible reflection,
- no absorption,
- constant real permittivity,
- normal incidence,
- simple Cartesian Raman tensor.

The optical fields should be approximately constant through the layer.

### Expected Result

The Raman amplitude should reduce to

```math
A_m \propto d\,\mathbf{e}_s^T R_m \mathbf{e}_i,
```

where `d` is the slab thickness.

The Raman strength should therefore be proportional to

```math
I_m \propto d^2\left|\mathbf{e}_s^T R_m \mathbf{e}_i\right|^2
```

for coherent depth integration.

If the code intentionally performs incoherent depth integration, the expected scaling is instead

```math
I_m \propto d\left|\mathbf{e}_s^T R_m \mathbf{e}_i\right|^2.
```

### Pass Criterion

The implementation must reproduce the expected thickness scaling for the chosen coherent or incoherent convention.

This test should be used to confirm and document exactly which convention the implementation uses.

---

## B2. Very Thin Slab Limit

### Objective

Check that the code approaches the uniform-field result when the slab thickness is much smaller than the optical wavelength.

### Test Condition

Use

```math
d \ll \lambda.
```

For example, for visible light with `lambda = 500 nm`, use layer thicknesses such as:

```text
0.1 nm, 0.5 nm, 1 nm, 2 nm
```

### Expected Result

The field should be nearly constant across the layer and the amplitude should satisfy

```math
A_m(d) \approx d\,
\mathbf{E}_s^{\mathrm{rec}}(0)^T
R_m
\mathbf{E}_i(0).
```

### Pass Criterion

For coherent integration:

```math
I_m(d) / d^2
```

should be approximately constant in the very thin slab limit.

For incoherent integration:

```math
I_m(d) / d
```

should be approximately constant.

---

## B3. Thickness Scaling with No Interference

### Objective

Verify that the depth integration is performed correctly when the optical fields do not vary significantly with depth.

### Test

Calculate Raman strength for a sequence of thicknesses:

```text
d, 2d, 4d, 8d
```

with weak reflection and negligible absorption.

### Expected Result

For coherent integration:

```math
I(2d) \approx 4 I(d).
```

For incoherent integration:

```math
I(2d) \approx 2 I(d).
```

### Pass Criterion

The observed scaling should match the intended integration convention.

---

## B4. Absorbing Thick Slab Saturation

### Objective

Check that absorption is treated with the correct sign and that the depth integration saturates for thick absorbing samples.

### Model Behaviour

If the incident field decays approximately as

```math
E_i(z) \sim e^{-\alpha_i z/2},
```

and the reciprocal scattered field as

```math
E_s^{\mathrm{rec}}(z) \sim e^{-\alpha_s z/2},
```

then the Raman amplitude should roughly behave as

```math
A(d) \propto
\int_0^d e^{-(\alpha_i+\alpha_s)z/2}\,dz
=
\frac{1-e^{-(\alpha_i+\alpha_s)d/2}}
{(\alpha_i+\alpha_s)/2}.
```

### Expected Result

As the slab becomes thick compared with the penetration depth, the Raman signal should approach a finite limiting value.

### Pass Criterion

Increasing the thickness beyond several optical penetration depths should not continue to increase the signal indefinitely.

Also check that increasing the imaginary part of the refractive index decreases or saturates the signal, rather than increasing it unphysically.

---

## B5. Scalar Isotropic Slab Compared with Independent Fresnel Theory

### Objective

Check the optical field part of the calculation independently of `pyGTM`.

### Test System

Use an isotropic film at normal incidence between two isotropic media.

For this case, the electromagnetic problem reduces to the standard scalar transfer-matrix/Fabry--Perot problem.

### Independent Benchmark

Write a small independent scalar transfer-matrix script to compute:

- the incident field profile `E_i(z)`,
- the reciprocal scattered field profile `E_s^rec(z)`,
- the Raman overlap integral.

At normal incidence, the field in the layer can be written as

```math
E_i(z) = A_i e^{ik_i z} + B_i e^{-ik_i z},
```

and similarly for the reciprocal scattered field:

```math
E_s^{\mathrm{rec}}(z) = A_s e^{ik_s z} + B_s e^{-ik_s z}.
```

The Raman amplitude is then a sum of analytic terms such as

```math
\int_0^d e^{i(k_i+k_s)z}\,dz
=
\frac{e^{i(k_i+k_s)d}-1}{i(k_i+k_s)}.
```

### Expected Result

The `pyGTM`-based calculation should agree with the independent scalar transfer-matrix calculation for isotropic normal-incidence cases.

### Pass Criterion

For non-resonant cases, aim for a relative error below

```text
1.0e-8 to 1.0e-6
```

depending on field normalisation and numerical integration accuracy.

---

# Part C: Layered-System Consistency Tests

These tests check whether the layered implementation behaves consistently under transformations that should not change the physical result.

---

## C1. Raman-Inactive Layer Test

### Objective

Check that layers with zero Raman tensor do not contribute directly to the Raman signal.

### Test

Construct a multilayer system with one or more layers having

```math
R_m = 0.
```

### Expected Result

These layers should contribute zero Raman amplitude.

They may still affect the optical fields in neighbouring Raman-active layers, but they should not generate Raman signal themselves.

### Pass Criterion

If all layers have `R_m = 0`, the total Raman strength should be exactly zero within numerical precision.

---

## C2. Layer Subdivision Invariance

### Objective

Check that artificial subdivision of a homogeneous layer does not change the result.

This is one of the most important regression tests for a layered Raman implementation.

### Test

Compute the Raman strength for a single Raman-active layer of thickness `d`, permittivity `epsilon`, and Raman tensor `R`.

Then replace the layer by two adjacent identical layers:

```math
d/2 + d/2.
```

Then repeat with 5, 10 and 20 identical sublayers.

### Expected Result

The Raman strength should be unchanged by subdivision.

### Pass Criterion

The relative difference should decrease with integration accuracy and should be small:

```text
relative error < 1.0e-8 to 1.0e-6
```

for well-conditioned cases.

### Errors This Test Can Reveal

- Incorrect layer indexing.
- Incorrect local `z` coordinate inside layers.
- Incorrect phase reference at interfaces.
- Failure to maintain field continuity at artificial boundaries.
- Numerical quadrature errors.

---

## C3. Origin-Shift Invariance

### Objective

Check that shifting the coordinate origin does not change physical intensities.

### Test

Represent the same physical stack using different choices of `z = 0`.

### Expected Result

The complex Raman amplitude may acquire an overall phase depending on convention, but the Raman strength should be invariant.

### Pass Criterion

The intensity should be unchanged within numerical tolerance.

---

## C4. Symmetric Stack Illumination Test

### Objective

Check left/right symmetry.

### Test System

Use a symmetric structure such as

```text
air / film / air
```

or

```text
air / spacer / Raman layer / spacer / air
```

with the Raman-active layer placed symmetrically.

Calculate the Raman signal for illumination from the left and from the right.

### Expected Result

For a symmetric stack, the Raman strength should be the same for illumination from either side, provided the incident and collection channels are transformed consistently.

### Pass Criterion

Left and right illumination should give the same strength within numerical tolerance.

---

## C5. Reversed Stack with Reversed Illumination

### Objective

Check consistency of layer order and propagation direction.

### Test

Take an asymmetric stack:

```text
medium A / layer 1 / Raman layer / layer 2 / medium B
```

Calculate the Raman strength for illumination from the left.

Then reverse the stack:

```text
medium B / layer 2 / Raman layer / layer 1 / medium A
```

and illuminate from the corresponding reversed side.

### Expected Result

The two calculations should describe the same physical situation and should give the same Raman strength, subject to the same source and detector normalisation.

### Pass Criterion

Agreement within numerical tolerance after accounting for field normalisation conventions.

---

## C6. Moving Raman Layer Through a Standing Wave

### Objective

Check whether the Raman signal follows the local optical field intensity/overlap in a multilayer structure.

### Test

Construct a transparent cavity or weak Fabry--Perot structure with a Raman-active thin layer placed at different positions inside the cavity.

### Expected Result

The Raman signal should oscillate as the Raman-active layer moves through nodes and antinodes of the optical overlap field

```math
\mathbf{E}_s^{\mathrm{rec}}(z)^T R \mathbf{E}_i(z).
```

### Pass Criterion

The signal should be largest near antinodes of the overlap field and smallest near nodes.

---

# Part D: Reciprocity and Reciprocal-Field Tests

These tests are especially important if the implementation uses the reciprocal-field or time-reversed-field method for the outgoing radiation.

---

## D1. Incident/Collection Reciprocity Test

### Objective

Check that the reciprocal-field method is implemented consistently.

### Test

1. Calculate the Raman amplitude for incidence from side A and collection into side B.
2. Swap source and detector channels in the reciprocal problem.
3. Compare the resulting amplitude or strength after accounting for normalisation.

### Expected Result

For reciprocal, non-magnetic media, the source/detector interchange should obey the expected reciprocity relation.

Depending on the exact field normalisation, the amplitudes may differ by a known factor, but the physical strength should be consistent.

### Pass Criterion

The two calculations should agree after applying the known normalisation convention.

---

## D2. No-Conjugation Versus Conjugation Test

### Objective

Check the complex-field convention in the Raman overlap integral.

In many Lorentz reciprocity formulations, the overlap uses an unconjugated transpose:

```math
\mathbf{E}_s^{\mathrm{rec},T} R \mathbf{E}_i,
```

rather than a Hermitian conjugate:

```math
\mathbf{E}_s^{\mathrm{rec},\dagger} R \mathbf{E}_i.
```

### Test

Use a lossy layer or a case with evanescent/complex fields.

Compare the implemented result with a known analytic or independent benchmark.

### Expected Result

The correct convention should reproduce reciprocity and the known limiting behaviour.

### Pass Criterion

This test should be documented carefully because the result may depend on the theoretical convention chosen in the implementation.

---

# Part E: Frequency-Dependence Tests

The incident and scattered optical fields are evaluated at different frequencies:

```math
\omega_L
```

and

```math
\omega_S = \omega_L - \omega_m.
```

These tests check that the implementation uses the correct frequency in each part of the calculation.

---

## E1. Zero Raman Shift Test

### Objective

Check the limiting case where the scattered and incident optical frequencies are identical.

### Test

Set

```math
\omega_m = 0,
```

so that

```math
\omega_S = \omega_L.
```

### Expected Result

The optical field solutions for the incident and reciprocal scattered frequencies should be identical apart from propagation direction and boundary conditions.

### Pass Criterion

The Raman strength should agree with the same-frequency benchmark.

---

## E2. Artificially Large Raman Shift Test

### Objective

Check that the scattered optical field is actually being evaluated at the scattered frequency, not accidentally at the incident frequency.

### Test

Use an artificially large Raman shift so that

```math
\omega_S
```

differs significantly from

```math
\omega_L.
```

Use a dispersive refractive index.

### Expected Result

The scattered reciprocal field profile should change, and therefore the Raman overlap should change.

### Pass Criterion

The result should differ from the zero-shift result in a way consistent with the change in optical phase thickness.

---

## E3. Non-Dispersive Control Test

### Objective

Separate optical field effects from explicit frequency prefactors.

### Test

Repeat E2 using a non-dispersive refractive index.

### Expected Result

If the refractive index is independent of frequency, the field-overlap part should change only through trivial changes in wavelength/phase, not through material dispersion.

### Pass Criterion

The behaviour should match the expected non-dispersive transfer-matrix result.

---

# Part F: Complex Permittivity and Complex Raman Tensor Tests

These tests check behaviour in absorbing and complex-valued systems.

---

## F1. Complex Permittivity Sign Convention

### Objective

Check that absorption has the correct sign.

### Test

Use a material with complex refractive index

```math
n = n' + i n''
```

or complex permittivity

```math
\epsilon = \epsilon' + i\epsilon''.
```

Increase `n''` or `epsilon''` gradually.

### Expected Result

The field should decay into the absorbing material.

The Raman signal from a thick slab should decrease or saturate as absorption increases.

### Pass Criterion

If increasing absorption increases the effective sampling depth or gives exponentially growing fields, the sign convention is wrong.

---

## F2. Complex Raman Tensor Phase Invariance in the Full Layered Code

### Objective

Extend the simple tensor phase test to the full layered implementation.

### Test

Apply a global phase to all Raman tensor components:

```math
R_m \rightarrow e^{i\theta}R_m.
```

### Expected Result

The Raman amplitude should acquire the same phase, but the strength should not change.

### Pass Criterion

The Raman strength should be invariant for several values of `theta`.

---

## F3. Complex Symmetric Tensor Test

### Objective

Check handling of complex Raman tensors that are symmetric but not Hermitian.

### Test

Use a tensor such as

```math
R =
\begin{pmatrix}
1+i & 0.2i & 0 \\
0.2i & 2-i & 0 \\
0 & 0 & 0.5+0.1i
\end{pmatrix}.
```

### Expected Result

The code should treat this as a symmetric Raman tensor, not as a Hermitian operator.

### Pass Criterion

Results should match a direct complex tensor contraction benchmark.

---

# Part G: Physical Validation Tests

These are not purely analytic tests, but they check whether the implementation gives physically sensible behaviour.

---

## G1. Single-Crystal Selection Rules

### Objective

Check the implementation against known single-crystal Raman selection rules.

### Suggested Approach

Use a simple high-symmetry crystal with known Raman tensors, or use artificial tensors with the same symmetry structure.

Examples:

- diagonal tensor mode,
- off-diagonal tensor mode,
- `zz`-only tensor mode,
- degenerate pair of modes.

### Expected Result

The calculated polarisation dependence should match the known selection rules in the appropriate scattering geometry.

---

## G2. Pure `zz` Raman Tensor at Normal and Oblique Incidence

### Objective

Check that the code correctly uses the full vector optical field, including possible `z` components.

### Test Tensor

Use

```math
R =
\begin{pmatrix}
0 & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & c
\end{pmatrix}.
```

### Expected Result

At normal incidence in an isotropic medium, transverse fields have no `z` component, so the mode should be absent.

At oblique incidence, especially for p-polarised light, the field may have a `z` component and the mode may become visible.

### Pass Criterion

The intensity should be zero at normal incidence and become non-zero only when the electromagnetic field has a `z` component.

---

## G3. Thin Film on Substrate: Interference Enhancement

### Objective

Check whether the implementation reproduces expected interference effects in a film on a substrate.

### Test

Use a Raman-active transparent film on a non-absorbing substrate and calculate the Raman strength as a function of film thickness.

### Expected Result

The Raman strength should oscillate with thickness due to Fabry--Perot interference.

The oscillation period should be related to the optical phase thickness of the film.

### Pass Criterion

The thickness dependence should agree with a scalar transfer-matrix benchmark for isotropic normal-incidence cases.

---

## G4. Absorbing Film on Substrate

### Objective

Check physically sensible behaviour in an absorbing film.

### Test

Repeat G3 with a finite imaginary part of the refractive index.

### Expected Result

The interference oscillations should be damped and the signal should approach a finite limit as the film becomes optically thick.

---

# Part H: Numerical Convergence Tests

These tests should be run for several representative systems.

---

## H1. Depth Integration Convergence

### Objective

Check convergence with the number of quadrature points used inside each layer.

### Test

For a representative slab or multilayer system, compute the Raman strength using increasing numbers of integration points:

```text
10, 20, 50, 100, 200, 500
```

### Expected Result

The result should converge smoothly.

### Pass Criterion

Choose a default integration density such that further refinement changes the result by less than a chosen tolerance, for example:

```text
relative change < 1.0e-6
```

for routine tests, or

```text
relative change < 1.0e-4
```

for difficult resonant cases.

---

## H2. Layer Subdivision Convergence

### Objective

Check that artificially subdividing a layer does not change the result.

This overlaps with C2, but should also be used as a convergence test.

### Test

Subdivide a homogeneous Raman-active layer into:

```text
1, 2, 5, 10, 20, 50
```

identical sublayers.

### Expected Result

The result should remain constant.

---

## H3. Resonant Fabry--Perot Case

### Objective

Check numerical stability near optical resonances.

### Test

Use a transparent film thickness near a Fabry--Perot resonance.

### Expected Result

The Raman signal may become large or rapidly varying, but it should converge with integration density and field-solver precision.

### Pass Criterion

No numerical instability or spurious discontinuity should appear as thickness is varied smoothly.

---

# Part I: Suggested Automated Test Suite

The following tests are recommended for automated regression testing.

---

## Fast Unit Tests

These should run every time the code is tested.

| Test ID | Name | Expected Runtime | Purpose |
|---|---|---:|---|
| A1 | Direct tensor contraction | Very fast | Checks basic amplitude |
| A2 | Rotation about `z` | Very fast | Checks rotation convention |
| A3 | Isotropic tensor invariance | Very fast | Checks orientation invariance |
| A4 | Off-diagonal selection rule | Very fast | Checks crossed polarisations |
| A5 | Complex phase invariance | Very fast | Checks complex tensor handling |
| B1 | Uniform-field slab | Fast | Checks full code in trivial optical limit |
| C1 | Zero Raman tensor | Fast | Checks inactive layers |
| C2 | Layer subdivision | Fast/moderate | Checks layer indexing and phase |

---

## Medium Tests

These can run in a more complete test suite.

| Test ID | Name | Purpose |
|---|---|---|
| B4 | Absorbing thick slab | Checks absorption sign and saturation |
| B5 | Scalar Fresnel comparison | Independent optical benchmark |
| C3 | Origin-shift invariance | Checks phase-reference consistency |
| C4 | Symmetric stack illumination | Checks left/right handling |
| D1 | Reciprocity | Checks reciprocal-field implementation |
| E1 | Zero Raman shift | Checks incident/scattered frequency handling |
| F2 | Full-code complex phase invariance | Checks complex full implementation |

---

## Slow or Validation Tests

These can be run before releases or major changes.

| Test ID | Name | Purpose |
|---|---|---|
| C6 | Moving Raman layer through standing wave | Checks spatial overlap physics |
| E2 | Large Raman shift | Checks frequency dependence |
| G3 | Thin film interference | Physical validation |
| G4 | Absorbing film on substrate | Physical validation |
| H3 | Resonant Fabry--Perot convergence | Numerical robustness |

---

# Part J: Recommended Development Order

The recommended order for implementing the verification tests is:

1. Implement tests A1--A5 for tensor algebra.
2. Implement B1 and B2 for the uniform-field and thin-slab limits.
3. Implement C1 and C2 for inactive layers and layer subdivision.
4. Implement B5 using an independent scalar Fresnel/transfer-matrix benchmark.
5. Implement B4 for absorption and thick-slab saturation.
6. Implement D1 for reciprocity.
7. Implement E1--E3 for incident/scattered frequency handling.
8. Implement F1--F3 for complex-valued tensors and permittivities.
9. Add G1--G4 as physical validation examples.
10. Add H1--H3 as systematic convergence tests.

---

# Part K: Particularly Important Tests

If only a small number of tests can be implemented initially, prioritise these:

## 1. Rotated Diagonal Raman Tensor

This checks:

- tensor rotation,
- polarisation contraction,
- lab/crystal frame convention,
- angular dependence.

Use

```math
R =
\begin{pmatrix}
2 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 0.5
\end{pmatrix}.
```

Expected values:

| `phi` | `I_xx` | `I_xy` |
|---:|---:|---:|
| `0°` | `4.00` | `0.00` |
| `45°` | `2.25` | `0.25` |
| `90°` | `1.00` | `0.00` |

## 2. Layer Subdivision Invariance

This checks:

- layer indexing,
- phase references,
- field continuity,
- depth integration.

A homogeneous Raman-active layer and the same layer split into many identical sublayers must give the same result.

## 3. Absorbing Thick Slab Saturation

This checks:

- sign of complex permittivity,
- physical decay of fields,
- finite sampling depth,
- correct treatment of absorption.

## 4. Independent Scalar Transfer-Matrix Benchmark

This checks:

- the optical fields from the layered solver,
- the Raman overlap integral,
- the frequency dependence of incident and scattered fields.

This is the best way to avoid testing the implementation only against itself.

---

# Part L: Documentation Recommendations

The implementation documentation should explicitly state:

1. Whether Raman amplitudes or intensities are integrated over depth.
2. Whether depth integration is coherent or incoherent.
3. The convention for rotating Raman tensors.
4. The convention for the reciprocal field.
5. Whether the Raman overlap uses transpose or Hermitian conjugate.
6. How incident and scattered optical frequencies are chosen.
7. How field normalisation is defined.
8. How layer coordinates and phase origins are defined.
9. How absorbing media are represented, including the sign convention for complex permittivity.

These conventions should be mentioned in the tests, because many apparent errors in Raman strength are actually convention mismatches.

---

# Part M: Summary

The most robust verification strategy is to move from simple to complex:

```text
Raman tensor algebra
        ↓
Uniform single-crystal slab
        ↓
Thin-slab and thickness-scaling limits
        ↓
Scalar Fresnel/transfer-matrix comparison
        ↓
Layer subdivision and origin invariance
        ↓
Reciprocity checks
        ↓
Complex permittivity and complex Raman tensors
        ↓
Realistic layered Raman validation examples
```

The two strongest early regression tests are:

1. the rotated diagonal Raman tensor angular-dependence test;
2. the layer subdivision invariance test.

Together, these catch many of the most likely errors in a layered crystal Raman implementation.
