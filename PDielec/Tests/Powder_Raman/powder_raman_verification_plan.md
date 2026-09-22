# Powder Raman Verification Plan

## 1. Purpose

This document proposes a structured verification plan for a powder Raman implementation in which the Raman response of small crystallites is modified by macroscopic electromagnetic local-field effects, particle shape, matrix permittivity, and orientational averaging.

The aim is not to validate every physical approximation at once. Instead, the tests should be built as a ladder:

1. Verify the Raman tensor algebra without local-field effects.
2. Verify orientation averaging against analytic powder Raman formulae.
3. Verify local-field and depolarisation tensor construction.
4. Verify limiting cases for spheres, isotropic particles, and matched particle/matrix permittivity.
5. Verify anisotropic particles and ellipsoids.
6. Verify phonon-frequency or denominator corrections, if included.
7. Verify numerical orientation sampling.
8. Verify consistency with existing PDielec infrared powder machinery where appropriate.

The most important principle is to create artificial test systems with known behaviour before attempting real materials.

---

## 2. Main quantities to test

The implementation likely contains several conceptually distinct operations. Each should be tested separately before the full powder Raman spectrum is tested.

The important quantities are:

- the microscopic or crystal Raman tensor for each mode, `R_m`;
- the orientation matrix, `G`, mapping crystal axes to laboratory or particle axes;
- the particle-shape depolarisation tensor, `L`;
- the particle permittivity tensor, `epsilon_i`;
- the external matrix permittivity, `epsilon_e`;
- the local-field tensor for incident light;
- the local-field tensor for scattered light;
- the effective or particle Raman tensor;
- the orientationally averaged VV and VH Raman strengths;
- any frequency-dependent correction to mode frequencies or oscillator denominators.

A useful general form for the particle Raman tensor is of the type

```text
R_particle = local_field_scattered^T * R_crystal_or_particle * local_field_incident
```

or a related expression depending on the exact convention used in the code.

The tests below are written so that they can be adapted to the exact implementation convention.

---

## 3. Verification ladder

The tests should be implemented in the following order.

1. Raman tensor algebra without local fields.
2. Analytic powder averaging for simple Raman tensors.
3. Isotropic sphere local-field tests.
4. Permittivity-matched particle/matrix tests.
5. Anisotropic sphere tests.
6. Ellipsoid depolarisation tensor tests.
7. Orientation sampling convergence tests.
8. Frequency-shift or phonon-denominator tests.
9. Full artificial powder spectra.
10. Real-material regression tests.

Each stage should have automated numerical tests with clear tolerances.

---

# Part A — Raman tensor and powder averaging tests without local-field effects

These tests should be run with all particle local-field corrections disabled, or with a fictitious case where the local-field tensors reduce to the identity.

## A1. Isotropic Raman tensor

### Setup

Use a single Raman tensor

```text
R = a I
```

where `I` is the identity matrix.

### Expected behaviour

The Raman response should be independent of crystal orientation.

For crossed polarisations, the Raman amplitude should vanish in the simple single-crystal geometry.

For powder averaging, the depolarisation ratio should be zero for a purely isotropic Raman tensor.

### Checks

- Rotate the tensor by many random orientations.
- Confirm that the parallel intensity is constant.
- Confirm that the crossed intensity is zero, within numerical tolerance.
- Confirm that orientation sampling does not change the result.

### Suggested tolerance

For double precision, relative errors should be close to machine precision, for example `1e-12` to `1e-10`, depending on the sampling method.

---

## A2. Diagonal anisotropic Raman tensor with analytic invariants

### Setup

Use

```text
R = diag(a, b, c)
```

with real values, for example

```text
a = 2
b = 1
c = 0.5
```

### Expected behaviour

The isotropic invariant is

```text
alpha = (a + b + c) / 3
```

The anisotropic invariant may be written in the common Placzek form

```text
gamma^2 = 0.5 * [ (a-b)^2 + (b-c)^2 + (c-a)^2 ]
```

for a symmetric diagonal tensor.

For ordinary powder Raman in the Placzek approximation, the common unnormalised parallel and perpendicular strengths are proportional to

```text
I_parallel     ~ 45 alpha^2 + 4 gamma^2
I_perpendicular ~ 3 gamma^2
```

The depolarisation ratio is

```text
rho = I_perpendicular / I_parallel
```

### Checks

- Compare numerical orientation averaging with analytic invariant formulae.
- Test several diagonal tensors.
- Include special cases:
  - `a = b = c`, isotropic tensor;
  - `a = b != c`, uniaxial tensor;
  - `a`, `b`, `c` all different.

### Pass criterion

Numerical orientation averages should agree with analytic values within the expected sampling error.

For deterministic dense quadrature, target relative errors below `1e-6`.

For Monte Carlo sampling, verify statistical convergence as `1/sqrt(N)`.

---

## A3. Off-diagonal Raman tensor

### Setup

Use a purely off-diagonal symmetric tensor, for example

```text
R = [[0, d, 0],
     [d, 0, 0],
     [0, 0, 0]]
```

### Expected behaviour

The isotropic invariant should be zero.

The powder spectrum should be entirely anisotropic.

The depolarisation ratio should approach the standard value for a purely anisotropic symmetric tensor:

```text
rho = 3 / 4
```

when the usual Placzek formulae are used:

```text
I_parallel     ~ 4 gamma^2
I_perpendicular ~ 3 gamma^2
```

### Checks

- Confirm zero isotropic invariant.
- Confirm the expected VV and VH ratio.
- Confirm that rotating the tensor before powder averaging does not change the average.

---

## A4. Antisymmetric Raman tensor test, if supported

Some formulations ignore antisymmetric Raman tensor components, while others keep them for resonance or magnetic cases.

### Setup

Use

```text
R = [[ 0, d, 0],
     [-d, 0, 0],
     [ 0, 0, 0]]
```

### Expected behaviour

The expected result depends on whether antisymmetric Raman tensors are supported.

### Checks

- If antisymmetric components are intentionally excluded, confirm that the code rejects them or symmetrises them explicitly.
- If antisymmetric components are supported, test the result against the chosen invariant convention.
- Document the convention clearly.

---

# Part B — Local-field tensor tests for spherical particles

The next tests verify the particle local-field machinery, initially for the simplest case: a spherical particle in an isotropic matrix.

## B1. Isotropic sphere in isotropic matrix

### Setup

Use scalar particle permittivity `epsilon_i` and scalar matrix permittivity `epsilon_e`.

For a sphere, the depolarisation tensor is

```text
L = I / 3
```

The internal local-field factor for a uniform external field is usually

```text
N = 3 epsilon_e / (epsilon_i + 2 epsilon_e)
```

or an equivalent tensor expression depending on convention.

### Expected behaviour

For an isotropic sphere, the incident and scattered local-field factors are scalar.

Therefore the effective Raman tensor should be the original Raman tensor multiplied by scalar local-field factors:

```text
R_eff ~ N_s * R * N_i
```

For scalar local fields this means

```text
R_eff ~ N_s N_i R
```

and the Raman intensity should scale as

```text
|N_s N_i|^2
```

relative to the uncorrected result.

### Checks

- Use `epsilon_i = epsilon_e`; verify `N = 1` and the result equals the uncorrected powder Raman result.
- Use real positive permittivities and compare with the analytic scalar expression.
- Use complex permittivities and compare with the complex scalar expression.
- Confirm that multiplying the Raman tensor by a global complex phase does not change the intensity.

### Suggested cases

```text
epsilon_i = 4.0
epsilon_e = 1.0
N = 3 / 6 = 0.5
```

If the same local-field factor is used at incident and scattered frequencies, the intensity should scale by

```text
|0.5 * 0.5|^2 = 1/16
```

relative to the uncorrected value.

---

## B2. Matched permittivity limit

### Setup

Set

```text
epsilon_i = epsilon_e
```

for both incident and scattered frequencies.

### Expected behaviour

The particle should become electromagnetically invisible to the matrix.

The local-field correction should reduce to the identity.

### Checks

- The effective Raman tensor equals the bare Raman tensor.
- No phonon-frequency correction due to dielectric contrast should occur.
- Raman intensities match the no-local-field powder average.
- The result should be independent of particle shape if all dielectric contrast terms vanish.

This is one of the most important tests.

---

## B3. High-permittivity contrast limit

### Setup

Use a real isotropic particle with large permittivity in air, for example

```text
epsilon_i = 100
epsilon_e = 1
```

### Expected behaviour

For a sphere,

```text
N = 3 / 102
```

so the internal optical field should be strongly suppressed.

The Raman intensity should be strongly reduced by the local-field factors.

### Checks

- Confirm monotonic decrease of Raman strength as `epsilon_i / epsilon_e` increases.
- Confirm that no unphysical divergence occurs away from resonant denominators.
- Confirm that the result is finite and well behaved for large positive real permittivity.

---

## B4. Surface-mode resonance test for negative permittivity

### Setup

For an isotropic sphere in an isotropic matrix, the local-field denominator contains

```text
epsilon_i + 2 epsilon_e
```

A resonance occurs near

```text
epsilon_i = -2 epsilon_e
```

if losses are small.

### Expected behaviour

The local-field factor should increase strongly near the surface-mode condition.

With a finite imaginary part of the permittivity, the response should remain finite.

### Checks

- Sweep real `epsilon_i` through `-2 epsilon_e` with a small imaginary component.
- Confirm that the response peaks near the expected condition.
- Increase the imaginary part and confirm that the peak broadens and decreases.
- Confirm that the code does not produce NaNs or infinities unless damping is exactly zero.

This test is useful if the Raman theory includes frequency-dependent local fields at optical or phonon frequencies.

---

# Part C — Anisotropic spherical particle tests

A sphere with anisotropic permittivity has a simple depolarisation tensor `L = I/3`, but the local-field tensor is not generally a scalar.

## C1. Diagonal anisotropic permittivity aligned with Raman tensor

### Setup

Use

```text
epsilon_i = diag(eps_x, eps_y, eps_z)
epsilon_e = scalar
L = I / 3
R = diag(a, b, c)
```

The local-field tensor should be diagonal in the same basis.

For a spherical particle, the local-field factors along principal axes are expected to have the form

```text
N_j = 3 epsilon_e / (epsilon_j + 2 epsilon_e)
```

for `j = x, y, z`, depending on the exact convention.

### Expected behaviour

The effective Raman tensor remains diagonal:

```text
R_eff_jj ~ N_s,j * R_jj * N_i,j
```

### Checks

- Compare each diagonal component with the analytic expression.
- Confirm that off-diagonal components remain zero.
- Confirm that if `eps_x = eps_y = eps_z`, the isotropic sphere result is recovered.

---

## C2. Rotated anisotropic permittivity

### Setup

Use a diagonal permittivity in the crystal frame, then rotate it with the same orientation matrix as the Raman tensor:

```text
epsilon_lab = G epsilon_crystal G^T
R_lab       = G R_crystal G^T
```

### Expected behaviour

For a spherical particle, rotating the entire crystal should rotate both `epsilon` and `R` consistently.

The powder-averaged result should be independent of the initial choice of crystal axes.

### Checks

- Compute the powder average using the original crystal axes.
- Apply a fixed arbitrary rotation to both `epsilon` and all Raman tensors.
- Repeat the powder average.
- The final VV and VH intensities should be unchanged.

This catches inconsistent rotation of `epsilon`, `R`, Born charges, or local-field tensors.

---

## C3. Complex symmetric anisotropic permittivity

### Setup

Use a complex symmetric tensor, for example

```text
epsilon = [[eps_x, delta, 0],
           [delta, eps_y, 0],
           [0,     0,     eps_z]]
```

where `eps_x`, `eps_y`, `eps_z`, and `delta` may be complex.

### Expected behaviour

The local-field tensor should be computed using the full tensor algebra rather than assuming real diagonal eigenvectors.

### Checks

- Confirm that the code does not assume Hermitian diagonalisation.
- Confirm that direct tensor inversion gives the same result as any internal implementation.
- Confirm that the result varies smoothly as the off-diagonal element `delta` is increased from zero.
- Confirm that the result reduces to the diagonal case when `delta = 0`.

---

# Part D — Ellipsoidal particle and depolarisation tensor tests

These tests verify the particle-shape machinery.

## D1. Spherical limit of ellipsoid

### Setup

Use an ellipsoid with semi-axes

```text
a = b = c
```

### Expected behaviour

The depolarisation tensor should be

```text
L = I / 3
```

### Checks

- Confirm that all three depolarisation factors are `1/3`.
- Confirm that the ellipsoid code reproduces the spherical particle Raman result.
- Confirm independence of orientation for isotropic `epsilon_i` and isotropic `R`.

---

## D2. Sum rule for depolarisation factors

### Setup

Use several ellipsoids:

- prolate: `a > b = c`;
- oblate: `a = b > c`;
- triaxial: `a != b != c`.

### Expected behaviour

For an ellipsoid in electrostatics, the principal depolarisation factors should satisfy

```text
L_x + L_y + L_z = 1
```

### Checks

- Confirm the sum rule numerically.
- Confirm that each `L_j` lies between 0 and 1 for ordinary ellipsoids.
- Confirm the correct limiting trends:
  - needle-like axis: `L` along long axis tends to 0;
  - disc normal axis: `L` along short axis tends to 1.

---

## D3. Isotropic particle in ellipsoid

### Setup

Use scalar `epsilon_i`, scalar `epsilon_e`, and ellipsoidal `L = diag(L_x, L_y, L_z)`.

### Expected behaviour

The local-field tensor should be diagonal in the ellipsoid principal-axis frame, with axis-dependent local-field factors.

A common form is

```text
N_j = epsilon_e / [epsilon_e + L_j (epsilon_i - epsilon_e)]
```

which reduces to the sphere expression when `L_j = 1/3`.

### Checks

- Compare with the analytic axis-by-axis expression.
- Confirm recovery of the spherical result when `L_x = L_y = L_z = 1/3`.
- Confirm that shape affects the Raman intensity even when the material itself is isotropic, unless `R` is also isotropic and the final powder average removes the anisotropy appropriately.

---

## D4. Rotation of particle shape relative to crystal axes

### Setup

For non-spherical particles, distinguish between:

1. crystal orientation relative to particle axes;
2. particle orientation relative to laboratory axes.

Use a triaxial ellipsoid and an anisotropic Raman tensor.

### Expected behaviour

If the implementation assumes that the particle shape axes are locked to the crystal axes, then rotating the crystal should rotate `R`, `epsilon`, and `L` together.

If the implementation allows independent particle-shape and crystal orientations, then both rotations must be sampled or specified.

### Checks

- Verify the documented convention.
- Confirm that rotating all tensors together does not change the powder average.
- Confirm that rotating only the Raman tensor relative to the depolarisation tensor does change the result for anisotropic cases.
- Confirm that the result is unchanged for isotropic `epsilon_i`, isotropic `R`, or spherical `L`, where appropriate.

This is a key test for non-spherical powder particles.

---

# Part E — Orientation averaging tests

These tests verify the numerical sampling used for random powder orientations.

## E1. Uniform SO(3) sampling test

### Purpose

A random powder average requires uniform sampling over 3D rotations, not merely over a polar angle or a solid angle direction.

### Setup

Generate orientations using the implementation's chosen method:

- Euler-angle quadrature;
- random quaternions;
- Sobol sequence mapped to quaternions;
- Lebedev directions plus rotation about the direction;
- other quadrature.

### Expected behaviour

For a uniform orientation average,

```text
<g_ij^2> = 1/3
```

for each direction-cosine element.

Also,

```text
<g_ij> = 0
```

for all elements.

### Checks

- Compute averages of direction cosines over the orientation set.
- Confirm they agree with isotropic expectations.
- Confirm that sampling over only two angles fails this test for quantities sensitive to rotation about the particle axis.

---

## E2. Convergence against analytic powder invariants

### Setup

Use the diagonal and off-diagonal Raman tensors from Part A.

### Expected behaviour

The numerical average should converge to the analytic powder-invariant result.

### Checks

- Run with increasing number of orientations `N`.
- Plot or record the error against `N`.
- For Monte Carlo, confirm approximately `1/sqrt(N)` convergence.
- For Sobol or deterministic sampling, confirm smoother and faster convergence.

### Suggested regression values

Use several fixed `N` values, for example:

```text
N = 16, 64, 256, 1024, 4096
```

Keep the same sequence deterministic so that tests are reproducible.

---

## E3. Rotation-invariance of powder average

### Setup

Choose an arbitrary fixed rotation `G0`.

Compute the powder average for a tensor set `{R_m, epsilon}`.

Then rotate every tensor by `G0` before performing the powder average:

```text
R_m'       = G0 R_m G0^T
epsilon'   = G0 epsilon G0^T
```

### Expected behaviour

The final powder-averaged VV and VH intensities should be unchanged.

### Checks

- Test for isotropic spheres.
- Test for anisotropic spheres.
- Test for ellipsoids where the ellipsoid axes are also rotated consistently.

This is one of the most important whole-implementation tests.

---

## E4. VV and VH polarisation definitions

### Setup

Use a powder with ordinary backscattering or right-angle scattering geometry, depending on the implementation.

Define incident polarisation and analysed scattered polarisations explicitly.

### Expected behaviour

For the standard Placzek powder result, the VV and VH intensities should reduce to the known invariant expressions when local fields are disabled.

### Checks

- Confirm that VV corresponds to parallel incident and analysed scattered polarisations.
- Confirm that VH corresponds to perpendicular analysed scattered polarisation.
- Confirm that unpolarised or averaged detection, if implemented, equals the intended combination of VV and VH.
- Confirm that the depolarisation ratio is correct for isotropic and purely anisotropic tensors.

---

# Part F — Frequency and phonon correction tests

These tests apply if the implementation includes frequency shifts or modifications to the dynamical matrix due to particle shape, Born charges, and dielectric boundary conditions.

## F1. No-contrast frequency correction test

### Setup

Set the matrix and particle optical/background permittivities equal:

```text
epsilon_i_infinity = epsilon_e
```

or the equivalent no-contrast condition in the implemented theory.

### Expected behaviour

The particle-induced correction to phonon frequencies should vanish.

The mode frequencies should remain equal to the input TO or reference frequencies.

### Checks

- Confirm zero or numerically negligible frequency shift.
- Confirm that Raman strengths still reduce to the ordinary local-field-free powder result when optical local-field corrections are also disabled or become identity.

---

## F2. Spherical isotropic polar mode test

### Setup

Use a simple artificial polar mode with scalar Born charge and scalar high-frequency permittivity.

For a sphere, the surface-mode condition is often associated with a denominator involving

```text
epsilon_i(omega) + 2 epsilon_e
```

or an equivalent depolarisation expression.

### Expected behaviour

The particle-modified polar-mode frequency should lie between the TO and LO frequencies for a simple polar material in a non-absorbing matrix, depending on the dielectric contrast.

For a sphere in vacuum, the surface optical mode is usually associated with

```text
epsilon_i(omega) = -2
```

for scalar permittivity.

### Checks

- Use a one-mode Lorentz dielectric function with known TO and LO frequencies.
- Solve the analytic surface-mode condition independently.
- Compare the code's frequency correction or spectral peak position with the analytic solution.
- Confirm that the result moves correctly when `epsilon_e` is changed.

---

## F3. Shape dependence of polar-mode correction

### Setup

Use an isotropic polar material in ellipsoids with different depolarisation factors.

For a principal axis `j`, a resonance condition commonly involves

```text
epsilon_e + L_j (epsilon_i(omega) - epsilon_e) = 0
```

or equivalently

```text
epsilon_i(omega) = -[(1 - L_j) / L_j] epsilon_e
```

for scalar permittivities.

### Expected behaviour

Changing particle shape should split or shift the surface-mode frequencies.

### Checks

- Compare sphere, prolate, and oblate particles.
- Confirm that the resonance condition follows the depolarisation factor.
- Confirm that the spherical result is recovered when `L_j = 1/3`.
- Confirm that in the weakly polar or zero-Born-charge limit the correction disappears.

---

## F4. Zero Born charge limit

### Setup

Set the mass-weighted Born charge tensor for a mode to zero.

### Expected behaviour

Any electrostatic phonon-frequency correction associated with long-range polar fields should vanish.

The Raman tensor may still be non-zero if the mode is Raman active.

### Checks

- Confirm no particle-induced frequency shift.
- Confirm Raman intensity remains if `R_m` is non-zero.
- Confirm infrared oscillator strength vanishes if the same mode is used in IR tests.

This test separates Raman activity from polar IR activity.

---

## F5. Degenerate-mode handling

### Setup

Use two or three modes with identical reference frequencies.

Assign Raman tensors and Born charge tensors that cause coupling under the particle correction.

### Expected behaviour

If the implemented frequency correction creates a mode-coupling matrix, degenerate modes may mix.

The final observable spectrum should be invariant under arbitrary rotations within the degenerate subspace.

### Checks

- Rotate the basis of degenerate modes with an arbitrary orthogonal matrix.
- Transform Raman tensors and Born charges consistently.
- Confirm that the total spectrum is unchanged.
- Confirm that individual mode labels may change, but summed intensity is invariant.

This is an important test if the implementation diagonalises a corrected dynamical matrix.

---

# Part G — Effective medium consistency tests

These tests are useful because PDielec already contains effective-medium treatments for infrared powder absorption.

## G1. IR absorption consistency for scalar sphere

### Setup

Use the same dielectric function and particle shape in the existing powder IR machinery and the new Raman local-field machinery.

### Expected behaviour

The local-field denominators controlling particle resonances should be consistent between IR and Raman calculations.

Raman and IR intensities need not be proportional, but resonance positions arising from the same electrostatic boundary condition should be consistent.

### Checks

- Use a simple one-oscillator isotropic material.
- Compute the IR absorption spectrum for small spheres.
- Compute the Raman local-field enhancement factor or Raman spectrum for the same particle.
- Confirm that local-field resonances occur at compatible frequencies.

---

## G2. Maxwell-Garnett dilute limit

### Setup

Use a very small volume fraction of spherical particles in a matrix.

### Expected behaviour

In the dilute limit, the effective-medium absorption should be dominated by the single-particle polarizability denominator.

The Raman local-field factor for an isolated particle should use the same single-particle denominator.

### Checks

- Compare resonance positions, not absolute intensities.
- Confirm convergence as volume fraction tends to zero.
- Confirm that increasing volume fraction changes the EMT result but not the isolated-particle Raman result, unless the Raman implementation explicitly includes interparticle interactions.

---

## G3. Bruggeman vs isolated-particle distinction

### Purpose

Avoid accidentally using an effective-medium result where an isolated-particle local-field factor is intended.

### Checks

- Run the same material with isolated-particle local fields and with EMT-derived quantities, if both options exist.
- Confirm that the two approaches agree only in appropriate limits.
- Document which model is used for the Raman calculation.

---

# Part H — Full artificial spectrum tests

These tests combine frequencies, line shapes, local fields, and orientation averaging.

## H1. Non-polar Raman-only material

### Setup

Create an artificial material with:

- several Raman-active modes;
- zero Born charges;
- real constant optical permittivity;
- no phonon-frequency correction.

### Expected behaviour

The powder Raman spectrum should be a sum of lines at the input frequencies.

Particle shape and matrix permittivity may affect intensities through optical local-field factors, but should not shift mode frequencies.

### Checks

- Confirm peak positions equal input frequencies.
- Confirm integrated intensities match analytic powder averages when local fields are disabled.
- Confirm local-field scaling for isotropic spheres.

---

## H2. Single polar Raman-active mode

### Setup

Use one polar mode with:

- non-zero Raman tensor;
- non-zero Born charge;
- Lorentzian dielectric response;
- spherical particle in isotropic matrix.

### Expected behaviour

Depending on the implemented theory, the Raman peak may occur at a particle-modified surface-mode frequency rather than the bulk TO frequency.

### Checks

- Confirm the no-contrast limit gives the reference frequency.
- Confirm the finite-contrast result agrees with the analytic scalar surface-mode condition.
- Confirm the shift changes with matrix permittivity.
- Confirm increasing damping broadens the feature but does not move it unphysically.

---

## H3. Multiple modes with overlapping lines

### Setup

Use two modes close in frequency with different Raman tensors and Born charges.

### Expected behaviour

The spectrum should be the sum of the mode contributions unless the theory explicitly includes mode mixing.

### Checks

- Confirm integrated area is approximately additive when lines overlap.
- Confirm changing one mode's Raman tensor only affects that mode unless mode coupling is active.
- Confirm that sorting modes by frequency does not change the final spectrum.

---

# Part I — Numerical stability and regression tests

## I1. Matrix inversion stability

### Setup

Choose permittivities close to local-field resonances, with small damping.

### Checks

- Confirm that matrix inversions are stable.
- Confirm useful warnings or graceful behaviour near singular denominators.
- Confirm no NaNs or infinities for finite damping.
- Confirm results converge as damping is varied.

---

## I2. Orientation-number convergence

### Setup

Run the same powder calculation with increasing orientation sample count.

### Checks

- Record VV and VH intensities.
- Record depolarisation ratios.
- Confirm monotonic or statistically sensible convergence.
- Store reference values for a fixed deterministic sampling grid.

---

## I3. Frequency-grid convergence

### Setup

Compute a spectrum using different frequency grid spacings and line widths.

### Checks

- Integrated intensity should be stable with grid refinement.
- Peak height should scale sensibly with line width.
- Peak position should not depend on grid spacing beyond interpolation error.

---

## I4. Invariance to mode ordering

### Setup

Shuffle the list of modes before computing the spectrum.

### Expected behaviour

The final spectrum should be unchanged.

### Checks

- Compare spectra point-by-point.
- Compare integrated intensities by mode, where available.

---

## I5. Units and scaling tests

### Setup

Scale Raman tensors, mode volumes, or normal-coordinate conventions in controlled ways.

### Expected behaviour

If

```text
R_m -> lambda R_m
```

then Raman intensity should scale as

```text
|lambda|^2
```

provided no other quantities are changed.

### Checks

- Scale `R_m` by 2 and confirm intensity scales by 4.
- Scale all Raman tensors by a complex phase and confirm intensity is unchanged.
- Confirm Bose population factors, if included, change the spectrum only through the intended frequency- and temperature-dependent prefactor.
- Confirm Stokes and anti-Stokes options, if present, use the correct thermal factors.

---

# Part J — Suggested automated test cases

The following tests are good candidates for a continuous integration test suite.

## J1. `test_isotropic_tensor_powder_average`

- Input: `R = a I`.
- Expected: VH intensity zero; VV independent of orientation.
- Purpose: basic Raman tensor and polarisation check.

## J2. `test_diagonal_tensor_powder_invariants`

- Input: `R = diag(2, 1, 0.5)`.
- Expected: numerical VV/VH match analytic Placzek formulae.
- Purpose: orientation averaging and invariant check.

## J3. `test_offdiagonal_tensor_depolarisation_ratio`

- Input: symmetric off-diagonal tensor.
- Expected: depolarisation ratio `rho = 0.75` for the usual non-resonant Placzek convention.
- Purpose: crossed-polarisation and anisotropic invariant check.

## J4. `test_sphere_matched_permittivity`

- Input: `epsilon_i = epsilon_e`.
- Expected: local-field tensor is identity; Raman result equals uncorrected result.
- Purpose: local-field denominator check.

## J5. `test_isotropic_sphere_scalar_local_field`

- Input: scalar `epsilon_i`, scalar `epsilon_e`, sphere.
- Expected: intensity scaled by `|N_i N_s|^2`.
- Purpose: analytic local-field benchmark.

## J6. `test_anisotropic_sphere_diagonal_local_field`

- Input: diagonal anisotropic `epsilon_i`, sphere.
- Expected: axis-by-axis analytic local-field factors.
- Purpose: tensor local-field test.

## J7. `test_ellipsoid_depolarisation_sum_rule`

- Input: several ellipsoid aspect ratios.
- Expected: `L_x + L_y + L_z = 1`.
- Purpose: shape-factor validation.

## J8. `test_layer_or_particle_rotation_invariance`

- Input: rotate all tensors by the same arbitrary rotation before powder averaging.
- Expected: final VV/VH results unchanged.
- Purpose: frame-convention validation.

## J9. `test_orientation_sampling_convergence`

- Input: anisotropic tensor with known analytic powder average.
- Expected: convergence to analytic value with increasing sample count.
- Purpose: sampling validation.

## J10. `test_zero_born_charge_no_frequency_shift`

- Input: Raman-active mode with zero Born charge.
- Expected: no particle-induced polar frequency shift.
- Purpose: separate Raman activity from polar electrostatic correction.

## J11. `test_isotropic_surface_mode_condition`

- Input: one Lorentz oscillator, spherical particle.
- Expected: resonance near `epsilon_i(omega) = -2 epsilon_e`.
- Purpose: polar frequency/local-field denominator validation.

## J12. `test_mode_ordering_invariance`

- Input: same set of modes in different orders.
- Expected: identical total spectrum.
- Purpose: spectrum assembly validation.

---

# Part K — Recommended reference artificial systems

## K1. Simple non-polar Raman material

Use three modes:

```text
mode 1 frequency = 100 cm^-1
R1 = diag(1, 1, 1)
Z1 = 0

mode 2 frequency = 200 cm^-1
R2 = diag(2, 1, 0.5)
Z2 = 0

mode 3 frequency = 300 cm^-1
R3 = [[0, 1, 0],
      [1, 0, 0],
      [0, 0, 0]]
Z3 = 0
```

Expected behaviour:

- mode 1 has zero depolarisation ratio;
- mode 2 has intermediate depolarisation ratio;
- mode 3 has depolarisation ratio near 0.75;
- no mode frequency shifts occur if Born charges are zero.

This is the best first full-spectrum test.

---

## K2. Isotropic polar sphere

Use one mode with a scalar Lorentz dielectric function:

```text
epsilon(omega) = epsilon_infinity * (omega_LO^2 - omega^2 - i gamma omega)
                 / (omega_TO^2 - omega^2 - i gamma omega)
```

For a sphere in a matrix, the surface-mode condition is approximately

```text
epsilon(omega) + 2 epsilon_e = 0
```

Expected behaviour:

- no damping gives a sharp resonance;
- finite damping gives a finite peak;
- increasing `epsilon_e` shifts the resonance;
- matched contrast should remove the particle-induced correction, depending on the exact implemented formulation.

---

## K3. Anisotropic Raman tensor with isotropic optical permittivity

Use anisotropic `R` but scalar `epsilon_i`.

Expected behaviour:

- local-field corrections scale all tensor components by the same factor;
- depolarisation ratio should be unchanged by scalar local-field factors;
- absolute intensities change by `|N_i N_s|^2`.

This test separates Raman anisotropy from dielectric anisotropy.

---

## K4. Isotropic Raman tensor with anisotropic optical permittivity

Use `R = a I` but anisotropic `epsilon_i`.

Expected behaviour:

- local-field tensors can make the effective particle Raman tensor anisotropic;
- powder depolarisation ratio may become non-zero if the local-field correction is anisotropic;
- if `epsilon_i` becomes isotropic, the depolarisation ratio returns to zero.

This is a sensitive test of local-field tensor placement.

---

# Part L — Practical validation workflow

A good practical sequence is:

1. Run all no-local-field Raman tensor tests.
2. Verify analytic powder invariants.
3. Turn on scalar local fields for isotropic spheres.
4. Verify matched-permittivity and scalar local-field scaling.
5. Turn on anisotropic permittivity for spherical particles.
6. Verify full tensor local-field behaviour.
7. Turn on ellipsoidal particle shapes.
8. Verify depolarisation factors and shape limits.
9. Add frequency-dependent dielectric functions.
10. Verify surface-mode resonance conditions.
11. Add Born-charge-based frequency corrections, if implemented.
12. Verify no-contrast and zero-Born-charge limits.
13. Run full artificial spectra.
14. Compare selected real-material cases with existing PDielec IR behaviour and with known qualitative Raman selection rules.

---

# Part M — Most important tests to implement first

If only a small number of tests can be implemented initially, prioritise these:

1. `R = a I` gives zero VH and orientation-independent VV.
2. `R = diag(2, 1, 0.5)` matches analytic powder invariants.
3. Symmetric off-diagonal `R_xy` gives depolarisation ratio near 0.75.
4. `epsilon_i = epsilon_e` gives identity local-field correction.
5. Isotropic sphere gives scalar local-field scaling `|N_i N_s|^2`.
6. Rotating all tensors before powder averaging leaves the result unchanged.
7. Splitting or reordering modes does not change the total spectrum.
8. Zero Born charge gives no polar frequency shift.
9. Isotropic polar sphere resonance occurs near `epsilon_i = -2 epsilon_e`.
10. Orientation sampling converges to the analytic powder average.

These tests together would give strong confidence that the powder Raman implementation is internally consistent.

---

# Part N — Common failure modes these tests should catch

## N1. Wrong rotation convention

Symptoms:

- angular dependence has the wrong phase;
- powder average changes after a global rotation;
- off-diagonal tensors give incorrect VV/VH ratios.

Likely cause:

```text
R_lab = G^T R G
```

used where

```text
R_lab = G R G^T
```

was intended, or vice versa.

---

## N2. Sampling only directions, not full orientations

Symptoms:

- spherical particles look correct;
- ellipsoids or anisotropic tensors give biased averages;
- results depend on arbitrary choice of rotation about the beam direction.

Likely cause:

Only two angular variables are sampled instead of full SO(3).

---

## N3. Incorrect local-field tensor placement

Symptoms:

- scalar tests pass, but anisotropic local-field tests fail;
- `R = a I` with anisotropic permittivity gives the wrong effective anisotropy;
- VV/VH ratios change incorrectly.

Likely cause:

Using

```text
N R N
```

instead of the intended transpose/order convention, or applying incident and scattered local-field tensors in the wrong order.

---

## N4. Incorrect complex conjugation

Symptoms:

- real permittivity tests pass;
- absorbing materials fail;
- results change incorrectly under global complex phase of `R`.

Likely cause:

Using Hermitian conjugation where the theory requires transpose, or vice versa.

---

## N5. Confusing optical local-field correction with phonon-frequency correction

Symptoms:

- non-polar modes shift frequency;
- matched permittivity still gives a frequency correction;
- zero Born charge modes shift.

Likely cause:

Applying the dielectric-boundary correction to all modes regardless of polar coupling.

---

## N6. EMT denominator used inconsistently

Symptoms:

- Raman resonances do not match isolated-particle IR resonances in the dilute limit;
- resonance positions depend incorrectly on volume fraction;
- isolated-particle and effective-medium models are mixed.

Likely cause:

Using Maxwell-Garnett or Bruggeman effective permittivity where a single-particle local-field tensor should be used.

---

# Part O — Suggested output diagnostics

For debugging and validation, the code should optionally print or return intermediate quantities:

- rotated Raman tensor for each sampled orientation;
- rotated permittivity tensor;
- depolarisation tensor in particle and lab frames;
- incident local-field tensor;
- scattered local-field tensor;
- effective particle Raman tensor;
- orientation-resolved VV and VH contributions;
- accumulated powder averages;
- mode frequency before and after correction;
- local-field denominator eigenvalues or condition numbers;
- Bose and frequency prefactors separately from tensor strengths.

These diagnostics make it much easier to locate mistakes when a full spectrum looks wrong.

---

# Part P — Final acceptance criteria

The powder Raman implementation can be considered well verified when the following are true:

1. Analytic powder averages are reproduced when local fields are disabled.
2. Isotropic sphere local-field factors reproduce scalar analytic results.
3. Matched particle/matrix permittivity removes local-field and frequency-shift effects.
4. Anisotropic tensor rotations are frame invariant after powder averaging.
5. Ellipsoid depolarisation factors obey known sum rules and limiting behaviour.
6. Orientation sampling converges reliably and reproducibly.
7. Zero Born charge modes do not receive polar frequency corrections.
8. Surface-mode or particle-mode frequencies agree with scalar analytic limits.
9. Complex permittivity tests are stable and obey phase-invariance checks.
10. Full artificial spectra are reproducible under mode reordering and orientation-grid refinement.

Passing these tests would give strong evidence that the implementation is correct at the level of tensor algebra, electrostatic local-field theory, orientation averaging, and spectrum assembly.
