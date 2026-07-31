# Modal-Pairs Final-State Design Note

## Motivation

The current modal-pairs implementation has accumulated several corrective rules:
NAC branch remapping, non-polar bypass, q-direction grouping, degenerate optical
q-channel grouping, and optional zero-q suppression.  Each rule addresses a real
physical issue, but the implementation still treats optical Berreman pairs as the
primary object and then decides how to combine them afterwards.

A more general implementation should start from the physical scattering channel
and only introduce modal-pair q labels when they correspond to physically
resolved phonon final states.  The internal Berreman decomposition is a useful
way to analyse phase matching, but it is not automatically a list of independent
Raman final states.

The central question is therefore:

```text
Which products of internal optical field components represent distinct phonon
final states, and which are only components of the same external-channel Raman
amplitude?
```

The design should answer that question explicitly rather than hiding it behind
geometry-dependent switches.

## Physical Starting Point

The experiment prepares an external incident photon channel and detects an
external outgoing photon channel.  The internal fields in a layered film contain
forward and backward Berreman components because Maxwell boundary conditions
produce reflections and standing-wave structure.  Those internal components are
not, by themselves, separately prepared or detected photon states.

For a finite coherent layer, the safest physical object is the external-channel
Raman amplitude:

```text
A = integral E_S_rec(z)^T R(q, branch) E_L(z) dz
```

where `E_L` and `E_S_rec` are the full internal fields associated with the
selected external incident and detected channels.  `E_S_rec` is the reciprocal detector
field: the field generated inside the stack by launching the detected analyser
channel backwards from the detector side at the Stokes frequency.  A modal
decomposition may be used to identify phase-matched parts of this amplitude,
but an incoherent sum over modal pairs requires an additional physical argument:
the phonon momentum/branch labels must be resolved final-state quantum numbers.

In the bulk or thick-layer limit, translational coherence can make phonon
momentum a good final-state label.  Then terms with different phonon wavevector
can become distinguishable and may be summed incoherently.  In a finite coherent
film, broken translational symmetry and boundary reflections relax strict
momentum conservation, so the coherent external-channel amplitude may remain the
more physical description.

### Reciprocal-Field Convention

The Raman source-overlap convention should use Lorentz reciprocity, not a
Hermitian inner product.  For time-harmonic complex amplitudes with reciprocal
non-magnetic media (`epsilon = epsilon^T`, `mu = mu^T`), the unconjugated
Lorentz reciprocity theorem gives

```text
integral_V (E_1 . J_2 - E_2 . J_1) dV
  = surface term.
```

Here "." is a bilinear dot product; it does not complex-conjugate either field.
Let problem 1 be the Raman-polarisation source inside the layer,

```text
J_R(z) = -i omega_S P_R(z),
P_R(z) = R_m(z) E_L(z),
```

and let problem 2 be the reciprocal detector problem: launch the detected
external analyser channel backwards from the detector side with the same field
normalisation used for optical collection.  Lorentz reciprocity then relates
the outgoing amplitude in the detected channel to

```text
A_m proportional to integral E_S_rec(z)^T R_m(z) E_L(z) dz.
```

The transpose is therefore ordinary transpose, not Hermitian conjugate.  This
is also the convention currently used by `LayeredRamanCalculator` and by the
Crystal Raman tensor-algebra and reciprocity tests.

Consequences:

- `E_S_rec` is not the complex conjugate of the physical scattered field.  It is
  the reciprocal solution of Maxwell's equations for the detector channel.
- Absorption is already contained in the complex reciprocal field returned by
  the optical solver.  Applying an additional complex conjugation would solve a
  different adjoint problem and generally breaks the bilinear reciprocity
  relation used here.
- The only final complex conjugation is the ordinary intensity operation
  `I = |A_m|^2` after the amplitude has been formed.
- For lossless symmetric test cases, the convention implies
  `A(p -> s) = A(s -> p)` for symmetric Raman tensors, because
  `E_s^T R E_p = E_p^T R^T E_s`.

The refactor should keep this convention unless a separate adjoint-GTM
derivation is introduced.  If such an adjoint formulation is ever added, it
must be a different named convention and must pass the same reciprocity tests.

### Berreman Pair Classification

When the full internal fields are decomposed into Berreman components, each
forward/backward Berreman pair `(i, j)` carries a net phonon wavevector
`q_ph(i,j) = k_L(i) - k_S(j)`.  For a backscattering geometry these pair types
are:

```text
(forward, forward)   -> q_ph ~ +2k  (LO-like)
(backward, backward) -> q_ph ~ -2k  (LO-like)
(forward, backward)  -> q_ph ~  0   (TO-like)
(backward, forward)  -> q_ph ~  0   (TO-like)
```

The cross-pairs `(forward, backward)` and `(backward, forward)` are
interference terms in the coherent external-channel amplitude — they are not
additional phonon final states.  Whether they need to be retained or may be
discarded depends on the coherence regime:

- In the **coherent finite film** regime, all four pair types contribute to the
  same external-channel amplitude and must be summed coherently before squaring.
  Discarding the cross-pairs would change the result; their contribution is a
  genuine part of the Raman amplitude, not an artefact.

- In the **incoherent depth** regime (thick layers integrated locally), phase
  averaging causes cross-pair amplitudes to cancel when integrated over depth.
  Only co-propagating pairs `(forward, forward)` and `(backward, backward)`
  survive the average.  Cross-pairs may therefore be discarded after the average.

- In the **bulk / phase-matched** limit, only the pairs whose phonon wavevector
  matches the externally imposed momentum transfer contribute.  For
  backscattering, this selects LO-like pairs; cross-pairs and the opposite-sign
  pair are suppressed by phase mismatch.

This classification is the physical rule that replaces all current geometry
switches.  The coherence regime is an explicit input; the pair classification
follows automatically.

### Wavevector Sign Convention

The pair classification above uses `k_S` as the physical scattered photon
wavevector, not necessarily the propagation direction of the reciprocal field
used to compute `E_S`.  This distinction must be made explicit in code.

For every optical component, the resolver should carry both:

```text
k_physical : photon wavevector in the actual scattering process
k_adjoint  : wavevector used in the reciprocal/adjoint field calculation
```

The phonon momentum transfer must use `k_physical`:

```text
q_ph = k_L_physical - k_S_physical
```

Changing this convention flips signs and can turn a nominal cross-pair into a
co-propagating pair.  The sign convention should therefore be covered by unit
tests for forward and backscattering at normal incidence.

## Target Formulation

For diagnostic or bulk-limit modal analysis, decompose the full external-channel
fields into internal optical components.  For a Raman-active layer, an incident
component `i`, a detected component `j`, and an original TO mode `m`, compute

```text
q_ph(i,j) = k_L(i) - k_S(j)
```

Then resolve `(m, q_ph, coherence regime, external channel)` into either:

- a physical phonon final state; or
- a coherent contribution to an already existing external-channel amplitude; or
- no separate contribution, if the component is only an artefact of an
  over-resolved internal decomposition.

When a separate phonon final state is justified, represent it as:

```text
FinalState(
    key,
    frequency_cm1,
    raman_tensor,
    linewidth_cm1,
    selected,
)
```

The intensity is then accumulated by physically justified final state, not
directly by optical pair:

```text
A_F = sum_{layers and optical pairs assigned to F} integral E_S_rec,j^T R_F E_L,i dz
I   = sum_F |A_F|^2
```

Distinct resolved final states are summed incoherently.  Contributions that are
not physically distinguishable final states are summed coherently before
squaring.

### Raman Tensor Assignment for Final States

For non-polar modes, `R_F` is q-independent and is the same for every pair
assigned to `F`.

For polar modes, `R_F` must be the NAC-modified Raman tensor evaluated at the
q direction that defines the final-state key.  All pairs assigned to the same
final state `F` use the same `R_F`, computed once at the q direction encoded in
`q_class_key`.  This is physically consistent provided the angular spread of
q directions within a single `q_class_key` bin is smaller than the angular
tolerance used to define the bin (see section on q_class_key below).

If two pairs assigned to the same final state have slightly different q
directions (e.g. at oblique incidence in a biaxial medium), the spread in NAC
frequency and tensor is a higher-order correction.  The implementation should
document this approximation rather than silently carrying it.  A diagnostic mode
may expose the per-pair tensor variation for validation.

### Polar Tensors in Coherent Finite Films and the Bulk Limit

The coherent finite-film regime is the most delicate case for polar modes.
Co-propagating internal pair products can be LO-like, while cross-pair products
can be near-zero-q and therefore TO-like.  If all of these products are coherent
components of one external-channel amplitude, a single Raman tensor `R_F` does
not describe the whole amplitude.

The coherent finite-film polar amplitude should therefore be written as a
coherent sum of pair-dependent NAC tensors:

```text
A_m = sum_{i,j} integral E_S_rec,j(z)^T R_m(q_ph(i,j)) E_L,i(z) dz
```

where `R_m(q_ph)` means:

- the TO Raman tensor for non-polar modes;
- the q-dependent NAC/EO Raman tensor for polar modes with finite `q_ph`;
- the TO-limit tensor for polar modes when `|q_ph|` is genuinely near zero.

This is not an incoherent modal-pair final-state sum.  It is one coherent
external-channel amplitude expanded in internal optical components.  The
ordinary intensity operation is applied only after the coherent pair sum:

```text
I_m = |A_m|^2.
```

For a laterally uniform layer, each pair contribution contains a longitudinal
phase factor of the form

```text
exp(i Delta k_ij z),  Delta k_ij = k_L,z(i) - k_S,z(j) - q_z
```

or the equivalent GTM sign convention.  Integrating over a thickness `L` gives
the finite-thickness phase-matching factor

```text
F_ij(L) = integral_0^L exp(i Delta k_ij z) dz
        = L exp(i Delta k_ij L/2) sinc(Delta k_ij L/2).
```

The bulk limit follows from this factor:

- If `Delta k_ij != 0`, `F_ij(L)/L -> 0` as `L -> infinity`; the pair is
  suppressed by phase mismatch.
- If `Delta k_ij = 0`, `F_ij(L)/L -> 1`; the pair survives as the bulk
  phase-matched contribution.

Therefore the coherent-film expression has the correct bulk limits without an
ad hoc switch:

- In external forward scattering, `q_ext -> 0`, so the surviving polar branch is
  TO-like.
- In external backscattering along the polar axis, `q_ext` is finite and the
  surviving polar branch is LO-like.
- Internal cross-pairs that produce near-zero `q_ph` in a backscattering field
  can contribute to a finite coherent film.  Under incoherent-depth
  integration they are local-field components of the same externally selected
  final state, rather than separate phonon final states.

The final-state resolver can still classify pair products, but for
`coherent_film` polar modes it should return coherent pair-resolved
contributions rather than forcing every pair into a resolved final-state bucket.
The `bulk_phase_matched` regime is then the asymptotic limit of the same
expression, not a different physical model.

## Final-State Resolver

The final-state resolver should be a separate, testable component.  It should
accept only the physical inputs needed to classify the phonon:

- original TO mode index;
- phonon wavevector `q_ph` in the lab frame;
- crystal-to-lab rotation for the layer;
- NAC ingredients and Raman tensors;
- mode-selection mask;
- external incident and detected scattering channels;
- coherence regime, such as coherent finite film, incoherent depth, or
  bulk/phase-matched limit.

It should return a classification, not just a frequency:

```text
ResolvedFinalState(FinalState)
CoherentExternalChannelContribution(label)
DiscardedInternalComponent(reason)
```

The resolver is where the physical criterion belongs.  It should not contain
ad hoc rules such as "backscatter suppresses q=0"; it should decide whether a
q-labelled phonon is actually a resolved final state for the chosen coherence
regime and external scattering channel.

## Coherence Regime Determination

The coherence regime is the primary input that controls how Berreman pairs map to
final states.  It must be an explicit user-visible setting, not inferred from
geometry heuristics.  Three regimes are sufficient for the planned use cases.

### coherent_film

All Berreman pairs within a layer are coherent with one another.  Cross-pairs
are coherent interference terms and must not be classified as independent final
states.  The resolver returns `CoherentExternalChannelContribution` for every
cross-pair product; the accumulator sums all contributions to a single
external-channel amplitude before squaring.

This is the physically correct regime for a thin film whose thickness is
comparable to the optical coherence length.  It is also the safest default when
the coherence length is unknown.

### incoherent_depth

The layer is divided into thin depth slices.  Within each slice the field
products are treated as in `coherent_film`.  Across slices, amplitudes are
summed as intensities (squared amplitudes), not as amplitudes.  Cross-pair
contributions at each slice carry a rapid phase factor `exp(i q_cross z)` that
averages to zero over the slice thickness when that thickness is large compared
to the coherence length `2pi / |q_cross|`.  The resolver may discard cross-pairs
whose phase-averaging factor is below a configurable threshold.

This regime is appropriate for bulk or thick-layer calculations where the phonon
coherence length is much shorter than the layer thickness.

### bulk_phase_matched

Only pairs whose phonon wavevector direction matches the externally imposed
momentum transfer `q_ext = k_inc - k_det` (evaluated from the true external
angles, not internal Berreman wavevectors) are retained.  All other pairs are
classified as `DiscardedInternalComponent`.  This regime recovers the standard
Porto selection rules in the bulk limit.

For finite thickness this selection should preferably emerge from an explicit
phase-matching/dephasing factor, not from an immediate hard discard.  A useful
intermediate form is:

```text
W_ij = sinc(Delta q_ij L / 2)^2
```

or the corresponding complex amplitude factor before squaring, where `Delta q`
measures mismatch from the externally imposed phonon momentum and `L` is the
coherent depth.  The hard `bulk_phase_matched` filter is the asymptotic
large-`L` limit of this weighting.  The implementation should keep this
distinction clear so that finite-thickness convergence can be tested rather
than imposed.

The external momentum transfer `q_ext` is computed by the `OpticalChannelResolver`
from the macroscopic scattering geometry and is passed to the resolver as an
input alongside the per-pair `q_ph`.

### Determining the Regime in Practice

The coherence regime should be selected by the user in the GUI or script,
defaulting to `coherent_film` for thin-film stacks and `incoherent_depth` for
bulk-like powder layers.  An automatic heuristic based on layer thickness vs
phonon coherence length may be provided as a convenience, but it must be
overridable and its threshold must be documented.

## Final-State Key

The key should encode physical indistinguishability, not implementation details.
A suitable starting point is:

```text
(
    layer_coherence_group,
    detected_channel_key,
    phonon_branch_key,
    frequency_key,
    q_class_key,
)
```

The details should follow these rules.

### Non-Polar Modes

Modes with negligible oscillator strength are q-independent.  Their final-state
key must not include `q_ph` or `q_hat`.  All optical paths that create the same
non-polar mode in the same coherence group should be summed coherently.

This replaces the current non-polar bypass with a general final-state rule.

A mode is classified as non-polar when its oscillator strength (the squared
norm of the Born-charge-weighted eigenvector) is below a threshold, suggested
initially as `1e-4` in units normalised by the largest oscillator strength in
the mode set.  This is a test parameter, not a settled physical constant.  The
threshold should be named, documented, and validated across materials before it
becomes a default.

### Polar Finite-q Modes

For polar modes with finite `|q_ph|`, NAC depends on the direction of `q_ph`.
If the coherence regime resolves phonon momentum, the final-state key should
include a q-direction class, not the q magnitude.

```text
q_class_key = discretised(q_hat_lab, angular_tolerance_deg)
```

The angular tolerance should default to **5 degrees** and be stored as a
named resolver parameter.  It should be specified as a physical angle, not as
a rounding of vector components.  Two phonon wavevectors are placed in the same
class when the angle between them is less than `angular_tolerance_deg`.

The 5 degree value is an initial validation setting, not a final default.  It
may be too coarse near strong anisotropy, avoided crossings, or rapidly varying
NAC/EO tensors.  Tests should include a tolerance-convergence check and a
diagnostic report of the maximum tensor/frequency spread within each q class.

The current GUI parameter `q_lab_hat` (already present in `SettingsTab`) maps
directly to this tolerance and should be renamed accordingly when the refactor
is complete.

### Polar q Near Zero

For polar modes with `|q_ph|` below threshold, the resolver must not apply a
geometry switch.  It must answer the physical question:

```text
Does this near-zero q arise from the external photon momentum transfer, or from
an internal standing-wave/reflection component of a finite-film field?
```

The concrete rule, following directly from the Berreman pair classification
above, is:

- Compute `q_ext = |k_inc - k_det|` from the **external** scattering geometry
  using the macroscopic wavevectors at the medium boundaries, not from internal
  Berreman components.
- If `|q_ext|` is below the same near-zero threshold, the zero-q polar mode is
  a genuine TO-like final state and is classified `ResolvedFinalState`.
- If `|q_ext|` is above threshold but a particular Berreman pair `(i,j)` happens
  to produce a near-zero `q_ph(i,j)` (i.e. it is a cross-pair in a backscattering
  geometry), the calculator checks the depth-integration regime:
  - In `coherent_film`: classify as `CoherentExternalChannelContribution` — the
    cross-pair is a coherent interference term, not an independent TO final state.
  - In `incoherent_depth`: resolve the phonon using `q_ext`, then add the
    pair's local intensity to that external final-state channel.  The internal
    `q_ph(i,j)` remains available for diagnostics and coherent phase matching,
    but must not suppress a valid external channel.

This rule replaces all current zero-q suppression switches.  The near-zero
threshold for `|q_ph|` should be `1e-3` in units of the laser wavevector
magnitude `|k_L|`.

The `1e-3 |k_L|` value is an initial numerical threshold.  It should be
validated by near-forward and near-backscatter convergence tests, and should
not be treated as a universal physical constant.

This is the main physics issue to settle.  The implementation should be able to
test both limits:

- infinite or bulk-like forward scattering should approach TO behaviour;
- infinite or bulk-like backscattering should approach LO behaviour for a polar
  mode with q along the polar axis;
- finite coherent films should be compared against the full-field coherent
  amplitude, not automatically against an incoherent modal-pair sum.

### NAC Branch Identity

NAC branch identity must be assigned by eigenvector overlap or subspace overlap,
never by sorted eigenvalue index.  When a NAC branch is mapped back to a TO mode,
the chosen frequency, Raman tensor, linewidth, and selection flag must all come
from that mapped branch.

The recommended procedure is:

1. Compute the NAC dynamical matrix at the q direction of the pair.
2. Diagonalise to obtain NAC eigenvectors `u_NAC`.
3. For each TO mode eigenvector `u_TO`, find the NAC branch `b` that maximises
   `|<u_NAC(b) | u_TO>|^2`.
4. If the maximum overlap is below `0.5`, emit a warning and fall back to the
   closest-frequency match, logging the fallback.
5. Use the matched NAC branch `b` for all downstream quantities.

The existing `dominant_to_by_nac` map already implements this; the refactor
should make it the only branch-mapping path, with no fallback to index equality.

### Degenerate Branches and Subspaces

Degenerate or nearly degenerate optical q channels should be combined before
final-state resolution when they represent the same optical propagation
subspace.  Degenerate or nearly degenerate phonon branches should be represented
by a branch subspace rather than an arbitrary single eigenvector where needed.

#### Optical Subspace Detection

Two Berreman forward (or two backward) modes are treated as a degenerate
subspace when their eigenvalue gap satisfies

```text
|kz_a - kz_b| / mean(|kz_a|, |kz_b|) < eps_optical
```

with `eps_optical = 1e-4` as default.  When this condition holds, replace the
two individual eigenvectors with the 2x2 projector

```text
P = u_a u_a^H + u_b u_b^H
```

and compute the Raman amplitude contribution from the subspace as the sum over
any orthonormal basis of `P`.  The result is then independent of any numerical
rotation within the degenerate subspace, which resolves the fragility at normal
incidence on uniaxial crystals where the two forward modes are exactly
degenerate.

The Euclidean projector above is appropriate for the synthetic FS5 tests, where
the test vectors are explicitly orthonormal.  Real Berreman eigenvectors from
GTM may not be orthonormal under the ordinary Euclidean inner product,
especially in anisotropic or absorbing media.  The implementation should either:

- orthonormalise the degenerate electric-field basis before forming the
  projector; or
- use the correct biorthogonal projector for the Berreman eigenproblem.

The chosen convention must be documented and tested.  A synthetic orthonormal
projector test is necessary but not sufficient for the full GTM implementation.

#### Phonon Branch Subspace Detection

Two phonon branches are treated as a degenerate subspace when their NAC
frequencies satisfy

```text
|omega_a - omega_b| / mean(omega_a, omega_b) < eps_phonon
```

with `eps_phonon = 1e-3` as default.  The Raman tensor contribution from the
subspace is the sum of the individual branch tensors, which is invariant to
unitary mixing within the subspace.

Both thresholds should be named constants, exposed as resolver parameters for
testing, but not user-facing GUI options.

## Proposed Code Structure

1. `OpticalChannelResolver`
   - Computes the full incident field `E_L(z)` and reciprocal detector field
     `E_S_rec(z)` for the selected external channels using GTM.
   - Optionally decomposes those fields into internal Berreman components for
     phase-matching analysis.
   - Computes the external momentum transfer `q_ext` from the macroscopic
     scattering geometry.
   - Provides optical wavevectors, `q_ext`, and full field profiles as outputs.
   - Does not know about NAC, Raman tensors, or mode selection.

2. `PhononFinalStateResolver`
   - Converts `(mode_idx, q_ph, layer, external_channel, coherence_regime)` to
     a final-state classification.
   - Owns q-independent mode handling, NAC branch mapping, q-direction
     tolerances, the polar/non-polar threshold, the near-zero-q rule, and the
     degenerate subspace logic.
   - Receives `q_ext` from the `OpticalChannelResolver` to distinguish genuine
     forward-scatter zero-q from cross-pair zero-q.
   - For coherent finite-film polar modes, may return pair-resolved coherent
     contributions rather than forcing every pair into a resolved final-state
     bucket.
   - Can be unit-tested without GTM field integrals.

3. `RamanAmplitudeAccumulator`
   - Receives classified amplitude records
     `(final_state_key, amplitude_contribution)`.
   - Sums amplitudes by final-state key.
   - Keeps coherent external-channel contributions coherent.
   - Applies incoherent sums only between physically resolved final states or
     explicitly incoherent layers/depth intervals.

This separation should remove most modal-pair special cases from the main
calculator loop.

## Implementation Phases

The refactor should proceed in phases to keep the existing `geometry` and
`dominant_mode` paths working throughout.  Those paths should not be touched
until Phase 4.

Backward compatibility with the current `modal_pairs` behaviour is not a
requirement.  The layered Raman implementation has not been released, so the
priority is a physically coherent design rather than preserving existing
modal-pair switches, diagnostics, or output conventions.  Temporary parallel
paths are useful for validation during the refactor, but they should not become
permanent compatibility layers.

### Phase 1: PhononFinalStateResolver in isolation

Implement `PhononFinalStateResolver` as a standalone class with no dependency on
GTM or the GUI.  Write unit tests covering:

- Non-polar mode classification (FS6, FS9).
- NAC branch mapping by overlap, including swapped-eigenvalue case (FS4).
- Degenerate optical subspace projector (FS5, synthetic).
- Near-zero-q rule under each coherence regime (FS14, synthetic).
- q_class_key binning with angular tolerance (FS7).
- Pair-dependent polar tensor classification and finite-thickness
  phase-matching weights (FS16, synthetic).

These tests should pass before any calculator code is changed.

### Phase 2: RamanAmplitudeAccumulator in isolation

Implement `RamanAmplitudeAccumulator` and test coherent vs incoherent summation
separately from field integrals.  Verify that:

- Amplitudes assigned to the same key are summed coherently.
- Amplitudes assigned to different keys contribute incoherently.
- Layer coherence group controls cross-layer accumulation (FS11).

### Phase 3: OpticalChannelResolver and wiring

Implement `OpticalChannelResolver` wrapping the existing GTM field calculation.
Add a parallel code path in `LayeredRamanCalculator` that calls the three new
components alongside the existing `modal_pairs` path and compares outputs for
known cases (FS1, FS2, FS3).  Both paths run; the new path is used only for
diagnostics at this stage.

### Phase 4: Replace modal_pairs

Once the parallel path produces correct results for all test matrix cases,
replace the `modal_pairs` code path with the new implementation.  Remove
mode-specific patches (non-polar bypass, zero-q suppression switch, ad hoc
branch-index assumptions) from the main loop.

### Phase 5: Cleanup

Remove the old `modal_pairs` implementation.  Rename GUI parameters as needed.
Update the `SettingsTab` to expose coherence regime selection.  Run the full
`pdmake tests` suite.

## Small Test Matrix

The tests below are intended to define the refactor target before implementation.
Tests marked **[synthetic]** require no DFT data and should be implemented as
pure-Python unit tests in Phase 1 or 2.  Tests marked **[ZnO]** require the
CASTEP ZnO dataset from `Examples/Crystal_Raman/`.  Tests marked
**[validation]** are slow convergence checks that are not run on every commit.

| ID | Case | Geometry | Expected Behaviour | Type |
|----|------|----------|--------------------|------|
| FS1 | c-plane ZnO `z(xx)-z` | normal-incidence backscatter | Spectrum contains A1(LO) and E2(high); A1(TO) is not the dominant polar line. | [ZnO] |
| FS2 | c-plane ZnO forward scattering | normal incidence, collection on substrate side | Polar A1 response is TO-dominated; LO branch is absent or strongly suppressed. | [ZnO] |
| FS3 | same external channel, different internal decomposition | coherent finite film | Changing the Berreman modal basis does not change the full external-channel Raman amplitude. | [synthetic] |
| FS4 | polar mode with branch reordering | synthetic swapped NAC eigenvalues | The TO mode maps to the correct NAC branch by overlap; no intensity appears on the wrong sorted branch. | [synthetic] |
| FS5 | degenerate optical q subspace | normal incidence, split Berreman basis | Rotating or numerically splitting the optical basis does not change the intensity. | [synthetic] |
| FS6 | non-polar E2-like mode | backscatter and oblique incidence | Intensity is invariant to q resolution because q is not a physical final-state label for this mode. | [ZnO] |
| FS7 | same q direction, different q magnitude | anisotropic oblique bulk-limit stack | Pairs with the same q direction and same phonon branch sum coherently when q is resolved only by NAC direction. | [synthetic] |
| FS8 | opposite q directions | bulk-limit backscatter with resolved phonon momentum | Opposite polar q directions are distinct final states unless symmetry explicitly makes them indistinguishable. | [synthetic] |
| FS9 | mode-selection mask with NAC reordering | reordered branches and disabled TO mode | Selection follows the mapped branch, not the raw sorted NAC index. | [synthetic] |
| FS10 | incoherent depth integration | thick layer/bulk-like limit | Final-state grouping remains the same; only depth accumulation changes from amplitude to local intensity. | [synthetic] |
| FS11 | coherent vs incoherent layer combination | subdivided identical layer | Final-state identity is unchanged; layer coherence policy alone controls cross-layer amplitude terms. | [synthetic] |
| FS12 | Porto compound polarisation | `z(x+y,x+y)-z` | Jones-vector excitation/detection uses the same final-state classification as pure `x` or `y` channels. | [ZnO] |
| FS13 | near-forward limit | external scattering angle tends to zero | Polar response tends continuously to TO behaviour when the external momentum transfer tends to zero. | [ZnO] |
| FS14 | near-backscatter internal zero-q component | coherent film with reflections | Internal zero-q modal products do not create an independent TO final state unless the coherence model says phonon q is resolved. | [synthetic] |
| FS15 | thickness convergence | finite coherent film to bulk-like slab | The modal/final-state result approaches the expected bulk LO/TO limits without changing ad hoc switches. | [validation] |
| FS16 | polar tensor bulk limit | synthetic pair phases with LO-like and TO-like tensors | The coherent pair-resolved formula tends to TO in forward scattering and LO in backscattering as thickness increases. | [synthetic]/[validation] |

### Synthetic Test Setup for FS3 and FS5

These two tests are the most important correctness invariants and require no
real phonon data.  The recommended setup is:

**FS3 — Basis invariance:**
1. Construct a single-layer stack with a known 3x3 Raman tensor `R` (e.g.
   diagonal, representing a non-degenerate polar mode).
2. Compute the GTM forward and backward Berreman eigenvectors for a chosen angle
   of incidence, giving a 4x4 eigenvector matrix `M`.
3. Apply a random unitary rotation `U` to the two forward-mode columns of `M`,
   producing a rotated basis `M'` that spans the same forward subspace.
4. Compute the external-channel Raman amplitude using both `M` and `M'`.
5. Assert that the full external-channel amplitudes agree to within `1e-10`.
   (They will only agree if the accumulator sums all pairs coherently before
   squaring.)
6. Repeat the same check for a unitary rotation of the reciprocal detector
   subspace.  The external-channel amplitude must be invariant to basis changes
   in both `E_L` and `E_S_rec`; testing only the incident side is not sufficient
   for the final implementation.

**FS5 — Degenerate optical subspace at normal incidence:**
1. Set up a uniaxial medium at normal incidence.  The two forward Berreman modes
   are exactly degenerate.
2. Parameterise the degenerate basis by an angle `theta` (rotation within the
   2D subspace).
3. Compute the Raman intensity for `theta = 0, 15, 30, 45` degrees.
4. Assert that the intensity is independent of `theta` to within `1e-8`.
5. Add a second implementation-level test where the reciprocal detector
   degenerate subspace is rotated.  The same basis-invariance requirement
   applies on the detector side.

**FS16 — Polar tensor bulk limit:**
1. Construct synthetic pair amplitudes with two LO-like co-propagating pairs and
   two TO-like cross-pairs.
2. Assign pair tensors `R_LO = R(q_ext)` and `R_TO = R(q=0)`.
3. Multiply each pair amplitude by the analytic finite-thickness factor
   `F_ij(L) = L exp(i Delta k_ij L/2) sinc(Delta k_ij L/2)`.
4. For forward scattering, set `q_ext = 0` and verify that the normalized
   large-`L` amplitude tends to the TO tensor contribution.
5. For backscattering, set `q_ext` along the polar axis and verify that the
   normalized large-`L` amplitude tends to the LO tensor contribution.
6. Verify that no hard zero-q suppression or geometry switch is needed; the
   limit follows from the phase-matching factor.

## Acceptance Criteria

- The final-state resolver can be tested without constructing a full GUI
  scenario.
- The main modal-pair loop contains no mode-specific patches such as "E2 bypass"
  or geometry-specific ad hoc frequency filtering.
- Backscatter and forward-scatter differences follow from external momentum
  transfer and the coherence regime, not from a hard-coded geometry policy.
- The implementation distinguishes internal optical field components from
  externally prepared/detected scattering channels.
- Direct calculator diagnostics can still expose individual final-state
  contributions for debugging.
- The default GUI behaviour corresponds to physical Porto/backscatter selection
  rules, while advanced scripts can request diagnostic modal decompositions
  explicitly.
- All synthetic tests (FS3–FS11, FS14) pass as pure-Python unit tests with no
  GUI or DFT data dependency.
- The `geometry` and `dominant_mode` code paths produce identical results before
  and after the refactor, verified by the existing `pdmake tests` suite.
- FS15 (thickness convergence) is tracked as a validation test, not a per-commit
  regression fixture, due to its computational cost.
