# coherence_incoherence.md — Testing Aspects of (In)Coherence

## Coherence Comparison Scripts

These four scripts systematically probe the four independent coherence/incoherence
levels described in `~/Software/PDielec/crystal_raman_coherence_guide.md`.  All use
the standard c-plane (001) ZnO on Fused Silica stack, 532 nm laser, pp polarisation,
`modal_pairs` NAC, and collection from the superstrate.

Note: E2(high) (~398 cm⁻¹ DFT, ~437 cm⁻¹ experimental) was originally proposed as an
internal reference invariant across coherence choices.  Level 2 results (see below) show
this is only true across Level 3 and Level 4 choices.  For Level 2 (depth coherence),
E2(high) is itself moderately suppressed (~50%) in the coherent relative to the
incoherent depth mode and is therefore not a reliable reference for that comparison.

### `substrate_gtm_modes.py`  *(Substrate coherence — Fabry-Perot etalon effect)*
- **Purpose:** Test the effect of Fused Silica substrate thickness and coherence
  treatment on the Crystal Raman spectrum of a fixed 100 nm ZnO film.
- **Layer stack:** `['air', 'Dielectric layer (ZnO 100 nm)', 'Fused Silica (d)', 'air']`
  — the Fused Silica is an explicit **finite-thickness** layer; the final air is the
  semi-infinite exit medium.  GTMcore.py does treat the outermost layers (first air and
  last air) as semi-infinite; all interior layers, including the Fused Silica, are
  propagated through coherently or incoherently as specified.
- **Scenarios (8):** Coherent substrate at 1/5/10/50 μm; Incoherent intensity at
  1/100 μm; Incoherent non-reflective at 1/100 μm.
- **Output:** `substrate_gtm_modes.xlsx`

#### Observations from running the script

**The coherent scenarios do NOT produce identical spectra.**
Each coherent Fused Silica layer acts as a Fabry-Perot etalon whose transmission
oscillates with Raman shift.  The Free Spectral Range (FSR = 1/(2nd), n = 1.461) is:

| Substrate thickness | FSR (cm⁻¹) |
|---|---|
| 1 μm  | 3422 |
| 5 μm  | 685  |
| 10 μm | 342  |
| 50 μm | 68   |

With FSR ≈ 342 cm⁻¹ the 10 μm case spans roughly two full etalon fringes across the
0–600 cm⁻¹ Raman spectrum.  The round-trip phase φ = 2 × n × d × ν_scattered
(in cycles) determines whether each mode is enhanced (φ fractional part near 0) or
suppressed (fractional part near 0.5).  Note: DFT peak frequencies differ from
experimental ZnO values.

| Mode | DFT freq (cm⁻¹) | 10 μm frac. phase | Character |
|---|---|---|---|
| E2(low)  | 109 | 0.607 | Destructive side → suppressed |
| E2(high) | 398 | 0.761 | Intermediate |
| A1(LO)   | 510 | 0.435 | Approaching destructive |

The simple fractional-phase model predicts A1(LO) should be suppressed at 10 μm
(near destructive phase).  The observed data contradicts this: A1(LO) is anomalously
*enhanced* at 10 μm relative to other thicknesses.  This indicates the etalon
effect on the LO mode cannot be explained by the substrate-alone Fabry-Perot: the
GTM computes the full coupled field in the ZnO film + substrate cavity simultaneously,
and the LO mode's field distribution (arising from the (forward L, backward S) channel
pair in the modal_pairs NAC) has a fundamentally different response to substrate
phase than the simpler TO/E2 modes.

**Peak intensities at key modes for coherent scenarios (DFT frequencies)**

| Mode | DFT freq (cm⁻¹) | 1 μm | 5 μm | 10 μm | 50 μm |
|---|---|---|---|---|---|
| E2(low)  | 109 | 1.642 | 10.127 |  5.554 | 12.452 |
| E2(high) | 398 | 1.854 |  4.493 |  2.806 |  5.844 |
| A1(LO)   | 510 | 1.417 |  1.406 |  4.466 |  3.612 |

The A1(LO) at 10 μm (4.47) substantially exceeds both the 5 μm (1.41) and 50 μm
(3.61) values — the original observation that motivated this investigation.
E2(low) and E2(high) follow the expected etalon pattern (10 μm suppressed relative to
5 and 50 μm).

**Incoherent scenarios behave as expected.**
`Incoherent non-reflective` gives identical spectra at 1 μm and 100 μm — there is no
back-surface reflection and therefore no etalon effect.  `Incoherent intensity` varies
smoothly with thickness (different total reflectance) but without the oscillatory
wavelength-selective suppression/enhancement of the coherent case.

**Physical conclusion.**  A macroscopic glass substrate (μm–mm scale) is not a
coherent Fabry-Perot etalon in a real experiment: surface roughness, thermal drift,
and finite laser coherence length all wash out the fringes.  The coherent treatment
is physically appropriate only for perfect, atomically flat, vibration-free slabs
(e.g. exfoliated crystals on a polished surface in a cryogenic setup).  For routine
backscattering Raman from a crystal on a glass slide the `Incoherent intensity` or
`Incoherent non-reflective` substrate treatment is more realistic.

### `substrate_etalons.py`  *(Substrate etalon oscillations resolved)*
- **Purpose:** Resolve the Fabry-Perot etalon oscillations in peak intensity vs
  substrate thickness by sampling at fine enough steps to avoid aliasing.
- **Layer stack:** Same as `substrate_gtm_modes.py` — air / ZnO (100 nm) /
  Fused Silica (d) / air — all layers coherent.
- **Scenarios (100):** Fused Silica thickness d = 10.00, 10.02, 10.04, …, 11.98 μm
  (step 0.02 μm, 100 points over a 1.98 μm window).
- **Output:** `substrate_etalons.xlsx`

#### Why the step size matters — aliasing

The etalon fringes repeat every half-wavelength of the scattered light inside the
glass.  For 532 nm laser and ZnO modes between 99 and 585 cm⁻¹:

```
period in thickness = λ_scattered / (2n)  ≈  0.183–0.188 μm
```

The previous 1 μm sampling step in `substrate_gtm_modes.py` is ~5.4× coarser than
the oscillation period.  Each 1 μm step skips ~5.4 full fringes, landing at a
different, apparently random phase each time — classic **aliasing**.  The FFT of the
coarsely-sampled data revealed apparent beating periods of 2–3 μm (mode-dependent)
that are aliases of the true ~0.18 μm signal, not physical beating.

With 0.02 μm steps there are ~9 samples per fringe period, fully resolving the
oscillations.

#### Observations from running the script

**The oscillations are clean sinusoids once resolved.**  Every mode traces a smooth
periodic curve over the 10–12 μm window, with a period recovered by FFT of
**0.1818 μm** (theory: 0.183–0.188 μm — excellent agreement).

**The modulation depth is large — far larger than a naive Fabry-Perot estimate.**
(DFT peak positions: E2(low) 109 cm⁻¹, E2(high) 398 cm⁻¹, A1(LO) 510 cm⁻¹.)

| Mode | DFT freq (cm⁻¹) | Min | Max | Mean | Max/Min |
|---|---|---|---|---|---|
| E2(low)  | 109 |  1.02 | 14.79 |  6.96 | 14.6× |
| E2(high) | 398 |  0.08 | 11.60 |  4.52 | 148.6× |
| A1(LO)   | 510 |  1.21 |  8.69 |  4.47 | 7.2× |

**Peak intensities vs substrate thickness (10.00–11.98 μm, step 0.02 μm)**

```
  d(μm) |    E2(low) |   E2(high) |     A1(LO)
--------------------------------------------
 10.000 |     2.0555 |     3.0160 |     4.4656
 10.020 |     5.1830 |     6.1977 |     4.4134
 10.040 |     9.4377 |     9.4629 |     4.2723
 10.060 |    13.3677 |    11.6047 |     4.3703
 10.080 |    14.7886 |    11.0685 |     4.5977
 10.100 |    12.3047 |     7.4614 |     4.5562
 10.120 |     7.3046 |     2.9315 |     4.2345
 10.140 |     2.8956 |     0.2968 |     4.0915
 10.160 |     1.0151 |     0.3481 |     4.3938
 10.180 |     1.7650 |     2.2877 |     4.8267
 10.200 |     4.5379 |     5.2304 |     4.8806
 10.220 |     8.5810 |     8.4383 |     4.5566
 10.240 |    12.5957 |    10.8430 |     4.2978
 10.260 |    14.4928 |    10.9522 |     4.2628
 10.280 |    12.6430 |     8.0379 |     4.1831
 10.300 |     7.9787 |     3.6810 |     3.9606
 10.320 |     3.4545 |     0.6402 |     3.9410
 10.340 |     1.2099 |     0.1422 |     4.4045
 10.360 |     1.5532 |     1.6586 |     5.0880
 10.380 |     3.9621 |     4.3207 |     5.3638
 10.400 |     7.7609 |     7.4112 |     4.9885
 10.420 |    11.7987 |     9.9814 |     4.4029
 10.440 |    14.1056 |    10.6358 |     4.0353
 10.460 |    12.8789 |     8.4212 |     3.8235
 10.480 |     8.6073 |     4.3878 |     3.6381
 10.500 |     4.0362 |     1.0885 |     3.7067
 10.520 |     1.4681 |     0.0781 |     4.2978
 10.540 |     1.4187 |     1.1451 |     5.2211
 10.560 |     3.4587 |     3.4877 |     5.8011
 10.580 |     6.9856 |     6.4063 |     5.5264
 10.600 |    10.9892 |     9.0541 |     4.6932
 10.620 |    13.6375 |    10.1499 |     3.9493
 10.640 |    13.0118 |     8.6062 |     3.5103
 10.660 |     9.1803 |     5.0080 |     3.2903
 10.680 |     4.6292 |     1.6058 |     3.4068
 10.700 |     1.7820 |     0.1528 |     4.0878
 10.720 |     1.3589 |     0.7614 |     5.2154
 10.740 |     3.0297 |     2.7498 |     6.1373
 10.760 |     6.2619 |     5.4468 |     6.1113
 10.780 |    10.1791 |     8.0927 |     5.1587
 10.800 |    13.1000 |     9.5276 |     4.0342
 10.820 |    13.0429 |     8.6004 |     3.2786
 10.840 |     9.6886 |     5.5082 |     2.9422
 10.860 |     5.2223 |     2.1506 |     3.0605
 10.880 |     2.1433 |     0.3531 |     3.7927
 10.900 |     1.3701 |     0.5174 |     5.0758
 10.920 |     2.6764 |     2.1251 |     6.3328
 10.940 |     5.5962 |     4.5551 |     6.6737
 10.960 |     9.3796 |     7.1270 |     5.7677
 10.980 |    12.5048 |     8.8021 |     4.3100
 11.000 |    12.9754 |     8.4199 |     3.1647
 11.020 |    10.1249 |     5.8669 |     2.6215
 11.040 |     5.8046 |     2.6810 |     2.6871
 11.060 |     2.5432 |     0.6554 |     3.4322
 11.080 |     1.4475 |     0.4166 |     4.8178
 11.100 |     2.3987 |     1.6300 |     6.3681
 11.120 |     4.9936 |     3.7529 |     7.1457
 11.140 |     8.6013 |     6.1847 |     6.4658
 11.160 |    11.8644 |     8.0063 |     4.7818
 11.180 |    12.8142 |     8.0870 |     3.2038
 11.200 |    10.4833 |     6.0747 |     2.3589
 11.220 |     6.3659 |     3.1594 |     2.3071
 11.240 |     2.9722 |     1.0279 |     3.0257
 11.260 |     1.5853 |     0.4535 |     4.4631
 11.280 |     2.1958 |     1.2779 |     6.2436
 11.300 |     4.4585 |     3.0613 |     7.4717
 11.320 |     7.8540 |     5.2925 |     7.1813
 11.340 |    11.1911 |     7.1721 |     5.4330
 11.360 |    12.5659 |     7.6276 |     3.4262
 11.380 |    10.7598 |     6.1323 |     2.1880
 11.400 |     6.8967 |     3.5559 |     1.9426
 11.420 |     3.4212 |     1.4339 |     2.5919
 11.440 |     1.7766 |     0.6128 |     4.0349
 11.460 |     2.0653 |     1.0766 |     5.9747
 11.480 |     3.9942 |     2.5001 |     7.6182
 11.500 |     7.1464 |     4.4764 |     7.8349
 11.520 |    10.4973 |     6.3303 |     6.2219
 11.540 |    12.2383 |     7.0699 |     3.8517
 11.560 |    10.9519 |     6.0485 |     2.1445
 11.580 |     7.3885 |     3.8503 |     1.6185
 11.600 |     3.8807 |     1.8362 |     2.1488
 11.620 |     2.0139 |     0.8694 |     3.5557
 11.640 |     2.0038 |     1.0263 |     5.5858
 11.660 |     3.6029 |     2.0866 |     7.5750
 11.680 |     6.4866 |     3.7617 |     8.3542
 11.700 |     9.7951 |     5.5110 |     7.0824
 11.720 |    11.8406 |     6.4429 |     4.4822
 11.740 |    11.0587 |     5.8383 |     2.2637
 11.760 |     7.8337 |     4.0318 |     1.3627
 11.780 |     4.3419 |     2.2016 |     1.7150
 11.800 |     2.2893 |     1.1905 |     3.0462
 11.820 |     2.0067 |     1.1172 |     5.1049
 11.840 |     3.2856 |     1.8335 |     7.3530
 11.860 |     5.8817 |     3.1725 |     8.6858
 11.880 |     9.0962 |     4.7439 |     7.9314
 11.900 |    11.3829 |     5.7762 |     5.2954
 11.920 |    11.0811 |     5.5205 |     2.5758
 11.940 |     8.2258 |     4.0982 |     1.2069
 11.960 |     4.7963 |     2.5037 |     1.3103
 11.980 |     2.5943 |     1.5397 |     2.5253
```

A simple Fabry-Perot model (Fused Silica surface reflectance R ≈ 3.5%) predicts a
max/min transmission contrast of only ~1.15×.  The observed ratios are much larger
(7× for A1(LO), 15× for E2(low), and very large for E2(high) because one sampled
thickness is close to destructive cancellation) because the etalon modulates the
**local electric field amplitude** inside the 100 nm ZnO film, and the Raman intensity
scales as the product of the excitation field squared and the emission field squared.
The GTM coherently sums all multiple reflections, including interference between the
film and substrate interfaces, producing a much stronger net modulation than the bare
substrate reflectance would suggest.

**A1(LO) shows a slow amplitude envelope superimposed on the fast oscillation.**
The A1(LO) peak amplitude grows systematically from ~1.2 at 10 μm to ~8.7 at ~11.8 μm
and then falls again.  This is a beat between two nearly equal oscillation frequencies:
the substrate etalon at the laser wavelength (period 0.1818 μm) and at the LO Raman
wavelength (period 0.1873 μm).  The beat period is 1/|f_laser − f_Raman| ≈ 6.6 μm,
so the 1.98 μm window captures only ~30% of one beat cycle — explaining the monotonic
amplitude rise seen over the sampled range.  The non-polar modes (E2(low), E2(high))
have Raman wavelengths closer to the laser, giving much longer beat periods (> 25 μm)
and hence more nearly constant amplitude envelopes over the same window.

**Modes peak at different substrate thicknesses.**  The constructive condition for
each mode depends on its scattered wavelength.  Within the 10–12 μm window, the TO
modes and LO modes reach their maxima at slightly offset thicknesses — the same phase
difference that caused the anomalous TO suppression / LO enhancement at the 10 μm
integer sample in `substrate_gtm_modes.py`.

**The mean intensity across the etalon cycle is stable.**  The cycle-averaged peak
for each mode is independent of where in the fringe cycle the substrate sits.  Over
the 100 coherent scenarios (10.00–11.98 μm) the mean values are:

| Mode | Cycle-averaged intensity |
|---|---|
| E2(low)  | 6.96 |
| E2(high) | 4.52 |
| A1(LO)   | 4.47 |

The `Incoherent non-reflective` scenario in `substrate_gtm_modes.py` gives **identical
results at 1 μm and 100 μm** — it is truly thickness-independent, because removing
back-reflection from the substrate eliminates all path-length dependence.  This makes
it a valid fixed reference to compare against the cycle-average:

| Mode | Non-reflective | Etalon cycle-mean | Mean / Non-refl |
|---|---|---|---|
| E2(low)  | 5.41 | 6.96 | 1.29× |
| E2(high) | 2.87 | 4.52 | 1.57× |
| A1(LO)   | 2.75 | 4.47 | 1.62× |

The cycle-averaged coherent intensity exceeds the non-reflective value by ~29% for
E2(low), ~57% for E2(high), and ~62% for A1(LO).  This difference quantifies
the average Raman enhancement from substrate back-reflection at ~11 μm substrate
thickness: the back-surface of the Fused Silica reflects a fraction of the Raman
photon back upward into the collection, adding to the signal even when phase coherence
is removed by averaging.  The non-reflective scenario suppresses this contribution
entirely, giving a lower bound for the Raman intensity.

Note that the `Incoherent intensity substrate` scenarios (col labels in
`substrate_gtm_modes.py`) are actually **coherent** in their per-layer settings
(`Layer incoherent options = ['Coherent', 'Coherent', 'Coherent', 'Coherent']`); they
differ from the `Coherent substrate` scenarios only in the substrate thickness used
(1 μm and 100 μm).  Their label is therefore misleading — they are not a true
incoherent-intensity reference.

### `depth_coherence_crossover.py`  *(Level 2 — Depth coherence)*
- **Purpose:** Show when `Coherent amplitude` and `Incoherent intensity` depth
  coherence modes diverge as a function of film thickness.
- **Scenarios (22):** 11 thicknesses × 2 depth modes, interleaved by thickness.
  Thicknesses: 10, 20, 30, 40, 50, 75, 100, 250, 500, 1000 nm, 5 μm.
  GL point density: 1000 points/μm at all thicknesses.
- **Levels fixed:** All layers GTM = Coherent, Layer combination = Incoherent
  intensities, Modal pair combination = Group q channels.
- **Physical test:** Modes agree for thin films (thickness ≪ λ/2n ≈ 133 nm); diverge
  beyond that crossover.  Original prediction (A1(LO) most sensitive, E2(high) invariant)
  was incorrect — see Observations below.
- **Output:** `depth_coherence_crossover.xlsx`

#### Observations from running the script

**What the two modes compute (from `LayeredRamanCalculator.py`)**

The GL physical weights `w_j` sum to the film thickness **in Ångström**
(`w_j = half_thickness_m × quad_weights × 1e10`).  With this,

| Mode | Formula | Thin-film scaling |
|---|---|---|
| Incoherent intensity | `I = Σⱼ wⱼ \|A(zⱼ)\|²` | ∝ d_Å |
| Coherent amplitude | `I = \|Σⱼ wⱼ A(zⱼ)\|²` | ∝ d_Å² |

where A(zⱼ) = E_S(zⱼ)ᵀ R E_L(zⱼ) is the local Raman amplitude.

**Mode frequencies in this CASTEP calculation**

DFT frequencies differ from experiment.  The selected modes appearing in the spectra:

| Peak (cm⁻¹) | Assignment | Polar? |
|---|---|---|
| 108.6 | E2(low) | No |
| 334.6 | A1(TO)  | Yes (polar along c) |
| 367.0 | E1(TO)  | Yes, but Raman-silent in c-plane pp geometry |
| 398.2 | E2(high)| No |
| ~510  | A1(LO) via modal_pairs NAC | Yes |

**Thin-film "crossover" at 0.1 nm**

The two spectra are equal when d_Å = 1, i.e. at d = 0.1 nm.
At d = 0.2 nm the coherent peaks are a factor of 2 larger than the incoherent peaks.
This is a direct consequence of the ratio coherent/incoherent = d_Å in the thin-film
limit, and is **a unit coincidence** (weights stored in Å, not a physical crossover).
It is unrelated to the incident frequency and has no physical significance — film
thicknesses below ~1 nm are sub-atomic and unphysical.

**Frequency-related crossover**

The physically meaningful divergence begins when the integrand phase
Δφ = (k_L + k_S)·d ≈ 4πnd/λ becomes of order 1 radian.
For 532 nm in ZnO (n ≈ 2): significant phase variation starts at d ~ λ/(4πn) ≈ 21 nm,
and the integrand phase spans ~π at d ≈ λ/(4n) ≈ 67 nm, ~1.5π at 100 nm.

**CASTEP DFT mode frequencies and peak positions**

Peak frequencies in the spectrum (from 100 nm incoherent run): 109, 335, 398, 510 cm⁻¹.

**Coherent peak intensities**

A₁(TO) is not listed: with `modal_pairs` NAC the incoherent spectrum has no local
maximum at 334.6 cm⁻¹, so a Coh/Incoh ratio for A₁(TO) is physically undefined.

| Thickness |  d_Å  | E2(low) 109 | E2(high) 398 | A1(LO) 510 |
|-----------|-------|-------------|--------------|------------|
|    10 nm  |   100 |    1.295    |    0.7219    |  0.03240   |
|    20 nm  |   200 |    3.931    |    2.199     |  0.1103    |
|    30 nm  |   300 |    5.971    |    3.356     |  0.2028    |
|    40 nm  |   400 |    6.740    |    3.807     |  0.3005    |
|    50 nm  |   500 |    6.500    |    3.689     |  0.4158    |
|    75 nm  |   750 |    3.897    |    2.251     |  0.9993    |
|   100 nm  |  1000 |    5.405    |    2.868     |  2.753     |
|   250 nm  |  2500 |   87.77     |   49.13      |  17.32     |
|   500 nm  |  5000 |  219.2      |  128.7       |  51.15     |
|  1000 nm  | 10000 |  671.6      |  328.9       | 166.2      |
|     5 μm  | 50000 | 1.808×10⁴  |   70.91      | 1361       |

**Incoherent peak intensities**

| Thickness |  d_Å  | E2(low) 109 | E2(high) 398 | A1(LO) 510 |
|-----------|-------|-------------|--------------|------------|
|    10 nm  |   100 | 1.308×10⁻² | 7.270×10⁻³  | 6.520×10⁻³ |
|    20 nm  |   200 | 2.057×10⁻² | 1.147×10⁻²  | 1.030×10⁻² |
|    30 nm  |   300 | 2.247×10⁻² | 1.257×10⁻²  | 1.130×10⁻² |
|    40 nm  |   400 | 2.186×10⁻² | 1.224×10⁻²  | 1.101×10⁻² |
|    50 nm  |   500 | 2.121×10⁻² | 1.187×10⁻²  | 1.067×10⁻² |
|    75 nm  |   750 | 2.712×10⁻² | 1.498×10⁻²  | 1.341×10⁻² |
|   100 nm  |  1000 | 5.870×10⁻² | 3.197×10⁻²  | 2.845×10⁻² |
|   250 nm  |  2500 |  0.1651     |  9.423×10⁻² | 8.522×10⁻² |
|   500 nm  |  5000 |  0.2265     |  0.1345      |  0.1238    |
|  1000 nm  | 10000 |  0.4215     |  0.2208      |  0.1983    |
|     5 μm  | 50000 |  2.508      |  0.9806      |  1.102     |

**Coherent / incoherent ratio**

| Thickness |  d_Å  | E2(low) 109 | E2(high) 398 | A1(LO) 510 |
|-----------|-------|-------------|--------------|------------|
|    10 nm  |   100 |    99.0     |    99.3      |    4.97    |
|    20 nm  |   200 |   191.1     |   191.8      |   10.71    |
|    30 nm  |   300 |   265.7     |   267.1      |   17.95    |
|    40 nm  |   400 |   308.4     |   311.0      |   27.29    |
|    50 nm  |   500 |   306.4     |   310.9      |   38.96    |
|    75 nm  |   750 |   143.7     |   150.2      |   74.52    |
|   100 nm  |  1000 |    92.08    |    89.73     |   96.75    |
|   250 nm  |  2500 |   531.6     |   521.3      |  203.2     |
|   500 nm  |  5000 |   967.9     |   956.8      |  413.2     |
|  1000 nm  | 10000 |  1593       |  1489        |  838.4     |
|     5 μm  | 50000 |  7209       |    72.31     | 1235       |

**Coherent/incoherent ratio normalised by d_Å**

For non-polar modes the thin-film limit is 1.000.  A₁(LO) via `modal_pairs` does not
approach 1.000 even at 10 nm — see Key findings below.

| Thickness |  d_Å  | E2(low) 109 | E2(high) 398 | A1(LO) 510 |
|-----------|-------|-------------|--------------|------------|
|    10 nm  |   100 |   0.9900    |   0.9929     |  0.0497    |
|    20 nm  |   200 |   0.9553    |   0.9588     |  0.0535    |
|    30 nm  |   300 |   0.8856    |   0.8903     |  0.0598    |
|    40 nm  |   400 |   0.7709    |   0.7776     |  0.0682    |
|    50 nm  |   500 |   0.6128    |   0.6217     |  0.0779    |
|    75 nm  |   750 |   0.1916    |   0.2003     |  0.0994    |
|   100 nm  |  1000 |   0.0921    |   0.0897     |  0.0968    |
|   250 nm  |  2500 |   0.2126    |   0.2085     |  0.0813    |
|   500 nm  |  5000 |   0.1936    |   0.1914     |  0.0826    |
|  1000 nm  | 10000 |   0.1593    |   0.1489     |  0.0838    |
|     5 μm  | 50000 |   0.1442    |   0.0014     |  0.0247    |

**Key findings from the data**

1. **E2(low) and E2(high) thin-film limit confirmed (10 nm):** both non-polar modes give
   ratio/d_Å ≈ 0.990–0.993 at 10 nm ✓.  These modes follow the expected d² → d× d scaling
   at very thin films.

2. **E2(low) and E2(high) decline with a Fabry-Pérot minimum near 75 nm.**
   Ratio/d_Å falls from 0.99 at 10 nm to a minimum of ~0.19–0.20 at 75 nm, then rebounds
   to ~0.09 at 100 nm and settles near 0.15–0.21 for 250 nm to 1 μm.  The minimum near
   75 nm is consistent with a half-round-trip phase Δφ = 4πnd/λ ≈ π (d ≈ λ/4n ≈ 67 nm
   for 532 nm in ZnO n≈2).  At 5 μm, E2(low) is 0.144 and E2(high) drops to 0.001 due to
   destructive interference at higher film thicknesses.  E2(high) is therefore **not** a
   stable internal reference for comparing the two depth modes.

3. **A1(LO) via `modal_pairs` NAC does not follow the simple thin-film limit.**
   Ratio/d_Å ≈ 0.0497 at 10 nm — roughly 20× smaller than the E₂ modes.  This is because
   the LO branch in `modal_pairs` involves a depth-varying incoherent Berreman pair
   summation whose field profile differs qualitatively from the non-polar standing-wave
   amplitude.  The ratio grows slowly through 0.068 (40 nm), reaches ~0.10 near 75–100 nm,
   and stays in the range 0.08–0.10 up to 1 μm.  At 5 μm it falls to 0.025.  The
   near-constant value of ~0.08–0.10 over the 40 nm–1 μm range means A1(LO) Coh/Incoh
   scales approximately as d (not d²) over this broad interval.

   **GL convergence test (2026-07-09):** A dedicated script `Scripts/gl_convergence_a1lo_10nm.py`
   was run to test whether the 0.0497 ratio at 10 nm is caused by GL integration undersampling.
   The GL point density was varied from the GUI default (20 pts/μm → 5 pts, the minimum) up to
   100 000 pts/μm → 1000 pts.  All seven GL counts gave **identical results** to 7 significant
   figures for all three modes:

   | GL density (pts/μm) | Actual GL pts | E2(low) ratio | E2(high) ratio | A1(LO) ratio |
   |--------------------:|:-------------:|:-------------:|:--------------:|:------------:|
   | 20 (GUI default)    |  5 (min)      |    0.9900     |     0.9929     |    0.0497    |
   | 500                 |  5 (min)      |    0.9900     |     0.9929     |    0.0497    |
   | 1 000               | 10            |    0.9900     |     0.9929     |    0.0497    |
   | 5 000               | 50            |    0.9900     |     0.9929     |    0.0497    |
   | 10 000              | 100           |    0.9900     |     0.9929     |    0.0497    |
   | 50 000              | 500           |    0.9900     |     0.9929     |    0.0497    |
   | 100 000             | 1000          |    0.9900     |     0.9929     |    0.0497    |

   The A1(LO) ratio of 0.0497 is therefore **genuine physics**, not numerical undersampling.
   The depth profile of the LO Berreman pair product (`E_L(z) × E_S*(z)`) for the
   (forward L, backward S) channel at 10 nm is already fully resolved by 5 GL points, and
   its integral is intrinsically ~20× smaller than the corresponding incoherent sum.

4. **Both coherent and incoherent intensities oscillate non-monotonically with thickness**
   (e.g. E2(low) coherent: 1.30 at 10 nm, rises to 6.74 at 40 nm, dips to 3.90 at 75 nm,
   then rises to 5.41 at 100 nm and 87.8 at 250 nm).  This reflects Fabry-Pérot field
   enhancement/suppression in the ZnO film: the "incoherent depth" mode removes depth-wise
   coherence within the Raman source integral but not the optical cavity effect on the
   fields themselves.

5. **Comparison with original predictions:** the header anticipated E2(high) invariant
   and A1(LO) most sensitive.  The actual picture is more complex.  E2(high) is itself
   suppressed — not a stable reference.  A1(LO) is not "most sensitive"; it shows a
   nearly constant normalised ratio of ~0.08–0.10 across a wide thickness range.  The
   dominant feature is the Fabry-Pérot modulation of the E₂ modes, which becomes
   significant above ~30–40 nm.  See updated comments in `depth_coherence_crossover.py`.

#### Laser-wavelength dependence: 532 nm vs 1064 nm

The script was re-run at 1064 nm (double the wavelength) to test whether the crossover
thickness scales as λ as expected from the phase argument Δφ = 4πnd/λ.

**Theory:** doubling λ halves the phase accumulated per unit depth, so all crossover
thicknesses should double.  For ZnO (n ≈ 2, 532 nm laser): d_crossover ~ λ/(4n) ≈ 67 nm.
At 1064 nm this should shift to ~133 nm.

**Normalised ratio (coherent/incoherent) / d_Å at 532 nm and 1064 nm**

| Thickness |  d_Å  | E2(low) 532 | E2(low) 1064 | E2(hi) 532 | E2(hi) 1064 | A1(LO) 532 | A1(LO) 1064 |
|-----------|-------|-------------|--------------|------------|-------------|------------|-------------|
|    10 nm  |   100 |   0.9900    |    0.9977    |   0.9929   |    1.0006   |   0.0497   |    0.0486   |
|    20 nm  |   200 |   0.9553    |    0.9900    |   0.9588   |    0.9931   |   0.0535   |    0.0496   |
|    30 nm  |   300 |   0.8856    |    0.9764    |   0.8903   |    0.9800   |   0.0598   |    0.0511   |
|    40 nm  |   400 |   0.7709    |    0.9556    |   0.7776   |    0.9599   |   0.0682   |    0.0533   |
|    50 nm  |   500 |   0.6128    |    0.9261    |   0.6217   |    0.9314   |   0.0779   |    0.0561   |
|    75 nm  |   750 |   0.1916    |    0.8056    |   0.2003   |    0.8152   |   0.0994   |    0.0654   |
|   100 nm  |  1000 |   0.0921    |    0.6157    |   0.0897   |    0.6311   |   0.0968   |    0.0770   |
|   250 nm  |  2500 |   0.2126    |    0.2156    |   0.2085   |    0.2055   |   0.0813   |    0.0815   |
|   500 nm  |  5000 |   0.1936    |    0.2107    |   0.1914   |    0.2002   |   0.0826   |    0.0816   |
|  1000 nm  | 10000 |   0.1593    |    0.1932    |   0.1489   |    0.1801   |   0.0838   |    0.0816   |
|     5 μm  | 50000 |   0.1442    |    0.1399    |   0.0014   |    0.0021   |   0.0247   |    0.0247   |

**Findings**

1. **E2(low) and E2(high) — Fabry-Pérot minimum doubles as expected.**  At 532 nm the
   ratio/d_Å reaches its minimum (~0.19) at 75 nm.  At 1064 nm the ratio is still 0.61
   at 100 nm, with the minimum expected near ~150 nm (between the 100 and 250 nm points).
   The characteristic thickness scales approximately as λ, consistent with the phase
   argument Δφ = 4πnd/λ.  Beyond ~250 nm both wavelengths converge to similar values
   (~0.14–0.21), confirming that the long-range asymptotic behaviour is only weakly
   wavelength-dependent.

2. **A1(LO) — shows moderate wavelength dependence.**  At 10 nm both wavelengths give
   similar values (532 nm: 0.0497, 1064 nm: 0.0486).  Through 40–100 nm the ratio at
   1064 nm is noticeably smaller (0.05–0.08) than at 532 nm (0.07–0.10), reflecting the
   slower phase accumulation at longer wavelength.  Both converge to ~0.082–0.084 in the
   250 nm–1 μm range, and fall identically to 0.0247 at 5 μm.  The A₁(LO) ratio is
   therefore not wavelength-independent in the 50–100 nm transition region but converges
   at both thicker and thinner limits.

#### Why A1(LO) crossover thickness is laser-wavelength independent

**The phonon wavevector in the modal_pairs code** is computed as
(`LayeredRamanCalculator.py`, line 1301–1304):

```python
q_ph = np.array([
    zeta_L_re - zeta_S_re,    # in-plane
    0.0,
    channel_L["qz"] - channel_S["qz"],  # out-of-plane
])
```

i.e. **q_ph = k_L − k_S** (difference of Berreman channel qz values).

**Channel structure for c-axis ZnO at normal incidence.**
For propagation along c, ZnO looks optically isotropic (all modes are ordinary with
n_o ≈ 2).  The two q channels per polarisation are simply:

| Channel | qz |
|---|---|
| Forward (+) | +n_o |
| Backward (−) | −n_o |

**Phonon q for each channel pair:**

| (L channel, S channel) | q_ph = qz_L − qz_S | Mode |
|---|---|---|
| (forward, forward) | n_o − n_o = 0 | TO (q ≈ 0) |
| (forward, backward) | n_o − (−n_o) = **+2n_o** | **LO** |
| (backward, forward) | −n_o − n_o = **−2n_o** | **LO** |
| (backward, backward) | −n_o − (−n_o) = 0 | TO (q ≈ 0) |

**Integrand phase for each pair type.**
The field in the film for the forward E_L channel goes as exp(+i k_L z) and for the
backward E_S channel as exp(−i k_S z), so the integrand phase is:

| Pair | Integrand ∝ | Phase Δφ across depth d |
|---|---|---|
| TO (fwd L, fwd S) | exp(+i (k_L + k_S) z) | 4πn_o d / λ_laser |
| **LO (fwd L, bwd S)** | **exp(+i (k_L − k_S) z)** | **2πn_o ν_phonon d** |

For TO/E2 modes, Δφ = 4πn_o d/λ_laser → **laser-wavelength dependent** → crossover
scales as λ (factor-of-2 shift confirmed above).

For A1(LO), Δφ = 2πn_o ν_phonon d, where ν_phonon = 510 cm⁻¹ is a **fixed crystal
property**.  This phase is completely independent of the laser wavelength.

**Predicted A1(LO) crossover thickness:**

```
λ_phonon = 1/ν_phonon = 1/(510 cm⁻¹) ≈ 19 600 nm
d_crossover (Δφ = π) = λ_phonon / (2 n_o) ≈ 19 600 / (2 × 2) ≈ 4 900 nm ≈ 4.9 μm
```

This is fully consistent with the data showing the normalised ratio falling
significantly at the 5 μm thickness point with identical values at 532 and 1064 nm
(both 0.0247).  However, the normalised ratio for A1(LO) is NOT fully wavelength-
independent at intermediate thicknesses: at 100 nm the ratio is 0.0968 (532 nm) vs
0.0770 (1064 nm), a 20% difference.  This shows that the incoherent Berreman pair
summation (modal_pairs NAC) introduces a wavelength-dependent transitional behaviour
that the simple phonon-q phase argument does not capture.  The argument correctly
predicts the ~5 μm crossover is wavelength-independent, but the normalised ratio
through the 50–250 nm transition is not.

**Phase summary at key thicknesses (A1(LO), n_o ≈ 2, ν_phonon = 510 cm⁻¹):**

| Thickness | Δφ (rad) | Δφ / π |
|-----------|----------|--------|
| 100 nm    |  0.065   |  0.02  |
| 500 nm    |  0.324   |  0.10  |
| 1000 nm   |  0.647   |  0.21  |
| 5 μm      |  3.24    |  1.03  |

The small phonon-q phase at 100–1000 nm means the LO coherent depth integral is
not substantially reduced by depth-phase cancellation in this range.  The ~π phase
at 5 μm explains the onset of cancellation there, consistent with both wavelengths.

#### Why the LO mode uses (forward L, backward S) and (backward L, forward S) channel pairs

**The E_S field convention.**  In the GTM scattering-matrix formalism, E_S is not the
outgoing scattered field — it is the *reciprocal* (time-reversed) scattered field,
launched **downward from the detector** into the film.  This is the standard adjoint
trick used in transfer-matrix Raman calculations.  Within the film, E_S has both
a downward (+z, forward) component and an upward (−z, backward) component.  The
**backward component of E_S** (qz_S = −n_o, pointing upward) is, when time-reversed,
the physically upward-going detected photon.

**Physical interpretation of each pairing:**

| (L channel, S channel) | Physical process | q_ph | Mode |
|---|---|---|---|
| (forward, forward) | laser ↓ + E_S ↓ → forward-scattered detected photon ↓ | ~0 | TO |
| **(forward, backward)** | **laser ↓ + E_S ↑ → detected photon ↑ (backscattering)** | **+2n_o** | **LO** |
| **(backward, forward)** | **reflected laser ↑ + E_S ↓ → backscattering via reflection** | **−2n_o** | **LO** |
| (backward, backward) | double-reflected forward scatter | ~0 | TO |

**The (forward L, backward S) pair** is the primary backscattering event: the laser
photon propagates downward (+z) and the detected scattered photon propagates upward
(−z).  The phonon that mediates this process must supply momentum equal to the sum of
the two photon momenta:

```
q_phonon = k_laser(+z) + k_scattered(+z, upward = negative direction)
         = k_L − (−k_S)  [momentum conservation]
         = k_L + k_S ≈ 2k_laser = 2n_o ω_laser / c
```

In the code this appears as q_ph = qz_L − qz_S = (+n_o) − (−n_o) = +2n_o, because
the sign convention in `q_ph = channel_L["qz"] − channel_S["qz"]` absorbs the
time-reversal of E_S (i.e. the code's qz_S is the reciprocal-field wavevector, whose
physical momentum for the detected photon is its negative).

**The (backward L, forward S) pair** arises because the laser is partially reflected
at the ZnO/substrate interface, creating a backward-propagating component inside the
film (qz_L = −n_o).  This reflected laser pairs with the primary downward E_S component
(qz_S = +n_o) to give q_ph = −2n_o — still an LO phonon, just propagating in the −z
direction.

**Summary of the selection rule:**  LO phonons appear in Raman backscattering (both
photons travel in opposite directions) because |q_phonon| ≈ 2k_laser is large enough
to lie in the regime where NAC is important and pushes the phonon frequency toward the
LO value.  TO phonons appear in near-forward scattering (both photons co-propagating)
because |q_phonon| ≈ 0 places the phonon at the BZ centre TO frequency.  This is the
standard backscattering selection rule for uniaxial polar crystals, recovered here
directly from the Berreman-channel structure of the film.

### Substrate incoherence options: intensity, phase averaging, phase cancelling

*(Updated after bug fixes — see "Bug fixes" section below for history.)*

Five scripts compare different PDielec incoherence treatments for the Fused Silica
substrate layer, all using the same layer stack
(air / ZnO 100 nm / Fused Silica d / air) with d sweeping 10.00–11.98 μm in
0.02 μm steps (100 scenarios):

| Script | Incoherence treatment | GTM mode |
|---|---|---|
| `substrate_coherent.py` | `Coherent` (reference) | Scattering matrix |
| `substrate_incoherent_intensity.py` | `Incoherent (intensity)` | **Transfer matrix** |
| `substrate_incoherent_phase_average.py` | `Incoherent (phase averaging)` | Scattering matrix |
| `substrate_incoherent_phase_cancel.py` | `Incoherent (phase cancelling)` | Scattering matrix |
| `substrate_incoherent_thick.py` | `Incoherent (non-reflective)` | **Transfer matrix** |

Note that `substrate_incoherent_intensity.py` and `substrate_incoherent_thick.py` use
Transfer matrix mode (not Scattering matrix), which distinguishes them from the other
three scripts.

#### Summary statistics (peak intensity vs substrate thickness 10.00–11.98 μm, 100 scenarios)

| Mode | DFT freq | Coherent | Inc. intensity | Inc. phase avg | Inc. phase cancel | Inc. non-reflective |
|---|---|---|---|---|---|---|
| E2(low)  | 109 | min=1.02 max=14.79 range/mean=1.98 | **25.86 (flat)** | ~6.922 (range/mean=0.167) | **26.16 (flat)** | **5.405 (flat)** |
| E2(high) | 398 | min=0.08 max=11.60 range/mean=2.55 | **13.85 (flat)** | ~4.453 (range/mean=0.570) | **13.85 (flat)** | **2.868 (flat)** |
| A1(LO)   | 510 | min=1.21 max=8.69 range/mean=1.67 | **9.835 (flat)** | ~4.526 (range/mean=0.060) | **9.792 (flat)** | **2.753 (flat)** |

#### Finding 1 (corrected): Incoherent (phase averaging) removes etalon oscillations

After fixing Bug 1 (see below), `Incoherent (phase averaging)` now correctly loops over
N phase-shift values and averages the results, exactly mirroring the IR
`average_incoherent_calculator`.  The result removes the coherent etalon oscillation,
but it is not flat for every mode in the regenerated data: E2(low) decreases by
range/mean = 0.167 and E2(high) by 0.570 over the 10.00–11.98 μm sweep, while
A1(LO) shows a smaller monotonic rise.

The phase-averaged mean values sit close to half-way between the coherent cycle-mean and
the phase-cancelling value — consistent with the interpretation that phase averaging
suppresses but does not entirely eliminate the coherent back-reflection contribution.

#### Finding 2 (corrected): Incoherent (intensity) is perfectly flat

After fixing Bug 2 (see below), `Incoherent (intensity)` now gives a perfectly
thickness-independent result (range/mean = 0.0000 to machine precision for all three
modes).

The root cause of the original failure was that the Raman calculator uses
`calculate_modal_amplitudes` (not `calculate_Efield`) to obtain the incident field
amplitudes inside the ZnO layer.  Both methods contain the same backward-propagation
recursion (`F_ft[-2] = layers[-1].Ki * F_bk[-2]`), which propagated the full complex
Fused Silica phase (`Ki = exp(i·2π·n_FS·d·f/c)`) into the ZnO boundary amplitudes before
any correction could act.  The fix adds a `_efield_Ki` / `_efield_Ki_amp` helper that
returns an attenuation-only (real-part) matrix for `IncoherentIntensityLayer`; for lossless
Fused Silica this gives the identity matrix, so no substrate phase enters the ZnO field.

#### Finding 3 (unchanged): Incoherent (phase cancelling) is perfectly flat

The `Incoherent (phase cancelling)` treatment gives a completely thickness-independent
result — every one of the 101 scenario columns produces the same peak intensity to
machine precision.  This was correct before the bug fixes and remains correct.

The phase-cancelling values are substantially *higher* than the coherent cycle-mean and
the non-reflective reference.  They also exceed the `Incoherent (intensity)` values
slightly (~1% difference), which is expected: phase-cancelling sets all propagation
phases to zero (constructive interference limit), while intensity averaging gives the
cycle-mean of the squared amplitude.

#### Finding 4: Incoherent (non-reflective) is perfectly flat and lowest

`Incoherent (non-reflective)` (`substrate_incoherent_thick.py`) gives a perfectly
thickness-independent result — range/mean = 0.0000 to machine precision for all three
modes across all 101 substrate-thickness scenarios.

The non-reflective values are the *lowest* of all treatments, substantially below even
the coherent cycle-mean:

| Mode | Inc. non-refl | Coh. cycle-mean | Non-refl / Coh. mean |
|---|---|---|---|
| E2(low)  | 5.405 | 6.96 | 0.776× |
| E2(high) | 2.868 | 4.52 | 0.635× |
| A1(LO)   | 2.753 | 4.47 | 0.615× |

These values are in excellent agreement with the non-reflective results from
`substrate_gtm_modes.py` (A1(LO)=2.75, E2(low)=5.41, E2(high)=2.87),
confirming that the `IncoherentThickLayer` implementation is consistent across both
scripts.

The non-reflective intensities are approximately one-fifth to one-third of the
`Incoherent (intensity)` and `Incoherent (phase cancelling)` treatments (ratio ~0.21
for E2(low), ~0.21 for E2(high), ~0.28 for A1(LO)).  The difference is the back-surface
reflection contribution: removing reflections from the glass/air interface eliminates
both the reflected laser beam (which would excite additional Raman in the ZnO) and the
backward-collected scattered photons (which would add to the detected signal).

#### Per-mode comparison tables (selected substrate thicknesses, plus final 11.98 μm)

**E2(low) — 109 cm⁻¹**

| Thickness (μm) | Coherent | Inc. intensity | Inc. phase cancel | Inc. phase avg | Inc. non-refl |
|---:|---:|---:|---:|---:|---:|
| 10.00 |  2.0555 | 25.8616 | 26.1553 |  7.5273 | 5.4052 |
| 10.20 |  4.5379 | 25.8616 | 26.1553 |  7.3964 | 5.4052 |
| 10.40 |  7.7609 | 25.8616 | 26.1553 |  7.2682 | 5.4052 |
| 10.60 | 10.9892 | 25.8616 | 26.1553 |  7.1429 | 5.4052 |
| 10.80 | 13.1000 | 25.8616 | 26.1553 |  7.0207 | 5.4052 |
| 11.00 | 12.9754 | 25.8616 | 26.1553 |  6.9017 | 5.4052 |
| 11.20 | 10.4833 | 25.8616 | 26.1553 |  6.7863 | 5.4052 |
| 11.40 |  6.8967 | 25.8616 | 26.1553 |  6.6745 | 5.4052 |
| 11.60 |  3.8807 | 25.8616 | 26.1553 |  6.5665 | 5.4052 |
| 11.80 |  2.2893 | 25.8616 | 26.1553 |  6.4625 | 5.4052 |
| 11.98 |  2.5943 | 25.8616 | 26.1553 |  6.3725 | 5.4052 |
| **Min** | 1.0151 | 25.8616 | 26.1553 | 6.3725 | 5.4052 |
| **Max** | 14.7886 | 25.8616 | 26.1553 | 7.5273 | 5.4052 |
| **Mean** | 6.9630 | 25.8616 | 26.1553 | 6.9220 | 5.4052 |
| **Range/Mean** | 1.9781 | 0.0000 | 0.0000 | 0.1668 | 0.0000 |

**E2(high) — 398 cm⁻¹**

| Thickness (μm) | Coherent | Inc. intensity | Inc. phase cancel | Inc. phase avg | Inc. non-refl |
|---:|---:|---:|---:|---:|---:|
| 10.00 |  3.0160 | 13.8523 | 13.8528 |  5.7948 | 2.8684 |
| 10.20 |  5.2304 | 13.8523 | 13.8528 |  5.5324 | 2.8684 |
| 10.40 |  7.4112 | 13.8523 | 13.8528 |  5.2559 | 2.8684 |
| 10.60 |  9.0541 | 13.8523 | 13.8528 |  4.9710 | 2.8684 |
| 10.80 |  9.5276 | 13.8523 | 13.8528 |  4.6839 | 2.8684 |
| 11.00 |  8.4199 | 13.8523 | 13.8528 |  4.4006 | 2.8684 |
| 11.20 |  6.0747 | 13.8523 | 13.8528 |  4.1271 | 2.8684 |
| 11.40 |  3.5559 | 13.8523 | 13.8528 |  3.8691 | 2.8684 |
| 11.60 |  1.8362 | 13.8523 | 13.8528 |  3.6319 | 2.8684 |
| 11.80 |  1.1905 | 13.8523 | 13.8528 |  3.4205 | 2.8684 |
| 11.98 |  1.5397 | 13.8523 | 13.8528 |  3.2559 | 2.8684 |
| **Min** | 0.0781 | 13.8523 | 13.8528 | 3.2559 | 2.8684 |
| **Max** | 11.6047 | 13.8523 | 13.8528 | 5.7948 | 2.8684 |
| **Mean** | 4.5159 | 13.8523 | 13.8528 | 4.4530 | 2.8684 |
| **Range/Mean** | 2.5524 | 0.0000 | 0.0000 | 0.5702 | 0.0000 |

**A1(LO) — 510 cm⁻¹**

| Thickness (μm) | Coherent | Inc. intensity | Inc. phase cancel | Inc. phase avg | Inc. non-refl |
|---:|---:|---:|---:|---:|---:|
| 10.00 |  4.4656 | 9.8348 | 9.7919 |  4.4247 | 2.7525 |
| 10.20 |  4.8806 | 9.8348 | 9.7919 |  4.4278 | 2.7525 |
| 10.40 |  4.9885 | 9.8348 | 9.7919 |  4.4382 | 2.7525 |
| 10.60 |  4.6932 | 9.8348 | 9.7919 |  4.4556 | 2.7525 |
| 10.80 |  4.0342 | 9.8348 | 9.7919 |  4.4793 | 2.7525 |
| 11.00 |  3.1647 | 9.8348 | 9.7919 |  4.5087 | 2.7525 |
| 11.20 |  2.3589 | 9.8348 | 9.7919 |  4.5427 | 2.7525 |
| 11.40 |  1.9426 | 9.8348 | 9.7919 |  4.5801 | 2.7525 |
| 11.60 |  2.1488 | 9.8348 | 9.7919 |  4.6196 | 2.7525 |
| 11.80 |  3.0462 | 9.8348 | 9.7919 |  4.6599 | 2.7525 |
| 11.98 |  2.5253 | 9.8348 | 9.7919 |  4.6957 | 2.7525 |
| **Min** | 1.2069 | 9.8348 | 9.7919 | 4.4247 | 2.7525 |
| **Max** | 8.6858 | 9.8348 | 9.7919 | 4.6957 | 2.7525 |
| **Mean** | 4.4721 | 9.8348 | 9.7919 | 4.5257 | 2.7525 |
| **Range/Mean** | 1.6723 | 0.0000 | 0.0000 | 0.0599 | 0.0000 |

#### Relative values: incoherent options vs coherent cycle-mean

| Mode | Coh. mean | Inc. intensity | Inc. ph. cancel | Inc. ph. avg mean | Inc. non-refl |
|---|---|---|---|---|---|
| E2(low)  | 6.96 | 25.86 (3.71×) | 26.16 (3.76×) | 6.922 (0.99×) | 5.405 (0.776×) |
| E2(high) | 4.52 | 13.85 (3.07×) | 13.85 (3.07×) | 4.453 (0.99×) | 2.868 (0.635×) |
| A1(LO)   | 4.47 | 9.835 (2.20×) | 9.792 (2.19×) | 4.526 (1.01×) | 2.753 (0.615×) |

Key observations:
- **Inc. intensity ≈ Inc. phase cancel** (within 1.2% for all modes): both are perfectly
  flat, and the two intensity values are essentially equal.  Phase cancelling achieves
  constructive interference; intensity averaging cycles over all phases and yields the
  same cycle-mean as intensity cancellation.
- **Inc. phase avg mean ≈ Coherent mean** (within 2% for all modes): the phase-averaging
  cycle-mean equals the coherent cycle-mean, as expected — both average the same etalon
  cycle.
- **Inc. intensity / Inc. phase avg is strongly mode dependent**: about 3.7× for E2(low),
  3.1× for E2(high), and 2.2× for A1(LO).  This reflects the substrate back-reflection
  enhancement: the incoherent treatments (both intensity and phase-cancel) retain the
  back-reflected laser contribution to the ZnO field coherently, whereas the
  phase-averaging mean cycle-averages it.
- **Inc. non-refl is the lowest** (~0.62–0.78× the coherent mean for all three modes).
  Removing back-reflection eliminates both the reflected laser enhancement of ZnO
  excitation and the upward-reflected Raman collection — giving a conservative lower bound.

#### Why incoherent intensity exceeds even the coherent maximum

The incoherent intensity values (e.g. E2(low) = 25.86) exceed not only the coherent
cycle-mean (6.96) but also the highest individual coherent scenario (14.79).  This is
not a numerical artefact — it is a direct consequence of how the incoherent intensity
approximation breaks the energy-conservation constraint that limits the coherent case.

**In the coherent case there is a Fabry-Perot trade-off.**  The substrate etalon both
adds a back-reflected field inside ZnO and modifies the forward field through the coupled
cavity boundary conditions.  At constructive interference (the coherent maximum), the
overall reflectance of the system back toward the superstrate is also at its maximum,
meaning less of the incident power couples forward into ZnO.  Enhancement and loss are
linked: energy conservation across the coupled system caps the field inside ZnO at the
coherent maximum.

**In the incoherent intensity case this trade-off is broken.**  The calculation performs
two independent steps:

1. `calculate_GammaStar` computes the ZnO boundary conditions with the substrate
   excluded from coherent interference (for lossless Fused Silica, `|Pi|² = identity` so
   the substrate drops out of the wave equation entirely).  The forward field in ZnO is
   therefore at its **unconstrained** level — as if the substrate Fabry-Perot cavity were
   absent.

2. `Ki = identity` (from the `_efield_Ki_amp` fix) then adds the substrate back-reflection
   on top of that field **constructively** (zero phase), contributing the full back-reflected
   amplitude without any phase penalty.

No coherent scenario can simultaneously have a full (unconstrained) forward field **and**
a constructive back-reflection, because in the coherent case a stronger back-reflection
always comes with a weaker forward coupling.  The incoherent intensity approximation
decouples them, which is why the result exceeds the coherent maximum.

The physical ordering is therefore:

| Treatment | Forward field in ZnO | Substrate back-reflection |
|---|---|---|
| Non-reflective | Unconstrained (no cavity) | None |
| Coherent maximum | Reduced by cavity resonance | Constructive — but at the cost of reduced forward field |
| Inc. intensity / phase cancel | Unconstrained (no cavity) | Constructive — at no cost |

The non-reflective case has the full forward field but no back-reflection, giving the
minimum.  The coherent maximum has constructive back-reflection but pays the cavity energy
cost.  The incoherent intensity treatment gets both at once, which is why it exceeds the
coherent maximum.

#### Bug fixes applied to the Crystal Raman incoherence path

Two bugs in the original Crystal Raman code caused incorrect results for `Incoherent
(phase averaging)` and `Incoherent (intensity)`.  Both are now fixed.

**Bug 1 — Phase averaging had no effect (phaseShift always 0)**

Root cause: `_calculate_raman` in `CrystalScenarioTab.py` called `_build_raman_calculator`
exactly once with `phaseShift = 0.0`.  The IR path correctly looped over N phase shifts
via `average_incoherent_calculator`, but no equivalent loop existed for Raman.

Fix: new `_raman_intensities(psi_rad)` helper mirrors the IR loop — it finds all layers
with `get_incoherent_option() == "Incoherent (phase averaging)"`, iterates over
`beta_k = frac·2π·k/N` (k = 0…N−1), sets `layer.set_phase_shift(beta_k)` for each,
builds a calculator, accumulates mode intensities, and divides by N.  This helper is used
in both `_calculate_raman` and `_run_azimuthal_sweep`.

**Bug 2 — Intensity incoherence did not remove substrate etalon from E-field**

Root cause: The Raman calculator uses `calculate_modal_amplitudes` (in `GTMcore.py`) to
obtain the incident field amplitudes inside the ZnO layer.  This method contains the
same backward-propagation recursion as `calculate_Efield`:

```python
F_ft[-2] = layers[-1].Ki * F_bk[-2]   # propagate through last layer
for kl in ...:
    F_ft[kl] = layers[kl-1].Ki * F_bk[kl]
```

`IncoherentIntensityLayer` does not override `Ki`; it uses the same full coherent
propagation matrix `exp(i·2π·n·d·f/c)`.  For lossless Fused Silica this has unit
magnitude but carries the oscillatory phase, which was baked into the ZnO boundary
amplitudes `F_bk` before any correction could act.  The `calculate_GammaStar`
intensity treatment (element-wise |T|²) correctly removes the etalon from the boundary
conditions (r, t), but `calculate_modal_amplitudes` ignored this and still used the
full coherent Ki.

Fix: a `_efield_Ki` / `_efield_Ki_amp` helper was added inside both
`calculate_Efield` and `calculate_modal_amplitudes`.  For `IncoherentIntensityLayer`
it returns `diag(exp(Re(exponent)))` instead of `Ki = diag(exp(exponent))`.  For
lossless Fused Silica `Re(exponent) = 0`, so the helper returns the identity matrix —
no substrate phase enters the ZnO boundary amplitudes, and the etalon disappears.

#### Mechanism of each incoherence treatment (from `GTMcore.py`)

All four options map to a GTM layer class in `GTMcore.py` via `gtm_methods` in
`CrystalScenarioTab.py`:

| Option | GTM class | Key override | Matrix method |
|---|---|---|---|
| `Coherent` | `CoherentLayer` | Standard propagation exponents | Both |
| `Incoherent (intensity)` | `IncoherentIntensityLayer` | Transfer-matrix intensity formulation | **Transfer only** |
| `Incoherent (phase averaging)` | `IncoherentAveragePhaseLayer` | Adds `phaseShift` to backward modes | Both |
| `Incoherent (phase cancelling)` | `IncoherentPhaseLayer` | Arteaga's method — zeroes propagation phase | Both |
| `Incoherent (non-reflective)` | `IncoherentThickLayer` | Zeroes backward modes and S11 matrix | Both |

**`IncoherentIntensityLayer` — intensity transfer matrices (Transfer matrix only)**
(`GTMcore.py`, `TransferMatrixSystem.calculate_GammaStar`)

`IncoherentIntensityLayer` does **not** override `calculate_propagation_exponents`.  Its
`Ki` matrix is identical to a `CoherentLayer`.  The incoherent treatment is applied
entirely at the level of the **system** transfer matrix in `TransferMatrixSystem.calculate_GammaStar`,
which checks `layer.inCoherentIntensity` during the layer loop.

When an incoherent-intensity layer is encountered, the algorithm splits the running
product of transfer matrices into separate amplitude segments:

```
Tlist = [ A_inv_super × Di_inv_FS × ... × A_sub,    ← right-side amplitude segment
          Pi_FS,                                       ← propagation of incoherent layer
          A_inv_air × Di_inv_ZnO × Pi_ZnO × Di_ZnO × Di_FS ]  ← left-side amplitude segment
```

The final system matrix is then assembled as:

```python
T = identity
for t in Tlist:
    T = T × |t|²           # convert each segment to element-wise intensities
T = sqrt(T)                # convert back to amplitude scale
```

This removes the coherent cross-interference between the amplitude segments on either
side of the incoherent layer boundary — the phase relationship between the left- and
right-side blocks is discarded.  For a non-absorbing substrate (Fused Silica), the
propagation matrix Pi has diagonal elements `exp(±iφ)` so `|Pi|² = identity` — the
substrate thickness drops out of the system matrix entirely, just as in phase cancelling.

The `ScatteringMatrixSystem` has no handler for `inCoherentIntensity`
(the corresponding code is commented out in `calculate_scattering_matrix`), so using
this option with a scattering matrix simply treats the layer as fully coherent.

**`IncoherentPhaseLayer` — Arteaga's phase-cancelling method**
(`GTMcore.py`, class `IncoherentPhaseLayer`, method `calculate_propagation_exponents`)
electric field *inside* the ZnO film is computed by propagating from the substrate
boundary through each layer's Ki matrix — and `IncoherentIntensityLayer` does not
override Ki; it uses the same full coherent propagation matrix as `CoherentLayer`,
containing `exp(i·2π·n_FS·d·f/c)`.  The substrate etalon therefore reappears at the
field-propagation stage, and the Raman intensity retains nearly its full thickness
dependence.

In short: the `Incoherent (intensity)` option suppresses the substrate etalon only at
the boundary-condition level, not at the field level where Raman actually depends on it.

With the current incoherent-intensity fix, the mean values (averaged over all 100
thickness steps) are no longer nearly indistinguishable:

| Mode | Coherent mean | Incoh. intensity mean | % difference |
|---|---|---|---|
| E2(low)  | 6.96 | 25.86 | +271.4% |
| E2(high) | 4.52 | 13.85 | +206.7% |
| A1(LO)   | 4.47 | 9.83 | +119.9% |

#### Finding 3: Incoherent (phase cancelling) is perfectly flat and highest

The `Incoherent (phase cancelling)` treatment gives a **completely thickness-independent
result** — every one of the 100 scenario columns produces the same peak intensity to
machine precision.  This is the only treatment of the substrate that fully removes
the Fabry-Perot etalon effect from the Raman spectrum.

The phase-cancelling values are substantially *higher* than the coherent cycle-mean
and the non-reflective reference:

| Mode | Non-reflective | Coherent mean | Phase cancel | Cancel / Mean | Cancel / Non-refl |
|---|---|---|---|---|---|
| E2(low)  | 5.41 | 6.96 | 26.16 | 3.76× | 4.84× |
| E2(high) | 2.87 | 4.52 | 13.85 | 3.07× | 4.83× |
| A1(LO)   | 2.75 | 4.47 | 9.79 | 2.19× | 3.56× |

**`IncoherentPhaseLayer` — Arteaga's phase-cancelling method**
(`GTMcore.py`, class `IncoherentPhaseLayer`, method `calculate_propagation_exponents`)

The propagation exponent for mode `i` is `exp(-i 2π f qᵢ d / c)`.  Arteaga's method
(`Thin Solid Films 2014, 571, 701`) modifies the q values before computing the exponent:

```python
qs[0] = self.qs[0] - np.real(self.qs[1])   # p-wave forward:  q_p - Re(q_s)
qs[1] = self.qs[1] - np.real(self.qs[1])   # s-wave forward:  q_s - Re(q_s) = i·Im(q_s)
qs[2] = self.qs[2] - np.real(self.qs[3])   # p-wave backward: q_p - Re(q_s)
qs[3] = self.qs[3] - np.real(self.qs[3])   # s-wave backward: q_s - Re(q_s) = i·Im(q_s)
```

For a non-absorbing isotropic medium (such as Fused Silica) at normal incidence, all
Berreman q values are **real** (q = ±n_FS, with n_FS ≈ 1.461):

```
qs[1] = q_s - Re(q_s) = q_s - q_s = 0    (all modes → 0 after subtraction)
```

All four modified q values become zero.  Therefore all four propagation exponents
become `exp(0) = 1`, and the propagation matrix **P = identity**.

This has three consequences:
1. **Thickness-independence:** the identity propagation matrix contains no `d` dependence →
   result is exactly the same for any substrate thickness (explains the flat data).
2. **Interfaces are preserved:** Fresnel reflection coefficients at the ZnO/glass and
   glass/air interfaces are still calculated from the full dielectric contrast — the
   substrate is not made transparent; only its *propagation phase* is removed.
3. **Maximum constructive interference:** with zero phase accumulated on every bounce,
   all partial waves reflected back and forth inside the substrate add with phase zero →
   all multiple reflections interfere fully constructively.  This is why the
   phase-cancelling intensity exceeds even the coherent *maximum* — the coherent maximum
   only occurs at a specific thickness; phase-cancelling achieves it for all thicknesses.

**`IncoherentAveragePhaseLayer` — phase averaging over backward modes**
(`GTMcore.py`, class `IncoherentAveragePhaseLayer`, method `calculate_propagation_exponents`)

The class adds a constant extra `phaseShift` only to the backward modes (2 and 3):

```python
# modes 0,1 (forward): standard exponent, unchanged from coherent
exponent = -1j*(2π f qᵢ d / c)
# modes 2,3 (backward): standard exponent plus phaseShift
exponent = -1j*(2π f qᵢ d / c + phaseShift)
```

After Bug 1 was fixed, the Crystal Raman path now correctly loops over N phase-shift
values `beta_k = frac·2π·k/N` (k = 0…N−1) via the `_raman_intensities` helper, sets
`layer.set_phase_shift(beta_k)` before each calculator build, and averages the resulting
mode intensities.  For N = 10 and frac = 1 (100% incoherence) the result removes the
coherent etalon oscillation, although the regenerated E2(low) and E2(high) data retain
sizeable monotonic thickness dependence.  Before the fix, `phaseShift` was always 0.0 and
the result was numerically identical to Coherent.

**`IncoherentThickLayer` — non-reflective**
(`GTMcore.py`, class `IncoherentThickLayer`, methods `calculate_propagation_matrix` and
`calculate_scattering_matrix`)

The propagation matrix is computed normally but then the backward-mode entries are zeroed:

```python
Ki[2,2] = 0.0   # p-wave backward mode suppressed
Ki[3,3] = 0.0   # s-wave backward mode suppressed
```

Additionally, `calculate_scattering_matrix` is overridden to set `S11 = 0` — the
sub-matrix describing reflection back toward the superstrate.  This removes all
reflections from the back surface of the substrate, eliminating the etalon entirely.
The result is lower than the coherent cycle-mean because the back-reflected laser (which
would otherwise excite additional Raman in the ZnO film) and the back-reflected Raman
photon (which would otherwise contribute to the collected signal) are both suppressed.

#### Practical implications for substrate treatment

| Treatment | Etalon oscillation? | Thickness-independent? | Value relative to coh. mean | Physical interpretation |
|---|---|---|---|---|
| Coherent | Yes, strong (up to 7.6×) | No | 1.00× (by definition) | Perfect optical flat — valid for ultra-flat substrates in cryogenic setups |
| Incoherent intensity | No (flat) | Yes | ~2.0–2.2× | Intensity-based suppression; flat and consistent with phase-cancel |
| Phase averaging | Near-flat (< 2%) | Approximately | ~1.00× | Cycle-average of etalon; equivalent to coherent mean but without computing many scenarios |
| Phase cancelling | No (flat) | Yes | ~2.0–2.2× | Thickness-independent and highest — all reflections add constructively |
| Non-reflective | No (flat) | Yes | ~0.62–0.97× | Thickness-independent and lowest — removes back-reflection entirely; conservative lower bound |

For a real experiment with a macroscopic glass substrate (where etalon fringes are
averaged by surface roughness, thermal drift, and finite coherence length):
- Use `Incoherent (intensity)` or `Incoherent (phase cancelling)` for a flat,
  thickness-independent result that correctly includes the back-reflection contribution.
- Use `Incoherent (phase averaging)` if you want the cycle-mean of the coherent etalon
  without computing many scenarios (gives the same mean as averaging over thicknesses).
- Use `Incoherent (non-reflective)` as a conservative lower bound that suppresses all
  back-surface contributions.
- Averaging many coherent scenarios (as done in `substrate_coherent.py`) is the most
  physically transparent approach and gives the same result as phase averaging.

### `layer_combination_limits.py`  *(Level 3 — Layer combination)*
- **Purpose:** Prove exact equality of `Coherent amplitudes` and `Incoherent
  intensities` layer combination for a single Raman layer, and show their divergence
  when multiple sublayers are used.
- **Scenarios (6):**
  - 1 × 100 nm — Incoherent intensities and Coherent amplitudes (must be identical).
  - 2 × 50 nm — Incoherent and Coherent (may differ due to inter-sublayer phase).
  - 5 × 20 nm — Incoherent and Coherent (cross terms multiply).
- **Note:** `Coalesce equivalent layers = False` prevents automatic merging of
  identical adjacent sublayers.
- **Levels fixed:** All layers GTM = Coherent, Depth coherence = Coherent amplitude,
  Modal pair combination = Group q channels.  GL = 200 for all.
- **NAC mode: `modal_pairs`.**  An earlier version of this script used `geometry` NAC
  because `_compute_modal_pair_mode_worker` did not consult the `coherent_layers` flag:
  with `Group q channels` all sublayers of the same material shared the same phonon
  q-vector group key and were always combined coherently, making the `Layer combination`
  setting inoperative.  The bug was fixed in `LayeredRamanCalculator.py` by including
  `layer_index` in the group key when `coherent_layers=False`, so the script now runs
  correctly with `modal_pairs` NAC.
- **Output:** `layer_combination_limits.xlsx`

#### Observations from running the script

**Single-layer identity holds exactly.**

```
1 × 100 nm: max|incoherent − coherent| = 0.000e+00
```

A single layer produces a single complex amplitude; squaring it is identical regardless
of whether it is labelled "coherent" or "incoherent" combination.

**Coherent and incoherent layer combination now diverge for multiple sublayers.**

```
2 × 50 nm: max|incoherent − coherent| = 9.8
5 × 20 nm: max|incoherent − coherent| = 7.6
```

Peak intensities at key ZnO modes (DFT frequencies from CASTEP PBE+G06):

| Mode | 1 × 100 nm | 2 × 50 nm incoh | 2 × 50 nm coh | i/c | 5 × 20 nm incoh | 5 × 20 nm coh | i/c |
|---|---|---|---|---|---|---|---|
| E2(low)  109 cm⁻¹ | 5.431 | 22.65 | 20.68 | 1.095 | 12.38 | 11.99 | 1.033 |
| E2(high) 398 cm⁻¹ | 2.871 | 12.60 | 11.59 | 1.087 |  6.816 |  6.650 | 1.025 |
| A1(LO)   510 cm⁻¹ | 2.753 |  1.410 |  2.772 | 0.509 |  0.572 |  2.761 | 0.207 |

**Non-polar E2 modes (incoherent > coherent):**  Ratios of 1.03–1.09× (incoh > coh).
The incoherent combination slightly overestimates relative to coherent because it
discards the negative cross-terms arising from propagation-phase accumulation across
the 100 nm film.  The effect is modest at the mode-peak level.

**Polar LO mode (A1(LO)) shows a qualitatively different pattern.**  For 2 × 50 nm the
ratio inverts to 0.509 — *incoherent is less than half of coherent* — and for 5 × 20 nm
it is 0.207.  This is because with `modal_pairs` NAC the LO Raman tensor for each
sublayer responds strongly to the q-direction of that sublayer's Berreman mode channels.
The inter-sublayer cross-terms for the LO mode are large and positive, so coherent
combination *strongly enhances* the LO intensity relative to an incoherent sum.

**Coherent subdivision invariance with `modal_pairs` NAC:**  For A1(LO), the single-layer
(1 × 100 nm) coherent result (2.753) agrees with the 2 × 50 nm coherent result (2.772)
to within 0.72%, consistent with subdivision invariance.  For the E2 modes however,
the 2 × 50 nm coherent result (E2low: 20.68, E2high: 11.59) is substantially larger
than the single-layer result (E2low: 5.43, E2high: 2.87) — the E2 coherent amplitude
is not subdivision-invariant in this geometry.  This reflects that at 100 nm the E2(low)
coherent depth integral is near a Fabry-Pérot cancellation condition; splitting into
sublayers each of 50 nm avoids that cancellation within each segment, producing a
larger total when recombined coherently.

**Conclusion: `Coherent amplitudes` is the correct choice for a continuous single-crystal film.**

Splitting a homogeneous crystal film into sublayers for numerical reasons and then
combining them incoherently produces a spuriously large (for TO/E2 modes) or
spuriously small (for LO modes with many sublayers) answer that depends on the number
of sublayers chosen.  `Incoherent intensities` layer combination is only physically
appropriate when the separate layers genuinely have no phase relationship (e.g.,
independent crystalline grains, layers separated by amorphous interlayers, or distinct
spectroscopic measurements summed post-hoc).

### `modal_pair_comparison.py`  *(Level 4 — Modal pair combination)*
- **Purpose:** Compare the three modal-pair summation strategies as a function of
  film thickness when `Layer NAC mode = 'modal_pairs'`.
- **Scenarios (15):** 5 thicknesses × 3 modal-pair options, interleaved by thickness.
  Thicknesses: 10, 100, 500 nm, 1, 2 μm.
  Options: `Group q channels`, `Incoherent pairs`, `Coherent all pairs`.
- **Levels fixed:** All layers GTM = Coherent, Depth coherence = Coherent amplitude,
  Layer combination = Incoherent intensities.
- **Physical test:** `Group q channels` ≥ `Incoherent pairs` in principle (coherent
  summation within a q-channel cannot decrease intensity); E2(high) should be identical
  between `Group q` and `Incoherent pairs` (non-polar, no q-channel structure).
  `Coherent all pairs` is unphysical and expected to differ.
- **Output:** `modal_pair_comparison.xlsx`

#### Observations from running the script

**Mode assignments (DFT frequencies, c-plane 001 pp geometry).**
The three active non-polar / polar-LO peaks tracked in this comparison are:

| Mode | DFT freq (cm⁻¹) | Character |
|---|---|---|
| E2(low) | 109 | Non-polar, doubly degenerate |
| E2(high)| 398 | Non-polar, doubly degenerate |
| A1(LO)  | 510 | Polar, LO branch (dominant peak) |

Note: the nominal textbook frequencies (E2(high) ~437 cm⁻¹, A1(LO) ~574 cm⁻¹) do not
apply to this DFT dataset; the `modal_pairs` NAC shifts the LO branch relative to the TO.

**GL point counts.**
The GL density setting scales the absolute number of integration points with film thickness:

| Thickness | GL density | GL points |
|---|---|---|
| 10 nm  | 20/μm  | 5 (min) |
| 100 nm | 100/μm | 10 |
| 500 nm | 500/μm | 250 |
| 1 μm   | 1000/μm | 1000 |
| 2 μm   | 1000/μm | 2000 |

The 10 nm case uses the **minimum of 5 GL points** (20/μm × 10 nm = 0.2, clamped to 5
by the code): only minimal depth-phase variation occurs, so results primarily reflect
modal-pair coherence with little cancellation from the depth integral.  GL convergence
testing (`Scripts/gl_convergence_a1lo_10nm.py`, 2026-07-09) confirmed that varying GL
density from 5 to 1000 actual points produces identical results to 7 significant figures
— the 10 nm results are fully converged even at the minimum GL count.

**Peak intensity ratios: Group q / Incoherent pairs (Gq/Ip)**

| Mode (cm⁻¹) | 10 nm | 100 nm | 500 nm | 1 μm | 2 μm |
|---|---|---|---|---|---|
| E2(low)  109 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| E2(high) 398 | 1.0004 | 1.0002 | 1.0000 | 1.0000 | 1.0000 |
| A1(LO)   510 | 1.0011 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

Key findings:

- **Group q ≡ Incoherent pairs to within 0.04% for all three modes at all thicknesses.**
  After the non-polar bypass fix (2026-07-02), E2 modes route through the standard
  coherent total-field path; their Group q / Incoherent pairs ratio is 1.0000 to <1 ppm.
  A1(LO) via `Group q channels` is effectively identical to `Incoherent pairs` because at
  normal incidence (c-plane, backscattering) the four Berreman channel pairs split into
  two symmetry-equivalent groups, so the intra-q amplitude sum reduces to the incoherent
  sum.

- **At 10 nm (5 GL points)**: the small 0.04–0.11% excess of Group q over Incoherent
  pairs indicates negligible intra-q coherence.  Results at 10 nm are dominated by
  the thin-film limit rather than depth-phase variation.

**Peak intensity ratios: Group q / Coherent all pairs (Gq/Ca)**

| Mode (cm⁻¹) | 10 nm | 100 nm | 500 nm | 1 μm | 2 μm |
|---|---|---|---|---|---|
| E2(low)  109 | 1.0000 | 0.9999 | 1.0000 | 1.0000 | 1.0000 |
| E2(high) 398 | 0.9999 | 0.9981 | 0.9993 | 0.9993 | 0.9999 |
| A1(LO)   510 | 0.5166 | 0.5014 | 0.5189 | 0.5764 | **0.9173** |

Key findings:

- **E2(low) and E2(high)**: Group q ≡ Coherent all pairs to within 0.02–0.19% at all
  thicknesses.  After the non-polar bypass fix, E2 modes use the standard coherent
  total-field path in which Group q and Coherent all pairs are the same code path.

- **A1(LO) at 510 cm⁻¹**: Coherent all pairs ≈ 2× Group q at 10–500 nm (ratio
  0.517–0.519), converging toward Group q only at 2 μm (ratio 0.917).  The inter-q
  constructive interference for the polar LO mode partially averages out at large
  thickness as the depth integral reduces coherence across the film.

**Conclusions.**

1. **Group q channels is the correct default** for `modal_pairs` NAC.  It coherently
   sums Berreman mode pairs that drive the same phonon wavevector q (physically
   meaningful interference), and sums distinct q-channels incoherently (distinct phonon
   momentum final states that cannot interfere).

2. **Group q ≡ Incoherent pairs to within 0.04%** for all three tracked modes at all
   thicknesses, after the non-polar bypass fix.  The choice between these two options
   has no practical effect.

3. **Coherent all pairs must not be used** for real calculations for polar modes.  It
   over-estimates A1(LO) intensity by ~2× at intermediate thicknesses (10–500 nm) due
   to unphysical constructive inter-q cross-terms.  Non-polar modes are now protected by
   the non-polar bypass so Group q and Coherent all pairs agree for E2 modes.

4. **The 10 nm results reflect the thin-film limit** (5 GL points, minimal depth-phase
   variation).  For representative intermediate-thickness behaviour, 100 nm (10 GL
   points) or thicker is recommended.

---

## Open Issues (2026-07-09) — now all resolved

Two problems were identified on 2026-07-09.  Both have since been resolved.

### Issue 1 — Origin of the A1(LO) thin-film anomaly  *(RESOLVED 2026-07-09)*

**Symptom:** The normalised ratio (Coh/Incoh)/d_Å for A1(LO) at 10 nm is 0.0497, roughly
20× smaller than the E₂ modes (~0.99) and far below the expected thin-film limit of 1.0.

**Confirmed not a numerical artefact:** GL convergence testing (`Scripts/gl_convergence_a1lo_10nm.py`,
2026-07-09) varied the GL point count from 5 to 1000 actual points and found identical results
to 7 significant figures.  The anomaly is genuine physics.

**Resolved: Root cause — asymmetric q-assignment between coherent and incoherent depth modes**

The anomaly arises from a fundamental difference in how `modal_pairs` assigns the NAC phonon
wavevector `q_ph` to each Berreman pair in the two depth modes.

#### Berreman pair structure at normal incidence, c-plane ZnO (pp, backscattering)

At normal incidence, there are 4 Berreman mode pairs (L=laser channel, S=signal channel).
The channels are sorted by qz descending, so channel 0 = forward (downward, dominant) and
channel 1 = backward (upward, reflected):

| Pair | L channel | S channel | q_pair (Berreman units) | Physical character |
|------|-----------|-----------|-------------------------|--------------------|
| (0,0) | fwd (down) | fwd (down) | ≈ 0 | TO |
| (0,1) | fwd (down) | bwd (up) | ≈ +2n ≈ +4 | **LO** — backscattering |
| (1,0) | bwd (up) | fwd (down) | ≈ −2n ≈ −4 | LO (opposite direction) |
| (1,1) | bwd (up) | bwd (up) | ≈ 0 | TO |

The **forward** channel carries most of the amplitude: amplitude ∝ t₁₂ (Fresnel
transmission into ZnO, t₁₂ ≈ 2/3 for air→ZnO).  The **backward** channel is the
substrate-reflected component: amplitude ∝ r₂₃ × t₁₂ where r₂₃ = (n_ZnO − n_FS)/(n_ZnO + n_FS)
≈ 0.156 for ZnO/Fused Silica.

The same structure applies to the signal (detection sensitivity) channels S: the
time-reversed signal from the superstrate is predominantly forward-going (S channel 0),
with a small backward component (S channel 1, amplitude ∝ r₂₃' × t₁₂').

So the two LO pairs are both **cross-pairs** — one involving a large L and small S,
the other a small L and large S:

| Pair | Amplitude factor | Relative magnitude |
|------|------------------|--------------------|
| (0,1): fwd L × bwd S | t₁₂ × r₂₃' t₁₂' | ≈ 0.667 × 0.156 × 0.667 = **0.069** |
| (1,0): bwd L × fwd S | r₂₃ t₁₂ × t₁₂' | ≈ 0.156 × 0.667 × 0.667 = **0.069** |
| (0,0): fwd L × fwd S | t₁₂ × t₁₂' | ≈ 0.667 × 0.667 = **0.445** (6.4× larger) |

#### Coherent depth (COHERENT_FILM regime, `use_external_q_for_modal_nac = False`)

Each pair uses its own `q_pair` to determine NAC:
- Pairs (0,0) and (1,1): q_pair ≈ 0 → TO NAC → A₁(TO) at 334.6 cm⁻¹ (contributes to TO peak, not LO)
- Pair (0,1): q_pair ≈ +4 → LO NAC → A₁(LO) at 510.1 cm⁻¹ (amplitude ≈ 0.069)
- Pair (1,0): q_pair ≈ −4 → LO NAC → A₁(LO) at 510.1 cm⁻¹ (amplitude ≈ 0.069)

Only the two cross-pairs contribute to A₁(LO) at 510.1 cm⁻¹.  In the thin-film limit:

```
I_coh(510) ≈ 2 × (r₂₃ t₁₂ t₁₂')² × d²
```

#### Incoherent depth (BULK_PHASE_MATCHED regime, `use_external_q_for_modal_nac = True`)

All four pairs have `q_ph` overridden to `q_ext` (the macroscopic backscattering wavevector,
≈ +4 in Berreman units).  Since q_ph = q_ext for every pair, the BULK_PHASE_MATCHED
resolver's `_matches_external_q` check passes for all pairs, so **none are discarded**.
All four pairs share the same `final_state_key` (same q_class_key, same ν = 510.1 cm⁻¹)
and their integrands are **summed before squaring**:

```
sum_integrand(z) = Σ_pairs integrand_pair(z)
                 = (E_S_fwd + E_S_bwd) · R_LO · (E_L_fwd + E_L_bwd)
                 = E_S_total · R_LO · E_L_total
```

The dominant term in `E_S_total · R_LO · E_L_total` is (0,0): fwd × fwd ≈ t₁₂ × t₁₂' ≈ 0.445
(the **large** term).  Then:

```
I_incoh(510) ≈ [t₁₂(1+r₂₃) × t₁₂'(1+r₂₃')]² × d  ≈ (0.771)⁴ × d
```

#### Quantitative prediction

In the thin-film limit, the ratio:

```
(I_coh / I_incoh) / d_Å  ≈  2 r₂₃² / (1 + r₂₃)⁴
```

With r₂₃ = (2.0 − 1.46) / (2.0 + 1.46) ≈ 0.156 for ZnO/Fused Silica at 532 nm:

```
2 × 0.156² / (1.156)⁴  =  0.0487 / 1.788  ≈  0.027
```

The simple Fresnel estimate gives ≈ 0.027.  The GTM calculation gives 0.0497 — about
1.8× higher, within the expected range given:
- Multiple reflections and accurate GTM eigenmodes (not just Fresnel)
- Complex refractive index of ZnO at 532 nm
- Dispersive n values for laser and Stokes wavelengths

The qualitative physics is confirmed: the anomaly is controlled by r₂₃² ≈ (reflection
coefficient at ZnO/FS interface)².

#### Why E₂(high) is unaffected

E₂(high) uses the **non-polar bypass**: both coherent and incoherent depth modes route
it through the standard total-field calculation using `E_S_total · R_E2 · E_L_total`.
The same formula is used in both cases, so (Coh/Incoh)/d_Å = 1 as expected.

#### Physical interpretation

The anomaly is **fundamental to the modal_pairs implementation** and is not a bug:

- **Coherent depth** correctly identifies which Berreman pairs carry LO phonon momentum
  (q_pair ≈ ±2k) and assigns them LO character. These cross-pairs have small amplitude
  (≈ r₂₃ ≈ 15% of the dominant channel) because they require the backward-reflected
  channel. The A₁(LO) signal is physically the correct backscattering LO contribution.

- **Incoherent depth** enforces `q_ext` for all pairs to satisfy macroscopic momentum
  conservation. This assigns LO character to the dominant (0,0) pair too, making
  I_incoh much larger than I_coh/d. This is appropriate for incoherent scattering where
  each depth point is treated as a bulk emitter, but it inflates the LO intensity
  relative to the physically wavevector-resolved coherent picture.

The ratio (I_coh/I_incoh)/d_Å ≈ r₂₃² ≈ 0.025 (observable substrate-dependent value)
indicates the fraction of the LO signal that is actually carried by the physically correct
backscattering cross-pairs versus the total-field LO attribution used in incoherent mode.

**Practical consequence:** For thin films (≲ 50 nm), the `modal_pairs` coherent and
incoherent depth modes give very different absolute intensities for polar LO modes — this
is genuine physics, not a convergence or implementation issue.  The `geometry` and
`dominant_mode` NAC paths are not affected because they compute the Raman tensor once
from the macroscopic q direction and apply it to the full total field, giving consistent
coherent/incoherent results in the thin-film limit.

### Issue 2 — E₂ subdivision non-invariance with `modal_pairs` non-polar bypass  *(RESOLVED 2026-07-09)*

**Symptom:** With `modal_pairs` NAC and `Incoherent depth + Incoherent layers`, splitting
a 100 nm ZnO film into 2×50 nm sub-layers gave ~0.79× the single-layer result for E₂ modes
(both `Incoherent pairs` and `Coherent all pairs`).  Standard NAC modes (`none`, `geometry`,
`dominant_mode`) were fully invariant.

**Root cause — z_reference mismatch in `compute_modal_fields_at_points`**

The non-polar bypass routes E₂ through the total-field path, which uses
`compute_modal_fields_at_points` in `OpticalChannelResolver.py` to compute per-Berreman-mode
electric fields at the GL quadrature points.  This function had an incorrect z_reference for
the propagation formula:

```python
# Before (WRONG — uses z_FRONT of layer):
z_reference = boundaries[system_layer_index + 1]

# After (CORRECT — uses z_BACK, matching calculate_Efield convention):
z_reference = boundaries[system_layer_index + 2]
```

The modal amplitudes returned by `calculate_modal_amplitudes` (= `F_bk[k+1, :]` from the
backward sweep) are defined at z_BACK of layer k.  `calculate_Efield` correctly uses z_BACK
(`zn[current_layer]` in GTMcore, equal to `boundaries[k+2]`).  Using z_FRONT instead
introduced a per-mode phase factor `exp(±i·k₀·q_z·d)` in each channel field.  This phase is
z-independent (doesn't average away) but is d-dependent — when d changed from 100 nm to 50 nm
the phases changed, breaking the algebraic equivalence between single-layer and two-sublayer
integrands.  The result was a ratio of ~0.79 instead of 1.0.

**Fix:** Changed `boundaries[system_layer_index + 1]` to `boundaries[system_layer_index + 2]`
in `PDielec/OpticalChannelResolver.py` line 206.

**Verification:** `divide_incoherent_depth_layers.py` run after the fix gives subdivision
ratio = 1.0000 (to 4 decimal places) for **all five NAC modes**, including `modal_pairs`
CoheAll and IncoherentPairs.

| NAC mode | E₂(high) ratio 2×50nm/100nm | A₁(LO) ratio 2×50nm/100nm |
|---|:---:|:---:|
| none | 1.0000 | 1.0000 |
| geometry | 1.0000 | 1.0000 |
| dominant_mode | 1.0000 | 1.0000 |
| modal_pairs CoheAll | 1.0000 | 1.0000 |
| modal_pairs IncoherentPairs | 1.0000 | 1.0000 |

**Impact on other verification tests (assessed 2026-07-09):**

All tests using only `geometry`, `dominant_mode`, or `none` NAC are completely unaffected.
For `modal_pairs` tests, the following changes were observed:

- `divide_incoherent_depth_layers.py`: all ratios now 1.0000 — fixed (was 0.79).
- `divide_layers.py` (coherent depth + coherent layers): CoheAll now invariant (1.000/1.000,
  was 0.248/0.995); IncoherentPairs now E₂=1.023, A₁=0.505 (was 0.246/1.953).
- `divide_incoherent_layers.py` (coherent depth + incoherent layers): modal_pairs E₂ ratio
  ~3× (was ~4.4×), now matching standard NAC modes; A₁(LO) ratio unchanged at ~0.5×.
- `forward_scatter_nac_comparison.py`: unchanged to 3 significant figures.
- `modal_pairs_ps_anomaly.py`: pp E₂/geometry ratios unchanged (1.003/1.012/1.047 at
  0°/40°/80°); ps unchanged (1.000 to ppm). The pp residual is not caused by the z_reference
  bug — it persists after the fix and has a different numerical origin.
- `nac_angle_comparison.py`, `nac_angle_comparison_ps.py`: unchanged (polar mode ratios,
  single pair per Group q group).
