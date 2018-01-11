Title: PDielec: Application Notes
author: John Kendrick
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom.
eMail: j.kendrick@leeds.ac.uk
author: Andrew Burnett
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom
email: a.d.burnett@leeds.ac.uk
bibliography: ./pdielec.bib

[INCLUDE="style"]

[TITLE]

[TOC]

Introduction
========

Several examples are given to illustrate applications of the package.  The calculations used to provide the data for the permittivities are sufficiently accurate to illustrate aspects of the theory. The examples are chosen to show the package being used with the QM packages CASTEP and VASP and with the MM package GULP.

MgO using CASTEP
================

Magnesium oxide is an isotropic medium, the initial unit cell and the space group symmetry ($Fm\overline{3}m$) were taken from the Inorganic Crystal Structure Database (ICSD) [@Hellenbrandt2015] reference number ICSD-52026 [@Tsirelson1998]. The primitive cell was optimized using CASTEP. Norm-conserving pseudo-potentials were used to represent the core electrons of magnesium and oxygen. An energy cutoff of 1000 eV was used with the PBE [@Kresse1999a] density functional and a k-point spacing for the Monkhorst-Pack grid of 0.04 Å^-1^. The primitive cell was optimized and a Density Functional Perturbation Theory (DFPT) calculation of the phonon spectrum at the gamma point was performed. The optimised lattice parameter was found to be 2.1234 Å, compared with the experimental value of 2.107 Å. Only 3 degenerate modes contribute to the permittivity. A summary of the results is presented in Table [#tab-mgo-properties].

~TableFigure {#tab-mgo-properties; caption: "Experimental and Calculated Properties of MgO"}
+------------------------------+--------------+------------|
| Property                     | Experimental | Calculated |
+------------------------------+:------------:+:----------:+
| Unit cell dimensions (Å)     | 2.107        |2.1234      |
| Space group                  | $Fm\bar{3}m$ |            |
| Optical permittivity         |              | 3.140       |
| Static permittivity          |              | 10.00       |
+------------------------------+--------------+------------+
| **TO phonons**                               |  **Frequency**<br>(cm^-1^) | **Intensity**<br>(cm^-1^ D^2^ Å^-2^ amu^-1^) |              | 17.1 (A)   |
| -----------------------------------------+--------------+------------+
| T                                        | 388.3        | 9.29    |
| -----------------------------------------+--------------+------------+
| **LO phonons **                          |              |            |
| -----------------------------------------+--------------+------------+
| T (001)                                 | 593.7        |        |
+------------------------------------------+--------------+------------+
~

Because MgO is isotropic with only a single frequency contributing to the permittivity, it makes a useful example application to illustrate several features of PDielec. The real and imaginary frequency dependent permittivities are shown in Figure [#fig-mgo-permittivity], where a damping factor ($\sigma$) of 10 cm^-1^ has been used. In the Figure the real permittivity at zero frequency corresponds to the static permittivity in Table 4, and at frequencies above the absorption at 388 cm^-1^ the permittivity tends to the optical permittivity as the frequency increases. The real permittivity has zero values at 388.3 and 693.7 cm^-1^ which are the TO and LO frequencies respectively.

~Figure {#fig-mgo-permittivity; Caption: "Permittivity of MgO"}
![#img-MgO-permittivity]
~ 

[#img-MgO-permittivity]: Figures/MgO_Permittivity.svg {width:auto; max-width:90%} 

~Figure {#fig-mgo-real; Caption: "Real and imaginary permittivities of a 10% volume fraction of MgO spheres in PTFE using the Maxwell-Garnett method"}
![#img-MgO-real]
~
[#img-MgO-real]: Figures/MgO_Real_Imaginary.svg {width:auto; max-width:90%}


Using the Maxwell-Garnett mixing rule, Figure [#fig-mgo-real] shows the calculated permittivities of a 10% volume fraction of MgO spheres in a supporting medium with a frequency independent permittivity of 2.0, which would be typical of a material such as PTFE. Due to the dilution effect the real component has shifted to a base line value close to 2, and the absorption, as indicated by the maximum in the imaginary component has shifted by about 150 cm^-1^ to 550 cm^-1^. 

The effect of volume fraction on the predicted molar absorption coefficient, using the Maxwell-Garnett mixing rule, is shown in Figure [#fig-mgo-vf]. The lowest volume fraction of MgO gives the largest shift of the absorption peak to high frequency. As the volume fraction increases the mixing rule predicts a broadening of the absorption, whilst the peak in the molar absorption coefficient moves to lower frequency. At the highest loading (f=0.9) the maximum absorption occurs quite close to the TO frequency. The Maxwell-Garnett mixing rule is regarded as being appropriate for low volume fractions and so should not be used for interpreting results in which higher volume fractions of absorbing media have been used [@Sihvola].

~Figure {#fig-mgo-vf; Caption: "Effect of volume fraction on the Maxwell-Garnett molar absorption coefficient of MgO spheres in PTFE"}
![#img-mgo-vf]
~
[#img-mgo-vf]: ./Figures/MgO_MG_Volume_Fraction.svg {width:auto; max-width:90%}

Figure 6 shows the same plot for the Bruggeman mixing rule. At low
volume fractions the Bruggeman mixing rule predicts a similar absorption
to the Maxwell-Garnett. Indeed as the volume fraction approaches zero
the two rules predict the same absorption characteristics. However, even
at the relatively low 1% loading , the Bruggeman mixing rule shows
additional broadening of the peak, the shape of the absorption peak has
lost its Lorentzian characteristic shape as can be seen clearly in
Figure 5 At 10% loading the Bruggeman predicted absorption is broad with
the peak shifted to lower wavenumber. This broadening increases with
increased loading until, at the higher loadings, the TO peak begins to
dominate the absorption.

~Figure {#fig-mgo-vf-bg; Caption: "Effect of volume fraction on the Bruggeman molar absorption coefficient of MgO spheres in PTFE" }
![#img-mgo-vf-bg]
~
[#img-mgo-vf-bg]: ./Figures/MgO_Bruggeman_Volume_Fraction.svg {width:auto; max-width:90%}

~Figure {#fig-mgo-varying-permittivity; Caption:"The Maxwell-Garnett molar absorption coefficients of spherical MgO particles, 1% volume fraction, embedded in media of varying permittivities" }
![#img-mgo-varying-permittivity]
~
[#img-mgo-varying-permittivity]: Figures/MgO_Varying_Permittivity.svg "MgO: Effect of permittivity according to Maxwell-Garnett"{ width:auto; max-width:90% }

Figure [#fig-mgo-varying-permittivity] shows the effect of varying the permittivity of the supporting medium. The calculations were performed on spherical MgO particles with a 1% volume fraction. The lowest permittivity is that of a vacuum (or air) and shows the highest shift of the absorption maximum to higher frequencies. Increasing the permittivity lowers the shift until it becomes quite small. A similar effect is seen for the Bruggeman mixing model. However, the absorption resulting for particles in a low dielectric medium is considerable broader than that seen in the Maxwell-Garnet case. This broadening reduces as the permittivity of the medium increases (see Figure [#fig-mgo-varying-permittivity-brug]).

~ Figure {#fig-mgo-varying-permittivity-brug; Caption:"The Bruggeman molar absorption coefficients of spherical MgO particles, 1% volume fraction, embedded in media of varying permittivites"}
![Figure07_MgO_varying_permittivity_bruggeman]
~
[Figure07_MgO_varying_permittivity_bruggeman]: Figures/MgO_Varying_Permittivity_Bruggeman.svg "MgO: Effect of permittivity according to Bruggeman" { width:auto; max-width:90% }


ZnO using VASP
==============

Zinc oxide crystallizes in space group $P6_3mc$ (wurtzite). All calculations were performed by VASP [@Hafner2008c] using projector augmented-wave PAW [@Kresse1999a] pseudo-potentials, the PBE [@Perdew1996a] density functional, an energy cutoff of 600 eV and a k-point resolution of approximately 0.1 Å^-1^. The initial unit cell was taken from the ICSD [@]Hellenbrandt2015] with code ICSD-26170 [@McNally2012a]. The unit cell and atom positions were optimized using
VASP and the permittivity was calculated using DFPT and the results reported in Table [#tab-zno-properties]. Only two of the bands showed any significant intensity, a doubly degenerate band (E) with a TO frequency of 372.1 cm^-1^ and a non-degenerate band (A) with a TO frequency of 350.0 cm^-1^. The LO frequency of the non-degenerate band is shifted to 502.0 cm^-1^ for a wave-vector with direction (001), whilst the degenerate modes are unaffected. In the case of the (010) direction the LO frequency of one of the E modes is shifted to 511.2 cm^-1^. It is known that ZnO can crystallize with a plate morphology [@Yamamoto1977] with the (001) surface dominant. Calculations of the molar absorption were performed for a sphere, plate and needle like shapes with the unique directions of the plate and the needle being normal to the (001) surface. A volume fraction of 1% was chosen for these calculations and the predicted molar absorption coefficients for the Maxwell-Garnett mixing rule is shown in Figure 9.

~TableFigure {#tab-zno-properties; caption: "Experimental and Calculated Properties of ZnO"}
+------------------------------------------+--------------+------------+
| Property                                 | Experimental | Calculated |
+------------------------------------------+:------------:+:----------:+
| Unit cell dimensions a,b (Å)             | 3.250        |3.295       |
| Unit cell dimensions a   (Å)             | 5.207        |5.285       |
| Space group                              | $P6_3mc$     |            |
| Optical permittivity   xx,yy               |              | 5.09       |
| Optical permittivity   zz                 |              | 6.00       |
| Static permittivity    xx, yy              |              | 10.83      |
| Static permittivity    zz                 |              | 11.67      |
| -----------------------------------------+--------------+------------+
| **TO phonons**                               |  **Frequency (cm^-1^)**| **Intensity (cm^-1^ D^2^ Å^-2^ amu^-1^)** |              | 17.1 (A)   |
| -----------------------------------------+--------------+------------+
| A                                        | 350.0        | 17.1       |
| E                                        | 372.1        | 16.4       |
| -----------------------------------------+--------------+------------+
| **LO phonons **                          |              |            |
| -----------------------------------------+--------------+------------+
| (001) A                                  | 502.0        |        |
| (010) E                                  | 511.2        |        |
+------------------------------------------+--------------+------------+
~

~ Figure {#fig-zno; caption:"The effect of shape on the Maxwell-Garnett molar absorption coefficient of 1% volume fraction ZnO in PTFE"}
![Figure08_ZnO]
~
[Figure08_ZnO]: Figures/ZnO-Maxwell_Garnett.svg "ZnO: Effect of shape on absorption" { width:auto; max-width:90% }

For the Maxwell-Garnett mixing rule the sphere morphology results in the two absorption peaks shifting from their TO positions to higher wavenumber by about 80 cm^-1^. The plate morphology results in one of the peaks moving to higher wavenumber by about 130 cm^-1^, whilst the other remains at the TO position. The Maxwell-Garnett results are in close accord with some experimental results by Yamamoto et al [@Yamamoto1977] who measured the infrared spectrum of ZnO smoke particles and observed peaks in the absorption at 380, 530 and 550 cm^-1^. Previous work [@Rendon1981;@Hayashi1977a] have also used effective medium theory to explain the observed spectrum.

Calcite using GULP
==================

Calcite is the most stable polymorph of calcium carbonate and the crystal structure belongs to the $R\overline{3}c$ space group. The force field and atomic structures used here are described in detail in work by Fisler [@Fisler2000]. Briefly, the oxygen ions are described using a core-shell model [@Dick1958]. The carbon - oxygen potential of the carbonate is taken to be a Morse potential and an additional 3 atom potential is used to maintain the O-C-O angle at 120^o^. The van der Waals interactions between non bonded atoms are taken to be Buckingham potentials and the charges on the calcium, carbon and oxygen ions are +2, +1.3435 and -1.1145 respectively. The shell charge of the oxygen ion is -2.133 and the spring constant for the core-shell interaction is 52.74 eV/Å^2^. The unit cell was optimized using the primitive unit cell and the full space group symmetry. The calculation of the phonon spectrum was performed without symmetry but still using the primitive cell of the lattice. A summary of the calculated properties is given in Table [#tab-calcite-properties].


~TableFigure {#tab-calcite-properties; Caption: "Experimental and Calculated Properties of Calcite" }

+------------------------------------------+--------------+------------+
| Property                                 | Experimental | Calculated |
+------------------------------------------+:------------:+:----------:+
| Unit cell dimensions a,b,c (Å)           | 6.375        |6.376       |
| Unit cell angles                         | 46.1         |46.0        |
| Space group                              | $R\overline{3}c$|            |
| Optical permittivity   xx,yy               |              | 1.91       |
| Optical permittivity   zz                 |              | 2.00       |
| Static permittivity    xx, yy              |              | 6.70       |
| Static permittivity    zz                 |              | 7.10       |
| -----------------------------------------+--------------+------------+
| **TO phonon symmetry**                               |  **Frequency (cm^-1^)**| **Intensity (cm^-1^ D^2^ Å^-2^ amu^-1^)** |              | 17.1 (A)   |
| -----------------------------------------+--------------+------------+
| E                                        | 114.8        | 2.39       |
| A                                        | 127.4        | 3.36       |
| A                                        | 249.3        | 1.23       |
| E                                        | 320.7        | 5.82       |
| A                                        | 338.1        | 4.14       |
| E                                        | 620.1        | 3.38       |
| A                                        | 732.0        | 26.89      |
| E                                        |1463.6        | 16.97      |
+------------------------------------------+--------------+------------+
~

Figure [#fig-calcite] shows the results of analysis of the results using PDielec. The damping parameter used in the calculation was a value of 5 cm^-1^. A 10% volume fraction was used with sphere and plate morphologies for the particles. The unique axis for the plate was taken to be the normal to the (211) surfaces in the primitive cell axes (or the {104} surfaces in the standard unit cell). Such surfaces define the rhombohedral faces commonly seen in calcite crystals.^50^ Figure 10 shows that the doubly degenerate TO absorption peak at 620 cm^-1^ is not significantly affected by spherical particles and there is a small shift to higher frequencies in the case of plate-like particles. The non-degenerate TO transition at 732 cm^-1^, which corresponds to motion of the carbon atom of the carbonate along the unique direction of the slab, shows a shift to 786 and 819 cm^-1^ for the sphere and plate respectively. The doubly degenerate peak at 1463 cm^-1^ is shifted to 1480 cm^-1^ by spherical particles and is split by plate-like particles with one component which shifts to 1491 cm^-1^ .

~ Figure {#fig-calcite; Caption:"Calculated Maxwell-Garnett absorption spectrum of 10% volume fraction of calcite in PTFE"}
![Figure09_calcite]
~
[Figure09_calcite]: Figures/Calcite_AP_Plate211_Sphere.svg "Calculated Maxwell-Garnett absorption" { width:auto; max-width:90% }

Fluoroapatite using VASP
========================

The line shapes of the infrared absorption of apatite and fluoroapatite were examined extensively by Balan *et al* [@Balan2008b]. Their calculations included the effect of crystallite habit on the spectrum and the results reported here are similar to their conclusions. The method used by Balan *et al*. is an infinitely dilute Maxwell-Garnett model, so the only difference between the methods used by them and those reported here using PDielec are the incorporation of the volume fraction into the theory and the use of an ellipsoidal shape for comparison with the other shapes.

All calculations were performed by VASP [@Hafner2008c] using projector augmented-wave PAW [@Kresse1999a] pseudo-potentials, the PBE [@Perdew1996a] density functional, an energy cutoff of 600 eV and a k-point resolution of approximately 0.1 Å^-1^. Table [#tab-fluoroapatite-properties] summarises the results of the calculations. Only the 3 highest frequency bands are reported and discussed. The TO intensity of the highest frequency band at 1038 cm^-1^ is low and will not be discussed further. The Bravais Friedel Donnay Harker (BFDH) [@Donnay1937] crystal habit of the optimized crystal is shown in Figure 11. The habit was calculated using the Mercury software package. [@Macrae2008]. The BFDH crystal habit is often used to give an idea of the likely important faces of a crystal. It uses only the crystal lattice and space group to determine the crystal morphology. Figure [#fig-fluoroapatite-morphology] shows that the {100} surfaces form a tube which are capped by the {011} surfaces. The effect of different particle shapes on the predicted spectrum is shown in Figure [#fig-fluoroapatite-absorption]. The calculations of the spectra were performed with a damping parameter (σ) of 2 cm^-1^. The ellipsoid was chosen to have an aspect ratio, a/b, of 2 and a principle axis along [001], which was compatible with the morphology predicted by the BDFH method. The two TO absorption frequencies at 981 and 986 cm^-1^ have A and E symmetry respectively. Spherical crystallites result in three absorption peaks at around 1000, 1010 and 1015 cm^-1^. Needle shaped crystallites leave the A symmetry TO absorption peak at 981 cm^-1^ unaffected, but shift and split the E symmetry TO peak to 1020 and 1046 cm^-1^. A plate morphology with (100) surfaces results in the A and one component of the E TO absorption peak remaining at the TO frequencies, with the other component of the E shifting 85 cm^-1^ to 1075 cm^-1^. The ellipsoidal morphology show three shifted peaks at 1000, 1018 and 1045 cm^-1^. These results are consistent with those of Balan *et al*. [@Balan2008b], who gave detailed results for hydroxyapatite.


~TableFigure {#tab-fluoroapatite-properties; Caption: "Experimental and Calculated Properties of Fluoroapatite" }

+------------------------------------------+--------------+------------+
| Property                                 | Experimental [@Hughes1989] | Calculated |
+------------------------------------------+:------------:+:----------:+
| Unit cell dimensions a,b (Å)             | 9.417       |9.447       |
| Unit cell dimensions c (Å)               | 6.875       |6.926     |
| Space group                              | $P6_3m$      |       |
| Optical permittivity   xx,yy               |              | 2.891      |
| Optical permittivity   zz                 |              | 2.894      |
| Static permittivity    xx, yy              |              | 12.081     |
| Static permittivity    zz                 |              | 8.841      |
| -----------------------------------------+--------------+------------+
| **TO phonon symmetry**  |  **Frequency (cm^-1^)**| **Intensity (cm^-1^ D^2^ Å^-2^ amu^-1^)** |
| -----------------------------------------+--------------+------------+
| A                                        | 981.8        | 112.6      |
| E                                        | 986.3        | 101.0      |
| E                                        | 1038.1       | 7.92       |
+------------------------------------------+--------------+------------+
~

~Figure {#fig-fluoroapatite-morphology; caption:"BFDH Morphology of fluoroapatite"}
![Figure10_Fluoroapatite_morphology]
~
[Figure10_Fluoroapatite_morphology]: Figures/Fluoroapatite_morphology.png "Figure10_Fluroapatite_morphology" { width:auto; max-width:90% }

~Figure {#fig-fluoroapatite-absorption; caption:"Calculated Maxwell-Garnett absorption spectra of 10% fluoroapatite in PTFE"}
![Figure11_Fluoroapatite_absorption]
~
[Figure11_Fluoroapatite_absorption]: Figures/Fluoroapatite_absorption.svg "Figure11_Fluoroapatite_absorption" { width:auto; max-width:90% }


L-aspartic Acid using CASTEP
============================

L-aspartic acid is a zwitterion in the solid state and so the shape of the particles used in the measurement of IR and THz spectra maybe important. The starting geometry for optimization of the unit cell and molecular structure of L-aspartic acid was taken from Derissen et al [@Derissen1968]. The PBE [@Perdew1996a] functional was used with a plane wave energy cutoff of 1000 eV and norm conserving pseudo-potentials. A dispersion correction using the Tkatchenko-Scheffler scheme [@Tkatchenko2009c] available in CASTEP was applied for both the geometry optimisation and the calculation of the phonon spectrum at the gamma point, with a value S~6~ scaling factor [@Juliano2015] of 1.0. A summary of the results of the calculations is shown in Table [#tab-aspartic-properties].


~TableFigure {#tab-aspartic-properties; Caption: "Experimental and calculated properties of L-aspartic acid" }

+------------------------------------------+--------------+------------+
| Property                                 | Experimental [@Derissen1968] | Calculated |
+------------------------------------------+:------------:+:----------:+
| Unit cell dimensions a (Å)             | 7.617      |7.597       |
| Unit cell dimensions b (Å)               | 6.982     |7.028    |
| Unit cell dimensions c (Å)               | 5.142      |5.113     |
| Unit cell angle $\beta$               | 99.84       |98.77   |
| Space group                              | $P2_1$      |       |
| Optical permittivity   xx               |              | 2.68      |
| Optical permittivity   yy               |              | 2.20      |
| Optical permittivity   zz                 |              | 2.56      |
| Static permittivity    xx              |              | 4.58     |
| Static permittivity    yy              |              | 3.65     |
| Static permittivity    zz                 |              | 3.65      |
| -----------------------------------------+--------------+------------+
| **TO phonons **  |  **Frequency (cm^-1^)**| **Intensity (cm^-1^ D^2^ Å^-2^ amu^-1^)** |
| -----------------------------------------+--------------+------------+
|                                      | 84.5     | 0.120    |
|                                        | 104.7        | 0.202     |
|                                        | 106.0       | 0.243     |
|                                        | 115.3       | 0.474     |
|                                        | 137.3        | 0.617    |
|                                        | 1290.0       | 55.0     |
|                                        | 2945.9        | 102.8     |
|                                        | 2947.3       | 48.2     |
|                                        | 3053.7       | 44.1    |
+------------------------------------------+--------------+------------+
~

The THz spectrum of L-aspartic acid has been reported by Juliano and Korter [@Juliano2015] in the frequency range 0-90 cm^-1^. The infrared spectrum has been reported and assigned by Lopez *et al* [@LopezNavarrete1994]. Figure [#fig-aspartic] shows the calculated absorption spectra for L-aspartic acid for three frequency ranges. The calculation of the spectra used the Maxwell-Garnett mixing rule with a 10% volume fraction of L-aspartic acid in PTFE and for comparison the TO mixing rule. A damping factor of 2 cm^-1^ was used. Spherical and a variety of plate-like inclusions were used to illustrate their effect on the absorption spectra. Figure [#fig-aspartic]a shows the frequency range from 60-130 cm^-1^ which is that covered by THz spectroscopy. The shifts observed for the different particle morphologies are not large, but the change in intensities is significant. The molecular motions associated with phonons at these frequencies tend to be whole molecule motion involving rotation. Figure [#fig-aspartic]b shows the frequency range from 1260-1340 cm^-1^. In this frequency range bending of the carboxylate anion contributes to the spectrum significantly. The three different plate morphologies show different and significant shifts in the TO absorption peak at 1290 cm^-1^. The spherical morphology shows a shift of around 25 cm^-1^ to higher wavenumber. Figure [#fig-aspartic]c shows the spectra in the frequency range 2900-3100 cm^-1^, which corresponds to the motion of O-H (below 2980 cm^-1^) and N-H (above 2980 cm^-1^) stretching. The effect of the different possible crystal morphologies is large with shifts to higher frequency of up to 50 cm^-1^. The spectra below 3000 cm^-1^ arises from two TO absorptions at 2946 and 2947 cm^-1^. Because the motions associated with each mode interact differently with the internal field within each crystal they give rise to different shifts producing more complex spectra.

~ Begin Figure {#fig-aspartic; caption: "Calculated Maxwell-Garnett absorption spectra of 10% volume fraction of L-aspartic acid in PTFE "}
~ Begin SubFigureRow {vertical-align: bottom}
~ SubFigure {#fig-aspartic-a; caption: "Frequency range 60-130 cm^-1^"}
![Aspartica]
~
~ SubFigure {#fig-aspartic-b; caption: "Frequency range 1260-1340 cm^-1^"}
![Asparticb]
~
~ SubFigure {#fig-aspartic-a; caption: "Frequency range 2900-3100 cm^-1^"}
![Asparticc]
~
~ End SubFigureRow
~ End Figure
 
[Aspartica]: Figures/Aspartic_thz_a.svg "Maxwell-Garnett absorption of L-aspartic acid" { width:auto; max-width:90%; max-height:90% }
[Asparticb]: Figures/Aspartic_thz_b.svg "Maxwell-Garnett absorption of L-aspartic acid" { width:auto; max-width:90%; max-height:90% }
[Asparticc]: Figures/Aspartic_thz_c.svg "Maxwell-Garnett absorption of L-aspartic acid" { width:auto; max-width:90%; max-height:90% }


MgO Example using Mie Scattering
================================

Figure 13 compares a Mie scattering calculation with the results from Maxwell-Garnett effective medium theory. The same data set was used for the CASTEP, MgO example. A volume fraction of 1% was used with a small sphere radius (0.1 μm) and a broadening of 5 cm^-1^ embedded in a matrix of PTFE. A power expansion in the size parameter of the Mie expressions
indicates that for small sizes of particles, the Mie and the Maxwell-Garnett methods should be the same. This is verified in Figure 13.

~ Figure {#fig-mgo-mie; caption: "Comparison of Mie and Maxwell methods. 1% volume fraction of MgO in PTFE, sphere radius of 0.1 μm and a broadening of 5 cm^-1^"}
![mgo-maxwell-mie]
~
[mgo-maxwell-mie]: Figures/MgO_Mie_Maxwell.svg "MgO-maxwell-mie" { width:90%; max-width:90% }

To better understand what makes particles large or small Table [#tab-mgo-mie-sizes] shows the dimensionless size parameter, $x$, as a function of wavenumber and of sphere radius. It has been assumed that the supporting medium is PTFE. Since the power expansion of the size parameters leads to terms which are quadratic in x, it should be expected that when x is less than about 0.1 μm, the particles can be considered small. It can be seen that particles less than 0.01 μm are small over the range of frequencies considered. But while 0.1 μm particles are small in the THz regime and low frequency infrared they should not be considered small over the more extended infrared frequencies. 1μm particles should be considered large for both the THz and the extended infrared.

~ TableFigure {#tab-mgo-mie-sizes; Caption: "Variation of size parameter with wavenumber and radius of sphere"}
|------------|----------------------|--------|--------|--------|+--------|
| Wavenumber | Radius of sphere (μm) |||||
|            |                      |        |        |        |         |
| (cm-1)     | 0.001                | 0.01   | 0.1    | 1      | 10      |
|:----------:|:--------------------:|:------:|:------:|:------:|:-------:|
| 100        | 0.0009               | 0.0089 | 0.0888 | 0.8884 | 8.8841  |
| 200        | 0.0018               | 0.0178 | 0.1777 | 1.7768 | 17.7682 |
| 300        | 0.0027               | 0.0267 | 0.2665 | 2.6652 | 26.6523 |
| 400        | 0.0036               | 0.0355 | 0.3554 | 3.5536 | 35.5364 |
| 500        | 0.0044               | 0.0444 | 0.4442 | 4.4420 | 44.4204 |
| 600        | 0.0053               | 0.0533 | 0.5330 | 5.3305 | 53.3045 |
| 700        | 0.0062               | 0.0622 | 0.6219 | 6.2189 | 62.1886 |
| 800        | 0.0071               | 0.0711 | 0.7107 | 7.1073 | 71.0727 |
| 900        | 0.0080               | 0.0800 | 0.7996 | 7.9957 | 79.9568 |
| 1000       | 0.0089               | 0.0888 | 0.8884 | 8.8841 | 88.8409 |
|------------|----------------------|--------|--------|--------|---------|
~

Figure [#fig-mgo-mie-sizes] shows how the Mie predictions change as the particle radius changes from 0.2 to 1.6 μm. As the particle size increases the peak above 500 cm^-1^ splits into two. One broader peak which moves to lower frequency as the particle size increases and the peak at about 550 cm^-1^ which looses intensity as the particle size increases. There is also the onset of absorption at 388 cm^-1^ which corresponds to the bulk TO modes.

~ Figure {#fig-mgo-mie-sizes; caption: "Variation in absorption calculated by the Mie method for different radii of spheres. 1% volume fraction of MgO in PTFE and a broadening of 5 cm^-1^"}
![mgo-mie-sizes]
~
[mgo-mie-sizes]: Figures/MgO_Mie_Sizes.svg "Mie calculations on MgO" { width:90%; max-width:90% }

Figure [#fig-mgo-mie-sizes-big] shows the effect of increasing the particle size further. More structure appears in the absorption, with increasing absorption around
the bulk TO frequency. Above 4.0 μm there is more low frequency structure
appearing, below 300 cm^-1^.

~ Figure { #fig-mgo-mie-sizes-big; caption: "Variation in absorption calculated by the Mie method for larger sphere radii. 1% volume fraction of MgO in PTFE and a broadening of 5 cm^-1^"}
![mgo-mie-sizes-big]
~
[mgo-mie-sizes-big]: Figures/MgO_Mie_Sizes_Big.svg "Mie calculations on MgO" { width:90%; max-width:90% }

ZnO Example using Mie Scattering
================================

ZnO is an anisotropic material, so the treatment described here using Mie scattering is an approximation. However, the permittivity constant tensor is diagonal due to the space group symmetry of the crystal. 

~Figure {#fig-zno-maxwell-garnett; caption: "ZnO spheres in PTFE using Maxwell Garnett"}
![zno-maxwell-garnett]
~
[zno-maxwell-garnett]: Figures/ZnO-Maxwell_Garnett.svg "ZnO spheres using Maxwell Garnett" { width:90%; max-width:90% }

~Figure {#fig-zno-mie; caption: "ZnO spheres in PTFE using Mie"}
![zno-mie]
~
[zno-mie]: Figures/ZnO_Mie.svg "ZnO spheres using Mie" { width:90%; max-width:90% }

Figures [#fig-zno-maxwell-garnett] and [#fig-zno-mie] compare the capabilities of the Maxwell-Garnett and Mie methods for describing volume fraction effects. Figure [#fig-zno-maxwell-garnett] shows that Maxwell-Garnett predicts a lowering of intensity and frequency of the high frequency peak as the volume fraction is increased. Figure [#fig-zno-mie] illustrates that Mie theory shows no effect of the change in volume fraction. This is to be expected as the theory assumes that each sphere is isolated and not affecting the other spheres around it. It should be pointed out that the figures are plotting molar absorption coefficients. The actual absorption would increase with volume fraction of ZnO.

Figure [#fig-zno-mie-sizes] shows that the variation of the Mie scattering with sphere radius follows a similar pattern to that observed in MgO, though slightly more complex. The initial peaks at about 440 and 460 cm^-1^ broaden and shift to lower frequencies as the particle size increases. Bulk bands around 350 and 372 cm^-1^ can be seen which grown in intensity and shift to lower frequencies as the particle size increases.

~Figure {#fig-zno-mie-sizes; caption: "Mie scattering of 10% volume fraction ZnO spheres in PTFE using a line broadening factor of 5 cm^-1^. "}
![zno-mie-sizes]
~
[zno-mie-sizes]: Figures/ZnO_Mie_Sizes.svg "zno-mie-sizes" { width:90%; max-width:90% }



[BIB]
