author: John Kendrick
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom.
eMail: j.kendrick@leeds.ac.uk
author: Andrew Burnett
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom
title: PDielec: Theory
email: a.d.burnett@leeds.ac.uk
Bibliography: ./pdielec.bib

[INCLUDE="style"]

<!-- Comment out some of the options -->
<!-- Csl Style: ieee -->
<!-- Math Mode: static -->
<!-- [INCLUDE="style"] -->
<!-- Colorizer: javascript -->
<!-- Doc class: [10pt]article -->

~ MathDefs
\newcommand{\water}{H_{2}O}
\newcommand{\tensor}[1]{\bar{\bar{#1}}}
\newcommand{\tensorbs}[1]{\bar{\bar{\bm{#1}}}}
\newcommand{\tensorbf}[1]{\bar{\bar{\bm{#1}}}}
\newcommand{\fieldbf}[1]{\bar{\bm{#1}}}
~

[TITLE]

[TOC]

# INTRODUCTION

The theory underlying the Python package PDielec is described.  PDielec calculates the infrared absorption characteristics of a crystalline material supported in a non-absorbing medium.  PDielec post processes solid state quantum mechanical and molecular mechanical calculations of the phonons and dielectric response of the crystalline material. The molecular and solid state quantum mechanical (QM) calculations of response properties such as the frequencies and intensities of infrared (IR) and terahertz (THz) radiation absorption has become generally available in many molecular and solid state computer programs. A common approach is to assume the harmonic approximation and calculate the mass weighted force constant matrix (for molecules) or the dynamical matrix at the gamma point (for periodic solids). Diagonalisation of the matrix gives the frequencies for absorption and the normal modes (molecules) or phonon displacements (periodic solids). The calculation of the absorption intensity for each mode requires the calculation of the change in dipole moment caused by the displacement of the atoms for that mode. For solids where there is a large separation of charge, there can be a large coupling between a phonon mode and the internal field within a particle resulting from its morphology. 

PDielec post processes the output of solid state quantum mechanical and molecular mechanics (MM) based codes such as VASP [@Hafner2008c], CASTEP [@Clark2005d], CRYSTAL [@Dovesi2014], Abinit [@Gonze2016] , Quantum Espresso [@Giannozzi2009], Phonopy [@Togo2015], and GULP [@Gale2003] to predict the infrared absorption of crystalline insulator materials whose crystal size is small compared with the wavelength of the absorbing radiation. The package is suited for the calculation of the complex, frequency dependent permittivity and its associated absorption of infrared radiation for a finely ground crystalline material dispersed in a low loss dielectric medium such KBr or Polytetrafluoroethylene (PTFE). A particular feature of the program is its ability to take into account the constant permittivity of the supporting medium and the particle shape of the material of interest through an effective medium theory. The paper outlines the theory used by the program and gives some examples of the application of the program for ionic and molecular materials.

# THEORY

Equations [#eq-beer1] and [#eq-beer2] describe Beer-Lambert's law [@Bertie2006] where $\alpha$ is the (decadic) absorption coefficient (usually given in cm^-1^), $I$ and $I_0$ are the intensities after and before absorption respectively and $d$ is the path length.

~ Equation {#eq-beer1}
\frac{I}{I_{0}} = 10^{- \alpha d}
~

~ Equation {#eq-beer2}
log\left( \frac{I}{I_{0}} \right) = -\alpha d 
~

It is common, especially in the chemistry community, when reporting infrared spectra to use a decadic molar absorption coefficient ($a$), which has units of L mol^-1^cm^-1^. The relationship between the absorption coefficient and the molar absorption coefficient [@Bertie2006] is;

~ Equation {#eq-aC}
\alpha = aC
~


where $C$ is the concentration of the absorbing species.

## Molecular Approach to Absorption Intensity

For molecules the transition intensity $I_k$ of the $k^{th}$ vibrational mode (calculated from the change in dipole moment along the mode displacement) can be converted to an integrated molar absorption coefficient, $A_k$, which can then be more readily compared with experiment. The theory for this is described by Wilson, Decius and Cross [@Wilson1955] and results in expressions such as those given in equations [#eq-Absorption1] and [#eq-Absorption2].  The first expression shows the relationship between the integrated molar absorption coefficient and the transition intensity and uses the number of molecules per unit volume ($N$), the velocity of light ($c$) and the degeneracy of the mode ($g_k$). The second expression shows the appropriate conversion factors if the units for the integrated molar absorption coefficient are L mol^‑1^cm^‑2^ (1 L mol^-1^cm^-2^ = 0.01 km mol^-1^) and the units for the transition intensity are D^2^ Å^-2^ amu^-1^, where D represents the Debye unit of dipole  moment and amu is an atomic mass unit. The factor log~e~10 arises due to the choice of a decadic Beer's law.

~ Equation {#eq-Absorption1}
A_k =  \frac{N\pi}{3c^{2}\log_e10}g_kI_k
~

~ Equation {#eq-Absorption2}
A_k =  \frac{Na\pi}{3000c^{2}2.302585}g_kI_k  = 4225.6I_k
~

The derivation of the above expressions assumes that the rotational levels are not quantised and that the vibrational levels are thermally occupied according to a Boltzmann distribution. In order to use the calculated molecular intensities to predict a spectrum it is usual to assume [@Wilson1955] that each transition is associated with a Lorentzian line shape with a full width at half maximum (FWHM) of $\sigma_k$. It is common, when reporting comparison between theoretical and experimental spectra, to assume that the line widths are the same for all modes [@Juliano2013; @Burnett2013].  Recent work on terahertz absorption in crystalline pentaerythritol tetranitrate (PETN) using molecular dynamics calculations [@Pereverzev2011b] in combination with the direct calculation of the cubic anharmonic couplings of the normal modes [@Pereverzev2011b] has shown that the FWHM of the intense absorptions may vary between 10 and 25 cm^-1^. Assuming a Lorentzian line shape, the molar absorption coefficient for the $k^{th}$ mode at wavenumber, ${\bar{\nu}}_{k}$, can be written as a function of frequency or wavenumber ($\bar{\nu}$);

~ Equation {#eq-lorentzian}
a_k(\bar{\nu}) = \frac{2A_k}{\pi}\frac{\sigma_k}{4\left( \bar{\nu} - {\bar{\nu}}_k \right)^{2} + \sigma_k^2} 
~

~ Equation {#eq-akmax}
a_k^{max} = \frac{2A_k}{\pi\sigma_k}
~

The maximum height of the Lorentzian, $a_{k}^{\max}$clearly depends upon the value of $\sigma_k$. As can be seen in Equation [#eq-integratedmolarintensity], the choice of normalisation for the Lorentzian means that integration of the molar absorption coefficient over wavenumber returns the integrated molar absorption coefficient and a sum over all the bands provides the total molar absorption coefficient $a^{mol}(\bar{\nu})$ as a function of wavenumber, calculated from the intensities of each band.  Equation [#eq-molarabsorption] shows the relationship between the absorption and the molar absorption coefficients. $C$ is the concentration usually expressed in mol L^-1^.
~ Equation {#eq-integratedmolarintensity}
A_k = \int{a_k(\bar{\nu})d\bar{\nu}}
~

~ Equation {#eq-amol}
a^{mol}(\bar{\nu}) = \sum_k{a_k(\bar{\nu})}
~

~ Equation {#eq-molarabsorption}
\alpha^{\text{mol}}(\bar{\nu}) = Ca^{\text{mol}}(\bar{\nu})
~

A comment should be made about the various units which can be used for these quantities. A common unit for the transition intensity is D^2^ Å^-2^ amu^-1^, another is km mol^-1^. However, it should be pointed out that strictly speaking the latter unit refers to the integrated molar absorption coefficient as defined above in Equation [#eq-integratedmolarintensity] and therefore relies on the assumptions made in its derivation. ( 1 D^2^ Å^-2^ amu^-1^ is equivalent to 42.256 km mol^-1^ ).

## Solid State Approach to Absorption Intensity

The optical properties of a solid are determined by its complex, frequency dependent relative permittivity,
$\boldsymbol{\varepsilon}(\bar{\nu})$, and in particular the real and imaginary components, $\boldsymbol{\kappa}$ and $\mathbf{n}$, of the complex refractive index, $\mathbf{N}{(\bar{\nu})}$ where;

~ Equation {#eq-refractiveindex1}
\tensorbf{N}({\bar{\nu}})^{2} = {\tensorbs{\varepsilon}}({\bar{\nu}})
~

~ Equation {#eq-refractiveindex2}
{\tensorbf{N}}\left( \bar{\nu} \right) = {\tensorbf{n}}\left( \bar{\nu} \right) + i{\tensorbs{\kappa}}\left( \bar{\nu} \right) \
~

The intensity of absorption is given by the imaginary component of the refractive index which for an isotropic material is [@VanDeHulst1981];

~ Equation {#eq-beerIntensity}
I = I_{0}e^{- \frac{4\pi\kappa(\bar{\nu})d}{\lambda}} \\                                 
I = I_{0}e^{- 4\pi\bar{\nu}\kappa(\bar{\nu})d} 
~

~ Equation {#eq-logIntensity}
- ln\left( \frac{I}{I_{0}} \right) = 4\pi\bar{\nu}\kappa(\bar{\nu})d \\              
- log\left( \frac{I}{I_{0}} \right) = 4\pi\bar{\nu}\kappa(\bar{\nu})d \cdot log(e)
~

Comparison with the definition of the absorption coefficient from Beer-Lambert's law, Equation [#eq-beer1], and using Equation [#eq-molarabsorption] gives;

~ Equation {#eq-alphasol}
\alpha^{sol}(\bar{\nu}) = 4\pi\bar{\nu}\kappa(\bar{\nu}) log(e) \
~

~ Equation {#eq-asol}
a^{sol}(\bar{\nu}) = \frac{\alpha^{\text{sol}}(\bar{\nu})}{C}
~

Since the refractive index is dimensionless, the absorption coefficient, $\alpha^\text{sol}$ is specified in cm^-1^. The superscripts 'sol,' for solid, and 'mol,' for molecular, are used here to distinguish between the two methods of calculating the absorption ($\alpha$) and molar absorption coefficients ($a$).  In the calculation of the imaginary component of the refractive index it is necessary to choose the solution which gives a positive value. This is consistent with the Kramers-Kronig relationship between the real and imaginary components [@Wooten1972].

In order to calculate the relationship between absorption and molar absorption coefficients it is necessary to know the concentration. For solid state calculations the required unit is; moles of unit cells per litre. One of the drawbacks of this molar absorption coefficient unit is that the number of molecules in a unit cell can change depending on whether a supercell, primitive or non primitive unit cell is being used. A more natural unit would be to use a mole of formula units, or a mole of molecules. In order to aid comparison between such calculations PDielec is able to calculate concentration in both moles of atoms and moles of molecules. However for the rest of this paper Equation [#eq-concentration] will be used, where $V$ is the volume of the unit cell, and therefore the concentration $C$ is moles of unit cell/litre.

~ Equation {#eq-concentration}
C = \frac{f \cdot 1000cm^{3}}{VN_{a}} 
~

The volume fraction, $f$, of the dielectric material in a supporting matrix of non-absorbing material is included in the expression for the concentration as it will be useful when the theory for mixtures is developed.

For a periodic system the permittivity tensor can be calculated as a sum over Lorentz oscillators, incorporating an imaginary loss component through the damping factor *σ~k~* [@Gonze1997].  The frequencies of the oscillators are the transverse optic (TO) phonon frequencies of the system.

~ Equation {#eq-permittivity}
{\tensorbs{\varepsilon}}(v) = {{\tensorbs{\varepsilon}}}_{\infty} + \frac{4\pi}{V}\sum_{k}\frac{{{\tensorbf{S}}}_{k}}{\nu_{k}^{2} - \nu^{2} - i\sigma_{k}\nu}
~

~ Equation {#eq-oscillatorstrength}
\tensorbf{S}_{k} = {\bar{\mathbf{Z}}}_{k}{\bar{\mathbf{Z}}}_{k}^\text{T}
~

~ Equation {#eq-borncharges}
\bar{\mathbf{Z}}_{k} =  \sum_{a}{\tensorbf{Z^{a}}}{\bar{\mathbf{U^{a}}}}_{k}
~

~ Equation {#eq-eigenvalues}
\tensorbf{D} {\bar{\mathbf{U}}}_{k} = \Lambda_{k}{\bar{\mathbf{U}}}_{k}\\
\nu_{k}^{2} = \Lambda_{k}
~

~ Equation {#eq-intensity}
I_{k} = tr\left( \tensorbf{S}_{k} \right) 
~

$V$ is the volume of the unit cell, ${\tensorbf{S}}_{k}$ is the dipole oscillator strength tensor for the $k_{th}$  transition, with a TO frequency of $\nu_k$ and $\tensor{\varepsilon}_{\infty}$ the optical permittivity tensor, which represents the electronic contribution to the permittivity. The intensity of a transition, $I_k$, is given by the trace of the oscillator strength tensor, Equation [#eq-intensity].  The damping factor $\sigma_{k}$ removes any discontinuities at the TO frequencies. Since the oscillator strengths and phonon frequencies can be calculated routinely in solid state quantum mechanical packages, the calculation of the frequency dependent complex permittivity using Equation [#eq-permittivity] is straightforward. 

In some cases, using Equations [#eq-oscillatorstrength] and [#eq-borncharges], PDielec calculates the oscillator strengths from the Born charge matrix for atom $a$ and its contributions to the $k_{th}$ phonon mode [@Gonze1997].  As shown in Equation [#eq-eigenvalues], at the $\Gamma$ point the $k_{th}$ phonon mode is described by the eigenvector, $\bar{\mathbf{U}}_{k},$ and eigenvalue, $\Lambda_{k}$, of the mass weighted, dynamical matrix, $\tensorbf{D}$, which is a 3Nx3N matrix, where N is the number of atoms in the unit cell. The eigenvalues are the squared frequencies of the phonon modes, Equation [#eq-eigenvalues].  The displacement of each atom in the $k_{th}$ is proportional to $m_{a}^{-1/2}$, where $m_a$ is the mass of atom $a$. The dynamical matrix has 3N eigenvectors and eigenvalues, of which three should be zero due to translational invariance. If there are any negative eigenvalues the system is unstable to some displacement and therefore not at an energy minimum.

For ionic systems it is common practice in solid state QM and MM programs to include a long wave-length, non-analytic correction to the mass weighted dynamical matrix at the $\Gamma$ point, which describes the coupling of the longitudinal optic (LO) modes to the induced field resulting from the vibration. This may be written for atoms $s$ and $t$ and their Cartesian components $\alpha$ and $\beta$ as [@Gonze1997];

~ Equation {#eq-LODynamicalMatrix}
\left( \tensorbf{D}^\text{LO}_{\mathbf{q \rightarrow 0}} \right)_{s,\alpha;t,\beta} = \left( \tensorbf{D} \right)_{s,\alpha;t,\beta}\mathbf{+}\frac{4\pi}{V\sqrt{M_{s}M_{t}}}\mathbf{\ }\frac{\left( {\bar{\mathbf{q}}}^{\text{T}}{{\mathbf{\ }\tensorbf{Z}}}_{s} \right)_{\alpha}\left( {\bar{\mathbf{q}}}^{\text{T}}{\tensorbf{Z}}_t \right)_\beta}{{\bar{\mathbf{q}}}^{\text{T}}\mathbf{\cdot}{{\tensorbs{\varepsilon}}}_{\infty} \cdot \bar{\mathbf{q}}}
~

The mass weighting has been incorporated through the mass of the atoms, $M_s$ and $M_t$.  The correction depends upon the direction, $\bar{\mathbf{q}}$, that the long wave-length limit is approached.  Diagonalisation of the corrected matrix gives the squared frequencies of N-1 LO modes and 2N-2 TO modes, Equation [#eq-eigenvalues].  In some of the examples given below the LO frequencies will be given for comparison with the TO frequencies.

## Effect of Particle Shape on Infrared Absorption

It has long been recognised that, especially for ionic materials, the local field within a crystal and its coupling with the transverse optical phonons has an important effect on the position and intensity of the absorption. Fröhlich [@Frohlich1948] was one of the first to point out that the frequency of absorption of a small ionic sphere embedded in a low dielectric medium is shifted to lie between the transverse and longitudinal optical frequencies of the material making up the sphere.

In the development of the theory used in PDielec an important assumption is that the particle size of the crystallites in the sample is small compared with the wavelength of light.  Using this approach Genzel and Martin [@Genzel1972a] were able to explain the observed infrared absorption of small spheres of MgO crystallites and the effect of the permittivity of the supporting medium on the spectrum.  Studies of the infrared absorption by small particles of α-Fe~2~O~3~ using an effective medium theory and an absorption/scattering theory [@Serna1987; @Iglesias1990] showed that in order to fit the experimental spectra it was necessary to adjust not only the damping factors in Equation [#eq-permittivity] but also the permittivity of the matrix and the volume fraction of the dielectric medium. The latter was used to account for aggregation effects as the volume fraction increased. It was also shown that effective medium theories were only applicable for particles smaller than the wavelength of light. For larger particles the
scattering from the particles becomes increasingly important.

More recently Balan and others in a series of papers [@Balan2010b;@Balan2008b;@Fourdrin2009] used density functional calculations together with an effective medium theory to calculate the infrared absorption of several minerals incorporating information about the crystallite shape.  In an experimental and theoretical study of irradiated kaolinite  [@Fourdrin2009] it was shown that exposure to radiation resulted in shifts in the infrared spectrum which could be accounted for by increasing the polarisability of the particles through
an increase in the optical permittivity tensor.

The underlying theory adopted by PDielec is based on similar premises to the work described above, namely that the dielectric response of small spherical, ellipsoidal, slab-like or needle-like crystallites randomly distributed in a non-absorbing medium such as PTFE, KBr or Nujol, is the same as that of an effective medium material whose frequency dependent dielectric response can be calculated from the frequency dependent permittivity tensor of the crystal (as calculated by solid state QM or MM calculations), the shape of the crystallites and the permittivity of the non-absorbing medium (taken to be a constant over the frequency range of interest).

The development of the theory reported here closely follows the work of Sihvola [@Sihvola].  It will be assumed that the inclusion particles, which may be non-isotropic, ellipsoidal (including spherical, needle-like and plate-like), are randomly orientated in an embedding, non-absorbing medium such as PTFE, KBr or Nujol. It should be emphasized that whilst PDielec can take account of particle shape, particle and matrix permittivity there are many additional aspects of infrared absorption which need to be considered when comparing calculated and experimental results. Most notable of these are; the coupling between phonons and mobile electrons or holes (so called phonon-polariton coupling)[@Ruggiero2015], the scattering which starts to dominate as the particles get larger [@Fourdrin2009] and the agglomeration of particles as the volume fraction increases.

### The polarisability of an isolated particle

Figure [#fig-polarisation] shows a schematic of the field and polarisation inside an inclusion with non-isotropic permittivity $\tensor{\varepsilon}_{i}$ embedded in a supporting medium with permittivity, $\varepsilon_e$.  The internal field within the inclusion is indicated by $\fieldbf{E}_{i}$ the external, applied field is indicated by $\fieldbf{E}_{e}$ and the induced polarisation in the inclusion is shown by $\fieldbf{P}$.

~ Figure {#fig-polarisation; caption: "Schematic showing the field and polarisation inside an inclusion with non-isotropic permittivity ${\tensorbf{\varepsilon}}_{i}$ embedded in a supporting medium with permittivity $\varepsilon_e$. The internal field within the inclusion is indicated by $\fieldbf{E}_i$, the external, applied field is indicated by $\fieldbf{E}_e$ and the induced polarisation in the inclusion is shown by $\fieldbf{P}$"}
![img-polarisation]
~
[img-polarisation]: ./Figures/Polarisation_Schematic.png "Polarisation and field schematic" {width:auto; max-width:90%}

The electric field internal to the inclusion gives rise to a polarisation density which is no longer necessarily aligned with the field because the material is non-isotropic. The polarisation density in the inclusion can be expressed as the tensor product of the permittivity contrast between the inclusion and the supporting medium and the (as yet unknown) internal field.

~ Equation {#eq-PolarisationDensity}
\fieldbf{P} = \left( \tensorbs{\varepsilon}_{i} - \varepsilon_{e}{\tensorbf{1}} \right){\fieldbf{E}}_{i}
~

For any ellipsoidal shape (including sphere, slab and needle) with volume *V*, the polarisation density throughout the particle is uniform and integrating over all space gives the field induced dipole moment of the inclusion, $\fieldbf{p}$.

~ Equation {#eq-polar1}
\fieldbf{p} = V\fieldbf{P} = V\left( \tensorbs{\varepsilon}_{i} - \varepsilon_{e}\tensorbf{1} \right)\fieldbf{E}_{i}
~

The dipole and the external field, ${\fieldbf{E}}_{e}$, are related by the polarisability tensor, $\tensorbs{\alpha}$.

~ Equation {#eq-polar2}
\fieldbf{p} = \tensorbf{\alpha}\fieldbf{E}_{e}
~

Equations [#eq-polar1] and [#eq-polar2] allow the determination of the polarisability, once the field internal to the inclusion has been expressed in terms of the shape of the inclusion and its permittivity. The polarisation within the inclusion gives rise to a depolarisation field, $\fieldbf{E}_{d}$, which depends on the shape of the inclusion through the symmetric and unit-trace depolarisation tensor, $\tensorbf{L}$.

~ Equation {#eq-DepolarisationField}
\fieldbf{E}_d = - \frac{1}{\varepsilon_e}\tensorbf{L} \cdot \fieldbf{P}
~

The internal field is the sum of the external field and the depolarisation field.

~ Equation {#eq-InternalField}
\fieldbf{E}_{i} = \fieldbf{E}_{e} + \fieldbf{E}_{d}
~

The depolarisation tensor acts as a projection or screening operator describing the effect of the geometry of the inclusion on the depolarisation field which results from its polarisation. For instance, in the case of a needle, only polarisation perpendicular to the needle axis contributes to the depolarizing field, whilst for a slab only polarization perpendicular to the slab face may contribute. Similarly for a sphere, all directions contribute and so the depolarisation matrix is diagonal with a value of 1/3 for each diagonal element, as the trace of the depolarisation tensor must be 1. It follows from Equations [#eq-PolarisationDensity], [#eq-DepolarisationField] and [#eq-InternalField] that;

~ Equation {#eq-InternalField2}
\fieldbf{E}_i = \fieldbf{E}_e - \frac{1}{\varepsilon_e} \tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_e{\tensorbf{1}} \right){\fieldbf{E}}_i
~

Rearrangement allows the internal field of the inclusion to be expressed in terms of the known permittivities, the shape of the inclusion and the external field.

~ Equation {#eq-InternalField3}
\fieldbf{E}_i \left( \varepsilon_e{\tensorbf{1}} + \tensorbf{L}  \left( \tensorbs{\varepsilon}_i - \varepsilon_e{\tensorbf{1}} \right) \right) = \varepsilon_e\fieldbf{E}_e
~

~ Equation {#eq-InternalField4}
\fieldbf{E}_i = \varepsilon_e{\fieldbf{E}}_e\left( \varepsilon_e\tensorbf{1} + \tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_e \tensorbf{1} \right) \right)^{- 1}
~

Substituting the internal field expression, [#eq-InternalField2], into Equation [#eq-polar1] for the dipole moment and requiring the dipole moments calculated using the polarisation density to equal those calculated from the polarisability allows the polarisability to be written;

~ Equation {#eq-polarisation}
\tensorbs{\alpha} = V\varepsilon_e\left( \tensorbs{\varepsilon}_i - \varepsilon_e\tensor{1} \right)\left( \varepsilon_e\tensorbf{1} + \tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_e\tensorbf{1} \right) \right)^{- 1}
~

Although it has not been specified explicitly the permittivity of the inclusion, and therefore the polarisability tensor, are frequency dependent through the oscillator strengths of each phonon mode contributing to the permittivity according to Equation [#eq-permittivity].  The calculation of the complex, frequency dependent polarisability tensor of the composite material is the key step in the calculation of its effective permittivity.

### The Effective permittivity of a mixture

To extend this approach to include the effect of a number of inclusions
we need to introduce the concept of an effective permittivity, $\tensorbf{\varepsilon}_{eff}$, which describes the behaviour of an average field, $\left\langle \fieldbf{E} \right\rangle$, where the angle brackets indicate an average over a volume of the composite material. It is required that the average electric flux density $\left\langle \fieldbf{D} \right\rangle$ is the same in the effective medium as in the composite medium;

~ Equation {#eq-averageFluxDensity}
\left\langle \fieldbf{D} \right\rangle = \tensorbs{\varepsilon}_{eff}\left\langle \fieldbf{E} \right\rangle = \varepsilon_e\left\langle \fieldbf{E} \right\rangle + \left\langle \fieldbf{P} \right\rangle
~

The averaging is necessary because the polarisation within a given inclusion has an effect on the field in other inclusions. The local field in the cavity left by a single inclusion embedded in the average polarisation field is given by;

~ Equation {#eq-localField}
\fieldbf{E}_L = \left\langle \fieldbf{E} \right\rangle + \frac{1}{\varepsilon_e}\tensorbf{L}\left\langle \fieldbf{P} \right\rangle
~

The local field 'excites' the inclusion resulting in a dipole moment $\mathbf{p}_{mix}$ that is related to the polarisation through the number density of inclusions, $n$, and through the polarisability of the inclusion, which is already known from Equation [#eq-polarisation].

~ Equation {#eq-PolarisationField}
\left\langle \fieldbf{P} \right\rangle = n{\fieldbf{p}}_{mix} = n\left\langle \tensorbs{\alpha}\fieldbf{E}_{L} \right\rangle
~

The angle brackets around the product of the polarisability and the local field indicate that it is necessary to average the polarisation according to the distribution of alignments of inclusions. In this work it will be assumed that the inclusions are randomly aligned.  Substituting the expression for the local field, Equation [#eq-localField], gives;

~ Equation {#eq-averagePolarisation}
\left\langle \fieldbf{P} \right\rangle = \left( \tensorbf{1}
- \frac{n\left\langle \tensorbf{\alpha}\tensorbf{L} \right\rangle} {\varepsilon_e} \right)^{- 1}
n\left\langle \tensorbs{\alpha} \right\rangle\left\langle \tensorbf{E} \right\rangle 
~

### Mixing rules

There are many mixing rules which have been proposed to describe the homogenization of composite materials and a lot of work has been done in comparing their accuracy. Here two methods will be used. The first and the most commonly used method is the Maxwell-Garnett mixing rule [@Sihvola].  Indeed this has been implied by the use of Equation [#eq-averageFluxDensity] to define the effective permittivity. The other commonly used method is the Bruggeman mixing rule [@Sihvola], which differs substantially in the way the two components of the system are treated. It is usually stated that the Maxwell-Garnet mixing rule is good for low volume fractions of the inclusion and the Bruggeman approach should be better for higher volume fractions [@Giordano2003a].  In addition to these mixing rules one other approach will be described, namely the Averaged Permittivity (AP) mixing rule, which calculates the absorption spectrum ignoring the effects of the internal field on the absorption and can therefore be used as an indicator of the shifts in frequency and intensity which have occurred owing to the effect of particle shape.

### Maxwell-Garnett mixing rule

The Maxwell-Garnett approach for treating the properties of heterogeneous mixtures assumes that the average field and the average flux density result from volume fraction weighted sums. Substituting Equation [#eq-averagePolarisation] into Equation [#eq-averageFluxDensity] gives the Maxwell-Garnett effective permittivity;

~ Equation {#eq-mg}
\tensorbs{\varepsilon}_{mg} = \tensorbf{1} + \left( \tensorbf{1} - \frac{n\left\langle {\tensorbs{\alpha}}{\tensorbf{L}} \right\rangle}{\varepsilon_e} \right)^{- 1}n\left\langle \tensorbs{\alpha} \right\rangle 
~

The fact that the polarisability tensor has a volume term in it (Equation [#eq-polarisation]) means that the terms in Equation [#eq-mg] containing $n\tensorbs{\alpha}$ depend on the volume fraction $f$. Although written as a tensor, because the assumption has been made that the inclusions are randomly orientated, the effective permittivity has to be diagonal with equal complex values. Since the polarisability is complex and frequency dependent the effective permittivity is also and its calculation using Equations [#eq-mg] and [#eq-polarisation] need to be calculated over the frequency range of interest.

### Bruggeman mixing rule

In the Maxwell-Garnett mixing formalism there is a distinction between the inclusion and the supporting medium which results in an asymmetry in the treatment of the two species in the mixture. Instead the Bruggeman mixing rule assumes that each species is polarized against the background of the effective medium and therefore the polarisation in the two components cancel each other out;

~ Equation {#eq-EqualPolarisation}
\left\langle \fieldbf{P}_1 \right\rangle + \left\langle \fieldbf{P}_2 \right\rangle = 0\
~

where the components are now labeled 1 and 2 rather than external and internal. The polarisation for species 1 and 2 with a number density of species represented by $n_1$  and $n_2$ can be obtained from the polarisability of the species (Equation [#eq-EqualPolarisation]);

~ Equation {#eq-polarisation1}
\left\langle \fieldbf{P}_1 \right\rangle = n_1 \left\langle \tensorbs{\alpha}_1 \right\rangle\fieldbf{E}
~

Substituting Equation [#eq-polarisation1] into Equation [#eq-EqualPolarisation] leads to the requirement that;

~ Equation {#eq-BruggemanPolarisability1}
n_1\left\langle \tensorbs{\alpha}_1 \right\rangle + n_2\left\langle \tensorbs{\alpha}_2 \right\rangle = 0
~

Taking Equation 19 and generalizing it for species $i$, (where $i$ is 1 or 2) embedded in an effective permittivity given by $\tensorbs{\varepsilon}_{br}$;

~ Equation {#eq-BruggemanPolarisability2}
\tensorbs{\alpha}_i = V_i \tensorbs{\varepsilon}_{br} \left( \tensorbs{\varepsilon}_i - \tensorbs{\varepsilon}_{br} \right)
\left( \tensorbs{\varepsilon}_{br} + \tensorbf{L} \left( \tensorbs{\varepsilon}_i - \tensorbs{\varepsilon}_{br} \right) \right)^{- 1}
~

Equation [#eq-bruggemanpolarisability1] has to be solved for$\tensorbs{\varepsilon}_{br}$. Since the systems considered here are isotropic with random inclusions, a solution has to be found for a complex value of the Bruggeman permittivity at each frequency considered.  An issue in the use of Equation [#eq-bruggemanpolarisability2] is that the same depolarisation matrix is being used for both species, which is clearly not always appropriate. The solution of Equation [#eq-BruggemanPolarisability1] can be achieved either by iteration or by casting the equation as a minimization problem. The iterative approach implemented in PDielec involves repeated application of Equation 29 until convergence [@Mackay2009].  The starting point for the iterations is taken as the Maxwell-Garnett solution for the first frequency and then the solution at the previous
frequency is used to start the iterations.

~ Equation {#eq-br}
\tensorbs{\varepsilon}_{br} = \frac{f_1{ {\tensorbs{\varepsilon}}}_1\left\lbrack {\tensorbf{1}} + {\tensorbf{L}}\left( {\tensorbs{\varepsilon}}_1 - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1} +
f_2\tensorbs{\varepsilon}_2\left\lbrack \tensorbf{1} + \tensorbf{L}\left( \tensorbs{\varepsilon}_{2} - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1}}
{f_{1}\left\lbrack \tensorbf{1} +{\tensorbf{L}}\left( \tensorbf{\varepsilon}_1 - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1} +
f_{2}\left\lbrack \tensorbf{1} + \tensorbf{L}\left( \tensorbs{\varepsilon}_2 - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1}} 
~
Although the Bruggeman permittivity is written here as a tensor, the polarisabilities in Equation [#eq-BruggemanPolarisability1] have to be averaged over the random orientation of the inclusions and therefore the homogenized material is isotropic with a single complex value for the diagonal tensor.  Also, as with the Maxwell-Garnett mixing rule, since the polarisability is complex and frequency dependent, the effective permittivity is also, and its calculation using Equation [#eq-mg] needs to be performed over the frequency range of interest.

The choice between using the Bruggeman or Maxwell-Garnett model is often governed by the assumption that the Maxwell-Garnett model works well at low concentrations and the Bruggeman model works better at higher concentrations. Work by Karkkainen *et al*. using a finite difference method for random mixtures of non-absorbing materials indicated that the Bruggeman approximation works best when there is some clustering of the inclusions and the Maxwell Garnett model works best when there is no clustering [@Karkkainen2001].

The Bruggeman solution has been shown to be unphysical in certain circumstances [@Jamaian2010].  In particular when the real components of the permittivity have different signs and when the absolute values of the real components are much larger than those of the imaginary components. Unfortunately, it may be that these conditions will apply to modelling infrared absorption. As a result only a few of the examples discussed below will include results using the Bruggeman mixing rule; the majority will use the Maxwell-Garnett mixing rule.

### Averaged-Permittivity mixing rule

It is useful to be able to compare the effective medium theories with the absorption predicted using no shape information, that is using only the TO frequencies.

~ Equation {#eq-ap}
\tensorbs{\varepsilon}_{TO} = f\left\langle \tensorbs{\varepsilon}_i \right\rangle + \left( 1 - f \right)\varepsilon_e
~

Equation [#eq-ap] defines an isotropic permittivity which can be used to calculate such an absorption coefficient. The angle brackets indicate an average of orientation. This mixing rule provides a useful comparison between the absorption calculated without any shape effects and that calculated including shape effects using the other mixing rules presented above. At low concentrations the peak positions of the AP mixing rule will be at the TO frequencies.

## Particle Size Effects

Meier and Wokaun [@Meier1983] outlined an approach to treating large (metal) spherical particles, where particle size is incorporating terms up to 3^rd^ order in the wave vector *k*. Using Equations [#eq-PolarisationDensity] and [#eq-InternalField] we can write;

~ Equation {#eq-PolarisationWithSize}
\fieldbf{P} = \left( \tensorbs{\varepsilon}_i - \varepsilon_e \tensorbf{1} \right)\left( \fieldbf{E}_e + \fieldbf{E}_d \right)
~

~ Equation {#eq-FieldWithSize}
  \begin{aligned}
\fieldbf{E}_d &= - \frac{\left( 1 - x^2 - i\frac{2}{3}x^3 \right)}{\varepsilon_e}{\tensorbf{L}} \fieldbf{P} \\
x &= ak = \frac{2\pi a}{\lambda} 
\end{aligned}
~

Here $a$ is the radius of thespherical particle and $x$ is the dimensionless 'size' of the particle with respect to the wavelength of the incident light. The first term relating the depolarisation field to the polarisability is the same as that used above in Equation [#eq-DepolarisationField]. The second term is a dynamic depolarisation term and the third term is a radiation damping correction to the electrostatic solution[@Meier1983].

A slightly different, but related, approach is presented by [@Sihvola].  Starting from Equation [#eq-InternalField2], a size dependent term $G(x)$ is introduced as indicated by the work of Peltoniemi [@Peltoniemi1996];

~ Equation {#eq-SizeTerms}
  \begin{aligned}
G\left( x \right) &= G_1\left( x \right) + G_2\left( x \right) \\
G_{1}\left( x \right) &= \ \frac{2}{3}\left\lbrack \left( 1 + ix\right)e^{- ix} - 1 \right\rbrack \\
G_{2}\left( x \right) &= \ \left ( 1 + ix - \frac{7}{15}x^{2} - i\frac{2}{15}x^{3} \right)e^{- ix}-1
 \end{aligned}
~

The modified equation for the relationship between the internal field and the external field based on equation [#eq-InternalField2] becomes; 

~ Equation {#eq-ModifiedInternalField}
\fieldbf{E}_i = \fieldbf{E}_e - \frac{\left( 1 - G\left( x \right) \right)}{\varepsilon_e}\tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_e\tensorbf{ 1} \right)\fieldbf{E}_i
~

This leads to a modified equation for the polarisability of spherical particles;

~ Equation {#eq-ModifiedPolarisability}
\tensorbs{\alpha}\left( x \right) = V\varepsilon_e\left( \tensorbs{\varepsilon}_i - 
\varepsilon_e\tensorbf{1} \right)\left( \varepsilon_e\mathbf{1} + 
\left( 1 - G\left( x\right) \right)\tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_{e}\mathbf{1} \right) \right)^{- 1}
~

Using the modified, sized dependent polarisability all the Bruggeman and Maxwell mixing schemes can be implemented in a way that incorporates size effects.  Generally speaking the on-set of changes in the calculated absorption is an indication that size effects are important and should be treated properly.

## Light Scattering by Spherical Particles using Mie Theory

For particles which are comparable in size to the wavelength of light, the theory developed by Mie and described fully by van de Hulst [@VanDeHulst1981] can be used. Unfortunately this theory is only applicable to spherical, isotropic particles where the separation between the particles is large compared with the wavelength of light.

PDielec implements a form of Mie theory using the Python library PyMieScatt [@Sumlin2018a]. In order to treat systems which are anisotropic, PDielec first of all transforms the real component of the permittivity so that the tensor is diagonal. Then the full  permittivity (real and imaginary) is transformed to this basis.

~ Equation {#eq-transform1}
  \begin{aligned}
& {\tensorbf{U}^{\text{T}}} \tensorbs{\varepsilon}^{real} \tensorbf{U} = \tensorbs{\varepsilon}^{diagonal} \\
& {\tensorbs{\varepsilon}}^{new} = \tensorbf{U}^\text{T} \tensorbs{\varepsilon}^{full}\tensorbf{U}
\end{aligned} 
~

The core of Mie theory is the calculation of efficiency factors S1(θ) and S2(θ) (page 125 in reference [@VanDeHulst1981].  These factors depend upon the size of the particle with relative to the wavelength of the light ($x$) and the refractive index of infrared active material ($m_1$) relative to that of the supporting medium ($m_2$). Both refractive indices are scalars, $m_1$ is complex as light is absorbed by the material and $m_2$ is assumed to be real.

For each of the three diagonal elements of $\tensorbs{\varepsilon}_{new}$ the relative refractive index, $m$ is calculated;

~ Equation {#eq-MieVariables}
  \begin{aligned}
m_1 &= \tensorbs{\varepsilon}_{w,w}^{new} \\
m &= \frac{m_{1}}{m_{2}} \\
x &= \frac{2\pi a}{\lambda} \\
\lambda &= \frac{\lambda_{\text{ vac}}}{m_{2}}
\end{aligned}
~

The size parameter depends upon the refractive index of the supporting medium because the wavelength of light is smaller than its vacuum value in a medium with a refractive index greater than one. The subscript $w$ indicates one of the three directions x, y or z in the basis which diagonalises the real permittivity tensor.

Using the routines available in PyMieScatt [@Sumlin2018a] the efficiency factor $S(0)$, where $S(0) = S1(0) = S2(0)$, is calculated from $m$ and $x$. The extinction (this includes absorption and scattering) along with retardation are reflected in the overall complex refractive index, page 129 in reference [@VanDeHulst1981];

~ Equation {#eq-MieVariable2}
  \begin{aligned}
\widetilde{m}_w &= m_2\left( 1 - iS\left( 0 \right)2\pi Nk^{- 3} \right) \\
N &= \frac{f}{V_{sphere}}  \\
V_{sphere} &= \frac{4}{3} \pi a^{3} \\
k &= \frac{2\pi}{\lambda}
\end{aligned}
~

Here $N$ is the number density of particles and $k$ is the wave-vector of the light in the surrounding medium. The subscript $w$ indicates that there are three values for the effective refractive index, depending on which diagonal value of the permittivity tensor is taken. Once the effective refractive index is known in each direction the effective permittivity can be calculated using the one third rule [@Stout2007];

~ Equation {#eq-effectiveri}
  \begin{aligned}
\widetilde{m}_{eff} &=  \frac{1}{3}\sum_w{\widetilde{m}}_w \\
\varepsilon_{eff} &= \widetilde{m}_{eff}^2
\end{aligned}
~

This approach to taking anisotropy into account when the embedded particles are anisotropic but randomly oriented is approximate, but has been shown to have a reasonably wide range of application [@Stout2007].

# IMPLEMENTATION

The above theory has been implemented in a Python 2/3 package which is available for download [@pdielec]. The package requires SCIPY [@scipy], NUMPY [@scipy}], PyYAML (py-yaml) [@pyyaml], PyMieScatt [@Sumlin2018a], xlswriter [@xlsxwriter] and if visualization of the predicted spectra is required MATPLOTLIB [@scipy]. The program is run from the command line. There are several command options and these are summarized below in Table 1. At the moment the package has interfaces to five solid state QM codes, VASP [@Hafner2008c], CASTEP [@Clark2005d], CRYSTAL14 [@Dovesi2014], Abinit [@Gonze2016], Quantum Espresso [@Giannozzi2009] and Phonopy [@Togo2015].  In addition an interface is available for GULP [@Gale2003] which is a force field based solid state code. Finally an interface has been written to and 'experiment' file format which allows the preparation of a user defined file specifying the permittivities and absorption frequencies. The origin of the dataset(s) used for processing is determined by a command line switch, -program. An outline of the interfaces to these codes is given here. The package used for the calculation is described by the --program option. In addition a file name is given which contains the output which will be processed by PDielec.

## VASP (-program vasp OUTCAR)
The name provided on the command line is an OUTCAR file. The OUTCAR is read by PDielec to determine the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. The VASP run can be a DFPT or numerical calculation of the response.

## CASTEP (-program castep seedname)
The name provided on the command line is the seedname for the calculation. The corresponding seedname.castep file in the current directory is read and processed to determine the unit cell, atomic masses, optical permittivity and born charge tensors. The normal modes and their frequencies are determined from the seedname.phonon file. The CASTEP run needs to be a DFPT (phonon+efield) task.

## CRYSTAL (-program crystal outputfilename)
The name on the command line is a file ending in .out, containing the output of a CRYSTAL14 run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory where PDielec is run from , it uses these files to calculate the Born charge tensors, frequencies and normal modes. The CRYSTAL calculation needs to be a frequency calculation (FREQCALC) with the infrared intensity (INTENS) selected. The default algorithm does not calculate the optical permittivity, so this needs to be provided on the command line. However, if the CPHF or CPKS algorithm is used for the frequency calculation, the optical permittivity is calculated and PDielec will automatically read it from the output file. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the PDielec package. Small differences in the calculated frequencies between the CRYSTAL program and PDielec have been observed. These have been found to be due to a slightly different method for symmetrising the 2^nd^ derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that PDielec should use the same symmetrisation as CRYSTAL14.

## ABINIT (-program abinit outputfilename)
The output file should come from a run containing three datasets. One to calculate the wavefunction at the optimized geometry, one to calculate the field perturbations and one to calculate the second derivatives. Examples of input files and output files are available with the distribution.

## QE (-program qe outputfilename)
The output file is the dynamical matrix file, specified by "filedyn" in a run of the quantum espresso phonon package. Examples of input and output files are given in the PDielec distribution

## PHONOPY (-program phonopy vasp OUTCAR.born)
Phonopy calculates the dynamical matrix through numerical differentiation. It has interfaces to several programs, although PDielec has only used the VASP interface. In principle other interfaces could be used. The second parameter for the --program directive is the PHONOPY interface that was used to calculate the forces. Typically these would be generated by performing;

        phonopy --d --dim="1 1 1"

to calculate the displacements in a set of POSCAR-* files. After running VASP a single point VASP calculation for each displacement. The FORCE\_SETS file can then be calculated using for example;

        phonopy --f DISP-*/vasprun.xml

where the DISP-\* directories are where the VASP calculation was performed. Finally a dynamical is written out using;

        phonopy --dim="1 1 1" --qpoints="0 0 0" --writedm

To calculate the infrared spectrum PDielec needs the Born charges for the atoms in the unit cell and these can be calculated using VASP and the optimized geometry of the unit cell. The OUTCAR file from this calculation can be copied to the current directory and renamed OUTCAR.born

## GULP (-program gulp outputfilename)
The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit cell. The key words; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it on the command line (see -optical and -optical\_tensor options below).

## EXPERIMENT (-program experiment\_filename)
This option has been added to allow the exploration of 'toy' problems. The file contains a minimum amount of data to allow a calculation to proceed. It is assumed that the systems will be isotropic as that make the input simpler. Calculations of the LO frequencies will not work with this option. An example input file is given here, which gives results very similar to that found for the Castep, MgO example;

```
species 2 # Define the atomic species, followed by the number of species
O 16.0    # Species name followed by its mass
Mg 24.3
lattice 2.12346 # Define the lattice, the parameter is a scaling factor
1 1 0           # The a vector
1 0 1           # The b vector
0 1 1           # The c vector
unitcell 2      # Define the unit cell contents and fractional coordinates
O 0.5 0.5 0.5   # Species fraca, fracb, fracc
Mg 0.0 0.0 0.0  # Species fraca, fracb, fracc
static          # Define the static permittivity
3.13969 0.0 0.0
0.0 3.13969 0.0
0.0 0.0 3.13969
frequencies 3       # Specify how many frequencies to read in
388.282 0.000073639 # The frequency in wavenumbers and oscillator strength
388.282 0.000073639
388.282 0.000073639
```

Examples of data sets for these packages are included with the distribution. The interface to these QM and MM codes reads information about the unit cell, the calculated normal modes and the Born charge matrices; from these the permittivity is calculated over the frequency range requested. The absorption and molar absorption coefficients can be plotted along with the real and imaginary permittivities. Optionally all the information can be written to a comma separated values (csv) file for direct importing into a spreadsheet. The program is run from the command line. There are several command options and these are summarized below. Some of the options may be repeated. The package needs a shape to be specified (sphere, needle, plate or ellipse). If no shape is specified on the command line a sphere is assumed.

## Program Options
The program options are summarised below, a ✔ indicates that the option can be specified more than once.  Where there is a default value, it value is shown.

* -program program
  : Program can be “abinit”,  “castep”, “crystal”, “gulp”, “qe”, “experiment” or “vasp” and specifies the program which generated the results to be analysed.
* -method maxwell (✔)
  : The method is given by the string and is either ‘ap’, ‘maxwell’, ‘bruggeman’ or ‘mie’.
* -sphere (✔)
  : The inclusion is a sphere, the default if no other shape is given.
* -needle h k l (✔)
  : The inclusion is a needle whose unique directionis given by the direction \[hkl\].
* -plate h k l (✔)
  : The inclusion is a plate whose  surface is defined by the Miller indices \(hkl\).   Note that needles and ellipsoid use  directions in crystal coordinates defined by \[hkl\].   For non-orthogonal lattices the normal to  the \(hkl\) is not necessarily the same as \[hkl\].
* -ellipse h k l z (✔)
  : The inclusion is an ellipsoid, whose  unique direction is given by \[hkl\],  z  specifies the eccentricity of the ellipsoid.
* -vf 0.1 (✔)
  : Specifies the volume fraction of inclusion.
* -mf 0.0 (✔)
  : Specifies a mass fraction from  which the volume fraction is calculated.   The calculation requires the density of the supporting matrix.
* -matrix ptfe 
  : The supporting matrix is defined by  the string.  Options are “ptfe”, “kbr”,  “ujol”, “air”, “vacuum”, “ldpe”, “mdpe”, “hdpe”.  If the matrix is given in this way both the  density and the permittivity of the supporting matrix are defined.  Alternatively the density and dielectric  options can be used.
* -density 2.2
  : Defines the density of the  supporting matrix.
* -dielectric 2.0
  : Defines the dielectric of the  supporting matrix.
* -LO h k l (✔)
  : The frequencies corresponding to the longitudinal optic modes with a k vector direction (h k l) are calculated  using Equations [#eq-LODynamicalMatrix].
* -LO_cart x y z (✔)
  : As above but for Cartesian  directions
* -sigma 5.0 
  : Specifies the damping factor, σ,  for all modes in cm^-1^, as used in Equation [#eq-permittivity].
* -mode_sigma k σ (✔)
  : The k’th mode is assigned a specific σ (cm-1).
* -vmin 0.0
  : Thee starting wavenumber (cm^-1^)  for the frequency range.
* -vmax 300.0
  : The final wavenumber (cm^-1^)  for the frequency range.
* -i 0.2
  : The increment used to cover  the frequency range (cm^-1^).
* -plot absorption (✔)
  : Plot types  are specified by the string and they can be  ‘absorption’, ‘molar_absorption’, ‘real’ or ‘imaginary’.
* -excel filename.xlsx
  : Writes the results to an excel  spread sheet with the specified name.
* -csv filename.csv
  : Output is sent to the specified comma separated  file.
* -csv_ext filename
  : Output is sent to 3 comma separated  files; filename_command.csv, filename_frequency.csv and filename_spectrum.csv.
* -print
  : Additional output is provided from  the QM or MM calculation.
* -ignore k (✔)
  : Ignore the kth mode (any mode below 5cm^-1^ is ignored automatically).
* -mode k (✔)
  : Only use the kth mode in the  calculation of the permittivity.
* -threshold 1.0e-10 5.0
  : The modes selected for inclusion in  the absorption calculation have to have an IR intensity greater than 1.0E-10 and a frequency greater than 5.0 cm^-1^.
* -eckart
  : The translational modes are  projected out of the hessian before diagonalisation.
* -hessian symm
  : This option can specify either “crystal” or  “symm”.  In the case of “crystal” the  hessian is symmetrised using the same algorithm as Crystal14.
* -optical z1 z2 z3
  : z1,z2 and z3 define the diagonal of  the optical permittivity tensor.
* -optical_tensor z1 z2 ..z9
  : z1,..9 define the full optical permittivity tensor.
* -masses average
  : Specified the mass definition to use, which can be either “program”, “average” or “isotopic”, meaning that the masses used in the calculation of the frequenciesare either taken from the QM program or are the average of the isotope abundances or are the most abundant isotope mass.
* -mass element mass (✔)
  : The atomic mass of the element is set to mass.  This can be used to explore the effect of isotope substitution on the calculated frequencies.
* -processors int
  : The number of processors to be used  in the calculation can be defined.  By  default all available processors are used.
* -molesof cells \[number_of_atoms_per_molecule\]
  : The default for the calculation of molar concentration is “cells”.  Other options are “atoms” or  “molecules”.  If ‘molecules’ is  specified then the number of atoms in a molecules must be provided.
* -size  0.00001 \[sigma\]
  : Modifies the polarisability (Eq.34) for spherical particles which incorporates the radius of the particle in microns(real).   It is also used to specify the dimension of the spherical particles for the Mie method.  It is also possible to specify a size distribution in which case the first number is the mean of the log distribution (in microns) and the second is its width (in log(microns)).
~

The shape options; ellipse, slab and needle, specify a unique axis \[hkl\] using the crystal axes of the unit cell. PDielec transforms these to a cartesian coordinate system using the unit cell lattice vectors. In the case of a slab morphology the unique direction is normal to the surface specified by its Miller indices \(hkl\. The definitions of the various depolarisation tensors are indicated in Table below.

~TableFigure {#tab-depolarisation; caption: "Definitions used for the depolarisation tensor";}

| Shape     | Depolarisation Tensor                                 |
| ----------| -----------------------------------------------------|
| Sphere    | $\tensorbf{L} =\frac{1}{3}\left( \fieldbf{V}_1 \fieldbf{V}_1^T + \fieldbf{V}_2 \fieldbf{V}_2^T+\fieldbf{V}_3 \fieldbf{V}_3^T \right)$                            |
| Slab      | $\tensorbf{L}=\fieldbf{V}_1 \fieldbf{V}_1^T$        |
| Needle    |$\tensorbf{L}=\frac{1}{2}\left( \fieldbf{V}_2\fieldbf{V}_2^T + \fieldbf{V}_3 \fieldbf{V}_3^T \right)$|
| Ellipsoid | $\tensorbf{L}=a\fieldbf{V}_1 \fieldbf{V}_1^T + b\fieldbf{V}_2 \fieldbf{V}_2^T + b\fieldbf{V}_3 \fieldbf{V}_3^T$                       |

~                                                    

The three directions defined by $\fieldbf{V}_1, \fieldbf{V}_2 \text{and} \fieldbf{V}_3$ are mutually orthogonal cartesian vectors calculated from \[hkl\] for an ellipse  or needle and \(hkl\) for a slab. In the case of a slab, needle or ellipsoid, ${\bar{V}}_{1}$ defines the unique direction and the other vectors are orthogonal to it. For the case of an ellipsoid, the parameters *a* and *b* in Table depend on the ratio, $z$, of the length of unique axis length over the length of an axis perpendicular to it [@Sihvola].

For z \> 1 the ellipsoid is prolate;
$e = \sqrt{1 - z^{- 2}},\ a = \frac{\left( 1 - e^{2} \right)}{2e^{3}}\left( \log\frac{1 + e}{1 - e} - 2e \right),\ b = \frac{1}{2}\left( 1 - a \right)$

For z \< 1 the ellipsoid is oblate
$e = \sqrt{z^{- 2} - 1},\ a = \frac{\left( 1 + e^{2} \right)}{e^{3}}\left( e - \operatorname{}e \right),\ b = \frac{1}{2}\left( 1 - a \right)$

From an experimental point of view it is often convenient to use a mass fraction rather than a volume fraction to indicate the amount of dielectrically active material present.  PDielec allows mass fractions to be specified instead of a volume fraction, but this requires that the density of the supporting matrix is known. For convenience the package has a small database of the common supporting materials shown in Table below.  These can be specified through the -matrix option. In the case that the properties of the support material are different the properties can be defined instead with the -dielectric and -density options. 

~TableFigure {#tab-matrix; caption: "Physical properties of matrix materials in PDielec";}
| Name   | Density | Permittivity | Description                 |
| ------ |: -------| ------------ | --------------------------- |
| ptfe   | 2.2     | 2.0          | polytetrafluorethylene      |
| air    | 0.0     | 1.0          | air                         |
| vacuum | 0.0     | 1.0          | vacuum                      |
| kbr    | 2.75    | 2.25         | potassium bromide           |
| nujol  | 0.838   | 2.155        | Nujol                       |
| hdpe   | 0.955   | 2.25         | high density polyethylene   |
| mdpe   | 0.933   | 2.25         | medium density polyethylene |
| ldpe   | 0.925   | 2.25         | low density polyethylene    |
~

The optical permittivity is normally calculated by the QM or MM program concerned. However, as this property reflects the electronic contribution to the permittivity at zero frequency, unless there is some treatment of electrons by the shell model, then in MM calculations the optical permittivity needs to be defined through the command line options -optical or -optical_tensor.  Unlike the other methods 'mie' method cannot work with particles of zero radius. All methods therefore use a default size of 10^-12^ μm. The Mie approach is only valid for dilute dispersions and for spherical particles. However, if other shapes are specified the Mie method will still be called and the results will be applicable to spheres of the specified size.

## Parallelization and threads

To improve the performance of the program python parallelization has been used to parallelize over the frequencies, shapes and methods. By default this parallelization spawns additional Python executables, depending on the number of cores available.

In addition to this form of parallelization the NUMPY library can use multi-threaded BLAS. NUMPY can be compiled with several different BLAS libraries, including; MKL (from Intel), OPENBLAS and ATLAS, or indeed no optimized BLAS library is necessary. To explore the BLAS version compiled with your version of NUMPY please use the test_numpy_1, test_numpy_2 and test_numpy_3 scripts included in the Python/ subdirectory. If you are using MKL or OPENBLAS, the number of threads being used needs to be set with the MKL_NUM_THREADS or OPENBLAS _NUM_THREADS environment variable (sometimes OMP_NUM_THREADS is also used). Use the test routines to determine the optimum for your system. Experience indicates that no performance benefit is obtained with more than two threads. 

In the case of Phonopy the dynamical matrix is read from a yaml file. This has been found to be very slow unless the C parser is used. If the C parser is not available a warning is issued and the program reverts back to the Python parser.

Finally the use of non-standard BLAS libraries seems to cause problems with the affinity settings for the multiprocessing. It has been noticed that the parallel processes can all end up executing on the same processor. In order to prevent this, before executing the pdielec and preader scripts it may be necessary to include;

        export OPENBLAS_MAIN_FREE=1

For some reason this also works if the MKL library is being used.

## Performance

PDielec has been written to make use of multiprocessor computers. On a 4 processor machine the speed-up is nearly linear up to 4 processors, as can be seen in the Figure below.

~ Figure {#fig-speedup; caption:"Speed-up on a four processor workstation"}
![img-speedup]
~
[img-speedup]: ./Figures/SpeedUp.png "Speed-up" {width:auto; max-width:90%}

## Example command line uses of PDielec

         pdielec -program vasp OUTCAR -method ap -method maxwell \
                 -sphere -plate 0 0 1 -needle 0 0 1 -LO 0 0 1

This performs a calculation using the Averaged-Permittivity and Maxwell-Garnett mixing rules for spherical particles, plate-like particles with a surface (001) and needle-like particles with a unique
direction lying along the \[001\] direction. The supporting matrix is taken to be PTFE and the default volume fraction (10%) is used. The results of a VASP calculation are stored in the current directory.  There is no absorption output from this command as neither the -plot nor the -csv options were specified. The output includes the calculation of the LO modes along the (001) direction.

        pdielec -program castep phonon -vmin 300 -vmax 800 \
                -sphere -dielectric 3 -vf 0.1 -vf 0.2 -sigma 10 -csv mgo.csv

This performs a calculation for spherical particles varying the frequency from 300 to 800 cm^‑1^, the permittivity of the supporting media is 3, two volume fractions are considered and a damping factor of 10 cm^-1^ is used. The results of a CASTEP calculation with the seed-name "phonon" are analysed and the results stored in mgo.csv for further analysis using a spreadsheet. In this example a Maxwell-Garnett mixing rule is used by default. If visual inspection of the results is required then

        pdielec -program castep phonon -vmin 300 -vmax 800 \
                -sphere -dielectric 3 -vf 0.1 -vf 0.2\
                -sigma 10 -csv mgo.csv -plot molar_absorption

will perform the same calculation but a graph showing the molar absorption coefficients will be displayed.

        pdielec -program gulp calcite.gout -matrix hdpe \
                -method ap -method maxwell -sphere -plate -1 -1 -2 \
                -vmax 2000 -mf 0.1 calcite.gout -csv calcite.csv

This command performs a calculation of the absorption spectrum resulting from a GULP calculation. The supporting matrix density and permittivity are those of high density polyethylene, the frequency range is 0 to 2000 cm^-1^, the mass fraction considered is 10%, the mixing rules used are Averaged-Permittivity and Maxwell-Garnett. Spheres and plates with the ($\bar{1}\bar{1}\bar{2})$ surface are considered.

        pdielec -program vasp OUTCAR -method mie -sphere \
                -mf 0.1 -size 0.1 -size 1.0 -size 1.0 -csv results.csv

This command performs a calculation of the absorption spectrum resulting from a VASP calculation using the Mie method for sphere with several particles sizes.

        pdielec -sphere -mf 0.1 -program experiment experiment.expt\
                -size 0.1 -size 1.0 -size 1.0 -excel results.xlsx

This command performs a calculation of the absorption spectrum resulting from the data stored in the experiment.expt file. Maxwell-Garnett calculations are performed with 3 different sized spheres and the results stored in a Excel file.

## Contents of the csv output file

If a csv output file is requested the file will contain the command used to perform the calculation. A brief summary is given of each active infrared mode; including the mode number, frequency, intensity, integrated molar absorption coefficient, its peak height (calculated from the intensity and damping factor) and the damping parameter used in the calculation. Following this is a table with a column for frequency followed by columns containing the real and imaginary permittivities, the absorption and molar absorption coefficients at each frequency.

[BIB]
