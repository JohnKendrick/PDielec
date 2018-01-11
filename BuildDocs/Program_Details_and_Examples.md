---
Title: PDielec The Calculation of Infrared and Terahertz Absorption for Powdered Crystals
Author: John Kendrick^1^\footnote{john@kendrick.me.uk}, Andrew Burnett^1^
Abstract: The Python package PDielec is described, which calculates the infrared absorption characteristics of a crystalline material supported in a non-absorbing medium. PDielec post processes solid state quantum mechanical and molecular mechanical calculations of the phonons and dielectric response of the crystalline material. Packages supported are; Abinit, Castep, Crystal14, Gulp, QuantumEspresso, Phonopy and VASP. Using an effective medium method, the package calculates the internal electric field arising from different particle morphologies and calculates the resulting shift in absorption frequency and intensity arising from the coupling between a phonon and the internal field. The theory of the approach is described, followed by a description of the implementation within PDielec. For the specific case of a low concentration of spherical particles, calculations based on Mie scattering allow the exploration of particle size effects. A section providing several examples of its application is given.
---
PDielec: The Calculation of Infrared and Terahertz Absorption for Powdered Crystals
===================================================================================

John Kendrick^1^ and Andrew Burnett^1^

^1^School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom.

$$
\label{eq:preamble}\notag \newcommand{\water}{H_{2}O}
\newcommand{\tensor}[1]{\overline{\overline{#1}}}
\newcommand{\tensorbs}[1]{\overline{\overline{\boldsymbol{#1}}}}
\newcommand{\tensorbf}[1]{\overline{\overline{\mathbf{#1}}}}
\newcommand{\fieldbf}[1]{\overline{\mathbf{#1}}}
$$

ABSTRACT
========

The Python package PDielec is described, which calculates the infrared absorption characteristics of a crystalline material supported in a non-absorbing medium. PDielec post processes solid state quantum mechanical and molecular mechanical calculations of the phonons and dielectric response of the crystalline material. Packages supported are; Abinit, Castep, Crystal14, Gulp, QuantumEspresso, Phonopy and VASP.
Using an effective medium method, the package calculates the internal electric field arising from different particle morphologies and calculates the resulting shift in absorption frequency and intensity arising from the coupling between a phonon and the internal field. The theory of the approach is described, followed by a description of the implementation within PDielec. For the specific case of a low concentration of spherical particles, calculations based on Mie scattering allow the exploration of particle size effects. A section providing several examples of its application is given.

INTRODUCTION
============

The molecular and solid state quantum mechanical (QM) calculations of response properties such as the frequencies and intensities of infrared (IR) and terahertz (THz) radiation absorption has become generally available in many molecular and solid state computer programs. A common approach is to assume the harmonic approximation and calculate the mass weighted force constant matrix (for molecules) or the dynamical matrix at the gamma point (for periodic solids). Diagonalisation of the matrix gives the frequencies for absorption and the normal modes (molecules) or phonon displacements (periodic solids). The calculation of the absorption intensity for each mode requires the calculation of the change in dipole moment caused by the displacement of the atoms for that mode. For solids where there is a large separation of charge, there can be a large coupling between a phonon mode and the internal field within a particle resulting from its morphology. 

This paper describes the PDielec package, which is written in Python and post processes the output of solid state quantum mechanical and molecular mechanics (MM) based codes such as VASP [@Hafner2008c], CASTEP [@Clark2005d], CRYSTAL [@Dovesi2014], Abinit [@Gonze2016] , Quantum Espresso [@Giannozzi2009], Phonopy [@Togo2015], and GULP [@Gale2003] to predict the infrared absorption of crystalline insulator materials whose crystal size is small compared with the wavelength of the absorbing radiation. The package is suited for the calculation of the complex, frequency dependent permittivity and its associated absorption of infrared radiation for a finely ground crystalline material dispersed in a low loss dielectric medium such KBr or Polytetrafluoroethylene (PTFE). A particular feature of the program is its ability to take into account the constant permittivity of the supporting medium and the particle shape of the material of interest through an effective medium theory. The paper outlines the theory used by the program and gives some examples of the application of the program for ionic and molecular materials.

THEORY
======

Equations $\eqref{eq:beer1}$ and $\eqref{eq:beer2}$ describe Beer-Lambert's law [@Bertie2006] where $\alpha$ is the (decadic) absorption coefficient (usually given in cm^-1^), $l$ and $l_0$ are the intensities after and before absorption respectively and $d$ is the path length.
$$
\frac{I}{I_{0}} = 10^{- \alpha d}\label{eq:beer1}
$$

$$
log\left( \frac{I}{I_{0}} \right) = -\alpha d \label{eq:beer2}
$$

It is common, especially in the chemistry community, when reporting infrared spectra to use a decadic molar absorption coefficient ($a$), which has units of Lmol^-1^cm^-1^. The relationship between the absorption coefficient and the molar absorption coefficient [@Berti2006] is;
$$
\alpha = aC
$$


where *C* is the concentration of the absorbing species.

Molecular Approach to Absorption Intensity
------------------------------------------

For molecules the transition intensity $l_k$ of the $k^{th}$ mode (calculated from the change in dipole moment along the mode displacement) can be converted to an integrated molar absorption coefficient, $A_k$, which can then be more readily compared with experiment. The theory for this is described by Wilson, Decius and Cross [@Wilson1955] and results in expressions such as the two equations below $\eqref{eq:Absorption1}$ and $\eqref{eq:Absorption2}$.  The first expression shows the relationship between the integrated molar absorption coefficient and the transition intensity and uses the number of molecules per unit volume ($N$), the velocity of light ($c$) and the degeneracy of the mode ($g_k$). The second expression shows the appropriate conversion factors if the units for the integrated molar absorption coefficient are L mol^‑1^cm^‑2^ (1 L mol^-1^cm^-2^ = 0.01 km mol^-1^) and the units for the transition intensity are D^2^ Å^-2^ amu^-1^, where D represents the Debye unit of dipole  moment and amu is an atomic mass unit. The factor log~e~10 arises due to the choice of a decadic Beer's law.
$$
A_k =  \frac{N\pi}{3c^{2}\log_e10}g_kI_k\label{eq:Absorption1}
$$

$$
A_k =  \frac{Na\pi}{3000c^{2}2.302585}g_kI_k  = 4225.6I_k \label{eq:Absorption2}
$$

The derivation of the above expressions assumes that the rotational levels are not quantised and that the vibrational levels are thermally occupied according to a Boltzmann distribution. In order to use the calculated molecular intensities to predict a spectrum it is usual to assume [@Wilson1955] that each transition is associated with a Lorentzian line shape with a full width at half maximum (FWHM) of $\sigma_k$. It is common, when reporting comparison between theoretical and experimental spectra, to assume that the line widths are the same for all modes [@Juliano2013,@Burnett2013].  Recent work on terahertz absorption in crystalline pentaerythritol tetranitrate (PETN) using molecular dynamics calculations [@Pereverzev2011] in combination with the direct calculation of the cubic anharmonic couplings of the normal modes [@Pereverzev2014] has shown that the FWHM of the intense absorptions may vary between 10 and 25 cm^-1^. Assuming a Lorentzian line shape, the molar absorption coefficient for the $k^{th}$ mode at wavenumber, ${\overline{\nu}}_{k}$, can be written as a function of frequency or wavenumber ($\overline{\nu}$);
$$
a_k(\overline{\nu}) = \frac{2A_k}{\pi}\frac{\sigma_k}{4\left( \overline{\nu} - {\overline{\nu}}_k \right)^{2} + \sigma_k^2} \label{eq:lorentzian}
$$

$$
a_k^{\max} = \frac{2A_k}{\pi\sigma_k}
$$

The maximum height of the Lorentzian, $a_{k}^{\max}$ clearly depends upon the value of $\sigma_k$. As can be seen in Equation $\eqref{eq:integratedmolarintensity}$, the choice of normalisation for the Lorentzian means that integration of the molar absorption coefficient over wavenumber returns the integrated molar absorption coefficient and a sum over all the bands provides the total molar absorption coefficient $a^{\text{mol}}(\overline{\nu})\$ as a function of wavenumber, calculated from the intensities of each band.  Equation $\eqref{eq:molarabsorption}$ shows the relationship between the absorption and the molar absorption coefficients. $$C$$ is the concentration usually expressed in mol L^-1^.
$$
A_k = \int{a_k(\overline{\nu})d\overline{\nu}} \label{eq:integratedmolarintensity}
$$

$$
a^{\text{mol}}(\overline{\nu}) = \sum_k{a_k(\overline{\nu})}
$$

$$
\alpha^{\text{mol}}(\overline{\nu}) = Ca^{\text{mol}}(\overline{\nu}) \label{eq:molarabsorption}
$$

A comment should be made about the various units which can be used for these quantities. A common unit for the transition intensity is D^2^ Å^-2^ amu^-1^, another is km mol^-1^. However, it should be pointed out that strictly speaking the latter unit refers to the integrated molar absorption coefficient as defined above in Equation $\eqref{eq:integratedmolarintensity}$ and therefore relies on the assumptions made in its derivation. ( 1 D^2^ Å^-2^ amu^-1^ is equivalent to 42.256 km mol^-1^ ).

Solid State Approach to Absorption Intensity
--------------------------------------------

The optical properties of a solid are determined by its complex, frequency dependent relative permittivity,
$\boldsymbol{\varepsilon}(\overline{\nu})$, and in particular the real and imaginary components, $\boldsymbol{\kappa}$ and $\mathbf{n}$, of the complex refractive index, $\mathbf{N}(\overline{\nu})$ where;
$$
\tensorbf{N}(\overline{\nu})^{2} = {\tensorbs{\varepsilon}}(\overline{\nu})
$$

$$
{\tensorbf{N}}\left( \overline{\nu} \right) = {\tensorbf{n}}\left( \overline{\nu} \right) + i{\tensorbs{\kappa}}\left( \overline{\nu} \right) \
$$

The intensity of absorption is given by the imaginary component of the refractive index which for an isotropic material is [VanDeHulst1981]^10^;
$$
I = I_{0}e^{- \frac{4\pi\kappa(\overline{\nu})d}{\lambda}} \\   
  I = I_{0}e^{- 4\pi\overline{\nu}\kappa(\overline{\nu})d} \
$$

$$
- ln\left( \frac{I}{I_{0}} \right) = 4\pi\overline{\nu}\kappa(\overline{\nu})d \\              
- log\left( \frac{I}{I_{0}} \right) = 4\pi\overline{\nu}\kappa(\overline{\nu})d \cdot log(e)
$$

Comparison with the definition of the absorption coefficient from Beer-Lambert's law, Equation $\eqref{eq:beer1}$, and using Equation $\eqref{eq:molarabsorption}$ gives;
$$
\alpha^{\text{sol}}(\overline{\nu}) = 4\pi\overline{\nu}\kappa(\overline{\nu}) \cdot log(e) \
$$

$$
a^{\text{sol}}(\overline{\nu}) = \frac{\alpha^{\text{sol}}(\overline{\nu})}{C}
$$

Since the refractive index is dimensionless, the absorption coefficient, $\alpha^\text{sol}$ is specified in cm^-1^. The superscripts 'sol,' for solid, and 'mol,' for molecular, are used here to distinguish between the two methods of calculating the absorption ($\alpha$) and molar absorption coefficients ($a$).  In the calculation of the imaginary component of the refractive index it is necessary to choose the solution which gives a positive value. This is consistent with the Kramers-Kronig relationship between the real and imaginary components [@Wooten1972].

In order to calculate the relationship between absorption and molar absorption coefficients it is necessary to know the concentration. For solid state calculations the required unit is; moles of unit cells per litre. One of the drawbacks of this molar absorption coefficient unit is that the number of molecules in a unit cell can change depending on whether a supercell, primitive or non primitive unit cell is being used. A more natural unit would be to use a mole of formula units, or a mole of molecules. In order to aid comparison between such calculations PDielec is able to calculate concentration in both moles of atoms and moles of molecules. However for the rest of this paper Equation $\eqref{eq:concentration}$ will be used, where *V* is the volume of the unit cell, and therefore the concentration $C$ is moles of unit cell/litre.
$$
C = \frac{f \cdot1000cm^{3}}{VN_{a}} \label{eq:concentration}
$$
The volume fraction, $f$, of the dielectric material in a supporting matrix of non-absorbing material is included in the expression for the concentration as it will be useful when the theory for mixtures is developed.

For a periodic system the permittivity tensor can be calculated as a sum over Lorentz oscillators, incorporating an imaginary loss component through the damping factor *σ~k~* [@Gonze1997].  The frequencies of the oscillators are the transverse optic (TO) phonon frequencies of the system.
$$
{\tensorbf{\varepsilon}}(v) = {{\tensorbf{\varepsilon}}}_{\infty} + \frac{4\pi}{V}\sum_{k}\frac{{{\tensorbf{S}}}_{k}}{\nu_{k}^{2} - \nu^{2} - i\sigma_{k}\nu}\label{eq:permittivity}
$$

$$
\tensorbf{S}_{k} = {\overline{\mathbf{Z}}}_{k}{\overline{\mathbf{Z}}}_{k}^{T}\label{eq:oscillatorstrength}
$$

$$
\overline{\mathbf{Z}}_{k} =  \sum_{a}{\tensorbf{Z^{a}}}{\overline{\mathbf{U^{a}}}}_{k}\label{eq:borncharges}
$$

$$
\tensorbf{D}{\overline{\mathbf{U}}}_{k} = \Lambda_{k}{\overline{\mathbf{U}}}_{k}\\
\nu_{k}^{2} = \Lambda_{k}\label{eq:eigenvalues}
$$

$$
I_{k} = tr\left( \tensorbf{S}_{k} \right) \label{eq:intensity}
$$

$V$ is the volume of the unit cell, ${\tensorbf{S}}_{k}$ is the dipole oscillator strength tensor for the $k_{th}$  transition, with a TO frequency of $\nu_k$ and $\tensor{\varepsilon}_{\infty}$ the optical permittivity tensor, which represents the electronic contribution to the permittivity. The intensity of a transition, $I_k$, is given by the trace of the oscillator strength tensor, Equation $\eqref{eq:intensity}$.  The damping factor $\sigma_{k}$ removes any discontinuities at the TO frequencies. Since the oscillator strengths and phonon frequencies can be calculated routinely in solid state quantum mechanical packages, the calculation of the frequency dependent complex permittivity using Equation $\eqref{eq:permittivity}$ is straightforward. In some cases, using Equations $$\eqref{eq:oscillatorstrength}$$ and $\eqref{eq:borncharges}$, PDielec calculates the oscillator strengths from the Born charge
matrix for atom $a$ and its contributions to the $k_{th}$ phonon mode [@Gonze1997].  As shown in Equation $\eqref{eq:eigenvalues}$, at the $\Gamma$ point the $k_{th}$ phonon mode is described by the eigenvector, $\overline{\mathbf{U}}_{k},$ and eigenvalue, $\Lambda_{k},\ $of the mass weighted, dynamical matrix, $\tensorbf{D}$, which is a 3Nx3N matrix, where N is the number of atoms in the unit cell. The eigenvalues are the squared frequencies of the phonon modes, Equation $\eqref{eq:eigenvalues}$.  The displacement of each atom in the $k_{th}$ is proportional to $m_{a}^{-1/2}$, where $m_a$ is the mass of atom $a$. The dynamical matrix has 3N eigenvectors and eigenvalues, of which three should be zero due to translational invariance. If there are any negative eigenvalues the system is unstable to some displacement and therefore not at an energy minimum.

For ionic systems it is common practice in solid state QM and MM programs to include a long wave-length, non-analytic correction to the mass weighted dynamical matrix at the $\Gamma$ point, which describes the coupling of the longitudinal optic (LO) modes to the induced field resulting from the vibration. This may be written for atoms $s$ and $t$ and their Cartesian components $\alpha$ and $\beta$ as [Gonze1997];
$$
\left( \tensorbf{D^{LO}}_{\mathbf{q \rightarrow 0}} \right)_{s,\alpha;t,\beta}\mathbf{=}\left( \tensorbf{D} \right)_{s,\alpha;t,\beta}\mathbf{+}\frac{4\pi}{V\sqrt{M_{s}M_{t}}}\mathbf{\ }\frac{\left( {\overline{\mathbf{q}}}^{\mathbf{T}}{{\mathbf{\ }\tensorbf{Z}}}_{s} \right)_{\alpha}\left( {\overline{\mathbf{q}}}^{\mathbf{T}}{{\mathbf{\ }\tensorbf{Z}}}_{t} \right)_{\beta}}{{\overline{\mathbf{q}}}^{\mathbf{T}}\mathbf{\cdot}{{\tensorbf{\varepsilon}}}_{\infty} \cdot \overline{\mathbf{q}}}
$$
The mass weighting has been incorporated through the mass of the atoms, $M_s$ and $M_t$.  The correction depends upon the direction, $\overline{\mathbf{q}}$**,** that the long wave-length limit is
approached.  Diagonalisation of the corrected matrix gives the squared frequencies of N-1 LO modes and 2N-2 TO modes, Equation $\eqref{eq:eigenvalues}$.  In some of the examples given below the LO frequencies will be given for comparison with the TO frequencies.

Effect of Particle Shape on Infrared Absorption
-----------------------------------------------

It has long been recognised that, especially for ionic materials, the local field within a crystal and its coupling with the transverse optical phonons has an important effect on the position and intensity of the absorption. Fröhlich [@Frohlich1948] was one of the first to point out that the frequency of absorption of a small ionic sphere embedded in a low dielectric medium is shifted to lie between the transverse and longitudinal optical frequencies of the material making up the sphere.

In the development of the theory used in PDielec an important assumption is that the particle size of the crystallites in the sample is small compared with the wavelength of light.  Using this approach Genzel and Martin [@Genzel1972a] were able to explain the observed infrared absorption of small spheres of MgO crystallites and the effect of the permittivity of the supporting medium on the spectrum.  Studies of the infrared absorption by small particles of α-Fe~2~O~3~ using an effective medium theory and an absorption/scattering theory [@Serna1987, @Iglesias1990] showed that in order to fit the experimental spectra it was necessary to adjust not only the damping factors in Equation $\eqref{eq:permittivity}$ but also the permittivity of the matrix and the volume fraction of the dielectric medium. The latter was used to account for aggregation effects as the volume fraction increased. It was also shown that effective medium theories were only applicable for particles smaller than the wavelength of light. For larger particles the
scattering from the particles becomes increasingly important.

More recently Balan and others in a series of papers [@Balan2010b,@Balan2008b,@Fourdrin2009]^20-24^ used density functional calculations together with an effective medium theory to calculate the infrared absorption of several minerals incorporating information about the crystallite shape.  In an experimental and theoretical study of irradiated kaolinite[@Fourdrin2009] it was shown that exposure to radiation resulted in shifts in the infrared spectrum which could be accounted for by increasing the polarisability of the particles through
an increase in the optical permittivity tensor.

The underlying theory adopted by PDielec is based on similar premises to the work described above, namely that the dielectric response of small spherical, ellipsoidal, slab-like or needle-like crystallites randomly distributed in a non-absorbing medium such as PTFE, KBr or Nujol, is the same as that of an effective medium material whose frequency dependent dielectric response can be calculated from the frequency dependent permittivity tensor of the crystal (as calculated by solid state QM or MM calculations), the shape of the crystallites and the permittivity of the non-absorbing medium (taken to be a constant over the frequency range of interest).

The development of the theory reported here closely follows the work of Sihvola [@Sihvola].  It will be assumed that the inclusion particles, which may be non-isotropic, ellipsoidal (including spherical, needle-like and plate-like), are randomly orientated in an embedding, non-absorbing medium such as PTFE, KBr or Nujol. It should be emphasized that whilst PDielec can take account of particle shape, particle and matrix permittivity there are many additional aspects of infrared absorption which need to be considered when comparing calculated and experimental results. Most notable of these are; the coupling between phonons and mobile electrons or holes (so called phonon-polariton coupling)[@Ruggiero2015], the scattering which starts to dominate as the particles get larger [@Fourdrin2009] and the agglomeration of particles as the volume fraction increases.

### The polarisability of an isolated particle

Figure $@fig:Polarisation$ shows a schematic of the field and polarisation inside an inclusion with non-isotropic permittivity $\tensor{\varepsilon}_{i}$ embedded in a supporting medium with permittivity, $\varepsilon_e$.  The internal field within the inclusion is indicated by $\fieldbf{E}_{i}$ the external, applied field is indicated by $\fieldbf{E}_{e}$ and the induced polarisation in the inclusion is shown by $\fieldbf{P}$.

![Schematic showing the field and polarisation inside an inclusion with non-isotropic permittivity ${\tensorbf{\varepsilon}}_{i}$ embedded in a supporting medium with permittivity$\varepsilon_e$. The internal field within the inclusion is indicated by $\fieldbf{E}_i$ , the external, applied field is indicated by $\fieldbf{E}_e$ and the induced polarisation in the inclusion is shown by $\fieldbf{P}$.](./Images/Figure01_Polarisation_Picture.tiff){#fig:Polarisation}

The electric field internal to the inclusion gives rise to a polarisation density which is no longer necessarily aligned with the field because the material is non-isotropic. The polarisation density in the inclusion can be expressed as the tensor product of the permittivity contrast between the inclusion and the supporting medium and the (as yet unknown) internal field.
$$
\fieldbf{P} = \left( \tensorbf{\varepsilon}_{i} - \varepsilon_{e}{\tensorbf{1}} \right){\fieldbf{E}}_{i}\label{eq:PolarisationDensity}
$$
For any ellipsoidal shape (including sphere, slab and needle) with volume *V*, the polarisation density throughout the particle is uniform and integrating over all space gives the field induced dipole moment of the inclusion, $\fieldbf{p}$.
$$
\fieldbf{p} = V\fieldbf{P} = V\left( \tensorbf{\varepsilon}_{i} - \varepsilon_{e}\tensorbf{1} \right)\fieldbf{E}_{i}\label{eq:polar1}
$$
The dipole and the external field, ${\fieldbf{E}}_{e}$, are related by the polarisability tensor, $\tensorbf{\alpha}$.
$$
\fieldbf{p} = \tensorbf{\alpha}\fieldbf{E}_{e}\label{eq:polar2}
$$
Equations $\eqref{eq:polar1}$ and $\eqref{eq:polar2}$ allow the determination of the polarisability, once
the field internal to the inclusion has been expressed in terms of the shape of the inclusion and its permittivity. The polarisation within the inclusion gives rise to a depolarisation field, $\fieldbf{E}_{d}$, which depends on the shape of the inclusion through the symmetric and unit-trace depolarisation tensor, $\tensorbf{L}$.
$$
\fieldbf{E}_d = - \frac{1}{\varepsilon_e}\tensorbf{L} \cdot \fieldbf{P}\label{eq:DepolarisationField}
$$
The internal field is the sum of the external field and the depolarisation field.
$$
\fieldbf{E}_{i} = \fieldbf{E}_{e} + \fieldbf{E}_{d}\label{eq:InternalField}
$$
The depolarisation tensor acts as a projection or screening operator describing the effect of the geometry of the inclusion on the depolarisation field which results from its polarisation. For instance, in the case of a needle, only polarisation perpendicular to the needle axis contributes to the depolarizing field, whilst for a slab only polarization perpendicular to the slab face may contribute. Similarly for a sphere, all directions contribute and so the depolarisation matrix is diagonal with a value of 1/3 for each diagonal element, as the trace of the depolarisation tensor must be 1. It follows from Equations $\eqref{eq:PolarisationDensity}$, $\eqref{eq:DepolarisationField}$ and $\eqref{eq:InternalField}$ that;
$$
\fieldbf{E}_i = \fieldbf{E}_e - \frac{1}{\varepsilon_e} \tensorbf{L} \left( \tensorbf{\varepsilon}_i - \varepsilon_e{\tensorbf{1}} \right){\fieldbf{E}}_i\label{eq:InternalField2}
$$
Rearrangement allows the internal field of the inclusion to be expressed in terms of the known permittivities, the shape of the inclusion and the external field.
$$
\fieldbf{E}_i \left( \varepsilon_e{\tensorbf{1}} + \tensorbf{L}  \left( \tensorbf{\varepsilon}_i - \varepsilon_e{\tensorbf{1}} \right) \right) = \varepsilon_e\fieldbf{E}_e
$$

$$
\fieldbf{E}_i = \varepsilon_e{\fieldbf{E}}_e\left( \varepsilon_e\tensorbf{1} + \tensorbf{L} \left( \tensorbf{\varepsilon}_i - \varepsilon_e \tensorbf{1} \right) \right)^{- 1}
$$

Substituting the internal field expression, $\eqref{eq:InternalField2}$, into Equation $\eqref{eq:polar1}$ for the dipole moment and requiring the dipole moments calculated using the polarisation density to equal those calculated from the polarisability allows the polarisability to be written;
$$
\tensorbs{\alpha} = V\varepsilon_e\left( \tensorbs{\varepsilon}_i - \varepsilon_e\tensor{1} \right)\left( \varepsilon_e\tensorbf{1} + \tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_e\tensorbf{1} \right) \right)^{- 1}\label{eq:polarisation}
$$
Although it has not been specified explicitly the permittivity of the inclusion, and therefore the polarisability tensor, are frequency dependent through the oscillator strengths of each phonon mode contributing to the permittivity according to Equation $\eqref{eq:permittivity}$.  The calculation of the complex, frequency dependent polarisability tensor of the composite material is the key step in the calculation of its effective permittivity.

### The Effective permittivity of a mixture

To extend this approach to include the effect of a number of inclusions
we need to introduce the concept of an effective permittivity, $\tensorbf{\varepsilon}_{eff}$, which describes the behaviour of an average field, $\left\langle \fieldbf{E} \right\rangle$, where the angle brackets indicate an average over a volume of the composite material. It is required that the average electric flux density $\left\langle \fieldbf{D} \right\rangle$ is the same in the effective medium as in the composite medium;
$$
\left\langle \fieldbf{D} \right\rangle = \tensorbf{\varepsilon}_{eff}\left\langle \fieldbf{E} \right\rangle = \varepsilon_e\left\langle \fieldbf{E} \right\rangle + \left\langle \fieldbf{P} \right\rangle\label{eq:averageFluxDensity}
$$
The averaging is necessary because the polarisation within a given inclusion has an effect on the field in other inclusions. The local field in the cavity left by a single inclusion embedded in the average polarisation field is given by;
$$
\fieldbf{E}_L = \left\langle \fieldbf{E} \right\rangle + \frac{1}{\varepsilon_e}\tensorbf{L}\left\langle \fieldbf{P} \right\rangle\label{eq:localField}
$$
The local field 'excites' the inclusion resulting in a dipole moment $\mathbf{p}_{mix} $that is related to the polarisation through the number density of inclusions, $n$, and through the polarisability of the inclusion, which is already known from Equation $\eqref{eq:polarisation}$.
$$
\left\langle \fieldbf{P} \right\rangle = n{\fieldbf{p}}_{mix} = n\left\langle \tensorbf{\alpha}\fieldbf{E}_{L} \right\rangle
$$
The angle brackets around the product of the polarisability and the local field indicate that it is necessary to average the polarisation according to the distribution of alignments of inclusions. In this work it will be assumed that the inclusions are randomly aligned.  Substituting the expression for the local field, Equation $\eqref{eq:localField}$, gives;
$$
\left\langle \fieldbf{P} \right\rangle = \left( \tensorbf{1}
- \frac{n\left\langle \tensorbf{\alpha}\tensorbf{L} \right\rangle} {\varepsilon_e} \right)^{- 1}
n\left\langle \tensorbs{\alpha} \right\rangle\left\langle \tensorbf{E} \right\rangle \label{eq:averagePolarisation}
$$

### Mixing rules

There are many mixing rules which have been proposed to describe the homogenization of composite materials and a lot of work has been done in comparing their accuracy. Here two methods will be used. The first and the most commonly used method is the Maxwell-Garnett mixing rule [@Sihvola].  Indeed this has been implied by the use of Equation $\eqref{eq:averageFluxDensity}$ to define the effective permittivity. The other commonly used method is the Bruggeman mixing rule [@Sihvola], which differs substantially in the way the two components of the system are treated. It is usually stated that the Maxwell-Garnet mixing rule is good for low volume fractions of the inclusion and the Bruggeman approach should be better for higher volume fractions [@Giordano2003a].  In addition to these mixing rules one other approach will be described, namely the Averaged Permittivity (AP) mixing rule, which calculates the absorption spectrum ignoring the effects of the internal field on the absorption and can therefore be used as an indicator of the shifts in frequency and intensity which have occurred owing to the effect of particle shape.

### Maxwell-Garnett mixing rule

The Maxwell-Garnett approach for treating the properties of heterogeneous mixtures assumes that the average field and the average flux density result from volume fraction weighted sums. Substituting Equation $\eqref{eq:averagePolarisation}$ into Equation $\eqref{eq:averageFluxDensity}$ gives the Maxwell-Garnett effective permittivity;
$$
\tensorbs{\varepsilon}_{mg} = \tensorbf{1} + \left( \tensorbf{1} - \frac{n\left\langle {\tensorbs{\alpha}}{\tensorbf{L}} \right\rangle}{\varepsilon_e} \right)^{- 1}n\left\langle \tensorbs{\alpha} \right\rangle \label{eq:mg}
$$
The fact that the polarisability tensor has a volume term in it (Equation $\eqref{eq:polarisation}$) means that the terms in Equation $\eqref{eq:mg}$ containing $n\tensorbs{\alpha}$ depend on the volume fraction *f*. Although written as a tensor, because the assumption has been made that the inclusions are randomly orientated, the effective permittivity has to be diagonal with equal complex values. Since the polarisability is complex and frequency dependent the effective permittivity is also and its calculation using Equations $\eqref{eq:mg}$ and $\eqref{eq:polarisation}$ need to be calculated over the frequency range of interest.

### Bruggeman mixing rule

In the Maxwell-Garnett mixing formalism there is a distinction between the inclusion and the supporting medium which results in an asymmetry in the treatment of the two species in the mixture. Instead the Bruggeman mixing rule assumes that each species is polarized against the background of the effective medium and therefore the polarisation in the two components cancel each other out;
$$
\left\langle \fieldbf{P}_1 \right\rangle + \left\langle \fieldbf{P}_2 \right\rangle = 0\label{eq:EqualPolarisation}
$$
where the components are now labeled 1 and 2 rather than external and internal. The polarisation for species 1 and 2 with a number density of species represented by $n_1$  and $n_2$ can be obtained from the polarisability of the species (Equation $\eqref{eq:EqualPolarisation}$);
$$
\left\langle \fieldbf{P}_1 \right\rangle = n_1 \left\langle \tensorbs{\alpha}_1 \right\rangle\fieldbf{E}\label{eq:polarisation1}
$$
Substituting Equation $\eqref{eq:polarisation1}$6 into Equation $\eqref{eq:EqualPolarisation}$ leads to the requirement that;
$$
n_1\left\langle \tensorbs{\alpha}_1 \right\rangle + n_2\left\langle \tensorbs{\alpha}_2 \right\rangle = 0 \label{eq:BruggemanPolarisability1}
$$
Taking Equation 19 and generalizing it for species *I*, (where *I* is 1 or 2) embedded in an effective permittivity given by $\tensorbs{\varepsilon}_{br}$;
$$
\tensorbs{\alpha}_i = V_i \tensorbs{\varepsilon}_{br} \left( \tensorbs{\varepsilon}_i - \tensorbs{\varepsilon}_{br} \right)
\left( \tensorbs{\varepsilon}_{br} + \tensorbf{L} \left( \tensorbs{\varepsilon}_i - \tensorbs{\varepsilon}_{br} \right) \right)^{- 1} \label{eq:BruggemanPolarisability2}
$$
Equation $\eqref{eq:BruggemanPolarisability1}$ has to be solved for$\tensorbs{\varepsilon}_{br}$. Since the systems considered here are isotropic with random inclusions, a solution has to be found for a complex value of the Bruggeman permittivity at each frequency considered.  An issue in the use of Equation $\eqref{eq:BruggemanPolarisability}$ is that the same depolarisation matrix is being used for both species, which is clearly not always appropriate. The solution of Equation $\eqref{eq:BruggemanPolarisability1}$ can be achieved either by iteration or by casting the equation as a minimization problem. The iterative approach implemented in PDielec involves repeated application of Equation 29 until convergence [@Mackay2009].  The starting point for the iterations is taken as the Maxwell-Garnett solution for the first frequency and then the solution at the previous
frequency is used to start the iterations.
$$
\tensorbs{\varepsilon}_{br} = \frac{f_1{ {\tensorbs{\varepsilon}}}_1\left\lbrack {\tensorbf{1}} + {\tensorbf{L}}\left( {\tensorbs{\varepsilon}}_1 - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1} +
f_2\tensorbs{\varepsilon}_2\left\lbrack \tensorbf{1} + \tensorbf{L}\left( \tensorbs{\varepsilon}_{2} - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1}}
{f_{1}\left\lbrack \tensorbf{1} +{\tensorbf{L}}\left( \tensorbf{\varepsilon}_1 - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1} +
f_{2}\left\lbrack \tensorbf{1} + \tensorbf{L}\left( \tensorbs{\varepsilon}_2 - \tensorbs{\varepsilon}_{br} \right) \right\rbrack^{- 1}} \label{eq:br}
$$
Although the Bruggeman permittivity is written here as a tensor, the polarisabilities in Equation $\eqref{eq:BruggemanPolarisability1}$ have to be averaged over the random orientation of the inclusions and therefore the homogenized material is isotropic with a single complex value for the diagonal tensor.  Also, as with the Maxwell-Garnett mixing rule, since the polarisability is complex and frequency dependent, the effective permittivity is also, and its calculation using Equation $\eqref{eq:mg}$ needs to be performed over the frequency range of interest.

The choice between using the Bruggeman or Maxwell-Garnett model is often governed by the assumption that the Maxwell-Garnett model works well at low concentrations and the Bruggeman model works better at higher concentrations. Work by Karkkainen *et al*. using a finite difference method for random mixtures of non-absorbing materials indicated that the Bruggeman approximation works best when there is some clustering of the inclusions and the Maxwell Garnett model works best when there is no clustering [@Karkkainen2001].

The Bruggeman solution has been shown to be unphysical in certain circumstances [@Jamaian2010].  In particular when the real components of the permittivity have different signs and when the absolute values of the real components are much larger than those of the imaginary components. Unfortunately, it may be that these conditions will apply to modelling infrared absorption. As a result only a few of the examples discussed below will include results using the Bruggeman mixing rule; the majority will use the Maxwell-Garnett mixing rule.

### Averaged-Permittivity mixing rule

It is useful to be able to compare the effective medium theories with the absorption predicted using no shape information, that is using only the TO frequencies.
$$
\tensorbs{\varepsilon}_{TO} = f\left\langle \tensorbs{\varepsilon}_i \right\rangle + \left( 1 - f \right)\varepsilon_e\label{eq:ap}
$$
Equation $\eqref{eq:ap}$ defines an isotropic permittivity which can be used to calculate such an absorption coefficient. The angle brackets indicate an average of orientation. This mixing rule provides a useful comparison between the absorption calculated without any shape effects and that calculated including shape effects using the other mixing rules presented above. At low concentrations the peak positions of the AP mixing rule will be at the TO frequencies.

Particle Size Effects
---------------------

Meier and Wokaun [@Meier1983] outlined an approach to treating large (metal) spherical particles, where particle size is incorporating terms up to 3^rd^ order in the wave vector *k*. Using Equations $\eqref{eq:PolarisationDensity}$ and $\eqref{eq:InternalField}$ we can write;
$$
\fieldbf{P} = \left( \tensorbs{\varepsilon}_i - \varepsilon_e \tensorbf{1} \right)\left( \fieldbf{E}_e + \fieldbf{E}_d \right)
$$

$$
\fieldbf{E}_d = - \frac{\left( 1 - x^2 - i\frac{2}{3}x^3 \right)}{\varepsilon_e}{\tensorbf{L}} \fieldbf{P} \\
x = ak = \frac{2\pi a}{\lambda} 
$$

Here $a$ is the radius of thespherical particle and $x$ is the dimensionless 'size' of the particle with respect to the wavelength of the incident light. The first term relating the depolarisation field to the polarisability is the same as that used above in Equation $\eqref{eq:DepolarisationField}$. The second term is a dynamic depolarisation term and the third term is a radiation damping correction to the electrostatic solution[@Meier1983].

A slightly different, but related, approach is presented by [@Sihvola].  Starting from Equation $\eqref{eq:InternalField2}$, a size dependent term $G(x)$ is introduced as indicated by the work of Peltoniemi [@Peltoniemi1996];
$$
G\left( x \right) = G_1\left( x \right) + G_2\left( x \right) \\
G_{1}\left( x \right) = \ \frac{2}{3}\left\lbrack \left( 1 + ix\right)e^{- ix} - 1 \right\rbrack \\
G_{2}\left( x \right) = \ \left ( 1 + ix - \frac{7}{15}x^{2} - i\frac{2}{15}x^{3} \right)e^{- ix}-1
$$
The modified equation for the relationship between the internal field and the external field based on equation $\eqref{eq:InternalField2}$ becomes; 
$$
\fieldbf{E}_i = \fieldbf{E}_e - \frac{\left( 1 - G\left( x \right) \right)}{\varepsilon_e}\tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_e\tensorbf{ 1} \right)\fieldbf{E}_i
$$
This leads to a modified equation for the polarisability of spherical particles;
$$
\tensorbs{\alpha}\left( x \right) = V\varepsilon_e\left( \tensorbs{\varepsilon}_i - 
\varepsilon_e\tensorbf{1} \right)\left( \varepsilon_e\mathbf{1} + 
\left( 1 - G\left( x\right) \right)\tensorbf{L} \left( \tensorbs{\varepsilon}_i - \varepsilon_{e}\mathbf{1} \right) \right)^{- 1} \label{eq:ModifiedPolarisability}
$$
Using the modified, sized dependent polarisability all the Bruggeman and Maxwell mixing schemes can be implemented in a way that incorporates size effects.  Generally speaking the on-set of changes in the calculated absorption is an indication that size effects are important and should be treated properly.

Light Scattering by Spherical Particles using Mie Theory
--------------------------------------------------------

For particles which are comparable in size to the wavelength of light, the theory developed by Mie and described fully by van de Hulst^14^ can be used. Unfortunately this theory is only applicable to spherical, isotropic particles where the separation between the particles is large compared with the wavelength of light.

PDielec implements a form of Mie theory using the Python library PyMieScatt [@Sumlin2018a]. In order to treat systems which are anisotropic, PDielec first of all transforms the real component of the permittivity so that the tensor is diagonal. Then the full  permittivity (real and imaginary) is transformed to this basis.
$$
\tensorbf{U}^T \tensorbs{\varepsilon}^{real} \tensorbs{U} = \tensorbs{\varepsilon}^{diagonal}\\
\tensorbs{\varepsilon}^{new} = \tensorbf{U}^{T} \tensorbs{\varepsilon}^{full}\tensorbf{U} 
$$
The core of Mie theory is the calculation of efficiency factors S1(θ) and S2(θ) (page 125 in reference [@VanDeHulst1981].  These factors depend upon the size of the particle with relative to the wavelength of the light ($x$) and the refractive index of infrared active material ($m_1$) relative to that of the supporting medium ($m_2$). Both refractive indices are scalars, $m_1$ is complex as light is absorbed by the material and $m_2$ is assumed to be real.

For each of the three diagonal elements of $\tensorbs{\varepsilon}_{new}$ the relative refractive index, $m$ is calculated;
$$
m_1 = \tensorbs{\varepsilon}_{w,w}^{new} \\
m = \frac{m_{1}}{m_{2}} \\
x = \frac{2\pi a}{\lambda} \\
\lambda = \frac{\lambda_{\text{ vac}}}{m_{2}}
$$
The size parameter depends upon the refractive index of the supporting medium because the wavelength of light is smaller than its vacuum value in a medium with a refractive index greater than one. The subscript $w$ indicates one of the three directions x, y or z in the basis which diagonalises the real permittivity tensor.

Using the routines available in PyMieScatt [@Sumlin2018a] the efficiency factor S(0) (S(0) = S1(0) = S2(0)) is calculated from $m$ and $x$. The extinction (this includes absorption and scattering) along with retardation are reflected in the overall complex refractive index, page 129 in reference [@VanDeHulst1981];
$$
\widetilde{m}_w = m_2\left( 1 - iS\left( 0 \right)2\pi Nk^{- 3} \right) \\
N = \frac{f}{V_{sphere}}  \\
V_{sphere} = \frac{4}{3} \pi a^{3} \\
k = \frac{2\pi}{\lambda}
$$
Here *N* is the number density of particles and *k* is the wave-vector of the light in the surrounding medium. The subscript *w* indicates that there are three values for the effective refractive index, depending on which diagonal value of the permittivity tensor is taken. Once the effective refractive index is known in each direction the effective permittivity can be calculated using the one third rule [@Stout2007];
$$
\widetilde{m}_{eff} =  \frac{1}{3}\sum_w{\widetilde{m}}_w \\
\varepsilon_{eff} = \widetilde{m}_{eff}^2
$$
This approach to taking anisotropy into account when the embedded particles are anisotropic but randomly oriented is approximate, but has been shown to have a reasonably wide range of application [@Stout2007].

IMPLEMENTATION
==============

The above theory has been implemented in a Python 2/3 package which is available for download [@pdielec]. ^35^ The package requires SCIPY [@scipy}], NUMPY [@scipy}],
PyYAML (py-yaml) [@pyyaml], PyMieScatt [@Sumlin2018a], xlswriter [@xlsxwriter] and if visualization of the predicted spectra is required MATPLOTLIB [@scipy]. The program is run from the command line. There are several command options and these are summarized below in Table 1. At the moment the package has interfaces to five solid state QM codes, VASP [@Hafner2008c], CASTEP [@Clark2005d], CRYSTAL14 [@Dovesi2014], Abinit [@Gonze2016], Quantum Espresso [@Giannozzi2009] and Phonopy [@Togo2015].  In addition an interface is available for GULP [@Gale2003] which is a force field based solid state code. Finally an interface has been written to and 'experiment' file format which allows the preparation of a user defined file specifying the permittivities and absorption frequencies. The origin of the dataset(s) used for processing is determined by a command line switch, -program. An outline of the interfaces to these codes is given here. The package used for the calculation is described by the --program option. In addition a file name is given which contains the output which will be processed by PDielec.

**VASP** -program vasp OUTCAR
The name provided on the command line is an OUTCAR file. The OUTCAR is read by PDielec to determine the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. The VASP run can be a DFPT or numerical calculation of the response.

**CASTEP** -program castep seedname
The name provided on the command line is the seedname for the calculation. The corresponding seedname.castep file in the current directory is read and processed to determine the unit cell, atomic masses, optical permittivity and born charge tensors. The normal modes and their frequencies are determined from the seedname.phonon file. The CASTEP run needs to be a DFPT (phonon+efield) task.

**CRYSTAL** -program crystal outputfilename
The name on the command line is a file ending in .out, containing the output of a CRYSTAL14 run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory where PDielec is run from , it uses these files to calculate the Born charge tensors, frequencies and normal modes. The CRYSTAL calculation needs to be a frequency calculation (FREQCALC) with the infrared intensity (INTENS) selected. The default algorithm does not calculate the optical permittivity, so this needs to be provided on the command line. However, if the CPHF or CPKS algorithm is used for the frequency calculation, the optical permittivity is calculated and PDielec will automatically read it from the output file. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the PDielec package. Small differences in the calculated frequencies between the CRYSTAL program and PDielec have been observed. These have been found to be due to a slightly different method for symmetrising the 2^nd^ derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that PDielec should use the same symmetrisation as
CRYSTAL14.

**ABINIT** -program abinit outputfilename
The output file should come from a run containing three datasets. One to calculate the wavefunction at the optimized geometry, one to calculate the field perturbations and one to calculate the second derivatives. Examples of input files and output files are available with the distribution.

**QE** -program qe outputfilename
The output file is the dynamical matrix file, specified by "filedyn" in a run of the quantum espresso phonon package. Examples of input and output files are given in the PDielec distribution

**PHONOPY** -program phonopy vasp OUTCAR.born
Phonopy calculates the dynamical matrix through numerical differentiation. It has interfaces to several programs, although PDielec has only used the VASP interface. In principle other interfaces could be used. The second parameter for the --program directive is the PHONOPY interface that was used to calculate the forces. Typically these would be generated by performing;
`phonopy --d --dim="1 1 1"`
to calculate the displacements in a set of POSCAR-\* files. After running VASP a single point VASP calculation for each displacement. The FORCE\_SETS file can then be calculated using for example;
`phonopy --f DISP-\*/vasprun.xml`
where the DISP-\* directories are where the VASP calculation was performed. Finally a dynamical is written out using;
`phonopy --dim="1 1 1" --qpoints="0 0 0" --writedm`
To calculate the infrared spectrum PDielec needs the Born charges for the atoms in the unit cell and these can be calculated using VASP and the optimized geometry of the unit cell. The OUTCAR file from this calculation can be copied to the current directory and renamed OUTCAR.born

**GULP** -program gulp outputfilename
The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit cell. The key words; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it on the command line (see -optical and -optical\_tensor options below).

**EXPERIMENT** -program experiment\_filename
This option has been added to allow the exploration of 'toy' problems. The file contains a minimum amount of data to allow a calculation to proceed. It is assumed that the systems will be isotropic as that make the input simpler. Calculations of the LO frequencies will not work with this option. An example input file is given here, which gives results very similar to that found for the Castep, MgO example;

`species 2 # Define the atomic species, followed by the number of species`
`O 16.0 # Species name followed by its mass`
`Mg 24.3`
`lattice 2.12346 # Define the lattice, the parameter is a scaling factor`
`1 1 0           # The a vector`
`1 0 1           # The b vector`
`0 1 1           # The c vector`
`unitcell 2 # Define the unit cell contents and fractional coordinates`
`O 0.5 0.5 0.5  # Species fraca, fracb, fracc`
`Mg 0.0 0.0 0.0 # Species fraca, fracb, fracc`
`static         # Define the static permittivity`
`3.13969 0.0 0.0`
`0.0 3.13969 0.0`
`0.0 0.0 3.13969`
`frequencies 3       # Specify how many frequencies to read in`
`388.282 0.000073639 # The frequency in wavenumbers and oscillator strength`
`388.282 0.000073639`
`388.282 0.000073639`

Examples of data sets for these packages are included with the distribution. The interface to these QM and MM codes reads information about the unit cell, the calculated normal modes and the Born charge matrices; from these the permittivity is calculated over the frequency range requested. The absorption and molar absorption coefficients can be plotted along with the real and imaginary permittivities. Optionally all the information can be written to a comma separated values (csv) file for direct importing into a spreadsheet. The program is run from the command line. There are several command options and these are summarized below in {@tbl:options}. Some of the options may be repeated. The package needs a shape to be specified (sphere, needle, plate or ellipse). If no shape is specified on the command line a sphere is assumed.

| Option                     | Default | Purpose                                  | R[^foot1] |
| :------------------------- | ------- | :--------------------------------------- | :-------: |
| -program string            |         | string can be “abinit”,  “castep”, “crystal”, “gulp”, “qe”, “experiment” or “vasp” and specifies the program which generated the results to be analysed |           |
| -method string             |         | The method is given by the string and is either ‘ap’, ‘maxwell’, ‘bruggeman’ or ‘mie’ |     ✔     |
| -sphere                    |         | The inclusion is a sphere, the default if no other shape is given. |           |
| -needle h k l              |         | The inclusion is a needle whose unique directionis given by the direction [hkl]. |     ✔     |
| -plate h k l               |         | The inclusion is a plate whose  surface is defined by the Miller indices (hkl).   Note that needles and ellipsoid use  directions in crystal coordinates defined by [hkl].   For non-orthogonal lattices the normal to  the (hkl) is not necessarily the same as [hkl]. |     ✔     |
| -ellipse h k l z           |         | The inclusion is an ellipsoid, whose  unique direction is given by [hkl],  z  specifies the eccentricity of the ellipsoid. |     ✔     |
| -vf real                   | 0.1     | real specifies the volume fraction       |     ✔     |
| -mf real                   | 0.0     | real specifies a mass fraction from  which the volume fraction is calculated.   The calculation requires the density of the supporting matrix. |     ✔     |
| -matrix string             | ptfe    | The supporting matrix is defined by  the string.  Options are “ptfe”, “kbr”,  “ujol”, “air”, “vacuum”, “ldpe”, “mdpe”, “hdpe”.  If the matrix is given in this way both the  density and the permittivity of the supporting matrix are defined.  Alternatively the density and dielectric  options can be used. |           |
| -density real              | 2.2     | real defines the density of the  supporting matrix |           |
| -dielectric real           | 2.0     | real defines the dielectric of the  supporting matrix |           |
| -LO h k l                  |         | The frequencies corresponding to the  longitudinal optic modes with a k vector direction (h k l) are calculated  using Equations 10 and 11 |     ✔     |
| -LO_cart x y z             |         | As above but for Cartesian  directions   |     ✔     |
| -sigma real                | 5.0     | real specifies the damping factor, σ,  for all modes in cm^-1^, as used in Equation 10a |           |
| -mode_sigma k z            |         | The k’th mode is assigned a specific  σ (cm-1) given by z. |           |
| -vmin real                 | 0.0     | real is the starting wavenumber (cm^-1^)  for the frequency range |           |
| -vmax real                 | 300.0   | real is the final wavenumber (cm^-1^)  for the frequency range |           |
| -i real                    | 0.2     | real is the increment used to cover  the frequency range (cm^-1^) |           |
| -plot sstring              |         | Plot types  are specified by the string and they can be  ‘absorption’, ‘molar_absorption’, ‘real’ or ‘imaginary’ |     ✔     |
| -excel s                   |         | Writes the results to an excel  spread sheet with the name  specified  by the string s. |           |
| -csv s                     |         | Output is sent to a comma separated  file specified by the string s. |           |
| -csv_ext string            |         | Output is sent to 3 comma separated  files; string_command.csv, string_frequency.csv and string_spectrum.csv |           |
| -print                     |         | Additional output is provided from  the QM or MM calculation |           |
| -ignore k                  |         | Ignore the kth mode (any mode below 5cm^-1^ is ignored automatically) |     ✔     |
| -mode k                    |         | Only use the kth mode in the  calculation of the permittivity |     ✔     |
| -threshold z1 z2           |         | The modes selected for inclusion in  the absorption calculation have to have an IR intensity greater than z1 and a  frequency greater than z2. By default z1 and z2 are 1.0e-10 and 5cm^-1^  respectively. |           |
| -eckart                    |         | The translational modes are  projected out of the hessian before diagonalisation |           |
| -hessian string            |         | string may be either “crystal” or  “symm”.  In the case of “crystal” the  hessian is symmetrised using the same algorithm as Crystal14.  “symm” is the default |           |
| -optical z1 z2 z3          |         | z1,z2 and z3 define the diagonal of  the optical permittivity tensor |           |
| -optical_tensor z1 z2 ..z9 |         | z1,..9 define the full optical permittivitytensor |           |
| -masses string             |         | string can be either “program”, “average” or“isotopic”, meaning that the masses used in the calculation of the frequenciesare either taken from the QM program or are the average of the isotopeabundances or are the most abundant isotope mass. |           |
| -mass string real          |         | The atomic mass of string is set to real.  This can be used to explore the effect ofisotope substitution on the calculated frequencies |     ✔     |
| -processors int            |         | The number of processors to be used  in the calculation can be defined.  By  default all available processors are used. |           |
| -molesof string [int]      |         | By default string is “cells”.  Other options are “atoms” or  “molecules”.  If ‘molecules’ is  specified then the number of atoms in a molecules must be provided |           |
| -size real [sigma]         |         | Used to modify the polarisability (Eq.34) for spherical particles which incorporates the radius of the particle inmicrons(real).   It is also used to specify the dimension of thespherical particles for the Mie method.  It is also possible to specify a sizedistribution in which case the first number is the mean of the log distribution(in microns) and the second is its width (in log(microns)) |           |
Table: PDielec command line options {#tbl:options}
[^foot1]: This column indicates if a command line option can be used more than once

The shape options; ellipse, slab and needle, specify a unique axis \[hkl\] using the crystal axes of the unit cell. PDielec transforms these to a cartesian coordinate system using the unit cell lattice vectors. In the case of a slab morphology the unique direction is normal to the surface specified by its Miller indices (hkl). The definitions of the various depolarisation tensors are indicated in Table {@tbl:depol} below.

| Shape     | Depolarisation Tensor                    |
| --------- | ---------------------------------------- |
| Sphere    | $$\tensorbf{L} =\frac{1}{3}\left( \fieldbf{V}_1 \fieldbf{V}_1^T + \fieldbf{V}_2 \fieldbf{V}_2^T + \fieldbf{V}_3 \fieldbf{V}_3^T \right)$$ |
| Slab      | $$\tensorbf{L}=\fieldbf{V}_1 \fieldbf{V}_1^T$$ |
| Needle    | $$\tensorbf{L}=\frac{1}{2}\left( \fieldbf{V}_2\fieldbf{V}_2^T + \fieldbf{V}_3 \fieldbf{V}_3^T \right)$$ |
| Ellipsoid | $$\tensorbf{L}=a\fieldbf{V}_1 \fieldbf{V}_1^T + b\fieldbf{V}_2 \fieldbf{V}_2^T + b\fieldbf{V}_3 \fieldbf{V}_3^T$$ |
Table: Definitions used of the depolarisation tensor.  {#tbl:depol}

The three directions defined by $\fieldbf{V}_1, \fieldbf{V}_2 \text{and} \fieldbf{V}_3$ are mutually orthogonal cartesian vectors calculated from \[hkl\] for an ellipse  or needle and (hkl) for a slab. In the case of a slab, needle or ellipsoid, ${\overline{V}}_{1}$ defines the unique direction and the other vectors are orthogonal to it. For the case of an ellipsoid, the parameters *a* and *b* in Table {@tbl:depol} depend on the ratio, $z$, of the length of unique axis length over the length of an axis perpendicular to it [@Sihvola].

For z \> 1 the ellipsoid is prolate;
$e = \sqrt{1 - z^{- 2}},\ a = \frac{\left( 1 - e^{2} \right)}{2e^{3}}\left( \log\frac{1 + e}{1 - e} - 2e \right),\ b = \frac{1}{2}\left( 1 - a \right)$

For z \< 1 the ellipsoid is oblate
$e = \sqrt{z^{- 2} - 1},\ a = \frac{\left( 1 + e^{2} \right)}{e^{3}}\left( e - \operatorname{}e \right),\ b = \frac{1}{2}\left( 1 - a \right)$

From an experimental point of view it is often convenient to use a mass fraction rather than a volume fraction to indicate the amount of dielectrically active material present.  PDielec allows mass fractions to be specified instead of a volume fraction, but this requires that the density of the supporting matrix is known. For convenience the package has a small database of the common supporting materials shown in Table {@tbl:matrix} below.  These can be specified through the -matrix option. In the case that the properties of the support material are different the properties can be defined instead with the -dielectric and -density options. 

| Name   | Density | Permittivity | Description                 |
| ------ | ------- | ------------ | --------------------------- |
| ptfe   | 2.2     | 2.0          | polytetrafluorethylene      |
| air    | 0.0     | 1.0          | air                         |
| vacuum | 0.0     | 1.0          | vacuum                      |
| kbr    | 2.75    | 2.25         | potassium bromide           |
| nujol  | 0.838   | 2.155        | Nujol                       |
| hdpe   | 0.955   | 2.25         | high density polyethylene   |
| mdpe   | 0.933   | 2.25         | medium density polyethylene |
| ldpe   | 0.925   | 2.25         | low density polyethylene    |
Table: Physical properties of matrix materials in PDielec {#tbl:matrix}

The optical permittivity is normally calculated by the QM or MM program concerned. However, as this property reflects the electronic contribution to the permittivity at zero frequency, unless there is some treatment of electrons by the shell model, then in MM calculations the optical permittivity needs to be defined through the command line options -optical or -optical_tensor.
Unlike the other methods 'mie' method cannot work with particles of zero radius. All methods therefore use a default size of 10^-12^ μm. The Mie approach is only valid for dilute dispersions and for spherical particles. However, if other shapes are specified the Mie method will still be called and the results will be applicable to spheres of the specified size.

Parallelization, threads and performance
----------------------------------------

To improve the performance of the program python parallelization has been used to parallelize over the frequencies, shapes and methods. By default this parallelization spawns additional Python executables, depending on the number of cores available.

In addition to this form of parallelization the NUMPY library can use multi-threaded BLAS. NUMPY can be compiled with several different BLAS libraries, including; MKL (from Intel), OPENBLAS and ATLAS, or indeed no optimized BLAS library is necessary. To explore the BLAS version compiled with your version of NUMPY please use the test_numpy_1, test\_numpy\_2 and test\_numpy\_3 scripts included in the Python/ subdirectory. If you are using MKL or OPENBLAS, the number of threads being used needs to be set with the MKL\_NUM\_THREADS or OPENBLAS\_NUM\_THREADS environment variable (sometimes OMP\_NUM\_THREADS is also used). Use the test routines to determine the optimum for your system. Experience indicates that no performance benefit is obtained with more than two threads. 

In the case of Phonopy the dynamical matrix is read from a yaml file. This has been found to be very slow unless the C parser is used. If the C parser is not available a warning is issued and the program reverts back to the Python parser.

Finally the use of non-standard BLAS libraries seems to cause problems with the affinity settings for the multiprocessing. It has been noticed that the parallel processes can all end up executing on the same processor. In order to prevent this, before executing the pdielec and preader scripts it may be necessary to include;

`export OPENBLAS\_MAIN\_FREE=1`

For some reason this also works if the MKL library is being used.

Example command line uses of PDielec
------------------------------------

`pdielec -program vasp OUTCAR -method ap -method maxwell \`
`              -sphere -plate 0 0 1 -needle 0 0 1 -LO 0 0 1`

This performs a calculation using the Averaged-Permittivity and Maxwell-Garnett mixing rules for spherical particles, plate-like particles with a surface (001) and needle-like particles with a unique
direction lying along the \[001\] direction. The supporting matrix is taken to be PTFE and the default volume fraction (10%) is used. The results of a VASP calculation are stored in the current directory.  There is no absorption output from this command as neither the -plot nor the -csv options were specified. The output includes the calculation of the LO modes along the (001) direction.

`pdielec -program castep phonon -vmin 300 -vmax 800 \`
`-sphere -dielectric 3 -vf 0.1 -vf 0.2 -sigma 10 -csv mgo.csv`
This performs a calculation for spherical particles varying the frequency from 300 to 800 cm^‑1^, the permittivity of the supporting media is 3, two volume fractions are considered and a damping factor of 10 cm^-1^ is used. The results of a CASTEP calculation with the seed-name "phonon" are analysed and the results stored in mgo.csv for further analysis using a spreadsheet. In this example a Maxwell-Garnett mixing rule is used by default. If visual inspection of the results is required then

`pdielec -program castep phonon -vmin 300 -vmax 800 \`
`-sphere -dielectric 3 -vf 0.1 -vf 0.2\`
`-sigma 10 -csv mgo.csv -plot molar_absorption`

will perform the same calculation but a graph showing the molar absorption coefficients will be displayed.

`pdielec -program gulp calcite.gout -matrix hdpe \`
`-method ap -method maxwell -sphere -plate -1 -1 -2 \`
`-vmax 2000 -mf 0.1 calcite.gout -csv calcite.csv`

This command performs a calculation of the absorption spectrum resulting from a GULP calculation. The supporting matrix density and permittivity are those of high density polyethylene, the frequency range is 0 to 2000 cm^-1^, the mass fraction considered is 10%, the mixing rules used are Averaged-Permittivity and Maxwell-Garnett. Spheres and plates with the ($\overline{1}\overline{1}\overline{2})$ surface are considered.

`pdielec -program vasp OUTCAR -method mie -sphere -mf 0.1 -size 0.1 -size 1.0 -size 1.0 -csv results.csv`

This command performs a calculation of the absorption spectrum resulting from a VASP calculation using the Mie method for sphere with several particles sizes.

`pdielec -sphere -mf 0.1 -program experiment experiment.expt -size 0.1 -size 1.0 -size 1.0 -excel results.xlsx`

This command performs a calculation of the absorption spectrum resulting from the data stored in the experiment.expt file. Maxwell-Garnett calculations are performed with 3 different sized spheres and the results stored in a Excel file.

Contents of the csv output file
-------------------------------

If a csv output file is requested the file will contain the command used to perform the calculation. A brief summary is given of each active infrared mode; including the mode number, frequency, intensity, integrated molar absorption coefficient, its peak height (calculated from the intensity and damping factor) and the damping parameter used in the calculation. Following this is a table with a column for frequency followed by columns containing the real and imaginary permittivities, the absorption and molar absorption coefficients at each frequency.

Parallelization
---------------

PDielec has been written to make use of multiprocessor computers. On a 4
processor machine the speed-up is nearly linear up to 4 processors, as
can be seen in the Figure below.

![](media/image2.png){width="3.9556988188976376in"
height="2.3819444444444446in"}

Figure : Speed up on a 4 processor machine

EXAMPLES
========

Several examples are given to illustrate applications of the package.
The calculations used to provide the data for the permittivities are
sufficiently accurate to illustrate aspects of the theory. The examples
are chosen to show the package being used with the QM packages CASTEP
and VASP and with the MM package GULP.

MgO using CASTEP
----------------

Magnesium oxide is an isotropic medium, the initial unit cell and the
space group symmetry ($Fm\overline{3}m$) were taken from the Inorganic
Crystal Structure Database (ICSD)^39^ reference number ICSD-52026.^40^
The primitive cell was optimized using CASTEP. Norm-conserving
pseudo-potentials were used to represent the core electrons of magnesium
and oxygen. An energy cutoff of 1000 eV was used with the PBE^41^
density functional and a k-point spacing for the Monkhorst-Pack grid of
0.04 Å^-1^. The primitive cell was optimized and a Density Functional
Perturbation Theory (DFPT) calculation of the phonon spectrum at the
gamma point was performed. The optimised lattice parameter was found to
be 2.1234 Å, compared with the experimental value of 2.107 Å. Only 3
degenerate modes contribute to the permittivity. A summary of the
results is presented in Table 4.

Table 4: Calculated Properties of MgO

  **Property**                      **Values(s)**         **Units**   
--------------------------------- --------------------- ----------- ----------------------
  Unit cell dimensions^a^           2.123 (2.107)         Å           
  Space group                       $$Fm\overline{3}m$$               
  Optical permittivity              3.14                              
  Static permittivity               10.0                              
  Phonon frequency (intensity)^b^   TO                    LO (001)    cm^-1^ (D/Å)^2^/amu)
                                    T 388.3 (9.29)        693.7       

^a^The experimental value is given in brackets\
^b^The intensities are given in brackets, T indicates a triply
degenerate mode

Because MgO is isotropic with only a single frequency contributing to
the permittivity, it makes a useful example application to illustrate
several features of PDielec. The real and imaginary frequency dependent
permittivities are shown in Figure 3, where a damping factor (σ) of 10
cm^-1^ has been used. In the Figure the real permittivity at zero
frequency corresponds to the static permittivity in Table 4, and at
frequencies above the absorption at 388 cm^-1^ the permittivity tends to
the optical permittivity as the frequency increases. The real
permittivity has zero values at 388.3 and 693.7 cm^-1^ which are the TO
and LO frequencies respectively.

![Figure2\_MgO\_permittivity.tiff](media/image3.tiff){width="5.195in"
height="4.07in"}

Figure : Permittivity of MgO

![Figure3\_Mgo\_Real\_imaginary.tiff](media/image4.tiff){width="5.18in"
height="4.07in"}

Figure : Real and Imaginary permittivities of a 10% volume fraction of
MgO spheres in PTFE, calculated using the Maxwell-Garnett method

Using the Maxwell-Garnett mixing rule, Figure 4 shows the calculated
permittivities of a 10% volume fraction of MgO spheres in a supporting
medium with a frequency independent permittivity of 2.0, which would be
typical of a material such as PTFE. Due to the dilution effect the real
component has shifted to a base line value close to 2, and the
absorption, as indicated by the maximum in the imaginary component has
shifted by about 150 cm^-1^ to 550 cm^-1^.

The effect of volume fraction on the predicted molar absorption
coefficient, using the Maxwell-Garnett mixing rule, is shown in Figure
5. The lowest volume fraction of MgO gives the largest shift of the
  absorption peak to high frequency. As the volume fraction increases the
  mixing rule predicts a broadening of the absorption, whilst the peak in
  the molar absorption coefficient moves to lower frequency. At the
  highest loading (f=0.9) the maximum absorption occurs quite close to the
  TO frequency. The Maxwell-Garnett mixing rule is regarded as being
  appropriate for low volume fractions and so should not be used for
  interpreting results in which higher volume fractions of absorbing media
  have been used.^25^

![Figure4\_Mgo\_MG\_volume\_fraction.tiff](media/image5.tiff){width="5.236666666666666in"
height="4.07in"}

Figure : Effect of volume fraction on the Maxwell-Garnett molar
absorption coefficient of MgO spheres in PTFE

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

![Figure5\_Mgo\_Bruggeman\_volume\_fraction.tiff](media/image6.tiff){width="5.236666666666666in"
height="4.07in"}

Figure : Effect of volume fraction on the Bruggeman molar absorption
coefficient of MgO spheres in PTFE

![Figure6\_MgO\_varying\_permittivity.tiff](media/image7.tiff){width="5.236666666666666in"
height="4.07in"}

Figure : The Maxwell-Garnett molar absorption coefficients of spherical
MgO particles, 1% volume fraction, embedded in media of varying
permittivities

Figure 7 shows the effect of varying the permittivity of the supporting
medium. The calculations were performed on spherical MgO particles with
a 1% volume fraction. The lowest permittivity is that of a vacuum (or
air) and shows the highest shift of the absorption maximum to higher
frequencies. Increasing the permittivity lowers the shift until it
becomes quite small. A similar effect is seen for the Bruggeman mixing
model. However, the absorption resulting for particles in a low
dielectric medium is considerable broader than that seen in the
Maxwell-Garnet case. This broadening reduces as the permittivity of the
medium increases (see Figure 8).

![Figure7\_MgO\_varying\_permittivity\_bruggeman.tiff](media/image8.tiff){width="5.236666666666666in"
height="4.07in"}

Figure 8: The Bruggeman molar absorption coefficients of spherical MgO
particles, 1% volume fraction, embedded in media of varying
permittivities

ZnO using VASP
--------------

Zinc oxide crystallizes in space group P6~3~mc (wurtzite). All
calculations were performed by VASP^1^ using projector augmented-wave
PAW^42^ pseudo-potentials, the PBE^41^ density functional, an energy
cutoff of 600 eV and a k-point resolution of approximately 0.1 Å^-1^.
The initial unit cell was taken from the ICSD^39^ with code
ICSD-26170.^43^ The unit cell and atom positions were optimized using
VASP and the permittivity was calculated using DFPT and the results
reported in Table 5. Only two of the bands showed any significant
intensity, a doubly degenerate band (E) with a TO frequency of 372.1
cm^-1^ and a non-degenerate band (A) with a TO frequency of 350.0
cm^-1^. The LO frequency of the non-degenerate band is shifted to 502.0
cm^-1^ for a wave-vector with direction (001), whilst the degenerate
modes are unaffected. In the case of the (010) direction the LO
frequency of one of the E modes is shifted to 511.2 cm^-1^. It is known
that ZnO can crystallize with a plate morphology^44^ with the (001)
surface dominant. Calculations of the molar absorption were performed
for a sphere, plate and needle like shapes with the unique directions of
the plate and the needle being normal to the (001) surface. A volume
fraction of 1% was chosen for these calculations and the predicted molar
absorption coefficients for the Maxwell-Garnett mixing rule is shown in
Figure 9.

Table 5: Calculated properties of ZnO

-------------------------------------------------------------------------------------------------------------------
  Property                          Values(s)                            Units                 
--------------------------------- ------------------------------------ ---------- ---------- ----------------------
  Unit cell dimensions^a^           a,b = 3.295(3.25) c = 5.285(5.207)   Å                     

  Space group                       P6~3~mc                                                    

  Optical permittivity^b^           5.09, 5.09, 6.0                                            

  Static permittivity^b^            10.83, 10.83, 11.67                                        

  Phonon frequency (intensity)^c^   TO                                   LO (001)   LO (010)   cm^-1^ (D/Å)^2^/amu)

                                    A 350.0 (17.1)\                      502.0      511.2      
                                    E 372.1 (16.4)                                             
-------------------------------------------------------------------------------------------------------------------

^a^The experimental values are given in brackets\
^b^Only the diagonal components are given\
^c^The intensities are given in brackets, E and A indicate a doubly and
non- degenerate mode respectively

![Figure8\_ZnO.tiff](media/image9.tiff){width="5.375in"
height="4.07in"}\
Figure : The effect of shape on the Maxwell-Garnett molar absorption
coefficient of 1% volume fraction ZnO in PTFE

For the Maxwell-Garnett mixing rule the sphere morphology results in the
two absorption peaks shifting from their TO positions to higher
wavenumber by about 80 cm^-1^. The plate morphology results in one of
the peaks moving to higher wavenumber by about 130 cm^-1^, whilst the
other remains at the TO position. The Maxwell-Garnett results are in
close accord with some experimental results by Yamamoto et al.^45^ who
measured the infrared spectrum of ZnO smoke particles and observed peaks
in the absorption at 380, 530 and 550 cm^-1^. Previous work^46,47^ has
also used effective medium theory to explain the observed spectrum.

Calcite using GULP
------------------

Calcite is the most stable polymorph of calcium carbonate and the
crystal structure belongs to the $R\overline{3}c$ space group. The force
field and atomic structures used here are described in detail in work by
Fisler *et al*.^48^ Briefly, the oxygen ions are described using a
core-shell model.^49^ The carbon - oxygen potential of the carbonate is
taken to be a Morse potential and an additional 3 atom potential is used
to maintain the O-C-O angle at 120^O^. The van der Waals interactions
between non bonded atoms are taken to be Buckingham potentials and the
charges on the calcium, carbon and oxygen ions are +2, +1.3435 and
-1.1145 respectively. The shell charge of the oxygen ion is -2.133 and
the spring constant for the core-shell interaction is 52.74 eV/Å^2^.

The unit cell was optimized using the primitive unit cell and the full
space group symmetry. The calculation of the phonon spectrum was
performed without symmetry but still using the primitive cell of the
lattice. A summary of the calculated properties is given in Table 6.

Table 6: Calculated properties of calcite

+---------------------------------+-----------------------+----------------------+
| Property                        | Values(s)             | Units                |
+=================================+=======================+======================+
| Primitive cell dimensions^a^    | a,b,c = 6.376 (6.375) | Å                    |
|                                 |                       |                      |
|                                 | α,β,γ = 46.0 (46.1)   | degrees              |
+---------------------------------+-----------------------+----------------------+
| Space group                     | $$R\overline{3}c$$    |                      |
+---------------------------------+-----------------------+----------------------+
| Optical permittivity^b^         | 1.91, 1.91, 2.0       |                      |
+---------------------------------+-----------------------+----------------------+
| Static permittivity^b^          | 6.7, 6.7, 7.1         |                      |
+---------------------------------+-----------------------+----------------------+
| Phonon frequency (intensity)^c^ | TO                    | cm^-1^ (D/Å)^2^/amu) |
+---------------------------------+-----------------------+----------------------+
|                                 | E 114.8 (2.39)        |                      |
|                                 |                       |                      |
|                                 | A 127.4 (3.36)        |                      |
|                                 |                       |                      |
|                                 | A 249.3 (1.23)        |                      |
|                                 |                       |                      |
|                                 | E 320.7 (5.82)        |                      |
|                                 |                       |                      |
|                                 | A 338.1 (4.14)        |                      |
|                                 |                       |                      |
|                                 | E 620.1 (3.38)        |                      |
|                                 |                       |                      |
|                                 | A 732.0 (26.89)       |                      |
|                                 |                       |                      |
|                                 | E 1463.6 (16.97)      |                      |
+---------------------------------+-----------------------+----------------------+

^a^The experimental values taken from reference^48^ are given in
brackets\
^b^Only the diagonal components are given\
^c^The intensities are given in brackets, A and E indicate a non- and
doubly- degenerate mode respectively.

Figure 10 shows the results of analysis of the results using PDielec.
The damping parameter used in the calculation was a value of 5 cm^-1^. A
10% volume fraction was used with sphere and plate morphologies for the
particles. The unique axis for the plate was taken to be the normal to
the (211) surfaces in the primitive cell axes (or the {104} surfaces in
the standard unit cell). Such surfaces define the rhombohedral faces
commonly seen in calcite crystals.^50^ Figure 10 shows that the doubly
degenerate TO absorption peak at 620 cm^-1^ is not significantly
affected by spherical particles and there is a small shift to higher
frequencies in the case of plate-like particles. The non-degenerate TO
transition at 732 cm^-1^, which corresponds to motion of the carbon atom
of the carbonate along the unique direction of the slab, shows a shift
to 786 and 819 cm^-1^ for the sphere and plate respectively. The doubly
degenerate peak at 1463 cm^-1^ is shifted to 1480 cm^-1^ by spherical
particles and is split by plate-like particles with one component which
shifts to 1491 cm^-1^ .

![Figure9\_calcite.tiff](media/image10.tiff){width="5.278333333333333in"
height="3.986666666666667in"}

Figure : Calculated Maxwell-Garnett absorption spectrum of 10% volume
fraction of calcite in PTFE

Fluoroapatite using VASP
------------------------

The line shapes of the infrared absorption of apatite and fluoroapatite
were examined extensively by Balan *et al*.^22^ Their calculations
included the effect of crystallite habit on the spectrum and the results
reported here are similar to their conclusions. The method used by Balan
*et al*. is an infinitely dilute Maxwell-Garnett model, so the only
difference between the methods used by them and those reported here
using PDielec are the incorporation of the volume fraction into the
theory and the use of an ellipsoidal shape for comparison with the other
shapes.

All calculations were performed by VASP^1^ using projector
augmented-wave PAW^42^ pseudo-potentials, the PBE^41^ density
functional, an energy cutoff of 600 eV and a k-point resolution of
approximately 0.1 Å^-1^. Table 7 summarises the results of the
calculations. Only the 3 highest frequency bands are reported and
discussed. The TO intensity of the highest frequency band at 1038 cm^-1^
is low and will not be discussed further. The Bravais Friedel Donnay
Harker (BFDH)^51^ crystal habit of the optimized crystal is shown in
Figure 11. The habit was calculated using the Mercury software
package.^52^ The BFDH crystal habit is often used to give an idea of the
likely important faces of a crystal. It uses only the crystal lattice
and space group to determine the crystal morphology. Figure 11 shows
that the {100} surfaces form a tube which are capped by the {011}
surfaces. The effect of different particle shapes on the predicted
spectrum is shown in Figure 12. The calculations of the spectra were
performed with a damping parameter (σ) of 2 cm^-1^. The ellipsoid was
chosen to have an aspect ratio, a/b, of 2 and a principle axis along
\[001\], which was compatible with the morphology predicted by the BDFH
method. The two TO absorption frequencies at 981 and 986 cm^-1^ have A
and E symmetry respectively. Spherical crystallites result in three
absorption peaks at around 1000, 1010 and 1015 cm^-1^. Needle shaped
crystallites leave the A symmetry TO absorption peak at 981 cm^-1^
unaffected, but shift and split the E symmetry TO peak to 1020 and 1046
cm^-1^. A plate morphology with (100) surfaces results in the A and one
component of the E TO absorption peak remaining at the TO frequencies,
with the other component of the E shifting 85 cm^-1^ to 1075 cm^-1^. The
ellipsoidal morphology show three shifted peaks at 1000, 1018 and 1045
cm^-1^. These results are consistent with those of Balan *et al*.^22^,
who gave detailed results for hydroxyapatite.

Table 7: Calculated properties of fluroapatite

+-----------------------+-----------------------+-----------------------+
| Property              | Values(s)             | Units                 |
+=======================+=======================+=======================+
| Primitive cell        | a,b = 9.447, c =      | Å                     |
| dimensions^a^         | 6.926 (9.417, 6.875)  |                       |
+-----------------------+-----------------------+-----------------------+
| Space group           | $$P6_{3}m$$           |                       |
+-----------------------+-----------------------+-----------------------+
| Optical               | 2.891, 2.891, 2.894   |                       |
| permittivity^b^       |                       |                       |
+-----------------------+-----------------------+-----------------------+
| Static                | 12.081, 12.081, 8.841 |                       |
| permittivity^b^       |                       |                       |
+-----------------------+-----------------------+-----------------------+
| Phonon frequency      | TO                    | cm^-1^ (D/Å)^2^/amu)  |
| (intensity)^c^        |                       |                       |
+-----------------------+-----------------------+-----------------------+
|                       | A 981.8 (112.6)       |                       |
|                       |                       |                       |
|                       | E 986.3 (101.0)       |                       |
|                       |                       |                       |
|                       | E 1038.1 (7.92)       |                       |
+-----------------------+-----------------------+-----------------------+

^a^The experimental values taken from Hughes et al.^53^ are given in
brackets\
^b^Only the diagonal components are given\
^c^The intensities are given in brackets, E and A indicate doubly and
non- degenerate modes respectively

![Figure10\_Fluroapatite\_morphology.tiff](media/image11.tiff){width="6.268055555555556in"
height="3.2020833333333334in"}

Figure : BDFH Morphology of fluoroapatite

![Figure11\_Fluoroapatite\_absorption.tiff](media/image12.tiff){width="5.361666666666666in"
height="3.986666666666667in"}

Figure : Calculated Maxwell-Garnett absorption spectra of 10%
fluoroapatite in PTFE

L-aspartic Acid using CASTEP
----------------------------

L-aspartic acid is a zwitterion in the solid state and so the shape of
the particles used in the measurement of IR and THz spectra maybe
important. The starting geometry for optimization of the unit cell and
molecular structure of L-aspartic acid was taken from Derissen et
al.^54^ The PBE^41^ functional was used with a plane wave energy cutoff
of 1000 eV and norm conserving pseudo-potentials. A dispersion
correction using the Tkatchenko-Scheffler scheme^55^ available in CASTEP
was applied for both the geometry optimisation and the calculation of
the phonon spectrum at the gamma point, with a value S~6~ scaling
factor^55^ of 1.0. A summary of the results of the calculations is shown
in Table 8.

Table 8: Calculated properties of L-aspartic Acid

-----------------------------------------------------------------------------
  Property                          Values(s)            Units
--------------------------------- -------------------- ----------------------
  Unit cell dimensions^a^           a = 7.597 (7.617)\   Å
                                    b = 7.028 (6.982)\   
                                    c = 5.113 (5.142)\   
                                    β=98.77 (99.84)      

  Space group                       $$P2_{1}$$           

  Optical permittivity^b^           2.68, 2.20, 2.56     

  Static permittivity^b^            4.58, 3.65, 3.65     

  Phonon frequency (intensity)^c^   TO                   cm^-1^ (D/Å)^2^/amu)

                                    84.5 (0.120)\        
                                    104.7 (0.202)\       
                                    106.0 (0.243)\       
                                    115.3 (0.474)\       
                                    137.3 (0.617)\       
                                    1290.0 (55.0)\       
                                    2945.9 (102.8)\      
                                    2947.3 (48.2)\       
                                    3053.7 (44.1)        
-----------------------------------------------------------------------------

^a^The experimental values are taken from Derissen *et al*.^54^ are
given in brackets\
^b^Only the diagonal components are given\
^c^Only selected transitions are tabulated. The intensities are given in
brackets.

The THz spectrum of L-aspartic acid has been reported by Juliano and
Korter^56^ in the frequency range 0-90 cm^-1^. The infrared spectrum has
been reported and assigned by Lopez *et al*.^57^ Figure 13 shows the
calculated absorption spectra for L-aspartic acid for three frequency
ranges. The calculation of the spectra used the Maxwell-Garnett mixing
rule with a 10% volume fraction of L-aspartic acid in PTFE and for
comparison the TO mixing rule. A damping factor of 2 cm^-1^ was used.
Spherical and a variety of plate-like inclusions were used to illustrate
their effect on the absorption spectra. Figure 13a shows the frequency
range from 60-130 cm^-1^ which is that covered by THz spectroscopy. The
shifts observed for the different particle morphologies are not large,
but the change in intensities is significant. The molecular motions
associated with phonons at these frequencies tend to be whole molecule
motion involving rotation. Figure 13b shows the frequency range from
1260-1340 cm^-1^. In this frequency range bending of the carboxylate
anion contributes to the spectrum significantly. The three different
plate morphologies show different and significant shifts in the TO
absorption peak at 1290 cm^-1^. The spherical morphology shows a shift
of around 25 cm^-1^ to higher wavenumber. Figure 13c shows the spectra
in the frequency range 2900-3100 cm^-1^, which corresponds to the motion
of O-H (below 2980 cm^-1^) and N-H (above 2980 cm^-1^) stretching. The
effect of the different possible crystal morphologies is large with
shifts to higher frequency of up to 50 cm^-1^. The spectra below 3000
cm^-1^ arises from two TO absorptions at 2946 and 2947 cm^-1^. Because
the motions associated with each mode interact differently with the
internal field within each crystal they give rise to different shifts
producing more complex spectra.

  ![Figure12a\_aspartic.tiff](media/image13.tiff){width="2.9880686789151354in" height="2.3622047244094486in"}   ![Figure12b\_aspartic.tiff](media/image14.tiff){width="3.1289720034995625in" height="2.3622047244094486in"}   ![Figure12c\_aspartic.tiff](media/image15.tiff){width="3.1289720034995625in" height="2.3622047244094486in"}
------------------------------------------------------------------------------------------------------------- ------------------------------------------------------------------------------------------------------------- -------------------------------------------------------------------------------------------------------------
  a\) Frequency range 60-130 cm^-1^                                                                             b\) Frequency range 1260-1340 cm^-1^                                                                          c\) Frequency range 2900-3100 cm^-1^

Figure 12: Calculated Maxwell-Garnett absorption spectra of 10% volume
fraction of L-aspartic acid in PTFE

MgO Example using Mie Scattering
--------------------------------

Figure 13 compares a Mie scattering calculation with the results from
Maxwell-Garnett effective medium theory. The same data set was used for
the CASTEP, MgO example. A volume fraction of 1% was used with a small
sphere radius (0.1 μ) and a broadening of 5 cm^-1^ embedded in a matrix
of PTFE. A power expansion in the size parameter of the Mie expressions
indicates that for small sizes of particles, the Mie and the
Maxwell-Garnett methods should be the same. This is verified in Figure
13.

![](media/image16.png){width="4.84375in" height="3.652083333333333in"}

Figure : Comparison of Mie and Maxwell methods. 1% volume fraction of
MgO in PTFE, sphere radius of 0.1 μ and a broadening of 5 cm^-1^

To better understand what makes particles large or small Table 9 shows
the dimensionless size parameter, *x*, as a function of wavenumber and
of sphere radius (equations 37) . It has been assumed that the
supporting medium is PTFE. Since the power expansion of the size
parameters leads to terms which are quadratic in x, it should be
expected that when x is less than about 0.1, the particles can be
considered small. It can be seen that particles less than 0.01μ are
small over the range of frequencies considered. But while 0.1μ particles
are small in the THz regime and low frequency infrared they should not
be considered small over the more extended infrared frequencies. 1μ
particles should be considered large for both the THz and the extended
infrared.

Figure 14 shows how the Mie predictions change as the particle radius
changes from 0.2 to 1.6μ. As the particle size increases the peak above
500 cm^-1^ splits into two. One broader peak which moves to lower
frequency as the particle size increases and the peak at about 550
cm^-1^ which looses intensity as the particle size increases. There is
also the onset of absorption at 388 cm^-1^which corresponds to the bulk
TO modes.

Table 9: Variation of size parameter with wavenumber and radius of
sphere

-------------------------------------------------------------------------
  Wavenumber\   Radius of sphere (μ)                              
  (cm^-1^)                                                        
------------- ---------------------- -------- -------- -------- ---------
                0.001                  0.01     0.1      1        10

  100           0.0009                 0.0089   0.0888   0.8884   8.8841

  200           0.0018                 0.0178   0.1777   1.7768   17.7682

  300           0.0027                 0.0267   0.2665   2.6652   26.6523

  400           0.0036                 0.0355   0.3554   3.5536   35.5364

  500           0.0044                 0.0444   0.4442   4.4420   44.4204

  600           0.0053                 0.0533   0.5330   5.3305   53.3045

  700           0.0062                 0.0622   0.6219   6.2189   62.1886

  800           0.0071                 0.0711   0.7107   7.1073   71.0727

  900           0.0080                 0.0800   0.7996   7.9957   79.9568

  1000          0.0089                 0.0888   0.8884   8.8841   88.8409
  -------------------------------------------------------------------------

![](media/image17.png){width="4.84375in" height="3.652083333333333in"}

Figure : Variation in absorption calculated by the Mie method for
different radii of spheres. 1% volume fraction of MgO in PTFE and a
broadening of 5 cm^-1^

Figure 15 shows the effect of increasing the particle size further. More
structure appears in the absorption, with increasing absorption around
the bulk TO frequency. Above 4.0μ there is more low frequency structure
appearing, below 300cm^-1^.

![](media/image18.png){width="4.809027777777778in" height="3.652083333333333in"}
================================================================================

Figure : Variation in absorption calculated by the Mie method for
different radii of spheres. 1% volume fraction of MgO in PTFE and a
broadening of 5 cm^-1^

ZnO Example using Mie Scattering
--------------------------------

ZnO is an anisotropic material, so the treatment described here using
Mie scattering is an approximation. However, the permittivity constant
tensor is diagonal due to the space group symmetry of the crystal.
Figure 16 compares the predicted absorption using Mie and
Maxwell-Garnett for different volume fractions of 0.1μ ZnO spheres
embedded in PTFE. A line broadening of 5 cm^-1^ was assumed.

Figures 16 and 17 compare the capabilities of the Maxwell-Garnett and
Mie methods for describing volume fraction effects. Figure 16 shows that
Maxwell-Garnett predicts a lowering of intensity and frequency of the
high frequency peak as the volume fraction is increased. Figure 17 shows
that Mie theory shows no effect of the change in volume fraction. This
is to be expected as the theory assumes that each sphere is isolated and
not affecting the other spheres around it. It should be pointed out that
the figures are plotting molar absorption coefficients. The actual
absorption would increase with volume fraction of ZnO (see equation 2).

+-----------------------------------+-----------------------------------+
| ![](media/image19.png){width="2.9 | ![](media/image20.png){width="2.9 |
| 3125in"                           | 3125in"                           |
| height="2.19375in"}               | height="2.19375in"}               |
|                                   |                                   |
| Figure : ZnO spheres in PTFE      | Figure : ZnO spheres in PTFE      |
| using Maxwell-Garnett             | using Mie                         |
+-----------------------------------+-----------------------------------+

Figure 18 shows that the variation of the Mie scattering with sphere
radius follows a similar pattern to that observed in MgO, though
slightly more complex. The initial peaks at about 440 and 460 cm-1
broaden and shift to lower frequencies as the particle size increases.
Bulk bands around 350 and 372 cm^-1^ can be seen which grown in
intensity and shift to lower frequencies as the particle size increases.

![](media/image21.png){width="4.886805555555555in" height="3.652083333333333in"}
================================================================================

Figure : Mie scattering of 10% volume fraction ZnO spheres in PTFE,
using a line broadening factor of 5 cm^-1^.

CONCLUSIONS
===========

The PDielec package has been described and examples given as to its
application in calculating the infrared absorption spectrum of a
dielectric material embedded in supporting matrix. The shape of the
crystallites can be taken into account by describing them as spheres,
plates, needles or ellipsoids. The package can calculate the dielectric
response of the effective medium as well as the infrared absorption as a
function of frequency. Several of the examples cover dielectric
materials which have been well studied, both experimentally and
theoretically and the results are in agreement with the previous work.
The package is written in Python and can be extended relatively
straightforwardly to interface with other packages. The results show the
sensitivity of the absorption spectrum to the particle morphology and
illustrate the complexity of interpreting IR and THz absorption spectra.

The PDielec package along with some example test cases for each QM or MM
package supported is available on GitHub.^35^ The data used to create
the figures and tables are openly available from the Leeds University
data repository.^58^

ACKNOWLEDGEMENTS
================

The authors would like to express their thanks to Professor Sihvola for
helpful comments on the manuscript. We would like to thank Lorenzo
Maschio who helped with the understanding of the symmetrisation of the
hessian needed for the Crystal14 interface.

ADB would also like to acknowledge financial support from both the EPSRC
(EP/I026657/1), and the DTRA (US) (HDTRA1-14-C-0013).

ACKNOWLEDGEMENTS
================

Citations of the use of this package should use reference.^59^

REFERENCES
==========

\[1\] J. Hafner, J. Comput. Chem., 2008, DOI:10.1002/jcc.21057.

\[2\] S. J. Clark, M. D. Segall, C. J. Pickard, P. J. Hasnip, M. I. J.
Probert, K. Refson, M. C. Payne, Zeitschrift für Krist., 2005,
DOI:10.1524/zkri.220.5.567.65075.

\[3\] R. Dovesi, R. Orlando, A. Erba, C. M. Zicovich-Wilson, B.
Civalleri, S. Casassa, L. Maschio, M. Ferrabone, M. De La Pierre, P.
D'Arco, Y. Noel, M. Causa, M. Rerat, B. Kirtman, Int. J. Quantum Chem.,
2014, DOI:10.1002/qua.24658.

\[4\] X.Gonze, F.Jollet, F.Abreu Araujo, D.Adams, B.Amadon,
T.Applencourt, C.Audouze, J.-M.Beuken, J.Bieder, A.Bokhanchuk,
E.Bousquet, F.Bruneval, D.Caliste, M.Cote, F.Dahm, F.Da Pieve,
M.Delaveau, M.Di Gennaro, B.Dorado, C.Espejo, G.Geneste, L.Genovese,
A.Gerossier, M.Giantomassi, Y.Gillet, D.R.Hamann, L.He, G.Jomard,
J.Laflamme Janssen, S.Le Roux, A.Levitt, A.Lherbier, F.Liu, I.Lukacevic,
A.Martin, C.Martins, M.J.T.Oliveira, S.Ponce, Y.Pouillon, T.Rangel,
G.-M.Rignanese, A.H.Romero, B.Rousseau, O.Rubel, A.A.Shukri,
M.Stankovski, M.Torrent, M.J.Van Setten, B.Van Troeye, M.J.Verstraete,
D.Waroquier, J.Wiktor, B.Xue, A.Zhou, J.W.Zwanziger. Computer Physics
Communications 205, 106 (2016).

\[5\] P. Giannozzi, et al J.Phys.:Condens.Matter, 21, 395502 (2009)
http://dx.doi.org/10.1088/0953-8984/21/39/395502 .

\[6\] Atsushi Togo and Isao Tanaka, Scr. Mater., 108, 1-5 (2015),
DOI:10.1016/j.scriptamat.2015.07.021

\[7\] J. D. Gale, A. L. Rohl, Mol. Simul., 2003,
DOI:10.1080/0892702031000104887.

\[8\] J. E. Bertie, Glossary of Terms used in Vibrational Spectroscopy,
Handbook of Vibrational Spectroscopy. John Wiley & Sons, Ltd,
Chichester, UK, 2006.

\[9\] E. B. Wilson, J. C. Decius, P. C. Cross, B. R. Sundheim, Molecular
Vibrations: The Theory of Infrared and Raman Vibrational Spectra,
Journal of The Electrochemical Society, 102. 235C, 1955.

\[10\] T. R. Juliano, T. M. Korter, J. Phys. Chem. A, 2013,
DOI:10.1021/jp407112w.

\[11\] A. D. Burnett, J. Kendrick, C. Russell, J. Christensen, J. E.
Cunningham, A. R. Pearson, E. H. Linfield, A. G. Davies., Anal. Chem.,
2013, DOI:10.1021/ac401657r.

\[12\] A. Pereverzev, T. D. Sewell, J. Chem. Phys., 2011,
DOI:10.1063/1.3518423.

\[13\] A. Pereverzev, T. D. Sewell, D. L. Thompson, J. Chem. Phys.,
2014, DOI:10.1063/1.4866896.

\[14\] H. C. Van De Hulst, in Light Scattering by Small Particles;
Dover, New York, 1981.

\[15\] F. Wooten, in Optical Properties of Solids; Academic Press, New
York, 1972.

\[16\] X. Gonze, C. Lee, Phys. Rev. B, 1997,
DOI:10.1103/PhysRevB.55.10355.

\[17\] H. Fröhlich, in Theory of Dielectrics; Oxford University Press,
Oxford, 1948.

\[18\] L. Genzel, T. P. Martin, Surf. Sci., 1973,
DOI:10.1016/0039-6028(73)90185-4.

\[19\] C. J. Serna, M. Ocafia, J. E. Iglesias, J. Phys. C Solid St.
Phys., 1987, 20, 473--484.

\[20\] J. E. Iglesias, M. Ocana, C. J. Serna, Appl.Spectrosc, 1990, 44,
418--426.

\[21\] E. Balan, S. Delattre, D. Roche, L. Segalen, G. Morin, M.
Guillaumet, M. Blanchard, M. Lazzeri, C. Brouder, E. K. H. Salje, Phys.
Chem. Miner., 2011, DOI:10.1007/s00269-010-0388-x.

\[22\] E. Balan, M. Blanchard, J.-F. Hochepied, M. Lazzeri, Phys. Chem.
Miner., 2008, DOI:10.1007/s00269-008-0221-y.

\[23\] E. Balan, M. Lazzeri, G. Morin, F. Mauri, Am. Mineral., 2006,
DOI:10.2138/am.2006.1922.

\[24\] C. Fourdrin, E. Balan, T. Allard, C. Boukari, G. Calas, Phys.
Chem. Miner., 2009, DOI:10.1007/s00269-008-0277-8.

\[25\] A. Sihvola, in Electromagnetic Mixing Formulas and Applications;
P. J. Clarricoats and E. V. Jull, Eds.; IET, The Institution of
Engineering and Technology, Michael Faraday House, Six Hills Way,
Stevenage SG1 2AY, UK, 1999.

\[26\] M. T. Ruggiero, T. Bardon, M. Strlič, P. F. Taday, T. M. Korter,
Phys. Chem. Chem. Phys., 2015, DOI:10.1039/C5CP01195G.

\[27\] S. Giordano, J. Electrostat., 2003,
DOI:10.1016/S0304-3886(02)00199-7.

\[28\] T. G. Mackay, A. Lakhtakia, Opt. Commun., 2009,
DOI:10.1016/j.optcom.2009.03.035.

\[29\] K. Karkkainen, A. Sihvola, K. Nikoskinen, IEEE Trans. Geosci.
Remote Sens., 2001, DOI:10.1109/36.921419.

\[30\] S. Jamaian, T. G. Mackay, J. Nanophotonics, 2010,
DOI:10.1117/1.3460908.

\[31\] M. Meier and A. Wokaun, Optics Letters, 1983, 8, 581-583

\[32\] I, Peltoniemi, J. Quant. Spectroc. Radiat. Transfer, 1996, 55,
637-647

\[new33\] B. J. Sumlin, W. R. Heinson, R. K. Chakrabarty, J. Quant.
Spectrosc. Radiat. Transf., 2018, DOI:10.1016/j.jqsrt.2017.10.012.

\[34\] B. Stout, M. Nevière, E. Popov, *J. Opt. Soc. Am. A*, **2007**,
*24*, 1120--1130.

\[35\] http://www.github.com/JohnKendrick/PDielec

\[36\] <http://wiki.scipy.org>

\[37\] <https://github.com/yaml>

\[38\] http://xlsxwriter.readthedocs.io/

\[38\] M. Hellenbrandt, Crystallogr. Rev., 2004,
DOI:10.1080/08893110410001664882.

\[39\] V. G. Tsirelson, A. S. Avilov, Y. A. Abramov, E. L. Belokoneva,
R. Kitaneh, D. Feil, Acta Crystallogr. Sect. B Struct. Sci., 1998,
DOI:10.1107/S0108768197008963.

\[40\] J. P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett., 1996,
DOI:10.1103/PhysRevLett.77.3865.

\[41\] G. Kresse, D. Joubert, Phys. Rev. B, 1999,
DOI:10.1103/PhysRevB.59.1758.

\[42\] S. C. Abrahams, J. L. Bernstein, Acta Crystallogr. Sect. B
Struct. Crystallogr. Cryst. Chem., 1969, DOI:10.1107/S0567740869003876.

\[43\] C. S. McNally, D. P. Turner, A. N. Kulak, F. C. Meldrum, G.
Hyett, Chem. Commun., 2012, DOI:10.1039/C2CC14468A.

\[44\] K. Yamamoto, C.-D. Tran, H. Shimizu, K. Abe, J. Phys. Soc. Japan,
1977, DOI:10.1143/JPSJ.42.587.

\[45\] J. L. Rendon, J. E. Iglesias, C. J. Serna, Opt. Pura Y Apl.,
1981, 14, 117--122.

\[46\] S. Hayashi, N. Nakamori, H. Kanamori, J. Phys. Soc. Japan, 1979,
DOI:10.1143/JPSJ.46.176.

\[47\] D. K. Fisler, J. D. Gale, T. Cygan, Randall, Am. Mineral., 2000,
DOI:10.2138/am-2000-0121.

\[48\] B. G. Dick, A. W. Overhauser, Phys. Rev., 1958,
DOI:10.1103/PhysRev.112.90.

\[49\] D. B. DeOliveira, R. A. Laursen, J. Am. Chem. Soc., 1997,
DOI:10.1021/ja972270w.

\[50\] J.D.H. Donnay, D. Harker. Am. Mineralogist, 1937, 22, 446

\[51\] C. F. Macrae, I. J. Bruno, J. A. Chisholm, P. R. Edgington, P.
McCabe, E. Pidcock, L. Rodriguez-Monge, R. Taylor, J. van de Streek, P.
A. Wood, J. Appl. Crystallogr., 2008, DOI:10.1107/S0021889807067908.

\[52\] J.M. Hughes, M. Cameron, K.D. Crowley Am. Mineralogist 1989, 74,
870--876

\[53\] J. L. Derissen, H. J. Endeman, A. F. Peerdeman, Acta Crystallogr.
Sect. B Struct. Crystallogr. Cryst. Chem., 1968,
DOI:10.1107/S0567740868004280.

\[54\] A. Tkatchenko, M. Scheffler, Phys. Rev. Lett., 2009,
DOI:10.1103/PhysRevLett.102.073005.

\[56\] T. R. Juliano, T. M. Korter, J. Phys. Chem. A, 2015,
DOI:10.1021/jp512359p.

\[57\] J. T. Lopez Navarrete, V. Hernandez, F. J. Ramirez, Biopolymers,
1994, DOI:10.1002/bip.360340810.

\[58\] John Kendrick and Andrew Burnett (2015): *Dataset associated
with* PDielec: The Calculation of Infrared and Terahertz Absorption for
Powdered Crystals. University of Leeds. <http://doi.org/10.5518/21>

\[59\] J. Kendrick, A. Burnett. J. Comput. Chem. 2016, 37, 1491--1504.
DOI: 10.1002/jcc.24344

[^foot1]: 