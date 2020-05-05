The PDielec Package
===================
Authors: John Kendrick and Andrew Burnett

The Python package, PDielec calculates the infrared absorption
characteristics of a crystalline material supported in a non absorbing
medium by post processesing the output of solid state quantum mechanical
and molecular mechanical calculations of the phonons or dielectric
response of the crystalline material. The package calculates the
internal electric field arising from different particle morphologies and
calculates the resulting shift in absorption frequency and intensity
arising from the coupling between a phonon and the internal field. The
theory of the approach has been published. Any use of the package should
cite; PDielec: The calculation of infrared and terahertz absorption for
powdered crystals John Kendrick and Andrew D. Burnett Journal of
Computational Chemistry 2016, vol 27, 1491-1504 DOI: 10.1002/jcc.24344

Introduction
============

The molecular and solid state quantum mechanical calculations of
response properties such as the frequencies and intensities of infrared
(IR) and terahertz (THz) radiation absorption has become generally
available in many molecular and solid state computer programs. A common
approach is to assume the harmonic approximation and calculate the mass
weighted force constant matrix (for molecules) or the dynamical matrix
at the gamma point (for periodic solids). Diagonalisation of the matrix
gives the frequencies for absorption and the normal modes (molecules) or
phonon displacements (periodic solids). The calculation of the
absorption intensity for each mode requires the calculation of the
change in dipole moment caused by the displacement of the atoms for that
mode. For solids where there is a large separation of charge, there can
be a large coupling between a phonon mode and the internal field within
a particle resulting from its morphology. The PDielec program is written
in Python and post processes the output of solid state quantum
mechanical (QM) and molecular mechanics (MM) based codes such as VASP,
CASTEP, CRYSTAL14, Abinit, QuantumEspresso, Phonopy and GULP to predict
the infrared absorption of crystalline insulator materials whose crystal
size is small compared with the wavelength of the absorbing radiation.
The package is suited for the calculation of the complex, frequency
dependent permittivity and its associated absorption of infrared
radiation for a finely ground crystalline material dispersed in a low
loss dielectric medium such KBr or PTFE. A particular feature of the
program is its ability to take into account the constant permittivity of
the supporting medium and the particle shape through an effective medium
theory. Calculation of ATR spectra has been included in v5.0. More
details as to the theory are available here;
<https://johnkendrick.github.io/PDielec/>

Installation
============

Full documentation of the program and its installation is available
here; <https://johnkendrick.github.io/PDielec/>
