.. include:: preamble.txt


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


============
Introduction
============

The Python package PDielec calculates the infrared absorption characteristics of a crystalline material supported in a non-absorbing medium. 
PDielec post processes the output of solid state quantum mechanical and molecular mechanics (MM) based codes such as VASP :cite:`Hafner2008c`, CASTEP :cite:`Clark2005d`, CRYSTAL :cite:`Dovesi2014`, Abinit :cite:`Gonze2016`, Quantum Espresso :cite:`Giannozzi2009`, Phonopy :cite:`Togo2015`, and GULP :cite:`Gale2003` to predict the infrared absorption of crystalline insulator materials whose crystal size is small compared with the wavelength of the absorbing radiation. The package is suited for the calculation of the complex, frequency dependent permittivity and its associated absorption of infrared radiation for a finely ground crystalline material dispersed in a low loss dielectric medium such KBr or Polytetrafluoroethylene (PTFE). A particular feature of the program is its ability to take into account the constant permittivity of the supporting medium and the particle shape of the material of interest through an effective medium theory. 

Using an effective medium method, the package calculates the internal electric field arising from different particle morphologies and calculates the resulting shift in absorption frequency and intensity arising from the coupling between a phonon and the internal field. The theory of the approach is described, followed by a description of the implementation within PDielec. For the specific case of a low concentration of spherical particles, calculations based on Mie scattering allow the exploration of particle size effects. A section providing several examples of its application is given.

