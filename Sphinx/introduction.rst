.. include:: preamble.txt


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, Phonopy, VASP, GULP, QE


============
Introduction
============

The Python package, PDielec, calculates the infrared absorption characteristics of a crystalline material supported in a non-absorbing medium and the absorption of single crystals (thin films or thick slabs) of absorbing material at different orientations. 
PDielec post-processes the output of solid-state quantum mechanical and molecular mechanics (MM) based codes such as VASP :cite:`Hafner2008c`, CASTEP :cite:`Clark2005d`, CRYSTAL :cite:`Dovesi2014`, Abinit :cite:`Gonze2016`, Quantum Espresso :cite:`Giannozzi2009`, Phonopy :cite:`Togo2015`, and GULP :cite:`Gale2003`.
The infrared absorption of crystalline insulator materials whose crystal size is small compared with the wavelength of the absorbing radiation is calculated using effective medium theory. The package is suited for the calculation of the complex, frequency-dependent permittivity and its associated absorption of infrared radiation for a finely ground crystalline material dispersed in a low-loss dielectric medium such as KBr or Polytetrafluoroethylene (PTFE). A particular feature of the program is its ability to take into account the constant permittivity of the supporting medium and the particle shape of the material of interest through an effective medium theory. 

Using an effective medium method, the package calculates the internal electric field arising from different particle morphologies and calculates the resulting shift in absorption frequency and intensity arising from the coupling between a phonon and the internal field. The theory of the approach is described, followed by a description of the implementation within PDielec. For the specific case of a low concentration of spherical particles, calculations based on Mie scattering allow the exploration of particle size effects. A section providing several examples of its application is given.

For spherical, isotropic particles, embedded in a supporting matrix the Mie method is available to calculate the absorption of infrared radiation.
The method uses the routines of PyMieScatt by B. Sumlin :cite:`Sumlin2018a,pymiescatt`.

The attenuated total reflection (ATR) spectrum can be calculated for a powdered material next to a high refractive index material such as diamond or germanium.

The calculation of the optical behaviour of a single crystal (thick slab or film) is performed using the pyGTM code available on GitHub :cite:`pygtm`.
This code implements a generalised transfer matrix method described by Passler et al. :cite:`Passler2020` and builds on previous publications :cite:`Passler2017,Passler2017a`.  
The agolorithms used by pyGTM have been extended to include the scattering matrix method for light transmission and reflectance :cite:`PyLlama2022`

An interface to the vibAnalysis package of Philipe Teixeira :cite:`TeixeiraVibAnal` is provided through the vibanalysis command.
The background to the theory of the analysis is given in a paper :cite:`Teixeira2019`.
The package performs Vibrational Mode Decompsition (VMD) using three computational approaches; Vibrational Mode Projection, Vibrational Mode Linear Decomposition, and Vibrational Mode Automatic Relevance Determination.
