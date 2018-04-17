.. include:: preamble.txt

.. contents::
   :local:


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


============
Applications
============

Several examples are given to illustrate applications of the package.  The calculations used to provide the data for the permittivities are sufficiently accurate to illustrate aspects of the theory. The examples are chosen to show the package being used with the QM packages CASTEP and VASP and with the MM package GULP.  All the calculations reported here have been performed with the command line interface to PDielec.

MgO using CASTEP
================

Magnesium oxide is an isotropic medium, the initial unit cell and the space group symmetry (:math:`Fm\overline{3}m`) were taken from the Inorganic Crystal Structure Database (ICSD) :cite:`Hellenbrandt2015` reference number ICSD-52026 :cite:`Tsirelson1998`. The primitive cell was optimized using CASTEP. Norm-conserving pseudo-potentials were used to represent the core electrons of magnesium and oxygen. An energy cutoff of 1000 eV was used with the PBE :cite:`Kresse1999a` density functional and a k-point spacing for the Monkhorst-Pack grid of 0.04 Å\ :superscript:`-1`. The primitive cell was optimized and a Density Functional Perturbation Theory (DFPT) calculation of the phonon spectrum at the gamma point was performed. The optimised lattice parameter was found to be 2.1234 Å, compared with the experimental value of 2.107 Å. Only 3 degenerate modes contribute to the permittivity. A summary of the results is presented in :numref:`tab-mgo-properties`.


.. table:: Experimental and Calculated Properties of MgO
   :name: tab-mgo-properties
   :widths: 1 1 1
   :column-dividers:   none none  none  
   :header-alignment: right right right 
   :column-alignment: right right right 

   +------------------------------------------+-----------------------------+--------------------------------------+
   | Property                                 | Experimental                | Calculated                           |
   +==========================================+=============================+======================================+
   | Unit cell dimensions (Å)                 | 2.107                       | 2.1234                               |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Space group                              | :math:`Fm\bar{3}m`          |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Optical permittivity                     |                             | 3.1400                               |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Static permittivity                      |                             | 10.000                               |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | **TO phonons**                           |  **Frequency** (|cm-1|)     | **Intensity** (|D2A-2amu-1|)         | 
   +------------------------------------------+-----------------------------+--------------------------------------+
   | T                                        | 388.3                       | 9.29                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | **LO phonons**                           |                             |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | T (001)                                  | 693.7                       |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+



Because MgO is isotropic with only a single frequency contributing to the permittivity, it makes a useful example application to illustrate several features of PDielec. The real and imaginary frequency dependent permittivities are shown in :numref:`fig-mgo-permittivity`, where a damping factor (:math:`\sigma`) of 10 |cm-1| has been used. In the Figure the real permittivity at zero frequency corresponds to the static permittivity in :numref:`tab-mgo-properties`, and at frequencies above the absorption at 388 |cm-1| the permittivity tends to the optical permittivity as the frequency increases. The real permittivity has zero values at 388.3 and 693.7 |cm-1| which are the TO and LO frequencies respectively.

.. _fig-mgo-permittivity:

.. figure:: ./_static/Figures/MgO_Permittivity.*
   :scale: 90%

   Permittivity of MgO

    
.. _fig-mgo-real:

.. figure:: ./_static/Figures/MgO_Real_Imaginary.*
   :scale: 90%

   Real and imaginary permittivities of a 10% volume fraction of MgO spheres in PTFE using the Maxwell-Garnett method
    
Using the Maxwell-Garnett mixing rule, :numref:`fig-mgo-real` shows the calculated permittivities of a 10% volume fraction of MgO spheres in a supporting medium with a frequency independent permittivity of 2.0, which would be typical of a material such as PTFE. Due to the dilution effect the real component has shifted to a base line value close to 2, and the absorption, as indicated by the maximum in the imaginary component has shifted by about 150 |cm-1| to 550 |cm-1|. 

The effect of volume fraction on the predicted molar absorption coefficient, using the Maxwell-Garnett mixing rule, is shown in :numref:`fig-mgo-vf`. The lowest volume fraction of MgO gives the largest shift of the absorption peak to high frequency. As the volume fraction increases the mixing rule predicts a broadening of the absorption, whilst the peak in the molar absorption coefficient moves to lower frequency. At the highest loading (f=0.9) the maximum absorption occurs quite close to the TO frequency. The Maxwell-Garnett mixing rule is regarded as being appropriate for low volume fractions and so should not be used for interpreting results in which higher volume fractions of absorbing media have been used :cite:`Sihvola`.

.. _fig-mgo-vf:

.. figure:: ./_static/Figures/MgO_MG_Volume_Fraction.*
   :scale: 90%

   Effect of volume fraction on the Maxwell-Garnett molar absorption coefficient of MgO spheres in PTFE

:numref:`fig-mgo-vf-bg` shows the same plot for the Bruggeman mixing rule. At low
volume fractions the Bruggeman mixing rule predicts a similar absorption
to the Maxwell-Garnett. Indeed as the volume fraction approaches zero
the two rules predict the same absorption characteristics. However, even
at the relatively low 1% loading, the Bruggeman mixing rule shows
additional broadening of the peak, the shape of the absorption peak has
lost its Lorentzian characteristic shape as can be seen clearly in
:numref:`fig-mgo-vf-bg`.  At 10% loading the Bruggeman predicted absorption is broad with
the peak shifted to lower wavenumber. This broadening increases with
increased loading until, at the higher loadings, the TO peak begins to
dominate the absorption.


.. _fig-mgo-vf-bg:

.. figure:: ./_static/Figures/MgO_Bruggeman_Volume_Fraction.*
   :scale: 90%

   Effect of volume fraction on the Bruggeman molar absorption coefficient of MgO spheres in PTFE

.. _fig-mgo-varying-permittivity:

.. figure:: ./_static/Figures/MgO_Varying_Permittivity.*
   :scale: 90%

   Effect of volume fraction on the Bruggeman molar absorption coefficient of MgO spheres in PTFE

:numref:`fig-mgo-varying-permittivity` shows the effect of varying the permittivity of the supporting medium. The calculations were performed on spherical MgO particles with a 1% volume fraction. The lowest permittivity is that of a vacuum (or air) and shows the highest shift of the absorption maximum to higher frequencies. Increasing the permittivity lowers the shift until it becomes quite small. A similar effect is seen for the Bruggeman mixing model. However, the absorption resulting for particles in a low dielectric medium is considerable broader than that seen in the Maxwell-Garnet case. This broadening reduces as the permittivity of the medium increases (see :numref:`fig-mgo-varying-permittivity-brug`).

.. _fig-mgo-varying-permittivity-brug:

.. figure:: ./_static/Figures/MgO_Varying_Permittivity_Bruggeman.*
   :scale: 90%

   The Bruggeman molar absorption coefficients of spherical MgO particles, 1% volume fraction, embedded in media of varying permittivites


ZnO using VASP
==============

Zinc oxide crystallizes in space group :math:`P6_3mc` (wurtzite). All calculations were performed by VASP :cite:`Hafner2008c` using projector augmented-wave PAW :cite:`Kresse1999a` pseudo-potentials, the PBE :cite:`Perdew1996a` density functional, an energy cutoff of 600 eV and a k-point resolution of approximately 0.1 Å\ :superscript:`-1`. The initial unit cell was taken from the ICSD :cite:`Hellenbrandt2015` with code ICSD-26170 :cite:`McNally2012a`. The unit cell and atom positions were optimized using
VASP and the permittivity was calculated using DFPT and the results reported in :numref:`tab-zno-properties`. Only two of the bands showed any significant intensity, a doubly degenerate band (E) with a TO frequency of 372.1 |cm-1| and a non-degenerate band (A) with a TO frequency of 350.0 |cm-1|. The LO frequency of the non-degenerate band is shifted to 502.0 |cm-1| for a wave-vector with direction (001), whilst the degenerate modes are unaffected. In the case of the (010) direction the LO frequency of one of the E modes is shifted to 511.2 |cm-1|. It is known that ZnO can crystallize with a plate morphology :cite:`Yamamoto1977` with the (001) surface dominant. Calculations of the molar absorption were performed for a sphere, plate and needle like shapes with the unique directions of the plate and the needle being normal to the (001) surface. A volume fraction of 1% was chosen for these calculations and the predicted molar absorption coefficients for the Maxwell-Garnett mixing rule is shown in :numref:`fig-zno`.


.. table:: Experimental and Calculated Properties of ZnO
   :name: tab-zno-properties
   :widths: 1 1 1
   :column-dividers:   none none  none  
   :header-alignment: right right right 
   :column-alignment: right right right 

   +------------------------------------------+-----------------------------+--------------------------------------+
   | Property                                 | Experimental                | Calculated                           |
   +==========================================+=============================+======================================+
   | Unit cell dimensions a,b (Å)             | 3.250                       | 3.295                                |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Unit cell dimensions a   (Å)             | 5.207                       | 5.285                                |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Space group                              | :math:`P6_3mc`              |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Optical permittivity   xx,yy             |                             | 5.09                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Optical permittivity   zz                |                             | 6.00                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Static permittivity    xx, yy            |                             | 10.83                                |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Static permittivity    zz                |                             | 11.67                                |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | **TO phonons**                           |  **Frequency (|cm-1|)**     | **Intensity** (|D2A-2amu-1|)         | 
   +------------------------------------------+-----------------------------+--------------------------------------+
   | A                                        | 350.0                       | 17.1                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | E                                        | 372.1                       | 16.4                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | **LO phonons**                           |                             |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | (001) A                                  | 502.0                       |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | (010) E                                  | 511.2                       |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+


.. _fig-zno:

.. figure:: ./_static/Figures/ZnO-Maxwell_Garnett.*
   :scale: 90%

   ZnO: Effect of shape on absorption

For the Maxwell-Garnett mixing rule the sphere morphology results in the two absorption peaks shifting from their TO positions to higher wavenumber by about 80 |cm-1|. The plate morphology results in one of the peaks moving to higher wavenumber by about 130 |cm-1|, whilst the other remains at the TO position. The Maxwell-Garnett results are in close accord with some experimental results by Yamamoto et al :cite:`Yamamoto1977` who measured the infrared spectrum of ZnO smoke particles and observed peaks in the absorption at 380, 530 and 550 |cm-1|. Previous work :cite:`Rendon1981,Hayashi1977a` have also used effective medium theory to explain the observed spectrum.

Calcite using GULP
==================

Calcite is the most stable polymorph of calcium carbonate and the crystal structure belongs to the :math:`R\overline{3}c` space group. The force field and atomic structures used here are described in detail in work by Fisler :cite:`Fisler2000`. Briefly, the oxygen ions are described using a core-shell model :cite:`Dick1958`. The carbon - oxygen potential of the carbonate is taken to be a Morse potential and an additional 3 atom potential is used to maintain the O-C-O angle at 120\ :superscript:`o`. The van der Waals interactions between non bonded atoms are taken to be Buckingham potentials and the charges on the calcium, carbon and oxygen ions are +2, +1.3435 and -1.1145 respectively. The shell charge of the oxygen ion is -2.133 and the spring constant for the core-shell interaction is 52.74 eV/Å\ :superscript:`2`. The unit cell was optimized using the primitive unit cell and the full space group symmetry. The calculation of the phonon spectrum was performed without symmetry but still using the primitive cell of the lattice. A summary of the calculated properties is given in :numref:`tab-calcite-properties`.


.. table:: Experimental and Calculated Properties of Calcite
   :name: tab-calcite-properties
   :widths: 1 1 1
   :column-dividers:   none none  none  
   :header-alignment: right right right 
   :column-alignment: right right right 

   +------------------------------------------+-----------------------------+--------------------------------------+
   | Property                                 | Experimental                | Calculated                           |
   +==========================================+=============================+======================================+
   | Unit cell dimensions a,b,c (Å)           | 6.375                       | 6.376                                |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Unit cell angles                         | 46.1                        | 46.1                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Space group                              | :math:`R\overline{3}c`      |                                      |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Optical permittivity   xx,yy             |                             | 1.91                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Optical permittivity   zz                |                             | 2.00                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Static permittivity    xx, yy            |                             | 6.70                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | Static permittivity    zz                |                             | 7.10                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | **TO phonons**                           |  **Frequency** (|cm-1|)     | **Intensity** (|D2A-2amu-1|)         | 
   +------------------------------------------+-----------------------------+--------------------------------------+
   | E                                        | 114.8                       | 2.39                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | A                                        | 127.4                       | 3.36                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | A                                        | 249.3                       | 1.23                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | E                                        | 320.7                       | 5.82                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | A                                        | 338.1                       | 4.14                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | E                                        | 620.1                       | 3.38                                 |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | A                                        | 732.0                       | 26.89                                |
   +------------------------------------------+-----------------------------+--------------------------------------+
   | E                                        | 1463.6                      | 16.97                                |
   +------------------------------------------+-----------------------------+--------------------------------------+


:numref:`fig-calcite` shows the results of analysis of the results using PDielec. The damping parameter used in the calculation was a value of 5 |cm-1|. A 10% volume fraction was used with sphere and plate morphologies for the particles. The unique axis for the plate was taken to be the normal to the (211) surfaces in the primitive cell axes (or the {104} surfaces in the standard unit cell). Such surfaces define the rhombohedral faces commonly seen in calcite crystals :cite:`DeOliveira1997`. :numref:`fig-calcite` shows that the doubly degenerate TO absorption peak at 620 |cm-1| is not significantly affected by spherical particles and there is a small shift to higher frequencies in the case of plate-like particles. The non-degenerate TO transition at 732 |cm-1|, which corresponds to motion of the carbon atom of the carbonate along the unique direction of the slab, shows a shift to 786 and 819 |cm-1| for the sphere and plate respectively. The doubly degenerate peak at 1463 |cm-1| is shifted to 1480 |cm-1| by spherical particles and is split by plate-like particles with one component which shifts to 1491 |cm-1| .


.. _fig-calcite:

.. figure:: ./_static/Figures/Calcite_AP_Plate211_Sphere.*
   :scale: 90%

   Calculated Maxwell-Garnett absorption spectrum of 10% volume fraction of calcite in PTFE

Fluoroapatite using VASP
========================

The line shapes of the infrared absorption of apatite and fluoroapatite were examined extensively by Balan *et al* :cite:`Balan2008b`. Their calculations included the effect of crystallite habit on the spectrum and the results reported here are similar to their conclusions. The method used by Balan *et al*. is an infinitely dilute Maxwell-Garnett model, so the only difference between the methods used by them and those reported here using PDielec are the incorporation of the volume fraction into the theory and the use of an ellipsoidal shape for comparison with the other shapes.

All calculations were performed by VASP :cite:`Hafner2008c` using projector augmented-wave PAW :cite:`Kresse1999a` pseudo-potentials, the PBE :cite:`Perdew1996a` density functional, an energy cutoff of 600 eV and a k-point resolution of approximately 0.1 Å\ :superscript:`-1`. :numref:`tab-fluoroapatite-properties` summarises the results of the calculations. Only the 3 highest frequency bands are reported and discussed. The TO intensity of the highest frequency band at 1038 |cm-1| is low and will not be discussed further. The Bravais Friedel Donnay Harker (BFDH) :cite:`Donnay1937` crystal habit of the optimized crystal is shown in :numref:`fig-fluoroapatite-morphology`. The habit was calculated using the Mercury software package. :cite:`Macrae2008`. The BFDH crystal habit is often used to give an idea of the likely important faces of a crystal. It uses only the crystal lattice and space group to determine the crystal morphology. :numref:`fig-fluoroapatite-morphology` shows that the {100} surfaces form a tube which are capped by the {011} surfaces. The effect of different particle shapes on the predicted spectrum is shown in :numref:`fig-fluoroapatite-absorption`. The calculations of the spectra were performed with a damping parameter (σ) of 2 |cm-1|. The ellipsoid was chosen to have an aspect ratio, a/b, of 2 and a principle axis along [001], which was compatible with the morphology predicted by the BDFH method. The two TO absorption frequencies at 981 and 986 |cm-1| have A and E symmetry respectively. Spherical crystallites result in three absorption peaks at around 1000, 1010 and 1015 |cm-1|. Needle shaped crystallites leave the A symmetry TO absorption peak at 981 |cm-1| unaffected, but shift and split the E symmetry TO peak to 1020 and 1046 |cm-1|. A plate morphology with (100) surfaces results in the A and one component of the E TO absorption peak remaining at the TO frequencies, with the other component of the E shifting 85 |cm-1| to 1075 |cm-1|. The ellipsoidal morphology show three shifted peaks at 1000, 1018 and 1045 |cm-1|. These results are consistent with those of Balan *et al*. :cite:`Balan2008b`, who gave detailed results for hydroxyapatite.


.. table:: Experimental and Calculated Properties of Fluoroapatite
   :name: tab-fluoroapatite-properties
   :widths: 1 1 1
   :column-dividers:   none none  none  
   :header-alignment: right right right 
   :column-alignment: right right right 

   +------------------------------------------+---------------------------------+------------------------------+
   | Property                                 | Experimental :cite:`Hughes1989` | Calculated                   |
   +==========================================+=================================+==============================+
   | Unit cell dimensions a,b (Å)             | 9.417                           | 9.447                        |
   +------------------------------------------+---------------------------------+------------------------------+
   | Unit cell dimensions c (Å)               | 6.875                           | 6.926                        |
   +------------------------------------------+---------------------------------+------------------------------+
   | Space group                              | :math:P6_3m`                    |                              |
   +------------------------------------------+---------------------------------+------------------------------+
   | Optical permittivity   xx,yy             |                                 | 2.891                        |
   +------------------------------------------+---------------------------------+------------------------------+
   | Optical permittivity   zz                |                                 | 2.894                        |
   +------------------------------------------+---------------------------------+------------------------------+
   | Static permittivity    xx, yy            |                                 | 12.081                       |
   +------------------------------------------+---------------------------------+------------------------------+
   | Static permittivity    zz                |                                 | 8.841                        |
   +------------------------------------------+---------------------------------+------------------------------+
   | **TO phonons**                           |  **Frequency** (|cm-1|)         | **Intensity** (|D2A-2amu-1|) | 
   +------------------------------------------+---------------------------------+------------------------------+
   | A                                        | 981.8                           | 112.6                        |
   +------------------------------------------+---------------------------------+------------------------------+
   | E                                        | 986.3                           | 101.0                        |
   +------------------------------------------+---------------------------------+------------------------------+
   | E                                        | 1038.1                          | 7.92                         |
   +------------------------------------------+---------------------------------+------------------------------+


.. _fig-fluoroapatite-morphology:

.. figure:: ./_static/Figures/Fluoroapatite_morphology.*
   :scale: 90%

   BFDH Morphology of fluoroapatite


.. _fig-fluoroapatite-absorption:

.. figure:: ./_static/Figures/Fluoroapatite_absorption.*
   :scale: 90%

   Calculated Maxwell-Garnett absorption spectra of 10% fluoroapatite in PTFE


L-aspartic Acid using CASTEP
============================

L-aspartic acid is a zwitterion in the solid state and so the shape of the particles used in the measurement of IR and THz spectra maybe important. The starting geometry for optimization of the unit cell and molecular structure of L-aspartic acid was taken from Derissen et al :cite:`Derissen1968`. The PBE :cite:`Perdew1996a` functional was used with a plane wave energy cutoff of 1000 eV and norm conserving pseudo-potentials. A dispersion correction using the Tkatchenko-Scheffler scheme :cite:`Tkatchenko2009c` available in CASTEP was applied for both the geometry optimisation and the calculation of the phonon spectrum at the gamma point, with a value S~6~ scaling factor :cite:`Juliano2015` of 1.0. A summary of the results of the calculations is shown in :numref:`tab-aspartic-properties`.


.. table:: Experimental and calculated properties of L-aspartic acid
   :name: tab-aspartic-properties
   :widths: 1 1 1
   :column-dividers:   none none  none  
   :header-alignment: right right right 
   :column-alignment: right right right 

   +------------------------------------------+-----------------------------------+------------------------------+
   | Property                                 | Experimental :cite:`Derissen1968` | Calculated                   |
   +==========================================+===================================+==============================+
   | Unit cell dimensions a (Å)               | 7.617                             | 7.597                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Unit cell dimensions b (Å)               | 6.982                             | 7.028                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Unit cell dimensions c (Å)               | 5.142                             | 5.113                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Unit cell angle :math:`\beta`            | 99.84                             | 98.77                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Space group                              | :math:`P2_1`                      |                              |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Optical permittivity   xx                |                                   | 2.68                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Optical permittivity   yy                |                                   | 2.20                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Optical permittivity   zz                |                                   | 2.56                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Static permittivity    xx                |                                   | 4.58                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Static permittivity    yy                |                                   | 3.65                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   | Static permittivity    zz                |                                   | 3.65                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   | **TO phonons**                           |  **Frequency** (|cm-1|)          | **Intensity** (|D2A-2amu-1|)  | 
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 84.5                              | 0.120                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 104.7                             | 0.202                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 106.0                             | 0.243                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 115.3                             | 0.474                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 137.3                             | 0.617                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 1290.0                            | 55.0                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 2945.9                            | 102.8                        |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 2947.3                            | 48.2                         |
   +------------------------------------------+-----------------------------------+------------------------------+
   |                                          | 3053.7                            | 44.1                         |
   +------------------------------------------+-----------------------------------+------------------------------+


The THz spectrum of L-aspartic acid has been reported by Juliano and Korter :cite:`Juliano2015` in the frequency range 0-90 |cm-1|. The infrared spectrum has been reported and assigned by Lopez *et al* :cite:`LopezNavarrete1994`. :numref:`fig-aspartica`, :numref:`fig-asparticb` and :numref:`fig-asparticc` shows the calculated absorption spectra for L-aspartic acid for three frequency ranges. The calculation of the spectra used the Maxwell-Garnett mixing rule with a 10% volume fraction of L-aspartic acid in PTFE and for comparison the TO mixing rule. A damping factor of 2 |cm-1| was used. Spherical and a variety of plate-like inclusions were used to illustrate their effect on the absorption spectra. :numref:`fig-aspartica` shows the frequency range from 60-130 |cm-1| which is that covered by THz spectroscopy. The shifts observed for the different particle morphologies are not large, but the change in intensities is significant. The molecular motions associated with phonons at these frequencies tend to be whole molecule motion involving rotation. :numref:`fig-asparticb` shows the frequency range from 1260-1340 |cm-1|. In this frequency range bending of the carboxylate anion contributes to the spectrum significantly. The three different plate morphologies show different and significant shifts in the TO absorption peak at 1290 |cm-1|. The spherical morphology shows a shift of around 25 |cm-1| to higher wavenumber. :numref:`fig-asparticc` shows the spectra in the frequency range 2900-3100 |cm-1|, which corresponds to the motion of O-H (below 2980 |cm-1|) and N-H (above 2980 |cm-1|) stretching. The effect of the different possible crystal morphologies is large with shifts to higher frequency of up to 50 |cm-1|. The spectra below 3000 |cm-1| arises from two TO absorptions at 2946 and 2947 |cm-1|. Because the motions associated with each mode interact differently with the internal field within each crystal they give rise to different shifts producing more complex spectra.

.. _fig-aspartica:

.. figure:: ./_static/Figures/Aspartic_thz_a.*
   :scale: 90%

   Aspartic Acid: Frequency range 60-130 |cm-1|

.. _fig-asparticb:

.. figure:: ./_static/Figures/Aspartic_thz_b.*
   :scale: 90%

   Aspartic Acid: Frequency range 1260-1340 |cm-1|

.. _fig-asparticc:

.. figure:: ./_static/Figures/Aspartic_thz_c.*
   :scale: 90%

   Aspartic Acid: Frequency range 2900-3100 |cm-1|


MgO Example using Mie Scattering
================================

:numref:`fig-mgo-mie` compares a Mie scattering calculation with the results from Maxwell-Garnett effective medium theory. The same data set was used for the CASTEP, MgO example. A volume fraction of 1% was used with a small sphere radius (0.1 μm) and a broadening of 5 |cm-1| embedded in a matrix of PTFE. A power expansion in the size parameter of the Mie expressions
indicates that for small sizes of particles, the Mie and the Maxwell-Garnett methods should be the same. This is verified in Figure 13.


.. _fig-mgo-mie:

.. figure:: ./_static/Figures/MgO_Mie_Maxwell.*
   :scale: 90%

   Comparison of Mie and Maxwell methods. 1% volume fraction of MgO in PTFE, sphere radius of 0.1 μm and a broadening of 5 |cm-1|

To better understand what makes particles large or small :numref:`tab-mgo-mie-sizes` shows the dimensionless size parameter, :math:`x`, as a function of wavenumber and of sphere radius. It has been assumed that the supporting medium is PTFE. Since the power expansion of the size parameters leads to terms which are quadratic in :math:`x`, it should be expected that when :math:`x` is less than about 0.1 μm, the particles can be considered small. It can be seen that particles less than 0.01 μm are small over the range of frequencies considered. But while 0.1 μm particles are small in the THz regime and low frequency infrared they should not be considered small over the more extended infrared frequencies. 1μm particles should be considered large for both the THz and the extended infrared.

.. table:: Variation of size parameter with wavenumber and radius of sphere
   :name: tab-mgo-mie-sizes
   :widths: 1 1 1 1 1 1
   :column-dividers:   none none none  none  none  none
   :header-alignment: right right right right right right
   :column-alignment: right right right right right right

   +------------+-----------------------+--------+--------+--------+---------+
   | Wavenumber | Radius of sphere (μm)                                      |
   +------------+-----------------------+--------+--------+--------+---------+
   | (cm-1)     | 0.001                 | 0.01   | 0.1    | 1      | 10      |
   +============+=======================+========+========+========+=========+
   | 100        | 0.0009                | 0.0089 | 0.0888 | 0.8884 | 8.8841  |
   +------------+-----------------------+--------+--------+--------+---------+
   | 200        | 0.0018                | 0.0178 | 0.1777 | 1.7768 | 17.7682 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 300        | 0.0027                | 0.0267 | 0.2665 | 2.6652 | 26.6523 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 400        | 0.0036                | 0.0355 | 0.3554 | 3.5536 | 35.5364 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 500        | 0.0044                | 0.0444 | 0.4442 | 4.4420 | 44.4204 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 600        | 0.0053                | 0.0533 | 0.5330 | 5.3305 | 53.3045 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 700        | 0.0062                | 0.0622 | 0.6219 | 6.2189 | 62.1886 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 800        | 0.0071                | 0.0711 | 0.7107 | 7.1073 | 71.0727 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 900        | 0.0080                | 0.0800 | 0.7996 | 7.9957 | 79.9568 |
   +------------+-----------------------+--------+--------+--------+---------+
   | 1000       | 0.0089                | 0.0888 | 0.8884 | 8.8841 | 88.8409 |
   +------------+-----------------------+--------+--------+--------+---------+


Figure :numref:`fig-mgo-mie-sizes` shows how the Mie predictions change as the particle radius changes from 0.2 to 1.6 μm. As the particle size increases the peak above 500 |cm-1| splits into two. One broader peak which moves to lower frequency as the particle size increases and the peak at about 550 |cm-1| which looses intensity as the particle size increases. There is also the onset of absorption at 388 |cm-1| which corresponds to the bulk TO modes.

.. _fig-mgo-mie-sizes:

.. figure:: ./_static/Figures/MgO_Mie_Sizes.*
   :scale: 90%

   Variation in absorption calculated by the Mie method for different radii of spheres. 1% volume fraction of MgO in PTFE and a broadening of 5 |cm-1|

Figure :numref:`fig-mgo-mie-sizes-big` shows the effect of increasing the particle size further. More structure appears in the absorption, with increasing absorption around
the bulk TO frequency. Above 4.0 μm there is more low frequency structure
appearing, below 300 |cm-1|.

.. _fig-mgo-mie-sizes-big:

.. figure:: ./_static/Figures/MgO_Mie_Sizes_Big.*
   :scale: 90%

   Variation in absorption calculated by the Mie method for larger sphere radii. 1% volume fraction of MgO in PTFE and a broadening of 5 |cm-1|


ZnO Example using Mie Scattering
================================

ZnO is an anisotropic material, so the treatment described here using Mie scattering is an approximation. However, the permittivity constant tensor is diagonal due to the space group symmetry of the crystal. 

.. _fig-zno-maxwell-garnett:

.. figure:: ./_static/Figures/ZnO-Maxwell_Garnett.*
   :scale: 90%

   ZnO spheres in PTFE using Maxwell Garnett

.. _fig-zno-mie:

.. figure:: ./_static/Figures/ZnO_Mie.*
   :scale: 90%

   ZnO spheres using Mie


Figures :numref:`fig-zno-maxwell-garnett` and :numref:`fig-zno-mie` compare the capabilities of the Maxwell-Garnett and Mie methods for describing volume fraction effects. Figure :numref:`fig-zno-maxwell-garnett` shows that Maxwell-Garnett predicts a lowering of intensity and frequency of the high frequency peak as the volume fraction is increased. Figure :numref:`fig-zno-mie` illustrates that Mie theory shows no effect of the change in volume fraction. This is to be expected as the theory assumes that each sphere is isolated and not affecting the other spheres around it. It should be pointed out that the figures are plotting molar absorption coefficients. The actual absorption would increase with volume fraction of ZnO.

Figure :numref:`fig-zno-mie-sizes` shows that the variation of the Mie scattering with sphere radius follows a similar pattern to that observed in MgO, though slightly more complex. The initial peaks at about 440 and 460 |cm-1| broaden and shift to lower frequencies as the particle size increases. Bulk bands around 350 and 372 |cm-1| can be seen which grown in intensity and shift to lower frequencies as the particle size increases.

.. _fig-zno-mie-sizes:

.. figure:: ./_static/Figures/ZnO_Mie_Sizes.*
   :scale: 90%

   ZnO spheres using Mie
   Mie scattering of 10% volume fraction ZnO spheres in PTFE using a line broadening factor of 5 |cm-1|. 
