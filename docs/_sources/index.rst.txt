.. include:: preamble.txt

:Authors:
    John Kendrick, Andrew Burnett  

:EMail:
    john@kendrick.me.uk, a.d.burnett@leeds.ac.uk  

.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP, Phonopy, QE

=======
PDielec
=======

This is the documentation for the PDielec package.  The package takes output from DFT calculations and calculates a material's infrared and terahertz response.  
For single crystals (thin films or slabs), a generalized transfer matrix or a scattering matrix method is used.  For powdered material in a support matrix, an effective medium theory is used which can take into account the effect of shape and particle size.

.. toctree::
   :caption: Software
   :maxdepth: 2

   introduction
   installation
   pdgui
   software

.. toctree::
   :caption: Theory
   :maxdepth: 2

   theory_powder
   theory_single_crystal
   analysis

.. toctree::
   :caption: Applications
   :maxdepth: 2

   application_notes_1
   application_notes_2
   application_notes_3

.. toctree::
   :caption: API
   :maxdepth: 3

   autoapi/index
   autoapi/PDielec/pdgui/index
   autoapi/PDielec/preader/index
   autoapi/PDielec/pdmake/index

.. toctree::
   :caption: References
   :maxdepth: 2

   zreferences


Search
======

* :ref:`search`
