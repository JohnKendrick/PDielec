.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP, FHI-Aims, Phonopy, QE


=========================
Other Software Components
=========================

preader
=======

preader is a 'helper' program that uses the underlying modules of PDielec to read output files and summarise the results of several MM/QM packages.  The program can be used to perform some straightforward calculations.
For instance, projection of any remaining centre-of-mass motion of the crystal can be performed to make sure that there are three zero frequencies.  Also, the masses used in the calculation of the dynamical matrix can be altered.
Unlike PDGui it is not necessary to have performed a full calculation of the dynamical matrix.  In the majority of cases, preader will read geometry optimisation runs.

Command options
---------------

Examples of data sets for these packages are included with the distribution and can be found in the Examples/'Package'/preader directory. The program is run from the command line. There are several command options and these are summarised below. Some options may be repeated.

  - \-program program_name
       | Program_name can be “abinit”,  "aims", “castep”, “crystal”, “gulp”, “qe”, “experiment”, "phonopy", “vasp”, “finite_field”, or “auto” and specifies the program which generated the results to be analysed
  - \-neutral
       | Impose neutrality on the Born charge matrices
  - \-spectroscopy infrared|raman
       | Select IR intensities (default) or Raman activities; Raman requires tensor data.

  - \-nocalculation
       | No calculations are performed.  This results in a single line of output with just information from the program.  If -eckart, -mass, -masses, -neutral or -crystal have -hessian crystal have been specified they will be ignored
  - \-eckart
       | Translational modes are  projected out of the hessian before diagonalisation
  - \-masses average
       | The element mass_definition can be either “program”, “average” or “isotopic”, meaning that the masses used in the calculation of the frequencies are either taken from the QM program or are the average of the isotope abundancies or are the most abundant isotope mass.
  - \-mass element mass
       | The atomic mass of the element is set to mass.  This can be used to explore the effect of isotope substitution on the calculated frequencies


Examples
--------

.. code-block:: bash

         preader -program vasp `find . -name OUTCAR` > results.csv

This reads all the VASP OUTCAR files in the current and any of its subdirectories and summarises the results to results.csv.

.. code-block:: bash

         preader -program castep -eckart `find . -name \*.castep` > results.csv

This reads all the CASTEP output files in the current and any of its subdirectories and summarises the results to results.csv.  For each file, the centre-of-mass motion of the crystal is projected.  The results file contains both the unprojected and the projected results.

There are examples of preader being used in the Examples/'Package'/preader subdirectories of the distribution of PDielec.

aims-pdielec-born
=================

To generate the Born charges from an FHI-Aims calculation of infrared vibrational frequencies and intensities, the aims-pdielec-born command reads the single point FHI-Aims calculations of forces and polarisation at the displaced geometries created by the get_vibrations.py script from the FHI-Aims distribution.
The aims-pdielec-born writes a single file BORN_PDIELEC containing the Born charges for all atoms in the unit cell.

In the example given below the FHI-Aims command ``python get_vibrations.py na2so42`` was used to generate the perturbations.

.. code-block:: bash

         aims-pdielec-born na2so42

The command line options are:

+----------------------------+--------------------------------------------------------------------------------------------------+
| Command line option        | Description                                                                                      |
+============================+==================================================================================================+
| name                       | "name" is the string used in the get_vibrations.py command to identify the set of displacements. |
+----------------------------+--------------------------------------------------------------------------------------------------+
| | -d delta                 | | Define the step used to generate the perturbed geometries                                      |
| | --delta delta            | | the default delta is 0.0025                                                                    |
+----------------------------+--------------------------------------------------------------------------------------------------+
| | -p dir                   | | Set the directory name holding the calculation of permittivity                                 |
| | --permittivity dir       | | the default is "Dielectric"                                                                    |
+----------------------------+--------------------------------------------------------------------------------------------------+
| | -h                       | | Show some help information                                                                     |
| | --help                   | |                                                                                                |
+----------------------------+--------------------------------------------------------------------------------------------------+

phonopy-pdielec-born
====================

The BORN file calculated by Phonopy holds the Born charges for the symmetry unique atoms in the unit-cell.
To allow PDielec to use these Born charges the phonopy-pdielec-born command reads the phonopy.yaml and BORN files and writes out the full set of charges to the terminal.
PDGui expects the Born charges to be in the file BORN_PDIELEC.  An example of the use of the script is given below.

.. code-block:: bash

       phonopy-pdielec-born  BORN_PDIELEC


The command line options are:

+----------------------------+----------------------------------------------------------+
| Command line option        | Description                                              |
+============================+==========================================================+
| | name                     | | Output the charges to a file called name               |
+----------------------------+----------------------------------------------------------+
| | -input name              | | Set the name of the Phonopy Born charge file           |
| | --input name             | | the default is "BORN"                                  |
+----------------------------+----------------------------------------------------------+
| | -output name             | | Send the output to this file                           |
| | --output name            | | the default is to send the output to standard out      |
+----------------------------+----------------------------------------------------------+
| | -yaml name               | | Set the name of the phonopy.yaml file                  |
| | --yaml name              | | the default is "phonopy.yaml"                          |
+----------------------------+----------------------------------------------------------+
| | -symprec float           | | Define the symmetry threshold                          |
| | --symprec float          | | the default 1.0E-5                                     |
+----------------------------+----------------------------------------------------------+
| | -nosym                   | | Do not use any symmetry                                |
+----------------------------+----------------------------------------------------------+
| | -h                       | | Show some help information                             |
| | --help                   | |                                                        |
+----------------------------+----------------------------------------------------------+



vibanalysis
===========

The vibanalysis command provides an interface the vibAnalysis package of Philipe Teixeira.  This package provides Vibrational Mode Decomposition of the phonon modes using a variety of methods including; Vibrational Mode Projection (VMP), Vibrational Mode Linear Decomposition (VMLD) and Vibrational Mode Automatic Relevance Determination (VMARD).  By default VMARD is used.

The package options can be determined using;

.. code-block:: bash

       vibanalysis OUTCAR
       vibanalysis phonon.castep

which will analyse a VASP and a CASTEP phonon calculation output to produce a summary in the files OUTCAR.nma and phonon.nma respectively.

pdmake
======
Use ``./pdmake`` from the repository root, or the installed ``pdmake``
command. Running it without arguments prints help. Run named test suites from
a checkout containing ``Examples/`` and ``PDielec/Tests/``.

.. list-table:: Commands and options
   :header-rows: 1
   :widths: 35 65

   * - Command or option
     - Behaviour
   * - ``tests (or test)``
     - Run the complete PDielec/Tests pytest suite first, including General, then the example regression suites.
   * - ``test-pytests``
     - Run the complete PDielec/Tests directory, equivalent to python -m pytest PDielec/Tests. New test directories are discovered automatically.
   * - ``test-pytest-powder_raman, test-pytest-crystal_raman, test-pytest-materials, test-pytest-calculator, test-pytest-unitcell, test-pytest-gtmcore, test-pytest-constants``
     - Run one configured pytest group.
   * - ``test-powder_ir, test-crystal_ir, test-powder_raman, test-crystal_raman, test-atr``
     - Run the corresponding example group.
   * - ``test-preader, test-p2cif, test-vibanalysis, test-pdgui``
     - Run the selected tool’s example suite.
   * - ``benchmarks``
     - Run benchmark examples without reference comparison.
   * - ``pypi``
     - Remove existing build artifacts and build source and wheel archives using python -m build. Does not upload; unavailable on Windows.
   * - ``clean``
     - Remove generated regression outputs below the root directory. Review the working tree first; unavailable on Windows.
   * - ``install``
     - Legacy Unix script-copying installation. Prefer the package installation described in the installation guide.
   * - ``scripts DIRECTORY``
     - Choose the destination for the legacy install action (default ~/bin).
   * - ``--root DIRECTORY``
     - Set the repository root used to locate tests and examples.
   * - ``--view command.pdmake``
     - Open a PDGui recipe interactively.
   * - ``--regenerate command.pdmake``
     - Replace reference outputs for the recipe; review differences first.
   * - ``--cpus N, --threads N, --threading``
     - Pass CPU count, thread count or threading selection to PDGui.
   * - ``--debug, --padding N``
     - Enable debugging or set the printed title width.
   * - ``--usesystem, -usesystem``
     - Use installed executables instead of local scripts.
   * - ``--directory, -directory``
     - Display directory names as test titles. This flag takes no value.

A recipe's first line is its title. Subsequent command lines invoke ``pdgui``,
``preader``, ``p2cif`` or ``vibanalysis`` with their arguments. Run recipes from
the directory containing their input files; for example::

    cd Examples/Powder_Raman/Castep
    pdmake command.pdmake
    pdmake --view command.pdmake

PDGui regression recipes load ``script.py`` and write ``results.xlsx``.
Use separate invocations for recipes and named suites when execution order
matters.

Excel regression comparisons include the *Settings* worksheet, including
Raman activity values and their unit headers, as well as the result sheets.
A checked sheet present in only one of the reference and generated workbooks
is reported as an error. Sheets absent from both are skipped, so infrared
examples do not require Raman result sheets. *Main* and *Scenarios* remain
excluded from the standard comparison. When running ``checkexcel`` directly,
use ``-settings`` to include Settings or ``-f`` to include all three metadata
sheets. Review reported differences before regenerating reference workbooks.


MM/QM Interfaces
================

The packages have interfaces to solid-state QM codes, VASP :cite:`Hafner2008c`, CASTEP :cite:`Clark2005d`, CRYSTAL14 :cite:`Dovesi2014`, Abinit :cite:`Gonze2016`, Quantum Espresso :cite:`Giannozzi2009`, FHI-Aims :cite:`Blum2009` and Phonopy :cite:`Togo2015`.  In addition, an interface is available for GULP :cite:`Gale2003` which is a force field based solid-state code. Finally, an interface has been written to an 'experiment' file format which allows the preparation of a user-defined file specifying the permittivities and absorption frequencies. The origin of the dataset(s) used for processing is determined by a command line switch, -program. An outline of the interfaces to these codes is given here.
The package used for the calculation is described by the -program option. In addition, a file name is given which contains the output to be processed by PDielec.

Raman input availability
------------------------

Phonon data alone do not supply Raman tensors. The following examples show
Raman-capable inputs; nonlinear susceptibility is additionally needed for an
electro-optic correction. Availability depends on what was calculated and saved.
The correction also needs the Born charges, optical permittivity, cell volume
and compatible phonon modes. A Raman activity alone cannot replace a Raman
tensor, and enabling the EO option cannot supply missing susceptibility data.
Companion files must describe the same structure, atom order and Cartesian
axes as the main calculation. See :doc:`CrystalRaman` for the EO conventions.

.. list-table:: Raman inputs and companion data
   :header-rows: 1
   :widths: 20 45 35

   * - Reader
     - Raman data
     - Nonlinear susceptibility
   * - CASTEP
     - Mode derivatives or atomic polar tensors in output; matching phonon file.
     - Read when present in output.
   * - CRYSTAL
     - ``TENS_RAMAN.DAT`` alongside the frequency output.
     - Susceptibility in output; optional ``CHI2.DAT`` takes precedence.
   * - ABINIT
     - Raman response in output; see ``Examples/Powder_Raman/AbInit``.
     - Read when present in output.
   * - Quantum ESPRESSO
     - Phonon log and/or ``tensors.xml``, with matching mode data.
     - Read from the log or XML when present.
   * - VASP
     - Optional ``Raman-Tensors.yaml`` alongside the electronic/phonon output.
     - Use a complete finite-field JSON dataset for the EO workflow.
   * - Finite field
     - JSON tensors or Cartesian dielectric derivatives with mode data.
     - Optional ``chi2`` in JSON; see the input-format reference below.

VASP
-----
The name provided on the command line is an OUTCAR file. The OUTCAR is read by PDielec to determine the unit-cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. The VASP run can be a DFPT or numerical calculation of the response.
pdgui is able to parse the OUTCAR file for the information it needs, but it is recommended to use vasprun.xml as this file gives greater precision to the required variables.
The repository includes VASP 5 and VASP 6 examples, including ``Examples/Vasp/Vasp6``. Raman calculations can supply the optional companion file ``Raman-Tensors.yaml``; see ``Examples/Powder_Raman/Vasp``.
This file contains mode Raman tensors and is required for Raman calculations
through this reader; it does not supply the nonlinear susceptibility for EO.
For a VASP finite-field calculation with EO, use the consolidated finite-field
JSON reader with both Raman data and ``chi2`` in the documented internal units.

CASTEP
-------
The name provided on the command line is the seedname for the calculation. The corresponding seedname.castep file in the current directory is read and processed to determine the unit-cell, atomic masses, optical permittivity and born charge tensors. The normal modes and their frequencies are determined from the seedname.phonon file. The CASTEP run needs to be a DFPT (phonon+efield) task.
The repository includes CASTEP 24.1 and 25.12 Raman examples. The reader supports mode polarizability derivatives and the newer atomic polar tensors, together with nonlinear susceptibility when present in the output.
Keep the matching ``.castep`` and ``.phonon`` files together. Raman tensor and
nonlinear susceptibility tables must be included in the ``.castep`` output;
no separate Raman or susceptibility companion file is needed.

CRYSTAL
-------
The name on the command line is a file ending in .out, containing the output of a CRYSTAL run. The contents of this file alone are sufficient to provide the unit-cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory containing the CRYSTAL output, it uses these files to calculate the Born charge tensors, frequencies and normal modes. The CRYSTAL calculation needs to be a frequency calculation (FREQCALC) with the infrared intensity (INTENS) selected. The default algorithm does not calculate the optical permittivity, so this needs to be provided on the command line. However, if the CPHF or CPKS algorithm is used for the frequency calculation, the optical permittivity is calculated and PDielec will automatically read it from the output file. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the PDielec package. Small differences in the calculated frequencies between the CRYSTAL program and PDielec have been observed. These have been found to be due to a slightly different method for symmetrising the 2\ :superscript:`nd` derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that PDielec should use the same symmetrisation as CRYSTAL.

For Raman calculations, place ``TENS_RAMAN.DAT`` alongside the output file.
This supplies the polarizability derivatives used to construct Raman tensors;
printed Raman activities alone are insufficient.

For electro-optic corrections, the reader recognises the
``FIRST HYPERPOLARIZABILITY (BETA) AND SECOND ELECTRIC SUSCEPTIBILITY (CHI(2))``
table in the main output. Alternatively, place that table in an optional
``CHI2.DAT`` beside the output, for example when the nonlinear response was
calculated separately. PDielec reads the ``d(MKS)`` column in pm/V and converts
it using :math:`\chi^{(2)}=2d`. A successfully read ``CHI2.DAT`` overrides the
susceptibility from the main output. It is redundant only when both sources
contain the same tensor.

``Examples/Crystal_Raman/Crystal23`` uses the susceptibility table embedded in
``opt_raman.out`` without a ``CHI2.DAT`` companion. Its reference spreadsheet
uses this embedded response for the EO correction.
``Examples/Powder_Raman/Crystal23`` retains ``CHI2.DAT`` as an example of the
optional companion format; its contents duplicate the embedded table, so it
can be omitted without changing the susceptibility.

Versions of CRYSTAL which are known to be compatible with the package are CRYSTAL14, 17, 21 & 23 later versions may be compatible also.

ABINIT
------
The output file should come from a run containing three datasets. One to calculate the wavefunction at the optimised geometry, one to calculate the field perturbations and one to calculate the second derivatives. Examples of input files and output files are available with the distribution.
The package should be compatible with AbInit 8 to 10.
For Raman and EO, use an output containing the atomic susceptibility
derivatives and the ``Non-linear optical susceptibility tensor d`` table,
as in ``Examples/Powder_Raman/AbInit/raman.abo``. The reader constructs Raman
tensors from these derivatives and the phonon modes. The example's DDB files
and ``analysis.abo`` are used for the separate ``anaddb`` analysis; they are
not automatically loaded as companions to ``raman.abo`` by PDielec.

QE
---
The output file is the dynamical matrix file, specified by "filedyn" in a run of the quantum espresso phonon package. Examples of input and output files are given in the PDielec distribution.  It is assumed that this files ends with a ".dynG" extension.
For Raman calculations, the reader supports tensors in the phonon log and
companion ``tensors.xml``. Keep matching dynamical-matrix and tensor files
together; see ``Examples/Powder_Raman/QE`` for a complete dataset.
Copy ``tensors.xml`` from the phonon calculation's ``_ph0/ph.save/`` directory
beside the input selected in PDielec. It is discovered automatically and its
Raman data take precedence over the lower-precision log tables. For EO, the
reader uses the log's ``Electro-optic tensor in cartesian axis`` block or the
XML ``ELOP_TNS`` data, preferring the XML when present. A dynamical-matrix file
alone does not supply these Raman and EO tensors; retain the tensor XML or
the matching phonon log containing them.
The latest version of QE for which the package has been tested is version 7.3.1.
pdgui is capable of reading the output log file of pwscf and ph, but it is recommended to use the pwscf.xml and the .dynG file written by ph and pwscf.

FHI-Aims
--------
The vibrational frequencies and intensities calculation uses FHI-Aims and the get_vibrations.py script which is distributed with it.
After determining the optimised cell, the script is used to generate perturbed single point calculations of the forces and the polarization which are then used to provide data for numerical differentiation to form the hessian and the IR intensities.
An example of the output from these calculations is provided in Examples/FHI-Aims/Na2SO42.
In this example, PDGui requires the hessian.na2so42.dat, masses.na2so42.dat, geometry.in and the Dielectric/ directory.
The Dielectric/ directory contains an FHI-Aims calculation of the electronic permittivity (:math:`\epsilon_{\infty}`). More details of the calculations are provided in the example directory.

The any file with the "in" or "dat" extension will cause PDGui to assume the DFT package used to create the files in that directory was the FHI-Aims package.  The current version of FHI-Aims which is supported is 250822.

PHONOPY
-------
Phonopy calculates the dynamical matrix through numerical differentiation. It has interfaces to several programs.
The second parameter for the --program directive is the PHONOPY interface that was used to calculate the forces. Typically these would be generated by performing;

        phonopy --d --dim="1 1 1"

to calculate the displacements in a set of POSCAR-* files. After running VASP a single point VASP calculation for each displacement. The FORCE\_SETS file can then be calculated using for example;

        phonopy --f DISP-\*/vasprun.xml

where the DISP-\* directories are where the VASP calculation was performed. Finally, a dynamical is written out using;

        phonopy --dim="1 1 1" --qpoints="0 0 0" --writedm

To calculate the infrared spectrum PDielec needs the Born charges for the atoms in the unit-cell and these can be calculated using the scripts provided by Phonopy.  The Phonopy scripts (for instance phonopy_vasp_born) to write a file called BORN, which contains the symmetry unique Born charges.
PDielec provides a script `phonopy-pdielec-born` which takes the BORN file and writes a file called BORN_PDIELEC containing the Born charges for all the atoms in the cell.

GULP
-----

The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit-cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit-cell. The keywords; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it.

.. _experimental-file-format:

Experimental File Format
========================
There is a way of reading in experimental or calculated permittivities and calculating the infrared optical behaviour of the material.  This uses the *experimental* file format which by default is assumed to be associated with files with *.exp* extension.
The available GUI features depend on the data supplied. Dielectric models support
optical spectra; the additional vibrational blocks described below supply normal
modes and Raman tensors for analysis, visualisation and Raman calculations.
Select *Experiment* in the Main Tab or use ``pdgui experiment material.exp``.

The general format of an experimental file is best shown by an example;::

    lattice
      5.027782
      0.866025403784439  -0.500000000000000   0.000000000000000
      0.000000000000000   1.000000000000000   0.000000000000000
      0.000000000000000   0.000000000000000   1.097684415337773
    species 2
      Si 0.0
       O 0.0
    unitcell 9
      Si 0.000000000000000   0.477290000000000   0.333333333333333
      Si 0.477290000000000   0.000000000000000   0.666666666666667
      Si 0.522710000000000   0.522710000000000   0.000000000000000
       O 0.160740000000000   0.745703000000000   0.537333670000000
       O 0.584963000000000   0.839260000000000   0.870667003333333
       O 0.415037000000000   0.254297000000000   0.795999663333333
       O 0.745703000000000   0.160740000000000   0.462666330000000
       O 0.839260000000000   0.584963000000000   0.129332996666667
       O 0.254297000000000   0.415037000000000   0.204000336666667
    epsinf
    2.296 0.0   0.0
    0.0   2.296 0.0
    0.0   0.0   2.334
    # model taken from Winta, C. J., Wolf, M., & Paarmann, A. (2019) Physical Review B, 99(14), 144308
    fpsq
    # 1.5K Data
    # Omega(TO) Gamma(TO) Omega(LO) Gamma(LO) (<2 in Table is taken to be 1.0)
    xx 6
    # E Phonon modes
      391.5     1.0        403.0    1.0
      454.0     2.6        510.5    1.0
      695.9     4.9        698.4    4.0
      797.2     4.8        810.0    4.3
     1063.7     6.1       1230.7    8.2
     1157.2     6.2       1154.9    6.1
    yy 6
      391.5     1.0        403.0    1.0
      454.0     2.6        510.5    1.0
      695.9     4.9        698.4    4.0
      797.2     4.8        810.0    4.3
     1063.7     6.1       1230.7    8.2
     1157.2     6.2       1154.9    6.1
    zz 4
    # A2 Phonon modes
      360.7     1.0        384.8    1.0
      497.9     3.1        553.6    2.8
      773.7     5.4        789.9    6.3
     1073.0     6.2       1238.7   12.4

The file starts with a definition of the lattice, the first number being a lattice constant followed by three lines specifying the a, b, and c directions of the unit-cell.

The next directive specifies the species and their masses, by default PDGui uses its own internal mass scheme, so unless otherwise required these masses will be overwritten.

The *unitcell* directive specifies the number of atoms in the unit-cell and their fractional coordinates.

The *epsinf* directive gives the values of the :math:`\epsilon_{\infty}` tensor

Finally, a model for the specification of the frequency-dependent permittivity is given.  In this case the four parameter semi-quantum model (FPSQ) is employed.  Each diagonal component of the tensor is specified by the number of contributions and then each contribution provides the TO frequency and its :math:`\gamma` followed by the same for the LO frequency (all frequencies are in |cm-1|. Full details of each permittivity model available are shown below.


Constant model (constant)
-------------------------

The constant model defines a frequency-independent permittivity.  The data for such a model is shown below.::

    constant
    2.0+0.1j 0.0      0.0
    0.0      2.0+0.1j 0.0
    0.0      0.0      2.0+0.1j

This would specify an isotropic permittivity with some absorption.




FPSQ model (fpsq)
-----------------

The FPSQ model defines a frequency-dependent permittivity using the Four Parameter Semi-Quantum model.  An example of the data for such a model is shown below.::

    fpsq
    xx 6
    # E Phonon modes
    # Omega(TO) Gamma(TO) Omega(LO) Gamma(LO)
      391.5     1.0        403.0    1.0
      454.0     2.6        510.5    1.0
      695.9     4.9        698.4    4.0
      797.2     4.8        810.0    4.3
     1063.7     6.1       1230.7    8.2
     1157.2     6.2       1154.9    6.1
    yy 6
    # E Phonon modes
    # Omega(TO) Gamma(TO) Omega(LO) Gamma(LO)
      391.5     1.0        403.0    1.0
      454.0     2.6        510.5    1.0
      695.9     4.9        698.4    4.0
      797.2     4.8        810.0    4.3
     1063.7     6.1       1230.7    8.2
     1157.2     6.2       1154.9    6.1
    zz 4
    # A2 Phonon modes
    # Omega(TO) Gamma(TO) Omega(LO) Gamma(LO)
      360.7     1.0        384.8    1.0
      497.9     3.1        553.6    2.8
      773.7     5.4        789.9    6.3
     1073.0     6.2       1238.7   12.4

The model only allows for a diagonal permittivity tensor and each component of the tensor specified requires the number of terms in the expansion to the specified.  Each component of the permittivity tensor is generated using the following formula;

.. _fpsq:
.. math::

    \epsilon (\omega )=\epsilon _{\infty}\prod_{j} \frac{\Omega^2_{LO_j}-\omega ^2-i\gamma _{LO_j}\omega }{\Omega^2_{TO_j}-\omega ^2-i\gamma _{TO_j}\omega}


Drude-Lorentz model (drude-lorentz)
-----------------------------------

The drude-lorentz model defines a frequency-dependent permittivity.  An example of the data for such a model describing MgO is shown below.::

    drude-lorentz
    xx 2
    # T Phonon modes
    # Frequency    Strength  Gamma
      413.7           1050.0  22.2
      652.2             90.0  64.6
    yy 2
      413.7           1050.0  22.2
      652.2             90.0  64.6
    zz 2
      413.7           1050.0  22.2
      652.2             90.0  64.6

The model only allows for a diagonal permittivity tensor and each component of the tensor specified requires the number of terms in the expansion to the specified.  Each component of the permittivity tensor is generated using the following formula;

.. _drude:
.. math::

    \epsilon (\omega )=\epsilon _{\infty} +  \sum_{j} \frac{ S_j^2 }{ \Omega^2_j-\omega^2 -i \sigma_j \omega }

Where :math:`S_j` is the strength of the phonon mode, :math:`\Omega_j` its frequency and :math:`\sigma_j` the width of the mode.



.. include:: input_formats.inc


The Format of the Materials' Database
=====================================

The materials' database stores some basic information about the supporting materials often used in THz and Infrared spectroscopy.
The file *PDielec/MaterialsDatabase.xlsx* can be opened and inspected to see the format of the entries.
The first sheet gives some information about the structure of the spreadsheet.
Each tab, except the first, represents a material and the tab name can be referred to in the *Support matrix* pulldown of the *PowderScenarioTab*.

In Columns G&H of a material sheet are, in arbitrary order: an entry, a density, a temperature; low and high frequency values, the number of points in a frequency-dependent entry and a reference.
Column G holds the directive and column H holds its value
For the case that the unit-cell of the material is required (for anisotropic materials for instance), the unit-cell can be provided either in terms of a, b, c, :math:`\alpha`, :math:`\beta` and :math:`\gamma` or as cell vectors.

The entry type may be one of: "Constant refractive index", "Constant permittivity", "Tabulated refractive index" or "Tabulated permittivity", "Drude-Lorentz" or "FPSQ".

For columns A to F the first row provides a label for the data in each column.  The label is not used by PDGui but can be used to describe the contents of the column.
The contents of these columns depend upon the entry type.

For tabulated entries, column A must have the frequency in |cm-1| and in ascending order.
Column B is not used by PDGui but because refractive indices are often tabulated with micron wavelengths, this column may be used to allow conversion to |cm-1|.

   Summary of columns used for different Entry modes

   +----------------------------+-------------+-----------------------------+-----------------------+-----------------------+-----------------------+-----------------------+
   | Entry                      | A:          |  B:                         |   C:                  |  D:                   |  E:                   |  F:                   |
   +============================+=============+=============================+===========+===========+=======================+===========+===========+=======================+
   | Constant permittivity      |             |                             | :math:`{\epsilon}_r`  | :math:`{\epsilon}_i`  |                       |                       |
   +----------------------------+-------------+-----------------------------+-----------------------+-----------------------+-----------------------+-----------------------+
   | Constant refractive index  |             |                             | :math:`\nu`           | :math:`\kappa`        |                       |                       |
   +----------------------------+-------------+-----------------------------+-----------------------+-----------------------+-----------------------+-----------------------+
   | Tabulated permittivity     | v( |cm-1| ) |                             | :math:`{\epsilon}_r`  | :math:`{\epsilon}_r`  |                       |                       |
   +----------------------------+-------------+-----------------------------+-----------------------+-----------------------+-----------------------+-----------------------+
   | Tabulated refractive index | v( |cm-1| ) |                             |  :math:`\nu`          |  :math:`\kappa`       |                       |                       |
   +----------------------------+-------------+-----------------------------+-----------------------+-----------------------+-----------------------+-----------------------+
   | FPSQ                       | xx/yy/zz    | :math:`{\epsilon}_{\infty}` | :math:`{\Omega}_{TO}` | :math:`{\gamma}_{TO}` | :math:`{\Omega}_{LO}` | :math:`{\gamma}_{LO}` |
   +----------------------------+-------------+-----------------------------+-----------------------+-----------------------+-----------------------+-----------------------+
   | Drude-Lorentz              | xx/yy/zz    | :math:`{\epsilon}_{\infty}` | :math:`{\Omega}_j`    | :math:`{S_j}`         | :math:`{\sigma_j}`    |                       |
   +----------------------------+-------------+-----------------------------+-----------------------+-----------------------+-----------------------+-----------------------+


   Summary of columns G and H

   +-----------------------+----------------------------------------+
   | Directive (G)         | Contents (H)                           |
   +=======================+========================================+
   | | Entry:              | One of;                                |
   | |                     |  - 'tablulated refractive index'       |
   | |                     |  - 'tablulated permittivity'           |
   | |                     |  - 'constant refractive index'         |
   | |                     |  - 'constant refractive index'         |
   | |                     |  - 'Drude-Lorentz'                     |
   | |                     |  - 'FPSQ'                              |
   +-----------------------+----------------------------------------+
   | Density:              | Density in g/ml                        |
   +-----------------------+----------------------------------------+
   | Temperature:          | Temperature of measurement             |
   +-----------------------+----------------------------------------+
   | Low:                  | The lowest valid frequency             |
   +-----------------------+----------------------------------------+
   | High:                 | The highest valid frequency            |
   +-----------------------+----------------------------------------+
   | NPoints:              | The number of points in a tabulation   |
   +-----------------------+----------------------------------------+
   | Reference:            | The reference for the data             |
   +-----------------------+----------------------------------------+
   | a_vector              | 3 numbers specifying the a-direction   |
   +-----------------------+----------------------------------------+
   | b_vector              | 3 numbers specifying the b-direction   |
   +-----------------------+----------------------------------------+
   | c_vector              | 3 numbers specifying the c-direction   |
   +-----------------------+----------------------------------------+
   | a:                    | The length of a-direction              |
   +-----------------------+----------------------------------------+
   | b:                    | The length of b-direction              |
   +-----------------------+----------------------------------------+
   | c:                    | The length of c-direction              |
   +-----------------------+----------------------------------------+
   | alpha:                | Cell angle alpha                       |
   +-----------------------+----------------------------------------+
   | beta:                 | Cell angle beta                        |
   +-----------------------+----------------------------------------+
   | gamma:                | Cell angle gamma                       |
   +-----------------------+----------------------------------------+


The Format of the Output Spreadsheet
====================================

The Excel spreadsheet, which can be written by the program, contains details of the system being analysed, a list of the scenario settings and tables of absorption, permittivity, reflection and transmission. It also includes mode frequencies, IR intensities and Raman activities when Raman data are available.
Result sheets are conditional on the scenario types present: powder IR/ATR,
crystal IR, powder Raman and crystal Raman each supply their corresponding
sheets. Spectrum frequencies are exported in cm⁻¹. Plotting-tab spectrum
renormalisation is a display operation and is not applied by the spreadsheet
exporter. The Settings sheet contains base mode data, rather than the optional
NAC/EO diagnostics displayed in the frequency table.

The spreadsheet is divided into different sheets;

-   **Main**
        | This stores the information concerning the file that is being analysed and the program associated with the output.  A list of the frequencies read in from the output file is given too.

-   **Settings**
        | The information here comes from the Settings Tab: masses, permittivities and the frequencies that will be used in subsequent calculations. These frequencies can include corrections to the dynamical matrix to project out translational modes, corrections to the Born charges to ensure they sum to zero and changes to the masses of the atoms. Each mode has an inclusion flag and a linewidth (Sigma), followed by its frequency, IR intensity, integrated molar absorption and absorption maximum. The Raman columns give total, parallel and perpendicular activities, labelled ``Raman R_total``, ``Raman R_parallel`` and ``Raman R_perp``. The activity values use the units selected by *Raman activity units* in the Settings Tab: Å⁴/amu for polarizability-volume activities or Å/amu for the internal convention. The column headers indicate the selected units. These activities describe the individual modes, rather than the Raman spectrum calculated for a particular scenario.

-   **Analysis**
        | The analysis of the vibrational modes into molecular, internal and external contributions is summarised here.

-   **Scenarios**
        | A complete list of the settings for every scenario used in the calculation is given.

-   **Powder Molar Absorption (cells)**
        | The molar absorption in moles of unit-cell is given as a function of frequency for each powder scenario.  The units are |Lmol-1cm-1|.

-   **Powder Absorption**
        | The absorption for each powder scenario is tabulated as a function frequency.  The units are |cm-1|

-   **Powder Real Permittivity**
        | The real component of the complex permittivity for each powder scenario is tabulated as a function frequency.

-   **Powder Imaginary Permittivity**
        | The imaginary component of the complex permittivity for each powder scenario is tabulated as a function frequency.

-   **Powder ATR Reflectance**
        | The extinction coefficient associated with the attenuated total reflectance for each powder scenario is tabulated as a function frequency.  See :ref:`ATR_theory` for more details.

-   **Crystal R_p**
        | The P polarised component of the reflectance for each single crystal scenario is tabulated as a function of frequency.

-   **Crystal R_s**
        | The S polarised component of the reflectance for each single crystal scenario is tabulated as a function of frequency.

-   **Crystal T_p**
        | The P polarised component of the transmittance for each single crystal scenario is tabulated as a function of frequency.

-   **Crystal T_s**
        | The S polarised component of the transmittance for each single crystal scenario is tabulated as a function of frequency.

-   **Real Crystal Permittivity**
        | The real components of the crystal permittivity tensor are tabulated as a function of frequency.

-   **Imag Crystal Permittivity**
        | The imaginary components of the crystal permittivity tensor are tabulated as a function of frequency.


-   **Powder Molar Absorption (mols)** and **Powder Molar Absorption (atoms)**
        | The selected molecular or atomic molar definition adds the corresponding sheet, in |Lmol-1cm-1|, alongside the unit-cell sheet.

-   **Crystal A_p** and **Crystal A_s**
        | Dimensionless p- and s-polarised absorptance for crystal infrared scenarios, alongside reflectance and transmittance.

-   **Powder Raman** and **Crystal Raman**
        | Raman spectra for the corresponding scenarios, in arbitrary intensity units, as a function of Raman shift in cm⁻¹. These include scenario-dependent geometry and line broadening; they are distinct from the Settings-sheet mode activities in Å⁴/amu or Å/amu.
