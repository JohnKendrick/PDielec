.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


=========================
Other Software Components
=========================

PReader
=======

preader is a 'helper' program which uses the underlying modules of PDielec to read output files and summarise the results of various MM/QM packages.  The program can be used to perform some straightforward calculations.  For instance projection of any remaining centre of mass motion of the crystal can be performed to make sure that there are three zero frequencies.  Also the masses used in the calculation of the dynamical matrix can be altered.
Unlike PDGui it is not necessary to have performed a full calculation of the dynamical matrix.  In the majority of cases preader will read geometry optimisation runs.

Command options
---------------

Examples of data sets for these packages are included with the distribution and can be found in the Examples/'Package'/preader directory. The program is run from the command line. There are several command options and these are summarised below. Some options may be repeated.

  - \-program program_name
       | Program_name can be “abinit”,  “castep”, “crystal”, “gulp”, “qe”, “experiment” or “vasp” and specifies the program which generated the results to be analysed
  - \-neutral
       | Impose neutrality on the Born charge matrices
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

::

         preader -program vasp `find . -name OUTCAR` > results.csv

This reads all the VASP OUTCAR files in the current and any of its subdirectories and summarises the results to results.csv. ::

         preader -program castep -eckart `find . -name \*.castep` > results.csv

This reads all the castep output files in the current and any of its subdirectories and summarises the results to results.csv.  For each file the centre of mass motion of the crystal is projected.  The results file contains both the unprojected and the projected results.

There are examples of preader being used in the Examples/'Package'/preader subdirectories of the distribution of PDielec.

Vibanalysis
===========

The vibanalysis command provides an interface the vibAnalysis package of Philipe Teixeira.  This package provides Vibrational Mode Decomposition of the phonon modes using a variety of methods including; Vibrational Mode Projection (VMP), Vibrational Mode Linear Decomposition (VMLD) and Vibrational Mode Automatic Relevance Determination (VMARD).  By default VMARD is used.

The package options can be determined using; 

::

       vibanalysis OUTCAR
       vibanalysis phonon.castep

which will analyse a VASP and a CASTEP phonon calculation output to produce a summary in the files OUTCAR.nma and phonon.nma respectively.

pdmake
======
The testing suite is run use the the pdmake command.  The command line options are summarised below.

   +-----------------------+----------------------------------------+
   | Command line option   | Description                            |
   +=======================+========================================+
   | test                  | Run all the tests                      |
   | tests                 | Run all the tests                      |
   +-----------------------+----------------------------------------+
   | test-preader          | Run preader tests                      |
   +-----------------------+----------------------------------------+
   | test-pdgui            | Run pdgui tests                        |
   +-----------------------+----------------------------------------+
   | test-p2cif            | Run p2cif tests                        |
   +-----------------------+----------------------------------------+
   | test-vibanalysis      | Run vibanalysis tests                  |
   +-----------------------+----------------------------------------+
   | benchmarks            | Run the benchmark suite (no testing)   |
   +-----------------------+----------------------------------------+
   | --regenerate          | Regenerate the reference data          |
   +-----------------------+----------------------------------------+
   | --root adirectory     | Set the root directory to adirectory   |
   +-----------------------+----------------------------------------+
   | --directory           | Print the directory name for the test  |
   +-----------------------+----------------------------------------+
   | --debug               | Switch on debugging                    |
   +-----------------------+----------------------------------------+
   | --view                | Use windowing with pdgui               |
   +-----------------------+----------------------------------------+
   | --padding 60          | Set the size of the title field to 60  |
   +-----------------------+----------------------------------------+
   | --usesystem           | Force the use of system executables    |
   +-----------------------+----------------------------------------+


The command can also be use in an examples directory.  There will be a file, command.pdmake, which contains the instructions to run the test.
The first line of the file is a title that is used to identify the test.
The second line is the command to be used in the test followed by the command parameters required.

Examples of the use of the command are;::

        pdmake                              # prints a help summary
        pdmake tests                        # runs all the tests
        pdmake command.pdmake               # runs the test in command.pdmake
        pdmake --debug command.pdmake       # run the test with debugging
        pdmake --view command.pdmake        # Bring up the GUI interface




MM/QM Interfaces
================

The packages have interfaces to five solid state QM codes, VASP :cite:`Hafner2008c`, CASTEP :cite:`Clark2005d`, CRYSTAL14 :cite:`Dovesi2014`, Abinit :cite:`Gonze2016`, Quantum Espresso :cite:`Giannozzi2009` and Phonopy :cite:`Togo2015`.  In addition an interface is available for GULP :cite:`Gale2003` which is a force field based solid state code. Finally an interface has been written to an 'experiment' file format which allows the preparation of a user defined file specifying the permittivities and absorption frequencies. The origin of the dataset(s) used for processing is determined by a command line switch, -program. An outline of the interfaces to these codes is given here. 
The package used for the calculation is described by the -program option. In addition, a file name is given which contains the output to be processed by PDielec.

VASP 
-----
The name provided on the command line is an OUTCAR file. The OUTCAR is read by PDielec to determine the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. The VASP run can be a DFPT or numerical calculation of the response.

CASTEP 
-------
The name provided on the command line is the seedname for the calculation. The corresponding seedname.castep file in the current directory is read and processed to determine the unit cell, atomic masses, optical permittivity and born charge tensors. The normal modes and their frequencies are determined from the seedname.phonon file. The CASTEP run needs to be a DFPT (phonon+efield) task.

CRYSTAL 
-------
The name on the command line is a file ending in .out, containing the output of a CRYSTAL14 run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory where PDielec is run from, it uses these files to calculate the Born charge tensors, frequencies and normal modes. The CRYSTAL calculation needs to be a frequency calculation (FREQCALC) with the infrared intensity (INTENS) selected. The default algorithm does not calculate the optical permittivity, so this needs to be provided on the command line. However, if the CPHF or CPKS algorithm is used for the frequency calculation, the optical permittivity is calculated and PDielec will automatically read it from the output file. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the PDielec package. Small differences in the calculated frequencies between the CRYSTAL program and PDielec have been observed. These have been found to be due to a slightly different method for symmetrising the 2\ :superscript:`nd` derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that PDielec should use the same symmetrisation as CRYSTAL14.

ABINIT
------
The output file should come from a run containing three datasets. One to calculate the wavefunction at the optimised geometry, one to calculate the field perturbations and one to calculate the second derivatives. Examples of input files and output files are available with the distribution.

QE 
---
The output file is the dynamical matrix file, specified by "filedyn" in a run of the quantum espresso phonon package. Examples of input and output files are given in the PDielec distribution

PHONOPY 
-------
Phonopy calculates the dynamical matrix through numerical differentiation. It has interfaces to several programs, although PDielec has only used the VASP interface. In principle other interfaces could be used. The second parameter for the --program directive is the PHONOPY interface that was used to calculate the forces. Typically these would be generated by performing;

        phonopy --d --dim="1 1 1"

to calculate the displacements in a set of POSCAR-* files. After running VASP a single point VASP calculation for each displacement. The FORCE\_SETS file can then be calculated using for example;

        phonopy --f DISP-\*/vasprun.xml

where the DISP-\* directories are where the VASP calculation was performed. Finally a dynamical is written out using;

        phonopy --dim="1 1 1" --qpoints="0 0 0" --writedm

To calculate the infrared spectrum PDielec needs the Born charges for the atoms in the unit cell and these can be calculated using VASP and the optimised geometry of the unit cell. The OUTCAR file from this calculation can be copied to the current directory and renamed OUTCAR.born

GULP 
-----

The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit cell. The key words; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it.

.. _experimental-file-format:

Experimental File Format
========================
There is a way of reading in experimental or calculated permittivities and calculating the infrared optical behaviour of the material.  The uses the *experimental* file format which by default is assumed to be associated with files with *.exp* extention. 
There are some significant limitations in providing the information in the manner.  Specifically as there is no Dynamical Matrix the normal modes are not known, so they will not be visible in the *Viewer tab*.   There may be little information about specific modes for the *Settings tab*.

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

Finally a model for the specification of the frequency dependent permittivity is given.  In this case the four parameter semi-quantum model (fpsq) is employed.  Each diagonal component of the tensor is specified by the number of contributions and then each contribution provides the TO frequency and its :math:`\gamma` followed by the same for the LO frequency (all frequencies are in cm-1. Full details of each permittivity model available are shown below.


Constant model (constant)
-------------------------

The constant model defines a frequency independent permittivity.  The data for such a model is shown below.::

    constant
    2.0+0.1j 0.0      0.0
    0.0      2.0+0.1j 0.0
    0.0      0.0      2.0+0.1j

This would specify an isotropic permittivity with some absorption.
    



FPSQ model (fpsq)
-----------------

The fpsq model defines a frequency dependent permittivity using the Four Parameter Semi-Quantum model.  An example of the data for such a model is shown below.::

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

.. math::
   :label: eq-fpsqb

    \epsilon (\omega )=\epsilon _{\infty}\prod_{j} \frac{\Omega^2_{LO_j}-\omega ^2-i\gamma _{LO_j}\omega }{\Omega^2_{TO_j}-\omega ^2-i\gamma _{TO_j}\omega}


Drude-Lorentz model (drude-lorentz)
-----------------------------------

The drude-lorentz model defines a frequency dependent permittivity.  An example of the data for such a model describing MgO is shown below.::

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

.. math::
   :label: eq-drude

    \epsilon (\omega )=\epsilon _{\infty}\sum{j} \frac{\Omega^2_{LO_j}-\omega ^2-i\gamma _{LO_j}\omega }{\Omega^2_{TO_j}-\omega ^2-i\gamma _{TO_j}\omega}



The Format of the Materials' Database
=====================================

The materials' database stores some basic information about the supporting materials often used in THz and Infrared spectroscopy.
The file *PDielec/MaterialsDatabase.xlsx* can be opened and inspected to see the format of the entries.
The first sheet gives some information about the structure of the spreadsheet.
Each tab, except the first, represents a material and the tab name can be referred to in the *Support matrix* pulldown of the *PowderScenarioTab*.

In Columns G&H of a material sheet are in arbitrary order: an entry, a density, a temperature; low and high frequency values, the number of points in a frequency dependent entry and a reference.
Column G holds the directive and column H holds its value
For the case that the unit cell of the material is required (for anisotropic materials for instance), the unit cell can be provided either in terms of a, b, c, :math:`\alpha`, :math:`\beta` and :math:`\gamma` or as cell vectors.

The entry type maybe one of: "Constant refractive index", "Constant permittivity", Tabulated refractive index' or "Tabulated permittivity", "Drude-Lorentz" or "FPSQ".

The density gives the material density in g/ml.

For columns A to F the first row provides a label for the data in each column.  The label is not used by PDGui.
The contents of these columns depends upon the entry type.

For Tabulated entries column A must have the frequency in |cm-1| and in ascending order.  
Column B is not used by PDGui but because refractive indices are often tabulated with micron wavelengths, this column may be used to allow conversion to |cm1|.

   Summary of columns used for different Entry modes

   +----------------------------+---------+-------+-----------+-----------+-----------+-----------+
   | Entry                      | A:      |  B:   |   C:      |  D:       |  E:       |  F:       |
   +============================+=========+=======+===========+===========+===========+===========+
   | Constant permittivity      |         |       | eps(real) | eps(imag) |           |           |
   +----------------------------+---------+-------+-----------+-----------+-----------+-----------+
   | Constant refractive index  |         |       | ri(real)  | ri(imag)  |           |           |
   +----------------------------+---------+-------+-----------+-----------+-----------+-----------+
   | Tabulated permittivity     | v(cm-1) |       | eps(real) | eps(imag) |           |           |
   +----------------------------+---------+-------+-----------+-----------+-----------+-----------+
   | Tabulated refractive index | v(cm-1) |       |  ri(real) |  ri(imag) |           |           |
   +----------------------------+---------+-------+-----------+-----------+-----------+-----------+
   | FPSQ                       | xx/yy/zz| epsinf|  Omega(TO)| Gamma(TO) | Omega(LO) | Gamma(LO) |
   +----------------------------+---------+-------+-----------+-----------+-----------+-----------+
   | Drude-Lorentz              | xx/yy/zz| epsinf|  Omega    | Strength  | Gamma     |           |
   +----------------------------+---------+-------+-----------+-----------+-----------+-----------+
  
 
   Summary of columns G and H

   +-----------------------+----------------------------------------+
   | Directive (G)         | Contents (H)                           |
   +=======================+========================================+
   | Entry:                | One of;                                |
   |                       | - 'tablulated refractive index'        |
   |                       | - 'tablulated permittivity'            |
   |                       | - 'constant refractive index'          |
   |                       | - 'constant refractive index'          |
   |                       | - 'Drude-Lorentz'                      |
   |                       | - 'FPSQ'                               |
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

The Excel spread sheet, which can be written by the program, contains details of the system being anaylysed, a list of the scenario settings and tables of absorption, permittivity, reflection and transmission.
The spreedsheet is divided into different sheets;

-   **Main**
        | This stores the information concerning the file that is being analysed and the program associated with the output.  A list of the frequencies read in from the output file is given too.

-   **Settings**
        | The information here comes from the settings tab.  Masses, permittivities and the frequencies which will be used in subsequent calculations.  These frequencies can include corrections to the dynamical matrix to project out translational modes, corrections to the Born charges to ensure they sum to zero and changes to the masses of the atoms.  Along with the frequencies the intensities, integrated absorption and absorption maximum are given.
  
-   **Analysis**
        | The analysis of the vibrational modes into molecular, internal and external contributions is summarised here.

-   **Scenarios**
        | A complete list of the settings for every scenario used in the calculation is given.

-   **Powder Molar Absorption (cells)**
        | The molar absorption in moles of cell is given as a function of frequency for each powder scenario.  The units are |Lmol-1cm-1|.

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



