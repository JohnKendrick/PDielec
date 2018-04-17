.. include:: preamble.txt

.. contents::
   :local:


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


===================
Software components
===================

There are three programs within the PDielec package.  *Pdielec* and *preader* are command line programs and *PDGui* is a graphical interface to the features of the PDielec package.


PDGui
=====


pdgui is provides a graphical front end to the capabilities of PDielec.
pdgui can be run from the command line without any parameters or the QM/MM program name followed by the file to be read can be provided.
For example; ::

    pdgui

or ::

    pdgui vasp OUTCAR

If a spreadsheet of results is required the name of the spreadsheet can be provided on the command line. ::

    pdgui vasp OUTCAR results.xlsx
    
The spreadsheet name must have an .xlsx extension.

The package requires the PyQt5 library.
After running the program the user sees a notebook interface with six tabs.

    *Main Tab* allows the user to specify the program and filename which are to be anaylsed.  It is also possible to specify here the spreadsheet name of results are going to be saved.

    *Settings Tab* is used for changing the settings within the package

    *Scenario Tab* specifies the material properties used to define the effective medium from which the  absorption characteristic are calculated.  This tab can also be used to change the algorithms used to calculate the effective medium.

    *Plotting Tab* shows the absorption or permittivity as a function of frequency.

    *Analysis Tab* shows the decomposition of the normal modes into molecular components.

    *Visualisation Tab* displays the normal modes as either arrows showing the atomic displacement or as an animation.


Main tab
--------

The *Main Tab* is used to pick the MM/QM package and the output file which will be analysed.

.. _fig-mainTab:

.. figure:: ./_static/Figures/MainTab.png 
   :scale: 80%

   The Main Tab

The QM/MM program can be chosen from the dropdown list.  The output file name can be input into the text box below it.  Any entry in this *Output file name* text box will be read when the 'Read the output file and start the calculation' button is pressed.  If the file is not valid a file chooser will pop-up and the user can select a file from that.

If an Excel spreadsheet file name given, it must have the extensions .xlsx.  The spreadsheet is written when the program exits.

Once the file has been specified and read the frequencies found in the calculation file are reported in the output text box.


Settings Tab
------------

The *Settings Tab* affects the calculation of the frequencies and their intensities.

.. _fig-settingsTab:

.. figure:: ./_static/Figures/SettingsTab.png 
   :scale: 80%

   The Settings Tab

The Eckart conditions control whether 3 modes at zero frequency due to translational invariance are projected out of the dynamical matrix.  By default the Eckart conditions are applied and the 3 modes are projected out.

There is an option to set the sum of the Born charge matrices to zero.  By default this is not applied.  

The atomic masses can be specified in a variety of ways.  They can be chosen from;
  - The average natural abundance
  - The masses used by the QM/MM program
  - The isotopic mass of the most abundant isotope

If necessary the mass of each element can be edited separately by clicking on the mass concerned and entering a new number.

The width (:math:`\sigma`) of all the absorptions can be set using the *Lorentzian width(sigma)*  spin box. 

Finally the optical permittivity at zero frequency is given in the *Optical permittivity* table.  In some cases it is necessary to enter the optical permittivity by hand.  This can be done by clicking each element in the table which needs changing and typing the new matrix element.

The output table at the bottom of the tab shows the calculated frequencies and their intensities.  Transitions which do not contribute to the Infrared absorption are greyed out.  These transitions will not be used in later calculations.  From this table it is possible to remove or add transitions to the later calculations and it is also possible to change the width of individual transitions.

Scenario Tabs
-------------

There can be more than one *Scenario Tab*.  Each one specifies a particular material, method or particle shape which will be used for the calculation of the effective medium.  The results of the calculations will be shown in the *Plotting Tab*.

.. _fig-scenarioTab:

.. figure:: ./_static/Figures/ScenarioTab.png 
   :scale: 80%

   The Scenario Tab

The support matrix into which the active dielectric material is dispersed can be selected from the *Support matrix* drop down menu.  The selected supporting material will change the density and permittivity shown in the respective text boxes.  The user can edit these values independently if necessary.  

The amount of dielectric material to be considered can be entered either as mass fraction (in percent) or as a volume fraction (in percent).  If the matrix support density is changed the calculated mass fraction will be updated.  It is assumed that the volume fraction has precedence.

The calculation of the effective medium can be performed using a variety of methods which can be chosen from the *Method* drop down menu.  If the *Mie* method is chosen the user can enter the particles radius (microns).  The *Particle sigma* specifies the width of a log-normal distribution.  If the width is 0.0 no sampling of the distribution is performed.

For effective medium theories other than the Mie method the particle shape can be specfied using the *Shape* pull down menu.  Possible shapes are *Sphere*, *Needle*, *Plate* and *Ellipsoid*.  For the cases of *Needle* and *Ellipsoid* the unique direction is specifed by a direction \[abc\] in lattice units.  In the case of *Plate* the unique direction is specifies as the normal to a plane (hkl) in reciprical lattice units.  If an *Ellipsood* shape is used the eccentricy factor can be specified in the *Ellipsoid a/b eccentricity* text box.

The plot in the *Plotting Tab* has a legend and the description used in the legend for each scenario can be set by filling out the text box associated with the *Scenario legend*.

Scenarios can be added or removed using the push buttons at the bottom of each *Scenario Tab*.  When a new scenario is created the settings are copied from the current scenario.  


Plotting Tab
------------

The *Plotting tab* controls and plots the absorption and permittivity as a function of frequency calculated for each scenario present in the pdgui notebook.


.. _fig-plottingTab:

.. figure:: ./_static/Figures/PlottingTab.png 
   :scale: 80%

   The Plotting Tab

The minimum and maximum frequencies can be specified along with the frequency increment.  An effective medium theory calculation will be performed for each scenario at every frequency between the minimum and maxmimum frequencies at the interval specified.

By default the calculation uses moles of unit cells to calculation the molar absorption coefficient.  This can be altered using the the 'Molar definition' pull down menu.  Options include *Unit cells*, *Atoms* and *Molecules*.  In the case of *Molecules* it is necessary to supply the number of atoms in a formula unit of the compound in the *Number of atoms per molecule* text box.

The title of the plot can be supplied by entering it into the *Plot title* text box and the frequency units used for the plot can also be changed from *wavenumber* to *THz*.

The molar absorption, the absorption and the real or imaginary permittivity can be plotted.  Once a plot has been requested the calculation progress is shown in the progress bar.  Some settings can be changed without the whole plot being recalculated.

Analysis tab
------------

The *Analysis tab* shows a breakdown of the phonon modes into molecular components.  The molecular structure of the unit cell is determined by the covalent radii of the atoms.  These can be speficied individually if needed.  The *Analysis tab* shows the number of molecules that have been found in the analysis.


.. _fig-analysisTab:

.. figure:: ./_static/Figures/AnalysisTab.png 
   :scale: 80%

   The Analysis Tab

The bar graph shows a break down of each normal mode in the chosen frequency range into either internal and external contributions or into molecular components.  The atom sizes and the molecular composition of the unit cell is displayed in the *3D viewer tab*

3D Viewer tab
-------------

The *3d viewer tab* shows the unit cell of the system using the molecular information and atomic sizes from the *analysis tab*.

The atomic displacement of each phonon can either be shown as arrows or as an animation.  The views and animations can be recorded in  .png and .mp4 files respectively.  If a .gif file is specified the animation is recorded but with reduced numbers of colours.


.. _fig-viewerTab:

.. figure:: ./_static/Figures/ViewerTab.png 
   :scale: 80%

   The 3D Viewer Tab

As well as being able to change the phonon mode being analysed.  The colours and many settings in the visualiser can be adjusted from the settings tab.


PDielec
=======

Pdielec is a command line interface to the capabilities of the PDielec package.

Options
-------

The program options are summarised below, a ✔ indicates that the option can be specified more than once.  Where there is a default value, it value is shown.

  -program program
    Program can be “abinit”,  “castep”, “crystal”, “gulp”, “qe”, “experiment” or “vasp” and specifies the program which generated the results to be analysed.
  -method maxwell (✔)
    The method is given by the string and is either ‘ap’, ‘maxwell’, ‘bruggeman’ or ‘mie’.
  -sphere (✔)
    The inclusion is a sphere, the default if no other shape is given.
  -needle h k l (✔)
    The inclusion is a needle whose unique directionis given by the direction \[hkl\].
  -plate h k l (✔)
    The inclusion is a plate whose  surface is defined by the Miller indices \(hkl\).   Note that needles and ellipsoid use  directions in crystal coordinates defined by \[hkl\].   For non-orthogonal lattices the normal to  the \(hkl\) is not necessarily the same as \[hkl\].
  -ellipse h k l z (✔)
    The inclusion is an ellipsoid, whose  unique direction is given by \[hkl\],  z  specifies the eccentricity of the ellipsoid.
  -vf 0.1 (✔)
    Specifies the volume fraction of inclusion.
  -mf 0.0 (✔)
    Specifies a mass fraction from  which the volume fraction is calculated.   The calculation requires the density of the supporting matrix.
  -matrix ptfe 
    The supporting matrix is defined by  the string.  Options are “ptfe”, “kbr”,  “ujol”, “air”, “vacuum”, “ldpe”, “mdpe”, “hdpe”.  If the matrix is given in this way both the  density and the permittivity of the supporting matrix are defined.  Alternatively the density and dielectric  options can be used.
  -density 2.2
    Defines the density of the  supporting matrix.
  -dielectric 2.0
    Defines the dielectric of the  supporting matrix.
  -LO h k l (✔)
    The frequencies corresponding to the longitudinal optic modes with a k vector direction (h k l) are calculated  using Equations :eq:`eq-LODynamicalMatrix`.
  -LO_cart x y z (✔)
    As above but for Cartesian  directions
  -sigma 5.0 
    Specifies the damping factor, σ,  for all modes in |cm-1|, as used in Equation :eq:`eq-permittivity`.
  -mode_sigma k σ (✔)
    The k’th mode is assigned a specific σ (|cm-1|).
  -vmin 0.0
    Thee starting wavenumber (|cm-1|)  for the frequency range.
  -vmax 300.0
    The final wavenumber (|cm-1|)  for the frequency range.
  -i 0.2
    The increment used to cover  the frequency range (|cm-1|).
  -plot absorption (✔)
    Plot types  are specified by the string and they can be  ‘absorption’, ‘molar_absorption’, ‘real’ or ‘imaginary’.
  -excel filename.xlsx
    Writes the results to an excel  spread sheet with the specified name.
  -csv filename.csv
    Output is sent to the specified comma separated  file.
  -csv_ext filename
    Output is sent to 3 comma separated  files; filename_command.csv, filename_frequency.csv and filename_spectrum.csv.
  \-print
    Additional output is provided from  the QM or MM calculation.
  -ignore k (✔)
    Ignore the kth mode (any mode below 5 |cm-1| is ignored automatically).
  -mode k (✔)
    Only use the kth mode in the  calculation of the permittivity.
  -threshold 1.0e-10 5.0
    The modes selected for inclusion in  the absorption calculation have to have an IR intensity greater than 1.0E-10 and a frequency greater than 5 |cm-1|.
  \-eckart
    The translational modes are  projected out of the hessian before diagonalisation.
  \-hessian symm
    This option can specify either “crystal” or  “symm”.  In the case of “crystal” the  hessian is symmetrised using the same algorithm as Crystal14.
  -optical z1 z2 z3
    z1,z2 and z3 define the diagonal of  the optical permittivity tensor.
  -optical_tensor z1 z2 ..z9
    z1,..9 define the full optical permittivity tensor.
  -masses average
    Specified the mass definition to use, which can be either “program”, “average” or “isotopic”, meaning that the masses used in the calculation of the frequenciesare either taken from the QM program or are the average of the isotope abundances or are the most abundant isotope mass.
  -mass element mass (✔)
    The atomic mass of the element is set to mass.  This can be used to explore the effect of isotope substitution on the calculated frequencies.
  -processors int
    The number of processors to be used  in the calculation can be defined.  By  default all available processors are used.
  -molesof cells \[number_of_atoms_per_molecule\]
    The default for the calculation of molar concentration is “cells”.  Other options are “atoms” or  “molecules”.  If ‘molecules’ is  specified then the number of atoms in a molecules must be provided.
  \-size  0.00001 [sigma]
    Modifies the polarisability (Eq.34) for spherical particles which incorporates the radius of the particle in microns(real).
    It is also used to specify the dimension of the spherical particles for the Mie method.
    It is also possible to specify a size distribution in which case the first number is the mean of the log distribution (in microns) and the second is its width (in log(microns)).

Support matrix
--------------
From an experimental point of view it is often convenient to use a mass fraction rather than a volume fraction to indicate the amount of dielectrically active material present.  PDielec allows mass fractions to be specified instead of a volume fraction, but this requires that the density of the supporting matrix is known. For convenience the package has a small database of the common supporting materials shown in Table below.  These can be specified through the -matrix option. In the case that the properties of the support material are different the properties can be defined instead with the -dielectric and -density options. 



.. table:: Physical properties of matrix materials
   :name: tab-matrix-properties
   :widths: 1 1 1 3
   :column-dividers:   none none  none  none
   :header-alignment: right right right left
   :column-alignment: right right right left

   +--------+---------+--------------+-----------------------------+
   | Name   | Density | Permittivity | Description                 |
   +========+=========+==============+=============================+
   | ptfe   | 2.200   | 2.000        | polytetrafluorethylene      |
   +--------+---------+--------------+-----------------------------+
   | air    | 0.000   | 1.000        | air                         |
   +--------+---------+--------------+-----------------------------+
   | vacuum | 0.000   | 1.000        | vacuum                      |
   +--------+---------+--------------+-----------------------------+
   | kbr    | 2.750   | 2.250        | potassium bromide           |
   +--------+---------+--------------+-----------------------------+
   | nujol  | 0.838   | 2.155        | Nujol                       |
   +--------+---------+--------------+-----------------------------+
   | hdpe   | 0.955   | 2.250        | high density polyethylene   |
   +--------+---------+--------------+-----------------------------+
   | mdpe   | 0.933   | 2.250        | medium density polyethylene |
   +--------+---------+--------------+-----------------------------+
   | ldpe   | 0.925   | 2.250        | low density polyethylene    |
   +--------+---------+--------------+-----------------------------+



Optical permittivity
--------------------

The optical permittivity is normally calculated by the QM or MM program concerned. However, as this property reflects the electronic contribution to the permittivity at zero frequency, unless there is some treatment of electrons by the shell model, then in MM calculations the optical permittivity needs to be defined through the command line options -optical or -optical_tensor.  Unlike the other methods 'mie' method cannot work with particles of zero radius. All methods therefore use a default size of 10\ :superscript:`-12` μm. The Mie approach is only valid for dilute dispersions and for spherical particles. However, if other shapes are specified the Mie method will still be called and the results will be applicable to spheres of the specified size.

Parallelization and threads
---------------------------


To improve the performance of the program python parallelization has been used to parallelize over the frequencies, shapes and methods. By default this parallelization spawns additional Python executables, depending on the number of cores available.

In addition to this form of parallelization the NUMPY library can use multi-threaded BLAS. NUMPY can be compiled with several different BLAS libraries, including; MKL (from Intel), OPENBLAS and ATLAS, or indeed no optimized BLAS library is necessary. To explore the BLAS version compiled with your version of NUMPY please use the test_numpy_1, test_numpy_2 and test_numpy_3 scripts included in the Python/ subdirectory. If you are using MKL or OPENBLAS, the number of threads being used needs to be set with the MKL_NUM_THREADS or OPENBLAS _NUM_THREADS environment variable (sometimes OMP_NUM_THREADS is also used). Use the test routines to determine the optimum for your system. Experience indicates that no performance benefit is obtained with more than two threads. 

In the case of Phonopy the dynamical matrix is read from a yaml file. This has been found to be very slow unless the C parser is used. If the C parser is not available a warning is issued and the program reverts back to the Python parser.

Finally the use of non-standard BLAS libraries seems to cause problems with the affinity settings for the multiprocessing. It has been noticed that the parallel processes can all end up executing on the same processor. In order to prevent this, before executing the pdielec and preader scripts it may be necessary to include;

        export OPENBLAS_MAIN_FREE=1

For some reason this also works if the MKL library is being used.

Performance
------------


PDielec has been written to make use of multiprocessor computers. On a 4 processor machine the speed-up is nearly linear up to 4 processors, as can be seen in the Figure below.

.. _fig-speedup:

.. figure:: ./_static/Figures/SpeedUp.png 
   :scale: 90%

   Speed-up on a four processor workstation

Example command line uses of PDielec
------------------------------------
::

         pdielec -program vasp OUTCAR -method ap -method maxwell \
                 -sphere -plate 0 0 1 -needle 0 0 1 -LO 0 0 1

This performs a calculation using the Averaged-Permittivity and Maxwell-Garnett mixing rules for spherical particles, plate-like particles with a surface (001) and needle-like particles with a unique
direction lying along the \[001\] direction. The supporting matrix is taken to be PTFE and the default volume fraction (10%) is used. The results of a VASP calculation are stored in the current directory.  There is no absorption output from this command as neither the -plot nor the -csv options were specified. The output includes the calculation of the LO modes along the (001) direction. ::

        pdielec -program castep phonon -vmin 300 -vmax 800 \
                -sphere -dielectric 3 -vf 0.1 -vf 0.2 -sigma 10 -csv mgo.csv

This performs a calculation for spherical particles varying the frequency from 300 to 800 |cm-1|, the permittivity of the supporting media is 3, two volume fractions are considered and a damping factor of 10 |cm-1| is used. The results of a CASTEP calculation with the seed-name "phonon" are analysed and the results stored in mgo.csv for further analysis using a spreadsheet. In this example a Maxwell-Garnett mixing rule is used by default. If visual inspection of the results is required then ::

        pdielec -program castep phonon -vmin 300 -vmax 800 \
                -sphere -dielectric 3 -vf 0.1 -vf 0.2\
                -sigma 10 -csv mgo.csv -plot molar_absorption

will perform the same calculation but a graph showing the molar absorption coefficients will be displayed. ::

        pdielec -program gulp calcite.gout -matrix hdpe \
                -method ap -method maxwell -sphere -plate -1 -1 -2 \
                -vmax 2000 -mf 0.1 calcite.gout -csv calcite.csv

This command performs a calculation of the absorption spectrum resulting from a GULP calculation. The supporting matrix density and permittivity are those of high density polyethylene, the frequency range is 0 to 2000 |cm-1|, the mass fraction considered is 10%, the mixing rules used are Averaged-Permittivity and Maxwell-Garnett. Spheres and plates with the (:math:`\bar{1}\bar{1}\bar{2})` surface are considered. ::

        pdielec -program vasp OUTCAR -method mie -sphere \
                -mf 0.1 -size 0.1 -size 1.0 -size 1.0 -csv results.csv

This command performs a calculation of the absorption spectrum resulting from a VASP calculation using the Mie method for sphere with several particles sizes. ::

        pdielec -sphere -mf 0.1 -program experiment experiment.expt\
                -size 0.1 -size 1.0 -size 1.0 -excel results.xlsx

This command performs a calculation of the absorption spectrum resulting from the data stored in the experiment.expt file. Maxwell-Garnett calculations are performed with 3 different sized spheres and the results stored in a Excel file.

Contents of the csv output file
-------------------------------

If a csv output file is requested the file will contain the command used to perform the calculation. A brief summary is given of each active infrared mode; including the mode number, frequency, intensity, integrated molar absorption coefficient, its peak height (calculated from the intensity and damping factor) and the damping parameter used in the calculation. Following this is a table with a column for frequency followed by columns containing the real and imaginary permittivities, the absorption and molar absorption coefficients at each frequency.


PReader
=======

preader is a 'helper' program which uses the underlying modules of PDielec to read output files and summarise the results of various MM/QM packages.  The program can be used to perform some straightforward calculations.  For instance projection of any remaining centre of mass motion of the crystal can be performed to make sure that there are three zero frequencies.  Also the masses used in the calculation of the dynamical matrix can be altered.
Unlike PDielec and PDGui it is not necessary to have performed a full calculation of the dynamical matrix.  In the majority of cases preader will read geometry optimisation runs.

Command options
---------------

Examples of data sets for these packages are included with the distribution and can be found in the Examples/'Package'/preader directory. The program is run from the command line. There are several command options and these are summarized below. Those options which may be repeated are indicated by a ✔ .

  \-program program
    Program can be “abinit”,  “castep”, “crystal”, “gulp”, “qe”, “experiment” or “vasp” and specifies the program which generated the results to be analysed
  \-neutral
    Impose neutrality on the Born charge matrices
  \-nocalculation
    No calculations are performed.  This results in a single line of output with just information from the program.  If -eckart, -mass, -masses, -neutral or -crystal have -hessian crystal have been specified they will be ignored
  \-eckart
    Translational modes are  projected out of the hessian before diagonalisation
  \-hessian symm
    The hessian can be symmetrised in one of two ways; either “crystal” or  “symm”.  In the case of “crystal” the  hessian is symmetrised using the same algorithm as Crystal14. 
  \-masses average
    The element mass_definition can be either “program”, “average” or “isotopic”, meaning that the masses used in the calculation of the frequencies are either taken from the QM program or are the average of the isotope abundances or are the most abundant isotope mass.
  \-mass element mass
    The atomic mass of the element is set to mass.  This can be used to explore the effect of isotope substitution on the calculated frequencies


Examples
--------

::

         preader -program vasp `find . -name OUTCAR` > results.csv

This reads all the VASP OUTCAR files in the current and any of its subdirectories and summarises the results to results.csv. ::

         preader -program castep -eckart `find . -name \*.castep` > results.csv

This reads all the castep output files in the current and any of its subdirectories and summarises the results to results.csv.  For each file the centre of mass motion of the crystal is projected.  The results file contains both the unprojected and the projected results.

There are examples of preader being used in the Examples/'Package'/preader subdirectories of the distribution of PDielec.


MM/QM Interfaces
================

The packages have interfaces to five solid state QM codes, VASP :cite:`Hafner2008c`, CASTEP :cite:`Clark2005d`, CRYSTAL14 :cite:`Dovesi2014`, Abinit :cite:`Gonze2016`, Quantum Espresso :cite:`Giannozzi2009` and Phonopy :cite:`Togo2015`.  In addition an interface is available for GULP :cite:`Gale2003` which is a force field based solid state code. Finally an interface has been written to and 'experiment' file format which allows the preparation of a user defined file specifying the permittivities and absorption frequencies. The origin of the dataset(s) used for processing is determined by a command line switch, -program. An outline of the interfaces to these codes is given here. The package used for the calculation is described by the --program option. In addition a file name is given which contains the output which will be processed by PDielec.

VASP 
-----
The name provided on the command line is an OUTCAR file. The OUTCAR is read by PDielec to determine the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. The VASP run can be a DFPT or numerical calculation of the response.

CASTEP 
-------
The name provided on the command line is the seedname for the calculation. The corresponding seedname.castep file in the current directory is read and processed to determine the unit cell, atomic masses, optical permittivity and born charge tensors. The normal modes and their frequencies are determined from the seedname.phonon file. The CASTEP run needs to be a DFPT (phonon+efield) task.

CRYSTAL 
-------
The name on the command line is a file ending in .out, containing the output of a CRYSTAL14 run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory where PDielec is run from , it uses these files to calculate the Born charge tensors, frequencies and normal modes. The CRYSTAL calculation needs to be a frequency calculation (FREQCALC) with the infrared intensity (INTENS) selected. The default algorithm does not calculate the optical permittivity, so this needs to be provided on the command line. However, if the CPHF or CPKS algorithm is used for the frequency calculation, the optical permittivity is calculated and PDielec will automatically read it from the output file. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the PDielec package. Small differences in the calculated frequencies between the CRYSTAL program and PDielec have been observed. These have been found to be due to a slightly different method for symmetrising the 2^nd^ derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that PDielec should use the same symmetrisation as CRYSTAL14.

ABINIT
------
The output file should come from a run containing three datasets. One to calculate the wavefunction at the optimized geometry, one to calculate the field perturbations and one to calculate the second derivatives. Examples of input files and output files are available with the distribution.

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

To calculate the infrared spectrum PDielec needs the Born charges for the atoms in the unit cell and these can be calculated using VASP and the optimized geometry of the unit cell. The OUTCAR file from this calculation can be copied to the current directory and renamed OUTCAR.born

GULP 
-----

The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit cell. The key words; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it on the command line (see -optical and -optical\_tensor options below).

EXPERIMENT 
----------

This option has been added to allow the exploration of 'toy' problems. The file contains a minimum amount of data to allow a calculation to proceed. It is assumed that the systems will be isotropic as that make the input simpler. Calculations of the LO frequencies will not work with this option. An example input file is given here, which gives results very similar to that found for the Castep, MgO example; ::

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

Examples of data sets for these packages are included with the distribution. The interface to these QM and MM codes reads information about the unit cell, the calculated normal modes and the Born charge matrices; from these the permittivity is calculated over the frequency range requested. The absorption and molar absorption coefficients can be plotted along with the real and imaginary permittivities. Optionally all the information can be written to a comma separated values (csv) file for direct importing into a spreadsheet. The program is run from the command line. There are several command options and these are summarized below. Some of the options may be repeated. The package needs a shape to be specified (sphere, needle, plate or ellipse). If no shape is specified on the command line a sphere is assumed.

