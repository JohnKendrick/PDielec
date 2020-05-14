.. include:: preamble.txt

.. contents::
   :local:


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


===================
Software components
===================

There are two programs within the PDielec package.  *Preader* is a command line program and *PDGui* a graphical interface to the features of the PDielec package.


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

There are several command line options which may be useful in running the package.

.. table:: PDGui command line options
   :name: command_line_options
   :widths: 1 3 
   :column-dividers:   none none  none  none
   :header-alignment: left left
   :column-alignment: left left

   +-----------------------+----------------------------------------+
   | Command line option   | Description                            |
   +=======================+========================================+
   | -debug                | Debugging mode, only for developers    |
   +-----------------------+----------------------------------------+
   | -script filename      | Specify a python script to run         |
   +-----------------------+----------------------------------------+
   | -exit                 | Exit the program                       |
   +-----------------------+----------------------------------------+
   | -nosplash             | Do not show the splash screen          |
   +-----------------------+----------------------------------------+
   | -threads              | Using threading for multiprocessing    |
   +-----------------------+----------------------------------------+
   | -cpus 2               | Using 2 processors for multiprocessing |
   +-----------------------+----------------------------------------+
   | -h                    | Print out help information             |
   +-----------------------+----------------------------------------+



The package requires the PyQt5 library.
After running the program the user sees a notebook interface with seven tabs.

    *Main Tab* allows the user to specify the program and filename which are to be anaylsed.  It is also possible to specify here the spreadsheet name of results are going to be saved.

    *Settings Tab* is used for changing the settings within the package

    *Scenario Tab* specifies the material properties used to define the effective medium from which the  absorption characteristic are calculated.  This tab can also be used to change the algorithms used to calculate the effective medium.

    *Plotting Tab* shows the absorption or permittivity as a function of frequency.

    *Analysis Tab* shows the decomposition of the normal modes into molecular components.

    *Visualisation Tab* displays the normal modes as either arrows showing the atomic displacement or as an animation.

    *Fitter Tab* allows a user to display an experimental spectrum and compare it with the calculated spectrum.


Main Tab
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
The supporting medium may have bubbles of air trapped in the matrix.  For the case that polyethylene spheres are used to make the sample, experimental information indicates a relationship between the size of the spheres and size of the air inclusions.


.. table:: Size of air inclusion (bubbles) in polyethylene
   :name: tab-air_bubbles
   :widths: 1 1 1 3
   :column-dividers:   none none  none  none
   :header-alignment: right right right left
   :column-alignment: right right right left

   +-------------------------------+--------------+----------------------+
   | PE Sphere Diameter (micron)   | Porosity (%) | Void radius (micron) |
   +===============================+==============+======================+
   | 60                            | 2.200        | 2.000                |
   +-------------------------------+--------------+----------------------+
   | 70                            | 0.000        | 1.000                |
   +-------------------------------+--------------+----------------------+
   | 360                           | 0.000        | 1.000                |
   +-------------------------------+--------------+----------------------+

The size of the air inclusions (bubbles) can be specified in the *Scenario Tab* along with the volume fraction.  If a non-zero value for the volume fraction of air voids is given the scattering due to the air inclusions is calculated before the effective medium calculation is performed.  The forward and backward scattering amplitudes for a single scatterer are calculated using the Mie method and then used to calculated the scattering using the Waterman-Truell approximation.


The amount of dielectric material to be considered can be entered either as mass fraction (in percent) or as a volume fraction (in percent).  If the matrix support density is changed the calculated mass fraction will be updated.  It is assumed that the volume fraction has precedence.  If an air void volume fraction is supplied then this has to be taken account of in the calculation of the volume and mass fractions.

The calculation of the effective medium can be performed using a variety of methods which can be chosen from the *Method* drop down menu.  If the *Mie* method is chosen the user can enter the particles radius (microns).  The *Particle sigma* specifies the width of a log-normal distribution.  If the width is 0.0 no sampling of the distribution is performed.

For effective medium theories other than the Mie method the particle shape can be specfied using the *Shape* pull down menu.  Possible shapes are *Sphere*, *Needle*, *Plate* and *Ellipsoid*.  For the cases of *Needle* and *Ellipsoid* the unique direction is specifed by a direction \[abc\] in lattice units.  In the case of *Plate* the unique direction is specifies as the normal to a plane (hkl) in reciprical lattice units.  If an *Ellipsood* shape is used the eccentricy factor can be specified in the *Ellipsoid a/b eccentricity* text box.

The plot in the *Plotting Tab* has a legend and the description used in the legend for each scenario can be set by filling out the text box associated with the *Scenario legend*.

Scenarios can be added or removed using the push buttons at the bottom of each *Scenario Tab*.  When a new scenario is created the settings are copied from the current scenario.  


Plotting Tab
------------

The *Plotting Tab* controls and plots the absorption and permittivity as a function of frequency calculated for each scenario present in the pdgui notebook.


.. _fig-plottingTab:

.. figure:: ./_static/Figures/PlottingTab.png 
   :scale: 80%

   The Plotting Tab

The minimum and maximum frequencies can be specified along with the frequency increment.  An effective medium theory calculation will be performed for each scenario at every frequency between the minimum and maxmimum frequencies at the interval specified.

By default the calculation uses moles of unit cells to calculation the molar absorption coefficient.  This can be altered using the the 'Molar definition' pull down menu.  Options include *Unit cells*, *Atoms* and *Molecules*.  In the case of *Molecules* it is necessary to supply the number of atoms in a formula unit of the compound in the *Number of atoms per molecule* text box.

The title of the plot can be supplied by entering it into the *Plot title* text box and the frequency units used for the plot can also be changed from *wavenumber* to *THz*.

The molar absorption, the absorption and the real or imaginary permittivity can be plotted.  Once a plot has been requested the calculation progress is shown in the progress bar.  Some settings can be changed without the whole plot being recalculated.

Analysis Tab
------------

The *Analysis Tab* shows a breakdown of the phonon modes into molecular components.  The molecular structure of the unit cell is determined by the covalent radii of the atoms.  These can be speficied individually if needed.  The *Analysis Tab* shows the number of molecules that have been found in the analysis.


.. _fig-analysisTab:

.. figure:: ./_static/Figures/AnalysisTab.png 
   :scale: 80%

   The Analysis Tab

The bar graph shows a break down of each normal mode in the chosen frequency range into either internal and external contributions or into molecular components.  The atom sizes and the molecular composition of the unit cell is displayed in the *3D viewer Tab*

3D Viewer Tab
-------------

The *3d viewer Tab* shows the unit cell of the system using the molecular information and atomic sizes from the *Analysis Tab*.

The atomic displacement of each phonon can either be shown as arrows or as an animation.  The views and animations can be recorded in  .png and .mp4 files respectively.  If a .gif file is specified the animation is recorded but with reduced numbers of colours.


.. _fig-viewerTab:

.. figure:: ./_static/Figures/ViewerTab.png 
   :scale: 80%

   The 3D Viewer Tab

As well as being able to change the phonon mode being analysed.  The colours and many settings in the visualiser can be adjusted from the settings tab.


Fitter Tab
-----------

The *Fitter Tab* imports an experimental spectrum.  The spectrum is stored in an Excel spreadsheet.  The spreadsheet should contain a single sheet with two columns.  The first column should be the frequency in cm-1 and the second should be the measured signal.  The signal can be molar absorption, absorption, real permittivity, imaginary permittivity or ATR absorbance.  Once imported the experimental spectrum can be compared with the calculated spectrum.  The tab shows the frequencies contributing to the spectrum and allows the Lorentzian widths of the transitions to be altered.  The frequency range used for the display is the same as that used in the *plotting tab*.
At the top of the Tab are some settings in a tabbed notebook.  The most important is the name of the Excel file containing the experimental spectrum.  In addition there are options to change the 'Plot type', include frequency scaling in any fitting, set the frequency scaling factor, set the number of iterations to be used when fitting, choose whether the plot should use indepent y-axes for the calculated and experimental spectra, set the method used to do the fitting and finally specifiy the spectral difference threshold.
The spectrum is shown at the bottom of the tab and is recalculated when the *Replot* or *Replot with frequency shift* buttons are pressed.  The data type stored in the experimental spreadsheet is defined by the *Plot and data type* setting.  One of the settings options is the ability to remove a baseline from the experimental spectrum.
If baseline removal is selected a Hodrick-Prescott filter is used. The value of the filter parameter is determined by the value in the *HP filter lambda* tab.  The actual value of lambda used in the filter is the value given in the tab raised by the power of 10.


.. _fig-fitter:

.. figure:: ./_static/Figures/FitterTab.png 
   :scale: 80%

   The Fitter Tab


After a replot the cross-correlation coefficient, current frequency scaling factor, the frequency shift needed to maximise that cross-correlation coefficient, and the root mean squared error between the calculated and experimental spectra are shown.  It is possible to apply a frequency scaling to the calculated spectrum. The *Replot with frequency shift* shows the calculated spectrum with the frequencies shifted by the shift calculated to maximise the cross-correlation coefficient.  Care must be taken in comparing the calculated and experimental spectra as different y-axes are used for each.

It is also possible to automatically adjust the Lorentzian width factors with the *Perform fitting* button.  However experience with this option shows that it maybe better to adjust the peak heights manually by altering the sigma values.  At the moment only 20 iterations are performed for each press of the button.  If requested the frequency scaling factor can be adjust too.

Two algorithms for performing the fitting are supported.  Either the cross-correlation coefficient can be maximised or the root mean squared error between the spectra can be minimised.  In the latter case the error is only calculated for signal strengths above the spectral difference threshold (which defaults to 0.05).


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

Scripting
---------
Nearly all aspects of the program can be accessed through scripting.  
The script is provided by the *-script* command line option.
A very useful feature is the ability to save a script which reflects the current state of the program.
By pressing *ctrlS* the user is able to save the settings to a file, which for this example will be called *settings.py*.
The complete set of tabs and their settings are saved in the file and the state of the program can be recreated using;

        pdgui -script settings.py

There are several things to be aware of when using this option.  
The file will contain parameters which have been calculated whilst reading the input file.  
For example the optical permittivities are calculated and stored in this file.  
If the file is then used on another DFT calculation then the optical permittivity for the settings file will be used and not the one which shoiuld have been calculated.  To use the script for other calculations but with the same scenarios the optical permittivities should be removed.

Parallelization and threads
---------------------------

To improve the performance of the program python parallelization has been used to parallelize over the frequencies, shapes and methods. By default this parallelization spawns additional Python executables, depending on the number of cores available.

In addition to this form of parallelization the NUMPY library can use multi-threaded BLAS. NUMPY can be compiled with several different BLAS libraries, including; MKL (from Intel), OPENBLAS and ATLAS, or indeed no optimized BLAS library is necessary. To explore the BLAS version compiled with your version of NUMPY please use the test_numpy_1, test_numpy_2 and test_numpy_3 scripts included in the Python/ subdirectory. If you are using MKL or OPENBLAS, the number of threads being used needs to be set with the MKL_NUM_THREADS or OPENBLAS _NUM_THREADS environment variable (sometimes OMP_NUM_THREADS is also used). Use the test routines to determine the optimum for your system. Experience indicates that no performance benefit is obtained with more than two threads. 

In the case of Phonopy the dynamical matrix is read from a yaml file. This has been found to be very slow unless the C parser is used. If the C parser is not available a warning is issued and the program reverts back to the Python parser.

Finally the use of non-standard BLAS libraries seems to cause problems with the affinity settings for the multiprocessing. It has been noticed that the parallel processes can all end up executing on the same processor. In order to prevent this, before executing the pdielec and preader scripts it may be necessary to include;

        export OPENBLAS_MAIN_FREE=1

For some reason this also works if the MKL library is being used.

There have been issues in running PDielec on Catalina, MacOS.  
This appears to be due to the multiprocessing features of Python not working as expected on this operating system.  
A work around is to use;

        pdgui -threads -cpu 1

This runs the program on a single processor using the threading library for multiprocessing.  
Sometimes it is more convenient to set these using the environment.  Two environment variables are read befor any command line options are processed;

        export PDIELEC_NUM_PROCESSORS=4

will run `pdgui` on 4 processors.  By default the number of processors is determined by interrogating the computer.  This acheived by setting the value of PDIELEC_NUM_PROCESSORS to 0.

        export PDIELEC_THREADING=TRUE

will run `pdgui` using the threading multiprocessing options.

Performance
------------

PDielec has been written to make use of multiprocessor computers. On a 4 processor machine using multiprocessing the speed-up is nearly linear up to 4 processors, as can be seen in the Figure below.  Threading does not appear to be an efficient way of running `PDielec` and it is only recommended if the default multiprocessing options is failing for some reason.

.. _fig-speedup:

.. figure:: ./_static/Figures/SpeedUp.png 
   :scale: 90%

   Speed-up on a four processor workstation


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
The name on the command line is a file ending in .out, containing the output of a CRYSTAL14 run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory where PDielec is run from, it uses these files to calculate the Born charge tensors, frequencies and normal modes. The CRYSTAL calculation needs to be a frequency calculation (FREQCALC) with the infrared intensity (INTENS) selected. The default algorithm does not calculate the optical permittivity, so this needs to be provided on the command line. However, if the CPHF or CPKS algorithm is used for the frequency calculation, the optical permittivity is calculated and PDielec will automatically read it from the output file. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the PDielec package. Small differences in the calculated frequencies between the CRYSTAL program and PDielec have been observed. These have been found to be due to a slightly different method for symmetrising the 2\ :superscript:`nd` derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that PDielec should use the same symmetrisation as CRYSTAL14.

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

The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit cell. The key words; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it.
