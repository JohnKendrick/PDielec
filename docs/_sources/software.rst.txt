.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


===================
Software components
===================

There are two programs within the PDielec package.  *preader* is a command line program, which allows the user to process and ectract information from QM/MM package output, and *pdgui* a graphical interface to the features of the PDielec library.


PDGui
=====


PDgui provides a graphical front end to the capabilities of PDielec.
PDgui can be run from the command line without any parameters or followed by the file to be read in.  PDGui makes a best guess at the DFT/MM program used to create the output.
It can also be used a command line tool without having to use the GUI, this needs a script file to be available.

For example; ::

    pdgui

or ::

    pdgui OUTCAR                                                # A Vasp/Phonpy calculation
    pdgui LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out    # A Crystal calculation
    pdgui run1.dynG                                             # A QE calculation
    pdgui run2.out                                              # An Abinit calculation
    pdgui run3.gout                                             # An Gulp calculation

If running a script which contains program and output file information, then::

    pdgui -script script.py

If a spreadsheet of results is required the name of the spreadsheet can be provided on the command line. ::

    pdgui OUTCAR -spreadsheet results.xlsx
    
The spreadsheet name must have an .xlsx extension.
To use as a command line tool an example would be; ::

    pdgui OUTCAR -spreadsheet results.xlsx -script script.py -nosplash -exit

There are several command line options which may be useful in running the package.

.. table:: PDGui command line options
   :name: command_line_options
   :widths: 1 3 
   :column-dividers:   single single  single 
   :header-alignment: left left
   :column-alignment: left left

   +-----------------------+----------------------------------------+
   | Command line option   | Description                            |
   +=======================+========================================+
   | -script filename      | Specify a python script to run         |
   +-----------------------+----------------------------------------+
   | -scenario type        | Changes the initialial scenario type,  |
   |                       | type can be 'powder' or 'crystal',     |
   |                       | the default is 'powder'                |
   +-----------------------+----------------------------------------+
   | -spreadsheet file     | An alternative way of defining the     |
   |                       | spreadsheet                            |
   +-----------------------+----------------------------------------+
   | -program name         | An alternative way of defining the     |
   |                       | program used to create the output      |
   +-----------------------+----------------------------------------+
   | -version              | Give information of the version no.    |
   +-----------------------+----------------------------------------+
   | -debug                | Debugging mode, only for developers    |
   +-----------------------+----------------------------------------+
   | -exit                 | Exit the program after running the     |
   |                       | script                                 |
   +-----------------------+----------------------------------------+
   | -nosplash             | Do not show the splash screen          |
   +-----------------------+----------------------------------------+
   | -cpus 2               | Using 2 processors for multiprocessing |
   |                       | or 2 threads if threading              |
   +-----------------------+----------------------------------------+
   | -threads              | Using threading for multiprocessing    |
   +-----------------------+----------------------------------------+
   | -h                    | Print out help information             |
   +-----------------------+----------------------------------------+



After running the program the user sees a notebook interface with seven tabs.

-   *Main Tab* allows the user to specify the program and filename which are to be anaylsed.  It is also possible to specify here the spreadsheet name to which results are saved to.  

-   *Settings Tab* is used for changing the settings within the package

-   *Scenario Tab* can be a powder scenario or a single crystal scenario.  In the case of a powder scenario the tab specifies the material properties used to define the effective medium from which the  absorption characteristic are calculated.  In the case of a single crystal the orientation of the the crystal can be defined and the angle of incidence of the radiation.

-   *Plotting Tab* shows the results of all the scenario tabs in one place.

-   *Analysis Tab* shows the decomposition of the normal modes into molecular components.

-   *Visualisation Tab* displays the normal modes as either arrows showing the atomic displacement or as an animation.

-   *Fitter Tab* allows a user to display an experimental spectrum and compare it with the calculated spectrum.


Main Tab
--------

The *Main Tab* is used to pick the MM/QM package and the output file which will be analysed.

.. _fig-mainTab:

.. figure:: ./_static/Figures/MainTab.png 
   :scale: 80%

   The Main Tab

The QM/MM program can be chosen from the dropdown list.  The output file name can be input into the text box below it, or it can be specified using the `File manager` button.  In addition to input files from QM/MM programs the program can read Python script files written by PDGui.  These should have a '.py' extension.

If an Excel spreadsheet file name given, it must have the extensions .xlsx.  The spreadsheet is written when the `Save results` button is pressed.

If the current program settings need to be stored they can be stored as a Python script file by pressing the `Save settings` button.  Any settings saved in the way can be re-instated by reading in the script using the `File manager` and selecting the `PDGui (*.py)` option in the file types selection.

Once the file has been specified and read, the unit cell is shown and the frequencies found in the calculation file are reported in the output text boxes.


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

When a new output file is read the atomic masses and the optical permittivity matrix are reset.

The output table at the bottom of the tab shows the calculated frequencies and their intensities.  Transitions which do not contribute to the Infrared absorption are greyed out.  These transitions will not be used in later calculations.  From this table it is possible to remove or add transitions to the later calculations and it is also possible to change the width of individual transitions.

Scenario Tabs
-------------
There are two sorts of scenario tab, one for powder calculations and one for single crystal calculations.  By default it is assumed that powder calculations are being performed.  There is a button ( `Switch to crystal scenario` or `Switch to powder scenario`) at the bottom of the tabs which allows the user to switch between the two.

Scenarios can be added or removed using the push buttons at the bottom of each *Scenario Tab*.  When a new scenario is created the settings are copied from the current scenario.  

There can be more than one *Scenario Tab*.  Each one specifies a particular material, method or particle shape which will be used for the calculation of the effective medium or single crystal.  The results of the calculations will be shown in the *Plotting Tab*.

Powder Scenarios
................

For each powder scenario it is assumed that the size of the particle embedded in the matrix is small compared with the wavelength of light.


.. _fig-scenarioTab-Powder:

.. figure:: ./_static/Figures/ScenarioTab-Powder.png 
   :scale: 80%

   The Scenario Tab for a Powder

The support matrix into which the active dielectric material is dispersed can be selected from the *Support matrix* drop down menu.  The selected supporting material will change the density and permittivity shown in the respective text boxes.  The user can edit these values independently if necessary.  
The supporting medium may have bubbles of air trapped in the matrix.  For the case that polyethylene spheres are used to make the sample, experimental information indicates a relationship between the size of the spheres and size of the air inclusions :cite:`Kaushik2014`.


.. table:: Size of air inclusion (bubbles) in polyethylene
   :name: tab-air_bubbles
   :widths: 1 1 1 3
   :column-dividers:   none none  none  none
   :header-alignment: right right right left
   :column-alignment: right right right left

   +-------------------------------+--------------+----------------------+
   | PE Sphere Diameter (micron)   | Porosity (%) | Void radius (micron) |
   +===============================+==============+======================+
   | 60                            | 24           | 24.0                 |
   +-------------------------------+--------------+----------------------+
   | 72                            | 25           | 28.0                 |
   +-------------------------------+--------------+----------------------+
   | 360                           | 44.5         | 90.0                 |
   +-------------------------------+--------------+----------------------+

The size of the air inclusions (bubbles) can be specified in the *Scenario Tab* along with the volume fraction.  If a non-zero value for the volume fraction of air voids is given, the scattering due to the air inclusions is calculated before the effective medium calculation is performed.
The forward and backward scattering amplitudes for a single scatterer are calculated using the Mie method and then used to calculated the scattering using the Waterman-Truell approximation.

The amount of dielectric material can be entered either as a mass fraction (in percent) or as a volume fraction (in percent).
If the matrix support density is changed the calculated mass fraction will be updated.
It is assumed that the volume fraction has precedence.
If an air void volume fraction is supplied then this is taken account of in the calculation of the volume and mass fractions.

The calculation of the effective medium can be performed using a variety of methods which can be chosen from the *Method* drop down menu.
If the *Mie* method is chosen the user can enter the particles radius (microns).
*Particle sigma* specifies the width of a log-normal distribution.
If the width is 0.0 no sampling of the distribution is performed.

For effective medium theories, the particle shape can be specfied using the *Shape* pull down menu.
Possible shapes are *Sphere*, *Needle*, *Plate* and *Ellipsoid*.
For the cases of *Needle* and *Ellipsoid* the unique direction is specifed by a direction \[abc\] in lattice units.
In the case of *Plate* the unique direction is specifies as the normal to a plane (hkl) in reciprical lattice units.
If an *Ellipsoid* shape is used the eccentricity factor can be specified in the *Ellipsoid a/b eccentricity* text box.

The plot in the *Plotting Tab* has a legend and the description used in the legend for each scenario can be set by filling out the text box associated with the *Scenario legend*.


Single Crystal Scenarios
........................

Single crystals are defined by a surface and can be either thin coherent films or thick slabs.  A thick slab assumes that the thickness is such that no radiation can be transmitted and all radiation is therefore reflected or absorbed.  As such only the reflectance has any meaning for thick slabs.



.. _fig-scenarioTab-SingleCrystal:

.. figure:: ./_static/Figures/ScenarioTab-SingleCrystal.png 
   :scale: 80%

   The Scenario Tab for a Single Crystal

The mode of calculation is selected at the top of the tab, followed by the definition of the crystal surface, the azimuthal angle and the angle incidence.


This information is sufficient to link the crystal axis system with the laboratory frame.  
The laboratory axes are defined with Z- being the normal to the surface (hkl), and the incident and reflected beams lie in the XZ- plane, which means the laboratory Y- axis is perpendicular to the XZ- plane. 
The incident radiation is assumed to be at an angle :math:`\theta` to the normal. 
p-polarised light lies in the XZ- plane and s-polarised light is parallel with the Y-axis.  The azimuthal angle :math:`\phi` defines the rotation of the crystal around the Z-axis.
To help understand the relationship between the crystal and laboratory co-ordinates the program outputs the direction of the lattice vectors in terms of the laboratory X-, Y- and Z- coordinates.

The single crystal scenario also allows the specification of permittivities for the incident and substrate (in the case of thick slab this is not needed).  Finally from thin films the thickness of the film in nanometres is specified.


Plotting Tab
------------

The *Plotting Tab* controls and plots the absorption and permittivity as a function of frequency calculated for each scenario present in the pdgui notebook.
The plots available depends on the availability of powder and single crystal scenarios.


.. _fig-plottingTab-Powder:

.. figure:: ./_static/Figures/PlottingTab-Powder.png 
   :scale: 80%

   The Plotting Tab Showing Powder Absorption

The minimum and maximum frequencies can be specified along with the frequency increment.  An effective medium theory calculation will be performed for each scenario at every frequency between the minimum and maxmimum frequencies at the interval specified.

By default the calculation uses moles of unit cells to calculation the molar absorption coefficient.  This can be altered using the the 'Molar definition' pull down menu.  Options include *Unit cells*, *Atoms* and *Molecules*.  In the case of *Molecules* it is necessary to supply the number of atoms in a formula unit of the compound in the *Number of atoms per molecule* text box.

The title of the plot can be supplied by entering it into the *Plot title* text box and the frequency units used for the plot can also be changed from *wavenumber* to *THz*.

For powder scenarios the molar absorption, the absorption and the real or imaginary permittivity can be plotted.  Once a plot has been requested the calculation progress is shown in the progress bar.  Some settings can be changed without the whole plot being recalculated.

For single crystal plots the crystal reflectance and transmittance for s and p polarised radiation can be plotted along with the absorptance which is defined in terms of the reflectance (R) and the transmittance (T) as :math:`1-R-T`.
For thick slabs only the transmittance is of any revelance.  For thin films the absorptance is useful as it removes some of the oscilations that occur in the transmittance or reflectance due to the film thickness.

.. _fig-plottingTab-SingleCrystal:

.. figure:: ./_static/Figures/PlottingTab-SingleCrystal.png 
   :scale: 80%

   The Plotting Tab showing Single Crystal results

For the case of a thick slab, the permittivity of the incident can be specified, but the substrate permittivity is not used.  For the case of a thin film both permittivities are used.

The title of the plot can be supplied by entering it into the *Plot title* text box and the frequency units used for the plot can also be changed from *wavenumber* to *THz*.  A default plot title is generated giving information about the type of plot being viewed.


Analysis Tab
------------

The *Analysis Tab* shows a breakdown of the phonon modes into molecular components or into internal and external modes.  The molecular structure of the unit cell is determined by the covalent radii of the atoms.  These can be speficied individually if needed.  The *Analysis Tab* shows the number of molecules that have been found in the analysis.


.. _fig-analysisTab:

.. figure:: ./_static/Figures/AnalysisTab.png 
   :scale: 80%

   The Analysis Tab

The bar graph shows a break down of each normal mode in the chosen frequency range into either internal and external contributions or into molecular components.  The atom sizes and the molecular composition of the unit cell is displayed in the *3D Viewer Tab*

3D Viewer Tab
-------------

The *3D Viewer Tab* shows the unit cell of the system using the molecular information and atomic sizes from the *Analysis Tab*.

The atomic displacement of each phonon can either be shown as arrows or as an animation.  The views and animations can be recorded in  .png and .mp4 files respectively.  If a .gif file is specified the animation is recorded but with reduced numbers of colours.  A series of cif files can also be generated which follows the motion of the mode being studied.


.. _fig-viewerTab:

.. figure:: ./_static/Figures/ViewerTab.png 
   :scale: 80%

   The 3D Viewer Tab

As well as being able to change the phonon mode being analysed.  The colours and many settings in the visualiser can be adjusted from the settings tab.  The unit cell can be shown as a supercell by altering the a, b, c parameters in the settings tab.  Other parameters affecting the display can also be modified.


Fitter Tab
-----------

The *Fitter Tab* imports an experimental spectrum.  The spectrum is stored in an Excel spreadsheet.  The spreadsheet should contain a single sheet with two columns.  The first column should be the frequency in |cm-1| and the second should be the measured signal.  The signal can be molar absorption, absorption, real permittivity, imaginary permittivity, ATR absorbance or any of the single crystal propertie.  Once imported the experimental spectrum can be compared with the calculated spectrum.  The tab shows the frequencies contributing to the spectrum and allows the Lorentzian widths of the transitions to be altered.  The frequency range used for the display is the same as that used in the *plotting tab*.
Indeed the calculated signal used for fitting is the same as that being shown in the *plotting tab*.
At the top of the tab are some settings in a tabbed notebook.  The most important is the name of the Excel file containing the experimental spectrum.  In addition there are options to change the 'Plot type', include frequency scaling in any fitting, set the frequency scaling factor, set the number of iterations to be used when fitting, choose whether the plot should use indepent y-axes for the calculated and experimental spectra, set the method used to do the fitting and finally specifiy the spectral difference threshold.
The spectrum is shown at the bottom of the tab and is recalculated when the *Replot* or *Replot with frequency shift* buttons are pressed.  The data type stored in the experimental spreadsheet is defined by the *Plot and data type* setting.  One of the settings options is the ability to remove a baseline from the experimental spectrum.
If baseline removal is selected a Hodrick-Prescott filter is used. The value of the filter parameter is determined by the value in the *HP filter lambda* tab.  The actual value of lambda used in the filter is the value given in the tab raised by the power of 10.


.. _fig-fitter:

.. figure:: ./_static/Figures/FitterTab.png 
   :scale: 80%

   The Fitter Tab


After a replot the cross-correlation coefficient, current frequency scaling factor, the frequency shift needed to maximise that cross-correlation coefficient, and the root mean squared error between the calculated and experimental spectra are shown.  It is possible to apply a frequency scaling to the calculated spectrum. The *Replot with frequency shift* shows the calculated spectrum with the frequencies shifted by the shift calculated to maximise the cross-correlation coefficient.  Care must be taken in comparing the calculated and experimental spectra as different y-axes are used for each.

It is also possible to automatically adjust the Lorentzian width factors with the *Perform fitting* button.  However experience with this option shows that it maybe better to adjust the peak heights manually by altering the sigma values.  At the moment only 20 iterations are performed for each press of the button.  If requested the frequency scaling factor can be adjusted too.

Two algorithms for performing the fitting are supported.  Either the cross-correlation coefficient can be maximised or the root mean squared error between the spectra can be minimised.  In the latter case the error is only calculated for signal strengths above the spectral difference threshold (which defaults to 0.05).

A summary of the optional settings available in the *Fitter Tab* is shown below;

.. table:: Fitter Tab options
   :name: fitter_tab_options
   :widths: 1 3 
   :column-dividers:   none none  none 
   :header-alignment: left left
   :column-alignment: left left

   +--------------------------+-----------------------------------------------------+
   | Fitter options           | Description                                         |
   +==========================+=====================================================+
   | Experimental spectrum    | Provide a file with the experimental                |
   |                          | data. A return brings up a file manager             |
   +--------------------------+-----------------------------------------------------+
   | Scenario                 | Allows the user to choose which                     |
   |                          | will be used in the fitting                         |
   +--------------------------+-----------------------------------------------------+
   | Frequency scaling factor | The calculated frequencies are scaled               |
   +--------------------------+-----------------------------------------------------+
   | No. of iterations        | Number of iterations for automatic                  |
   |                          | fitting                                             |
   +--------------------------+-----------------------------------------------------+
   | Independent y-axes       | The calculated spectrum uses the left               |
   |                          | axis and the experimental the right                 |
   +--------------------------+-----------------------------------------------------+
   | Fitting type             | Refers to the method used for automatically fitting |
   |                          |     - minimise x-correlation                        |
   |                          |     - minimise spectral difference                  |
   +--------------------------+-----------------------------------------------------+
   | Optimise scaling         | Tick if frequency scaling is to be included in the  |
   |                          | automatic optimisation                              |
   +--------------------------+-----------------------------------------------------+
   | Spectral difference      | Used if fitting type is minimise spectral difference|
   | threshold                |                                                     |
   +--------------------------+-----------------------------------------------------+
   | Baseline removal?        | Tick to request that the Hodrick-Prescott (HP)      |
   |                          | method is used to remove a sloping baseline         |
   |                          | in the experimental data                            |
   +--------------------------+-----------------------------------------------------+
   | HP filter Lambda         | Provide the lambda value in the HP method           |
   +--------------------------+-----------------------------------------------------+

Scripting
---------
Nearly all aspects of the program can be accessed through scripting.  
The script is provided by the *-script* command line option or by choosing the PDGui (\*.py) file type from the file manager on the *Main tab*.
A very useful feature is the ability to save a script which reflects the current state of the program.  This is available from the *Main tab*.

There are several things to be aware of when using this option to read in settings from one calculation into another.  
The file will contain parameters which have been calculated whilst reading the input file.  
For example the optical permittivities are calculated and stored in this file.  
If the file is then used on another DFT calculation then the optical permittivity for the settings file will be used and not the one which should have been calculated.  This is also true if masses have been defined in once calculation, they may not be appropriate for a different molecule. To use the script for other calculations but with the same scenarios the optical permittivities should be removed and any mass specifications removed. 

Parallelisation and threads
---------------------------

To improve the performance of the program python parallelisation has been used to parallelise over the frequencies, shapes and methods. By default this parallelisation spawns additional Python executables, depending on the number of cores available.

In addition to this form of parallelisation the NUMPY library can use multi-threaded BLAS. NUMPY can be compiled with several different BLAS libraries, including; MKL (from Intel), OPENBLAS and ATLAS, or indeed no optimised BLAS library is necessary. To explore the BLAS version compiled with your version of NUMPY please use the test_numpy_1, test_numpy_2 and test_numpy_3 scripts included in the Python/ subdirectory. If you are using MKL or OPENBLAS, the number of threads being used needs to be set with the MKL_NUM_THREADS or OPENBLAS _NUM_THREADS environment variable (sometimes OMP_NUM_THREADS is also used). Use the test routines to determine the optimum for your system. Experience indicates that no performance benefit is obtained with more than two threads. 

In the case of Phonopy the dynamical matrix is read from a yaml file. This has been found to be very slow unless the C parser is used. If the C parser is not available a warning is issued and the program reverts back to the Python parser.

Finally the use of non-standard BLAS libraries seems to cause problems with the affinity settings for the multiprocessing. It has been noticed that the parallel processes can all end up executing on the same processor. In order to prevent this, before executing the pdielec and preader scripts it may be necessary to include;

        export OPENBLAS_MAIN_FREE=1

For some reason this also works if the MKL library is being used.

There have been issues in running PDielec on Catalina, MacOS.  
This appears to be due to the multiprocessing features of Python not working as expected on this operating system.  
A work around is to use;

        pdgui -threads -cpu 1

This runs the program on a single processor using the threading library for multiprocessing.  
Sometimes it is more convenient to set these using the environment.  Two environment variables are read before any command line options are processed;

        export PDIELEC_NUM_PROCESSORS=4

will run `pdgui` on 4 processors.  By default the number of processors is determined by interrogating the computer.  This is acheived by setting the value of PDIELEC_NUM_PROCESSORS to 0.

        export PDIELEC_THREADING=TRUE

will run `pdgui` using the threading multiprocessing options.

.. _performance:

Performance
------------

PDielec has been written to make use of multiprocessor computers. On an 8 processor (16 thread) machine using multiprocessing the speed-up soon slows down after three processors, as there is a lot of the code which is not parallelised .  Threading does not appear to be an efficient way of running `PDielec` and it is only recommended if the default multiprocessing options is failing for some reason.
The speedup is calculated by running the benchmark suite (pdmake benchmarks), which consists of 30 different examples of using the code and is therefore representative of the type of calculations performed.

.. _fig-speedup:

.. figure:: ./_static/Figures/SpeedUp.png 
   :scale: 90%

   Speed-up on a eight processor workstation

Using the same suite of benchmarks on an 8 processor Intel core i9-9900K running at 3.6GHz took 92.6s.

On a Windows laptop with an Intel Skylake i7-650, which has 2 processors running at 2.5GHz, the benchmark ran in 432s seconds.  On the same machine but running the Linux operating system the benchmark suite took 358s, making the impact of Windows on performance a factor of 1.2 times slower.


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


The Format of the Excel Spreadsheet
===================================

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



