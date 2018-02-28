author: John Kendrick
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom.
eMail: j.kendrick@leeds.ac.uk
author: Andrew Burnett
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom
title: PDGUI
email: a.d.burnett@leeds.ac.uk
Bibliography: ./pdielec.bib

[INCLUDE="style"]

<!-- Comment out some of the options -->
<!-- Csl Style: ieee -->
<!-- Math Mode: static -->
<!-- [INCLUDE="style"] -->
<!-- Colorizer: javascript -->
<!-- Doc class: [10pt]article -->

~ MathDefs
\newcommand{\water}{H_{2}O}
\newcommand{\tensor}[1]{\bar{\bar{#1}}}
\newcommand{\tensorbs}[1]{\bar{\bar{\bm{#1}}}}
\newcommand{\tensorbf}[1]{\bar{\bar{\bm{#1}}}}
\newcommand{\fieldbf}[1]{\bar{\bm{#1}}}
~

[TITLE]

[TOC]


# INTRODUCTION

pdgui is provides a graphical front end to the capabilities of PDielec.

# Usage

pdgui can be run from the command line without any parameters or the QM/MM program name followed by the file to be read can be provided.
For example;

    pdgui

or

    pdgui vasp OUTCAR

The package requires the PyQt5 library.
After running the program the user sees a notebook interface with four tabs.

The *Main Tab* allows the user to specify the program and filename which are to be anaylsed

The *Settings Tab* is used for changing the settings within the package

The *Scenario Tab* specifies the material properties used to define the effective medium from which the  absorption characteristic are calculated.  This tab can also be used to change the algorithms used to calculate the effective medium.

The *Plotting Tab* shows the absorption or permittivity as a function of frequency.


## Main tab

The *Main Tab* is used to pick the MM/QM package and the output file which will be analysed.

~ Figure { #fig-mainTab; caption: "The Main Tab"; page-align:inplace }
![img-mainTab]
~
[img-mainTab]: Figures/MainTab.png { width:80%; }

The QM/MM program can be chosen from the dropdown list.  The output file name can be input into the text box below it.  Any entry in this *Output file name* text box will be read when the user types an 'Enter' command.  If the file is not valid a file chooser will pop-up and the user can select a file from that.

Once the file has been specified and read the frequencies found in the calculation file are reported in the output text box.


## Settings Tab

The *Settings Tab* affects the calculation of the frequencies and their intensities.

~ Figure { #fig-settingsTab; caption: "The Settings Tab"; page-align:inplace }
![img-settingsTab]
~
[img-settingsTab]: Figures/SettingsTab.png { width:80%; }

The Eckart conditions control whether 3 modes at zero frequency due to translational invariance are projected out of the dynamical matrix.  By default the Eckart conditions are applied and the 3 modes are projected out.

There is an option to set the sum of the Born charge matrices to zero.  By default this is not applied.  

The atomic masses can be specified in a variety of way.  They can be chosen from;
. The average natural abundance
. The masses used by the QM/MM program
. The isotopic mass of the most abundant isotope

In addition the mass of each element can be edited separately by clicking on the mass concerned and entering a new number.

The width (sigma) of all the absorptions can be set using the *Lorentzian width(sigma)*  spin box. 

Finally the optical permittivity at zero frequency is given in the *Optical permittivity* table.  In some cases it is necessary to enter the optical permittivity by hand.  This can be done by clicking each element in the table which needs changing and typing the new matrix element.

The output table at the bottom of the tab shows the calculated frequencies and their intensities.  Transitions which do not contribute to the Infrared absorption are greyed out.  These transitions will not be used in later calculations.  From this table it is possible to remove or add transitions to the later calculations and it is also possible to change the width of individual transitions.

## Scenario Tabs

There can be more than one *Scenario Tab".  Each one specifies a particular material, method or particle shape which will be used in the *Plotting Tab*.

~ Figure { #fig-scenarioTab; caption: "The Scenario Tab"; page-align:inplace }
![img-scenarioTab]
~
[img-scenarioTab]: Figures/ScenarioTab.png { width:80%; }

The support matrix into which the active dielectric material is dispersed can be selected from the *Support matrix* drop down menu.  The selected supporting material will change the density and permittivity shown in the respective text boxes.  The user can edit these values independently if necessary.  

The amount of dielectric material to be considered can be entered either as mass fraction (in percent) or as a volume fraction (in percent).  If the matrix support density is changed the calculated mass fraction will be updated.  It is assumed that the volume fraction has precedence.

The calculation of the effective medium can be performed using a variety of methods which can be chosen from the *Method* drop down menu.  If the *Mie* method is chosen the user can enter the particles radius (microns).  The *Particle sigma* specifies the width of the LogNormal distribution.  If the width is 0.0 no sampling of the distribution is performed.

For effective medium theories other than the Mie method the particle shape can be specfied using the *Shape* pull down menu.  Possible shapes are *Sphere*, *Needle*, *Plate* and *Ellipsoid*.  For the cases of *Needle* and *Ellipsoid* the unique direction is specifed by a direction [abc] in lattice units.  In the case of *Plate* the unique direction is specifies as the normal to a plane (hkl) in reciprical lattice units.  If an *Ellipsood* shape is used the eccentricy factor can be specified in the *Ellipsoid a/b eccentricity* text box.

The plot in the *Plotting Tab* has a legend and the description used in the legend for each scenario can be set by filling out the text box associated with the *Scenario legend*.

Scenarios can be added or removed using the push buttons at the bottom of each *Scenario Tab*.  When a new scenario is created the settings are copied from the current scenario.  


## Plotting Tab

The *plotting tab* controls and plots the absorption and permittivity as a function of frequency calculated for each scenario present in the pdgui notebook.

~ Figure { #fig-PlottingTab; caption: "The Plotting Tab"; page-align:inplace }
![img-PlottingTab]
~
[img-PlottingTab]: Figures/PlottingTab.png { width:80%; }

The minimum and maximum frequencies can be specified along with the frequency increment.  An effective medium theory calculation will be performed for each scenario at every frequency between the minimum and maxmimum frequencies at the interval specified.

By default the calculation uses moles of unit cells to calculation the molar absorption coefficient.  This can be altered using the the 'Molar definition' pull down menu.  Options include *Unit cells*, *Atoms* and *Molecules*.  In the case of *Molecules* it is necessary to supply the number of atoms in a formula unit of the compound in the *Number of atoms per molecule* text box.

The title of the plot can be supplied by entering it into the *Plot title* text box and the frequency units used for the plot can also be changed from *wavenumber* to *THz*.

The molar absorption, the absorption and the real or imaginary permittivity can be plotted.  Once a plot has been requested the calculation progress is shown in the progress bar.  Some settings can be changed with the whole plot being recalculated.


# INSTALLATION

pdgui requires the PyQT5 libraries in addition to the libraries used by PDielec, these will be installed automatically by matplotlib.

Installation on Microsoft Windows 10 was succesful using Anaconda with Python 3.6

    conda create -n pdielec36 python=3.6
    conda activate pdielec36
    conda install numpy
    conda install scipy
    conda install matplotlib
    pip install termcolor
    pip install xlsxwriter
    pip install pyyaml

Installation on Linux requires pyqt5 modules.  If they have not been installed already with maplotlib.  If you don't have sudo access you can install the module in your user directory;


    sudo pip install pyqt5

or

    pip install --user pyqt5

