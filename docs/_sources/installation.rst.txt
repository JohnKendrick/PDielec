============
Installation
============

.. contents::
   :local:


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


Anaconda
========

Anaconda is a Python package manager and is useful when trying out different versions of Python as there is a very convenient environment manager.  Anaconda can be installed as a user program and so root access to the computer is not needed.
The software can be downloaded from https://conda.io/miniconda.html or https://www.anaconda.com/download/.
Unless you are interested in the full package, miniconda installation should be sufficient for installing the rest of the modules needed for the PDielec installation.
A miniconda is the recommended environment for installing PDielec.  On windows it is necessary to run the commands in an Anaconda command terminal, which will be added to the user's menu when miniconda is installed.

PDielec is available on the conda-forge channel and installation using Anaconda can be performed on Linux, Windows and Mac systems.
However, there are no example files installed, only the executables.  
Example files for each of the DFT packages supported can be downloaded from the Git repository.

The cond-forge package has a Python 3 environment included in the package.  After installation of miniconda or anaconda, a new evironment should be created in which to install the package::

   conda create  --name pdielec
   conda activate pdielec
   conda config --add channels conda-forge
   conda install -c conda-forge pdielec
   pdgui


If at a later stage you want to update PDielec to the latest version from conda-forge, you should be able to update the environment in the following way::

   conda activate pdielec
   conda update pdielec

This only works if the conda-forge channel has been added to the channel list for the environment.


PyPi
====

If you really dont want to use Anaconda, PDielec is available on pypi.org and can be installed using pip.::

   pip install --user pdielec
   pdgui


GitHub - Linux
==============

The package is available on GitHub and can be downloaded from https://github.com/JohnKendrick/PDielec.
cd to a directory where PDielec will be installed and use git to clone a copy of the program.  I use a 'Software/' directory in my home directory to store programs, so for me the commands to obtain PDielec would look like this.::

  cd ~/Software
  git clone https://github.com/JohnKendrick/PDielec.git

This will create a directory \~/Software/PDielec.  To make the program generally available to the user;::

  cd ~/Software/PDielec
  make install

This copies the commands and modules to the \~/bin directory which should be included in your PATH environment variable.::

 export PATH=~/bin;$PATH


As part of the installation you will have to install a number of Python packages into your enviroment.  The full list of packages is as follows;::

   pyqt5
   psutil
   numpy
   scipy
   dill
   multiprocess
   matplotlib
   pyopengl
   openpyxl
   pyyaml
   xlsxwriter
   pymiescatt
   termcolor
   imageio
   imageio-ffmpeg
   sklearn

GitHub - Windows
================

An Anaconda installation is now the recommended way of installing to a Windows machine.  See the section on Anaconda above.  If the user still wishes to proceed with a local installation based purely on the GitHub releases, see below.


Standalone installation
-----------------------

For older releases there is a Windows 10 installation script available on the GitHub at https://github.com/JohnKendrick/PDielec/releases is no longer supported, but an old copy of the package will remain available there.
Releases after version 4.0 have a Windows installable script which installs all the software required to run PDielec, PDGui and PReader.  The installation includes Python 3, and all the modules that are required.  The downloadable executable is very large so be be patient.  The imageio library uses ffmpeg to create videos.  ffmpeg is not distributed with the installation script.  Instead, when creating a video of the phonon motion, PDGui checks for the availability of ffmpeg.  If it is not present then it is downloaded automatically.

Installation from repository
----------------------------

This Windows installation is only needed if the standalone version is found to be too cumbersome.
A Windows 10 installation from the git repository which works for users without administrator rights involves a few steps but gives an installation which can run all the test cases.  In the following instructions replace 'yourusername' with your user name.

Install git
...........

First of all install a Windows version of git from www.git-scm.com. A 64 bit version of Windows 10 will be assumed for the following instructions.

* From the download page download and run the 64 bit Git for Windows setup.
* During the installation install the Quick Launch and Desktop icons as these make using the program easier.
* If you are not familar with the vi or vim editor, it is probably best to use the Nano editor, although if you wish you can install Notepad++ and use that as the default editor.
* In the section concerning the PATH environment I would recommend the last option "Use git and optional Unix tools from the Windows Command prompt".  This option will mean that some Windows commands eg. find and sort will be replaced by the Unix commands, so be careful.
* Leave the https certificate choice as the default, namely the OpenSSL library.
* Line endings are best left to the default setting of Windows-style for checkout and Unix-style for check-in.
* The Console I use is the MinTTY console it has a larger scrolling buffer than the Windows console.
* Under the "Configuring extra options" I leave everything as the default.

I have seen a few hickups in the installation of Git.  Occasionally I have to do the installation twice and occasionally I am left with the Setup Installing window indicating that I should wait, when in fact the installation has completed.  When this happens I kill the setup process with the task manager

Install make
............

Once you have Git installed open the "Git bash" shell and create the following directories in your home directory;::

  mkdir bin
  mkdir Software

From the web url http://www.equation.com/servlet/equation.cmd?fa=make download the 64-bit version of make.exe and copy it into your bin/ directory.  You should find this directory in c:\users\yourusername\bin in the file manager.

Install Python
..............

* From https://www.python.org/downloads/windows/ download and run the Windows x86-64 executable installer for the latest Python 3.x version
* Uncheck install launcher for all users
* Check "Add Python 3.x to PATH"
* Click on "Install now" button
* Check installation ran OK by running the Idle python environment

Open a git bash console and type;::


   pip install pyqt5
   pip install psutil
   pip install numpy
   pip install scipy
   pip install dill
   pip install multiprocess
   pip install matplotlib
   pip install pyopengl
   pip install openpyxl
   pip install pyyaml
   pip install xlsxwriter
   pip install pymiescatt
   pip install termcolor
   pip install imageio
   pip install imageio-ffmpeg
   pip install sklearn

Install PDielec
...............

Open a git bash console and type;::

  cd Software
  git clone -c core.symlinks=true https://github.com/JohnKendrick/PDielec.git

This should create a directory in Software called PDielec.  The "-c core.symlinks=true" means that the commands; pdgui and preader are treated as windows symlinks to their equivalent .py file.

Testing PDielec
...............

Open a git bash console.  If you have installed Python using anaconda then you need to 'source activate' the environment you have established before typing;::

  cd Software/PDielec
  make test_preader
  make test_pdgui

Installing PDielec to run in any git bash console
.................................................

Open a git bash console and type;::

  cd Software/PDielec
  export SCRIPTS=~/bin
  make install

Updating PDielec from the git repository
........................................

Open a git bash console and type;::

  cd Software/PDielec
  git pull

GitHub - PDielec directory structure
====================================

* PDielec/ is the home directory and contains the `pdgui` and `preader` commands
* PDielec/PDielec holds the source for the modules used by the pdielec and preader commands
* PDielec/PDielec/GUI holds the python code for PDGui
* PDielec/Examples a set of examples are available for Abinit, Crystal14, CASTEP, GULP, Phonopy, Mie and VASP.  Each example directory holds the input files to the QM/MM program and the relevant output files which are post processed by PDielec.  For each program there is also a preader directory which holds test output for the preader command.
* PDielec/Sphinx holds the documentation as restructured text documents (.rst).  Sphinx can be used to build the documentation in either html or pdf format.
* PDielec/docs holds the final html documentation.

.. _Examples:

Examples
========

Each example directory has the relevant input data sets use to run the QM/MM program and the output files from that run which are post-processed by PDielec.  There is a file `script.py` which contains an of PDGui usage line and which has been used to create the reference output file `results.ref.xlsx`.  The example can be run interactively::

 make pdgui-view

The output can be compared with the reference data to see if the program is working correctly.  The checkexcell command can be used to do this automatically.

The main Examples/ directory also has a  Makefile file which can be used to verify the correct working of the package.  Simply by typing `make tests` in the Examples directory each example will be run automatically and the output compared with the reference files.  To remove the intermediate files after running the tests automatically, type `make clean`.

A benchmark can be run for comparison of the performance of PDielec on different platforms by typing; ::

    make benchmark

This runs a range of calculations on different systems and provides a real world view of the performance.  An indication of the likely performance of the program is given in the :ref:`performance` section of the documentation.

A summary of the different examples and their purpose is shown below;

.. table:: Summary of the Examples available in the Examples/ directory
   :widths: 2 1 2 8
   :header-alignment: center center center center
   :column-alignment: left center left left
   :column-wrapping: false false false true
   :column-dividers: none single single single none

   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Directory                 | Program          | Molecule      | Description                                                                                                                                     |
   +===========================+==================+===============+=================================================================================================================================================+
   | ATR/AlAs                  | AbInit           | AlAs          | Maxwell Garnett calculation of the ATR spectrum of an ellipsoid along   [001].  The incident angle varies from   0 to 80 degrees.               |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | ATR/Na2SO42               | Vasp             | Na2(SO4)2     | Maxwel-Garnett calculation of the ATR spectrum, change the S polarisation   component from 0 to 100%                                            |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | ATR/Na2SO42_fit           | Vasp             | Na2(SO4)2     | Maxwel-Garnett calculation of the ATR spectrum, example of fitting the spectrum to experiment                                                   |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | AbInit/AlAs               | AbInit           | AlAs          | Average permittivity and Maxwell Garnett calculation of sphere, plate and   ellipsoid                                                           |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | AbInit/BaTiO3             | AbInit           | BaTiO3        | Average permittivity and Maxwell-Garnett calculations of sphere, plate   and ellipsoid, using average isotope masses                            |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | AbInit/BaTiO3-phonana     | AbInit           | BaTiO3        | Average permittivity and Maxwell-Garnett calculations of sphere, plate   and ellipsoid, using program defined masses                            |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | AbInit/Na2SO42            | AbInit           | Na2(SO4)2     | Average permittivity and Maxwell-Garnett calculations of Na2(SO4)2,   sphere, plate and ellipsoid, using program defined masses                 |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Castep/AsparticAcid       | Castep           | Aspartic Acid | Average permittivity and Maxwell-Garnett calculations of sphere, plate   and ellipsoid, using program defined masses                            |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Castep/Bubbles            | Castep           | MgO           | Maxwell-Garnett calculation showing the effect of air bubbles at 24%   volume fraction and 30 micron radius                                     |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Castep/Castep17           | Castep           | beta-Lactose  | Castep 17,Maxwell-Garnett sphere and plates with 3 surfaces                                                                                     |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Castep/Isoleucine         | Castep           | Isoleucine    | Maxwell-Garnett sphere                                                                                                                          |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Castep/MgO                | Castep           | MgO           | Comparison of MG, Bruggeman and AP methods changing shapes and volume   fractions                                                               |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Castep/Na2SO42            | Castep           | Na2(SO4)2     | Comparison of MG and Bruggeman, for needle, ellipsoid and plate shapes                                                                          |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Crystal/Leucine           | Crystal          | Leuscine      | Comparison of MG, plates and ellipsoids                                                                                                         |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Crystal/Na2SO42           | Crystal          | Na2(SO4)2     | Comparison of MG for needle, ellipsoid and plate shapes                                                                                         |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Crystal/Na2SO42_C17       | Crystal          | Na2(SO4)2     | Comparison of MG for needle, ellipsoid and plate shapes, reading output   from Crystal 17                                                       |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Crystal/Quartz            | Crystal          | Quartz        | Comparison of MG for needle, ellipsoid and plate shapes                                                                                         |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Crystal/ZnO/CPHF          | Crystal          | ZnO           | Coupled Hartree-Fock, Maxwell-Garnett Sphere, Needle and Plate                                                                                  |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Crystal/ZnO/Default       | Crystal          | ZnO           | Default Crystal calculation of IR spectrum, Maxwell-Garnett Sphere,   Needle and Plate                                                          |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Crystal/ZnO/NoEckart      | Crystal          | ZnO           | As above, but no Eckart projection in Crystal,  Maxwell-Garnett Sphere, Needle and Plate                                                        |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Experiment/Forsterite     | Experiment       | Forsterite    | Single crystal calculations of a thick slab, for a, b and c axis   alignments with polarisation direction.    Uses FPSQ model for permittivity. |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Experiment/Mayerhofer     | Experiment       | Toy model     | Example of a Drude Lorentz model permittivity                                                                                                   |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Experiment/constant       | Experiment       | Constant      | Example of a constant permittivity with loss                                                                                                    |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Experiment/drude-lorentz  | Experiment       | MgO           | A Drude-Lorentz model for MgO, varying the angle of incidence                                                                                   |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Experiment/fpsq           | Experiment       | Quartz        | An FPSQ model for Quartz, showing polarisation on along different axes   and different incident angles.                                         |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Experiment/interpolation  | Experiment       | Quartz        | Under development                                                                                                                               |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Gulp/Na2SO42              | Gulp             | Na2(SO4)2     | Maxwell-Garnett and Bruggeman on needle, ellipsoid and plate                                                                                    |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Gulp/calcite              | Gulp             | Calcite       | Maxwell-Garnett method on Sphere and Plate                                                                                                      |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Mie/MgO                   | Castep           | MgO           | Mie method with varying volume fractions and sphere sizes                                                                                       |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Mie/MgO_lognormal         | Castep           | MgO           | Mie method with varying volume fractions and sphere size distributions                                                                          |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Phonopy/Na2SO42           | Phonopy          | Na2(SO4)2     | Maxwell-Garnett and Bruggeman method for needle, ellipsoid and plate   shapes, with varying volume fractions                                    |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Phonopy/ZnO               | Phonopy          | ZnO           | Maxwell-Garnett and Bruggeman method for needle, ellipsoid and plate   shapes                                                                   |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | QE/Cocaine                | Quantum Espresso | Cocaine       | Maxwel-Garnett sphere                                                                                                                           |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | QE/Na2SO42                | Quantum Espresso | Na2(SO4)2     | Maxwell-Garnett and Bruggeman on needle, ellipsoid and plate                                                                                    |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | QE/ZnO                    | Quantum Espresso | ZnO           | Maxwell-Garnett and Bruggeman on needle, ellipsoid and plate                                                                                    |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | SingleCrystal/Bi2Se3      | Vasp             | Bi2Se3        | Single crystal example of thick slab, angle of incedence varies from 0 to   90                                                                  |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | SingleCrystal/Bi2Se3_film | Vasp             | Bi2Se3        | Single crystal example of thin film, angle of incedence varies from 0 to   90                                                                   |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | SizeEffects/BaTiO3        | Abinit           | BaTiO3        | Exploration of size effects in Bruggeman effective medium theory                                                                                |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | SizeEffects/MgO           | Castep           | MgO           | Exploration of size effects in Bruggeman and Maxwell-Garnett effective   medium theories                                                        |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | SizeEffects/ZnO           | Vasp             | ZnO           | Exploration of size effects in Maxwell-Garnett effective medium theory                                                                          |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Vasp/F-Apatite            | Vasp             | F-Apatite     | Maxwell-Garnett, sphere plates and needles                                                                                                      |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Vasp/Na2SO42              | Vasp             | Na2(SO4)2     | Maxwell-Garnett and Bruggeman, needle, plate and needle                                                                                         |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
   | Vasp/ZnO                  | Vasp             | ZnO           | Maxwell-Garnett and Bruggeman, needle, plate and needle, mass fraction                                                                          |
   +---------------------------+------------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------------------+
