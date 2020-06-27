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

PDielec is available on the conda-forge channel and can be installed as show in the following couple of examples.

Installation using Anaconda can be performed on Linux, Windows and Mac systems.
However, there are no example files installed, only the executables.  
Example files for each of the DFT packages supported can be downloaded from the Git repository.

The cond-forge package has a python 3 environment included in the package::

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

   numpy
   scipy
   pyopengl
   matplotlib
   pyyaml
   termcolor
   xlsxwriter
   openpyxl
   pyyaml
   imageio
   imageio-ffmpeg

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
* In the section concerning the PATH environment I would recommend the last option "Use git and optional Unix tools from the Windows Command prompt".  This option will mean that some Windows commands eg. find and sort will be replace by the Unix commands.
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

  pip install numpy
  pip install scipy
  pip install pyopengl
  pip install matplotlib
  pip install PyYaml
  pip install XlsxWriter
  pip install openpyxl
  pip install termcolor
  pip install imageio
  pip install imageio-ffmpeg
  pdgui


Install PDielec
...............

Open a git bash console and type;::

  cd Software
  git clone -c core.symlinks=true https://github.com/JohnKendrick/PDielec.git

This should create a directory in Software called PDielec.  The "-c core.symlinks=true" means that the commands; pdielec, pdgui and preader are treated as windows symlinks to their equivalent .py file.

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

* PDielec/ is the home directory and contains the `pdielec` and `preader` commands
* PDielec/Python holds the source for the modules used by the pdielec and preader commands
* PDielec/Python/GUI holds the python code for PDGui
* PDielec/Examples a set of examples are available for Abinit, Crystal14, CASTEP, GULP, Phonopy, Mie and VASP.  Each example directory holds the input files to the QM/MM program and the relevant output files which are post processed by PDielec.  For each program there is also a preader directory which holds test output for the preader command.
* PDielec/Sphinx holds the documentation as restructured text documents (.rst).  Sphinx can be used to build the documentation in either html or pdf format.
* PDielec/docs holds the final html documentation.

Examples
========

Each example directory has the relevant input data sets use to run the QM/MM program and the output files from that run which are post-processed by PDielec.  There is a file `script.py` which contains an of PDGui usage line and which has been used to create the reference output file `results.ref.xlsx`.  The example can be run interactively::

 make pdgui-view

The output can be compared with the reference data to see if the program is working correctly.

The main Examples/ directory also has a  Makefile file which can be used to verify the correct working of the package.  Simply by typing `make` in the Examples directory each example will be run automatically and the output compared with the reference files.  To remove the intermediate files after running the tests automatically, type `make clean`.

