============
Installation
============

.. contents::
   :local:


Linux
=====

The package is available on GitHub and can be downloaded from https://github.com/JohnKendrick/PDielec.
cd to a directory where PDielec will be installed and use git to clone a copy of the program.  I use a 'Software/' directory in my home directory to store programs, so for me the commands to obtain PDielec would look like this.::

  cd ~/Software
  git clone https://github.com/JohnKendrick/PDielec.git

This will create a directory \~/Software/PDielec.  To make the program generally available to the user;::

  cd ~/Software/PDielec
  make install

This copies the commands and modules to the \~/bin directory which should be included in your PATH environment variable.::

 export PATH=~/bin;$PATH


Anaconda
--------

Anaconda is a Python package manager and is useful when trying out different versions of Python as there is a very convenient environment manager.  Anaconda can be installed as a user program and so root access to the computer is not needed.  The software can be downloaded from https://conda.io/miniconda.html or https://www.anaconda.com/download/.  Unless you are interested in the full package the miniconda installation should be sufficient for installing the rest of the modules needed for the PDielec installation.  

As an example here is an installation using a Python 2.7 environment.  Using 'pip install' means that the anaconda build system is bypassed.ld system in Anaconda normally makes sure only compatible software versions are installed, but this seems to cause problems in this case. ::


   conda create  --name pdielec2 python=2
   conda install --name pdielec2 scipy
   source activate pdielec2
   pip install matplotlib
   pip install termcolor
   pip install xlsxwriter
   pip install openpyxl
   pip install pyyaml
   pip install imageio

Here is an example of creating a Python 3 environment in Anaconda::

   conda create  --name pdielec3 python=3
   conda install --name pdielec3 -c conda-forge scipy
   source activate pdielec3
   pip install matplotlib
   pip install termcolor
   pip install xlsxwriter
   pip install openpyxl
   pip install pyyaml
   pip install imageio


Pip (only needed if not using Anaconda)
---------------------------------------

If you are working in an environment with Python already installed then it is possible to use pip to add the additional modules that are needed.  If the user does not have administrator rights on the machine then the use of pip install --user will install the python modules in a user's directory.  If the user does have administrator rights and wants to install the packages for all users then the '--user' option can be omitted.::


   pip install --user numpy
   pip install --user scipy
   pip install --user matplotlib
   pip install --user PyYaml
   pip install --user XlsxWriter
   pip install --user openpyxl
   pip install --user termcolor
   pip install imageio


Windows
=======

Standalone installation
-----------------------

A Windows 10 installation script is available on the GitHub at
https://github.com/JohnKendrick/PDielec/releases .
Releases after version 4.0 have a Windows installable script which installs all the software required to run PDielec, PDGui and PReader.  The installation includes Python 3, and all the modules that are required.  The downloadable executable is very large so be be patient.  The imagio library uses ffmpeg to create videos.  ffmpeg is not distributed with the installation script.  Instead, when creating a video of the phonon motion, PDGui checks for the availability of ffmpeg.  If it is not present then it is downloaded automatically.

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

Anaconda
........

Anaconda is a Python package manager and is useful when trying out different versions of Python as there is a very convenient environment manager.  Anaconda can be installed as a user program and so root access to the computer is not needed.  The software can be downloaded from https://conda.io/miniconda.html or https://www.anaconda.com/download/.  Unless you are interested in the full package the miniconda installation should be sufficient for installing the rest of the modules needed for the PDielec installation.

As an example here is an installation using a Python 2.7 environment.  Using 'pip install' means that the anaconda build system is bypassed.ld system in Anaconda normally makes sure only compatible software versions are installed, but this seems to cause problems in this case.  ::

   conda create  --name pdielec2 python=2
   conda install --name pdielec2 scipy
   source activate pdielec2
   pip install matplotlib
   pip install termcolor
   pip install xlsxwriter
   pip install openpyxl
   pip install pyyaml
   pip install imageio

Here is an example of creating a Python 3 environment in Anaconda::

   conda create  --name pdielec3 python=3
   conda install --name pdielec3 -c conda-forge scipy
   source activate pdielec3
   pip install matplotlib
   pip install termcolor
   pip install xlsxwriter
   pip install openpyxl
   pip install pyyaml
   pip install imageio


Install Python (only needed if not using anaconda)
..................................................

* From https://www.python.org/downloads/windows/ download and run the Windows x86-64 executable installer for the latest Python 3.6 version
* Uncheck install launcher for all users
* Check "Add Python 3.6 to PATH"
* Click on "Install now" button
* Check installation ran OK by running the Idle python environment

Open a git bash console and type;::

  pip install numpy
  pip install scipy
  pip install matplotlib
  pip install PyYaml
  pip install XlsxWriter
  pip install openpyxl
  pip install termcolor
  pip install imageio


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
  make test_pdielec

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

PDielec directory structure
===========================

* PDielec/ is the home directory and contains the `pdielec` and `preader` commands
* PDielec/Python holds the source for the modules used by the pdielec and preader commands
* PDielec/Python/PyMieScatt holds a customised version of the PyMieScatt code downloaded from https://github.com/bsumlin/PyMieScatt.  Thanks to Benjamin Sumlin for making this available
* PDielec/Python/GUI holds the python code for PDGui
* PDielec/Examples a set of examples are available for Abinit, Crystal14, CASTEP, GULP, Phonopy, Mie and VASP.  Each example directory holds the input files to the QM/MM program and the relevant output files which are post processed by PDielec.  For each program there is also a preader directory which holds test output for the preader command.
* PDielec/Sphinx holds the documentation as restructured text documents (.rst).  Sphinx can be used to build the documentation in either html or pdf format.
* PDielec/docs holds the final html documentation.

Examples
========

Each example directory has the relevant input data sets use to run the QM/MM program and the output files from that run which are post-processed by PDielec.  There is a file `command.sh` which contains a typical example of a PDielec command line and which has been used to create the reference output files `command.ref.out` and `csvfile.ref.csv`.  The example can be run::

 bash command.sh

The output can be compared with the reference data to see if the program is working correctly.

The main Examples/ directory also has a  Makefile file which can be used to verify the correct working of the package.  Simply by typing `make` in the Examples directory each example will be run automatically and the output compared with the reference files.  To remove the intermediate files after running the tests automatically, type `make clean`.

