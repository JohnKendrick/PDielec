#!/usr/bin/env python
#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""A Makefile style system for testing and installing PDGui and PDielec.

This script accepts various command-line arguments to perform different actions such as running tests, benchmarks, cleaning build directories, and installing scripts.

Command line options
--------------------
    - `test` or `tests`: Runs all test cases.
    - `test-singlecrystal`: Runs tests for single crystal analysis.
    - `test-preader`: Runs tests for the preader component.
    - `test-p2cif`: Runs tests for the p2cif component.
    - `test-pdgui`: Runs tests for the PDGUI component.
    - `test-vibanalysis`: Runs tests for vibrational analysis.
    - `benchmarks`: Executes benchmark tests.
    - `pypi`: Runs commands for PyPI distribution.
    - `clean`: Cleans up directories.
    - `install`: Installs necessary scripts.
    - `scripts` or `-scripts`: Specifies the directory for scripts (expects a directory path next).
    - `-root` or `--root`: Specifies the root directory to operate in.
    - `-debug` or `--debug` or `-d`: Enables debug mode.
    - `-usesystem`: Use system installed binaries rather than local ones.
    - `-view` or `--view` or `-v`: Enables viewing mode.
    - `-regenerate`: Forces regeneration of test data.
    - `-padding`: Sets padding for some operation (expects an integer value next).
    - `-directory` or `--directory`: Sets the directory for some operation (expects a directory path next).
    - `<filename>.pdmake`: Specifies a PDMakefile to execute.
"""

import contextlib
import functools
import os
import subprocess
import sys
import time
from shutil import copyfile

import numpy as np

from PDielec.checkcsv import main as main_checkcsv
from PDielec.checkexcel import main as main_checkexcel

# Force a flush on print, needed for Windows
print = functools.partial(print, flush=True)

settings = {}
settings["padding"] = 50
settings["title"] = "title"

debug = False
viewing = False
useLocal = True

test_preader = [
	"Castep/preader",
	"Vasp/preader",
	"Phonopy/preader",
	"Gulp/preader",
	"Crystal/preader",
	"AbInit/preader",
	"QE/preader",
    ]

test_p2cif = [
	"P2Cif/",
    ]

test_pdgui = [
	"ATR/AlAs",
	"ATR/Na2SO42",
	"ATR/Na2SO42_fit",
	"AbInit/AlAs",
	"AbInit/BaTiO3",
	"AbInit/Na2SO42",
	"Castep/AsparticAcid",
	"Castep/Bubbles",
	"Castep/Castep17",
	"Castep/Isoleucine",
	"Castep/MgO",
	"Castep/Na2SO42",
	"Crystal/Leucine",
	"Crystal/Na2SO42",
	"Crystal/Na2SO42_C17",
	"Crystal/Quartz",
	"Crystal/ZnO/CPHF",
	"Crystal/ZnO/Default",
	"Crystal/ZnO/NoEckart",
	"Experiment/AlN",
	"Experiment/Forsterite",
	"Experiment/constant",
	"Experiment/drude-lorentz",
	"Experiment/fpsq",
	"Experiment/interpolation",
	"Experiment/Mayerhofer",
	"Experiment/Sapphire",
	"Gulp/Na2SO42",
	"Gulp/calcite",
	"Mie/MgO",
	"Mie/MgO_lognormal",
	"Phonopy/Al2O3/Primitive_cell",
	"Phonopy/Al2O3/Standard_cell",
	"Phonopy/ZnO",
	"Phonopy/Na2SO42",
	"Phonopy/QE",
	"Phonopy/Vasp",
	"QE/Cocaine",
	"QE/Na2SO42",
	"QE/Na2SO42-v7",
	"QE/Urea",
	"QE/ZnO",
	"SingleCrystal/Bi2Se3",
	"SingleCrystal/Bi2Se3_film",
	"SingleCrystal/L-alanine",
	"SingleCrystal/MgO",
	"SingleCrystal/ScatteringMatrix",
	"SizeEffects/BaTiO3",
	"SizeEffects/MgO",
	"SizeEffects/ZnO",
	"Vasp/F-Apatite",
	"Vasp/Na2SO42",
	"Vasp/Na2SO42_v",
	"Vasp/Vasp6",
	"Vasp/Urea",
	"Vasp/ZnO",
    ]

test_singlecrystal = [
	"Experiment/fpsq",
	"Experiment/Mayerhofer",
	"SingleCrystal/Bi2Se3",
	"SingleCrystal/Bi2Se3_film",
	"SingleCrystal/L-alanine",
	"SingleCrystal/MgO",
	"SingleCrystal/ScatteringMatrix",
    ]

test_vibanalysis = [
	"VibAnalysis/AsparticAcid",
	"VibAnalysis/BaTiO3",
	"VibAnalysis/Cocaine",
	"VibAnalysis/F-Apatite",
	"VibAnalysis/Isoleucine",
	"VibAnalysis/Leucine",
    ]

benchmarks = [
	"Castep/AsparticAcid",
	"Castep/Isoleucine",
	"Castep/MgO",
	"Castep/Bubbles",
	"Castep/Na2SO42",
	"Vasp/F-Apatite",
	"Vasp/ZnO",
	"Vasp/Na2SO42",
	"Gulp/calcite",
	"Gulp/Na2SO42",
	"Crystal/Na2SO42",
	"Crystal/Leucine",
	"Crystal/Quartz",
	"Crystal/ZnO/Default",
	"Crystal/ZnO/CPHF",
	"Crystal/ZnO/NoEckart",
	"AbInit/AlAs",
	"AbInit/BaTiO3",
	"AbInit/Na2SO42",
	"QE/ZnO",
	"QE/Na2SO42",
	"QE/Cocaine",
	"Phonopy/ZnO",
	"Phonopy/Na2SO42",
	"SizeEffects/BaTiO3",
	"SizeEffects/MgO",
	"SizeEffects/ZnO",
	"Mie/MgO",
	"Mie/MgO_lognormal",
	"ATR/AlAs",
	"ATR/Na2SO42",
    ]

def usage():
    """Print usage instructions for the command line tool.

    This function outputs instructions to `sys.stderr` for using the command line tool, including available commands and options. It also handles platform-specific instructions for Linux systems.

    SystemExit
        This function will terminate the program after displaying the usage message.
    """    
    print("pdmake:", file=sys.stderr)
    print("  test" , file=sys.stderr)
    print("  tests" , file=sys.stderr)
    print("           run all the tests               " , file=sys.stderr)
    print("  test-preader" , file=sys.stderr)
    print("           run all the preader tests" , file=sys.stderr)
    print("  test-pdgui" , file=sys.stderr)
    print("           run all the pdgui tests" , file=sys.stderr)
    print("  test-singlecrystal" , file=sys.stderr)
    print("           run all the singlecrystal tests" , file=sys.stderr)
    print("  test-p2cif" , file=sys.stderr)
    print("           run all the p2cif tests" , file=sys.stderr)
    print("  test-vibanalysis" , file=sys.stderr)
    print("           run all the vibanalysis tests" , file=sys.stderr)
    print("  benchmarks      " , file=sys.stderr)
    print("           run the benchmarks   " , file=sys.stderr)
    print("  --regenerate      " , file=sys.stderr)
    print("           regenerate all the test reference data" , file=sys.stderr)
    print("  --root directory  " , file=sys.stderr)
    print("           sets the root directory for the installation" , file=sys.stderr)
    print("  --directory", file=sys.stderr)
    print("           Use the directory when listing the tests" , file=sys.stderr)
    print("  --padding 60", file=sys.stderr)
    print("           change the width of the output title field" , file=sys.stderr)
    print("  --view  " , file=sys.stderr)
    print("           run pdgui with windowing on" , file=sys.stderr)
    print("  --usesystem " , file=sys.stderr)
    print("           Force the use of system executables" , file=sys.stderr)
    print("  --debug " , file=sys.stderr)
    print("           run with debug options switched on" , file=sys.stderr)
    if os.name != "nt":
        print("  pypi      " , file=sys.stderr)
        print("           Only use on linux installations...          " , file=sys.stderr)
        print("           prepare for pypi uploading" , file=sys.stderr)
        print("  clean" , file=sys.stderr)
        print("           Only use on linux installations...          " , file=sys.stderr)
        print("           clean up" , file=sys.stderr)
        print("  install" , file=sys.stderr)
        print("           Only use on linux installations...          " , file=sys.stderr)
        print("           install the executables to the scripts directory" , file=sys.stderr)
        print("  --scripts directory" , file=sys.stderr)
        print("           Only use on linux installations...          " , file=sys.stderr)
        print("           Set the scripts directory to be used for installation" , file=sys.stderr)
        print("           the default is ~/bin                       " , file=sys.stderr)
    sys.exit()

def checkLocalExecutables():
    """Check to see if the rootDirectory has the executables we need.

    Parameters
    ----------
    rootDirectory : str
        The path to the root directory to check.

    Returns
    -------
    bool
        True if the executables are found, False otherwise.

    """
    global rootDirectory
    os.chdir(rootDirectory)
    result = True
    if not os.path.isfile("preader"):
        result =  False
    if not os.path.isfile("p2cif"):
        result =  False
    if not os.path.isfile("pdgui"):
        result =  False
    if not os.path.isfile("vibanalysis"):
        result =  False
    return result

def install(scripts):
    """Install a set of scripts into a specified directory.

    Parameters
    ----------
    scripts : str
        The directory where the scripts will be copied to.

    Notes
    -----
    This function copies a predefined list of scripts and directories into 
    the specified `scripts` directory. It utilizes shell commands to perform the 
    copy operation and to create necessary subdirectories. Global variable 
    `rootDirectory` should be set before calling this function as it is used to 
    determine the source directory of the scripts. The installation progress is 
    printed to the console.

    Raises
    ------
    This function does not explicitly raise exceptions, but subprocess.run 
    may raise exceptions related to the execution of shell commands.

    """    
    print("Performing installation of scripts into ",scripts)
    global rootDirectory
    mychdir(rootDirectory)
    print("Installing preader",end="")
    subprocess.run(f"cp -P preader {scripts}",shell=True)
    print(",p1reader",end="")
    subprocess.run(f"cp -P p1reader {scripts}",shell=True)
    print(",pdgui",end="")
    subprocess.run(f"cp -P pdgui {scripts}",shell=True)
    print(",p2cif",end="")
    subprocess.run(f"cp -P p2cif {scripts}",shell=True)
    print(",pdcompare",end="")
    subprocess.run(f"cp -P pdcompare {scripts}",shell=True)
    print(",graphdatagenerator",end="")
    subprocess.run(f"cp -P graphdatagenerator {scripts}",shell=True)
    print(",vibanalysis",end="")
    subprocess.run(f"cp -P vibanalysis {scripts}",shell=True)
    print(",pdmake",end="")
    subprocess.run(f"cp -P pdmake {scripts}",shell=True)
    print(",python scripts")
    subprocess.run(f"mkdir -p {scripts}/PDielec",shell=True)
    subprocess.run(f"mkdir -p {scripts}/PDielec/GUI",shell=True)
    subprocess.run(f"cp -r PDielec/*.py {scripts}/PDielec",shell=True)
    subprocess.run(f"cp -r PDielec/GUI/*.py {scripts}/PDielec/GUI/",shell=True)
    subprocess.run(f"cp -r PDielec/GUI/*.png {scripts}/PDielec/GUI/",shell=True)
    print("Finished installation of scripts into ",scripts)

@contextlib.contextmanager
def redirect(file):
    # capture all outputs to a log file while still printing it
    """Redirect the standard output to a specified file.

    This function is a context manager that redirects `sys.stdout` to a file,
    optionally suppressing stderr output if not in debug mode. It is useful for
    capturing the output of a block of code to a file.

    Parameters
    ----------
    file : str
        The path to the file where the stdout will be redirected.

    Yields
    ------
    _io.TextIOWrapper
        A writable stream associated with the log file.

    Notes
    -----
    `sys.stderr` is redirected to `os.devnull` if `debug` is `False`, effectively
    suppressing any error messages unless in debug mode. Ensure that debug status is
    properly set before using this context manager.

    Upon exiting the `with` block, stdout and stderr are restored to their original
    state.

    Raises
    ------
    AttributeError
        If an attribute not supported by the `sys.stdout` is being accessed through
        the Logger class.

    """    
    class Logger:
        def __init__(self, file):
            fd = open(file,'w')             # noqa SIM115
            self.terminal = sys.stdout
            self.log = fd

        def write(self, message):
            global debug
            if debug:
                self.terminal.write(message)
            self.log.write(message)

        def __getattr__(self, attr):
            return getattr(self.terminal, attr)

    logger = Logger(file)

    sys.stdout = logger
    if not debug:
        sys.stderr = open(os.devnull,'w')    # noqa SIM115
    try:
        yield logger.log
    finally:
        sys.stdout = sys.__stdout__
        if not debug:
            sys.stderr = sys.__stderr__

def readNmaFile(file):
    """Read and extract data from an NMA file.

    Parameters
    ----------
    file : str
        The path to the NMA file to be read.

    Returns
    -------
    tuple
        A tuple containing three elements:
        - The sum of all frequencies excluding the first three.
        - The sum of all R squared values excluding the first three.
        - The sum of all explained variance values excluding the first three.

    Notes
    -----
    This function reads an NMA (Normal Mode Analysis) file, looking for lines that contain
    'Mode', 'R**2', or 'Explained'. For lines starting with 'Mode', it extracts the
    frequency values, excluding the first three modes which are typically zero or negligible.
    For 'R**2' and 'Explained' lines, it extracts their respective values under the same conditions.
    All extracted values are then summed up and returned.

    Raises
    ------
    FileNotFoundError
        If the file specified by the 'file' parameter does not exist.

    """    
    nmodes = 0
    frequencies = []
    rsquared = []
    explained = []
    with open(file) as fd:
        line = fd.readline()
        while line:
            splits = line.split()
            if len(splits) > 0:
                if splits[0] == "Mode":
                    if splits[1] != "too":
                        nmodes += 1
                        frequencies.append(float(splits[2]))
                elif splits[0] == "R**2":
                    rsquared.append(float(splits[2]))
                elif splits[0] == "Explained":
                    explained.append(float(splits[3]))
                # end if
            # end if
            line = fd.readline()
        #end while
    #end with open
    f = np.array(frequencies[3:])
    r = np.array(rsquared[3:])
    e = np.array(explained[3:])
    return np.sum(f),np.sum(r),np.sum(e)

def compareNmaFiles(file1,file2):
    """Compare two NMA files for differences within a tolerance.

    Parameters
    ----------
    file1 : str
        The path to the first NMA file.
    file2 : str
        The path to the second NMA file.

    Returns
    -------
    int
        The number of differences detected in frequency, redox, and energy 
        values between the two NMA files.

    Notes
    -----
    This function reads two NMA files using the `readNmaFile` function, 
    which returns the frequency (f), redox (r), and energy (e) values 
    for each file. The comparison is done based on a tolerance of 1e-06, 
    and the total number of differences in these values are counted and returned.

    """    
    f1,r1,e1 = readNmaFile(file1)
    f2,r2,e2 = readNmaFile(file2)
    nerrors = 0
    if abs(f1 - f2) > 1.0e-6:
        nerrors += 1
    if abs(r1 - r2) > 1.0e-6:
        nerrors += 1
    if abs(e1 - e2) > 1.0e-6:
        nerrors += 1
    return nerrors

def compareFiles(file1,file2):
    """Compare two text files for differences.

    This function opens and reads two text files line by line and then compares them word by word,
    ignoring certain punctuation marks. If any differences are found in either numeric values
    (floating-point comparison) or text strings, they are counted as errors. Numeric differences
    are considered errors if the relative difference between them is more than 0.0001, provided
    both numbers are not insignificantly small. The total count of such errors is returned.
    Optionally, if the global variable 'debug' is set to True, details of each difference are printed.

    Parameters
    ----------
    file1 : str
        The path to the first file to be compared.
    file2 : str
        The path to the second file to be compared.

    Returns
    -------
    int
        The total number of errors (differences) found.

    Notes
    -----
    - The function makes use of a global variable 'debug' which, if set to True,
      enables printing detailed information about each found difference directly to the console.
    - It's assumed that both files can be opened and read successfully. The function 
      does not explicitly handle I/O errors, such as a file not existing.
    - Special characters '(', ')', '%', '/', and '\\' are ignored during the comparison.
    - For numerical comparison, the function employs a relative difference check to
      determine if two numbers should be considered different, taking into account
      floating-point inaccuracies.

    """    
    global debug
    with open(file1) as fd1, open(file2) as fd2:
        lines1 = fd1.readlines()
        lines2 = fd2.readlines()
        nerrors = 0
        for line1,line2 in  zip(lines1,lines2):
            for word1,word2 in zip(line1.split(),line2.split()):
                word1 = word1.replace("(","")
                word1 = word1.replace(")","")
                word1 = word1.replace("%","")
                word1 = word1.replace("/","")
                word1 = word1.replace("\\","")
                word2 = word2.replace("(","")
                word2 = word2.replace(")","")
                word2 = word2.replace("%","")
                word2 = word1.replace("/","")
                try:
                    float1 = float(word1)
                    float2 = float(word2)
                    if abs(float1)+abs(float2) > 1.e-12 and 2*abs(float1-float2)/(abs(float1)+abs(float2)) > 1.0e-4:
                        if debug:
                            print("Float difference",float1,float2)
                            print("Line 1",line1)
                            print("Line 2",line2)
                        nerrors += 1
                except Exception:
                    if word1 != word2:
                        if debug:
                            print("Word difference",word1,word2)
                            print("Line 1",line1)
                            print("Line 2",line2)
                        nerrors += 1
    return nerrors


def runP2CifTest(title, instructions, regenerate):
    """Run a test on p2cif.

    Parameters
    ----------
    title : str
        The title of the test
    instructions : str
        The command line parameters.
    regenerate : bool
        Set to True if the reference file is overwritten.

    """
    global rootDirectory
    global useLocal
    print(title,end="",flush=True)
    if useLocal:
        sys.argv = ["python"]
        sys.argv.append(os.path.join(rootDirectory,"p2cif"))
    else:
        sys.argv = ["p2cif"]
    sys.argv.extend(instructions)
    outputfile = "all.cif"
    if os.path.exists("all.cif"):
        os.remove("all.cif")
    if regenerate:
        outputfile = "all.ref.cif"
    with open(outputfile,"w") as stdout:
        if debug:
            subprocess.run(sys.argv,stdout=stdout)
        else:
            subprocess.run(sys.argv,stdout=stdout,stderr=subprocess.DEVNULL)
        # end if debug
    # end with open
    # If not doing a regeneration perform a check
    if not regenerate:
        nerrors = compareFiles("all.cif", "all.ref.cif")
        if nerrors > 0:
            print(f" {nerrors} ERRORS:")
        else:
            print(" OK:")
        # end if
    else:
        print(" Regenerated:")
    # end if
    return

def runVibAnalysis(title, instructions, regenerate):
    """Run a vibanalysis test.

    Parameters
    ----------
    title : str
        The title of the test
    instructions : str
        The command line parameters.
    regenerate : bool
        Set to True if the reference file is overwritten.

    """
    global rootDirectory
    global useLocal
    print(title,end="",flush=True)
    if useLocal:
        sys.argv = ["python"]
        sys.argv.append(os.path.join(rootDirectory,"vibanalysis"))
    else:
        sys.argv = ["vibanalysis"]
    sys.argv.extend(instructions)
    filename = None
    for option in instructions:
        if not option.startswith("-"):
            filename = option
        # end if
    # end for
    if filename is None:
        print("Error in runVibAnalysis: there is no filename specified")
    header = os.path.splitext(filename)[0]
    nmafile = header+".nma"
    reffile = header+".nma.ref"
    if debug:
        subprocess.run(sys.argv)
    else:
        subprocess.run(sys.argv,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    # If not doing a regeneration perform a check
    if not regenerate:
        nerrors = compareNmaFiles(nmafile, reffile)
        if nerrors > 0:
            print(f" {nerrors} ERRORS:")
        else:
            print(" OK:")
        # end if
    else:
        copyfile(nmafile,reffile)
        print(" Regenerated:")
    # end if
    return

def runPreaderTest(title, instructions, regenerate):
    """Run a preader test.

    Parameters
    ----------
    title : str
        The title of the test
    instructions : str
        The command line parameters.
    regenerate : bool, optional
        Set to true if the reference file is overwritten. Default is False.

    """
    global debug
    global rootDirectory
    global useLocal
    global debug
    if debug:
        print("runPreaderTest:",title,instructions)
    print(title,end="",flush=True)
    if useLocal:
        sys.argv = ["python"]
        sys.argv.append(os.path.join(rootDirectory,"preader"))
    else:
        sys.argv = ["preader"]
    sys.argv.extend(instructions)
    if os.path.exists("command.csv"):
        os.remove("command.csv")
    outputfile = "command.csv"
    if regenerate:
        outputfile = "command.ref.csv"
    with open(outputfile,"w") as fd:
        if debug:
            result = subprocess.run(sys.argv,stdout=fd)
        else:
            result = subprocess.run(sys.argv,stdout=fd,stderr=subprocess.DEVNULL)
    # end with
    # If not doing a regeneration perform a check
    if not regenerate:
        sys.argv = ["checkcsv", "command.ref.csv","command.csv"]
        with redirect(os.devnull):
            result = main_checkcsv()
        nerrors,keep_line_number,keep_word1,keep_word2,max_percentage_error = result
        if nerrors > 0:
            print(" ERRORS:"+f"LARGEST ON LINE {nerrors} - max %error={keep_line_number}")
        else:
            print(" OK:"+f" - max %error={max_percentage_error}")
        # end if
    else:
        print(" Regenerated:")
    # end if
    return

def runPDGuiTest(title, instructions, regenerate, benchmarks=False):
    """Run a pdgui test.

    Parameters
    ----------
    title : str
        The title of the test
    instructions : str
        The command line parameters.
    regenerate : bool
        Set to true if the reference file is overwritten.
    benchmarks : boolean
        Set to true if benchmarks are to be run

    """
    global viewing
    global start_time
    global debug
    global rootDirectory
    global useLocal
    global settings
    print(title,end="",flush=True)
    if useLocal:
        sys.argv = ["python"]
        sys.argv.append(os.path.join(rootDirectory,"pdgui"))
    else:
        sys.argv = ["pdgui"]
    if os.path.exists("results.xlsx"):
        os.remove("results.xlsx")
    if os.path.exists("results.csv"):
        os.remove("results.csv")
    # Append any threading or cpu information
    sys.argv.extend(settings["PDGui tokens"])
    sys.argv.extend(instructions)
    if debug:
        result = subprocess.run(sys.argv)
    else:
        result = subprocess.run(sys.argv,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    if not viewing:
        # If not doing a regeneration perform a check
        if not regenerate and not benchmarks:
            sys.argv = ["checkexcel", "results.ref.xlsx","results.xlsx"]
            with redirect(os.devnull):
                result = main_checkexcel()
            nerrors,row,col,sheet,file1,value1,file2,value2,max_percentage_error = result
            if nerrors > 0:
                print(f"{nerrors} ERRORS:"+f"{sheet}@{row},{col} - max %error={max_percentage_error}")
            else:
                print(" OK:"+f" - max %error={max_percentage_error}")
            # end if
        elif benchmarks:
            end_time = time.time()
            elapsed_time = end_time - start_time
            start_time = end_time
            print(" OK:"+f" - elapsed time {elapsed_time:.3f}s")
        elif not benchmarks:
            # If we asked for a benchmarking then don't do a regenerate
            copyfile("results.xlsx","results.ref.xlsx")
            print(" Regenerated:")
        # end if not regenerate
    #end if not viewing
    return

def changePadding(all):
    """Change the padding based on the maximum length of elements in a list or the predefined setting.

    Parameters
    ----------
    all : list
        A list containing elements whose lengths are to be compared.

    Returns
    -------
    int
        The maximum length found in the list or the predefined padding value, whichever is greater.

    Notes
    -----
    This function assumes the presence of a global dictionary named `settings` with a key 'padding'
    that stores the default padding value. If the maximum length of an element in the list is
    less than 'settings['padding']', the function returns 'settings['padding']'. Otherwise,
    it returns the maximum length found.

    """    
    global settings
    maxlen = 0
    for d in all:
        maxlen = max(maxlen,len(d))
    return max(maxlen, settings["padding"])

def runTests(testlist, testType, regenerate):
    """Run the tests given in the directories stored in testlist.

    Parameters
    ----------
    testlist : list of str
        The directories where the tests are stored.
    testType : {'pdgui', 'preader', 'p2cif', 'benchmark'}
        The type of tests to run.
    regenerate : bool
        If True, the test data will be regenerated.

    """
    global start_time
    global debug
    mychdir(os.path.join(rootDirectory,"Examples"))
    homedir = os.getcwd()
    convertTestType={}
    convertTestType["pdgui"] = "PDGui tests"
    convertTestType["vibanalysis"] = "VibAnalysis tests"
    convertTestType["preader"] = "PReader tests"
    convertTestType["p2cif"] = "P2Cif tests"
    convertTestType["benchmarks"] = "BenchMarks"
    convertTestType["singlecrystal"] = "Single Crystals"
    print("")
    print("--------------------------------------------------")
    print(convertTestType[testType],"starting" )
    start_time = time.time()
    test_start_time = start_time
    if debug:
        print("RunTests: homedir=",homedir)
    pdmakefile = "command.pdmake"
    for directory in testlist:
        if debug:
            print("RunTest: directory=",directory)
        if not os.path.isdir(directory):
            print("Error: command needs to be executed in the PDielec home directory")
            print(f"       current directory is {homedir}")
            print(f"       required directory is {directory}")
            sys.exit()
        if testType == "benchmarks":
            runPdMakefile(directory,pdmakefile,regenerate,benchmarks=True)
        else:
            runPdMakefile(directory,pdmakefile,regenerate)
    # end for
    elapsed_time = time.time() - test_start_time
    print("--------------------------------------------------")
    print(convertTestType[testType],f"completed in {elapsed_time:.3f}s")

def readPdMakefile(directory,filename):
    # The first word is the dictionary key, the rest is its value
    """Read a custom makefile and return the title and instructions.

    Parameters
    ----------
    directory : str
        The directory in which the makefile is located. If the directory name is
        provided and the `settings['title']` is not 'title', the directory name is used as the title.
        If the directory is not provided (empty string), the current working directory is used as the title.
    filename : str
        The name of the makefile to be read.

    Returns
    -------
    tuple
        A tuple returning a dictionary of instructions
        - title (str): The title for the makefile, determined based on the `settings['title']` value or the directory name, or the current working directory.
        - instructions (dict): A dictionary where each key is an instruction extracted from the makefile, and the value is a list of arguments for that instruction.

    Raises
    ------
    FileNotFoundError
        If the specified file cannot be found or read.

    Notes
    -----
    - This function uses two global variables:
        - `settings`: A dictionary that contains various settings, including 'title' and 'padding'.
          'title' is used to determine if the title should be extracted from the makefile's first line,
          and 'padding' determines how much padding the title should have when displayed.
        - `debug`: A boolean flag used to control debugging output. If `True`, additional information
          about the function's operation will be printed to the standard output.

    - The function removes line continuation characters and treats the continued lines as a single line.

    - Lines starting with `#` are treated as comments and ignored.

    - It's assumed that `os` module has been imported for `os.getcwd()` to work, though it is not explicitly
      stated in the function definition.

    """    
    global settings
    global debug
    instructions = {}
    if debug:
        print("ReadPdMakefile:",directory,filename)
    with open(filename) as fd:
        line = fd.readline()[:-1]
        if settings["title"] == "title":
            title = line.ljust(settings["padding"])
        elif len(directory)>1:
            title = directory.ljust(settings["padding"])
        else:
            title = os.getcwd().ljust(settings["padding"])
        line = fd.readline()[:-1]
        while line:
            if len(line)>0 and not line.startswith("#"):
                full_line = line.replace("\\","")
                while line.endswith("\\"):
                    line = fd.readline()[:-1]
                    full_line += line.replace("\\","")
                splits = full_line.split()
                key = splits[0]
                value = splits[1:]
                instructions[key] = value
            # end if
            line = fd.readline()[:-1]
        # end while
    #end open
    return title,instructions

def mychdir(directory):
    """Local version of os.chdir with some testing.

    Notes
    -----
    This is a local implementation of os.chdir which includes
    additional testing or functionality tailored to specific needs.

    """
    if directory != "":
        if os.path.isdir(directory):
            os.chdir(directory)
        else:
            print("mychdir: Error directory does not exist",directory)
            print("mychdir: Current directory is          ",os.getcwd())
        # end if
    # end if
    return

def runPdMakefile(directory,pdmakefile,regenerate,benchmarks=False):
    """Run a specific pdMakefile.

    Parameters
    ----------
    directory : str
        Directory where the pdmakefile will be run
    pdmakefile : str
        Name of the pdmakefile
    regenerate : bool
        True if this is a regeneration of data
    benchmarks : bool (Option default is false)
        True if this is a benchmark

    Notes
    -----
    The routine reads the pdmakefile and determines the type of pdmakefile and the instructions
    Then the instructions are processed according to the context (preader, pdgui, p2cif or vibanalysis)

    """
    global debug
    homedir = os.getcwd()
    if debug:
        print("RunPdMakefile: directory =",directory)
        print("RunPdMakefile: pdmakefile =",pdmakefile)
        print("RunPdMakefile: cwd =",homedir)
    mychdir(directory)
    title,instructions = readPdMakefile(directory,pdmakefile)
    if debug:
        print("RunPdMakefile: title =",title)
        print("RunPdMakefile: instructions =",instructions)
    for key in instructions:
        parameters = instructions[key]
        if key == "preader":
            runPreaderTest(title,parameters,regenerate)
        elif key == "vibanalysis":
            runVibAnalysis(title,parameters,regenerate)
        elif key == "p2cif":
            runP2CifTest(title,parameters,regenerate)
        elif key == "pdgui":
            if viewing:
                parameters.extend(["-script", "script.py"])
            elif benchmarks:
                parameters.extend(["-nosplash", "-exit", "-script", "script.py"])
            else:
                parameters.extend(["-nosplash", "results.xlsx", "-exit", "-script", "script.py"])
            if debug:
                parameters.append("-d")
            runPDGuiTest(title,parameters,regenerate,benchmarks=benchmarks)
    mychdir(homedir)

def runClean():
    """Clean old results from the Examples directory.

    This function navigates to the root directory and removes specific file types related to old results. Supported file types for deletion include .xlsx, .csv, .cif, and .nma files. This operation is not supported on Windows platforms and will return immediately if attempted.

    Notes
    -----
    - The function checks the operating system and will not execute on Windows systems (`os.name == 'nt'`).
    - It navigates to a directory specified by a global variable `rootDirectory`.
    - Uses `subprocess.run` with `shell=True` to execute shell commands for deleting files. This operation involves running `find` commands to locate and remove files with specified extensions.
    - It assumes that the `rootDirectory` global variable and necessary imports (`os`, `subprocess`) are correctly defined elsewhere.

    Raises
    ------
    None

    """    
    if os.name == "nt":
        print("Unable to clean installation in Windows")
        return
    print("Cleaning old results from the Examples directory")
    os.chdir(rootDirectory)
    subprocess.run("find . -name results.xlsx -exec rm -f {} \\;",shell=True)
    subprocess.run("find . -name results.csv -exec rm -f {} \\;",shell=True)
    subprocess.run("find . -name command.csv -exec rm -f {} \\;",shell=True)
    subprocess.run("find . -name all.cif -exec rm -f {} \\;",shell=True)
    subprocess.run("find . -name \\*.nma -exec rm -f {} \\;",shell=True)
    subprocess.run("find . -name \\*.pyc -exec rm -f {} \\;",shell=True)
    print("Cleaning complete")

def runPyPi():
    """Create PyPi distribution files, restricted by OS.

    This function changes the current working directory to `rootDirectory` 
    and runs a series of shell commands that clean up previous distribution files and 
    create new PyPi distribution files. The function is not designed to run on Windows OS.
    """    
    if os.name == "nt":
        print("Unable to create PyPi installation in Windows")
        return
    print("Creating PyPi distribution files")
    os.chdir(rootDirectory)
    # deprecated subprocess.run('rm -rf build dist PDielec.egg-info; python setup.py sdist bdist_wheel',shell=True)
    # Following command needs an installed copy of setuptools
    # pip install setuptools
    subprocess.run("rm -rf build dist PDielec.egg-info; python -m build ",shell=True)

def testForRootDirectory(path):
    """Test if a given path contains a subdirectory named 'Examples'.

    This function checks if the given path has a direct subdirectory
    called 'Examples' and returns a boolean indicating the result.

    Parameters
    ----------
    path : str
        The directory path to test for the existence of a 'Examples' subdirectory.

    Returns
    -------
    bool
        Returns True if the 'Examples' subdirectory exists within the given path,
        otherwise returns False.

    """    
    test = os.path.join(path,"Examples")
    return os.path.isdir(test)

def findRootDirectory(start):
    """Find the root directory starting from start.

    Parameters
    ----------
    start : str
        The starting directory path.

    Returns
    -------
    str
        The root directory path.

    """
    lastpath = start
    path = os.path.join(start,"Examples")
    while path != lastpath:
        lastpath = path
        path,tail=os.path.split(path)
        if testForRootDirectory(path):
            return path
    print("Unable to locate root directory")
    sys.exit()

def main():
    """Process command line arguments and execute actions (main routine).

    This script accepts various command-line arguments to perform different actions such as running tests, benchmarks, cleaning build directories, and installing scripts.

    Command line options
    --------------------
    - `test` or `tests`: Runs all test cases.
    - `test-singlecrystal`: Runs tests for single crystal analysis.
    - `test-preader`: Runs tests for the preader component.
    - `test-p2cif`: Runs tests for the p2cif component.
    - `test-pdgui`: Runs tests for the PDGUI component.
    - `test-vibanalysis`: Runs tests for vibrational analysis.
    - `benchmarks`: Executes benchmark tests.
    - `pypi`: Runs commands for PyPI distribution.
    - `clean`: Cleans up directories.
    - `install`: Installs necessary scripts.
    - `scripts` or `-scripts`: Specifies the directory for scripts (expects a directory path next).
    - `-root` or `--root`: Specifies the root directory to operate in.
    - `-debug` or `--debug` or `-d`: Enables debug mode.
    - `-usesystem`: Use system installed binaries rather than local ones.
    - `-view` or `--view` or `-v`: Enables viewing mode.
    - `-regenerate`: Forces regeneration of test data.
    - `-padding`: Sets padding for some operation (expects an integer value next).
    - `-directory` or `--directory`: Sets the directory for some operation (expects a directory path next).
    - `<filename>.pdmake`: Specifies a PDMakefile to execute.

    Raises
    ------
    SystemExit
        If an invalid command line argument is provided or certain conditions are not met (e.g., invalid root directory, command not executed in the PDielec home directory).

    Notes
    -----
    - It uses global variables for configuration options like `debug`, `viewing`, `settings`, `rootDirectory`, and `useLocal`.
    - Requires the `os`, `sys`, and specific project function imports to function.
    - The script changes the current working directory multiple times, which affects relative path usage.

    """    
    global debug
    global viewing
    global settings
    global rootDirectory
    global useLocal
    # Start processing the directories
    originalDirectory = os.getcwd()
    if len(sys.argv) <= 1 :
        usage()
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    regenerate = False
    scriptsDirectory="~/bin"
    rootDirectory = originalDirectory
    rootDirectory = findRootDirectory(rootDirectory)
    actions = []
    pdmakefiles = []
    pdgui_tokens = []
    # Loop over the tokens on the command line
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token in ( "test", "tests" ):
            actions.append("test all")
        elif token == "test-singlecrystal":
            actions.append("test singlecrystal")
        elif token == "test-preader":
            actions.append("test preader")
        elif token == "test-p2cif":
            actions.append("test p2cif")
        elif token == "test-pdgui":
            actions.append("test pdgui")
        elif token == "test-vibanalysis":
            actions.append("test vibanalysis")
        elif token == "benchmarks":
            actions.append("run benchmarks")
        elif token in ( "-root", "--root" ):
            itoken += 1
            rootDirectory = tokens[itoken]
        elif token in ( "-debug", "--debug", "-d" ):
            debug = True
        elif token == "-usesystem":
            useLocal = False
        elif token in ( "-view", "--view", "-v" ):
            viewing = True
        elif token in ( "-regenerate", "--regenerate" ):
            regenerate = True
        elif token in ( "-padding", "--padding" ):
            itoken += 1
            settings["padding"] = int(tokens[itoken])
        elif token in ( "-directory", "--directory" ):
            itoken += 1
            settings["title"] = "directory"
        elif token == "pypi":
            actions.append("run pypi")
        elif token == "clean":
            actions.append("run clean")
        elif token == "install":
            actions.append("install")
        elif token in ( "scripts", "-scripts" ):
            itoken += 1
            scriptsDirectory = tokens[itoken]
        elif token.endswith(".pdmake"):
            pdmakefiles.append(token)
        elif token in ( "-cpus", "--cpus" ):
            itoken += 1
            pdgui_tokens.append("-cpus")
            pdgui_tokens.append(tokens[itoken])
        elif token in ( "-threads", "--threads" ):
            itoken += 1
            pdgui_tokens.append("-threads")
            pdgui_tokens.append(tokens[itoken])
        elif token in ( "-threading", "--threading" ):
            pdgui_tokens.append("-threading")
        else:
            usage()
            sys.exit()
        # endif
    # end while
    settings["PDGui tokens"] = pdgui_tokens
    #
    # Change to the rootDirectory
    if not os.path.isdir(rootDirectory):
        print("Error: The root directory is not valid",rootDirectory)
        sys.exit()
    mychdir(rootDirectory)
    #
    # Check to see if the required executables are present
    # If the are then use them
    # 
    if useLocal:
        useLocal = checkLocalExecutables()
    # Change directory to the Examples directory
    if not os.path.isdir("Examples"):
        print("Error: command needs to be executed in the PDielec home directory")
        sys.exit()
    #
    # Handle any pdmake files there might be
    #
    for pdmakefile in pdmakefiles:
        directory = originalDirectory
        runPdMakefile(directory,pdmakefile,regenerate)
    #
    # Now loop over the actions
    #
    mychdir("Examples")
    if settings["title"] != "title":
        test_all = []
        test_all.extend(test_p2cif)
        test_all.extend(test_preader)
        test_all.extend(test_pdgui)
        test_all.extend(test_vibanalysis)
        settings["padding"] = changePadding(test_all)
    for action in actions:
        if action == "test all":
            runTests(test_p2cif      ,"p2cif"      ,regenerate)
            runTests(test_preader    ,"preader"    ,regenerate)
            runTests(test_pdgui      ,"pdgui"      ,regenerate)
            runTests(test_vibanalysis,"vibanalysis",regenerate)
        elif action == "test preader":
            runTests(test_preader    ,"preader"    ,regenerate)
        elif action == "test p2cif":
            runTests(test_p2cif      ,"p2cif"      ,regenerate)
        elif action == "test pdgui":
            runTests(test_pdgui      ,"pdgui"      ,regenerate)
        elif action == "test singlecrystal":
            runTests(test_singlecrystal      ,"singlecrystal"      ,regenerate)
        elif action == "test vibanalysis":
            runTests(test_vibanalysis,"vibanalysis",regenerate)
        elif action == "run benchmarks":
            runTests(benchmarks,"benchmarks",regenerate)
        elif action == "run pypi":
            runPyPi()
        elif action == "run clean":
            runClean()
        elif action == "install":
            install(scriptsDirectory)
        else:
            usage()
            sys.exit()
        #end if
    #end for

# end of def main

if __name__ == "__main__":
    main(sys)
