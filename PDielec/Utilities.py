#!/usr/bin/python
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
"""Utility Functions.

A set of utility functions that may be used anywhere in the package.
"""

import os
import sys

from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.ExperimentOutputReader import ExperimentOutputReader
from PDielec.GulpOutputReader import GulpOutputReader
from PDielec.PhonopyOutputReader import PhonopyOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.VaspOutputReader import VaspOutputReader


def printsp(name,matrix):
    """Print 4x4 matrices or 4 vectors.

    Parameters
    ----------
    name : string
        The name of the matrix
    matrix : 4x4 or 4 array
        The matrix to be printed

    Returns
    -------
    None

    Notes
    -----
    None

    """
    print("")
    print(name)
    if len(matrix.shape) == 1:
        columns = matrix.shape[0]
        string = ""
        for i in range(columns):
            string = string + f"{matrix[i]:+.5f}" + "   "
        print(string)
    else:
        rows = matrix.shape[0]
        columns = matrix.shape[1]
        for j in range(rows):
            string = ""
            for i in range(columns):
                string = string + f"{matrix[j,i]:+.5f}" + "   "
            print(string)
    return

def find_program_from_name( filename ):
    # Determine the program to use from the file name being used
    """Determine the simulation program from a given filename.

    Parameters
    ----------
    filename : str
        The complete path (absolute or relative) to a file.

    Returns
    -------
    (str1,str2)
        A tuple of 2 strings, the first, str1, is the program name that was used to calculate the frequencies.
        In the case of phonopy, the second string is the QM program that calculates the forces

    Notes
    -----
    This function examines the file extension and, in some cases, the presence of specific files in the same directory, to determine the associated simulation program. It recognizes files from several popular materials simulation programs, such as 'phonopy', 'gulp', 'vasp', and others.

    Examples
    --------
    ::

        program = find_program_from_name('./data/structure.castep')
        print(program)
        # Output: ("castep","")
    
        program = find_program_from_name('path/to/simulation/OUTCAR')
        print(program)
        # Output: ("phonopy","vasp")
        # Depends on the presence of 'phonopy.yaml'

    """    
    head,tail = os.path.split(filename)
    root,ext = os.path.splitext(tail)
    head = "./" if head == "" else head+"/"
    if tail == "OUTCAR":
        if os.path.isfile(head+"phonopy.yaml"):
            return "phonopy","vasp"
        else:
            return "vasp",""
    if ext == ".gout":
        return "gulp",""
    if ext == ".born":
        return "phonopy","vasp"
    if ext == ".castep":
            if os.path.isfile(head+"phonopy.yaml"):
                return "phonopy","castep"
            else:
                return "castep",""
    if ext ==  ".out":
        if os.path.isfile(head+root+".files"):
            if os.path.isfile(head+"phonopy.yaml"):
                return "phonopy","abinit"
            else:
                return "abinit",""
        elif os.path.isfile(head+root+".dynG"):
            if os.path.isfile(head+"phonopy.yaml"):
                return "phonopy","quantum espresso"
            else:
                return "quantum espresso",""
        elif os.path.isfile(head+"phonopy.yaml"):
            return "phonopy","crystal"
        else:
            return "crystal",""
    if ext ==  ".dynG":
        if os.path.isfile(head+"phonopy.yaml"):
            return "phonopy","quantum espresso"
        else:
            return "quantum espresso",""
    if ext ==  ".exp":
        return "experiment",""
    if ext ==  ".py":
        return "pdgui",""
    return "",""

def get_reader( name, program, qmprogram):
    """Get the appropriate output reader based on the simulation program and, if specified, the quantum mechanical program.

    This function is designed to create an output reader object for various simulation programs (like CASTEP, VASP, etc.) and, for phonopy simulations, it can additionally create a quantum mechanical output reader based on the specified quantum mechanical program.

    Parameters
    ----------
    name : str
        The primary filepath or name associated with the output file(s).
    program : str
        The name of the simulation program. Supported values are 'castep', 'vasp', 'gulp', 'crystal', 'abinit', 'qe', and 'phonopy'.
    qmprogram : str
        The name of the quantum mechanical program used in conjunction with phonopy. Supported values are the same as for `program`, but its relevance is exclusive to when `program` is set to 'phonopy'.

    Returns
    -------
    object
        An instance of the appropriate output reader class based on the input `program` and, if applicable, `qmprogram`.

    Notes
    -----
    For phonopy simulations, both phonon and quantum mechanical output files are considered. If the `qmprogram` is not compatible or not specified, the function will exit.

    Examples
    --------
    >>> reader = get_reader("output.log", "castep", "")
    >>> reader = get_reader("output", "phonopy", "vasp")

    """    
    program = program.lower()
    qmprogram = qmprogram.lower()
    fulldirname = name
    head,tail = os.path.split(fulldirname)
    root,ext = os.path.splitext(tail)
    if program == "castep":
        names = [ name ]
        reader = CastepOutputReader( names )
    elif program == "vasp":
        name1 = name
        name2 = os.path.join(head,"KPOINTS")
        names = [ name1, name2 ]
        reader = VaspOutputReader( names )
    elif program == "gulp":
        names = [ name ]
        reader = GulpOutputReader( names )
    elif program == "crystal":
        names = [ name ]
        reader = CrystalOutputReader( names )
    elif program == "abinit":
        names = [ name ]
        reader = AbinitOutputReader( names )
    elif program == "quantum espresso":
        tail1 = root+".dynG"
        tail2 = root+".log"
        tail3 = root+".out"
        name1 = os.path.join(head,tail2)
        name2 = os.path.join(head,tail3)
        name3 = os.path.join(head,tail)
        # We want the dynG entry last rounding causes problems otherwise
        name4 = os.path.join(head,tail1)
        names = []
        for n in [ name1, name2, name3, name4 ]:
            if os.path.isfile(n) and n not in names:
                names.append(n)
        reader = QEOutputReader( names )
    elif program == "phonopy":
        # The order is important
        pname1 = os.path.join(head,"qpoints.yaml")
        pname2 = os.path.join(head,"phonopy.yaml")
        # Which QM program was used by PHONOPY?
        if qmprogram == "castep":
            print("Error in qmreader",qmprogram)
            sys.exit()
            qmreader = CastepOutputReader(names)
        elif qmprogram == "vasp":
            vname1 = name
            vname2 = os.path.join(head,"KPOINTS")
            pnames = [ pname1, pname2 ]
            vnames = [ vname1, vname2 ]
            pnames.extend(vnames)
            names = pnames
            qmreader = VaspOutputReader(vnames)
        elif qmprogram == "gulp":
            print("Error in qmreader",qmprogram)
            sys.exit()
            qmreader = GulpOutputReader(names)
        elif qmprogram == "crystal":
            names = [ name ]
            qmreader = CrystalOutputReader(names)
            print("I got here 4")
        elif qmprogram == "abinit":
            print("Error in qmreader",qmprogram)
            sys.exit()
            qmreader = AbinitOutputReader(names)
        elif qmprogram == "quantum espresso":
            print("Error in qmreader",qmprogram)
            sys.exit()
            qmreader = QEOutputReader(names)
        # The QM reader is used to get info about the QM calculation
        reader = PhonopyOutputReader(pnames,qmreader)
    else:
        print("Program name not recognized",program,file=sys.stderr)
        sys.exit()
    return reader


def pdgui_get_reader(program,names,qmprogram):
    """Get the appropriate reader based on the provided program and file names.

    This function determines the correct output reader to use based on the program specified and the presence of certain files
    in the provided directory or file names list. It supports various quantum mechanics and molecular dynamics
    programs as well as experiment data.

    Parameters
    ----------
    program : str
        The name of the program or an empty string to attempt automatic detection. Supported programs include 'vasp',
        'castep', 'gulp', 'crystal', 'abinit', 'qe', 'phonopy', and 'experiment'.
    names : list of str
        A list containing either directory names or file names to be read or analysed.
    qmprogram : str
        Specifies the quantum mechanics program used in conjunction with 'phonopy'. It is only relevant if 'program'
        is set to 'phonopy'. Supported programs include 'castep', 'vasp', 'gulp', 'crystal', 'abinit', and 'qe'.

    Returns
    -------
    object
        Returns an instance of the correct output reader based on the specified program and available files. If the
        necessary files are not found or if the program is not recognized, it returns None.

    Notes
    -----
    - The readers returned are capable of parsing output files from their respective programs and extracting relevant data.
    - For 'phonopy', additional handling is included to combine phonopy output with quantum mechanics results from another
      specified program.
    - This function also handles directory or file existence checks and will print errors if required files are missing.

    Examples
    --------
    To get a VASP output reader for a directory containing VASP output files:

    >>> reader = pdgui_get_reader('', ['path/to/vasp/output'], '')
    >>> print(reader)

    For phonopy with VASP as the quantum mechanics program:

    >>> qm_reader = pdgui_get_reader('phonopy', ['path/to/phonopy/files'], 'vasp')
    >>> print(qm_reader)

    """    
    reader = None
    program = program.lower()
    qmprogram = qmprogram.lower()
    if program == "":
        #  This is the old behaviour.  It copes with VASP, CASTEP and Crystal
        #  If names[0] is a directory then we will use a vaspoutputreader
        #  Otherwise it is a seedname for castep, or a gulp output file, or a crystal output file
        if os.path.isdir(names[0]):
            print(f"Analysing VASP directory: {names[0]} ")
            outcarfile = os.path.join(names[0], "OUTCAR")
            kpointsfile = os.path.join(names[0], "KPOINTS")
            if not os.path.isfile(outcarfile):
                print("Error: NO OUTCAR FILE IN DIRECTORY")
                reader = None
                return None
            reader = VaspOutputReader( [outcarfile, kpointsfile] )
        elif names[0].find("OUTCAR") >= 0 and os.path.isfile("OUTCAR"):
            reader = VaspOutputReader(names)
        elif names[0].find(".gout") >= 0 and os.path.isfile(names[0]):
            reader = GulpOutputReader(names)
        elif names[0].find(".out") >= 0 and os.path.isfile(names[0]):
            reader = CrystalOutputReader(names)
        elif names[0].find(".castep") >= 0 and os.path.isfile(names[0]):
            reader = CastepOutputReader(names)
        elif os.path.isfile(names[0]+".castep") and os.path.isfile(names[0]+".castep"):
            reader = CastepOutputReader([names[0]+".castep"])
        else:
            print("No valid file name has been found on the command line")
            print("Try using the -program option to specify the")
            print("files which will be read")
            reader = None
    else:
        # New Specification of Program used to define the input files
        # Abinit and QE need a couple of files to be specified
        #
        # First Check that the file(s) we requested are there
        #
        checkfiles = []
        if program == "castep":
            if names[0].find(".castep") >= 0:
                seedname, ext = os.path.splitext(names[0])
            else:
                seedname = names[0]
            checkfiles.append(seedname+".castep")
        elif program == "phonopy":
            # We only have a VASP / Phonopy interface
            # Creat a list of phonopy files
            pnames = []
            head,tail = os.path.split(names[0])
            pnames.append(os.path.join(head,"qpoints.yaml"))
            pnames.append(os.path.join(head,"phonopy.yaml"))
            # Creat a list of VASP files NB.  They all have to be in the same directory
            vnames = names
            pnames.extend(vnames)
            checkfiles = pnames
        else:
            checkfiles = names
        for f in checkfiles:
            if not os.path.isfile(f):
                print(f"Output files created by program: {program}")
                print(f"Error: file not available: {f}")
                return None
        # The files requested are available so read them
        if program == "castep":
            reader = CastepOutputReader(names)
        elif program == "vasp":
            reader = VaspOutputReader(names)
        elif program == "gulp":
            reader = GulpOutputReader(names)
        elif program == "crystal":
            reader = CrystalOutputReader(names)
        elif program == "abinit":
            reader = AbinitOutputReader(names)
        elif program == "quantum espresso":
            reader = QEOutputReader(names)
        elif program == "phonopy":
            # Which QM program was used by PHONOPY?
            if qmprogram == "castep":
                qmreader = CastepOutputReader(vnames)
            elif qmprogram == "vasp":
                qmreader = VaspOutputReader(vnames)
            elif qmprogram == "gulp":
                qmreader = GulpOutputReader(vnames)
            elif qmprogram == "crystal":
                qmreader = CrystalOutputReader(vnames)
            elif qmprogram == "abinit":
                qmreader = AbinitOutputReader(vnames)
            elif qmprogram == "quantum espresso":
                qmreader = QEOutputReader(vnames)
            # The QM reader is used to get info about the QM calculation
            reader = PhonopyOutputReader(pnames,qmreader)
        elif program == "experiment":
            reader = ExperimentOutputReader(names)
        # endif
    # end if
    return reader

class Debug:
    """A class aimed at providing a structured way to include debug messages in code.

    Methods
    -------
    print(args, level=0)
        Prints debug messages based on the debug level provided as argument compared to the object's debug level.
    state()
        Returns the current state (enabled or disabled) of debugging.

    Notes
    -----
    The `print` method provides a flexible way to include additional information along with the base debug message, allowing for a detailed and adjustable debugging output.

    """

    def __init__(self,debug,text,level=0):
        """Initialize an instance with debug status, text, and optional level.

        Parameters
        ----------
        debug : bool
            A boolean indicating whether debugging is enabled.
        text : str
            The text associated with the instance.
        level : int, optional
            The level of the instance, by default 0.

        """        
        self.debug = debug
        self.text  = text
        self.level = level
        return

    def print(self,*args,level=0):
        """Print message if debugging level allows.

        Parameters
        ----------
        args
            Variable length argument list for the message to be printed.
        level : int, optional
            The level of the message that determines if it gets printed or not,
            based on the instance's `level`. Default value is 0.

        Notes
        -----
        This method will only print the message if the instance's `debug` flag is True
        and the provided `level` is less than or equal to the instance's `level`.

        """        
        if self.debug and level <= self.level:
            print(self.text,*args)
        return

    def state(self,):
        """Get the debug state.

        Parameters
        ----------
        None

        Returns
        -------
        The current debug state.

        """        
        return self.debug

