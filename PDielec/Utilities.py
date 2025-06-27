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
    str
        The program name that was used to calculate the frequencies.

    Notes
    -----
    This function examines the file extension and, in some cases, the presence of specific files in the same directory, to determine the associated simulation program. It recognizes files from several popular materials simulation programs, such as 'phonopy', 'gulp', 'vasp', and others.

    Examples
    --------
    ::

        program = find_program_from_name('./data/structure.castep')
        print(program)
        # Output: "castep"
    
        program = find_program_from_name('path/to/simulation/phonopy.yaml')
        print(program)
        # Output: "phonopy"

    """    
    head,tail = os.path.split(filename)
    root,ext = os.path.splitext(tail)
    head_root = os.path.join(head,root)
    if ext == ".dynG":
        return "quantum espresso"
    if tail == "OUTCAR":
        return "vasp"
    if ext ==  ".abo":
            return "abinit"
    if ext ==  ".exp":
        return "experiment"
    if ext ==  ".py":
        return "pdgui"
    if ext ==  ".yaml":
        return "phonopy"
    if ext == ".gout":
        return "gulp"
    if ext == ".castep":
        return "castep"
    if ext ==  ".out":
        if os.path.isfile(head_root+".files"):
            return "abinit"
        elif os.path.isfile(head_root+".dynG"):
            return "quantum espresso"
        else:
            return "crystal"
    if ext ==  ".log":
        if os.path.isfile(head_root+".files"):
            return "abinit",""
        elif os.path.isfile(head_root+".dynG"):
            return "quantum espresso",""
        else:
            return "crystal"
    if os.path.isfile(os.path.join(head,"vasprun.xml")):
        return "vasp"
    if os.path.isfile(os.path.join(head,"pwscf.xml")):
        return "quantum espresso"
    return ""

def get_reader( name, program, debug=False):
    """Get the appropriate output reader based on the simulation program and, if specified, the quantum mechanical program.

    This function is designed to create an output reader object for various simulation programs (like CASTEP, VASP, etc.) and, for phonopy simulations, it can additionally create a quantum mechanical output reader based on the specified quantum mechanical program.

    Parameters
    ----------
    name : str
        The primary filepath or name associated with the output file(s).
    program : str
        The name of the simulation program. Supported values are 'castep', 'vasp', 'gulp', 'crystal', 'abinit', 'qe', and 'phonopy'.
    debug : boolean
        Optional.  If true print debug information.  Default is false

    Returns
    -------
    object
        An instance of the appropriate output reader class based on the input `program`

    Examples
    --------
    >>> reader = get_reader("output.log", "castep")
    >>> reader = get_reader("output", "phonopy")

    """    
    program = program.lower()
    fulldirname = name
    head,tail = os.path.split(fulldirname)
    root,ext = os.path.splitext(tail)
    if debug:
        print("get_reader:  program = ",program)
        print("get_reader:  fulldirname = ",fulldirname)
        print("get_reader:  head = ",head)
        print("get_reader:  tail = ",tail)
        print("get_reader:  root = ",root)
        print("get_reader:  ext = ",ext)
    if program == "castep":
        names = [ name ]
        reader = CastepOutputReader( names )
    elif program == "vasp":
        names = list()
        name1 = name
        if "OUTCAR" in name:
            name2 = os.path.join(head,"KPOINTS")
            names = [ name1, name2 ]
        else:
            names = [ name1 ]
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
        pwscf_name = os.path.join(head,"pwscf.xml")
        tensors_name = os.path.join(head,"tensors.xml")
        if ext in [".xml",".dynG"] and os.path.isfile(pwscf_name) and os.path.isfile(tensors_name):
            # Only handle xml but add the dynG file
            tail3 = root+".dynG"               # The order is important
            name1 = pwscf_name
            name2 = tensors_name
            name3 = os.path.join(head,tail3)
            names = []
            for n in [ name1, name2, name3 ]:
                if os.path.isfile(n):
                    names.append(n)
        else:
            # try and read the log file or out file based on the dynG file name
            tail1 = root+".log"
            tail2 = root+".out"
            tail3 = root+".dynG"               # The order is important
            name1 = os.path.join(head,tail1)
            name2 = os.path.join(head,tail2)
            name3 = os.path.join(head,tail3)
            names = []
            for n in [ name1, name2, name3 ]:
                if os.path.isfile(n):
                    names.append(n)
        if debug:
            print("get_reader:  names = ",names,flush=True)
        reader = QEOutputReader( names )
    elif program == "phonopy":
        # The order is important
        pnames = []
        pnames.append( os.path.join(head,"qpoints.yaml") )
        pnames.append( os.path.join(head,"phonopy.yaml") )
        pnames.append( os.path.join(head,"BORN_PDIELEC") )
        reader = PhonopyOutputReader(pnames)
    elif program == "experiment":
        names = [ name ]
        reader = ExperimentOutputReader(names)
    else:
        print("Program name not recognized",program,file=sys.stderr)
        sys.exit()
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
            print(self.text,*args,flush=True)
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

