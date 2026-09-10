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

import glob
import logging
import os
import sys

from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.AimsOutputReader import AimsOutputReader
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.ExperimentOutputReader import ExperimentOutputReader
from PDielec.FiniteFieldOutputReader import FiniteFieldOutputReader
from PDielec.GulpOutputReader import GulpOutputReader
from PDielec.PhonopyOutputReader import PhonopyOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.VaspOutputReader import VaspOutputReader

logger = logging.getLogger(__name__)



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
    logger.debug("")
    logger.debug(name)
    if len(matrix.shape) == 1:
        columns = matrix.shape[0]
        string = ""
        for i in range(columns):
            string = string + f"{matrix[i]:+.5f}" + "   "
        logger.debug(string)
    else:
        rows = matrix.shape[0]
        columns = matrix.shape[1]
        for j in range(rows):
            string = ""
            for i in range(columns):
                string = string + f"{matrix[j,i]:+.5f}" + "   "
            logger.debug(string)
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
    This function examines the file extension and, in some cases, the presence of specific files in the same directory,
    to determine the associated simulation program. It recognizes files from several popular materials simulation
    programs, such as 'phonopy', 'gulp', 'vasp', and others.

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
    if ext == ".json":
        try:
            import json
            with open(filename) as fd:
                schema = json.load(fd).get("schema", "")
            if schema.startswith("pdielect-") and ("raman" in schema or "r_epsilon" in schema):
                return "finite_field"
        except Exception:
            pass
    if ext ==  ".out":
        if os.path.isfile(head_root+".files"):
            return "abinit"
        elif os.path.isfile(head_root+".dynG") or glob.glob(os.path.join(head if head else ".", "*.dynG")):
            return "quantum espresso"
        else:
            return "crystal"
    if ext ==  ".log":
        if os.path.isfile(head_root+".files"):
            return "abinit"
        elif os.path.isfile(head_root+".dynG") or glob.glob(os.path.join(head if head else ".", "*.dynG")):
            return "quantum espresso"
        else:
            return "crystal"
    if ext ==  ".dat":
        return "aims"
    if tail ==  "aims.out":
        return "aims"
    if ext == ".xml":
        if tail == "vasprun.xml":
            return "vasp"
        # Peek inside to check for a QE output or tensors file
        try:
            import xml.etree.ElementTree as ET
            xml_tree = ET.parse(filename)
            xml_root = xml_tree.getroot()
            xml_tag = xml_root.tag
            # pwscf / SCF output: root tag contains 'espresso'
            if "espresso" in xml_tag:
                return "quantum espresso"
            # tensors.xml: root is 'Root' with an EF_TENSORS child
            if xml_tag == "Root" and xml_root.find("EF_TENSORS") is not None:
                return "quantum espresso"
        except Exception:
            pass
    if os.path.isfile(os.path.join(head,"vasprun.xml")):
        return "vasp"
    if os.path.isfile(os.path.join(head,"pwscf.xml")):
        return "quantum espresso"
    if os.path.isfile(os.path.join(head,"geometry.in")):
        return "aims"
    if os.path.isfile(os.path.join(head,"control.in")):
        return "aims"
    if os.path.isfile(os.path.join(head,"aims.out")):
        return "aims"
    return ""

def get_reader( name, program):
    """Get the appropriate output reader based on the simulation program and, if specified, the quantum mechanical program.

    This function is designed to create an output reader object for various simulation programs (like CASTEP, VASP,
    etc.) and, for phonopy simulations, it can additionally create a quantum mechanical output reader based on the
    specified quantum mechanical program.

    Parameters
    ----------
    name : str
        The primary filepath or name associated with the output file(s).
    program : str
        The name of the simulation program. Supported values are 'castep', 'vasp', 'gulp', 'crystal', 'abinit',
        'qe', and 'phonopy'.

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
    logger.debug(f"get_reader:  program = {program}")
    logger.debug(f"get_reader:  fulldirname = {fulldirname}")
    logger.debug(f"get_reader:  head = {head}")
    logger.debug(f"get_reader:  tail = {tail}")
    logger.debug(f"get_reader:  root = {root}")
    logger.debug(f"get_reader:  ext = {ext}")
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
        tensors_name = os.path.join(head, "tensors.xml")
        if ext == ".xml":
            # User passed an xml file directly (e.g. ZnO.xml or pwscf.xml)
            # Use it as the primary source, then add tensors.xml and matching dynG
            names = [name]
            if os.path.isfile(tensors_name) and os.path.abspath(name) != os.path.abspath(tensors_name):
                names.append(tensors_name)
            # Look for a dynG whose root matches the xml filename (e.g. ZnO.dynG for ZnO.xml)
            dynG_name = os.path.join(head, root + ".dynG")
            if os.path.isfile(dynG_name):
                names.append(dynG_name)
            else:
                # The dynG may have a different name (set via fildyn= in the ph.x input).
                # If there is exactly one *.dynG file in the same directory, use it.
                search_dir = head if head else "."
                dynG_candidates = glob.glob(os.path.join(search_dir, "*.dynG"))
                if len(dynG_candidates) == 1:
                    names.append(dynG_candidates[0])
                    logger.debug(f"get_reader: auto-selected dynG {dynG_candidates[0]}")
        else:
            # Try to read the log file, out file, and/or dynG file derived from the given name
            tail1 = root+".log"
            tail2 = root+".out"
            tail3 = root+".dynG"               # The order is important
            name1 = os.path.join(head,tail1)
            name2 = os.path.join(head,tail2)
            name3 = os.path.join(head,tail3)
            names = []
            # Add a matching SCF xml file first if present (e.g. ZnO.xml for ZnO.dynG)
            xml_name = os.path.join(head, root + ".xml")
            if os.path.isfile(xml_name) and os.path.abspath(xml_name) != os.path.abspath(tensors_name):
                names.append(xml_name)
            for n in [ name1, name2, name3 ]:
                if os.path.isfile(n):
                    names.append(n)
            # If no exact-match dynG was found, look for any single *.dynG in the directory
            # (handles cases where fildyn= in ph.x gave the dynG a different prefix)
            if name3 not in names:
                search_dir = head if head else "."
                dynG_candidates = glob.glob(os.path.join(search_dir, "*.dynG"))
                if len(dynG_candidates) == 1:
                    names.append(dynG_candidates[0])
                    logger.debug(f"get_reader: auto-selected dynG {dynG_candidates[0]}")
            # Also read tensors.xml if present (higher-precision Raman tensors)
            if os.path.isfile(tensors_name):
                names.append(tensors_name)
        logger.debug(f"get_reader:  names = {names}")
        reader = QEOutputReader( names )
    elif program == "phonopy":
        # The order is important
        pnames = []
        pnames.append( os.path.join(head,"qpoints.yaml") )
        pnames.append( os.path.join(head,"phonopy.yaml") )
        pnames.append( os.path.join(head,"BORN_PDIELEC") )
        reader = PhonopyOutputReader(pnames)
    elif program == "aims":
        if root.endswith(".dat"):
            identifier = root.split(".")[1]
        else:
            all_files = glob.glob(os.path.join(head,"hessian.*.dat"))
            if len(all_files) < 1:
                identifier="unkown_and_ignore"
            else:
                identifier = all_files[-1]
                _head,_tail = os.path.split(identifier)
                identifier = _tail.split(".")[1]
        # The order is important
        pnames = []
        aimsfile = os.path.join(head,"aims.out")
        if not os.path.exists(aimsfile):
            aimsfile = aimsfile+"unkown_and_ignore"
        pnames.append( aimsfile )
        pnames.append( os.path.join(head,"masses."+identifier+".dat") )
        pnames.append( os.path.join(head,"geometry.in") )
        pnames.append( os.path.join(head,"hessian."+identifier+".dat") )
        bornfile = os.path.join(head,"BORN_PDIELEC")
        if not os.path.exists(bornfile):
            bornfile = bornfile+"unkown_and_ignore"
        pnames.append(bornfile)
        reader = AimsOutputReader(pnames)
    elif program == "experiment":
        names = [ name ]
        reader = ExperimentOutputReader(names)
    elif program == "finite_field":
        names = [ name ]
        reader = FiniteFieldOutputReader(names)
    else:
        logger.debug(f"Program name not recognized {program}")
        sys.exit()
    return reader
