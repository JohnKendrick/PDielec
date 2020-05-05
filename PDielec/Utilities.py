#!/usr/bin/python
#
# Copyright 2015 John Kendrick
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
# You should have received a copy of the MIT License
# along with this program, if not see https://opensource.org/licenses/MIT
#
"""Generic reader of output files.  Actual reader should inherit from this class"""
from __future__ import print_function
import math
import os
import sys
import numpy as np

from PDielec.VaspOutputReader import VaspOutputReader
from PDielec.PhonopyOutputReader import PhonopyOutputReader
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.GulpOutputReader import GulpOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.ExperimentOutputReader import ExperimentOutputReader



def get_reader(program,names,qmprogram):
    reader = None
    program = program.lower()
    qmprogram = qmprogram.lower()
    #print("get_reader",names,program)
    if program == "":
        #  This is the old behaviour.  It copes with VASP, CASTEP and Crystal
        #  If names[0] is a directory then we will use a vaspoutputreader
        #  Otherwise it is a seedname for castep, or a gulp output file, or a crystal output file
        if os.path.isdir(names[0]):
            print('Analysing VASP directory: {} '.format(names[0]))
            outcarfile = os.path.join(names[0], "OUTCAR")
            kpointsfile = os.path.join(names[0], "KPOINTS")
            if not os.path.isfile(outcarfile):
                print("Error: NO OUTCAR FILE IN DIRECTORY")
                reader = None
                return
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
            print('No valid file name has been found on the command line')
            print('Try using the -program option to specify the')
            print('files which will be read')
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
            pnames.append(os.path.join(head,'qpoints.yaml'))
            pnames.append(os.path.join(head,'phonopy.yaml'))
            # Creat a list of VASP files NB.  They all have to be in the same directory
            vnames = names
            pnames.extend(vnames)
            checkfiles = pnames
        else:
            checkfiles = names
        #jk print("")
        #jk print("Program used to perform the phonon calculation was: {}".format(program))
        for f in checkfiles:
            #jk print("The file containing the output is: {}".format(f))
            if not os.path.isfile(f):
                print("Output files created by program: {}".format(program))
                print("Error: file not available: {}".format(f))
                reader = None
                return reader
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
        elif program == "qe":
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
            elif qmprogram == "qe":
                qmreader = QEOutputReader(vnames)
            # The QM reader is used to get info about the QM calculation
            reader = PhonopyOutputReader(pnames,qmreader)
        elif program == "experiment":
            reader = ExperimentOutputReader(names)
        # endif
    # end if
    return reader

class Debug():
    def __init__(self,debug,text):
        self.debug = debug
        self.text  = text
        return
    def print(self,*args):

        if self.debug:
            print(self.text,*args)
        return

