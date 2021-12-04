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
import os

from PDielec.VaspOutputReader import VaspOutputReader
from PDielec.PhonopyOutputReader import PhonopyOutputReader
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.GulpOutputReader import GulpOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.ExperimentOutputReader import ExperimentOutputReader

def find_program_from_name( filename ):
    # Determine the program to use from the file name being used
    head,tail = os.path.split(filename)
    root,ext = os.path.splitext(tail)
    if head == '':
        head = './'
    else:
        head = head+'/'
    if tail == 'OUTCAR':
        if os.path.isfile(head+'phonopy.yaml'):
            return 'phonopy'
        else:
            return 'vasp'
    if ext == '.gout':
        return 'gulp'
    if ext == '.born':
        return 'phonopy'
    if ext == '.castep':
        return 'castep'
    if ext ==  '.out':
        if os.path.isfile(head+root+'.files'):
            return 'abinit'
        elif os.path.isfile(head+root+'.dynG'):
            return 'qe'
        else:
            return 'crystal'
    if ext ==  '.dynG':
        return 'qe'
    if ext ==  '.exp':
        return 'experiment'
    if ext ==  '.py':
        return 'pdgui'
    return ''

def get_reader( name, program, qmprogram):
    fulldirname = name
    head,tail = os.path.split(fulldirname)
    root,ext = os.path.splitext(tail)
    if program == "castep":
        names = [ name ]
        reader = CastepOutputReader( names )
    elif program == "vasp":
        name1 = name
        name2 = os.path.join(head,'KPOINTS')
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
    elif program == "qe":
        tail1 = root+'.dynG'
        tail2 = root+'.log'
        tail3 = root+'.out'
        name1 = os.path.join(head,tail2)
        name2 = os.path.join(head,tail3)
        name3 = os.path.join(head,tail)
        # We want the dynG entry last rounding causes problems otherwise
        name4 = os.path.join(head,tail1)
        names = []
        for n in [ name1, name2, name3, name4 ]:
            if os.path.isfile(n):
                if not n in names:
                    names.append(n)
        reader = QEOutputReader( names )
    elif program == "phonopy":
        # The order is important
        pname1 = os.path.join(head,'qpoints.yaml')
        pname2 = os.path.join(head,'phonopy.yaml')
        # Only works for VASP at the moment
        vname1 = name
        vname2 = os.path.join(head,'KPOINTS')
        pnames = [ pname1, pname2 ]
        vnames = [ vname1, vname2 ]
        pnames.extend(vnames)
        names = pnames
        # Which QM program was used by PHONOPY?
        if qmprogram == "castep":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = CastepOutputReader(names)
        elif qmprogram == "vasp":
            qmreader = VaspOutputReader(vnames)
        elif qmprogram == "gulp":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = GulpOutputReader(names)
        elif qmprogram == "crystal":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = CrystalOutputReader(names)
        elif qmprogram == "abinit":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = AbinitOutputReader(names)
        elif qmprogram == "qe":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = QEOutputReader(names)
        # The QM reader is used to get info about the QM calculation
        reader = PhonopyOutputReader(pnames,qmreader)
    else:
        print('Program name not recognized',program,file=sys.stderr)
        exit()
    return reader


def pdgui_get_reader(program,names,qmprogram):
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
    def __init__(self,debug,text,level=0):
        self.debug = debug
        self.text  = text
        self.level = level
        return

    def print(self,*args,level=0):
        if self.debug:
            if level <= self.level:
                print(self.text,*args)
        return

