#!/usr/bin/env python
"""Read the contents of a directory containing DFT output and create a cif file of the structure"""
from __future__ import print_function
import string
import re
import numpy as np
import os, sys
from PDielec.Constants import amu, PI, avogadro_si, wavenumber, angstrom, isotope_masses, average_masses
from PDielec.VaspOutputReader import VaspOutputReader
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.GulpOutputReader import GulpOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.PhonopyOutputReader import PhonopyOutputReader
from multiprocessing import Pool, cpu_count
import PDielec.Calculator as Calculator

def set_affinity_on_worker():
    """When a new worker process is created, the affinity is set to all CPUs"""
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK for the time being this is simply commented out, but might be useful at some point
    #os.system("taskset -p 0xff %d > /dev/null" % os.getpid())

def read_a_file( calling_parameters):
    name, program, qmprogram, debug = calling_parameters
    fulldirname = name
    head,tail = os.path.split(fulldirname)
    root,ext = os.path.splitext(tail)
    if program == "castep":
        names = [ name ]
        reader = CastepOutputReader( names )
    elif program == "vasp":
        name1 = os.path.join(head,'OUTCAR')
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
        vname1 = os.path.join(head,'OUTCAR')
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

    print('  Analysing ',names, file=sys.stderr)
    reader.debug = debug
    reader.read_output()
    return name,reader.unit_cells[-1]

def main():
    # Start processing the directories
    if len(sys.argv) <= 1 :
        print('p2cif -program program filenames .....', file=sys.stderr)
        print('  \"program\" must be one of \"abinit\", \"castep\", \"crystal\", \"gulp\"       ', file=sys.stderr)
        print('           \"phonopy\", \"qe\", \"vasp\"                                         ', file=sys.stderr)
        print('           If phonopy is used it must be followed by the QM package              ', file=sys.stderr)
        print('  -debug   to switch on more debug information                                   ', file=sys.stderr)
        exit()

    files = []
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    program = ''
    qmprogram = ''
    debug = False
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token == "-debug":
            debug = True
        elif token == "-program":
            itoken += 1
            program = tokens[itoken]
            if program == 'phonopy':
                itoken += 1
                qmprogram = tokens[itoken]
        else:
            files.append(token)

    if len(program) < 1:
        print('Please use -program to define the package used to generate the output files',file=sys.stderr)
        exit()

    if not program in ['abinit','castep','crystal','gulp','qe','vasp','phonopy']:
        print('Program is not recognised: ',program,file=sys.stderr)
        exit()

    if program == 'phonopy':
        if not qmprogram in ['abinit','castep','crystal','gulp','qe','vasp']:
            print('Phonopy QM program is not recognised: ',qmprogram,file=sys.stderr)
            exit()
        print('  QM program used by Phonopy is: ',qmprogram,file=sys.stderr)

    print('  Program is ',program,file=sys.stderr)

    for f in files:
        if not os.path.isfile(f):
            print('Error file requested for analysis does not exist',f,file=sys.stderr)
            exit()

    #
    # Create a pool of processors to handle reading the files
    #
    p = Pool(initializer=set_affinity_on_worker)
    # Create a tuple list of calling parameters
    calling_parameters = []
    files.sort()
    for name in files:
        calling_parameters.append( (name, program, qmprogram, debug) )
    # Calculate the results in parallel
    results_map_object = p.map_async(read_a_file,calling_parameters)
    results_map_object.wait()
    results = results_map_object.get()
    # Convert the results into a dictionary
    results_dictionary = {}
    for name,cell in results:
        cell.write_cif(filename=name,file_=sys.stdout)
    #
    exit()
# end of def main

if __name__ == "__main__":
    main(sys)
