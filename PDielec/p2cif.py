#!/usr/bin/env python
"""Read the contents of a directory containing DFT output and create a cif file of the structure"""
from __future__ import print_function
import os, sys
import psutil
from PDielec.VaspOutputReader import VaspOutputReader
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.GulpOutputReader import GulpOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.PhonopyOutputReader import PhonopyOutputReader
import PDielec.Utilities as Utilities
from multiprocessing import Pool

def set_affinity_on_worker():
    """When a new worker process is created, the affinity is set to all CPUs"""
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK for the time being this is simply commented out, but might be useful at some point
    #os.system("taskset -p 0xff %d > /dev/null" % os.getpid())

def read_a_file( calling_parameters):
    name, program, qmprogram, debug = calling_parameters
    reader = Utilities.get_reader(name,program,qmprogram)
    reader.debug = debug
    reader.read_output()
    return name,reader.unit_cells[-1]

def main():
    # Start processing the directories
    if len(sys.argv) <= 1 :
        print('p2cif -program program filenames .....', file=sys.stderr)
        print('  \"program\" must be one of \"abinit\", \"castep\", \"crystal\", \"gulp\"       ', file=sys.stderr)
        print('           \"phonopy\", \"qe\", \"vasp\", \"experiment\", \"auto\"               ', file=sys.stderr)
        print('           The default is auto, so the program tries to guess the package from   ', file=sys.stderr)
        print('           the contents of the directory.  However this is not fool-proof!       ', file=sys.stderr)
        print('           If phonopy is used it must be followed by the QM package              ', file=sys.stderr)
        print('           in auto mode if the file was created by a phonopy VASP is assumed     ', file=sys.stderr)
        print('  -debug   to switch on more debug information                                   ', file=sys.stderr)
        exit()

    files = []
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    program = 'auto'
    qmprogram = 'vasp'
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

    if not program in ['abinit','castep','crystal','gulp','qe','vasp','phonopy','experiment','auto']:
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
    number_of_processors = psutil.cpu_count(logical=False)
    p = Pool(number_of_processors,initializer=set_affinity_on_worker)
    # Create a tuple list of calling parameters
    calling_parameters = []
    files.sort()
    for name in files:
        prog = program
        if program == 'auto':
            prog = Utilities.find_program_from_name(name)
        calling_parameters.append( (name, prog, qmprogram, debug) )
    # Calculate the results in parallel
    results_map_object = p.map_async(read_a_file,calling_parameters)
    results_map_object.wait()
    results = results_map_object.get()
    for name,cell in results:
        cell.write_cif(filename=name,file_=sys.stdout)
    #
    exit()
# end of def main

if __name__ == "__main__":
    main(sys)
