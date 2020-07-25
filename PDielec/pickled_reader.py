#!/usr/bin/env python
"""Read the contents of a directory containing DFT output and create a csv style file of information"""
from __future__ import print_function
import string
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
import dill
import PDielec.__init__
version = PDielec.__init__.__version__

def print_help():
    print('pickled_reader filenames', file=sys.stderr)
    print('  Read in a pickled (actually using dill to pickle the object) reader            ', file=sys.stderr)
    print('  The pickled file should have been created using preader -pickle                ', file=sys.stderr)
    print('  Version ',version,file=sys.stderr)
    exit()

def main():
    # Start processing the directories
    if len(sys.argv) <= 1 :
        print_help()
    # Read in the pickled reader objects from the dump file
    readers = []
    picklefile = sys.argv[1]
    with open(picklefile,'rb') as f:
        readers.append(dill.load(f))
    for reader in readers:
        reader.print_info()
        print('LAST CELL')
        lastcell = reader.unit_cells[-1]
        lastcell.print_info()
    # End of for loop over readers
# end of def main

if __name__ == "__main__":
    main()
