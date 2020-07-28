#!/usr/bin/env python
"""Read the contents of a directory containing DFT output and create a csv style file of information"""
from __future__ import print_function
import sys
import dill as pickle
import PDielec.__init__
version = PDielec.__init__.__version__

def print_help():
    print('pickled_reader filenames', file=sys.stderr)
    print('  Read in a pickled (actually using dill to pickle the object) reader            ', file=sys.stderr)
    print('  The pickled file should have been created using preader -pickle                ', file=sys.stderr)
    print('  Version ',version,file=sys.stderr)
    exit()

def main():
    #
    # Print out the help file if there is nothing else on the command line
    #
    if len(sys.argv) <= 1 :
        print_help()
    #
    # Read in the pickled reader objects from the dump file
    #
    picklefile = sys.argv[1]
    #
    # store each reader in a list
    #
    readers = []
    #
    # Open the pickled file as binary and for reading only
    # keep reading until we reach an end of file
    #
    with open(picklefile,'rb') as f:
        try:
            while True:
                readers.append(pickle.load(f))
        except EOFError:
            pass
    #
    print('Read in {} readers'.format(len(readers)))
    #
    # Loop over the readers and print out some information - assign a variable
    #
    for reader in readers:
        print('NEW READER type={}, file={}'.format(reader.type,reader.names[0]))
        reader.print_info()
        print('LAST CELL')
        lastcell = reader.unit_cells[-1]
        lastcell.print_info()
    # End of for loop over readers
# end of def main

if __name__ == "__main__":
    main()
