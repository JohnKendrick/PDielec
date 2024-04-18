#!/usr/bin/env python
"""
Read the contents of a directory containing DFT output and create a csv style file of information

The MIT License (MIT)

Copyright (c) 2024 John Kendrick

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import sys
import dill as pickle
import PDielec.__init__
version = PDielec.__init__.__version__

def print_help():
    """
    Print help information to standard error and exit the program.

    This function prints out instructions for using a command-line tool that reads in a pickled file. The file should have been created with a specific flag ('-pickle') using the 'preader' command. It also outputs the version of the program.

    No parameters are accepted or required by this function.

    The function uses `sys.stderr` to print messages, ensuring that they are output to the standard error stream rather than the standard output. This is useful for separating actual program output from error or help messages.

    The function explicitly calls `exit()`, which terminates the program execution. This is used here to halt further execution after displaying the help information.
    """    
    print('pickled_reader filenames', file=sys.stderr)
    print('  Read in a pickled (actually using dill to pickle the object) reader            ', file=sys.stderr)
    print('  The pickled file should have been created using preader -pickle                ', file=sys.stderr)
    print('  Version ',version,file=sys.stderr)
    exit()

def main():
    """
    Main program entry point.

     his function reads a serialized object from a pickle file specified as the first command line argument. It deserializes objects from the file until it reaches the end of the file and then iterates through these objects, printing their details and the details of their last unit cell.

    The expected structure for the serialized objects (readers) should have attributes `type`, `names`, and `unit_cells`, and must have a `print()` method. The `unit_cells` attribute is expected to be a list, where each element has a `print()` method as well.

    Parameters
    ----------
    This function does not accept any parameters directly but expects at least one argument from the command line:

    1. The path to the pickle file as the first command line argument.

    Raises
    ------
    FileNotFoundError
        If the specified pickle file does not exist.

    EOFError
        If there is an issue with reading the file to its end. This exception is caught and handled within the function.

    Any exception related to deserialization issues with `pickle.load()` is also caught and should be handled as per the application's needs.

    Note
    ----
    The function requires the 'sys' and 'pickle' modules to be imported.
    """    
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
        reader.print()
        print('LAST CELL')
        lastcell = reader.unit_cells[-1]
        lastcell.print()
    # End of for loop over readers
# end of def main

if __name__ == "__main__":
    main()
