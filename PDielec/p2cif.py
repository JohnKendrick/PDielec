#!/usr/bin/env python
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
"""Read the contents of a directory containing DFT output and create a cif file of the structure.

This script accepts command-line arguments to define the program used for generating output files, enables debug mode, and processes specified files. If no arguments are provided, it prints usage information. It sorts input files, processes them using a pool of workers equal to the number of logical processors available, and writes the results to standard output.

Command-line arguments:
    - `-debug`: Enable debug mode for more verbose output.
    - `-program <program>`: Specify the program used to generate the output files. The program must be one of: 'abinit', 'castep', 'crystal', 'gulp', 'phonopy', 'qe', 'vasp', 'experiment', 'auto'. For 'phonopy', a quantum mechanics program must be specified immediately after.
    - Files to process: Specify one or more filenames after the arguments.

This function uses multiprocessing to parallelize file processing, with the number of processes equal to the physical cores available on the machine.

Examples
--------
Command-line usage examples:

- Process files in automatic mode:
    ```
    python script_name.py file1 file2 ...
    ```

- Specify a program and enable debug mode:
    ```
    python script_name.py -debug -program gulp file1 file2 ...
    ```

"""
import os
import sys
from multiprocessing import Pool

import psutil

from PDielec import Utilities


def set_affinity_on_worker():
    """When a new worker process is created, the affinity is set to all CPUs.

    Notes
    -----
    None.

    """
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK for the time being this is simply commented out, but might be useful at some point
    #os.system("taskset -p 0xff %d > /dev/null" % os.getpid())

def read_a_file( calling_parameters):
    """Read data from a file and process it using specified reader utilities.

    Parameters
    ----------
    calling_parameters : tuple
        A tuple containing the parameters for the file reading operation. Expected to contain
        - name (str): The name of the file to be read.
        - program (str): The type of program or file format to use for reading.
        - debug (bool): A flag indicating whether debugging is enabled.

    Returns
    -------
    tuple
        A tuple containing:
        - name (str): The name of the file read.
        - unit_cell: The last unit cell read from the file.

    Notes
    -----
    This function requires that a `Utilities` class with a `get_reader` method is available 
    in the scope where this function is used. The `get_reader` method should return an object 
    with a `read_output` method and a `unit_cells` attribute that is a list.

    """    
    name, program, debug = calling_parameters
    reader = Utilities.get_reader(name,program)
    reader.debug = debug
    reader.read_output()
    return name,reader.get_unit_cell()

def main():
    # Start processing the directories
    """Process and convert files.

    This script accepts command-line arguments to define the program used for generating output files, enables debug mode, and processes specified files. If no arguments are provided, it prints usage information. It sorts input files, processes them using a pool of workers equal to the number of logical processors available, and writes the results to standard output.

    Notes
    -----
    Command-line arguments:
    - `-debug`: Enable debug mode for more verbose output.
    - `-program <program>`: Specify the program used to generate the output files. The program must be one of: 'abinit', 'castep', 'crystal', 'gulp', 'phonopy', 'qe', 'vasp', 'experiment', 'auto'. For 'phonopy', a quantum mechanics program must be specified immediately after.
    - Files to process: Specify one or more filenames after the arguments.

    This function uses multiprocessing to parallelize file processing, with the number of processes equal to the physical cores available on the machine.

    Examples
    --------
    Command-line usage examples:

    - Process files in automatic mode:
        ```
        python script_name.py file1 file2 ...
        ```

    - Specify a program and enable debug mode:
        ```
        python script_name.py -debug -program gulp file1 file2 ...
        ```

    """    
    if len(sys.argv) <= 1 :
        print("p2cif -program program filenames .....", file=sys.stderr)
        print('  "program" must be one of "abinit", "castep", "crystal", "gulp"       ', file=sys.stderr)
        print('           "phonopy", "qe", "vasp", "experiment", "auto"               ', file=sys.stderr)
        print("           The default is auto, so the program tries to guess the package from   ", file=sys.stderr)
        print("           the contents of the directory.  However this is not fool-proof!       ", file=sys.stderr)
        print("  -debug   to switch on more debug information                                   ", file=sys.stderr)
        sys.exit()

    files = []
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    program = "auto"
    debug = False
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token == "-debug":
            debug = True
        elif token == "-program":
            itoken += 1
            program = tokens[itoken]
        else:
            files.append(token)

    if program == "qe":
        program = "quantum espresso"

    if len(program) < 1:
        print("Please use -program to define the package used to generate the output files",file=sys.stderr)
        sys.exit()

    if program not in ["abinit","castep","crystal","gulp","quantum espresso","vasp","phonopy","experiment","auto"]:
        print("Program is not recognised: ",program,file=sys.stderr)
        sys.exit()

    print("  Program is ",program,file=sys.stderr)

    for f in files:
        if not os.path.isfile(f):
            print("Error file requested for analysis does not exist",f,file=sys.stderr)
            sys.exit()

    #
    # Create a pool of processors to handle reading the files
    #
    number_of_processors = psutil.cpu_count(logical=False)
    p = Pool(number_of_processors)
    # Create a tuple list of calling parameters
    calling_parameters = []
    files.sort()
    for name in files:
        prog = program
        if program == "auto":
            prog = Utilities.find_program_from_name(name)
        calling_parameters.append( (name, prog, debug) )
    # Calculate the results in parallel
    results_map_object = p.map_async(read_a_file,calling_parameters)
    results_map_object.wait()
    results = results_map_object.get()
    for _name,cell in results:
        cell.write_cif()
    p.close()
    p.join()
    return
# end of def main

if __name__ == "__main__":
    main()
