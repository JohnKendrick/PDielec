#!/usr/bin/env python
"""
Read the contents of a directory containing DFT output and create a cif file of the structure

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
    """
    When a new worker process is created, the affinity is set to all CPUs.

    Notes
    -----
    None.
    """
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK for the time being this is simply commented out, but might be useful at some point
    #os.system("taskset -p 0xff %d > /dev/null" % os.getpid())

def read_a_file( calling_parameters):
    """
    Read data from a file and process it using specified reader utilities.

    Parameters
    ----------
    calling_parameters : tuple
        A tuple containing the parameters for the file reading operation. Expected to contain
        - name (str): The name of the file to be read.
        - program (str): The type of program or file format to use for reading.
        - qmprogram (str): Specific quantum mechanics program details, if applicable.
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
    name, program, qmprogram, debug = calling_parameters
    reader = Utilities.get_reader(name,program,qmprogram)
    reader.debug = debug
    reader.read_output()
    return name,reader.get_unit_cell()

def main():
    # Start processing the directories
    """
    Main entry point for processing and converting files.

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
    p = Pool(number_of_processors)
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
    p.close()
    p.join()
    return
# end of def main

if __name__ == "__main__":
    main()
