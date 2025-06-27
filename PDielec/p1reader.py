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
"""Read the contents of a single file containing DFT output and create a csv style file of information.

Parse command-line arguments to set program configuration and process a specified file.

This function does not have a formal return but exits with different messages and statuses based on the command-line inputs. It processes various command-line arguments to configure the behavior of a software utility, specifically targeting file analysis with support for multiple programs. It also supports debug mode and handles command-line requests for help and version information.

Parameters
----------
This function does not accept parameters directly through its definition. Instead, it uses command-line arguments parsed from `sys.argv`:

- `-debug` : Enable debug mode.
- `-help` : Print help message and exit.
- `-version` : Print software version and exit.
- `-program` : Specify the program to use for analysis. When set to 'phonopy', an additional quantum mechanics program argument is required.
- Any other argument is considered as the filename for analysis.

"""

import os
import sys

import PDielec.__init__
from PDielec import Utilities

version = PDielec.__init__.__version__

def print_help():
    """Print the help message for the p1reader program to stderr.

    This function displays the usage of the `p1reader` command, including accepted programs, flags, and additional
    optional arguments.

    Notes
    -----
    This function directly prints to `sys.stderr` and exits the program, it does not return any value.
    It also prints the current version of the PDielec library being used.

    Examples
    --------
    To display the help message, simply call the function without any arguments:

    >>> print_help()

    """    
    print("p1reader -program program [-version] filename", file=sys.stderr)
    print('  "program" must be one of "abinit", "castep", "crystal", "gulp"       ', file=sys.stderr)
    print('           "phonopy", "qe", "vasp", "experiment", "auto"               ', file=sys.stderr)
    print("           The default is auto, so the program tries to guess the package from   ", file=sys.stderr)
    print("           the contents of the directory.  However this is not fool-proof!       ", file=sys.stderr)
    print("  -debug   to switch on more debug information                                   ", file=sys.stderr)
    print("  -version print the version of PDielec library being used                       ", file=sys.stderr)
    print("  Version ",version,file=sys.stderr)
    sys.exit()


def main():
    # Start processing the directories
    """Parse command-line arguments to set program configuration and process a specified file.

    This function does not have a formal return but exits with different messages and statuses based on the command-line inputs. It processes various command-line arguments to configure the behavior of a software utility, specifically targeting file analysis with support for multiple programs. It also supports debug mode and handles command-line requests for help and version information.

    Parameters
    ----------
    This function does not accept parameters directly through its definition. Instead, it uses command-line arguments parsed from `sys.argv`:
    - `-debug` : Enable debug mode.
    - `-help` : Print help message and exit.
    - `-version` : Print software version and exit.
    - `-program` : Specify the program to use for analysis. When set to 'phonopy', an additional quantum mechanics program argument is required.
    - Any other argument is considered as the filename for analysis.

    Raises
    ------
    SystemExit
        This function may call `exit()` which stops script execution and can be considered as raising `SystemExit` in various scenarios:
        - If less than one argument is provided.
        - If the `-version` argument is provided.
        - If no program is specified or an unrecognized program is given.
        - If the specified file does not exist.
        - When critical configuration issues are encountered.

    Notes
    -----
    - It is assumed that `print_help`, `Utilities.find_program_from_name`, and `Utilities.get_reader` are defined elsewhere.
    - The function makes extensive use of `sys.argv` for argument parsing and `sys.stderr` for error and status messaging.
    - The function supports a flexible addition of program types and respective validation.
    - The handling of `-debug`, `-help`, `-version`, and program-specific arguments provides an example of a rudimentary command-line interface.

    """    
    if len(sys.argv) <= 1 :
        print_help()
    filename = ""
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    program = "auto"
    debug = False
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        token = token.replace("--","-")
        if token == "-debug":
            debug = True
        elif token == "-help":
            print_help()
        elif token == "-version":
            print("  Version ",version,file=sys.stderr)
            sys.exit()
        elif token == "-program":
            itoken += 1
            program = tokens[itoken]
        else:
            filename = tokens[itoken]

    if len(program) < 1:
        print("Please give a filename to be read in",file=sys.stderr)
        sys.exit()

    if program == "qe":
        program = "quantum espresso"

    if program not in ["auto","abinit","castep","crystal","gulp","quantum espresso","vasp","phonopy","experiment"]:
        print("Program is not recognised: ",program,file=sys.stderr)
        sys.exit()

    print("  Program is ",program,file=sys.stderr)

    if not os.path.isfile(filename):
        print("Error file requested for analysis does not exist",filename,file=sys.stderr)
        sys.exit()
    #
    # If no program information was given try and work out what package created the outputfile
    #
    if program == "auto":
        program = Utilities.find_program_from_name(filename)
    #
    # Print out what we are doing
    #
    print(f"  Analysing {filename} generated by {program}",file=sys.stderr)
    #
    # Get the reader from the filename and package used to create it
    #
    reader = Utilities.get_reader(filename, program)
    #
    # applying before reading the file debug 
    #
    reader.debug = debug
    #
    # Now read the output file
    #
    reader.read_output()
    #
    # Test to make sure we have a functioning reader
    #
    reader.print()
    return

if __name__ == "__main__":
    main()
