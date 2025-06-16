#!/usr/bin/env python
#
# Copyright 2025 John Kendrick & Andrew Burnett
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
"""PDielec interface program to Phonopy.

The program reads the phonopy.yaml file and the BORN file.
The BORN file must be created using the standard Phonopy tools
The program writes out the static permittivity and the Born charges in a format similar to Phonopy's.   But all atoms in the unit cell are in the file.  Not just the symmetry unique ones.

phonopy-pdielec-born BORN-PDIELEC
    This command writes out the permitivity and the Born charges to the file "BORN-PDIELEC"

"""

import os
import sys
import phonopy
from contextlib import nullcontext
from phonopy.file_IO import parse_BORN


def show_usage():
    """Show the usage message."""
    print("")
    print("USAGE:")
    print("phonopy-pdielec-born [--symprec 1.0E-5] [-input BORN] [-nosym] [-output BORN_PDIELEC] [-yaml phonopy.yaml] [outputfilename]")
    print(" -symprec Sets the precision used in determing the symmetry (default is 1.0E-5)")
    print(" -input   Sets the input Born charges file (Defaults is BORN)")
    print(" -output  Sets the output Born charges file")
    print(" -yaml    phonopy.yaml Sets the output Born charges file (Defaults is phonopy.yaml)")
    print(" -nosym   Stops symmetry being used in parsing the BORN file")
    print("")
    print(" If no output file is given the output goes to stdout       ")


def main():
    """Define the main routine for the phonopy-pdielec-born command.

    The routine makes extensive use of Phonopy's API
    It reads the phonopy.yaml file to determine the primitive cell
    then it reads the BORN file, before writing out the permitivity and Born charges.
    The output file is formatted in a similar way to the BORN file, but
    all the atoms/ions are present, not just the symmetry unique ones.

    Options:
    -symprec 1.0E-5       Sets the precision used in determing the symmetry (default is 1.0E-5)
    -input   BORN         Sets the input Born charges file (Defaults is BORN)
    -output  BORN_PDIELEC Sets the output Born charges file 
    -yaml    phonopy.yaml Sets the output Born charges file (Defaults is phonopy.yaml)
    -nosym                Stops symmetry being used in parsing the BORN file
    """
    symprec = 1.e-5
    inputfilename = "BORN"
    outputfilename = None
    symmetry = True
    tokens = sys.argv[1:]
    ntokens = len(tokens)
    itoken = 0
    yamlfilename = "phonopy.yaml"
    while itoken < ntokens:
        token = tokens[itoken]
        if "-symprec" in token:
            itoken += 1
            symprec = float(tokens[itoken])
        elif "-input" in token:
            itoken += 1
            inputfilename = tokens[itoken]
        elif "-nosym" in token:
            symmetry = False
        elif "-output" in token:
            itoken += 1
            outputfilename = tokens[itoken]
        elif "-yaml" in token:
            itoken += 1
            yamlfilename = tokens[itoken]
        elif "-h" in token:
            itoken += 1
            show_usage()
            sys.exit()
        else:
            outputfilename = token
        itoken += 1
    print("Phonopy BORN charges are being converted to PDielec format",file=sys.stderr)
    print(f"The input BORN charges are from: {inputfilename}",file=sys.stderr)
    print(f"The phonopy.yaml is            : {yamlfilename}",file=sys.stderr)
    print(f"The symmetry precision is      : {symprec}",file=sys.stderr)
    print(f"The final output is in         : {outputfilename}",file=sys.stderr)
    if outputfilename and os.path.exists(outputfilename):
        print("")
        print("Error output file already exists",file=sys.stderr)
        sys.exit()
    if not os.path.exists(yamlfilename):
        print("")
        print("Error phonopy.yaml does not exist",file=sys.stderr)
        sys.exit()
    phonon = phonopy.load(yamlfilename)
    primitive = phonon.primitive
    nac = parse_BORN(primitive, symprec=symprec, is_symmetry=symmetry, filename=inputfilename)
    borns = nac["born"]
    epsilon = nac["dielectric"]
    text = "# epsilon and Z* of atoms PDielec/PDGui file format "
    lines = [text]
    lines.append(("%13.8f " * 9) % tuple(epsilon.flatten()))

    with open(outputfilename, "w") if outputfilename else nullcontext(sys.stdout) as fd:
        for z in borns:
            lines.append(("%13.8f " * 9) % tuple(z.flatten()))
        print("\n".join(lines),file=fd)
    print("Command finished")
    
    if __name__ == "__main__":
        main()
