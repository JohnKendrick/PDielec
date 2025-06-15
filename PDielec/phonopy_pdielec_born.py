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
from phonopy.file_IO import parse_BORN


def show_usage():
    """Show the usage message."""
    print("")
    print("USAGE:")
    print("phonopy-pdielec-born [--symprec 1.0E-5] [-input BORN] [-nosym] outputfilename")
    print("")
    print(" -symprec Sets the precision used in determing the symmetry (default is 1.0E-5)")
    print(" -input   Sets the input BORN charges file (Defaults is BORN)")
    print(" -nosym   Stops symmetry being used in parsing the BORN file")

def main():
    """Define the main routine for the phonopy-pdielec-born command.

    The routine makes extensive use of Phonopy's API
    It reads the phonopy.yaml file to determine the primitive cell
    then it reads the BORN file, before writing out the permitivity and Born charges.
    The output file is formatted in a similar way to the BORN file, but
    all the atoms/ions are present, not just the symmetry unique ones.

    Options:
    -symprec 1.0E-5       Sets the precision used in determing the symmetry (default is 1.0E-5)
    -input   BORN         Sets the input BORN charges file (Defaults is BORN)
    -nosym                Stops symmetry being used in parsing the BORN file
    """
    if len(sys.argv) <= 1:
        show_usage()
        sys.exit()

    symprec = 1.e-5
    inputfilename = "BORN"
    outputfilename = None
    symmetry = True
    tokens = sys.argv[1:]
    ntokens = len(tokens)
    itoken = 0
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
        else:
            outputfilename = token
        itoken += 1
    print("Phonopy BORN charges are being converted to PDielec format")
    print(f"The input BORN charges are from: {inputfilename}")
    print(f"The symmetry precision is      : {symprec}")
    print(f"The final output is in         : {outputfilename}")
    if outputfilename is None:
        print("")
        print("No output file has been given")
        show_usage()
        sys.exit()
    if os.path.exists(outputfilename):
        print("")
        print("Error output file already exists")
        sys.exit()
    phonon = phonopy.load("phonopy.yaml")
    primitive = phonon.primitive
    nac = parse_BORN(primitive, symprec=symprec, is_symmetry=symmetry, filename=inputfilename)
    borns = nac["born"]
    epsilon = nac["dielectric"]
    text = "# epsilon and Z* of atoms PDielec/PDGui file format "
    lines = [text]
    lines.append(("%13.8f " * 9) % tuple(epsilon.flatten()))
    with open(outputfilename,"w") as fd:
        for z in borns:
            lines.append(("%13.8f " * 9) % tuple(z.flatten()))
        print("\n".join(lines),file=fd)
    print("Command finished")
    
    if __name__ == "__main__":
        main()
