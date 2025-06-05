#!/usr/bin/python
#
# Copyright 2024 John Kendrick
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
# You should have received a copy of the MIT License
# along with this program, if not see https://opensource.org/licenses/MIT
#
import sys, os
home_directory = os.path.join('..','..','..')
examples_directory = os.path.join(home_directory,'Examples')
sys.path.insert(0,home_directory)

import numpy                      as np
# Helper routines
from PDielec.HelperRoutines   import getMaterial
from PDielec.HelperRoutines   import calculatePowderSpectrum
from PDielec.HelperRoutines   import calculateSingleCrystalSpectrum
import PDielec.Utilities      as Utilities
import PDielec.Calculator     as Calculator
from PDielec.Constants  import angs2bohr, amu, wavenumber

def longitudinal_frequencies(file):
    '''Calculate longitudinal frequencies'''
    # Set up the reader for this file
    program = Utilities.find_program_from_name(file)
    reader = Utilities.get_reader(file,program)
    # Initiate reading of the file
    reader.read_output()
    # Extract the information we need for calculating the longitudinal modes
    frequencies  = wavenumber * np.array(reader.frequencies)
    normal_modes = np.array(reader.mass_weighted_normal_modes)
    born_charges = np.array(reader.born_charges)
    masses       = amu * np.array(reader.masses)
    epsilon_inf  = np.array(reader.zerof_optical_dielectric)
    volume       = angs2bohr*angs2bohr*angs2bohr*reader.volume
    # Create a list over which we calculate the longitudinal modes
    qlist = []
    qlist.append( [1,0,0] )
    qlist.append( [0,1,0] )
    qlist.append( [0,0,1] )
    lo_vs = Calculator.longitudinal_modes(frequencies, normal_modes, born_charges, masses, epsilon_inf, volume, qlist, reader)
    print('            TO  ',end='')
    for q in qlist:
        print(q,end='  ')
    print(' ')
    for i,to in enumerate(frequencies):
        print('{:3d} {:10.5f}'.format(i,to/wavenumber),end=' ')
        for  lo in lo_vs:
            print('{:10.5f}'.format(lo[i]/wavenumber),end=' ')
        print(' ')

def main():
    longitudinal_frequencies(os.path.join(examples_directory,'Castep','MgO','phonon.castep'))
    return

if __name__ == '__main__':
    main()
