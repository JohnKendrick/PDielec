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
import numpy as np
from PDielec.UnitCell import UnitCell
from PDielec.HelperRoutines   import getMaterial

def main():
    # Read in asparticAcid calculation by CASTEP
    asparticAcid = getMaterial(os.path.join(examples_directory,'Castep','AsparticAcid','phonon.castep'))
    uc = asparticAcid.getCell()
    # bonds etc, are only defined if a new cell is created of whole molecules
    number_of_molecules = uc.calculate_molecular_contents()
    print('Number of whole molecules {}'.format(number_of_molecules))
    print('CIF file')
    uc.write_cif('test.cif')
    density = uc.getDensity('cm')
    #
    bonds,bondlengths = uc.get_bonds()
    print('BONDS')
    print('=====')
    for (i,j),bondlength in zip(bonds,bondlengths):
        print('Bond {}-{} = {} Angstrom'.format(i,j,bondlength))
    #
    angles,bondangles = uc.get_bond_angles()
    print('BOND ANGLES')
    print('===========')
    for (i,j,k),bondangle in zip(angles,bondangles):
        print('Bond angles {}-{}-{} = {} degrees'.format(i,j,k,bondangle))
    #
    torsions,torsionangles = uc.get_torsions()
    print('TORSION ANGLES')
    print('==============')
    for (i,j,k,l),torsionangle in zip(torsions,torsionangles):
        print('Torsion angles {}-{}-{}-{} = {} degrees'.format(i,j,k,l,torsionangle))
    return

if __name__ == '__main__':
    main()
