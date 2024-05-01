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

import sys
sys.path.insert(0, "../../..")
import numpy as np
from PDielec.UnitCell import UnitCell
from PDielec.HelperRoutines   import getMaterial

def main():
    # Read in asparticAcid calculation by CASTEP
    asparticAcid = getMaterial('AsparticAcid/phonon.castep')
    uc = asparticAcid.getCell()
    # bonds etc, are only defined if a new cell is created of whole molecules
    newcell,number_of_molecules,old_atom_order = uc.calculate_molecular_contents()
    print('Number of whole molecules {}'.format(number_of_molecules))
    print('CIF file')
    newcell.write_cif('test.cif')
    density = newcell.getDensity('cm')
    #
    bonds,bondlengths = newcell.get_bonds()
    print('BONDS')
    print('=====')
    for (i,j),bondlength in zip(bonds,bondlengths):
        print('Bond {}-{} = {} Angstrom'.format(i,j,bondlength))
    #
    angles,bondangles = newcell.get_bond_angles()
    print('BOND ANGLES')
    print('===========')
    for (i,j,k),bondangle in zip(angles,bondangles):
        print('Bond angles {}-{}-{} = {} degrees'.format(i,j,k,bondangle))
    #
    torsions,torsionangles = newcell.get_torsions()
    print('TORSION ANGLES')
    print('==============')
    for (i,j,k,l),torsionangle in zip(torsions,torsionangles):
        print('Torsion angles {}-{}-{}-{} = {} degrees'.format(i,j,k,l,torsionangle))
    return

if __name__ == '__main__':
    main()
