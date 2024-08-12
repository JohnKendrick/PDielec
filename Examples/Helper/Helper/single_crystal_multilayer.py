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
# Utility routines
from PDielec.GUI.SingleCrystalLayer import SingleCrystalLayer

def singleCrystalTest():
    '''Test single crystal calculation'''
    # Define the materials
    air1  = getMaterial('air')
    ptfe = getMaterial('ptfe')
    kbr  = getMaterial('kbr')
    sapphire  = getMaterial('Sapphire')
    air2  = getMaterial('air')
    # Prepare the layers
    layers = []
    # Add a substrate
    layers.append(SingleCrystalLayer(air1,thickness=1.0,thicknessUnit='um'))
    layers.append(SingleCrystalLayer(kbr,thickness=1.0,thicknessUnit='um'))
    layers.append(SingleCrystalLayer(ptfe,thickness=1.0,thicknessUnit='um'))
    layers.append(SingleCrystalLayer(sapphire,hkl=[0,1,0],thickness=1.0,thicknessUnit='um'))
    layers.append(SingleCrystalLayer(air2,thickness=1.0,thicknessUnit='um'))
    frequencies_cm1 = np.array([200.0])
    incident_angle = 80.0
    global_azimuthal_angle = 0.0
    (r_sm, t_sm, a_sm) = calculateSingleCrystalSpectrum(frequencies_cm1,layers,incident_angle, global_azimuthal_angle, method='Scattering matrix')
    print('')
    print('Scattering matrix')
    print('R',r_sm[0], r_sm[1])
    print('T',t_sm[0], t_sm[1])
    print('A',a_sm[0], a_sm[1])
    print('S',r_sm[0]+t_sm[0]+a_sm[0], r_sm[1]+t_sm[1]+a_sm[1])
    (r_tm, t_tm, a_tm) = calculateSingleCrystalSpectrum(frequencies_cm1,layers,incident_angle, global_azimuthal_angle, method='Transfer matrix')
    print('')
    print('Transfer matrix')
    print('R',r_tm[0], r_tm[1])
    print('T',t_tm[0], t_tm[1])
    print('A',a_tm[0], a_tm[1])
    print('S',r_tm[0]+t_tm[0]+a_tm[0], r_tm[1]+t_tm[1]+a_tm[1])

def main():
    singleCrystalTest()
    return

if __name__ == '__main__':
    main()
