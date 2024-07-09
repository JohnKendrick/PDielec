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
import numpy                      as np
# Helper routines
from PDielec.HelperRoutines   import getMaterial
from PDielec.HelperRoutines   import calculatePowderSpectrum
from PDielec.HelperRoutines   import calculateSingleCrystalSpectrum
# Utility routines
from PDielec.GUI.SingleCrystalLayer import SingleCrystalLayer

def powderTest():
    '''Test powder calculation'''
    frequencies_cm1 = np.arange( 0, 200, 0.2 )
    matrix = getMaterial('ptfe')
    dielectric = getMaterial('Sapphire')
    method = 'Maxwell-Garnett' 
    shape = 'Sphere'
    volume_fraction = 0.1
    absorption,permittivity = calculatePowderSpectrum(frequencies_cm1,dielectric, matrix, volume_fraction)
    print('')
    print('Powder absorption results')
    print('      freq', ' absorption (cm-1)')
    for f,a in zip(frequencies_cm1,absorption):
        print(f'{f:10.2f}{a:12.5f}')
    
def singleCrystalTest():
    '''Test single crystal calculation'''
    # Define the materials
    air = getMaterial('air')
    ptfe = getMaterial('ptfe')
    Sapphire = getMaterial('Sapphire')
    # Prepare the layers
    layers = []
    # Add a substrate
    layers.append(SingleCrystalLayer(ptfe,thickness=1.0,thicknessUnit='um'))
     # Add the dielectric layer
    layers.append(SingleCrystalLayer(Sapphire,hkl=[0,0,1],azimuthal=0.0,thickness=1.0,thicknessUnit='um',
                                     incoherentOption='Coherent'))
     # Add the superstrate
    layers.append(SingleCrystalLayer(air,thickness=1.0,thicknessUnit='um'))
    frequencies_cm1 = np.arange( 0, 200, 0.2 )
    incident_angle = 80.0
    global_azimuthal_angle = 0.0
    (reflectance, transmittance, absorptance) = calculateSingleCrystalSpectrum(frequencies_cm1,layers,incident_angle, global_azimuthal_angle, method='Scattering matrix')
    rps = reflectance[0]
    rss = reflectance[1]
    tps = transmittance[0]
    tss = transmittance[1]
    bps = absorptance[0]
    bss = absorptance[1]
    # Print single crystal results
    print('')
    print('Single Crystal Results')
    print('      freq', '      r_p', '      r_s', '      t_p', '      t_s', '      a_p', '      a_s')
    for rs,rp,ts,tp,bs,bp,f in zip(rss,rps,tss,tps,bss,bps,frequencies_cm1):
        print(f'{f:10.2f}{rp:10.5f}{rs:10.5f}{tp:10.5f}{ts:10.5f}{bp:10.5f}{bs:10.5f}')

def main():
    singleCrystalTest()
    powderTest()
    return

if __name__ == '__main__':
    main()
