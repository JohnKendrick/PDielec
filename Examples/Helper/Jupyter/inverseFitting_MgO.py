#!/usr/bin/env python
# coding: utf-8

# # Inverse fitting of permittivities
# Start by importing some useful routines

import numpy                      as np
import sys
import os
sys.path.insert(0, "../../..")     # This should allow the PDielec modules to be found
from PDielec.HelperRoutines   import getMaterial
from PDielec.HelperRoutines   import calculatePowderSpectrum
from PDielec.GUI.FitterTab    import read_experimental_file
from PDielec.Materials        import Tabulated
import matplotlib.pyplot as plt

def mg_mixing(eps0,eps,f):
    '''Calculate the effective mixing permittivity using Maxwell-Garnett
       eps0  is the permittivity of the support
       eps   is the permittivity of the spherical inclusion
       f     is the volume fraction of inclusions'''
    mg = eps0+3*f*eps0*(eps-eps0)/(eps+2*eps0-f*(eps-eps0))
    return mg

def inverse_mg_mixing(eps0,mg,f):
    '''Calculate the crystalline inclusion permittivity which 
       gives the required effective permittivity 
       eps0  is the permittivity of the support
       mg    is the effective permittivity of the mixture
       f     is the volume fraction of inclusions'''
    eps = eps0*(2*eps0*f-2*eps0+f*mg+2*mg)/(2*eps0*f+eps0+f*mg-mg)
    return eps
#
# Read in a permittivity that will be fitted to.
# The Na2SO42/ directory has a file experimental_permittivity.xlsx
# This contains the real and imaginary permittivity from a calculation on Na2(SO4)2
# This will be used to simulate our 'experimental' powder spectrum
# 

file = os.path.join('Na2SO42','experimental_permittivity.xlsx')
file = os.path.join('MgO','experimental_permittivity.xlsx')
#
# Read the real and imaginary components separately and then combine them
#
frequencies,eps_r = read_experimental_file(file,
                                           sheet_name='Powder Real Permittivity',
                                          frequency_column=2,
                                          spectrum_column=3)
frequencies,eps_i = read_experimental_file(file,sheet_name='Powder Imaginary Permittivity',
                                          frequency_column=2,
                                          spectrum_column=3)
reference_eps_1 = np.array(eps_r) + 1j*np.array(eps_i)
#
# The _1 indicates this is a scalar, other conventions for epsilon are _3 for the 3 diagonal components
# For the input to the tabulation routines _6 is used indicating a list with xx, yy, zz, xy, xz, yz components
# The full tensor is indicated with _9
#
matrix = getMaterial('ptfe')
vf = 0.1
#
# Set the electronic and the phonon contributions to the crystal permittivity
#
crystal_eps_1 = []
check_mg_permittivity_1 = []
matrix_permittivity_function = matrix.getPermittivityFunction()
for v,mg_permittivity in zip(frequencies,reference_eps_1):
    matrix_permittivity = matrix_permittivity_function(v)
    crystal_eps_1.append(inverse_mg_mixing(matrix_permittivity,mg_permittivity,vf))
for v,crystal_eps in zip(frequencies,crystal_eps_1):
    matrix_permittivity = matrix_permittivity_function(v)
    check_mg_permittivity_1.append( mg_mixing(matrix_permittivity,crystal_eps,vf))
crystal_eps_3 = [ [ eps, eps, eps ] for eps in crystal_eps_1]
crystal_eps_3 = np.array(crystal_eps_3)
crystal = Tabulated("Guess", vs_cm1=frequencies, permittivities=crystal_eps_3.T, density=2.0)
crystal_epsilon_infinity = np.zeros( (3,3) )
crystal.setEpsilonInfinity(crystal_epsilon_infinity)
# Calculate the effective medium theory permittivity
absorption,emt_eps_1 = calculatePowderSpectrum(frequencies,crystal,matrix,vf)
plt.figure(figsize = (10,5)) # set the size of the figure
plt.plot(frequencies, np.real(reference_eps_1),label='reference') 
plt.plot(frequencies, np.real(emt_eps_1),label='emt') 
plt.plot(frequencies, np.real(check_mg_permittivity_1),label='check') 
#plt.plot(frequencies, np.real(crystal_eps_3.T[0]),label='crystal_xx') 
plt.xlabel('Frequency (cm-1)')
plt.ylabel('Permittivity')
plt.legend(), plt.show()

plt.figure(figsize = (10,5)) # set the size of the figure
plt.plot(frequencies, np.imag(reference_eps_1),label='reference')
plt.plot(frequencies, np.imag(emt_eps_1),label='emt')
plt.plot(frequencies, np.imag(check_mg_permittivity_1),label='check') 
#plt.plot(frequencies, np.imag(crystal_eps_3.T[0]),label='crystal_xx')
plt.xlabel('Frequency (cm-1)')
plt.ylabel('Permittivity')
plt.legend(), plt.show()

