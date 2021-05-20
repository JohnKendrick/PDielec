#!/usr/bin/python
#
# Copyright 2015 John Kendrick
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

"""
DielectricFunction provides an interface to different mechanisms for providing dielectric information
to pdielec and pdgui
"""

from __future__ import print_function
import numpy as np
import sys
from PDielec.Constants import wavenumber


class DielectricFunction:
    """Provide an interface to different dielectric functions"""
    possible_epsTypes = ['constant','file','interpolate','dft','fpsq']
    possible_units = ['cm-1','microns','mu','nm','thz','hz']
    def __init__(self, epsType=None, volume=None, filename=None, parameters=None, units='cm-1' ):
        """
        Create a dielectricFunction
        epsType      - specifies the type of function
        filename   - used to read an experimental or calculated dielectric constant tensor
        parameters - used to calculate a dielectric constant from DFT calculations
        volume     - volume in angstrom
        """
        if  epsType not in self.possible_epsTypes:
            print('Catastrophic error in DielectricFunction: type not recognised', epsType)
            exit()
        if units not in self.possible_units:
            print('Catastrophic error in DielectricFunction: units not recognised', units)
            exit()
        self.epsType              = epsType
        self.filename             = filename
        self.parameters           = parameters
        self.frequency_units      = units
        self.epsilon_infinity     = None
        self.volume_angs          = volume
        if epsType == 'file':
            self.readFile()

    def setVolume(self,vol):
        """Set the volume for dielectric calculations"""
        self.volume_angs     = vol
        return

    def setEpsilonInfinity(self,eps):
        """ Set epsilon infinity for dielectric calculations"""
        self.epsilon_infinity     = eps
        return

    def function(self):
        """Return the function to report the dielectric"""
        return  self.calculate

    def readFile(self):
        """
        Read a dielectric function from a file
        """
        print('Reading file ',self.filename)
        exit()
        return
    def convert(self, f, units='cm-1'):
        """
        Convert frequency units to cm-1
        """
        if units == 'cm-1':
            return f                            # return cm-1
        elif units == 'microns' or units == 'mu':
            return 1.0 / (f * 1.0E-6 *1.0E2 )   # convert microns to cm-1
        elif units == 'nm':
            return 1.0 / (f * 1.0E-9 *1.0E2 )   # convert nm to cm-1
        elif units == 'thz':
            return f * 0.0299792458             # convert thz to cm-1
        elif units == 'hz':
            return f * 0.0299792458 * 1.012     # convert thz to cm-1
        else:
            print('Units in DielectricFunction not recognised:',units)
            exit()

    def calculate(self, frequency):
        """
        Calculate the dielectric constant at this frequency
        """
        frequency_cm1 = self.convert(frequency,self.frequency_units)
        if self.epsType == None:
            print('Error in DielectricFunction type is not set')
            exit()
        elif self.epsType == 'constant':
            return self.parameters
        elif self.epsType == 'filename':
            return  self.epsFromFile(frequency_cm1)
        elif self.epsType == 'interpolation':
            return  self.epsFromInterpolation(frequency_cm1)
        elif self.epsType == 'dft':
            return  self.epsFromDFT(frequency_cm1)
        elif self.epsType == 'lorentz':
            return  self.epsFromLorentz(frequency_cm1)
        elif self.epsType == 'fpsq':
            return  self.epsFromFpsq(frequency_cm1)
        else:
            print('Error in DielectricFunction type is not recognised: ', self.epsType)
            exit()
        return

    def epsFromFile(self, f):
        """
        Calculate the dielectric constant from file data
        """
        print('epsFromFile: ',f)
        return

    def epsFromInterpolation(self, f):
        """
        Calculate the dielectric constant from interpolation data
        """
        print('epsFromInterpolation: ',f)
        return

    def epsFromDFT(self, f):
        """
        Calculate the dielectric constant from DFT parameters
        Used internally by the code as it expects input parameters to be in atomic units
        parameter is a list of the following
                     mode_list,
                     mode_frequencies,           (au)
                     mode_sigmas,                (au)
                     mode_oscillator_strengths,  (au)
                     crystal_volume,             (au)
                     epsilon_inf,
                     drude,                      (true or false)
                     drude_plasma,               (au)
                     drude_sigma                 (au)
        """
        mode_list, mode_frequencies, mode_sigmas, mode_oscillator_strengths, crystal_volume, epsilon_inf, drude, drude_plasma, drude_sigma = self.parameters
        v_cm1 = self.convert(f,self.frequency_units)
        v_au = v_cm1 * wavenumber
        eps = self.dielectriContributionsFromModes(v_au, mode_list, mode_frequencies, mode_sigmas, mode_oscillator_strengths, crystal_volume)
        if drude:
            eps = eps + self.dielectricContributionFromDrude(vau, drude_plasma, drude_sigma, crystal_volume)
        return eps + epsilon_inf

    def epsFromIsoLorentz(self, f):
        """
        Setup a Lorentz model for the permittivity for an isotropic system
        """
        f_cm1 = self.convert(f,self.frequency_units)
        f_au  = wavenumber*f_cm1
        eps = np.zeros((3, 3), dtype=complex)
        for mode in self.parameters:
            v_cm, strength, sigma_cm = mode
            v_au     = v_cm * wavenuumber
            sigma_au = sigmas_cm * wavenumber
            eps += dielectric + strength / np.complex((v_au*v_au - f_au*f_au), -sigma_au*f_au)
        return eps * (4.0*np.pi/volume) + self.epsilon_infinity

    def epsFromFpsq(self, f):
        """
        Setup a fpsq model for the permittivity
        """
        #jk print('Setting up fpsq ',self.parameters)
        f_cm1 = self.convert(f,self.frequency_units)
        f_au  = wavenumber*f_cm1
        eps = self.epsilon_infinity
        for xyz,contribution in enumerate(self.parameters):
            for omega_to, gamma_to, omega_lo, gamma_lo  in parameters:
                eps[xyz] *= (omega_to**2 - f_cm1**2 - np.complex(0,gamm_lo)*f_cm1)/(omega_to**2 - f_cm1**2 - np.complex(0,gamm_lo)*f_cm1)
        return eps 


    def dielectriContributionsFromDrude(self, f, frequency, sigma, volume):
        """Calculate the dielectric function for a set of a Drude oscillator
        f is the frequency of the dielectric response in au
        frequency(au), sigmas(au) are the plasma frequency and the width
        The output from this calculation is a complex dielectric tensor
        """
        dielectric = np.zeros((3, 3), dtype=np.complex)
        unit = identity(3,dtype=np.complex)
        # Avoid a divide by zero if f is small
        if f <= 1.0e-8:
            f = 1.0e-8
        # Assume that the drude contribution is isotropic
        dielectric = dielectric - unit * frequency*frequency / np.complex(-f*f, -sigma*f)
        return dielectric * (4.0*np.pi/volume)


    def dielectriContributionsFromModes(self, f, modes, frequencies, sigmas, strengths, volume):
        """Calculate the dielectric function for a set of modes with their own widths and strengths
           f is the frequency of the dielectric response in au
           modes are a list of the modes
           frequencies(au), sigmas(au) and strengths(au) are fairly obvious
           The output from this calculation is a complex dielectric tensor
        """
        dielectric = np.zeros((3, 3), dtype=complex)
        for mode in modes:
            v = frequencies[mode]
            sigma = sigmas[mode]
            strength = strengths[mode].astype(complex)
            dielectric = dielectric + strength / np.complex((v*v - f*f), -sigma*f)
        return dielectric * (4.0*np.pi/volume)


