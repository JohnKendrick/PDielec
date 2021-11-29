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
from PDielec.Constants import wavenumber, angstrom, speed_light_si
from scipy import interpolate


class DielectricFunction:
    """Provide an interface to different dielectric functions"""
    possible_epsTypes = ['constant','constant_tensor','file','interpolate_3','interpolate_6','dft','fpsq','drude-lorentz']
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
        if self.volume_angs:
            self.volume_au = volume*angstrom*angstrom*angstrom
        if epsType == 'file':
            self._readFile()
        elif epsType == 'interpolate_6':
            self._generate6Interpolators()
        elif epsType == 'interpolate_3':
            self._generate3Interpolators()
        elif epsType == 'interpolate':
            self._generate3Interpolators()

    def setUnits(self,units):
        """Set the volume for dielectric calculations"""
        if units not in self.possible_units:
            print('Catastrophic error in DielectricFunction: units not recognised', units)
            exit()
        self.frequency_units      = units
        return

    def setVolume(self,volume):
        """Set the volume for dielectric calculations"""
        self.volume_angs     = volume
        self.volume_au       = volume*angstrom*angstrom*angstrom
        return

    def setEpsilonInfinity(self,eps):
        """ Set epsilon infinity for dielectric calculations"""
        self.epsilon_infinity     = eps
        return

    def function(self):
        """Return the function to report the dielectric"""
        return  self.calculate

    def _readFile(self):
        """
        Read a dielectric function from a file
        """
        print('Reading file ',self.filename)
        exit()
        return

    def _generate6Interpolators(self):
        # If the interpolators do not exist we have to create them from the parameters
        vs = []
        rxxs = []
        ryys = []
        rzzs = []
        ixxs = []
        iyys = []
        izzs = []
        rxys = []
        rxzs = []
        ryzs = []
        ixys = []
        ixzs = []
        iyzs = []
        for v_cm1, epsrxx, epsixx, epsryy, epsiyy, epsrzz, epsizz, epsrxy, epsixy, epsrxz, epsixz, epsryz, epsiyz in self.parameters:
            vs.append(v_cm1)
            epsixx = min(0.0,epsixx)
            epsiyy = min(0.0,epsiyy)
            epsizz = min(0.0,epsizz)
            epsixy = min(0.0,epsixy)
            epsixz = min(0.0,epsixz)
            epsiyz = min(0.0,epsiyz)
            rxxs.append(epsrxx)
            ryys.append(epsryy)
            rzzs.append(epsrzz)
            rxys.append(epsrzz)
            rxzs.append(epsrzz)
            ryzs.append(epsrzz)
            ixxs.append(epsixx)
            iyys.append(epsiyy)
            izzs.append(epsizz)
            ixys.append(epsizz)
            ixzs.append(epsizz)
            iyzs.append(epsizz)
        self.interpolaterxx = interpolate.InterpolatedUnivariateSpline(vs,np.array(rxxs))
        self.interpolateixx = interpolate.InterpolatedUnivariateSpline(vs,np.array(ixxs))
        self.interpolateryy = interpolate.InterpolatedUnivariateSpline(vs,np.array(ryys))
        self.interpolateiyy = interpolate.InterpolatedUnivariateSpline(vs,np.array(iyys))
        self.interpolaterzz = interpolate.InterpolatedUnivariateSpline(vs,np.array(rzzs))
        self.interpolateizz = interpolate.InterpolatedUnivariateSpline(vs,np.array(izzs))
        self.interpolaterxy = interpolate.InterpolatedUnivariateSpline(vs,np.array(rxys))
        self.interpolateixy = interpolate.InterpolatedUnivariateSpline(vs,np.array(ixys))
        self.interpolaterxz = interpolate.InterpolatedUnivariateSpline(vs,np.array(rxzs))
        self.interpolateixz = interpolate.InterpolatedUnivariateSpline(vs,np.array(ixzs))
        self.interpolateryz = interpolate.InterpolatedUnivariateSpline(vs,np.array(ryzs))
        self.interpolateiyz = interpolate.InterpolatedUnivariateSpline(vs,np.array(iyzs))
        self.number_of_interpolators = 6

    def _generate3Interpolators(self):
        # If the interpolators do not exist we have to create them from the parameters
        vs = []
        rxxs = []
        ryys = []
        rzzs = []
        ixxs = []
        iyys = []
        izzs = []
        for v_cm1, epsrxx, epsixx, epsryy, epsiyy, epsrzz, epsizz in self.parameters:
            vs.append(v_cm1)
            epsixx = min(0.0,epsixx)
            epsiyy = min(0.0,epsiyy)
            epsizz = min(0.0,epsizz)
            rxxs.append(epsrxx)
            ryys.append(epsryy)
            rzzs.append(epsrzz)
            ixxs.append(epsixx)
            iyys.append(epsiyy)
            izzs.append(epsizz)
        self.interpolaterxx = interpolate.InterpolatedUnivariateSpline(vs,np.array(rxxs))
        self.interpolateixx = interpolate.InterpolatedUnivariateSpline(vs,np.array(ixxs))
        self.interpolateryy = interpolate.InterpolatedUnivariateSpline(vs,np.array(ryys))
        self.interpolateiyy = interpolate.InterpolatedUnivariateSpline(vs,np.array(iyys))
        self.interpolaterzz = interpolate.InterpolatedUnivariateSpline(vs,np.array(rzzs))
        self.interpolateizz = interpolate.InterpolatedUnivariateSpline(vs,np.array(izzs))
        self.number_of_interpolators = 3

    def _convert(self, f):
        """
        Convert the incoming frequency units to cm-1
        Internally it is assumed that everything is in cm-1
        """
        units = self.frequency_units
        if units == 'cm-1':
            result = f                            # return cm-1
        elif units == 'microns' or units == 'mu':
            result = 1.0 / (f * 1.0E-6 *1.0E2 )           # convert microns to cm-1
        elif units == 'nm':
            result = 1.0 / (f * 1.0E-9 *1.0E2 )           # convert nm to cm-1
        elif units == 'thz':
            result = f*1.0e12 / (speed_light_si *1.0e2)   # convert hz to cm-1
        elif units == 'hz':
            result = f / (speed_light_si *1.0e2)          # convert hz to cm-1
        else:
            print('Units in DielectricFunction not recognised:',units)
            exit()
        return result

    def calculate(self, frequency):
        """
        Calculate the dielectric constant at this frequency
        The frequency units are defined by the convert function
        """
        frequency_cm1 = self._convert(frequency)
        if self.epsType == None:
            print('Error in DielectricFunction type is not set')
            exit()
        elif self.epsType == 'constant':
            return self.parameters
        elif self.epsType == 'constant_tensor':
            return self.parameters*np.eye(3)
        elif self.epsType == 'filename':
            return  self._epsFromFile(frequency_cm1)
        elif self.epsType == 'interpolate_6':
            return  self._epsFromInterpolate6(frequency_cm1)
        elif self.epsType == 'interpolate_3':
            return  self._epsFromInterpolate3(frequency_cm1)
        elif self.epsType == 'interpolate':
            return  self._epsFromInterpolate3(frequency_cm1)
        elif self.epsType == 'dft':
            return  self._epsFromDFT(frequency_cm1)
        elif self.epsType == 'drude-lorentz':
            return  self._epsFromDrudeLorentz(frequency_cm1)
        elif self.epsType == 'fpsq':
            return  self._epsFromFpsq(frequency_cm1)
        else:
            print('Error in DielectricFunction type is not recognised: ', self.epsType)
            exit()
        return

    def _epsFromFile(self, v_cm1):
        """
        Calculate the dielectric constant from file data
        """
        print('_epsFromFile: ',v_cm1)
        return

    def _epsFromInterpolate3(self, v_cm1):
        """
        Calculate the dielectric constant from interpolation data
        """
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolaterxx(v_cm1),self.interpolateixx(v_cm1))
        eps[1,1] = complex(self.interpolateryy(v_cm1),self.interpolateiyy(v_cm1))
        eps[2,2] = complex(self.interpolaterzz(v_cm1),self.interpolateizz(v_cm1))
        return eps + self.epsilon_infinity

    def _epsFromInterpolate6(self, v_cm1):
        """
        Calculate the dielectric constant from interpolation data
        """
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolaterxx(v_cm1),self.interpolateixx(v_cm1))
        eps[1,1] = complex(self.interpolateryy(v_cm1),self.interpolateiyy(v_cm1))
        eps[2,2] = complex(self.interpolaterzz(v_cm1),self.interpolateizz(v_cm1))
        eps[0,1] = complex(self.interpolaterxy(v_cm1),self.interpolateixy(v_cm1))
        eps[0,2] = complex(self.interpolaterxz(v_cm1),self.interpolateixz(v_cm1))
        eps[1,2] = complex(self.interpolateryz(v_cm1),self.interpolateiyz(v_cm1))
        eps[1,0] = eps[0,1]
        eps[2,0] = eps[0,2]
        eps[2,1] = eps[1,2]
        return eps + self.epsilon_infinity

    def _epsFromDFT(self, v_cm1):
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
        mode_list, mode_frequencies, mode_sigmas, mode_oscillator_strengths, crystal_volume, drude, drude_plasma, drude_sigma = self.parameters
        v_au = v_cm1 * wavenumber
        eps = self.dielectriContributionsFromModes(v_au, mode_list, mode_frequencies, mode_sigmas, mode_oscillator_strengths, crystal_volume)
        if drude:
            eps = eps + self.dielectricContributionFromDrude(vau, drude_plasma, drude_sigma, crystal_volume)
        return eps + self.epsilon_infinity

    def _epsFromDrudeLorentz(self, f_cm1):
        """
        Setup a Drude-Lorentz model for a diagonal permittivity 
        """
        eps = np.zeros( (3,3), dtype=complex )
        for xyz,contribution in enumerate(self.parameters):
            for mode in contribution:
                v_cm1, strength_cm1, sigma_cm1 = mode
                eps[xyz,xyz] += strength_cm1 * strength_cm1 / np.complex((v_cm1*v_cm1 - f_cm1*f_cm1), -sigma_cm1*f_cm1)
        return eps + self.epsilon_infinity

    def _epsFromFpsq(self, f_cm1):
        """
        Setup a fpsq model for the permittivity
        """
        #jk print('Setting up fpsq ',self.parameters)
        eps = np.array(self.epsilon_infinity,dtype=complex)
        for xyz,contribution in enumerate(self.parameters):
            for omega_to, gamma_to, omega_lo, gamma_lo  in contribution:
                contribution = (omega_lo**2 - f_cm1**2 - np.complex(0,gamma_lo)*f_cm1)/(omega_to**2 - f_cm1**2 - np.complex(0,gamma_to)*f_cm1)
                # make sure the imaginary component is not negative!
                real_contribution = np.real(contribution)
                imag_contribution = max(-1.0E-12,np.imag(contribution))
                eps[xyz,xyz] *= np.complex(real_contribution,imag_contribution)
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


