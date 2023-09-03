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
    possible_epsTypes = ['dft','fpsq','drude-lorentz']
    possible_units = ['cm-1','microns','mu','nm','thz','hz']
    def __init__(self, units='cm-1' ):
        """
        Create a dielectricFunction
        frequency_units        units of frequency (default is cm-1)
        volume_angs            volume in angstrom
        epsilon_infinity       epsilon infinity
        """
        if units not in self.possible_units:
            print('Catastrophic error in DielectricFunction: units not recognised', units)
            exit()
        self.frequency_units      = units
        self.epsilon_infinity     = np.zeros( (3,3) )
        self.volume_angs          = None
        if self.volume_angs:
            self.volume_au = volume*angstrom*angstrom*angstrom

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
        self.epsilon_infinity     = np.array(eps)
        return

    def function(self):
        """Return the function to report the dielectric"""
        return  self.calculate

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


    def dielectriContributionsFromDrude(self, f, frequency, sigma, volume):
        """Calculate the dielectric function for a set of a Drude oscillator
        f is the frequency of the dielectric response in au
        frequency(au), sigmas(au) are the plasma frequency and the width
        The output from this calculation is a complex dielectric tensor
        """
        dielectric = np.zeros((3, 3), dtype=complex)
        unit = identity(3,dtype=complex)
        # Avoid a divide by zero if f is small
        if f <= 1.0e-8:
            f = 1.0e-8
        # Assume that the drude contribution is isotropic
        dielectric = dielectric - unit * frequency*frequency / complex(-f*f, -sigma*f)
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
            dielectric = dielectric + strength / complex((v*v - f*f), -sigma*f)
        return dielectric * (4.0*np.pi/volume)


class ConstantTensor(DielectricFunction):
    def __init__(self, value, units='cm-1'):
        """ A simple constant dielectric function that returns a complex constant
            The value can be a 3x3 tensor or a scalar
            This routine returns a diagonal tensor
        """
        self.value = value
        DielectricFunction.__init__(self,units=units)

    def calculate(self,v):
        return self.value * np.eye(3) + self.epsilon_infinity

class ConstantScalar(DielectricFunction):
    def __init__(self, value, units='cm-1'):
        """ A simple constant dielectric function that returns a complex constant
            The value can be a 3x3 tensor or a scalar
            This routine returns a scalar vale
        """
        self.value = value
        DielectricFunction.__init__(self,units=units)

    def calculate(self,v):
        return self.value

class Tabulate1(DielectricFunction):
    def __init__(self, vs_cm1, permittivities,units='cm-1'):
        """ A dielectric function that interpolates a 1d vector of values
            The initialisation creates all the interpolators
            The results is a 3x3 tensor
        """
        vs = np.array(vs_cm1)
        eps = np.array(permittivities)
        eps_r = np.real(eps)
        eps_i = np.imag(eps)
        self.interpolater = interpolate.InterpolatedUnivariateSpline(vs,eps_r)
        self.interpolatei = interpolate.InterpolatedUnivariateSpline(vs,eps_i)
        DielectricFunction.__init__(self,units=units)

    def calculate(self,v):
        """
        Calculate the dielectric constant from interpolation data
        """
        v_cm1 = self._convert(v)
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater(v_cm1),self.interpolatei(v_cm1))
        eps[1,1] = eps[0,0]
        eps[2,2] = eps[0,0]
        return eps + self.epsilon_infinity

class Tabulate3(DielectricFunction):
    def __init__(self, vs_cm1, epsxx, epsyy, epszz,units='cm-1'):
        """ A dielectric function that interpolates the diagonal components of the permittivity
            The initialisation creates all the interpolators
            The results is a 3x3 tensor
        """
        vs = np.array(vs_cm1)
        eps = [np.array(epsxx), np.array(epsyy), np.array(epszz)]
        epsr = np.zeros(3)
        epsi = np.zeros(3)
        self.interpolater = []
        self.interpolatei = []
        for i in enum(eps):
            epsr[i] = np.real(eps)
            epsi[i] = np.imag(eps)
            self.interpolater.append(interpolate.InterpolatedUnivariateSpline(vs,epsr[i]))
            self.interpolatei.append(interpolate.InterpolatedUnivariateSpline(vs,epsi[i]))
        DielectricFunction.__init__(self,units=units)

    def calculate(self,v):
        """
        Calculate the dielectric constant from interpolation data
        """
        v_cm1 = self._convert(v)
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater[0](v_cm1),self.interpolatei[0](v_cm1))
        eps[1,1] = complex(self.interpolater[1](v_cm1),self.interpolatei[1](v_cm1))
        eps[2,2] = complex(self.interpolater[2](v_cm1),self.interpolatei[2](v_cm1))
        return eps + self.epsilon_infinity

class Tabulate6(DielectricFunction):
    def __init__(self, vs_cm1, epsxx, epsyy, epszz, epsxy, epsxz, epsyz,units='cm-1'):
        """ A dielectric function that interpolates all components of the permittivity
            The initialisation creates all the interpolators
            The results is a 3x3 tensor
        """
        vs = np.array(vs_cm1)
        eps = [np.array(epsxx), np.array(epsyy), np.array(epszz), np.array(epsxy), np.array(epsxz), np.array(epsyz)]
        epsr = np.zeros(3)
        epsi = np.zeros(3)
        self.interpolater = []
        self.interpolatei = []
        for i in enum(eps):
            epsr[i] = np.real(eps)
            epsi[i] = np.imag(eps)
            self.interpolater.append(interpolate.InterpolatedUnivariateSpline(vs,epsr[i]))
            self.interpolatei.append(interpolate.InterpolatedUnivariateSpline(vs,epsi[i]))
        DielectricFunction.__init__(self,units=units)

    def calculate(self,v):
        """
        Calculate the dielectric constant from interpolation data
        """
        v_cm1 = self._convert(v)
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater[0](v_cm1),self.interpolatei[0](v_cm1))
        eps[1,1] = complex(self.interpolater[1](v_cm1),self.interpolatei[1](v_cm1))
        eps[2,2] = complex(self.interpolater[2](v_cm1),self.interpolatei[2](v_cm1))
        eps[1,0] = complex(self.interpolater[3](v_cm1),self.interpolatei[3](v_cm1))
        eps[2,0] = complex(self.interpolater[4](v_cm1),self.interpolatei[4](v_cm1))
        eps[2,1] = complex(self.interpolater[5](v_cm1),self.interpolatei[5](v_cm1))
        eps[0,1] = eps[1,0]
        eps[0,2] = eps[2,0]
        eps[1,2] = eps[2,1]
        return eps + self.epsilon_infinity

class DFT(DielectricFunction):
    def __init__(self,mode_list, mode_frequencies, mode_sigmas, mode_oscillator_strengths, crystal_volume, drude, drude_plasma, drude_sigma,units='cm-1'):
        """
        Calculate the dielectric constant from DFT parameters
        Used internally by the code as it expects input parameters to be in atomic units
        parameter is a list of the following
                     mode_list,
                     mode_frequencies,           (au)
                     mode_sigmas,                (au)
                     mode_oscillator_strengths,  (au)
                     crystal_volume,             (au)
                     drude,                      (true or false)
                     drude_plasma,               (au)
                     drude_sigma                 (au)
        """
        self.mode_list = mode_list
        self.mode_frequencies = mode_frequencies
        self.mode_sigmas = mode_sigmas
        self.mode_oscillator_strengths = mode_oscillator_strengths
        self.crystal_volume = crystal_volume
        self.drude = drude
        self.drude_plasma = drude_plasma
        self.drude_sigma = drude_sigma
        DielectricFunction.__init__(self,units=units)

    def calculate(self,v):
        """
        Calculate the permittivity from DFT data
        """
        v_cm1 = self._convert(v)    # Shouldn't be needed
        v_au = v_cm1 * wavenumber
        eps = self.dielectriContributionsFromModes(v_au, self.mode_list, self.mode_frequencies, 
                                                   self.mode_sigmas, self.mode_oscillator_strengths, 
                                                   self.crystal_volume)
        if self.drude:
            eps = eps + self.dielectricContributionFromDrude(vau, self.drude_plasma, self.drude_sigma, self.crystal_volume)
        return eps + self.epsilon_infinity

class DrudeLorentz(DielectricFunction):
    def __init__(self, vs_cm1, strengths_cm1, sigmas_cm1,units='cm-1'):
        """Use a Drude Lorentz mode to calculate the permittivity
           vs_cm1 are the frequencies of the absorptions
           strengths_cm1 are the strengths
           sigmas are their widths
           If (3,n) arrays are passed through then each array is for xx, yy, zz
           If an (n) array is passed through then the same information is used for xx, yy and zz
        """
        rhombic = False
        if isinstance(vs_cm1[0], list):
            rhombic = True
        if rhombic:
            self.vs_cm1        = vs_cm1
            self.strengths_cm1 = strengths_cm1
            self.sigmas_cm1    = sigmas_cm1
        else:
            self.vs_cm1        = [vs_cm1, vs_cm1, vs_cm1]
            self.strengths_cm1 = [strengths_cm1, strengths_cm1, strengths_cm1]
            self.sigmas_cm1    = [sigmas_cm1, sigmas_cm1, sigmas_cm1]
        DielectricFunction.__init__(self,units=units)
        return

    def calculate(self, v):
        """
        Setup a Drude-Lorentz model for a diagonal permittivity 
        """
        f_cm1 = self._convert(v)
        eps = np.zeros( (3,3), dtype=complex )
        for xyz, (vs, strengths, sigmas) in enumerate(zip(self.vs_cm1,self.strengths_cm1, self.sigmas_cm1)):
            for v, strength, sigma in zip(vs, strengths,sigmas):
                eps[xyz,xyz] += strength * strength / complex((v*v - f_cm1*f_cm1), -sigma*f_cm1)
        return eps + self.epsilon_infinity

class FPSQ(DielectricFunction):
    def __init__(self, omega_tos, gamma_tos, omega_los, gamma_los,units='cm-1'):
        """Use a FPSQ mode to calculate the permittivity
           If (3,n) arrays are passed through then each array is for xx, yy, zz
           If an (n) array is passed through then the same information is used for xx, yy and zz
        """
        rhombic = False
        if isinstance(omega_tos[0], list):
            rhombic = True
        if rhombic:
            self.omega_tos = omega_tos
            self.gamma_tos = gamma_tos
            self.omega_los = omega_los
            self.gamma_los = gamma_los
        else:
            self.omega_tos = [omega_tos, omega_tos, omega_tos]
            self.gamma_tos = [gamma_tos, gamma_tos, gamma_tos]
            self.omega_los = [omega_los, omega_los, omega_los]
            self.gamma_los = [gamma_los, gamma_los, gamma_los]
        DielectricFunction.__init__(self,units=units)
        return

    def calculate(self, v):
        """
        Setup a FPSQ for a diagonal permittivity 
        """
        f_cm1 = self._convert(v)
        eps = np.array(self.epsilon_infinity,dtype=complex)
        for xyz,(omega_tos,gamma_tos,omega_los,gamma_los) in enumerate(zip(self.omega_tos,self.gamma_tos,self.omega_los, self.gamma_los)):
            for omega_to, gamma_to, omega_lo, gamma_lo  in zip(omega_tos,gamma_tos,omega_los,gamma_los):
                contribution = (omega_lo**2 - f_cm1**2 - complex(0,gamma_lo)*f_cm1)/(omega_to**2 - f_cm1**2 - complex(0,gamma_to)*f_cm1)
                # make sure the imaginary component is not negative!
                real_contribution = np.real(contribution)
                imag_contribution = max(-1.0E-12,np.imag(contribution))
                eps[xyz,xyz] *= complex(real_contribution,imag_contribution)
        return eps 

