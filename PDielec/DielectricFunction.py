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

'''
DielectricFunction provides an interface to different mechanisms for providing dielectric information
to pdielec and pdgui
'''

from __future__ import print_function
import numpy as np
import sys
from PDielec.Constants import wavenumber, angstrom, speed_light_si
import PDielec.Calculator as Calculator
from scipy import interpolate


class DielectricFunction:
    '''
    Provide a base class to different dielectric functions

    Attributes
    ----------
    frequency_units : str       
                           can be cm-1, microns, mu,nm, thz or hz (default is cm-1)
                           supplied frequencies are converyed to cm-1
    epsilon_infinity: 3x3 tensor (numpy array)
                           Optical permittivity
    self.volume: float
    self.volume_au: float  
                           Volume of crystal unit cells 
    self.vs_cm1 : np.array of floats
                           The frequencies of any tabulation
    self.isScalarFunction  True for a scalar function

    Methods
    -------
    setUnits(units)         Changes the frequency units
    setVolume(volume)       Change the unit cell volume 
    setEpsilonInfinity(eps) Change epsilon infinity
    function()              Return the calculate function for this object
    calculate(v)            Return the value of the permittivy at frequency v
    dielectriContributionsFromDrude(self, f, frequency, sigma, volume):
                            Calculates a permittivity from a Drude Model
    dielectriContributionsFromModes(self, f, modes, frequencies, sigmas, strengths, volume):
                            Calculates a permittivity from a Drude-Lorentz model
    _convert()              Converts frequency units
    '''
    possible_epsTypes = ['dft','fpsq','drude-lorentz','sellmeier']
    possible_units = ['cm-1','microns','mu','nm','thz','hz']

    def __init__(self, units='cm-1' ):
        '''
        Parameters
        ----------
        units : str
             The units of frequency: cm-1,microns,mu,nm,thz or hz (default is cm-1)
        '''
        if units not in self.possible_units:
            print('Catastrophic error in DielectricFunction: units not recognised', units)
            exit()
        self.isScalarFunction             = None
        self.vs_cm1               = 0
        self.frequency_units      = units
        self.epsilon_infinity     = np.zeros( (3,3) )
        self.volume_angs          = None
        if self.volume_angs:
            self.volume_au = volume*angstrom*angstrom*angstrom

    def isTensor(self):
        '''Returns true if the dielectric function returns a tensor'''
        return not self.isScalarFunction

    def isScalar(self):
        '''Returns true if the dielectric function returns a scalar value'''
        return self.isScalarFunction

    def setUnits(self,units):
        '''Set the units of frequency to be used in calls to calculate()

        Parameters
        ----------
        units : str
             The units of frequency: cm-1,microns,mu,nm,thz or hz (default is cm-1)
        '''
        if units not in self.possible_units:
            print('Catastrophic error in DielectricFunction: units not recognised', units)
            exit()
        self.frequency_units      = units
        return

    def getUnits(self):
        '''get units of frequency to be used in calls

        '''
        return self.frequency_units

    def setVolume(self,volume):
        '''Set the volume for dielectric calculations

        Parameters
        ----------
        volume : float
             The volume in angs^3
        '''
        self.volume_angs     = volume
        self.volume_au       = volume*angstrom*angstrom*angstrom
        return

    def setEpsilonInfinity(self,eps):
        ''' 
        Set epsilon infinity for dielectric calculations

        Parameters
        ----------
        eps : float
             The value of eps_infinity
        '''
        self.epsilon_infinity     = np.array(eps)
        return

    def function(self):
        '''
        Return
        ------
        Return the function to used to calculate the permittivity
        '''
        return  self.calculate

    def _convert(self, f):
        ''' 
        Internal routine to convert a frequency in to cm-1

        Internally it is assumed that everything is in cm-1.  So
        when a frequency is provided it may have to be converted.

        Parameters
        ----------
        f : float
            The frequency

        Result
        ------
        The frequency in cm-1
        '''
        units = self.frequency_units
        if units == 'cm-1':
            result = f                            # return cm-1
        elif units == 'microns' or units == 'mu':
            result = 1.0 / (f * 1.0E-6 *1.0E2 )           # convert microns to cm-1
        elif units == 'nm':
            result = 1.0 / (f * 1.0E-9 *1.0E2 )           # convert nm to cm-1
        elif units == 'thz':
            result = f*1.0e12 / (speed_light_si *1.0e2)   # convert thz to cm-1
        elif units == 'hz':
            result = f / (speed_light_si *1.0e2)          # convert hz to cm-1
        else:
            print('Units in DielectricFunction not recognised:',units)
            exit()
        return result

    def getLowestFrequency(self):
        ''' 
        Return the lowest tabulated frequency in cm-1

        Result
        ------
        The lowest tabulated frequency in cm-1
        '''
        lowestFrequency = np.min(self.vs_cm1)
        return lowestFrequency


    def getHighestFrequency(self):
        ''' 
        Return the highest tabulated frequency in cm-1

        Result
        ------
        The highest tabulated frequency in cm-1
        '''
        highestFrequency = np.max(self.vs_cm1)
        return highestFrequency


    def dielectriContributionsFromDrude(self, f, frequency, sigma, volume):
        ''' 
        Calculate the dielectric function for a set of a Drude oscillator

        Parameters
        ----------
        f : float
            the frequency of the dielectric response
        frequency : float
            The frequency of the drude oscillator
        sigma : float
            The width of the oscillator

        Returns
        -------
        A complex scalar dielectric
        '''
        dielectric = np.zeros((3, 3), dtype=complex)
        unit = identity(3,dtype=complex)
        # Avoid a divide by zero if f is small
        if f <= 1.0e-8:
            f = 1.0e-8
        # Assume that the drude contribution is isotropic
        dielectric = dielectric - unit * frequency*frequency / complex(-f*f, -sigma*f)
        return dielectric * (4.0*np.pi/volume)


    def dielectriContributionsFromModes(self, f, modes, frequencies, sigmas, strengths, volume):
        ''' 
        Calculate the dielectric function for a set of modes

        Parameters
        ----------
        mode : list of ints
             A list of the modes to be considered
        frequencies : list of floats
             The frequency of each mode
        sigmas : list of floats
            The width of each mode
        sigmas : list of floats
            The strengths of each mode

        Returns
        -------
        A complex scalar dielectric
        '''
        dielectric = np.zeros((3, 3), dtype=complex)
        for mode in modes:
            v = frequencies[mode]
            sigma = sigmas[mode]
            strength = strengths[mode].astype(complex)
            dielectric = dielectric + strength / complex((v*v - f*f), -sigma*f)
        return dielectric * (4.0*np.pi/volume)

    def print(self, v1_cm1, v2_cm1, v_inc=1,diagonal_only=True,file=None):
        '''Print the permittivity over a range of frequencies'''
        fd = sys.stdout
        closeFlag = False
        if file is not None:
            closeFlag = True
            fd = open(file,'w')
        c = ','
        if diagonal_only:
            print('f,eps00.real,epss11.real,eps22.real,eps00.imag,eps11.imag,eps22,imag',file=fd)
        else:
            print('f,eps00,epss11,eps22,eps01,eps02,,eps12',file=fd)
        units = self.getUnits()
        self.setUnits('cm-1')
        for v in np.arange(v1_cm1, v2_cm1, v_inc):
            epsilon = self.calculate(v)
            if diagonal_only:
                print(v,c,epsilon[0,0].real,c,epsilon[1,1].real,c,epsilon[2,2].real,c,epsilon[0,0].imag,c,epsilon[1,1].imag,c,epsilon[2,2].imag,file=fd)
            else:
                print(v,c,epsilon[0,0],c,epsilon[1,1],c,epsilon[2,2],c,epsilon[0,1],c,epsilon[0,2],c,epsilon[1,2],file=fd)
        self.setUnits(units)
        if closeFlag:
            fd.close()
        return

class ConstantTensor(DielectricFunction):
    '''
    A class for a complex constant tensor dielectric function

    Inherits from DielectricFunction
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, value, units='cm-1'):
        ''' 
            Parameters
            ----------
            value : 3x3 numpy array 
                    Initialise the permittivity with this tensor
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.value = value
        self.isScalarFunction = False

    def calculate(self,v):
        ''' 
            Parameters
            ----------
            v : float
                The frequency in internal units

           Returns
           -------
           The permittivity at frequency v as a 3x3 tensor
        '''
        return self.value * np.eye(3) + self.epsilon_infinity

class ConstantScalar(DielectricFunction):
    '''
    A class for a complex constant scalar dielectric function

    Inherits from DielectricFunction
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, value, units='cm-1'):
        ''' 
            Parameters
            ----------
            value : 3x3 numpy array 
                    Initialise the permittivity with this tensor
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.value = value
        self.isScalarFunction = True

    def calculate(self,v):
        '''
            Parameters
            ----------
            v : float
                The frequency in internal units

           Returns
           -------
           The permittivity at frequency v scalar
        '''
        return self.value

class TabulateScalar(DielectricFunction):
    '''
    A class for a complex tabulated scalar dielectric function

    Inherits from DielectricFunction
    Represents an isotropic material dielectric function
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, vs_cm1, permittivities,units='cm-1'):
        ''' 
            Parameters
            ----------
            vs_cm1 : list of floats
                     frequencies in cm-1 of the following permittivities
            permittivities : list of complex
                     The permittivities at each of the frequencies
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.isScalarFunction = True
        self.vs_cm1 = np.array(vs_cm1)
        eps = np.array(permittivities)
        eps_r = np.real(eps)
        eps_i = np.imag(eps)
        self.interpolater = interpolate.CubicSpline(self.vs_cm1,eps_r)
        self.interpolatei = interpolate.CubicSpline(self.vs_cm1,eps_i)

    def calculate(self,v):
        '''
            Parameters
            ----------
            v : float
                The frequency in internal units

           Returns
           -------
           The permittivity at frequency v as a scalar
        '''
        v_cm1 = self._convert(v)
        eps = complex(self.interpolater(v_cm1),self.interpolatei(v_cm1))
        return eps

class Tabulate1(DielectricFunction):
    '''
    A class for a complex constant tensor dielectric function

    Inherits from DielectricFunction
    Represents an isotropic dielectric function
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, vs_cm1, permittivities,units='cm-1'):
        ''' 
            Parameters
            ----------
            vs_cm1 : list of floats
                     frequencies in cm-1 of the following permittivities
            permittivities : list of complex
                     The permittivities at each of the frequencies
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.isScalarFunction = False
        self.vs_cm1 = np.array(vs_cm1)
        eps = np.array(permittivities)
        eps_r = np.real(eps)
        eps_i = np.imag(eps)
        self.interpolater = interpolate.CubicSpline(self.vs_cm1,eps_r)
        self.interpolatei = interpolate.CubicSpline(self.vs_cm1,eps_i)

    def calculate(self,v):
        '''
            Parameters
            ----------
            v : float
                The frequency in internal units

           Returns
           -------
           The isotropic permittivity at frequency v as a 3x3 tensor
        '''
        v_cm1 = self._convert(v)
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater(v_cm1),self.interpolatei(v_cm1))
        eps[1,1] = eps[0,0]
        eps[2,2] = eps[0,0]
        return eps + self.epsilon_infinity

class Tabulate3(DielectricFunction):
    '''
    A class for a complex constant tensor dielectric function

    Inherits from DielectricFunction
    Represents an orthrhombic dielectric function
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, vs_cm1, epsxx, epsyy, epszz,units='cm-1'):
        ''' 
            Parameters
            ----------
            vs_cm1 : list of floats
                     frequencies in cm-1 of the following permittivities
            epsxx : list of complex
                     The eps(xx) permittivities at each of the frequencies
            epsyy : list of complex
                     The eps(yy) permittivities at each of the frequencies
            epszz : list of complex
                     The eps(zz) permittivities at each of the frequencies
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.isScalarFunction = False
        self.vs_cm1 = np.array(vs_cm1)
        eps = [np.array(epsxx), np.array(epsyy), np.array(epszz)]
        epsr = np.zeros(3)
        epsi = np.zeros(3)
        self.interpolater = []
        self.interpolatei = []
        for i in enum(eps):
            epsr[i] = np.real(eps)
            epsi[i] = np.imag(eps)
            self.interpolater.append(interpolate.CubicSpline(self.vs_cm1,epsr[i]))
            self.interpolatei.append(interpolate.CubicSpline(self.vs_cm1,epsi[i]))

    def calculate(self,v):
        '''
            Parameters
            ----------
            v : float
                The frequency in internal units

           Returns
           -------
           The permittivity at frequency v as a diagonal 3x3 tensor
        '''
        v_cm1 = self._convert(v)
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater[0](v_cm1),self.interpolatei[0](v_cm1))
        eps[1,1] = complex(self.interpolater[1](v_cm1),self.interpolatei[1](v_cm1))
        eps[2,2] = complex(self.interpolater[2](v_cm1),self.interpolatei[2](v_cm1))
        return eps + self.epsilon_infinity

class Tabulate6(DielectricFunction):
    '''
    A class for a complex constant tensor dielectric function

    Inherits from DielectricFunction
    Represents an non-isotropic dielectric function
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, vs_cm1, epsxx, epsyy, epszz, epsxy, epsxz, epsyz,units='cm-1'):
        ''' 
            Parameters
            ----------
            vs_cm1 : list of floats
                    The frequencies of the tabulated dielectric 
            epsxx : list of complex
                     The eps(xx) permittivities at each of the frequencies
            epsyy : list of complex
                     The eps(yy) permittivities at each of the frequencies
            epszz : list of complex
                     The eps(zz) permittivities at each of the frequencies
            epsxy : list of complex
                     The eps(xy) permittivities at each of the frequencies
            epsxz : list of complex
                     The eps(xz) permittivities at each of the frequencies
            epsyz : list of complex
                     The eps(yz) permittivities at each of the frequencies
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.isScalarFunction = False
        self.vs_cm1 = np.array(vs_cm1)
        eps = [np.array(epsxx), np.array(epsyy), np.array(epszz), np.array(epsxy), np.array(epsxz), np.array(epsyz)]
        epsr = np.zeros(3)
        epsi = np.zeros(3)
        self.interpolater = []
        self.interpolatei = []
        for i in enum(eps):
            epsr[i] = np.real(eps)
            epsi[i] = np.imag(eps)
            self.interpolater.append(interpolate.CubicSpline(self.vs_cm1,epsr[i]))
            self.interpolatei.append(interpolate.CubicSpline(self.vs_cm1,epsi[i]))

    def calculate(self,v):
        '''
            Parameters
            ----------
            v : float
                The frequency in the units set by setUnit()

           Returns
           -------
           The full symmetric permittivity tensor at frequency v as a 3x3 tensor
        '''
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
    '''
    A class for a complex constant tensor dielectric function

    Inherits from DielectricFunction
    Suitable for use in calculating the permittivity from DFT calculations
    Internal units for this function are atomic units
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self,mode_list, mode_frequencies, mode_sigmas, mode_oscillator_strengths, crystal_volume, drude, drude_plasma, drude_sigma):
        ''' 
            Parameters
            ----------
            mode_list : list of ints
                    A list of active modes
            mode_frequencies : list of floats
                    All mode frequencies
            mode_sigmas : list of floats
                    All mode frequencies
            mode_oscillator_strengths : list of floats
                    All mode oscillator strengths
            crystal_volume,             
                    The crystal volume 
            drude : boolean
                    True if a drude oscillator is to be included                      
            drude_plasma : float               
                     Drude plasma frequency
            drude_sigma  : float
                     Drude plasma width
            All parameters are given in atomic units
        '''
        DielectricFunction.__init__(self)
        self.isScalarFunction = False
        self.mode_list = mode_list
        self.mode_frequencies = mode_frequencies
        self.mode_sigmas = mode_sigmas
        self.mode_oscillator_strengths = mode_oscillator_strengths
        self.crystal_volume = crystal_volume
        self.drude = drude
        self.drude_plasma = drude_plasma
        self.drude_sigma = drude_sigma

    def calculate(self,v):
        '''
            Parameters
            ----------
            v : float
                The frequency in cm-1

           Returns
           -------
           The permittivity at frequency v_cm1 as a 3x3 tensor
        '''
        # For this class we force the use of atomic units and ignore internal units
        v_cm1 = self._convert(v)
        v_au = v_cm1 * wavenumber
        eps = self.dielectriContributionsFromModes(v_au, self.mode_list, self.mode_frequencies, 
                                                   self.mode_sigmas, self.mode_oscillator_strengths, 
                                                   self.crystal_volume)
        if self.drude:
            eps = eps + self.dielectricContributionFromDrude(v_au, self.drude_plasma, self.drude_sigma, self.crystal_volume)
        return eps + self.epsilon_infinity

class DrudeLorentz(DielectricFunction):
    '''
    A class for a complex constant tensor dielectric function

    Inherits from DielectricFunction
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, vs_cm1, strengths_cm1, sigmas_cm1,units='cm-1'):
        ''' 
            Parameters
            ----------
            vs_cm1 : list of ints
                    Frequencies of the strengths and sigmas in cm-1
            strengths_cm1 : list of floats shape is either [] or [[xxs][yys][zzs]]
                    strengths of isotropic contribution ([])
                    strengths of orthonormal contributions ([[xxs][yys][zzs]])
            sigmas_cm1 : list of floats shape is either [] or [[xxs][yys][zzs]]
                    sigmas of isotropic contribution ([])
                    sigmas of orthonormal contributions ([[xxs][yys][zzs]])
        '''
        DielectricFunction.__init__(self,units=units)
        self.isScalarFunction = False
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
        return

    def calculate(self, v):
        '''
            Parameters
            ----------
            v : float
                The frequency in the units set by setUnits()

           Returns
           -------
           The permittivity at frequency v as a diagonal 3x3 tensor
        '''
        f_cm1 = self._convert(v)
        eps = np.zeros( (3,3), dtype=complex )
        for xyz, (vs, strengths, sigmas) in enumerate(zip(self.vs_cm1,self.strengths_cm1, self.sigmas_cm1)):
            for v, strength, sigma in zip(vs, strengths,sigmas):
                eps[xyz,xyz] += strength * strength / complex((v*v - f_cm1*f_cm1), -sigma*f_cm1)
        return eps + self.epsilon_infinity

class FPSQ(DielectricFunction):
    '''
    A complex constant tensor FPSQ dielectric function

    Inherits from DielectricFunction
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, omega_tos, gamma_tos, omega_los, gamma_los,units='cm-1'):
        ''' 
            Parameters
            ----------
            omega_tos : list of floats
                    Omega (to) in the FPSQ equation
            gamma_tos : list of floats
                    Gamma (to) in the FPSQ equation
            omega_los : list of floats
                    Omega (lo) in the FPSQ equation
            gamma_los : list of floats
                    Gamma (lo) in the FPSQ equation
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.isScalarFunction = False
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
        return

    def calculate(self, v):
        '''
            Parameters
            ----------
            v : float
                The frequency in units specified by setUnit()

           Returns
           -------
           The diagonal permittivity 3x3 tensor at frequency v
        '''
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

class Sellmeier(DielectricFunction):
    '''
    A complex constant tensor Sellmeier dielectric function

    Inherits from DielectricFunction
    Provides a calculate() function to return the permittivity
    '''
    def __init__(self, Bs, Cs, units='cm-1'):
        ''' 
            Parameters
            ----------
            Bs : list of floats
                    B parameter for Sellmeier equation (microns^2)
            Cs : list of floats
                    C parameter for Sellmeier equation (microns^2)
            units : string 
                    Define the units of frequency (default is cm-1)
        '''
        DielectricFunction.__init__(self,units=units)
        self.isScalarFunction = True
        self.Bs = Bs
        self.Cs = Cs
        return

    def calculate(self, v):
        '''
            Parameters
            ----------
            v : float
                The frequency in units specified by setUnit()

           Returns
           -------
           The diagonal permittivity 3x3 tensor at frequency v

           The Sellmeier parameters are in microns^2
        '''
        f_cm1 = self._convert(v)
        # Convert to wavelength in microns
        wavelength_mu = f_cm1 * 1.0E+6 *1.0E-2            # convert cm-1 to microns
        wavelength2 = wavelength_mu*wavelength_mu
        n2 = 1.0
        for B,C in zip(self.Bs,self.Cs):
            n2 += B*wavelength2/(wavelength2-C) 
        return n2 

