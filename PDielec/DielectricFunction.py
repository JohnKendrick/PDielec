#!/usr/bin/python
#
# Copyright 2024 John Kendrick & Andrew Burnett
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
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""DielectricFunction.

Provides an interface to different mechanisms for providing dielectric information

- Base class :class:`DielectricFunction`

    - Subclass :class:`ConstantTensor`
    - Subclass :class:`ConstantScalar`
    - Subclass :class:`TabulateScalar`
    - Subclass :class:`Tabulate1`
    - Subclass :class:`Tabulate3`
    - Subclass :class:`Tabulate6`
    - Subclass :class:`DFT`
    - Subclass :class:`DrudeLorentz`
    - Subclass :class:`FPSQ`
    - Subclass :class:`Sellmeier`
"""

import sys
from contextlib import nullcontext

import numpy as np
from scipy import interpolate

from PDielec.Constants import angstrom, wavenumber


class DielectricFunction:
    """A base class to different dielectric functions.

    There is a close relationship between the dielectric function and a material.
    The Material module allows data to be read in for different dielectric models.
    Once the data has been read in a DielectricFunction is created which is appropriate for the data.
    The tensor and scalar versions of the functions behave differently with respect to 
    epsilon infinity: tensor functions add epsilon infinity to the calculated permittivity, scalar functions do not.

    Attributes
    ----------
    epsilon_infinity : ndarray
        A 3x3 tensor (numpy array) representing the optical permittivity.
    volume : float
        Volume of crystal unit cells.
    volume_au : float
        Volume of crystal unit cells in atomic units.
    vs_cm1 : ndarray
        An array of floats representing the frequencies of any tabulation.
    isScalarFunction : bool
        True if this function is a scalar function, False otherwise.

    Methods
    -------
    setVolume(volume)
        Change the unit cell volume.
    setEpsilonInfinity(eps)
        Change epsilon infinity.
    function()
        Return the calculate function for this object.
    calculate(v)
        Return the value of the permittivity at frequency v.
    setFrequencies(frequencies)
        Set frequencies 
    getFrequencies()
        Get frequencies 
    setSigmas(sigmas)
        Set sigmas 
    getSigmas()
        Get sigmas 
    setOscillatorStrengths(strengths)
        Set strengths 
    getOscillatorStrengths()
        Get strengths 
    dielectricContributionsFromDrude(f, frequency, sigma, volume)
        Calculates permittivity from a Drude Model.
    dielectricContributionsFromModes(f, modes, frequencies, sigmas, strengths, volume)
        Calculates permittivity from a Drude-Lorentz model.

    """

    possible_epsTypes = ["dft","fpsq","drude-lorentz","sellmeier"]

    def __init__(self ):
        """Initialise and instance of the function."""
        self.isScalarFunction             = None
        self.vs_cm1               = 0
        self.epsilon_infinity     = np.zeros( (3,3) )

    def setFrequencies(self, frequencies):
        """Set the frequencies of a Lorentzian permittivity.

        The default is to print an error message unless the DielectricFunction is a Lorentzian

        Parameters
        ----------
        frequencies : 1d array of floats
            The frequencies in cm-1

        Returns
        -------
        None

        """
        print("Error is DielectricFunction.setFrequencies: unable to set frequencies of a non-Lorentzian function")
        return 


    def getFrequencies(self):
        """Return the frequencies of a Lorentzian permittivity.

        The default is to return None, only return an array if the permittivity is Lorentzian

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        return


    def setOscillatorStrengths(self, strengths):
        """Set the oscillator strengths of a Lorentzian permittivity.

        The default is to print an error message unless the DielectricFunction is a Lorentzian

        Parameters
        ----------
        strengths : 1d array of floats
            The oscillator strengths in cm-1

        Returns
        -------
        None

        """
        print("Error is DielectricFunction.setOscillatorStrengths: unable to set oscillator strengths of a non-Lorentzian function")
        return 


    def getOscillatorStrengths(self):
        """Return the oscillator strengths of a Lorentzian permittivity.

        The default is to return None, only return an array if the permittivity is Lorentzian

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        return

    def setSigmas(self, sigmas):
        """Set the sigma parameters of a Lorentzian permittivity.

        The default is to print an error message unless the DielectricFunction is a Lorentzian

        Parameters
        ----------
        sigmas : 1d array of floats
            The sigma parameters in cm-1

        Returns
        -------
        None

        """
        print("Error is DielectricFunction.setSigmas: unable to set sigmas of a non-Lorentzian function")
        return 


    def getSigmas(self):
        """Return the sigma parameters of a Lorentzian permittivity.

        The default is to return None, only return an array if the permittivity is Lorentzian

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        return


    def isTensor(self):
        """Return true if the dielectric function is a tensor.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the dielectric function returns a tensor, otherwise False.

        """
        return not self.isScalarFunction

    def isScalar(self):
        """Return true if the dielectric function is a scalar.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the dielectric function returns a scalar value, False otherwise.

        """
        return self.isScalarFunction

    def setVolume(self,volume):
        """Set the volume for dielectric calculations.

        Parameters
        ----------
        volume : float
            The volume in angs^3.

        Returns
        -------
        None

        """
        self.volume_angs     = volume
        self.volume_au       = volume*angstrom*angstrom*angstrom
        return

    def setEpsilonInfinity(self,eps):
        """Set epsilon infinity for dielectric calculations.

        Parameters
        ----------
        eps : tensor (3x3)
            The value of eps_infinity.

        Returns
        -------
        None

        """
        self.epsilon_infinity     = np.array(eps)
        return

    def function(self):
        """Return the function used to calculate the permittivity.

        Parameters
        ----------
        None

        Returns
        -------
        a function

        """
        return  self.calculate

    def getLowestFrequency(self):
        """Return the lowest tabulated frequency in cm-1.

        Parameters
        ----------
        None

        Returns
        -------
        The lowest tabulated frequency in cm-1.

        """
        return np.min(self.vs_cm1)

    def getHighestFrequency(self):
        """Return the highest tabulated frequency in cm-1.

        Parameters
        ----------
        None

        Returns
        -------
        float
            The highest tabulated frequency in cm-1.

        """
        return np.max(self.vs_cm1)

    def dielectriContributionsFromDrude(self, f, frequency, sigma, volume):
        """Calculate the dielectric function for a set of a Drude oscillator.

        The units of the parameters must be self-consistent.  Atomic units is a good choice.

        Parameters
        ----------
        f : float
            The frequency of the dielectric response.
        frequency : float
            The frequency of the drude oscillator.
        sigma : float
            The width of the oscillator.
        volume : float
            The volume of the cell

        Returns
        -------
        tensor (3x3)
            A complex scalar dielectric.

        """
        dielectric = np.zeros((3, 3), dtype=complex)
        unit = np.eye(3,dtype=complex)
        # Avoid a divide by zero if f is small
        f = max(1.0e-8, f)
        # Assume that the drude contribution is isotropic
        dielectric = dielectric - unit * frequency*frequency / complex(-f*f, -sigma*f)
        return dielectric * (4.0*np.pi/volume)


    def dielectriContributionsFromModes(self, f, modes, frequencies, sigmas, strengths, volume):
        """Calculate the dielectric function for a set of modes.

        Parameters
        ----------
        f : float
            The frequency considered
        modes : list of int
            A list of the modes to be considered.
        frequencies : list of float
            The frequency of each mode.
        sigmas : list of float
            The width of each mode.
        strengths : list of float
            The strengths of each mode.
        volume : float
            The volume of the cell

        Returns
        -------
        tensor (3x3)
            A complex scalar representing the dielectric function.

        """
        dielectric = np.zeros((3, 3), dtype=complex)
        for mode in modes:
            v = frequencies[mode]
            sigma = sigmas[mode]
            strength = strengths[mode].astype(complex)
            dielectric = dielectric + strength / complex((v*v - f*f), -sigma*f)
        return dielectric * (4.0*np.pi/volume)

    def print(self, v1_cm1, v2_cm1, v_inc=1,diagonal_only=True,file=None):
        """Print the permittivity over a range of frequencies to a file.

        Parameters
        ----------
        v1_cm1 : float
            Starting frequency
        v2_cm1 : float
            Finishing frequency
        v_inc : float, optional
            Incrementing frequency, default is 1
        diagonal_only: bool,optional
            If true only the diagonal values are printed
        file : str, optional
            Filename for printing.  If 'None' then stdout is used

        Returns
        -------
        None

        """
        c = ","
        with open(file, "w") if file else nullcontext(sys.stdout) as fd:
            if diagonal_only:
                print("f,eps00.real,epss11.real,eps22.real,eps00.imag,eps11.imag,eps22,imag",file=fd)
            else:
                print("f,eps00,epss11,eps22,eps01,eps02,,eps12",file=fd)
            for v in np.arange(v1_cm1, v2_cm1, v_inc):
                epsilon = self.calculate(v)
                if diagonal_only:
                    print(v,c,epsilon[0,0].real,c,epsilon[1,1].real,c,epsilon[2,2].real,c,epsilon[0,0].imag,c,epsilon[1,1].imag,c,epsilon[2,2].imag,file=fd)
                else:
                    print(v,c,epsilon[0,0],c,epsilon[1,1],c,epsilon[2,2],c,epsilon[0,1],c,epsilon[0,2],c,epsilon[1,2],file=fd)
        return

class ConstantTensor(DielectricFunction):
    """A complex constant tensor permittivity.

    Inherits from DielectricFunction.
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is used to calculate the total permittivity

    Methods
    -------
    calculate()
        Calculates and returns the permittivity.

    Examples
    --------
    .. code-block:: python

        permittivity = np.array([ [ 2.0, 0.0, 0.0 ] , [ 0.0, 3.0, 0.0 ], [ 0.0, 0.0, 4.0 ] ])
        df = ConstantTensor(permittivity)
        newMaterial = Material('New material name',permittivityObject=df,cell=UnitCell(2, 2, 2))

    """

    def __init__(self, value):
        """Initialise ConstantTensor.

        Parameters
        ----------
        value : 3x3 numpy array
            Initialise the permittivity with this tensor

        """
        DielectricFunction.__init__(self)
        self.value = value
        self.isScalarFunction = False

    def calculate(self,v):
        """Calculate the permittivity at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in internal units.

        Returns
        -------
        ndarray
            The permittivity at frequency v as a 3x3 tensor.

        """
        return self.value * np.eye(3) + self.epsilon_infinity

class ConstantScalar(DielectricFunction):
    """A complex constant scalar permittivity.

    Inherits from DielectricFunction. Provides a calculate() function to return
    the permittivity.  Epsilon infinity is not used for this function.

    Methods
    -------
    calculate()
        Calculates and returns the permittivity of the dielectric function.

    """

    def __init__(self, value):
        """Intialise ConstantScalar dielectric function.

        Parameters
        ----------
        value : float
            Initialise the permittivity with this tensor

        """
        DielectricFunction.__init__(self)
        self.value = value
        self.isScalarFunction = True

    def calculate(self,v):
        """Calculate the permittivity at a specified frequency.

        Parameters
        ----------
        v : float
            The frequency at which to calculate the permittivity.

        Returns
        -------
        float
            The permittivity at frequency v.

        """
        return self.value

class TabulateScalar(DielectricFunction):
    """A complex tabulated scalar permittivity.

    Inherits from DielectricFunction. 
    Represents an isotropic material dielectric function. 
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is not used in the calculation.

    Methods
    -------
    calculate()
        Compute the permittivity of the material.

    Notes
    -----
    None

    Examples
    --------
    None

    """

    def __init__(self, vs_cm1, permittivities):
        """Intialise TabulateScalar dielectric function.

        Parameters
        ----------
        vs_cm1 : list of float
            Frequencies in cm-1 of the following permittivities.
        permittivities : list of complex
            The permittivities at each of the frequencies.

        """
        DielectricFunction.__init__(self)
        self.isScalarFunction = True
        self.vs_cm1 = np.array(vs_cm1)
        eps = np.array(permittivities)
        eps_r = np.real(eps)
        eps_i = np.imag(eps)
        self.interpolater = interpolate.CubicSpline(self.vs_cm1,eps_r)
        self.interpolatei = interpolate.CubicSpline(self.vs_cm1,eps_i)

    def calculate(self,v):
        """Calculate the permittivity at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in internal units.

        Returns
        -------
        float
            The permittivity at frequency v as a scalar.

        """
        return complex(self.interpolater(v),self.interpolatei(v))

class Tabulate1(DielectricFunction):
    """A complex tabulated isotropic tensor permittivity.

    Inherits from DielectricFunction.
    Represents an isotropic dielectric function.
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is included in the total permittivity

    Methods
    -------
    calculate()
        Returns the permittivity of the dielectric function.

    """

    def __init__(self, vs_cm1, permittivities):
        """Intialise Tabulate1 dielectric function.

        Parameters
        ----------
        vs_cm1 : list of float
            Frequencies in cm-1 of the following permittivities.
        permittivities : list of complex
            The permittivities at each of the frequencies.

        """
        DielectricFunction.__init__(self)
        self.isScalarFunction = False
        self.vs_cm1 = np.array(vs_cm1)
        eps = np.array(permittivities)
        eps_r = np.real(eps)
        eps_i = np.imag(eps)
        self.interpolater = interpolate.CubicSpline(self.vs_cm1,eps_r)
        self.interpolatei = interpolate.CubicSpline(self.vs_cm1,eps_i)

    def calculate(self,v):
        """Calculate the isotropic permittivity at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in internal units.

        Returns
        -------
        numpy.ndarray
            The isotropic permittivity at frequency v as a 3x3 tensor.

        """
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater(v),self.interpolatei(v))
        eps[1,1] = eps[0,0]
        eps[2,2] = eps[0,0]
        return eps + self.epsilon_infinity

class Tabulate3(DielectricFunction):
    """A complex tabulated diagonal tensor permittivity.

    Inherits from DielectricFunction.
    Represents an orthorhombic dielectric function.
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is included in the total permittivity

    Methods
    -------
    calculate()
        Return the permittivity of the dielectric function.

    """

    def __init__(self, vs_cm1, epsxx, epsyy, epszz):
        """Intialise a Tabulate3 dielectric function object.

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

        """
        DielectricFunction.__init__(self)
        self.isScalarFunction = False
        self.vs_cm1 = np.array(vs_cm1)
        epsd = [np.array(epsxx), np.array(epsyy), np.array(epszz)]
        self.interpolater = []
        self.interpolatei = []
        for eps in epsd:
            self.interpolater.append(interpolate.CubicSpline(self.vs_cm1,np.real(eps)))
            self.interpolatei.append(interpolate.CubicSpline(self.vs_cm1,np.imag(eps)))

    def calculate(self,v):
        """Calculate the permittivity at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in internal units.

        Returns
        -------
        numpy.ndarray
            The permittivity at frequency v as a diagonal 3x3 tensor.

        """
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater[0](v),self.interpolatei[0](v))
        eps[1,1] = complex(self.interpolater[1](v),self.interpolatei[1](v))
        eps[2,2] = complex(self.interpolater[2](v),self.interpolatei[2](v))
        return eps + self.epsilon_infinity

class Tabulate6(DielectricFunction):
    """A complex tabulated, anisotropic tensor permittivity.

    Inherits from DielectricFunction.
    Represents an non-isotropic dielectric function.
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is included in the total permittivity

    Methods
    -------
    calculate()
        Returns the permittivity of the dielectric.

    """

    def __init__(self, vs_cm1, epsxx, epsyy, epszz, epsxy, epsxz, epsyz):
        """Intialise a Tabulate6 dielectric function object.

        Parameters
        ----------
        vs_cm1 : list of float
            The frequencies of the tabulated dielectric.
        epsxx : list of complex
            The eps(xx) permittivities at each of the frequencies.
        epsyy : list of complex
            The eps(yy) permittivities at each of the frequencies.
        epszz : list of complex
            The eps(zz) permittivities at each of the frequencies.
        epsxy : list of complex
            The eps(xy) permittivities at each of the frequencies.
        epsxz : list of complex
            The eps(xz) permittivities at each of the frequencies.
        epsyz : list of complex
            The eps(yz) permittivities at each of the frequencies.

        """
        DielectricFunction.__init__(self)
        self.isScalarFunction = False
        self.vs_cm1 = np.array(vs_cm1)
        eps6 = [np.array(epsxx), np.array(epsyy), np.array(epszz), np.array(epsxy), np.array(epsxz), np.array(epsyz)]
        self.interpolater = []
        self.interpolatei = []
        for eps in eps6:
            epsr = np.real(eps)
            epsi = np.imag(eps)
            self.interpolater.append(interpolate.CubicSpline(self.vs_cm1,epsr))
            self.interpolatei.append(interpolate.CubicSpline(self.vs_cm1,epsi))

    def calculate(self,v):
        """Calculate the symmetric permittivity tensor at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in cm-1.

        Returns
        -------
        numpy.ndarray
            The full permittivity tensor at frequency v as a 3x3 tensor.

        """
        eps = np.zeros( (3,3), dtype=complex )
        eps[0,0] = complex(self.interpolater[0](v),self.interpolatei[0](v))
        eps[1,1] = complex(self.interpolater[1](v),self.interpolatei[1](v))
        eps[2,2] = complex(self.interpolater[2](v),self.interpolatei[2](v))
        eps[1,0] = complex(self.interpolater[3](v),self.interpolatei[3](v))
        eps[2,0] = complex(self.interpolater[4](v),self.interpolatei[4](v))
        eps[2,1] = complex(self.interpolater[5](v),self.interpolatei[5](v))
        eps[0,1] = eps[1,0]
        eps[0,2] = eps[2,0]
        eps[1,2] = eps[2,1]
        return eps + self.epsilon_infinity

class DFT(DielectricFunction):
    """A complex tensor permittivity from a DFT calculation.

    Inherits from DielectricFunction.
    Suitable for use in calculating the permittivity from DFT calculations.
    Internal units for this function are atomic units.
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is included in the total permittivity

    Methods
    -------
    calculate()
        Returns the permittivity based on the DFT calculations.

    Notes
    -----
    The internal units for all calculations within this class are atomic units.

    """

    def __init__(self,mode_list, mode_frequencies, mode_sigmas, mode_oscillator_strengths, crystal_volume, drude, drude_plasma, drude_sigma):
        """Intialise a DFT dielectric function.

        Parameters
        ----------
        mode_list : list of int
            A list of active modes.
        mode_frequencies : list of float
            All mode frequencies.
        mode_sigmas : list of float
            All mode frequencies.
        mode_oscillator_strengths : list of float
            All mode oscillator strengths.
        crystal_volume : float
            The crystal volume.
        drude : bool
            True if a drude oscillator is to be included.
        drude_plasma : float
            Drude plasma frequency.
        drude_sigma : float
            Drude plasma width.
        All parameters are given in atomic units.

        """
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
        """Calculate the permittivity at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in cm-1.

        Returns
        -------
        ndarray
            The permittivity at frequency v as a 3x3 tensor.

        """
        # For this class we force the use of atomic units and ignore internal units
        v_au = v * wavenumber
        eps = self.dielectriContributionsFromModes(v_au, self.mode_list, self.mode_frequencies, 
                                                   self.mode_sigmas, self.mode_oscillator_strengths, 
                                                   self.crystal_volume)
        if self.drude:
            eps = eps + self.dielectricContributionFromDrude(v_au, self.drude_plasma, self.drude_sigma, self.crystal_volume)
        return eps + self.epsilon_infinity

    def setFrequencies(self,frequencies):
        """Set the frequencies of a DFT permittivity.

        The frequencies are converted to atomic units

        Parameters
        ----------
        frequencies : 1d array of floats
            The frequencies in cm-1

        Returns
        -------
        None

        """
        frequencies = np.array(frequencies)
        self.mode_frequencies = wavenumber * frequencies
        return 


    def getFrequencies(self):
        """Return the frequencies of a DFT permittivity.

        The frequencies are in units of cm-1

        Parameters
        ----------
        None

        Returns
        -------
        1d array of floats
            The frequencies in cm-1

        """
        return self.mode_frequencies / wavenumber

    def setOscillatorStrengths(self,strengths):
        """Set the oscillator strengths of a DFT permittivity.

        The oscillator strengths are converted to atomic units

        Parameters
        ----------
        strengths : 3x3 array for each frequency
            The oscillator strengths parameters in cm-1

        Returns
        -------
        None

        """
        strengths = np.array(strengths)
        self.mode_oscillator_strengths = wavenumber * strengths
        return 


    def getOscillatorStrengths(self):
        """Return the oscillator strengths of a DFT permittivity.

        The returned strengths are in units of cm-1

        Parameters
        ----------
        None

        Returns
        -------
        3x3 array of floats for each frequency
            The oscillator strengths in cm-1

        """
        return self.mode_oscillator_strengths / wavenumber

    def setSigmas(self,sigmas):
        """Set the sigma parameters of a DFT permittivity.

        The sigmas parameter is converted to atomic units

        Parameters
        ----------
        sigmas : 1d array of floats
            The sigma parameters in cm-1

        Returns
        -------
        None

        """
        sigmas = np.array(sigmas)
        self.mode_sigmas = wavenumber * sigmas
        return 


    def getSigmas(self):
        """Return the sigma parameters of a DFT permittivity.

        The returned sigmas are in units of cm-1

        Parameters
        ----------
        None

        Returns
        -------
        1d array of floats
            The sigma parameters of  Drude-Lorentzian permittivity

        """
        return self.mode_sigmas / wavenumber


class DrudeLorentz(DielectricFunction):
    """A complex tensor Drude-Lorentz permittivity.

    Inherits from DielectricFunction.
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is included in the total permittivity

    Methods
    -------
    calculate()
        Calculates and returns the permittivity.

    Examples
    --------
    None

    """

    def __init__(self, vs_cm1, strengths_cm1, sigmas_cm1):
        """Initialise a DrudeLorenz dielectric function.

        Parameters
        ----------
        vs_cm1 : list of ints
            Frequencies of the strengths and sigmas in cm-1.
        strengths_cm1 : list of floats
            strengths of isotropic contribution ([]).
            strengths of orthonormal contributions ([[xxs][yys][zzs]]), shape is either [] or [[xxs][yys][zzs]].
        sigmas_cm1 : list of floats
            Sigmas of isotropic contribution ([]).
            Sigmas of orthonormal contributions ([[xxs][yys][zzs]]), shape is either [] or [[xxs][yys][zzs]].

        """
        DielectricFunction.__init__(self)
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

    def calculate(self, f_cm1):
        """Calculate the permittivity at a given frequency.

        Parameters
        ----------
        f_cm1 : float
            The frequency in cm-1

        Returns
        -------
        tensor
            The permittivity at frequency v as a diagonal 3x3 tensor.

        """
        eps = np.zeros( (3,3), dtype=complex )
        for xyz, (vs, strengths, sigmas) in enumerate(zip(self.vs_cm1,self.strengths_cm1, self.sigmas_cm1)):
            for v, strength, sigma in zip(vs, strengths,sigmas):
                eps[xyz,xyz] += strength * strength / complex((v*v - f_cm1*f_cm1), -sigma*f_cm1)
        return eps + self.epsilon_infinity

    def setFrequencies(self,frequencies):
        """Set the frequencies of a Drude-Lorentz permittivity.

        The frequencies are in cm-1

        Parameters
        ----------
        frequencies : 1d array of floats
            The frequencies in cm-1

        Returns
        -------
        None

        """
        frequencies = np.array(frequencies)
        self.vs_cm1 = frequencies
        return 


    def getFrequencies(self):
        """Return the frequencies of a Drude-Lorentz permittivity.

        The frequencies are in units of cm-1

        Parameters
        ----------
        None

        Returns
        -------
        1d array of floats
            The frequencies in cm-1

        """
        return self.vs_cm1

    def setOscillatorStrengths(self,strengths):
        """Set the oscillator strengths of a Drude-Lorentz permittivity.

        The oscillator strengths are in cm-1

        Parameters
        ----------
        strengths : 1d array of floats
            The oscillator strengths parameters in cm-1

        Returns
        -------
        None

        """
        strengths = np.array(strengths)
        self.strengths = strengths
        return 

    def getOscillatorStrengths(self):
        """Return the oscillator strengths of a Drude-Lorentz permittivity.

        The returned strengths are in units of cm-1

        Parameters
        ----------
        None

        Returns
        -------
        1d array of floats
            The oscillator strengths in cm-1

        """
        return self.strengths

    def setSigmas(self,sigmas):
        """Set the sigma parameters of a Drude-Lorentz permittivity.

        The default is to print an error message unless the DielectricFunction is a Lorentzian

        Parameters
        ----------
        sigmas : 1d array of floats
            The sigma parameters in cm-1

        Returns
        -------
        None

        """
        self.sigmas_cm1 = sigmas
        return 


    def getSigmas(self):
        """Return the sigma parameters of a Lorentzian permittivity.

        The returned sigmas are in units of cm-1

        Parameters
        ----------
        None

        Returns
        -------
        1d array of floats

        """
        return self.sigmas_cm1


class FPSQ(DielectricFunction):
    """A complex tensor FPSQ dielectric function.

    Inherits from DielectricFunction.
    Provides a calculate() function to return the permittivity.
    Epsilon infinity is included in the total permittivity
    """

    def __init__(self, omega_tos, gamma_tos, omega_los, gamma_los):
        """Initialise an FPS dielectric function object.

        Parameters
        ----------
        omega_tos : list of floats
            Omega (to) in the FPSQ equation.
        gamma_tos : list of floats
            Gamma (to) in the FPSQ equation.
        omega_los : list of floats
            Omega (lo) in the FPSQ equation.
        gamma_los : list of floats
            Gamma (lo) in the FPSQ equation.

        """
        DielectricFunction.__init__(self)
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
        """Calculate the diagonal permittivity 3x3 tensor at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in cm-1

        Returns
        -------
        numpy.ndarray
            The diagonal permittivity 3x3 tensor at frequency v.

        """
        f_cm1 = v
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
    """A real tensor Sellmeier permittivity.

    Inherits from DielectricFunction. Provides a calculate() function to return the permittivity.
    Epsilon infinity is not included in the total permittivity

    """

    def __init__(self, Bs, Cs):
        """Initialise an Sellmeier dielectric function object.

        Parameters
        ----------
        Bs : list of floats
            B parameter for Sellmeier equation (microns^2)
        Cs : list of floats
            C parameter for Sellmeier equation (microns^2)

        """
        DielectricFunction.__init__(self)
        self.isScalarFunction = True
        self.Bs = Bs
        self.Cs = Cs
        return

    def calculate(self, v):
        """Calculate the diagonal permittivity 3x3 tensor at a given frequency.

        Parameters
        ----------
        v : float
            The frequency in cm-1.

        Returns
        -------
        ndarray
            The diagonal permittivity 3x3 tensor at frequency v.

        Notes
        -----
        The Sellmeier parameters are in microns^2.

        """
        f_cm1 = v
        # Convert to wavelength in microns
        wavelength_mu = f_cm1 * 1.0E+6 *1.0E-2            # convert cm-1 to microns
        wavelength2 = wavelength_mu*wavelength_mu
        n2 = 1.0
        for B,C in zip(self.Bs,self.Cs):
            n2 += B*wavelength2/(wavelength2-C) 
        return n2 

