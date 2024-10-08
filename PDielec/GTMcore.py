"""The pyGTM module.

It has been heavily modified by John Kendrick for inclusion in the PDielec library
Quite a lot has been removed as the fields spatial distribution was not needed
Thanks to the authors of the original code for all their hard work

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) Mathieu Jeannin 2019 2020 <math.jeannin@free.fr>.

This module implements the generalized 4x4 transfer matrix (GTM) method
poposed in `Passler, N. C. and Paarmann, A., JOSA B 34, 2128 (2017)
<http://doi.org/10.1364/JOSAB.34.002128>`_
and corrected in
`JOSA B 36, 3246 (2019) <http://doi.org/10.1364/JOSAB.36.003246>`_,
as well as the layer-resolved absorption proposed in
`Passler, Jeannin and Paarman <https://arxiv.org/abs/2002.03832>`_.
This code uses inputs from D. Dietze's FSRStools library
https://github.com/ddietze/FSRStools

Please cite the relevant associated publications if you use this code.

Author:
    - Mathieu Jeannin math.jeannin@free.fr (permanent)

Affiliations:
    - Laboratoire de Physique de l'Ecole Normale Superieure (2019)
    - Centre de Nanosciences et Nanotechnologies (2020-2021)

Layers are represented by the :py:class:`Layer` class that holds all parameters
describing the optical properties of a single layer.
The optical system is assembled using the :py:class:`System` class.
  
           
**Change log:**

    *01-10-2024*:
        - Moving back to 63 bit arithmetic as scattering formalism is stable

    *01-01-2023*: 
        
        - Major changes introduced by John Kendrick to make the package compatible with PDielec and PDGUI.
        - Allowed only full 3x3 tensors for the permittivity
        - Changed from 64 to 128 bit arithmetic where possible
        - Introduced a scattering matrix formalism
        - Added treatment of incoherence
    
    *15-10-2021*:
        
        - Fixed rounding error bug in lag.eig() causing the program to crash randomly for negligibly small imaginary parts of the wavevectors
        - Corrected a sign error in gamma32 that lead to field discontinuities

    *23 June 2021*:
        - integrated the code into pdielec and pdgui

    *19-03-2020*:

        - Adapted the code to compute the layer-resolved absorption as proposed
          by Passler et al. (https://arxiv.org/abs/2002.03832), using
          :py:func:`System.calculate_Poynting_Absorption_vs_z`.

        - Include the correct calculation of intensity transmission coefficients
          in :py:func:`System.calculate_r_t`.
          **This BREAKS compatibility** with the previous definition of the function.

        - Corrected bugs in :py:func:`System.calculate_Efield`
          and added magnetic field option

        - Adapted :py:func:`System.calculate_Efield` to allow hand-defined,
          irregular grid and a shorthand to compute only at layers interfaces.
          Regular grid with fixed resolution is left as an option.

    *20-09-2019*:
        - Added functions in the :py:class:`System` class to compute in-plane
          wavevector of guided modes and dispersion relation for such guided surface modes.
          This is *highly prospective* as it depends on the robustness of the minimization
          procedure (or the lack of thereoff)
"""
######## general utilities

import copy
import sys

import numpy as np
import scipy.linalg as lag

from PDielec.Constants import epsilon_0_si, speed_light_si

#jk c_const = 299792458 # m/s
#jk eps0 = 8.854e-12 ## vacuum permittivity
#
# JK using PDielec constants
#
c_const = speed_light_si
eps0 = epsilon_0_si            # Doesn't appear to be used

def vacuum_eps(f):
    """Vacuum permittivity function.

    Parameters
    ----------
    f : float or 1D-array
        Frequency (in Hz)

    Returns
    -------
    eps : complex or 1D-array of complex
        Complex value of the vacuum permittivity (1.0 + 0.0j)

    """
    try:
        return np.ones(len(f))
    except Exception:
        return 1.0 + 0.0j
    
def exact_inv_2x2(M):
    """Calculate the inverse of 2x2 complex matrix, M.

    Parameters
    ----------
    M : complex 2x2 array
        The matrix of which to find the inverse.

    Returns
    -------
    ndarray
        The inverse of M as a 2x2 complex numpy array.

    """
    assert M.shape == (2, 2)
    m11, m21, m12, m22 = M.flatten()
    determinant = m11*m22 - m21*m12
    return np.array([[ m22, -m21],
                     [-m12,  m11]])/determinant

    
def exact_inv_3x3(M):
    """Calculate the inverse of a 3x3 complex matrix M.

    Parameters
    ----------
    M : complex 3x3 array
        The matrix to be inverted.

    Returns
    -------
    numpy.ndarray
        The inverse of M as a 3x3 complex numpy array.

    """
    assert M.shape == (3, 3)
    m1, m2, m3, m4, m5, m6, m7, m8, m9 = M.flatten()
    determinant = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9  
    return np.array([[m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5],
                     [m6*m7-m4*m9, m1*m9-m3*m7, m3*m4-m1*m6],
                     [m4*m8-m5*m7, m2*m7-m1*m8, m1*m5-m2*m4]])/determinant

    
def exact_inv_4x4(M):
    """Compute the 'exact' inverse of a 4x4 matrix using the analytical result.

    Parameters
    ----------
    M : array-like, shape (4, 4)
        Matrix to be inverted, consisting of float or complex numbers.

    Returns
    -------
    array-like, shape (4, 4), complex
        Inverse of this matrix or Moore-Penrose approximation if matrix cannot be inverted.

    Notes
    -----
    This should give a higher precision and speed at a reduced noise, following D.Dietze's implementation in FSRStools. 

    See Also
    --------
    http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html : For mathematical background related to this computation.

    """
    assert M.shape == (4, 4)

    # the following equations use algebraic indexing; transpose input matrix to get indexing right
    A = M.T
    detA = A[0, 0] * A[1, 1] * A[2, 2] * A[3, 3] + A[0, 0] * A[1, 2] * A[2, 3] * A[3, 1] + A[0, 0] * A[1, 3] * A[2, 1] * A[3, 2]
    detA = detA + A[0, 1] * A[1, 0] * A[2, 3] * A[3, 2] + A[0, 1] * A[1, 2] * A[2, 0] * A[3, 3] + A[0, 1] * A[1, 3] * A[2, 2] * A[3, 0]
    detA = detA + A[0, 2] * A[1, 0] * A[2, 1] * A[3, 3] + A[0, 2] * A[1, 1] * A[2, 3] * A[3, 0] + A[0, 2] * A[1, 3] * A[2, 0] * A[3, 1]
    detA = detA + A[0, 3] * A[1, 0] * A[2, 2] * A[3, 1] + A[0, 3] * A[1, 1] * A[2, 0] * A[3, 2] + A[0, 3] * A[1, 2] * A[2, 1] * A[3, 0]

    detA = detA - A[0, 0] * A[1, 1] * A[2, 3] * A[3, 2] - A[0, 0] * A[1, 2] * A[2, 1] * A[3, 3] - A[0, 0] * A[1, 3] * A[2, 2] * A[3, 1]
    detA = detA - A[0, 1] * A[1, 0] * A[2, 2] * A[3, 3] - A[0, 1] * A[1, 2] * A[2, 3] * A[3, 0] - A[0, 1] * A[1, 3] * A[2, 0] * A[3, 2]
    detA = detA - A[0, 2] * A[1, 0] * A[2, 3] * A[3, 1] - A[0, 2] * A[1, 1] * A[2, 0] * A[3, 3] - A[0, 2] * A[1, 3] * A[2, 1] * A[3, 0]
    detA = detA - A[0, 3] * A[1, 0] * A[2, 1] * A[3, 2] - A[0, 3] * A[1, 1] * A[2, 2] * A[3, 0] - A[0, 3] * A[1, 2] * A[2, 0] * A[3, 1]


    if detA == 0:
        try:
            print("Warning 4x4 inversion problem 1")
            result = np.cdouble(lag.pinv(np.cdouble(M)))
        except Exception:
            print("Warning 4x4 inversion problem 2")
            result = np.cdouble(lag.pinv(np.cdouble(M)))
        return result

    B = np.zeros(A.shape, dtype=np.cdouble)
    B[0, 0] = A[1, 1] * A[2, 2] * A[3, 3] + A[1, 2] * A[2, 3] * A[3, 1] + A[1, 3] * A[2, 1] * A[3, 2] - A[1, 1] * A[2, 3] * A[3, 2] - A[1, 2] * A[2, 1] * A[3, 3] - A[1, 3] * A[2, 2] * A[3, 1]
    B[0, 1] = A[0, 1] * A[2, 3] * A[3, 2] + A[0, 2] * A[2, 1] * A[3, 3] + A[0, 3] * A[2, 2] * A[3, 1] - A[0, 1] * A[2, 2] * A[3, 3] - A[0, 2] * A[2, 3] * A[3, 1] - A[0, 3] * A[2, 1] * A[3, 2]
    B[0, 2] = A[0, 1] * A[1, 2] * A[3, 3] + A[0, 2] * A[1, 3] * A[3, 1] + A[0, 3] * A[1, 1] * A[3, 2] - A[0, 1] * A[1, 3] * A[3, 2] - A[0, 2] * A[1, 1] * A[3, 3] - A[0, 3] * A[1, 2] * A[3, 1]
    B[0, 3] = A[0, 1] * A[1, 3] * A[2, 2] + A[0, 2] * A[1, 1] * A[2, 3] + A[0, 3] * A[1, 2] * A[2, 1] - A[0, 1] * A[1, 2] * A[2, 3] - A[0, 2] * A[1, 3] * A[2, 1] - A[0, 3] * A[1, 1] * A[2, 2]

    B[1, 0] = A[1, 0] * A[2, 3] * A[3, 2] + A[1, 2] * A[2, 0] * A[3, 3] + A[1, 3] * A[2, 2] * A[3, 0] - A[1, 0] * A[2, 2] * A[3, 3] - A[1, 2] * A[2, 3] * A[3, 0] - A[1, 3] * A[2, 0] * A[3, 2]
    B[1, 1] = A[0, 0] * A[2, 2] * A[3, 3] + A[0, 2] * A[2, 3] * A[3, 0] + A[0, 3] * A[2, 0] * A[3, 2] - A[0, 0] * A[2, 3] * A[3, 2] - A[0, 2] * A[2, 0] * A[3, 3] - A[0, 3] * A[2, 2] * A[3, 0]
    B[1, 2] = A[0, 0] * A[1, 3] * A[3, 2] + A[0, 2] * A[1, 0] * A[3, 3] + A[0, 3] * A[1, 2] * A[3, 0] - A[0, 0] * A[1, 2] * A[3, 3] - A[0, 2] * A[1, 3] * A[3, 0] - A[0, 3] * A[1, 0] * A[3, 2]
    B[1, 3] = A[0, 0] * A[1, 2] * A[2, 3] + A[0, 2] * A[1, 3] * A[2, 0] + A[0, 3] * A[1, 0] * A[2, 2] - A[0, 0] * A[1, 3] * A[2, 2] - A[0, 2] * A[1, 0] * A[2, 3] - A[0, 3] * A[1, 2] * A[2, 0]

    B[2, 0] = A[1, 0] * A[2, 1] * A[3, 3] + A[1, 1] * A[2, 3] * A[3, 0] + A[1, 3] * A[2, 0] * A[3, 1] - A[1, 0] * A[2, 3] * A[3, 1] - A[1, 1] * A[2, 0] * A[3, 3] - A[1, 3] * A[2, 1] * A[3, 0]
    B[2, 1] = A[0, 0] * A[2, 3] * A[3, 1] + A[0, 1] * A[2, 0] * A[3, 3] + A[0, 3] * A[2, 1] * A[3, 0] - A[0, 0] * A[2, 1] * A[3, 3] - A[0, 1] * A[2, 3] * A[3, 0] - A[0, 3] * A[2, 0] * A[3, 1]
    B[2, 2] = A[0, 0] * A[1, 1] * A[3, 3] + A[0, 1] * A[1, 3] * A[3, 0] + A[0, 3] * A[1, 0] * A[3, 1] - A[0, 0] * A[1, 3] * A[3, 1] - A[0, 1] * A[1, 0] * A[3, 3] - A[0, 3] * A[1, 1] * A[3, 0]
    B[2, 3] = A[0, 0] * A[1, 3] * A[2, 1] + A[0, 1] * A[1, 0] * A[2, 3] + A[0, 3] * A[1, 1] * A[2, 0] - A[0, 0] * A[1, 1] * A[2, 3] - A[0, 1] * A[1, 3] * A[2, 0] - A[0, 3] * A[1, 0] * A[2, 1]

    B[3, 0] = A[1, 0] * A[2, 2] * A[3, 1] + A[1, 1] * A[2, 0] * A[3, 2] + A[1, 2] * A[2, 1] * A[3, 0] - A[1, 0] * A[2, 1] * A[3, 2] - A[1, 1] * A[2, 2] * A[3, 0] - A[1, 2] * A[2, 0] * A[3, 1]
    B[3, 1] = A[0, 0] * A[2, 1] * A[3, 2] + A[0, 1] * A[2, 2] * A[3, 0] + A[0, 2] * A[2, 0] * A[3, 1] - A[0, 0] * A[2, 2] * A[3, 1] - A[0, 1] * A[2, 0] * A[3, 2] - A[0, 2] * A[2, 1] * A[3, 0]
    B[3, 2] = A[0, 0] * A[1, 2] * A[3, 1] + A[0, 1] * A[1, 0] * A[3, 2] + A[0, 2] * A[1, 1] * A[3, 0] - A[0, 0] * A[1, 1] * A[3, 2] - A[0, 1] * A[1, 2] * A[3, 0] - A[0, 2] * A[1, 0] * A[3, 1]
    B[3, 3] = A[0, 0] * A[1, 1] * A[2, 2] + A[0, 1] * A[1, 2] * A[2, 0] + A[0, 2] * A[1, 0] * A[2, 1] - A[0, 0] * A[1, 2] * A[2, 1] - A[0, 1] * A[1, 0] * A[2, 2] - A[0, 2] * A[1, 1] * A[2, 0]

    return B.T / detA


###########################
#### The Layer Class   ####
###########################
class Layer:
    """Layer class. An instance is a single layer.
    
    The inherited layer classes wich are used in PDielec are:

        - :class:`CoherentLayer`
        - :class:`SemiInfiniteLayer`
        - :class:`IncoherentIntensityLayer`
        - :class:`IncoherentAveragePhaseLayer`
        - :class:`IncoherentPhaseLayer`
        - :class:`IncoherentThickLayer`

    Attributes
    ----------
    thickness : float
        Thickness of the layer in meters (m).
    epsilon : complex function
        Function `epsilon(frequency)` representing the full dielectric constant tensor.
    theta : float
        Euler angle theta (colatitude) in radians (rad).
    phi : float
        Euler angle phi in radians (rad).
    psi : float
        Euler angle psi in radians (rad).

    Notes
    -----
    If instantiated with default values, it generates a 1 µm thick layer of air.
    Properties can be checked/changed dynamically using the corresponding get/set methods.

    """

    # Define class variables
    qsd_thr = 1e-10                  # threshold for wavevector comparison
    zero_thr = 1e-10                 # threshold for eigenvalue comparison to zero
    jk_shift = sys.float_info.min    # The smallest number to be used to shift a division by zero

    def __init__(self, thickness=1.0e-6, epsilon=None, theta=0, phi=0, psi=0, exponent_threshold=700):

        ## epsilon is a 3x3 matrix of permittivity at a given frequency
        """Initialize the class instance with default parameters or specified properties.

        Parameters
        ----------
        thickness : float, optional
            Thickness of the material layer, specified in meters. Default is 1e-06.
        epsilon : np.ndarray or None, optional
            The dielectric tensor of the material. If `None`, an identity tensor is used. Default is None.
        theta : float, optional
            The Euler angle theta (rotation around z-axis) in radians. Default is 0.
        phi : float, optional
            The Euler angle phi (rotation around x-axis) in radians. Default is 0.
        psi : float, optional
            The Euler angle psi (rotation around z-axis again) in radians. Defaults to 0.
        exponent_threshold : float, optional
            The threshold for the exponent under which the calculations are considered safe from overflow errors. Default is 700.

        Attributes
        ----------
        epsilon : np.ndarray
            Dielectric tensor of the material.
        mu : float
            The magnetic permeability of the material. Initialized to 1.0.
        exponent_errors : int
            Counts the number of exponent errors encountered during calculations.
        largest_exponent : float
            Tracks the largest exponent encountered that is still within a safe range.
        exponent_threshold : float
            Threshold for safe exponent calculations.
        coherent : bool
            Flag indicating if coherent calculations should be performed.
        M : np.ndarray
            Placeholder matrix for calculations.
        a : np.ndarray
            Placeholder matrix for calculations.
        Delta : np.ndarray
            Placeholder matrix for calculations.
        qs : np.ndarray
            Placeholder vector for calculations.
        Py : np.ndarray
            Placeholder matrix for calculations.
        gamma : np.ndarray
            Placeholder matrix for calculations.
        Ai : np.ndarray
            Placeholder matrix for calculations.
        Ki : np.ndarray
            Placeholder matrix for calculations.
        Ti : np.ndarray
            Placeholder matrix for calculations.
        Berreman : np.ndarray
            Placeholder matrix used for Berreman's calculations.
        useBerreman : bool
            Flag indicating whether Berreman's method is used.
        euler : np.ndarray
            Euler rotation matrix used for orienting the material's tensor properties.
        epsilon_tensor_function : callable or None
            Function to compute the dielectric tensor if provided.
        propagation_exponents : np.ndarray
            Placeholder array for propagation exponent calculations.

        Raises
        ------
        ValueError
            When the input `epsilon` is not compatible in dimensions or data type.

        """        
        self.epsilon = np.identity(3, dtype=np.cdouble)
        self.mu = 1.0 ### mu=1 for now

        ## Added to handle overflow 
        self.exponent_errors = 0
        self.largest_exponent = 0.0
        self.exponent_threshold = exponent_threshold
        self.coherent = True

        ### initialization of all important quantities
        self.M = np.zeros((6, 6), dtype=np.cdouble)                  # constitutive relations
        self.a = np.zeros((6, 6), dtype=np.cdouble) 
        self.Delta = np.zeros((4, 4), dtype=np.cdouble) 
        self.qs = np.zeros(4, dtype=np.cdouble)                      # out of plane wavevector
        self.Py = np.zeros((3,4), dtype=np.cdouble)                  # Poynting vector
        self.gamma = np.zeros((4, 3), dtype=np.cdouble)
        self.Ai = np.zeros((4, 4), dtype=np.cdouble) 
        self.Ki = np.zeros((4, 4), dtype=np.cdouble) 
        self.Ti = np.zeros((4, 4), dtype=np.cdouble)                 # Layer transfer matrix
        self.Berreman = np.zeros((4,3), dtype=np.cdouble)            # Stores the Berreman modes
        self.useBerreman = False          # Boolean to replace Xu's eigenvectors by Berreman's in case of Birefringence
        self.euler = np.identity(3, dtype=np.cdouble)                # rotation matrix
        self.epsilon_tensor_function = None                              # Added by JK to store epsilon
        self.propagation_exponents = np.zeros((4), dtype=np.cdouble) # Stores the raw propagation exponents
        
        # Store initialiser variables
        self.set_thickness(thickness)          # set the thickness, 1um by default
        self.set_epsilon(epsilon)              # set epsilon, vacuum by default (JK change)
        self.set_euler(theta, phi, psi)        # set orientation of crystal axis w/ respect to the lab frame

    def isCoherent(self):
        """Return True if the layer is a coherent layer.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the layer is coherent, False if not.

        """
        return self.coherent


    def set_thickness(self, thickness):
        """Set the layer thickness.

        Parameters
        ----------
        thickness : float
            The layer thickness (in m).

        Returns
        -------
        None

        """
        self.thick = thickness

    def set_epsilon(self, epsilon_function):
        """Set the dielectric function.

        Parameters
        ----------
        epsilon_function : complex function
            A function representing the complex dielectric function.

        Returns
        -------
        None

        """
        self.epsilon_tensor_function = epsilon_function

    def calculate_epsilon(self, f_in):
        """Set the value of epsilon in the (rotated) lab frame.

        Parameters
        ----------
        f_in : float
            frequency (in Hz)

        Returns
        -------
        None

        Notes
        -----
        The values are set according to the epsilon_fi (i=1..3) functions
        defined using the :py:func:`set_epsilon` method, at the given frequency f.
        The rotation with respect to the lab frame is computed using the Euler angles.

        Use explicitly if you *don't* use the :py:func:`Layer.update` function!
        Modification by JK to allow the use of a full dielectric tensor.

        """
        f = f_in / (c_const * 100.0)
        epsilon_xstal = self.epsilon_tensor_function(f)
        if not isinstance(epsilon_xstal,np.ndarray):
            epsilon_xstal = np.eye(3,dtype=np.cdouble)*epsilon_xstal
        #JK Changed the order of the transformation to a transformation of basis vectors
        # self.epsilon = np.matmul(self.euler_inverse, np.matmul(epsilon_xstal,self.euler))
        # rather than a transformation of the coordinates active->passive
        self.epsilon = np.matmul(self.euler, np.matmul(epsilon_xstal,self.euler_inverse))
        return self.epsilon


    def set_euler(self,theta,phi,psi):
        """Set the values for the Euler rotations angles.

        Parameters
        ----------
        theta : float
            Euler angle theta (colatitude) in rad.
        phi : float
            Euler angle phi in rad.
        psi : float
            Euler angle psi in rad.

        Returns
        -------
        None

        """
        self.theta = theta
        self.phi = phi
        self.psi = psi
        # euler matrix for rotation of dielectric tensor
        self.euler[0, 0] = np.cos(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.sin(psi)
        self.euler[0, 1] = -np.sin(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.cos(psi)
        self.euler[0, 2] = np.sin(theta) * np.sin(phi)
        self.euler[1, 0] = np.cos(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.sin(psi)
        self.euler[1, 1] = -np.sin(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.cos(psi)
        self.euler[1, 2] = -np.sin(theta) * np.cos(phi)
        self.euler[2, 0] = np.sin(theta) * np.sin(psi)
        self.euler[2, 1] = np.sin(theta) * np.cos(psi)
        self.euler[2, 2] = np.cos(theta)
        #JK Added the inverse calculation here so it is only done once
        #JK self.euler_inverse = np.cdouble(lag.pinv(np.cdouble(self.euler)))
        self.euler_inverse = exact_inv_3x3(self.euler)

    def calculate_matrices(self, zeta):
        """Calculate the principal matrices necessary for the GTM algorithm.

        Parameters
        ----------
        zeta : complex
            In-plane reduced wavevector kx/k0 in the system.

        Returns
        -------
        None

        Notes
        -----
        Note that zeta is conserved through the whole system and set externally
        using the angle of incidence and `System.superstrate.epsilon[0,0]` value.

        Requires prior execution of :py:func:`calculate_epsilon`.

        """
        ## Constitutive matrix (see e.g. eqn (4))
        self.M[0:3, 0:3] = np.cdouble(self.epsilon.copy())
        self.M[3:6, 3:6] = np.cdouble(self.mu*np.identity(3))

        ## from eqn (10)
        b = self.M[2,2]*self.M[5,5] - self.M[2,5]*self.M[5,2]

        ## a matrix from eqn (9)
        self.a[2,0] = (self.M[5,0]*self.M[2,5] - self.M[2,0]*self.M[5,5])/b
        self.a[2,1] = ((self.M[5,1]-zeta)*self.M[2,5] - self.M[2,1]*self.M[5,5])/b
        self.a[2,3] = (self.M[5,3]*self.M[2,5] - self.M[2,3]*self.M[5,5])/b
        self.a[2,4] = (self.M[5,4]*self.M[2,5] - (self.M[2,4]+zeta)*self.M[5,5])/b
        self.a[5,0] = (self.M[5,2]*self.M[2,0] - self.M[2,2]*self.M[5,0])/b
        self.a[5,1] = (self.M[5,2]*self.M[2,1] - self.M[2,2]*(self.M[5,1]-zeta))/b
        self.a[5,3] = (self.M[5,2]*self.M[2,3] - self.M[2,2]*self.M[5,3])/b
        self.a[5,4] = (self.M[5,2]*(self.M[2,4]+zeta) - self.M[2,2]*self.M[5,4])/b

        ## S Matrix (Don't know where it comes from since Delta is just S re-ordered)
        ## Note that after this only Delta is used
        S = np.zeros((4, 4), dtype=np.cdouble) ##
        S[0,0] = self.M[0,0] + self.M[0,2]*self.a[2,0] + self.M[0,5]*self.a[5,0]
        S[0,1] = self.M[0,1] + self.M[0,2]*self.a[2,1] + self.M[0,5]*self.a[5,1]
        S[0,2] = self.M[0,3] + self.M[0,2]*self.a[2,3] + self.M[0,5]*self.a[5,3]
        S[0,3] = self.M[0,4] + self.M[0,2]*self.a[2,4] + self.M[0,5]*self.a[5,4]
        S[1,0] = self.M[1,0] + self.M[1,2]*self.a[2,0] + (self.M[1,5]-zeta)*self.a[5,0]
        S[1,1] = self.M[1,1] + self.M[1,2]*self.a[2,1] + (self.M[1,5]-zeta)*self.a[5,1]
        S[1,2] = self.M[1,3] + self.M[1,2]*self.a[2,3] + (self.M[1,5]-zeta)*self.a[5,3]
        S[1,3] = self.M[1,4] + self.M[1,2]*self.a[2,4] + (self.M[1,5]-zeta)*self.a[5,4]
        S[2,0] = self.M[3,0] + self.M[3,2]*self.a[2,0] + self.M[3,5]*self.a[5,0]
        S[2,1] = self.M[3,1] + self.M[3,2]*self.a[2,1] + self.M[3,5]*self.a[5,1]
        S[2,2] = self.M[3,3] + self.M[3,2]*self.a[2,3] + self.M[3,5]*self.a[5,3]
        S[2,3] = self.M[3,4] + self.M[3,2]*self.a[2,4] + self.M[3,5]*self.a[5,4]
        S[3,0] = self.M[4,0] + (self.M[4,2]+zeta)*self.a[2,0] + self.M[4,5]*self.a[5,0]
        S[3,1] = self.M[4,1] + (self.M[4,2]+zeta)*self.a[2,1] + self.M[4,5]*self.a[5,1]
        S[3,2] = self.M[4,3] + (self.M[4,2]+zeta)*self.a[2,3] + self.M[4,5]*self.a[5,3]
        S[3,3] = self.M[4,4] + (self.M[4,2]+zeta)*self.a[2,4] + self.M[4,5]*self.a[5,4]


        ## Delta Matrix from eqn (8)
        self.Delta[0,0] = S[3,0]
        self.Delta[0,1] = S[3,3]
        self.Delta[0,2] = S[3,1]
        self.Delta[0,3] = - S[3,2]
        self.Delta[1,0] = S[0,0]
        self.Delta[1,1] = S[0,3]
        self.Delta[1,2] = S[0,1]
        self.Delta[1,3] = - S[0,2]
        self.Delta[2,0] = -S[2,0]
        self.Delta[2,1] = -S[2,3]
        self.Delta[2,2] = -S[2,1]
        self.Delta[2,3] = S[2,2]
        self.Delta[3,0] = S[1,0]
        self.Delta[3,1] = S[1,3]
        self.Delta[3,2] = S[1,1]
        self.Delta[3,3] = -S[1,2]

    def calculate_q(self):
        """Calculate the 4 out-of-plane wavevectors for the current layer.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        From this we also get the Poynting vectors. Wavevectors are sorted according
        to (trans-p, trans-s, refl-p, refl-s). Birefringence is determined according
        to a threshold value `qsd_thr` set at the beginning of the script.

        """
        Delta_loc = np.zeros((4,4), dtype=np.cdouble)
        transmode = np.zeros((2), dtype=int)
        reflmode = np.zeros((2), dtype=int)
        
        Delta_loc = np.cdouble(self.Delta.copy())
        ## eigenvals // eigenvects as of eqn (11)
        qsunsorted, psiunsorted = lag.eig(Delta_loc)
        ##### remove extremely small real/imaginary parts that are due to numerical inaccuracy
        for km in range(4):
            if (np.abs(np.imag(qsunsorted[km])) > 0) and (np.abs(np.imag(qsunsorted[km])) < Layer.zero_thr):
                qsunsorted[km] = np.real(qsunsorted[km]) + 0.0j
            if (np.abs(np.real(qsunsorted[km])) > 0) and (np.abs(np.real(qsunsorted[km])) < Layer.zero_thr):
                qsunsorted[km] = 0.0 + 1.0j*np.imag(qsunsorted[km])
        for comp in range(4):
            if (np.abs(np.real(psiunsorted[km][comp]))>0) and (np.abs(np.real(psiunsorted[km][comp])) < Layer.zero_thr):
                psiunsorted[km][comp] = 0.0 + 1.0j*np.imag(psiunsorted[km][comp])
            if (np.abs(np.imag(psiunsorted[km][comp]))>0) and (np.abs(np.imag(psiunsorted[km][comp])) < Layer.zero_thr):
                psiunsorted[km][comp] = np.real(psiunsorted[km][comp]) + 0.0j

                
        Berreman_unsorted = np.zeros((4,3), dtype=np.cdouble)
        
        kt = 0 
        kr = 0
        ## sort berremann qi's according to (12)
        ##JK small modification to fix problems when qs are nearly real
        ##JK if any(np.abs(np.imag(qsunsorted))):
        if any(np.abs(np.imag(qsunsorted))>1.0E-16):
            for km in range(0,4):
                if np.imag(qsunsorted[km])>=0 :
                    transmode[kt] = km
                    kt = kt + 1
                else:
                    reflmode[kr] = km
                    kr = kr +1
        else:
            for km in range(0,4):
                ##JK if np.real(qsunsorted[km])>0:
                if np.real(qsunsorted[km])>=0 and kt < 2  :
                    transmode[kt] = km
                    kt = kt + 1
                else:
                    reflmode[kr] = km
                    kr = kr +1
        ## Calculate the Poynting vector for each Psi using (16-18)
        for km in range(0,4):
            Ex = psiunsorted[0,km]
            Ey = psiunsorted[2,km]
            Hx = -psiunsorted[3,km]
            Hy = psiunsorted[1,km]
            ## from eqn (17)
            Ez = self.a[2,0]*Ex + self.a[2,1]*Ey + self.a[2,3]*Hx + self.a[2,4]*Hy
            # from eqn (18)
            Hz = self.a[5,0]*Ex + self.a[5,1]*Ey + self.a[5,3]*Hx + self.a[5,4]*Hy
            ## and from (16)
            self.Py[0,km] = Ey*Hz-Ez*Hy
            self.Py[1,km] = Ez*Hx-Ex*Hz
            self.Py[2,km] = Ex*Hy-Ey*Hx
            ### Berreman modes (unsorted) in case they are needed later (birefringence)
            Berreman_unsorted[km,0] = Ex
            Berreman_unsorted[km,1] = Ey
            Berreman_unsorted[km,2] = Ez
        ## check Cp using either the Poynting vector for birefringent
        ## materials or the electric field vector for non-birefringent
        ## media to sort the modes       
        
        ## first calculate Cp for transmitted waves
        Cp_t1 = np.abs(self.Py[0,transmode[0]])**2/(np.abs(self.Py[0,transmode[0]])**2+np.abs(self.Py[1,transmode[0]])**2+Layer.jk_shift)
        Cp_t2 = np.abs(self.Py[0,transmode[1]])**2/(np.abs(self.Py[0,transmode[1]])**2+np.abs(self.Py[1,transmode[1]])**2+Layer.jk_shift)

        if np.abs(Cp_t1-Cp_t2) > Layer.qsd_thr: ## birefringence
            self.useBerreman = True ## sets _useBerreman fo the calculation of gamma matrix below
            if Cp_t2>Cp_t1:
                transmode = np.flip(transmode,0) ## flip the two values
            ## then calculate for reflected waves if necessary
            Cp_r1 = np.abs(self.Py[0,reflmode[1]])**2/(np.abs(self.Py[0,reflmode[1]])**2+np.abs(self.Py[1,reflmode[1]])**2+Layer.jk_shift)
            Cp_r2 = np.abs(self.Py[0,reflmode[0]])**2/(np.abs(self.Py[0,reflmode[0]])**2+np.abs(self.Py[1,reflmode[0]])**2+Layer.jk_shift)
            if Cp_r1>Cp_r2:
                reflmode = np.flip(reflmode,0) ## flip the two values

        else:     ### No birefringence, use the Electric field s-pol/p-pol
            Cp_te1 = np.abs(psiunsorted[0,transmode[1]])**2/(np.abs(psiunsorted[0,transmode[1]])**2+np.abs(psiunsorted[2,transmode[1]])**2+Layer.jk_shift)
            Cp_te2 = np.abs(psiunsorted[0,transmode[0]])**2/(np.abs(psiunsorted[0,transmode[0]])**2+np.abs(psiunsorted[2,transmode[0]])**2+Layer.jk_shift)
            if Cp_te1>Cp_te2:
                transmode = np.flip(transmode,0) ## flip the two values
            Cp_re1 = np.abs(psiunsorted[0,reflmode[1]])**2/(np.abs(psiunsorted[0,reflmode[1]])**2+np.abs(psiunsorted[2,reflmode[1]])**2+Layer.jk_shift)
            Cp_re2 = np.abs(psiunsorted[0,reflmode[0]])**2/(np.abs(psiunsorted[0,reflmode[0]])**2+np.abs(psiunsorted[2,reflmode[0]])**2+Layer.jk_shift)
            if Cp_re1>Cp_re2:
                reflmode = np.flip(reflmode,0) ## flip the two values

        ## finaly store the sorted version
        ####### q is (trans-p, trans-s, refl-p, refl-s)
        self.qs[0] = qsunsorted[transmode[0]]
        self.qs[1] = qsunsorted[transmode[1]]
        self.qs[2] = qsunsorted[reflmode[0]]
        self.qs[3] = qsunsorted[reflmode[1]]
        Py_temp = self.Py.copy()
        self.Py[:,0] = Py_temp[:,transmode[0]]
        self.Py[:,1] = Py_temp[:,transmode[1]]
        self.Py[:,2] = Py_temp[:,reflmode[0]]
        self.Py[:,3] = Py_temp[:,reflmode[1]]
        ### Store the (sorted) Berreman modes
        self.Berreman[0] = Berreman_unsorted[transmode[0],:]
        self.Berreman[1] = Berreman_unsorted[transmode[1],:]
        self.Berreman[2] = Berreman_unsorted[reflmode[0],:]
        self.Berreman[3] = Berreman_unsorted[reflmode[1],:]
        return
        
    def calculate_gamma(self, zeta):
        """Calculate the gamma matrix.

        Parameters
        ----------
        zeta : complex
            In-plane reduced wavevector kx/k0.

        Returns
        -------
        None

        """
        ### this whole function is eqn (20)
        self.gamma[0,0] = 1.0 + 0.0j
        self.gamma[1,1] = 1.0 + 0.0j
        self.gamma[3,1] = 1.0 + 0.0j
        self.gamma[2,0] = -1.0 + 0.0j
        
        ### convenience definition of the repetitive factor
        mu_eps33_zeta2 = (self.mu*self.epsilon[2,2]-zeta**2)

        if np.abs(self.qs[0]-self.qs[1])<Layer.qsd_thr:
            gamma12 = 0.0 + 0.0j
            
            gamma13 = -(self.mu*self.epsilon[2,0]+zeta*self.qs[0])/mu_eps33_zeta2
            
            gamma21 = 0.0 + 0.0j
            
            gamma23 = -self.mu*self.epsilon[2,1]/mu_eps33_zeta2
        
        else:
            gamma12_num = self.mu*self.epsilon[1,2]*(self.mu*self.epsilon[2,0]+zeta*self.qs[0])
            gamma12_num = gamma12_num - self.mu*self.epsilon[1,0]*mu_eps33_zeta2
            gamma12_denom = mu_eps33_zeta2*(self.mu*self.epsilon[1,1]-zeta**2-self.qs[0]**2)
            gamma12_denom = gamma12_denom - self.mu**2*self.epsilon[1,2]*self.epsilon[2,1]
            gamma12 = gamma12_num/gamma12_denom
            if np.isnan(gamma12):
                gamma12 = 0.0 + 0.0j
            
            gamma13 = -(self.mu*self.epsilon[2,0]+zeta*self.qs[0])
            gamma13 = gamma13-self.mu*self.epsilon[2,1]*gamma12 
            gamma13 = gamma13/mu_eps33_zeta2
            
            if np.isnan(gamma13):
                gamma13 = -(self.mu*self.epsilon[2,0]+zeta*self.qs[0])
                gamma13 = gamma13/mu_eps33_zeta2

            gamma21_num = self.mu*self.epsilon[2,1]*(self.mu*self.epsilon[0,2]+zeta*self.qs[1])
            gamma21_num = gamma21_num-self.mu*self.epsilon[0,1]*mu_eps33_zeta2
            gamma21_denom = mu_eps33_zeta2*(self.mu*self.epsilon[0,0]-self.qs[1]**2)
            gamma21_denom = gamma21_denom-(self.mu*self.epsilon[0,2]+zeta*self.qs[1])*(self.mu*self.epsilon[2,0]+zeta*self.qs[1])
            gamma21 = gamma21_num/gamma21_denom
            if np.isnan(gamma21):
                gamma21 = 0.0+0.0j
                
            gamma23 = -(self.mu*self.epsilon[2,0] +zeta*self.qs[1])*gamma21-self.mu*self.epsilon[2,1]
            gamma23 = gamma23/mu_eps33_zeta2
            if np.isnan(gamma23):
                gamma23 = -self.mu*self.epsilon[2,1]/mu_eps33_zeta2

        if np.abs(self.qs[2]-self.qs[3])<Layer.qsd_thr:
            gamma32 = 0.0 + 0.0j
            gamma33 = (self.mu*self.epsilon[2,0]+zeta*self.qs[2])/mu_eps33_zeta2
            gamma41 = 0.0 + 0.0j
            gamma43 = -self.mu*self.epsilon[2,1]/mu_eps33_zeta2
        
        else:
            gamma32_num = self.mu*self.epsilon[1,0]*mu_eps33_zeta2
            gamma32_num = gamma32_num-self.mu*self.epsilon[1,2]*(self.mu*self.epsilon[2,0]+zeta*self.qs[2])
            gamma32_denom = mu_eps33_zeta2*(self.mu*self.epsilon[1,1]-zeta**2-self.qs[2]**2)
            gamma32_denom = gamma32_denom-self.mu**2*self.epsilon[1,2]*self.epsilon[2,1]
            gamma32 = gamma32_num/gamma32_denom
            if np.isnan(gamma32):
                gamma32 = 0.0 + 0.0j
            
            gamma33 = self.mu*self.epsilon[2,0] + zeta*self.qs[2]
            gamma33 = gamma33 + self.mu*self.epsilon[2,1]*gamma32 
            gamma33 = gamma33/mu_eps33_zeta2
            if np.isnan(gamma33):
                gamma33 = (self.mu*self.epsilon[2,0] + zeta*self.qs[2])/mu_eps33_zeta2

            gamma41_num = self.mu*self.epsilon[2,1]*(self.mu*self.epsilon[0,2]+zeta*self.qs[3])
            gamma41_num = gamma41_num - self.mu*self.epsilon[0,1]*mu_eps33_zeta2
            gamma41_denom = mu_eps33_zeta2*(self.mu*self.epsilon[0,0]-self.qs[3]**2)
            gamma41_denom = gamma41_denom - (self.mu*self.epsilon[0,2]+zeta*self.qs[3])*(self.mu*self.epsilon[2,0]+zeta*self.qs[3])
            gamma41 = gamma41_num/gamma41_denom
            if np.isnan(gamma41):
                gamma41 = 0.0 + 0.0j
                
            gamma43 = -(self.mu*self.epsilon[2,0]+zeta*self.qs[3])*gamma41
            gamma43 = gamma43-self.mu*self.epsilon[2,1]
            gamma43 = gamma43/mu_eps33_zeta2
            if np.isnan(gamma43):
                gamma43 = -self.mu*self.epsilon[2,1]/mu_eps33_zeta2
        
        ### gamma field vectors should be normalized to avoid any birefringence problems
        # use double square bracket notation to ensure correct array shape
        gamma1 = np.array([[self.gamma[0,0], gamma12, gamma13]],dtype=np.cdouble)
        gamma2 = np.array([[gamma21, self.gamma[1,1], gamma23]],dtype=np.cdouble)
        gamma3 = np.array([[self.gamma[2,0], gamma32, gamma33]],dtype=np.cdouble)
        gamma4 = np.array([[gamma41, self.gamma[3,1], gamma43]],dtype=np.cdouble)
        
        #### Regular case, no birefringence, we keep the Xu fields
        self.gamma[0,:] = gamma1/np.sqrt(np.sum(gamma1**2))      #lag.norm(gamma1)
        self.gamma[1,:] = gamma2/np.sqrt(np.sum(gamma2**2))      #lag.norm(gamma2)
        self.gamma[2,:] = gamma3/np.sqrt(np.sum(gamma3**2))      #lag.norm(gamma3)
        self.gamma[3,:] = gamma4/np.sqrt(np.sum(gamma4**2))      #lag.norm(gamma4)

        #### In case of birefringence, use Berreman fields
        for ki in range(4): 
            ### normalize them first
            #self.Berreman[ki] = self.Berreman[ki]/lag.norm(self.Berreman[ki])
            self.Berreman[ki] = self.Berreman[ki]/np.sqrt(np.sum(self.Berreman[ki]**2))
        if self.useBerreman:
            self.gamma = self.Berreman
        return

    def calculate_transfer_matrix(self, f, zeta):
        """Compute the transfer matrix of the whole layer :math:`T_i=A_iP_iA_i^{-1}`.

        Parameters
        ----------
        f : float
            Frequency (in Hz).
        zeta : complex
            Reduced in-plane wavevector kx/k0.

        Returns
        -------
        None

        """
        self.calculate_ai(zeta)
        Aim1 = exact_inv_4x4(self.Ai.copy())
        self.calculate_propagation_matrix(f)
        self.Ti = np.matmul(self.Ai,np.matmul(self.Ki,Aim1))
        return

    def calculate_ai(self,zeta):
        """Calculate A_i.

        Boundary matrix :math:`A_i` of the layer.

        Parameters
        ----------
        zeta : float
            Reduced in-plane wavevector kx/k0.

        """
        ## eqn(22)
        self.Ai[0,:] = self.gamma[:,0].copy()
        self.Ai[1,:] = self.gamma[:,1].copy()
        self.Ai[2,:] = (self.qs*self.gamma[:,0]-zeta*self.gamma[:,2])/self.mu
        self.Ai[3,:] = self.qs*self.gamma[:,1]/self.mu
        return 
   

    def update_tm(self, f, zeta):
        """Shortcut to recalculate all layer properties.

        Appropriate for a transfer matrix method

        Parameters
        ----------
        f : float
            frequency (in Hz)
        zeta : complex
            reduced in-plane wavevector kx/k0

        Returns
        -------
        Ai : 4x4-array
            Boundary matrix :math:`A_i` of the layer
        Ki : 4x4-array
            Propagation matrix :math:`K_i` of the layer
        Ai_inv : 4x4-array
            Inverse of the :math:`A_i` matrix
        Ti : 4x4-array
            Transfer matrix of the whole layer

        """
        self.calculate_epsilon(f)
        self.calculate_matrices(zeta)
        self.calculate_q()
        self.calculate_gamma(zeta)
        self.calculate_transfer_matrix(f, zeta)
        Ai_inv = exact_inv_4x4(self.Ai.copy())
        return[self.Ai.copy(), self.Ki.copy(), Ai_inv.copy(), self.Ti.copy()]

    def update_sm(self, f, zeta):
        """Shortcut to recalculate all layer properties for scattering method.

        Appropriate for a scattering matrix method.
        This avoids the calculation of the exponential; it just calculates the exponents.

        Parameters
        ----------
        f : float
            frequency (in Hz)
        zeta : complex
               reduced in-plane wavevector kx/k0

        Calculates
        ----------
        Ai : array_like, shape (4, 4)
            Boundary matrix :math:`A_i` of the layer.
        Ki : array_like, shape (4, 4)
            Propagation matrix :math:`K_i` of the layer.
        Ai_inv : array_like, shape (4, 4)
            Inverse of the :math:`A_i` matrix.
        Ti : array_like, shape (4, 4)
            Transfer matrix of the whole layer.

        """
        self.calculate_epsilon(f)
        self.calculate_matrices(zeta)
        self.calculate_q()
        self.calculate_gamma(zeta)
        self.calculate_ai(zeta)
        self.calculate_propagation_exponents(f)
        return

    def calculate_propagation_exponents(self,f):
        """Calculate the propagation exponents.

        Parameters
        ----------
        f : float
            The frequency

        """
        for ii in range(4):
            self.propagation_exponents[ii] = np.cdouble(-1.0j*(2.0*np.pi*f*self.qs[ii]*self.thick)/c_const)
        return

    def calculate_propagation_matrix(self,f):
        """Calculate the matrix Ki (or Pi depending on the paper).

        Parameters
        ----------
        f : float
            Frequency

        Returns
        -------
        Ki : ndarray
            A 4x4-array representing the boundary matrix :math:`K_i` of the layer.

        """
        Ki = np.zeros( (4,4), dtype=np.cdouble)
        self.calculate_propagation_exponents(f)
        for ii in range(4):
            exponent = self.propagation_exponents[ii]
            if np.abs(exponent) > self.exponent_threshold:
                self.largest_exponent = max(self.largest_exponent,np.abs(exponent))
                exponent = exponent/np.abs(exponent)*self.exponent_threshold
                self.exponent_errors += 1
            Ki[ii,ii] = np.exp(exponent)
        self.Ki = Ki
        return  self.Ki
        

###############################
### The Layer Sub Classes   ###
###############################
class CoherentLayer(Layer):
    """Define a coherent layer inherits from Layer class."""

    def __init__(self, layer, theta=0, phi=0, psi=0, exponent_threshold=700):
        """Initialise an instance of a coherent layer.

        Parameters
        ----------
        layer : Layer object
            Use this to get the thickness
        theta : float, optional
            Theta value for layer (default is 0)
        phi : float, optional
            Phi value for layer (default is 0)
        psi : float, optional
            Psi value for layer (default is 0)
        exponent_threshold : float, optional
            Prevent calculation of large exponents (default is 700)

        """
        thickness = layer.getThicknessInMetres()
        epsilon = layer.getPermittivityFunction()
        Layer.__init__(self, thickness, epsilon, theta, phi, psi, exponent_threshold)
        self.coherent = True
        self.inCoherentIntensity = False
        self.inCoherentPhase = False
        self.inCoherentAveragePhase = False
        self.inCoherentThick = False
        self.SMatrix = None            # A holder for the scattering matrix for the layer
        return

    def calculate_scattering_matrix(self,b):
        """Calculate the scattering matrix of this layer with layer b.

        This method calculates the scattering matrix based on equations 18-21 in the PyLama paper.
        However, it's important to note that the definition of a scattering matrix used here 
        differs from that in the paper.

        Parameters
        ----------
        b : a pyGTM layer class
            The second layer for which the scattering matrix will be calculated

        Returns
        -------
        None

        Notes
        -----
        Explain any important details about the implementation, external references like the PyLama paper,
        and any differences in definitions or approaches.

        """
        # This is layer 'a'
        a = self
        # The forward an back propagators are defined here
        # Note that the exponent is calculated from the raw exponent, no
        # test is being made to avoid overflow
        forward_exponents =      a.propagation_exponents
        reverse_exponents = -1.0*a.propagation_exponents
        K_f = np.eye(4,dtype=np.cdouble)
        K_f[0,0] = np.exp(reverse_exponents[0])
        K_f[1,1] = np.exp(reverse_exponents[1])
        # Define the backward propagators
        K_b = np.eye(4,dtype=np.cdouble)
        K_b[2,2] = np.exp(forward_exponents[2])
        K_b[3,3] = np.exp(forward_exponents[3])
        # Change rows and columns so rows are the eigenvectors
        aa = a.Ai.T
        ab = b.Ai.T
        #if self.inCoherentIntensity:
        #    K_f = np.sign(K_f)*(np.absolute(K_f)**2)
        #    K_b = np.sign(K_b)*(np.absolute(K_b)**2)
        #    aa = np.sign(aa)*(np.absolute(aa)**2)
        #    ab = np.sign(ab)*(np.absolute(ab)**2)
        # build the in and out A matrices row by row
        # transform back so the columns are the eigenvectors
        A_out = np.array([ aa[0], aa[1], -ab[2], -ab[3] ]).T
        A_in  = np.array([ ab[0], ab[1], -aa[2], -aa[3] ]).T
        # Calculate the scattering matrix
        # Avoid doing 2 inversions on K_b by using the forward exponents
        # Taken from the paper: S = exact_inv_4x4(K_b).dot(exact_inv_4x4(A_in)).dot(A_out).dot(K_f)
        S = K_b.dot(exact_inv_4x4(A_in)).dot(A_out).dot(K_f)
        #if self.inCoherentIntensity:
        #    S = np.sqrt(S)
        # So far S-matrix is defined here as follows
        #  [ E_{i+1}^f ]   = S_{i,i+1} [ E_{i}^f   ]
        #  [ E_{i}^b   ]               [ E_{i+1}^b ]
        # But the following definition is much more usual
        #  [ E_{i}^b   ]   = S_{i,i+1} [ E_{i}^f   ]
        #  [ E_{i+1}^f ]               [ E_{i+1}^b ]
        brows = S[0:2,:].copy()
        arows = S[2:4,:].copy()
        S[0:2,:] = arows
        S[2:4,:] = brows
        # Finally store the S-matrix and return
        self.SMatrix = SMatrix(S)
        return

class SemiInfiniteLayer(CoherentLayer):
    """Define a semi-infinite layer, inherits from CoherentLayer class."""

    def __init__(self, layer, theta=0, phi=0, psi=0, exponent_threshold=700):
        """Initialise an instance of a semi-infinite layer.

        Parameters
        ----------
        layer : Layer object
            Use this to get the thickness
        theta : float, optional
            Theta value for layer (default is 0)
        phi : float, optional
            Phi value for layer (default is 0)
        psi : float, optional
            Psi value for layer (default is 0)
        exponent_threshold : float, optional
            Prevent calculation of large exponents (default is 700)

        """
        thickness = layer.getThicknessInMetres()
        epsilon = layer.getPermittivityFunction()
        Layer.__init__(self, thickness, epsilon, theta, phi, psi, exponent_threshold)
        self.coherent = True
        self.inCoherentIntensity = False
        self.inCoherentPhase = False
        self.inCoherentAveragePhase = False
        self.inCoherentThick = False
        self.SMatrix = None            # A holder for the scattering matrix
        return

    def calculate_propagation_exponents(self,f):
        """Calculate the propagation exponents.

        For the semi-infinite case the propagation exponents are set to 0

        Parameters
        ----------
        f : float
            The frequency (ignored for this object)

        Returns
        -------
        None

        """
        self.propagation_exponents = np.zeros( (4), dtype=np.cdouble)
        return

    def calculate_propagation_matrix(self,f):
        """Calculate the matrix Ki (or Pi) depending on the paper.

        This routine makes the material transparent. It is therefore suitable for 
        a semi-infinite layer.

        Parameters
        ----------
        f : float
            Frequency

        Returns
        -------
        Ki : tensor (3x3)

        """
        self.calculate_propagation_exponents(f)
        return np.eye( 4, dtype=np.cdouble)
        
class IncoherentIntensityLayer(CoherentLayer):
    """Define an incoherent layer using intensity transfer matrices.

    Inherits from the CoherentLayer class.
    """

    def __init__(self, layer, theta=0, phi=0, psi=0, exponent_threshold=700):
        """Initialise an instance of an incoherent layer based on intensities.

        Parameters
        ----------
        layer : Layer object
            Use this to get the thickness
        theta : float, optional
            Theta value for layer (default is 0)
        phi : float, optional
            Phi value for layer (default is 0)
        psi : float, optional
            Psi value for layer (default is 0)
        exponent_threshold : float, optional
            Prevent calculation of large exponents (default is 700)

        """
        thickness = layer.getThicknessInMetres()
        epsilon = layer.getPermittivityFunction()
        Layer.__init__(self, thickness, epsilon, theta, phi, psi, exponent_threshold)
        self.coherent = False
        self.inCoherentIntensity = True
        self.inCoherentPhase = False
        self.inCoherentAveragePhase = False
        self.inCoherentThick = False
        return

class IncoherentAveragePhaseLayer(CoherentLayer):
    """Define an incoherent layer using an average phase in the propagation matrix.

    Inherits from CoherentLayer
    """

    def __init__(self, layer, percentage_incoherence=100, number_of_samples=4, theta=0, phi=0, psi=0, exponent_threshold=700):
        """Initialise an instance of an incoherent layer, which uses phase averaging.

        Parameters
        ----------
        layer : Layer object
            Use this to get the thickness
        theta : float, optional
            Theta value for layer (default is 0)
        phi : float, optional
            Phi value for layer (default is 0)
        psi : float, optional
            Psi value for layer (default is 0)
        exponent_threshold : float, optional
            Prevent calculation of large exponents (default is 700)
        percentage_incoherence : float
            The percentage of incoherence to be incorporated into the calculation
        number_of_samples : int, optional
            The number of samples to be used in averages (default is 4)

        """
        thickness = layer.getThicknessInMetres()
        epsilon = layer.getPermittivityFunction()
        Layer.__init__(self, thickness, epsilon, theta, phi, psi, exponent_threshold)
        self.coherent = False
        self.inCoherentIntensity = False
        self.inCoherentPhase = False
        self.inCoherentAveragePhase = True
        self.inCoherentThick = False
        self.phaseShift = layer.getPhaseShift()
        return

    def calculate_propagation_exponents(self,f):
        """Calculate the matrix Ki (or Pi depending on the paper).

        A phase shift is included in the calculation.
        The phase shift can vary between 0 and :math:`\\pi`.

        Parameters
        ----------
        f : float
            Frequency

        Returns
        -------
        None

        """
        for ii in range(4):
            if ii < 2:
                exponent = np.cdouble(-1.0j*(2.0*np.pi*f*self.qs[ii]*self.thick/c_const))
            else:
                exponent = np.cdouble(-1.0j*(2.0*np.pi*f*self.qs[ii]*self.thick/c_const + self.phaseShift))
            self.propagation_exponents[ii] = exponent
        return

class IncoherentPhaseLayer(CoherentLayer):
    """Define an incoherent layer using Arteaga's modification of the phase in the propagation matrix.

    Inherits from CoherentLayer
    """

    def __init__(self, layer, theta=0, phi=0, psi=0, exponent_threshold=700):
        """Initialise an instance of an incoherent layer, which uses Arteaga's method.

        Parameters
        ----------
        layer : Layer object
            Use this to get the thickness
        theta : float, optional
            Theta value for layer (default is 0)
        phi : float, optional
            Phi value for layer (default is 0)
        psi : float, optional
            Psi value for layer (default is 0)
        exponent_threshold : float, optional
            Prevent calculation of large exponents (default is 700)

        """
        thickness = layer.getThicknessInMetres()
        epsilon = layer.getPermittivityFunction()
        Layer.__init__(self, thickness, epsilon, theta, phi, psi, exponent_threshold)
        self.coherent = False
        self.inCoherentIntensity = False
        self.inCoherentPhase = True
        self.inCoherentAveragePhase = False
        self.inCoherentThick = False
        return

    def calculate_propagation_exponents(self,f):
        """Calculate the matrix Ki (or Pi) depending on the paper.

        This routine is taken from Arteaga et al., Thin Solid Films 2014, 571, 701-705.
        The propagation matrix is suitable for incoherent films. In Passler, the ordering is:
        Arteaga:
        0 p/o -> p ->
        1 s/e -> p <-
        2 p/o <- s ->
        3 s/e <- s <-

        Parameters
        ----------
        f : float
            Frequency

        Returns
        -------
        None

        """
        qs = np.zeros( 4, dtype=np.cdouble)
        # p-wave ->
        qs[0] = self.qs[0] - np.real(self.qs[1])    # Subtract the real part of s ->
        # s-wave ->
        qs[1] = self.qs[1] - np.real(self.qs[1])    # Subtract the real part of s ->
        # p-wave <-
        qs[2] = self.qs[2] - np.real(self.qs[3])    # Subtract the real part of s <-
        # s-wave <-
        qs[3] = self.qs[3] - np.real(self.qs[3])    # Subtract the real part of s <-
        for ii in range(4):
            exponent = np.cdouble(-1.0j*(2.0*np.pi*f*qs[ii]*self.thick)/c_const)
            self.propagation_exponents[ii] = exponent
        return

class IncoherentThickLayer(CoherentLayer):
    """Define an incoherent layer using the thick slab approximation.

    Inherits from CoherentLayer
    """

    def __init__(self, layer, theta=0, phi=0, psi=0, exponent_threshold=700):
        """Initialise an instance of an incoherent layer, which uses the thick slab approximation.

        Parameters
        ----------
        layer : Layer object
            Use this to get the thickness
        theta : float, optional
            Theta value for layer (default is 0)
        phi : float, optional
            Phi value for layer (default is 0)
        psi : float, optional
            Psi value for layer (default is 0)
        exponent_threshold : float, optional
            Prevent calculation of large exponents (default is 700)

        """
        thickness = layer.getThicknessInMetres()
        epsilon = layer.getPermittivityFunction()
        Layer.__init__(self, thickness, epsilon, theta, phi, psi, exponent_threshold)
        self.coherent = False
        self.inCoherentIntensity = False
        self.inCoherentPhase = False
        self.inCoherentAveragePhase = False
        self.inCoherentThick = True
        return

    def calculate_propagation_matrix(self,f):
        """Routine to calculate the matrix Ki using the thick slab approximation.

        Parameters
        ----------
        f : float
            Frequency

        Returns
        -------
        None

        """
        Ki = super().calculate_propagation_matrix(f)
        Ki[2,2] = 0.0
        Ki[3,3] = 0.0
        return  Ki

    def calculate_scattering_matrix(self,b):
        """Calculate the scattering matrix of this layer with layer b.

        This method calculates the scattering matrix based on equations 18-21 in the PyLama paper.
        However, it's important to note that the definition of a scattering matrix used here 
        differs from that in the paper.

        Parameters
        ----------
        b : pyGTM layer
            The second layer

        Returns
        -------
        None

        """
        # First calculate the scattering matrix as usual
        # NB. the propagation matrix has zero's in the s/p backward scattering
        #     but the scattering matrix calculation uses the raw propagation exponents
        super().calculate_scattering_matrix(b)
        # Now manually set the backwards reflections for both s and p to zero
        self.SMatrix.S11 = np.zeros( (2,2) )
        self.SMatrix.calculateS()
        return

###########################
#### The System Class  ####
###########################


class System:
    """System class. An instance is an optical system with substrate, superstrate, and layers.
  
     Derived classes used by PDielec are:

        - :class:`TransferMatrixSystem`
        - :class:`ScatteringMatrixSystem`

    Attributes
    ----------
    theta : float
        Angle of incidence, in radians.
    substrate : Layer
        The substrate layer. Defaults to vacuum (an empty layer instance).
    superstrate : Layer
        The superstrate layer, defaults to vacuum (an empty layer instance).
    layers : list of Layer
        List of the layers in the system.

    Notes
    -----
    Layers can be added and removed (not inserted).

    The whole system's transfer matrix is computed using :py:func:`calculate_GammaStar`, which calls :py:func:`Layer.update` for each layer. General reflection and transmission coefficient functions are given; they require the prior execution of :py:func:`calculate_GammaStar`. The electric fields can be visualized in the case of an incident plane wave using :py:func:`calculate_Efield`.

    """

    def __init__(self, substrate=None, superstrate=None, layers=None):
        """Initialise the system for a transfer matrix calculation.

        Parameters
        ----------
        substrate : pyGTM layer, optional
            The substrate should be semi-infinite, default is None
        superstrate : pyGTM layer, optional
            The superstrate must be semi-infinite, default is None
        layers : list of pyGTM layers
             default is an empty list

        """
        self.layers=[]
        if layers is not None:
            self.layers=layers

        ## system transfer matrix
        self.Gamma = np.zeros((4,4), dtype=np.cdouble)
        self.GammaStar = np.zeros((4,4), dtype=np.cdouble)

        if substrate is not None:
            self.substrate = substrate
        else:
            self.substrate=Layer() ## should default to 1µm of vacuum
        if superstrate is not None:
            self.superstrate = superstrate
        else:
            self.superstrate=Layer() ## should default to 1µm of vacuum

    def set_substrate(self,sub):
        """Set the substrate.

        Parameters
        ----------
        sub : Layer
            Instance of the Layer class, serving as the substrate.

        Returns
        -------
        None

        """
        self.substrate=sub

    def set_superstrate(self,sup):
        """Set the superstrate.

        Parameters
        ----------
        sup : Layer
            Instance of the layer class, superstrate

        Returns
        -------
        None

        """
        self.superstrate=sup

    def get_all_layers(self):
        """Return the list of all layers in the system.

        Returns
        -------
        self.layers : list
            A list of all layers.

        Returns
        -------
        None

        """
        return self.layers

    def get_layer(self,pos):
        """Get the layer at a given position.

        Parameters
        ----------
        pos : int
            Position in the stack.

        Returns
        -------
        self.layers[pos] : Layer
            The layer at the position `pos`.

        """
        return self.layers[pos]

    def get_superstrate(self):
        """Return the System's superstrate.

        Parameters
        ----------
        None

        Returns
        -------
        self.superstrate : Layer
            The system superstrate.

        """
        return self.superstrate

    def get_substrate(self):
        """Return the System's substrate.

        Parameters
        ----------
        None

        Returns
        -------
        self.substrate : Layer
            The system substrate.

        """
        return self.substrate

    def add_layer(self,layer):
        """Add a layer instance.

        Parameters
        ----------
        layer : Layer
            The layer to be added on the stack.

        Returns
        -------
        None

        Notes
        -----
        The layers are added in *from superstrate to substrate* order.
        Light is incident *from the superstrate*.

        Note that this function adds a reference to `layer` to the list.
        If you are adding the same layer several times, be aware that if you
        change something for one of them, it changes for all of them.

        """
        self.layers.append(layer)

    def del_layer(self,pos):
        """Remove a layer at given position. Does nothing for invalid position.

        Parameters
        ----------
        pos : int
            Index of layer to be removed.

        Returns
        -------
        None

        """
        if pos >= 0 and pos < len(self.layers):
            self.layers.pop(pos)
        else:
            print("Wrong position given. No layer deleted")

    def initialize_sys(self, f):
        """Set the values of epsilon at the given frequency in all the layers.

        Parameters
        ----------
        f : float
            Frequency (Hz)

        Returns
        -------
        None

        Notes
        -----
        This function allows to define the in-plane wavevector (:math:`zeta`)
        outside of the class, and thus to explore also guided modes of the system.

        """
        self.superstrate.calculate_epsilon(f)
        self.substrate.calculate_epsilon(f)
        for li in self.layers:
            li.calculate_epsilon(f)

    def calculate_GammaStar(self,f, zeta_sys):
        """Calculate the whole system's transfer matrix.

        If one layer is incoherent then intensities rather than amplitudes are used

        Parameters
        ----------
        f : float
            Frequency (Hz)
        zeta_sys : complex
            In-plane wavevector kx/k0

        Returns
        -------
        GammaStar : complex ndarray
            4x4 complex matrix representing the system transfer matrix :math:`\\Gamma^{*}`

        """
        A_super, K_super, A_inv_super, T_super = self.superstrate.update_tm(f, zeta_sys)
        A_sub, K_sub, A_inv_sub, T_sub = self.substrate.update_tm(f, zeta_sys)
        Delta1234 = np.array([[1,0,0,0],
                              [0,0,1,0],
                              [0,1,0,0],
                              [0,0,0,1]],dtype=np.cdouble)
        incoherent = False
        Tlist = []
        # Initialise T with the substrate Transfer Matrix
        T = A_sub
        for layer in reversed(self.layers):
            Di, Pi, Di_inv, Ti = layer.update_tm(f, zeta_sys)
            if layer.inCoherentIntensity:
                # Deal with an incoherent layer
                # Finish any exisiting T matrix with Di_inv and add as a new matrix in Tlist
                # Then add Pi to Ttot
                # Then start a new T matrix with Di
                T = np.matmul(Di_inv,T)
                Tlist.append(T)
                Tlist.append(Pi)
                T = np.copy(Di)
                incoherent = True
            else:
                # coherent, so just multiply the T matrices
                T = np.matmul(Ti,T)
        # Complete with the superstrate inverse dynamic matrix
        T = np.matmul(A_inv_super,T)
        Tlist.append(T)
        if incoherent:
            # If there was incoherence then process the whole list 
            # by calculating amplitudes everywhere
            T =  np.identity(4, dtype=np.cdouble)
            for t in Tlist:
                T = np.matmul(T,np.absolute(t)**2)
            T = np.sqrt(T)
        # If there was no incoherence then T is the total transfer matrix
        Gamma = T
        GammaStar = np.matmul(exact_inv_4x4(Delta1234),np.matmul(Gamma,Delta1234))
        self.Gamma = Gamma.copy()
        self.GammaStar = GammaStar.copy()
        return self.GammaStar.copy()

    def calculate_r_t(self, zeta_sys):
        """Calculate various field and intensity reflection and transmission coefficients, as well as the 4-valued vector of transmitted field.

        Parameters
        ----------
        zeta_sys : complex
            Incident in-plane wavevector

        Returns
        -------
        r_out : array_like, shape (4,)
            Complex *field* reflection coefficients `r_out=([rpp,rps,rss,rsp])`
        R_out : array_like, shape (4,)
            Real *intensity* reflection coefficients `R_out=([Rpp,Rss,Rsp,Tps])`
        t_out : array_like, shape (4,)
            Complex *field* transmission coefficients `t=([tpp, tps, tsp, tss])`
        T_out : array_like, shape (4,)
            Real *intensity* transmission coefficients `T_out=([Tp,Ts])` (mode-inselective)

        Notes
        -----
        **IMPORTANT**

        - As of version 19-03-2020: All intensity coefficients are now well defined. Transmission is defined mode-independently. It could be defined mode-dependently for non-birefringent substrates in future versions. The new definition of this function **BREAKS compatibility** with the previous one.

        - As of version 13-09-2019: Note that the field reflectivity and transmission coefficients r and t are well defined. The intensity reflection coefficient is also correct. However, the intensity transmission coefficients T are ill-defined so far. This will be corrected upon future publication of the correct intensity coefficients.

        Note also the different ordering of the coefficients, for consistency with Passler's matlab code.

        """
        shift = sys.float_info.min
        #
        # Can we renormalise to get rid of numerical problems?
        #
        # Work with a copy of GammaStar
        largest = np.max(np.abs(self.GammaStar))
        GammaStar = self.GammaStar.copy()
        GammaStar = GammaStar / largest

        # common denominator for all coefficients
        Denom = GammaStar[0,0]*GammaStar[2,2]-GammaStar[0,2]*GammaStar[2,0]

        if Denom == 0:
            Denom += shift

        # field reflection coefficients
        rpp = GammaStar[1,0]*GammaStar[2,2]-GammaStar[1,2]*GammaStar[2,0]
        rpp = rpp/Denom

        rss = GammaStar[0,0]*GammaStar[3,2]-GammaStar[3,0]*GammaStar[0,2]
        rss = rss/Denom

        rps = GammaStar[3,0]*GammaStar[2,2]-GammaStar[3,2]*GammaStar[2,0]
        rps = rps/Denom

        rsp = GammaStar[0,0]*GammaStar[1,2]-GammaStar[1,0]*GammaStar[0,2]
        rsp = rsp/Denom

        # Try and remove any outliers by making sure the intensity is always < 1

        # Intensity reflection coefficients are just square moduli
        Rpp = np.abs(rpp)**2
        Rss = np.abs(rss)**2
        Rps = np.abs(rps)**2
        Rsp = np.abs(rsp)**2
        r_out = np.array([rpp,rps,rss,rsp]) ## order matching Passler Matlab code
        R_out = np.array([Rpp,Rss,Rsp,Rps]) ## order matching Passler Matlab code

        # field transmission coefficients
        t_out = np.zeros(4, dtype=np.cdouble)
        tpp = (GammaStar[2,2]/Denom)/largest
        tss = (GammaStar[0,0]/Denom)/largest
        tps =-(GammaStar[2,0]/Denom)/largest
        tsp =-(GammaStar[0,2]/Denom)/largest

        #CHECK the ORDERING! JK
        t_out = np.array([tpp, tps, tsp, tss],dtype=np.cdouble)

        #### Intensity transmission requires Poynting vector analysis
        ## N.B: could be done mode-dependentely later
        ## start with the superstrate
        ## Incident fields are either p or s polarized
        ksup = np.zeros((4,3), dtype=np.cdouble) ## wavevector in superstrate
        ksup[:,0] = zeta_sys
        for ii, qi in enumerate(self.superstrate.qs):
            ksup[ii,2] = qi
        ksup = ksup/c_const     ## omega simplifies in the H field formula
        Einc_pin = self.superstrate.gamma[0,:] ## p-pol incident electric field
        Einc_sin = self.superstrate.gamma[1,:] ## s-pol incident electric field
        ## Poynting vector in superstrate (incident, p-in and s-in)
        Sinc_pin = 0.5*np.real(np.cross(Einc_pin,np.conj(np.cross(ksup[0,:],Einc_pin))))
        Sinc_sin = 0.5*np.real(np.cross(Einc_sin,np.conj(np.cross(ksup[1,:],Einc_sin))))

        ### Substrate Poynting vector
        ## Outgoing fields (eqn 17)
        Eout_pin = t_out[0]*self.substrate.gamma[0,:]+t_out[1]*self.substrate.gamma[1,:] #p-in, p or s out
        Eout_sin = t_out[2]*self.substrate.gamma[0,:]+t_out[3]*self.substrate.gamma[1,:] #s-in, p or s out
        ksub = np.zeros((4,3), dtype=np.cdouble)
        ksub[:,0] = zeta_sys
        for ii, qi in enumerate(self.substrate.qs):
            ksub[ii,2] = qi
        ksub = ksub/c_const ## omega simplifies in the H field formula

        ###########################
        ## outgoing Poynting vectors, 2 formulations
        Sout_pin = 0.5*np.real(np.cross(Eout_pin,np.conj(np.cross(ksub[0,:],Eout_pin))))
        Sout_sin = 0.5*np.real(np.cross(Eout_sin,np.conj(np.cross(ksub[1,:],Eout_sin))))
        ### Intensity transmission coefficients are only the z-component of S !
        T_pp = (Sout_pin[2]/Sinc_pin[2]) ## z-component only
        T_ss = (Sout_sin[2]/Sinc_sin[2]) ## z-component only

        T_out = np.array([T_pp, T_ss])

        return r_out, R_out, t_out, T_out

    def overflowErrors(self):
        """Return the total number of overflow errors encountered.

        Parameters
        ----------
        None

        Returns
        -------
        int
            The total number of overflow errors encountered.

        """
        count = 0
        largest_exponent = 0.0
        for layer in self.layers:
            count = count + layer.exponent_errors
            largest_exponent = max(layer.largest_exponent,largest_exponent)
        return count,largest_exponent

class TransferMatrixSystem(System):
    """Define a system of layers which is described by a Transfer Matrix method.
    
    This is basically a place holder for the original pyGTM methods and it is used
    to distinguish itself from a system which uses the scattering matrix method
    """

    def __init__(self, substrate=None, superstrate=None, layers=None):
        """Initialize a new instance of the system with optional substrate, superstrate, and layers.

        Use for transfer matrix method calculations using the original pyGTM methods

        Parameters
        ----------
        substrate : Layer, optional
            The substrate layer used in the system. The default is None.
        superstrate : Layer, optional
            The superstrate layer used in the system. The default is None.
        layers : list of layers, optional
            A list of layer objects used in the system. The default is an empty list.

        """        
        System.__init__(self, substrate, superstrate, layers)
        return

class ScatteringMatrixSystem(System):
    """Define a system of layers which is described by a Scattering Matrix method.
    
    This implements additional routines needed for the scattering matrix method

    """

    def __init__(self, substrate=None, superstrate=None, layers=None):
        """Initialize a new instance of the System class suitable for scattering matrix calculations.

        Parameters
        ----------
        substrate : object, optional
            The substrate component of the system. The default is None.
        superstrate : object, optional
            The superstrate component of the system. The default is None.
        layers : list, optional
            A list representing the layers in the system. Defaults to an empty list.

        """        
        System.__init__(self, substrate, superstrate, layers)
        return

    def calculate_scattering_matrices(self):
        """Calculate the scattering matrices in every layer.

        Loop through all the layers calculating the scattering matrix of this layer with the next.
        Care is taken when dealing with the superstrate and the substrate

        Parameters
        ----------
        None

        Returns
        -------
        None

        Note
        ----
        The scattering matrix is stored as an attribute 
        If the layer system is changed in anyway all the scattering matrices must be recomputed

        """
        # First the superstrate with the first layer in the device
        if len(self.layers) < 1:
            self.superstrate.calculate_scattering_matrix(self.substrate)
        else:
            self.superstrate.calculate_scattering_matrix(self.layers[0])
            #Utilities.printsp('SM super-0 ',self.superstrate.SMatrix.S)
        # Next every layer in the devices with the next layer, but not the last layer in the device
        for index,layer in enumerate(self.layers[:-1]):
            layer.calculate_scattering_matrix(self.layers[index+1])
            #Utilities.printsp('SM {}-{}'.format(index,index+1),layer.SMatrix.S)
        # Now calculate the scattering matrix between the last layer in the device and the substrate
        if len(self.layers) > 0:
            self.layers[-1].calculate_scattering_matrix(self.substrate)
            #Utilities.printsp('SM last-substrate',self.layers[-1].SMatrix.S)
        # Note that the substrate does not have a scattering matrix
        return

    def calculate_GammaStar(self,f, zeta_sys):
        """Calculate the whole system's scattering matrix.

        This replaces the equivalent routine used in the transfer matrix method.

        Parameters
        ----------
        f : float
            Frequency (Hz)
        zeta_sys : complex
            In-plane wavevector kx/k0

        Returns
        -------
        Stotal : ndarray
            System scattering matrix :math:`S`

        """
        # Use the super class to calculate GammaStar
        #super().calculate_GammaStar(f,zeta_sys)

        self.superstrate.update_sm(f, zeta_sys)
        self.substrate.update_sm(f, zeta_sys)
        for layer in self.layers:
            layer.update_sm(f, zeta_sys)
        # Ask for calculation of scattering matrices in all layers
        self.calculate_scattering_matrices()
        # Loop through all the layers to calculate the total scattering matrix
        self.Stotal = copy.copy(self.superstrate.SMatrix)
        for layer in self.layers:
            self.Stotal = self.Stotal.redheffer(layer.SMatrix)
        #self.Stotal = self.Stotal.redheffer(self.substrate.SMatrix)
        return self.Stotal

    def calculate_r_t(self, zeta_sys):
        """Calculate the reflectance and transmittance coefficients using scattering matrix information.

        Calculate various field and intensity reflection and transmission coefficients, as well as the 4-valued vector of transmitted field.

        Parameters
        ----------
        zeta_sys : complex
            Incident in-plane wavevector

        Returns
        -------
        r_out : array_like, shape (4,)
            Complex *field* reflection coefficients (r_out=([rpp,rps,rss,rsp])).
        R_out : array_like, shape (4,)
            Real *intensity* reflection coefficients (R_out=([Rpp,Rss,Rsp,Tps])).
        t_out : array_like, shape (4,)
            Complex *field* transmission coefficients (t=([tpp, tps, tsp, tss])).
        T_out : array_like, shape (4,)
            Real *intensity* transmission coefficients (T_out=([Tp,Ts])) (mode-inselective).

        """
        # Extract the scattering matrix
        S = self.Stotal.S

        # field reflection coefficients
        rpp = S[0,0]
        rss = S[1,1]
        rps = S[0,1]
        rsp = S[1,0]

        # Intensity reflection coefficients are just square moduli
        Rpp = np.abs(rpp)**2
        Rss = np.abs(rss)**2
        Rps = np.abs(rps)**2
        Rsp = np.abs(rsp)**2
        r_out = np.array([rpp,rps,rss,rsp]) ## order matching Passler Matlab code
        R_out = np.array([Rpp,Rss,Rsp,Rps]) ## order matching Passler Matlab code

        # field transmission coefficients
        tpp = S[2,0]
        tss = S[3,1]
        tps = S[2,1]
        tsp = S[3,0]
        t_out = np.array([tpp, tps, tsp, tss],dtype=np.cdouble)

        # Intensity transmission using Poynting vector analysis
        # start with the superstrate
        # Incident fields are either p or s polarized
        ksup = np.zeros((4,3), dtype=np.cdouble) ## wavevector in superstrate
        ksup[:,0] = zeta_sys
        for ii, qi in enumerate(self.superstrate.qs):
            ksup[ii,2] = qi
        ksup = ksup/c_const                     # omega simplifies in the H field formula
        # Incident electric field
        Einc_pin = self.superstrate.gamma[0,:]  # p-pol incident electric field
        Einc_sin = self.superstrate.gamma[1,:]  # s-pol incident electric field

        # Poynting vector in superstrate (incident, p-in and s-in)
        Sinc_pin = 0.5*np.real(np.cross(Einc_pin,np.conj(np.cross(ksup[0,:],Einc_pin))))
        Sinc_sin = 0.5*np.real(np.cross(Einc_sin,np.conj(np.cross(ksup[1,:],Einc_sin))))

        # Substrate Poynting vector
        Eout_pin = t_out[0]*self.substrate.gamma[0,:]+t_out[1]*self.substrate.gamma[1,:] #p-in, p or s out
        Eout_sin = t_out[2]*self.substrate.gamma[0,:]+t_out[3]*self.substrate.gamma[1,:] #s-in, p or s out
        ksub = np.zeros((4,3), dtype=np.cdouble)
        ksub[:,0] = zeta_sys
        for ii, qi in enumerate(self.substrate.qs):
            ksub[ii,2] = qi
        ksub = ksub/c_const # omega simplifies in the H field formula

        # outgoing Poynting vectors, 2 formulations
        Sout_pin = 0.5*np.real(np.cross(Eout_pin,np.conj(np.cross(ksub[0,:],Eout_pin))))
        Sout_sin = 0.5*np.real(np.cross(Eout_sin,np.conj(np.cross(ksub[1,:],Eout_sin))))

        # Intensity transmission coefficients are only the z-component of S !
        T_pp = (Sout_pin[2]/Sinc_pin[2]) ## z-component only
        T_ss = (Sout_sin[2]/Sinc_sin[2]) ## z-component only

        T_out = np.array([T_pp, T_ss])

        return r_out, R_out, t_out, T_out

class SMatrix:
    """A class for storing and manipulating scattering matrices.

    S-matrices are square matrices used in physics and engineering to describe the 
    scattering and reflection of waves or particles. They are especially useful in the
    field of microwave engineering, quantum mechanics, and optical physics.

    Attributes
    ----------
    S11 : np.ndarray
        The top-left 2x2 submatrix of the S-matrix, representing input to input scattering.
    S22 : np.ndarray
        The bottom-right 2x2 submatrix of the S-matrix, representing output to output scattering.
    S21 : np.ndarray
        The bottom-left 2x2 submatrix of the S-matrix, representing input to output transmission.
    S12 : np.ndarray
        The top-right 2x2 submatrix of the S-matrix, representing output to input transmission.
    S : np.ndarray
        A 4x4 composite S-matrix combining S11, S22, S21, and S12.

    Methods
    -------
    unitMatrix()
        Initializes the S-matrix to a unit matrix with appropriate S11, S22, S21, and S12.
    redheffer(b)
        Performs the Redheffer star product, a specialized matrix multiplication for S-matrices, 
        with another SMatrix instance `b`.
    calculateS()
        Reconstructs the composite 4x4 S-matrix from its constituent submatrices (S11, S22, S21, S12).

    Examples
    --------
    >>> import numpy as np
    >>> from smatrix import SMatrix  # Assuming this class is saved in smatrix.py
    >>> s_matrix_data = np.array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15]])
    >>> smatrix = SMatrix(S=s_matrix_data)
    >>> print(smatrix.S)  # This displays the initialized S-matrix

    Notes
    -----
    The Redheffer star product is particularly useful for cascading two-port networks in
    microwave engineering and quantum physics. It allows for the efficient calculation of
    the overall scattering matrix of the combined system.

    """

    def __init__(self, S=None):
        """Initialise a class to store a scattering matrix and to allow multiplication of matrices.

        Parameters
        ----------
        S : tensor (4x4)
            A np.darray 4x4 tensor used to initialise the scattering matrix

        """
        if S is None:
            self.unitMatrix()
        else:
            self.S11 = S[0:2,0:2]
            self.S22 = S[2:4,2:4]
            self.S21 = S[2:4,0:2]
            self.S12 = S[0:2,2:4]
            self.calculateS()
        return

    def unitMatrix(self):
        """Reset the scattering parameters to form a unit matrix.

        Resets the scattering parameters S11, S22, S21, and S12 of the instance to represent a unit matrix, where S11 and S22 are set to zero matrices and S21 and S12 are set to identity matrices. After the reset, calculates the scattering parameters through `calculateS`.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method is typically used to initialize or reset the scattering parameters in a two-port network representation to a known state before performing further operations or calculations.

        """        
        self.S11 = np.zeros( (2,2) )
        self.S22 = np.zeros( (2,2) )
        self.S21 = np.eye( 2 )
        self.S12 = np.eye( 2 )
        self.calculateS()
        return

    def redheffer(self,b):
        """Calculate the product of S-matrices.

        Parameters
        ----------
        b : ScatteringMatrix
           The second scattering matrix to be used the the multiplication

        Returns
        -------
        Sab : ScatteringMatrix
            The product of Sa and Sb

        """
        a = self
        # If b is none then the S-matrix was not initialised (could be the substrate)
        if b is None:
            b = SMatrix()
        C = np.eye(2,dtype=np.cdouble) - ( b.S11.dot(a.S22) ) 
        D = np.eye(2,dtype=np.cdouble) - ( a.S22.dot(b.S11) ) 
        C = exact_inv_2x2(C)
        D = exact_inv_2x2(D)
        Sab = SMatrix()
        Sab.S11 = a.S11+a.S12.dot(C).dot(b.S11).dot(a.S21)
        Sab.S12 = a.S12.dot(C).dot(b.S12)
        Sab.S21 = b.S21.dot(D).dot(a.S21)
        Sab.S22 = b.S22+b.S21.dot(D).dot(a.S22).dot(b.S12)
        Sab.calculateS()
        return Sab

    def calculateS(self):
        """Calculate and assign the scattering parameter matrix, S.

        Fills a 4 by 4 matrix, `S`, with the sub-matrices `S11`, `S22`, `S12`, and `S21` representing the scattering parameters of a two-port network. The `S` matrix combines these sub-matrices as follows:

        - `S11` is placed in the top-left quadrant.
        - `S22` is placed in the bottom-right quadrant.
        - `S12` is placed in the top-right quadrant.
        - `S21` is placed in the bottom-left quadrant.

        Parameters
        ----------
        None

        Returns
        -------
        S : ndarray, shape (4, 4), dtype=np.cdouble
            The complete scattering parameter matrix of the two-port network.

        Notes
        -----
        The scattering parameters (`S11`, `S22`, `S12`, `S21`) must be defined as attributes of the instance before calling this method.

        """        
        self.S = np.zeros( (4,4),dtype=np.cdouble )
        self.S[0:2,0:2] = self.S11
        self.S[2:4,2:4] = self.S22
        self.S[0:2,2:4] = self.S12
        self.S[2:4,0:2] = self.S21
        return self.S
      
