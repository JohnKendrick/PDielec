# This file is part of the pyGTM module.
# It has been heavily modified by John Kendrick for inclusion in the pdielec library
# Quite a lot has been removed as the fields spatial distribution was not needed
# Thanks to the authors of the original code for all their hard work
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) Mathieu Jeannin 2019 2020 <math.jeannin@free.fr>.

"""
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
    
    *15-10-2021*:
        
        - Fixed rounding error bug in lag.eig() causing the program to crash randomly
        for negligibly small imaginary parts of the wavevectors
        
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

import numpy as np
import scipy.linalg as lag
import sys

c_const = 299792458 # m/s
eps0 = 8.854e-12 ## vacuum permittivity
qsd_thr = 1e-10 ### threshold for wavevector comparison
zero_thr = 1e-10 ### threshold for eigenvalue comparison to zero

def vacuum_eps(f):
    """
    Vacuum permittivity function

    Parameters
    ----------
    f: float or 1D-array
       frequency (in Hz)

    Returns
    -------
    eps : complex or 1D-array of complex
        Complex value of the vacuum permittivity (1.0 + 0.0j)
    """
    try:
        return np.ones(len(f))
    except:
        return 1.0 + 0.0j
    
    
def exact_inv(M):
    """Compute the 'exact' inverse of a 4x4 matrix using the analytical result.

    Parameters
    ----------
    M : 4X4 array (float or complex)
      Matrix to be inverted

    Returns
    -------
    out : 4X4 array (complex)
        Inverse of this matrix or Moore-Penrose approximation if matrix cannot be inverted.

    Notes
    -----
    This should give a higher precision and speed at a reduced noise.
    From D.Dietze code https://github.com/ddietze/FSRStools

    .. seealso:: http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html

    Modified to work with clongdouble by JK

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
            print('Warning 4x4 inversion problem 1')
            result = np.clongdouble(lag.pinv(np.cdouble(M)))
        except:
            print('Warning 4x4 inversion problem 2')
            result = np.clongdouble(lag.pinv(np.complex64(M)))
        return result

    B = np.zeros(A.shape, dtype=np.clongdouble)
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

    out = B.T / detA
    return out


#%%
###########################
#### The Layer Class   ####
###########################
class Layer:
    """
    Layer class. An instance is a single layer:

    Attributes
    -----------
    thickness : float
              thickness of the layer in m
    epsilon1 : complex function
             function epsilon(frequency) for the first axis. If none, defaults to vacuum.
    epsilon2 : complex function
             function epsilon(frequency) for the second axis. If none, defaults to epsilon1.
    epsilon3 : complex function
             function epsilon(frequency) for the third axis. If none, defaults to epsilon1.
    epsilon  : complex function   ADDED by JK
             function epsilon(frequency) the full dielectric constant tensor, if not none all other epsilons are ignored
    theta : float
          Euler angle theta (colatitude) in rad
    phi : float
        Euler angle phi in rad
    psi : float
        Euler angle psi in rad

    Notes
    -----
    If instanciated with defaults values, it generates a 1um thick layer of air.
    Properties can be checked/changed dynamically using the corresponding get/set methods.
    """

    def __init__(self, thickness=1.0e-6, epsilon1=None, epsilon2=None, epsilon3=None,  epsilon=None,
                                         theta=0, phi=0, psi=0):

        ## epsilon is a 3x3 matrix of permittivity at a given frequency
        self.epsilon = np.identity(3, dtype=np.clongdouble)
        self.mu = 1.0 ### mu=1 for now

        ### initialization of all important quantities
        self.M = np.zeros((6, 6), dtype=np.clongdouble) ## constitutive relations
        self.a = np.zeros((6, 6), dtype=np.clongdouble) ##
        self.S = np.zeros((4, 4), dtype=np.clongdouble) ##
        self.Delta = np.zeros((4, 4), dtype=np.clongdouble) ##
        self.qs = np.zeros(4, dtype=np.clongdouble) ## out of plane wavevector
        self.Py = np.zeros((3,4), dtype=np.clongdouble) ## Poyting vector
        self.gamma = np.zeros((4, 3), dtype=np.clongdouble) ##
        self.Ai = np.zeros((4, 4), dtype=np.clongdouble) ##
        self.Ki = np.zeros((4, 4), dtype=np.clongdouble) ##
        self.Ti = np.zeros((4, 4), dtype=np.clongdouble) ## Layer transfer matrix
        self.Berreman = np.zeros((4,3), dtype=np.clongdouble) ## Stores the Berreman modes, used for birefringent layers
        self.useBerreman = False ### Boolean to replace Xu's eigenvectors by Berreman's in case of Birefringence
        
        self.euler = np.identity(3, dtype=np.clongdouble) ## rotation matrix
        self.epsilon_tensor_function = None              ## Added by JK

        self.set_thickness(thickness) ## set the thickness, 1um by default
        self.set_epsilon(epsilon1, epsilon2, epsilon3, epsilon) # set epsilon, vacuum by default (JK change)
        self.set_euler(theta, phi, psi) ## set orientation of crystal axis w/ respect to the lab frame


    def set_thickness(self, thickness):
        """
        Sets the layer thickness

        Parameters
        ----------
        thickness : float
                  the layer thickness (in m)

        Returns
        -------
               None
        """
        self.thick = thickness

    def set_epsilon(self, epsilon1=vacuum_eps, epsilon2=None, epsilon3=None, epsilon=None):
        """
        Sets the dielectric functions for the three main axis.

        Parameters
        -----------
        epsilon1 : complex function
                 function epsilon(frequency) for the first axis. If none, defaults to :py:func:`vacuum_eps`
        epsilon2 : complex function
                 function epsilon(frequency) for the second axis. If none, defaults to epsilon1.
        epsilon3 : complex function
                 function epsilon(frequency) for the third axis. If none, defaults to epsilon1.
        epsilon  : complex function
                 JK modification to allow the full dielectric matrix to be supplied.  If specified
                 all other definitions are overridden
        func epsilon1: function returning the first (xx) component of the complex permittivity tensor in the crystal frame.

        Returns
        -------
               None

        Notes
        ------
        Each *epsilon_i* function returns the dielectric constant along axis i as
        a function of the frequency f in Hz.

        If no function is given for epsilon1, it defaults to :py:func:`vacuum_eps` (1.0 everywhere).
        epsilon2 and epsilon3 default to epsilon1: if None, a homogeneous material is assumed
        JK Modification to allow a fully complex tensor given by epsilon
        JK if epsilon is not specified it behaves as before
        JK if epsilon is given the epsilon_tensor function is defined by epsilon
        """
        if epsilon == None:
            if epsilon1==None:
                self.epsilon1_f = vacuum_eps
            else:
                self.epsilon1_f = epsilon1

            if epsilon2 == None:
                self.epsilon2_f = self.epsilon1_f
            else:
                self.epsilon2_f = epsilon2
            if epsilon3 == None:
                self.epsilon3_f = self.epsilon1_f
            else:
                self.epsilon3_f = epsilon3
        else:
            self.epsilon1_f = None
            self.epsilon2_f = None
            self.epsilon3_f = None
            self.epsilon_tensor_function = epsilon

    def calculate_epsilon(self, f):
        """
        Sets the value of epsilon in the (rotated) lab frame.

        Parameters
        ----------
        f : float
            frequency (in Hz)
        Returns
        -------
            None

        Notes
        ------
        The values are set according to the epsilon_fi (i=1..3) functions
        defined using the :py:func:`set_epsilon` method, at the given frequency f.
        The rotation with respect to the lab frame is computed using the Euler angles.

        Use only explicitely if you *don't* use the :py:func:`Layer.update` function!
        Modification by JK to allow the use of a full dielectric tensor
        """
        if self.epsilon_tensor_function == None:
            epsilon_xstal = np.zeros((3,3), dtype=np.clongdouble)
            epsilon_xstal[0,0] = self.epsilon1_f(f)
            epsilon_xstal[1,1] = self.epsilon2_f(f)
            epsilon_xstal[2,2] = self.epsilon3_f(f)
        else:
            epsilon_xstal = self.epsilon_tensor_function(f)
        # self.epsilon = np.matmul(self.euler_inverse, np.matmul(epsilon_xstal,self.euler))
        #JK Changed the order of the transformation to a transformation of basis vectors
        # rather than a transformation of the coordinates active->passive
        self.epsilon = np.matmul(self.euler, np.matmul(epsilon_xstal,self.euler_inverse))
        return self.epsilon


    def set_euler(self,theta,phi,psi):
        """
        Sets the values for the Euler rotations angles.

        Parameters
        ----------
        theta : float
              Euler angle theta (colatitude) in rad
        phi : float
            Euler angle phi in rad
        psi : float
            Euler angle psi in rad

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
        #JK self.euler_inverse = np.clongdouble(lag.pinv(np.cdouble(self.euler)))
        self.euler_inverse = self.invert(self.euler)

    def invert(self,m):
        """Calculate the inverse of m

        Parameters
        ----------
        m : complex 3x3 array

        Returns
        -------
        the inverse as a 3x3 complex np array
        """
        m1, m2, m3, m4, m5, m6, m7, m8, m9 = m.flatten()
        determinant = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9  
        return np.array([[m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5],
                         [m6*m7-m4*m9, m1*m9-m3*m7, m3*m4-m1*m6],
                         [m4*m8-m5*m7, m2*m7-m1*m8, m1*m5-m2*m4]])/determinant


    def calculate_matrices(self, zeta):
        """
        Calculate the principal matrices necessary for the GTM algorithm.

        Parameters
        ----------
        zeta : complex
             In-plane reduced wavevector kx/k0 in the system.

        Returns
        -------
             None

        Notes
        -----
        Note that zeta is conserved through the whole system and set externaly
        using the angle of incidence and `System.superstrate.epsilon[0,0]` value

        Requires prior execution of :py:func:`calculate_epsilon`

        """
        ## Constitutive matrix (see e.g. eqn (4))
        self.M[0:3, 0:3] = np.clongdouble(self.epsilon.copy())
        self.M[3:6, 3:6] = np.clongdouble(self.mu*np.identity(3))

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
        self.S[0,0] = self.M[0,0] + self.M[0,2]*self.a[2,0] + self.M[0,5]*self.a[5,0];
        self.S[0,1] = self.M[0,1] + self.M[0,2]*self.a[2,1] + self.M[0,5]*self.a[5,1];
        self.S[0,2] = self.M[0,3] + self.M[0,2]*self.a[2,3] + self.M[0,5]*self.a[5,3];
        self.S[0,3] = self.M[0,4] + self.M[0,2]*self.a[2,4] + self.M[0,5]*self.a[5,4];
        self.S[1,0] = self.M[1,0] + self.M[1,2]*self.a[2,0] + (self.M[1,5]-zeta)*self.a[5,0];
        self.S[1,1] = self.M[1,1] + self.M[1,2]*self.a[2,1] + (self.M[1,5]-zeta)*self.a[5,1];
        self.S[1,2] = self.M[1,3] + self.M[1,2]*self.a[2,3] + (self.M[1,5]-zeta)*self.a[5,3];
        self.S[1,3] = self.M[1,4] + self.M[1,2]*self.a[2,4] + (self.M[1,5]-zeta)*self.a[5,4];
        self.S[2,0] = self.M[3,0] + self.M[3,2]*self.a[2,0] + self.M[3,5]*self.a[5,0];
        self.S[2,1] = self.M[3,1] + self.M[3,2]*self.a[2,1] + self.M[3,5]*self.a[5,1];
        self.S[2,2] = self.M[3,3] + self.M[3,2]*self.a[2,3] + self.M[3,5]*self.a[5,3];
        self.S[2,3] = self.M[3,4] + self.M[3,2]*self.a[2,4] + self.M[3,5]*self.a[5,4];
        self.S[3,0] = self.M[4,0] + (self.M[4,2]+zeta)*self.a[2,0] + self.M[4,5]*self.a[5,0];
        self.S[3,1] = self.M[4,1] + (self.M[4,2]+zeta)*self.a[2,1] + self.M[4,5]*self.a[5,1];
        self.S[3,2] = self.M[4,3] + (self.M[4,2]+zeta)*self.a[2,3] + self.M[4,5]*self.a[5,3];
        self.S[3,3] = self.M[4,4] + (self.M[4,2]+zeta)*self.a[2,4] + self.M[4,5]*self.a[5,4];


        ## Delta Matrix from eqn (8)
        self.Delta[0,0] = self.S[3,0]
        self.Delta[0,1] = self.S[3,3]
        self.Delta[0,2] = self.S[3,1]
        self.Delta[0,3] = - self.S[3,2]
        self.Delta[1,0] = self.S[0,0]
        self.Delta[1,1] = self.S[0,3]
        self.Delta[1,2] = self.S[0,1]
        self.Delta[1,3] = - self.S[0,2]
        self.Delta[2,0] = -self.S[2,0]
        self.Delta[2,1] = -self.S[2,3]
        self.Delta[2,2] = -self.S[2,1]
        self.Delta[2,3] = self.S[2,2]
        self.Delta[3,0] = self.S[1,0]
        self.Delta[3,1] = self.S[1,3]
        self.Delta[3,2] = self.S[1,1]
        self.Delta[3,3] = -self.S[1,2]

    def calculate_q(self):
        """
        Calculates the 4 out-of-plane wavevectors for the current layer.

        Returns
        -------
        None


        Notes
        -----
        From this we also get the Poynting vectors.
        Wavevectors are sorted according to (trans-p, trans-s, refl-p, refl-s)
        Birefringence is determined according to a threshold value `qsd_thr`
        set at the beginning of the script.
        """
        Delta_loc = np.zeros((4,4), dtype=np.cdouble)
        transmode = np.zeros((2), dtype=int)
        reflmode = np.zeros((2), dtype=int)
        
        Delta_loc = np.cdouble(self.Delta.copy())
        ## eigenvals // eigenvects as of eqn (11)
        qsunsorted, psiunsorted = lag.eig(Delta_loc)
        ##### remove extremely small real/imaginary parts that are due to numerical inaccuracy
        for km in range(4):
            if (np.abs(np.imag(qsunsorted[km])) > 0) and (np.abs(np.imag(qsunsorted[km])) < zero_thr):
                qsunsorted[km] = np.real(qsunsorted[km]) + 0.0j
            if (np.abs(np.real(qsunsorted[km])) > 0) and (np.abs(np.real(qsunsorted[km])) < zero_thr):
                qsunsorted[km] = 0.0 + 1.0j*np.imag(qsunsorted[km])
        for comp in range(4):
            if (np.abs(np.real(psiunsorted[km][comp]))>0) and (np.abs(np.real(psiunsorted[km][comp])) < zero_thr):
                psiunsorted[km][comp] = 0.0 + 1.0j*np.imag(psiunsorted[km][comp])
            if (np.abs(np.imag(psiunsorted[km][comp]))>0) and (np.abs(np.imag(psiunsorted[km][comp])) < zero_thr):
                psiunsorted[km][comp] = np.real(psiunsorted[km][comp]) + 0.0j

                
        Berreman_unsorted = np.zeros((4,3), dtype=np.clongdouble)
        
        kt = 0 
        kr = 0;
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
        shift = sys.float_info.min
        Cp_t1 = np.abs(self.Py[0,transmode[0]])**2/(np.abs(self.Py[0,transmode[0]])**2+np.abs(self.Py[1,transmode[0]])**2+shift)
        Cp_t2 = np.abs(self.Py[0,transmode[1]])**2/(np.abs(self.Py[0,transmode[1]])**2+np.abs(self.Py[1,transmode[1]])**2+shift)

        if np.abs(Cp_t1-Cp_t2) > qsd_thr: ## birefringence
            self._useBerreman = True ## sets _useBerreman fo the calculation of gamma matrix below
            if Cp_t2>Cp_t1:
                transmode = np.flip(transmode,0) ## flip the two values
            ## then calculate for reflected waves if necessary
            Cp_r1 = np.abs(self.Py[0,reflmode[1]])**2/(np.abs(self.Py[0,reflmode[1]])**2+np.abs(self.Py[1,reflmode[1]])**2+shift)
            Cp_r2 = np.abs(self.Py[0,reflmode[0]])**2/(np.abs(self.Py[0,reflmode[0]])**2+np.abs(self.Py[1,reflmode[0]])**2+shift)
            if Cp_r1>Cp_r2:
                reflmode = np.flip(reflmode,0) ## flip the two values

        else:     ### No birefringence, use the Electric field s-pol/p-pol
            Cp_te1 = np.abs(psiunsorted[0,transmode[1]])**2/(np.abs(psiunsorted[0,transmode[1]])**2+np.abs(psiunsorted[2,transmode[1]])**2+shift)
            Cp_te2 = np.abs(psiunsorted[0,transmode[0]])**2/(np.abs(psiunsorted[0,transmode[0]])**2+np.abs(psiunsorted[2,transmode[0]])**2+shift)
            if Cp_te1>Cp_te2:
                transmode = np.flip(transmode,0) ## flip the two values
            Cp_re1 = np.abs(psiunsorted[0,reflmode[1]])**2/(np.abs(psiunsorted[0,reflmode[1]])**2+np.abs(psiunsorted[2,reflmode[1]])**2+shift)
            Cp_re2 = np.abs(psiunsorted[0,reflmode[0]])**2/(np.abs(psiunsorted[0,reflmode[0]])**2+np.abs(psiunsorted[2,reflmode[0]])**2+shift)
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
        
    def calculate_gamma(self, zeta):
        """
        Calculate the gamma matrix
        
        Parameters
        ----------
        zeta : complex
             in-plane reduced wavevector kx/k0

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

        if np.abs(self.qs[0]-self.qs[1])<qsd_thr:
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

        if np.abs(self.qs[2]-self.qs[3])<qsd_thr:
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
        gamma1 = np.array([[self.gamma[0,0], gamma12, gamma13]],dtype=np.clongdouble)
        gamma2 = np.array([[gamma21, self.gamma[1,1], gamma23]],dtype=np.clongdouble)
        gamma3 = np.array([[self.gamma[2,0], gamma32, gamma33]],dtype=np.clongdouble)
        gamma4 = np.array([[gamma41, self.gamma[3,1], gamma43]],dtype=np.clongdouble)
        
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
            print('replaced gamma by Berreman')
            self.gamma = self.Berreman
        
        
    def calculate_transfer_matrix_mp(self, f, zeta):
        """
        Compute the transfer matrix of the whole layer :math:`T_i=A_iP_iA_i^{-1}`
        Uses arbitrary precision maths

        Parameters
        ----------
        f : float
            frequency (in Hz)
        zeta : complex
               reduced in-plane wavevector kx/k0
        Returns
        -------
        None

        """
        ## eqn(22)
        #from mpmath import *
        #mp.dps = 20
        #self.Ai[0,:] = self.gamma[:,0].copy()
        #self.Ai[1,:] = self.gamma[:,1].copy()
        
        Ai[0,:] = self.gamma[:,0].copy()
        Ai[1,:] = self.gamma[:,1].copy()

        self.Ai[2,:] = (self.qs*self.gamma[:,0]-zeta*self.gamma[:,2])/self.mu
        self.Ai[3,:] = self.qs*self.gamma[:,1]/self.mu

        for ii in range(4):
            ## looks a lot like eqn (25). Why is K not Pi ?
            exponent = np.clongdouble(-1.0j*(2.0*np.pi*f*self.qs[ii]*self.thick)/c_const)
            if np.abs(exponent) > exponent_threshold:
                exponent = exponent/np.abs(exponent)*exponent_threshold
            self.Ki[ii,ii] = np.nan_to_num(np.exp(exponent))
            #  JK original line
            # self.Ki[ii,ii] = np.exp(-1.0j*(2.0*np.pi*f*self.qs[ii]*self.thick)/c_const)

        Aim1 = exact_inv(self.Ai.copy())
        ## eqn (26)
        self.Ti = np.matmul(self.Ai,np.matmul(self.Ki,Aim1))

    def calculate_transfer_matrix(self, f, zeta):
        """
        Compute the transfer matrix of the whole layer :math:`T_i=A_iP_iA_i^{-1}`

        Parameters
        ----------
        f : float
            frequency (in Hz)
        zeta : complex
               reduced in-plane wavevector kx/k0
        Returns
        -------
        None

        """
        exponent_threshold = 400.0
        ## eqn(22)
        self.Ai[0,:] = self.gamma[:,0].copy()
        self.Ai[1,:] = self.gamma[:,1].copy()

        self.Ai[2,:] = (self.qs*self.gamma[:,0]-zeta*self.gamma[:,2])/self.mu
        self.Ai[3,:] = self.qs*self.gamma[:,1]/self.mu
   
        for ii in range(4):
            ## looks a lot like eqn (25). Why is K not Pi ?
            exponent = np.clongdouble(-1.0j*(2.0*np.pi*f*self.qs[ii]*self.thick)/c_const)
            if np.abs(exponent) > exponent_threshold:
                exponent = exponent/np.abs(exponent)*exponent_threshold
            self.Ki[ii,ii] = np.nan_to_num(np.exp(exponent))
            #  JK original line
            # self.Ki[ii,ii] = np.exp(-1.0j*(2.0*np.pi*f*self.qs[ii]*self.thick)/c_const)
        Aim1 = exact_inv(self.Ai.copy())
        ## eqn (26)
        self.Ti = np.matmul(self.Ai,np.matmul(self.Ki,Aim1))

    def update(self, f, zeta):
        """Shortcut to recalculate all layer properties.

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
        Ai_inv = exact_inv(self.Ai.copy())

        return[self.Ai.copy(), self.Ki.copy(), Ai_inv.copy(), self.Ti.copy()]



#%%
###########################
#### The System Class  ####
###########################


class System:
    """
    System class. An instance is an optical system with substrate, superstrate and layers.

    Attributes
    ----------
    theta : float
            Angle of incidence, in radians
    substrate : Layer
            The substrate layer. Defaults to vacuum (empty layer instance)
    superstrate : Layer
            The superstrate layer, defaults to vacuum (empty layer instance)
    layers : list of layers
            list of the layers in the system

    Notes
    -----
    Layers can be added and removed (not inserted).

    The whole system's transfer matrix is computed using :py:func:`calculate_GammaStar`,
    which calls :py:func:`Layer.update` for each layer.
    General reflection and transmission coefficient functions are given, they require prior
    execution of :py:func:`calculate_GammaStar`.
    The electric fields can be visualized in the case of incident plane wave
    using :py:func:`calculate_Efield`

    """
    def __init__(self, substrate=None, superstrate=None, layers=[]):#,
                 #theta=0.0, phi=0.0, psi=0.0):

        self.layers=[]
        if len(layers)>0:
            self.layers=layers

        ## system transfer matrix
        self.Gamma = np.zeros((4,4), dtype=np.clongdouble)
        self.GammaStar = np.zeros((4,4), dtype=np.clongdouble)

        if substrate is not None:
            self.substrate = substrate
        else:
            self.substrate=Layer() ## should default to 1µm of vacuum
        if superstrate is not None:
            self.superstrate = superstrate
        else:
            self.superstrate=Layer() ## should default to 1µm of vacuum

    def set_substrate(self,sub):
        """Sets the substrate

        Parameters
        ----------
        sub : Layer
            Instance of the layer class, substrate
        Returns
        -------
        None
        """
        self.substrate=sub

    def set_superstrate(self,sup):
        """Set the superstrate

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
        """Returns the list of all layers in the system

        Returns
        -------
        l : list
            list of all layers
        """
        return self.layers

    def get_layer(self,pos):
        """Get the layer at a given position

        Parameters
        ----------
        pos : int
            position in the stack

        Returns
        -------
        L : Layer
            The layer at the position `pos`
        """
        return self.layers[pos]

    def get_superstrate(self):
        """Returns the System's superstrate

        Returns
        -------
        L : Layer
            The system superstrate
        """
        return self.superstrate

    def get_substrate(self):
        """Returns the System's substrate

        Returns
        -------
        L : Layer
            The system substrate
        """
        return self.substrate

    def add_layer(self,layer):
        """Add a layer instance.

        Parameters
        -----------
        layer : Layer
                The layer to be added on the stack

        Returns
        -------
        None

        Notes
        -----
        The layers are added *from superstrate to substrate* order.
        Light is incident *from the superstrate*.

        Note thate this function adds a reference to L to the list.
        If you are adding the same layer several times, be aware that if you
        change something for one of them, it changes all of them.
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
            print('Wrong position given. No layer deleted')

    def initialize_sys(self, f):
        """Sets the values of epsilon at the given frequency in all the layers.

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

    def calculate_incoherent_GammaStar(self,f, zeta_sys):
        """
        Calculate the whole system's transfer matrix.
        But instead of using amplitudes, use intensities
        This routine is specific for a single layer system with sub and superstrate

        Parameters
        -----------
        f : float
            Frequency (Hz)
        zeta_sys : complex
            In-plane wavevector kx/k0

        Returns
        -------
        GammaStar: 4x4 complex matrix
                   System transfer matrix :math:`\Gamma^{*}`
        """

        nan = np.nan_to_num

        Di_super, Pi_super, Di_inv_super, T_super = self.superstrate.update(f, zeta_sys)
        Di_sub, Pi_sub, Di_inv_sub, T_sub = self.substrate.update(f, zeta_sys)

        Delta1234 = np.array([[1,0,0,0],
                              [0,0,1,0],
                              [0,1,0,0],
                              [0,0,0,1]],dtype=np.clongdouble)


        Gamma = np.zeros(4, dtype=np.clongdouble)
        GammaStar = np.zeros(4, dtype=np.clongdouble)
        Dip1 = Di_sub
        Tloc = np.identity(4, dtype=np.clongdouble)
        rescale2 = 1.0
        for layer in reversed(self.layers):
            Di, Pi, Di_inv, Ti = layer.update(f, zeta_sys)
            rescale = np.max(np.abs(Pi))
            rescale2 = rescale*rescale2
            Pi = Pi/rescale
            Tloc = nan(np.matmul( nan(np.absolute( (nan(np.matmul(Di_inv,Dip1))))**2),Tloc))
            Tloc = nan(np.matmul( nan(np.absolute(Pi)**2),Tloc))
            Dip1 = Di
        Tloc = nan(np.matmul( nan(np.absolute(nan(np.matmul(Di_inv_super,Dip1)))**2),Tloc))
        Gamma = nan(np.sqrt(Tloc))
        Gamma = Gamma * rescale2
        GammaStar = np.matmul(exact_inv(Delta1234),np.matmul(Gamma,Delta1234))
        self.Gamma = Gamma.copy()
        self.GammaStar = GammaStar.copy()
        return self.GammaStar.copy()


    def calculate_GammaStar(self,f, zeta_sys):
        """
        Calculate the whole system's transfer matrix.

        Parameters
        -----------
        f : float
            Frequency (Hz)
        zeta_sys : complex
            In-plane wavevector kx/k0

        Returns
        -------
        GammaStar: 4x4 complex matrix
                   System transfer matrix :math:`\Gamma^{*}`
        """
        Ai_super, Ki_super, Ai_inv_super, T_super = self.superstrate.update(f, zeta_sys)
        Ai_sub, Ki_sub, Ai_inv_sub, T_sub = self.substrate.update(f, zeta_sys)

        Delta1234 = np.array([[1,0,0,0],
                              [0,0,1,0],
                              [0,1,0,0],
                              [0,0,0,1]],dtype=np.clongdouble)


        Gamma = np.zeros(4, dtype=np.clongdouble)
        GammaStar = np.zeros(4, dtype=np.clongdouble)
        Tloc = np.identity(4, dtype=np.clongdouble)

        for ii in range(len(self.layers))[::-1]:
            Ai, Ki, Ai_inv, T_ii = self.layers[ii].update(f, zeta_sys)
            Tloc = np.matmul(T_ii,Tloc)

        Gamma = np.matmul(Ai_inv_super,np.matmul(Tloc,Ai_sub))
        GammaStar = np.matmul(exact_inv(Delta1234),np.matmul(Gamma,Delta1234))

        self.Gamma = Gamma.copy()
        self.GammaStar = GammaStar.copy()
        return self.GammaStar.copy()

    def scale_intensity(self, r):
        """r is the complex reflectance field, attempt to make the intensity sensible"""
        if np.abs(r) > 1:
              r = r / np.abs(r) / np.abs(r)
        return r
    
    def calculate_r_t(self, zeta_sys):
        """ Calculate various field and intensity reflection and transmission coefficients, as well as the 4-valued vector of transmitted field.

        Parameters
        -----------
        zeta_sys : complex
            Incident in-plane wavevector
        Returns
        -------
        r_out : len(4)-array
                Complex *field* reflection coefficients r_out=([rpp,rps,rss,rsp])
        R_out : len(4)-array
                Real *intensity* reflection coefficients R_out=([Rpp,Rss,Rsp,Tps])
        t_out : len(4)-array
                Complex *field* transmition coefficients t=([tpp, tps, tsp, tss])
        T_out : len(4)-array
                Real *intensity* transmition coefficients T_out=([Tp,Ts]) (mode-inselective)

        Notes
        -----
        **IMPORTANT**
        ..version 19-03-2020:
        All intensity coefficients are now well defined. Transmission is defined
        mode-independently. It could be defined mode-dependently for non-birefringent
        substrates in future versions.
        The new definition of this function **BREAKS compatibility** with the previous
        one.

        ..version 13-09-2019:
        Note that the field reflectivity and transmission coefficients
        r and t are well defined. The intensity reflection coefficient is also correct.
        However, the intensity transmission coefficients T are ill-defined so far.
        This will be corrected upon future publication of the correct intensity coefficients.

        Note also the different ordering of the coefficients, for consistency w/ Passler's matlab code

        """
        nan = np.nan_to_num
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
        rpp = nan(GammaStar[1,0]*GammaStar[2,2])-nan(GammaStar[1,2]*GammaStar[2,0])
        rpp = np.nan_to_num(rpp/Denom)

        rss = nan(GammaStar[0,0]*GammaStar[3,2])-nan(GammaStar[3,0]*GammaStar[0,2])
        rss = np.nan_to_num(rss/Denom)

        rps = nan(GammaStar[3, 0]*GammaStar[2,2])-nan(GammaStar[3,2]*GammaStar[2,0])
        rps = np.nan_to_num(rps/Denom)

        rsp = nan(GammaStar[0,0]*GammaStar[1,2])-nan(GammaStar[1,0]*GammaStar[0,2])
        rsp = np.nan_to_num(rsp/Denom)

        # Try and remove any outliers by making sure the intensity is always < 1
        rpp = self.scale_intensity(rpp)
        rss = self.scale_intensity(rss)
        rps = self.scale_intensity(rps)
        rsp = self.scale_intensity(rsp)

        # Intensity reflection coefficients are just square moduli
        Rpp = np.abs(rpp)**2
        Rss = np.abs(rss)**2
        Rps = np.abs(rps)**2
        Rsp = np.abs(rsp)**2
        r_out = np.array([rpp,rps,rss,rsp]) ## order matching Passler Matlab code
        R_out = np.array([Rpp,Rss,Rsp,Rps]) ## order matching Passler Matlab code

        # field transmission coefficients
        #t_field = np.zeros(4, dtype=np.complex128)
        t_out = np.zeros(4, dtype=np.clongdouble)
        tpp = np.nan_to_num( (GammaStar[2,2]/Denom)/largest)
        tss = np.nan_to_num( (GammaStar[0,0]/Denom)/largest)
        tps = np.nan_to_num(-(GammaStar[2,0]/Denom)/largest)
        tsp = np.nan_to_num(-(GammaStar[0,2]/Denom)/largest)

        # Try and remove any outliers by making sure the intensity is always < 1
        tpp = self.scale_intensity(tpp)
        tss = self.scale_intensity(tss)
        tps = self.scale_intensity(tps)
        tsp = self.scale_intensity(tsp)
        t_out = np.array([tpp, tps, tsp, tss],dtype=np.clongdouble)
        #t_field = np.array([tpp, tps, tsp, tss])

        #### Intensity transmission requires Poynting vector analysis
        ## N.B: could be done mode-dependentely later
        ## start with the superstrate
        ## Incident fields are either p or s polarized
        ksup = np.zeros((4,3), dtype=np.clongdouble) ## wavevector in superstrate
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
        ksub = np.zeros((4,3), dtype=np.clongdouble)
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
        T_pp = self.scale_intensity(T_pp)
        T_ss = self.scale_intensity(T_ss)

        T_out = np.array([T_pp, T_ss])

        return r_out, R_out, t_out, T_out

