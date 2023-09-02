# This file is part of the pyGTM module.
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
Example file for the pyGTM package

It reproduces the main result of section IV in Passler, Jeannin and Paarman
https://arxiv.org/abs/2002.03832
fig 2

It also demonstrates how to use the euler angles for the layers
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import GTMcore as GTM
import Permittivities as mat
from matplotlib import rcParams, rcParamsDefault

c_const = 299792458.0 # m/s

# setup Simulation

## Parameters
f_cm = np.linspace(100.0, 400.0, 100) ## frequency range (cm-1->Hz)
f_sp = 100.0*c_const*f_cm ## frequency range (cm-1->Hz)
## **careful** angle of incidence is theta_in, 
## it is **not** the euler anle theta for the layers
theta_in = np.deg2rad(80.0)

## Note that Phi in the text is actually Euler angle psi
psi_v = np.deg2rad(0.0) # to get the last value
phi = 0.0
# the principal vertical axis of the layer is still z so theta=0
theta = 0.0

## physical dimensions
region1_thickness = 1000e-6
region2_thickness = 1000e-6
crystal_thickness = 10.0e-6
crystal_dielectric_function = lambda x: np.array([ [2.0+0.1j, 0.0, 0.0], [0.0, 2.0+0.1j, 0.0], [0.0, 0.0, 2.0+0.1j] ])
region1_dielectric_function = lambda x: np.array([ [1.1, 0.0, 0.0], [0.0, 1.1, 0.0], [0.0, 0.0, 1.1] ])
region2_dielectric_function = lambda x: np.array([ [1.5, 0.0, 0.0], [0.0, 1.5, 0.0], [0.0, 0.0, 1.5] ])

## substrate, superstrate and layers, angles are left for later
region1  = GTM.Layer(thickness=region1_thickness, epsilon=region1_dielectric_function)
crystal  = GTM.Layer(thickness=crystal_thickness, epsilon=crystal_dielectric_function)
region2  = GTM.Layer(thickness=region2_thickness, epsilon=region2_dielectric_function)
## setup the system
S = GTM.System()
S.set_superstrate(region1)
S.set_substrate(region2)
S.add_layer(crystal)


#
R_s = np.zeros((len(f_sp)))
R_p = np.zeros((len(f_sp)))
T_s = np.zeros((len(f_sp)))
T_p = np.zeros((len(f_sp)))

## set the layers orientation
S.substrate.set_euler(theta=theta, phi=phi, psi=psi_v)
S.superstrate.set_euler(theta=theta, phi=phi, psi=0.0) # this guy stands still
for L in S.layers:
    L.set_euler(theta=theta, phi=phi, psi=psi_v)
# loop for frequency
for ii, fi in enumerate(f_sp):
    S.initialize_sys(fi) # sets the epsilons w/ correct euler rotation
    zeta_sys = np.sin(theta_in)*np.sqrt(S.superstrate.epsilon[0,0]) # in-plane wavevector
    S.calculate_GammaStar(fi, zeta_sys) # main computation (discard output)
    r, R_loc, t, T_loc = S.calculate_r_t(zeta_sys) # calculate reflectivity
    R_p[ii] = R_loc[0] # p-pol only
    R_s[ii] = R_loc[1] # s-pol only
    T_p[ii] = T_loc[0] # p-pol only
    T_s[ii] = T_loc[1] # s-pol only
    zplot, E_out, H_out, zn_plot = S.calculate_Efield(fi, zeta_sys, magnetic=True)
    zplot, E_out, H_out, zn_plot = S.calculate_Efield(fi, zeta_sys, magnetic=True)
    S_loc, A_loc = S.calculate_Poynting_Absorption_vs_z(zplot, E_out, H_out, R_loc)

# Make pretty plot
for v, rp,rs,tp,ts in zip(f_cm,R_p,R_s,T_p,T_s):
    print(v,rs,rp,ts,tp)

