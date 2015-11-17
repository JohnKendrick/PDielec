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
"""Define some useful constants"""
import math 

#define some useful basic constants
PI = 2.0 * math.acos( 0.0 )
inv4pi = 1.0 / ( 4.0 * PI )
speed_light_si = 299792458.0
mu_0_si = 4.0*PI*1e-7
epsilon_0_si = 1.0/(mu_0_si*speed_light_si*speed_light_si)
planck_si = 6.62606896e-34
elementary_charge_si = 1.602176487e-19
electron_mass_si = 9.10938215e-31
proton_mass_si = 1.672621637e-27
electron_gyromagnetic_ratio_si = 1.76085977e11
avogadro_si = 6.02214179e23
molar_gas_si = 8.314472
hbar_si = planck_si/(2.0*PI)
fine_structure_si = elementary_charge_si*elementary_charge_si / ( 4*PI*epsilon_0_si*hbar_si*speed_light_si)
boltzmann_si = molar_gas_si/avogadro_si
amu_si = 1.0E-3 / avogadro_si
# Derived units
coulomb = 1.0 / elementary_charge_si
metre =  electron_mass_si*speed_light_si*fine_structure_si/hbar_si
angstrom = 1.E-10 * metre
amu   =    1.E-3 / avogadro_si / electron_mass_si
debye_si = 1.0E-21 / speed_light_si
debye = debye_si * coulomb * metre
d2byamuang2 = debye*debye / (angstrom * angstrom * amu)
joule = 1.0 / ( fine_structure_si*fine_structure_si * electron_mass_si * speed_light_si*speed_light_si)
hertz = planck_si*joule
wavenumber = hertz*speed_light_si*100.0
strength_conv  =  1.0 / d2byamuang2
