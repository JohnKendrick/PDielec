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

"""Hold unit cell information and its associated calculated properties"""

import numpy as np


class UnitCell:
    """Hold unit cell information and its associated calculated properties"""
    def __init__(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None):
        self.fractional_coordinates = []
        self.xyz_coordinates = []
        self.lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        if gamma is not None:
            self.lattice = self.convert_abc_to_unitcell(a, b, c, alpha, beta, gamma)
        else:
            self.lattice[0] = a
            self.lattice[1] = b
            self.lattice[2] = c
        self._calculate_reciprocal_lattice(self.lattice)

    def convert_unitcell_to_abc(self):
        """Convert a unit cell to the equivalent a, b, c, alpha, beta, gamma designation"""
        a = np.sqrt(np.dot(self.lattice[0], self.lattice[0]))
        b = np.sqrt(np.dot(self.lattice[1], self.lattice[1]))
        c = np.sqrt(np.dot(self.lattice[2], self.lattice[2]))
        gamma = np.arccos(np.dot(self.lattice[0], self.lattice[1]) / (a*b))
        beta = np.arccos(np.dot(self.lattice[0], self.lattice[2]) / (a*c))
        alpha = np.arccos(np.dot(self.lattice[1], self.lattice[2]) / (b*c))
        return a, b, c, np.degrees(alpha), np.degrees(beta), np.degrees(gamma)

    def convert_abc_to_unitcell(self, a, b, c, alpha_degs, beta_degs, gamma_degs):
        """Convert a,b,c,alpha,beta,gamma to a unit cell"""
        # This is castep convention, need to check it works with vasp
        alpha = np.radians(alpha_degs)
        beta = np.radians(beta_degs)
        gamma = np.radians(gamma_degs)
        lattice = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        lattice[2] = [0.0, 0.0, c]
        lattice[1] = [0.0, b*np.sin(np.radians(alpha)), b*np.cos(alpha)]
        z = a * np.cos(beta)
        y = a * (np.cos(gamma) - np.cos(alpha)*np.cos(beta))/np.sin(alpha)
        x = a * np.sqrt(1.0 - np.cos(beta)**2 - (y/a)**2)
        lattice[0] = [x, y, z]
        self._calculate_reciprocal_lattice(lattice)
        return self.lattice

    def _calculate_reciprocal_lattice(self, lattice):
        """Calculate the reciprocal lattice"""
        self.lattice = np.array(lattice)
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = self.convert_unitcell_to_abc()
        self.volume = self.calculate_volume()
        self.inverse_lattice = np.linalg.inv(self.lattice)
        self.reciprocal_lattice = self.inverse_lattice

    def calculate_volume(self):
        """Calculate the volume"""
        volume = np.abs(np.dot(self.lattice[0], np.cross(self.lattice[1], self.lattice[2])))
        return volume

    def convert_xyz_to_abc(self, xyz):
        """Convert xyz coordinates to abc lattice coordinates"""
        xyz = np.array(xyz)
        abc = np.dot(xyz, self.reciprocal_lattice)
        return abc

    def convert_hkl_to_xyz(self, hkl):
        """Convert hkl coordinates to xyz coordinates"""
        hkl = np.array(hkl)
        xyz = np.dot(hkl, self.reciprocal_lattice.T)
        return xyz

    def convert_abc_to_xyz(self, abc):
        """Convert abc coordinates to xyz coordinates"""
        abc = np.array(abc)
        xyz = np.dot(abc, self.lattice)
        return xyz

    def convert_hkl_to_xyz2(self, hkl):
        """Convert hkl coordinates to xyz coordinates using a different method"""
        # Use the hkl miller indices to calculate the normal to a plane in cartesian space
        p = []
        if not hkl[0] == 0:
            p.append(np.array(self.lattice[0] / hkl[0]))
        if not hkl[1] == 0:
            p.append(np.array(self.lattice[1] / hkl[1]))
        if not hkl[2] == 0:
            p.append(np.array(self.lattice[2] / hkl[2]))
        if hkl[0] == 0:
            p.append(self.lattice[0] + p[0])
        if hkl[1] == 0:
            p.append(self.lattice[1] + p[0])
        if hkl[2] == 0:
            p.append(self.lattice[2] + p[0])
        p12 = p[1] - p[0]
        p13 = p[2] - p[0]
        normal = np.cross(p12, p13)
        norm = np.linalg.norm(normal)
        if norm < 1.0e-8:
            print("Error in unit cell, calculation of normal to miller index failed")
            exit(1)
        normal = normal / norm
        return normal

    def set_fractional_coordinates(self, coords):
        """Set the fractional coordinates and calculate the xyz coordinates"""
        self.fractional_coordinates = coords
        self.xyz_coordinates = []
        for abc in coords:
            xyz = self.convert_abc_to_xyz(abc)
            self.xyz_coordinates.append(xyz)
        return

    def set_xyz_coordinates(self, coords):
        """Set the xyz coordinates and calculate the fractional coordinates"""
        self.xyz_coordinates = coords
        self.fractional_coordinates = []
        for xyz in coords:
            abc = self.convert_xyz_to_abc(xyz)
            self.fractional_coordinates.append(abc)
        return
