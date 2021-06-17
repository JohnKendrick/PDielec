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

from __future__ import print_function
from collections import deque
import numpy as np


class SuperCell:
    """Hold super-cell information; the unit cell and the super-cell size"""
    def __init__(self, aUnitCell,anImageSpecifier):
        self.unitCell       = aUnitCell         # the unit cell on which everything is based
        self.imageList      = None              # a list of the images of the unit cell
        self.imageSpecifier = anImageSpecifier  # the specification used to create the imageList
        self.XYZ            = None              # Cartesian coordinates
        self.UVW            = None              # Displacements
        self.createImages(anImageSpecifier)

    def createImages(self,anImageSpecifier):
        """Create a list of images extending by ia, jb, kc units  in the a, b, c directions.
           The image specifier is a list of the form [ia, jb, kc]"""
        ia,jb,kc = anImageSpecifier
        if ia < 1 or jb < 1 or kc < 1 :
            print("Error in createSuperCell: parameters must be 1 or more")
            return
        self.imageSpecifier = anImageSpecifier
        self.imageList = []
        for i in range(ia):
            for j in range(jb):
                for k in range(kc):
                    self.imageList.append( (i,j,k) )
        #jk print('CreateImages: ',self.imageList)

    def getElementNames(self):
        """calculate and return the element names of the supercell"""
        names = []
        cell_names = self.unitCell.element_names
        for  i,j,k in self.imageList:
            names.extend( cell_names )
        #end for i,j,k
        #jk print('getElementNames: ',names)
        return names

    def calculateNormalModes(self,modes):
        """Calculate and return the normal modes of the supercell
        Normal modes are initially an np array of nmodes, 3*nions
        The number of modes stays the same but the number of ions
        increases according to the number of images"""
        nmodes,nions3 = np.shape(modes)
        mult = len(self.imageList)
        normal_modes = np.zeros( (nmodes, mult*nions3) )
        for mode_index,mode in enumerate(modes):
            pos = 0
            for  i,j,k in self.imageList:
                for m in mode:
                  normal_modes[mode_index,pos] = m
                  pos += 1
                #end for m in mode
            #end for i,j,k
        #end for i,mode 
        #jk print('calculateNormalModes: ',modes[0])
        #jk print('calculateNormalModes: ',normal_modes[0])
        return normal_modes

    def calculateXYZ(self):
        """Calculate and return the cartesian cordinates of the supercell"""
        self.XYZ = []
        coords = []
        for  i,j,k in self.imageList:
            #jk print('calculateXYZ:', i,j,k)
            for fx,fy,fz in self.unitCell.fractional_coordinates:
                fx += i
                fy += j
                fz += k
                #jk print('calculateXYZ:', fx,fy,fz)
                coords.append( (fx,fy,fz) )
            # end for fx,fy,fz
        #end for i,j,k
        self.XYZ = self.unitCell.convert_abc_to_xyz(coords)
        #jk print('calculateXYZ: ',coords)
        #jk print('calculateXYZ: ',self.XYZ)
        return self.XYZ

    def calculateBonds(self):
        """Calculate the bonding in the supercell """
        self.bonds = []
        natoms = self.unitCell.nions
        bonds = self.unitCell.bonds
        cell = 0
        for i,j,k in self.imageList:
            for na,nb in bonds:
              na += cell*natoms
              nb += cell*natoms
              self.bonds.append( (na,nb) )
            # end for na,nb
            cell += 1
        #end for i,j,k
        #jk print('calculateBonds: ',self.bonds)
        return self.bonds


    def getBoundingBox(self,originABC):
        # Return a box with 12 edges which represent the supercel in cartesian space
        #  7---6
        # /|  /|
        # 3---2|
        # |4--|5
        # 0---1/
        # Corners in abc space
        # shift the origin to originABC
        i,j,k = self.imageSpecifier
        i = float(i)
        j = float(j)
        k = float(k)
        corners = []
        corners.append( np.array( [0.0,0.0,0.0] ) )
        corners.append( np.array( [  i,0.0,0.0] ) )
        corners.append( np.array( [  i,  j,0.0] ) )
        corners.append( np.array( [0.0,  j,0.0] ) )
        corners.append( np.array( [0.0,0.0,  k] ) )
        corners.append( np.array( [  i,0.0,  k] ) )
        corners.append( np.array( [  i,  j,  k] ) )
        corners.append( np.array( [0.0,  j,  k] ) )
        centre = np.array ( [0.0, 0.0, 0.0] )
        for corner in corners:
            centre += corner
        centre = centre / 8.0
        corners = [ corner - centre + originABC for corner  in corners ]
        # Now calculate the corners in cartesian space
        #jk corners_xyz = [ self.UnitCell.convert_abc_to_xyz(corner) for corner in corners ]
        corners_xyz = self.unitCell.convert_abc_to_xyz(corners)
        # Now calculate the edges each edge is a tuple with a beginning coordinate and an end coordinate
        edges = []
        edges.append( (corners_xyz[1] , corners_xyz[0]) )
        edges.append( (corners_xyz[2] , corners_xyz[1]) )
        edges.append( (corners_xyz[3] , corners_xyz[2]) )
        edges.append( (corners_xyz[3] , corners_xyz[0]) )
        edges.append( (corners_xyz[4] , corners_xyz[0]) )
        edges.append( (corners_xyz[5] , corners_xyz[1]) )
        edges.append( (corners_xyz[6] , corners_xyz[2]) )
        edges.append( (corners_xyz[7] , corners_xyz[3]) )
        edges.append( (corners_xyz[5] , corners_xyz[4]) )
        edges.append( (corners_xyz[6] , corners_xyz[5]) )
        edges.append( (corners_xyz[7] , corners_xyz[6]) )
        edges.append( (corners_xyz[7] , corners_xyz[4]) )
        return corners_xyz,edges

    def calculateCentreOfBox(self):
        """Not strictly the centre of mass but it will do.  We calculate the centre of the box in cartesian coordinates"""
        i,j,k = self.imageSpecifier
        i = 0.5*float(i)
        j = 0.5*float(j)
        k = 0.5*float(k)
        ijk = [ i,j,k ]
        centre = self.unitCell.convert_abc_to_xyz(ijk)
        return centre


    def calculateCentreOfMass(self,units='xyz'):
        # Calculate the centre of mass 
        # The centre of mass can be returned in units of 'xyz' space or 'abc' space
        # if units='all' a tuple of (mass,cm_xyz,cm_abc) is returned
        mass = 0.0
        cm_fractional = np.zeros(3)
        unit_mass, unit_cm_xyz, unit_cm_fractional = self.unitCell.calculateCentreOfMass(units=all)
        for  i,j,k in self.imageList:
            cm_fractional[0] += i + unit_cm_fractional [0]
            cm_fractional[1] += j + unit_cm_fractional [1]
            cm_fractional[2] += k + unit_cm_fractional [2]
            mass += unit_mass
        #end for i,j,k
        cm_fractional = cm_fractional/len(self.imageList)
        cm_xyz = self.unitCell.convert_abc_to_xyz(cm_fractional)
        if units == 'xyz':
            return cm_xyz
        elif units == 'mass':
            return mass
        elif units == 'abc':
            return cm_fractional
        else:
            return mass, cm_xyz, cm_fractional


