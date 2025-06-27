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
"""The SuperCell module."""

import numpy as np

from PDielec.Plotter import print_ints, print_reals


class SuperCell:
    """Hold super-cell information; the unit cell and the super-cell size.

    Parameters
    ----------
    aUnitCell : :class:`~PDielec.UnitCell.UnitCell`
        The unit cell associated with this instance.
    anImageSpecifier : list or tuple
        The list or tuples of 3 numbers used to generate images.

    Attributes
    ----------
    unitCell : :class:`~PDielec.UnitCell.UnitCell`
        The unit cell associated with this instance. (see :class:`~PDielec.UnitCell.UnitCell`)
    imageList : NoneType or list
        Initialized as None, can later hold a list of images of the cell
    imageSpecifier : list or tuple
        A list or tuple of 3 numbers which indicate the number of images of the cell to be created
    XYZ : NoneType or array_like
        Initialized as None, can later hold XYZ coordinates. (Cartesian coordinates)

    Notes
    -----
    This method automatically invokes the createImages method at the end of initialization.

    See Also
    --------
    :class:`~PDielec.UnitCell.UnitCell`

    """

    def __init__(self, aUnitCell,anImageSpecifier):
        """Initialize the supercell with a unit cell and an image specifier.

        Parameters
        ----------
        aUnitCell : :class:`~PDielec.UnitCell.UnitCell`
            The unit cell associated with this instance.
        anImageSpecifier : list or tuple
            The list or tuples of 3 numbers used to generate images.

        Attributes
        ----------
        unitCell : :class:`~PDielec.UnitCell.UnitCell`
            The unit cell associated with this instance.
        imageList : NoneType or list
            Initialized as None, can later hold a list of images of the cell
        imageSpecifier : list or tuple
            A list or tuple of 3 numbers which indicate the number of images of the cell to be created
        XYZ : NoneType or array_like
            Initialized as None, can later hold XYZ coordinates. (Cartesian coordinates)

        Notes
        -----
        This method automatically invokes the createImages method at the end of initialization.

        """
        self.unitCell       = aUnitCell         # the unit cell on which everything is based
        self.imageList      = None              # a list of the images of the unit cell
        self.imageSpecifier = anImageSpecifier  # the specification used to create the imageList
        self.XYZ            = None              # Cartesian coordinates of the super-cell
        self.bonds          = []                # A list of bonds in the super-cell
        self.createImages(anImageSpecifier)
        self.calculateXYZ()

    def print(self):
        """Print the contents of the SuperCell.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Example
        --------
        ```
        # Create a 3x2x1 super-cell
        import PDielec.SuperCell as SC
        import PDielec.UnitCell as UC
        uc=UC.UnitCell(1.0,2.0,3.0,90.0,90.0,90.0)
        imageSpecifier = [3,2,1]
        sc=SC.SuperCell(uc,imageSpecifier)
        sc.print()
        ```

        """
        print("SuperCell")
        print_ints("Image_specifier",self.imageSpecifier)
        self.unitCell.print()
        print("SuperCell Image list")
        for l in self.imageList:
            print_ints("",l)
        corners,edges,labels = self.getBoundingBox()
        print("SuperCell Corners")
        for corner in corners:
            print_reals("",corner,format="{:12.6f}")
        if self.XYZ is not None:
            print()
            print("SuperCell XYZ (Coordinates)")
            for xyz in self.XYZ:
                print_reals("",xyz,format="{:12.6f}")
        if len(self.bonds) > 0:
            print()
            print("SuperCell Bond list")
            for bond in self.bonds:
                print_ints("",bond)

    def createImages(self,anImageSpecifier):
        """Create a list of images extending by ia, jb, kc units in the a, b, c directions. The image specifier is a list of the form [ia, jb, kc].

        Parameters
        ----------
        anImageSpecifier : list (ia, jb, kc)
            ia an integer specifying extension in the a direction.
            jb an integer specifying extension in the b direction.
            kc an integer specifying extension in the c direction.

        Returns
        -------
        list
            A list of created images.

        """
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

    def getElementNames(self):
        """Calculate and return the element names of the supercell.

        Parameters
        ----------
        None

        Returns
        -------
        list
            The element names of the supercell.

        """
        names = []
        cell_names = self.unitCell.element_names
        for  _i,_j,_k in self.imageList:
            names.extend( cell_names )
        #end for i,j,k
        return names

    def calculateNormalModes(self,modes):
        """Calculate and return the normal modes of the supercell.

        Normal modes are initially an np array of nmodes, 3*nions.
        The number of modes stays the same but the number of ions
        increases according to the number of images.

        Parameters
        ----------
        modes : a list of floats
            A list of the modes to be used in the calculation

        Returns
        -------
        np.ndarray
            Array representing the normal modes of the supercell. Initially, this is an array
            of shape (nmodes, 3*nions), where `nmodes` is the number of modes and `nions` is the
            number of ions per mode. The shape may vary as the number of ions increases with the
            number of images.

        """
        nmodes,nions,_ = np.shape(modes)
        nions3 = nions*3
        # Make a copy of the modes, transforming [nmodes, nions, 3] to [ nmodes, nions3 ]
        modes_copy = np.zeros( (nmodes, nions3) )
        for imode,mode in enumerate(modes):
            pos = 0
            modes_copy[imode] = mode.flatten()
        # Now make multiple copies of the normal modes for each image.
        mult = len(self.imageList)
        normal_modes = np.zeros( (nmodes, mult*nions3) )
        for mode_index,mode in enumerate(modes_copy):
            pos = 0
            for  _i,_j,_k in self.imageList:
                for m in mode:
                  normal_modes[mode_index,pos] = m
                  pos += 1
                #end for m in mode
            #end for i,j,k
        #end for i,mode 
        return normal_modes

    def calculateXYZ(self):
        """Calculate and return the cartesian coordinates of the supercell.

        The fractional cell coordinates of the unit cell are used to calculate the XYZ coordinates of the super-cell

        Parameters
        ----------
        None

        Returns
        -------
        output : array
            The Cartesian coordinates of the supercell.

        """
        self.XYZ = []
        coords = []
        for  i,j,k in self.imageList:
            for fx,fy,fz in self.unitCell.fractional_coordinates:
                fx += i
                fy += j
                fz += k
                coords.append( (fx,fy,fz) )
            # end for fx,fy,fz
        #end for i,j,k
        self.XYZ = self.unitCell.convert_abc_to_xyz(coords)
        return self.XYZ

    def calculateBonds(self):
        """Calculate the bonding in the supercell.

        Parameters
        ----------
        None

        Returns
        -------
        list
            A list of the bonds as pairs of atom indexes

        Notes
        -----
        The bonding in the unit-cell is assumed to be valid whole-molecule bonding
        the bonds are not recalculated, just extended into the super-cell

        """
        self.bonds = []
        natoms = self.unitCell.nions
        bonds = self.unitCell.bonds
        for cell in range(len(self.imageList)):
            for na,nb in bonds:
              na += cell*natoms
              nb += cell*natoms
              self.bonds.append( (na,nb) )
            # end for na,nb
        #end for i,j,k
        return self.bonds


    def getBoundingBox(self,originABC=None):
        """Calculate the bounding box of an object in Cartesian coordinates based on its image specifications and an origin.

        Return a box with 8 corners and 12 edges which represent the supercel in cartesian space
        - Corners in abc space
        - shift the origin to originABC

        Parameters
        ----------
        originABC : array-like, optional
            The origin coordinates (A, B, C) to which the bounding box will be translated.
            defaults to 0,0,0

        Returns
        -------
        tuple
            A tuple containing three elements:
    
            - corners_xyz : list
              A list of corner points in XYZ coordinates specifying the bounding box, shifted to the provided origin.
    
            - edges : list
              A list of tuples, each representing an edge of the bounding box defined by a pair of XYZ coordinates.

            - labels : list of strings
              A list of strings labelling the bounding box a,b,c,o


        Notes
        -----
        The function calculates the center of the bounding box, re-centers the corners to the provided origin, and then builds edges between these corners in XYZ space.

        """
        if originABC is None:
            originABC = [0.0, 0.0, 0.0]
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
        labels = [ "o", "a", "", "b", "c", "", "", "" ]
        centre = np.array ( [0.0, 0.0, 0.0] )
        for corner in corners:
            centre += corner
        centre = centre / 8.0
        corners = [ corner - centre + originABC for corner  in corners ]
        # Now calculate the corners in cartesian space
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
        return corners_xyz,edges,labels

    def calculateCentreOfBox(self):
        """Calculate the centre of the super-cell in cartesian coordinates.

        This is the centre of the box, not the centre of mass

        Parameters
        ----------
        None

        Returns
        -------
        vector
            The centre of the images in the box

        """
        i,j,k = self.imageSpecifier
        i = 0.5*float(i)
        j = 0.5*float(j)
        k = 0.5*float(k)
        ijk = [ i,j,k ]
        return self.unitCell.convert_abc_to_xyz(ijk)


    def calculateCentreOfMass(self,output="xyz"):
        """Calculate the center of mass for a molecular structure.

        Parameters
        ----------
        output : str, optional
            The output. Options are 'xyz' for mass and Cartesian coordinates,
            'mass' for the total mass only, 'abc' for mass and fractional coordinates, or any other string
            to return a tuple with all three values (mass, Cartesian coordinates, fractional coordinates).
            The default value is 'xyz'.

        Returns
        -------
        numpy.ndarray or float or tuple
            Depending on the value of `units`:
            If 'xyz', returns a tuple with mass and a numpy.ndarray of Cartesian coordinates of the center of mass.
            If 'mass', returns the total mass as a float.
            If 'abc', returns tuple with mass and a numpy.ndarray of fractional coordinates of the center of mass.
            For any other value of `output`, returns a tuple containing the total mass (float), Cartesian
            coordinates (numpy.ndarray), and fractional coordinates (numpy.ndarray) of the center of mass.

        """
        # Calculate the centre of mass 
        # The centre of mass can be returned as 'xyz' space or 'abc' space
        # if output='all' a tuple of (mass,cm_xyz,cm_abc) is returned
        # if output='mass' a mass is returned
        mass = 0.0
        cm_fractional = np.zeros(3)
        unit_mass, unit_cm_xyz, unit_cm_fractional = self.unitCell.calculateCentreOfMass(output=all)
        for  i,j,k in self.imageList:
            cm_fractional[0] += i + unit_cm_fractional [0]
            cm_fractional[1] += j + unit_cm_fractional [1]
            cm_fractional[2] += k + unit_cm_fractional [2]
            mass += unit_mass
        #end for i,j,k
        cm_fractional = cm_fractional/len(self.imageList)
        cm_xyz = self.unitCell.convert_abc_to_xyz(cm_fractional)
        return {
            "xyz"  : cm_xyz,
            "mass" : mass,
            "abc"  : cm_fractional,
            "all"  : (mass, cm_xyz, cm_fractional ),
        }.get(output, (mass, cm_xyz, cm_fractional))
