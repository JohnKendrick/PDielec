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
"""Hold unit cell information and its associated calculated properties."""

import math
import sys
from contextlib import nullcontext

import numpy as np
import spglib

from PDielec.Calculator import (
    calculate_angle,
    calculate_distance,
    calculate_torsion,
    cleanup_symbol,
)
from PDielec.Constants import atomic_number_to_element, avogadro_si, covalent_radii, element_to_atomic_number
from PDielec.Plotter import print_ints, print_reals, print_strings


def convert_length_units(value, units_in, units_out):
    """"Convert between different length units.

    The 'internal' unit is taken to be the Angstrom so units are relative to the Angstrom
    The unit strings are made lowercase, so case should be irrelevant

    Parameters
    ----------
    value : float or list of floats or a numpy array
        The value(s) for which the conversion is to be made.
    units_in : str
        The units of the input value(s). Can be one of '
    units_out : str
        The units of the output value(s). Must be one of 'a.u. au bohr ang angs angstrom nm um mm cm m'

    Returns
    -------
    scalar or numpy array
        The converted value(s) in the output units specified.

    Notes
    -----
    The input can be either a scalar value, a list or a numpy array of values. The function will return the converted value(s) in the output units specified.

    """
    # the conversion dictionary has a value that converts the key unit to angstroms
    angstroms = { "a.u."       : 0.5291772, 
                  "au"         : 0.5291772, 
                  "bohr"       : 0.5291772, 
                  "ang"        : 1.0000000,
                  "angs"       : 1.0000000,
                  "angstrom"   : 1.0000000,
                  "angstroms"  : 1.0000000,
                  "nm"         : 1.0E1    ,
                  "um"         : 1.0E4    ,
                  "mm"         : 1.0E7    ,
                  "cm"         : 1.0E8    ,
                  "m"          : 1.0E10   ,
                 }
    units_in  = units_in.lower()
    units_out = units_out.lower()
    if units_in == units_out:
        return value
    # Deal with a possible list
    wasList = False
    if isinstance(value, list):
        wasList = True
        value = np.array(value)
    # Convert the input unit angstrom
    scale = angstroms[units_in]
    value = scale * value     
    # convert the internal value from angstrom to the output unit
    scale = angstroms[units_out]
    value = value / scale
    if wasList:
        value = value.tolist()
    return value


class UnitCell:
    """Hold unit cell information and its associated calculated properties.

    Initialize the class instance with optional lattice parameters and calculate the reciprocal lattice.

    Parameters
    ----------
    a, b, c : float or array, optional
        Lattice vectors or cell lengths. If not specified, they default to [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], and [0.0, 0.0, 1.0], respectively.
    alpha, beta, gamma : float, optional
        Lattice angles (in degrees). These are only used if all three angles are specified, otherwise, the default lattice (orthorhombic) is used.
    units : str
        An optional unit such as 'a.u., au bohr angs angstrom Angs Angstrom or nm'  The default is Angstrom.
        The internal unit is always the angstrom

    Notes
    -----
    This constructor initializes an instance with empty lists for fractional coordinates, xyz coordinates, element names, atom labels, bonds, and molecules. It also initializes an empty list for atomic masses and zero for total mass. If the angles alpha, beta, and gamma are provided, it attempts to convert the provided lattice parameters (a, b, c, alpha, beta, gamma) into a 3x3 lattice matrix. If not, it directly assigns a, b, and c as lattice vectors. Finally, it calculates and sets the reciprocal lattice for the instance.

    Examples
    --------
    ::
    
        a = [2.853604, -1.647529, 0.0]
        b = [0.0,       3.295058, 0.0]
        c = [0.0,       0.0,      5.284824]
        cell = UnitCell( a, b, c )
        cell.set_element_names(['Zn', 'Zn', 'O', 'O'])
        coords = [ [ 0.333333,    0.666667,    0.000900 ]
                   [ 0.666667,    0.333333,    0.500900 ]
                   [ 0.333333,    0.666667,    0.381600 ]
                   [ 0.666667,    0.333333,    0.881600 ] ]
        cell.set_fractional_coordinates(coords)
        cell.print()

    """

    def __init__(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None, units="Angstrom"):
        """Initialize the class instance with optional lattice parameters and calculate the reciprocal lattice.

        Parameters
        ----------
        a, b, c : float or array, optional
            Lattice vectors or cell lengths. If not specified, they default to [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], and [0.0, 0.0, 1.0], respectively.
        alpha, beta, gamma : float, optional
            Lattice angles (in degrees). These are only used if all three angles are specified, otherwise, the default lattice (orthorhombic) is used.
        units : str
            An optional unit such as 'a.u., au bohr angs angstrom Angs Angstrom or nm'  The default is Angstrom.
            The internal unit is always the angstrom

        Notes
        -----
        This constructor initializes an instance with empty lists for fractional coordinates, xyz coordinates, element names, atom labels, bonds, and molecules. It also initializes an empty list for atomic masses and zero for total mass. If the angles alpha, beta, and gamma are provided, it attempts to convert the provided lattice parameters (a, b, c, alpha, beta, gamma) into a 3x3 lattice matrix. If not, it directly assigns a, b, and c as lattice vectors. Finally, it calculates and sets the reciprocal lattice for the instance.

        Examples
        --------
        ::
        
            a = [2.853604, -1.647529, 0.0]
            b = [0.0,       3.295058, 0.0]
            c = [0.0,       0.0,      5.284824]
            cell = UnitCell( a, b, c )
            cell.set_element_names(['Zn', 'Zn', 'O', 'O'])
            coords = [ [ 0.333333,    0.666667,    0.000900 ]
                       [ 0.666667,    0.333333,    0.500900 ]
                       [ 0.333333,    0.666667,    0.381600 ]
                       [ 0.666667,    0.333333,    0.881600 ] ]
            cell.set_fractional_coordinates(coords)
            cell.print()

        """        
        self.fractional_coordinates = []
        self.xyz_coordinates = []
        self.element_names = []
        self.atom_labels = []
        self.bonds = []
        self.nions = 0
        self.molecules = []
        self.atomic_masses = []
        self.centres_of_mass = []
        self.total_mass = 0.0
        self.lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        # The internal units will be angstrom
        self.units = "Angstrom"
        if  None not in [a, b, c, alpha, beta, gamma]:
            abc = convert_length_units( [a, b, c],units,self.units )
            a, b, c = tuple(abc)
            self.lattice = self.convert_abc_to_unitcell(a, b, c, alpha, beta, gamma)
        elif None not in [a, b, c]:
            abc = convert_length_units( [a, b, c],units,self.units )
            a, b, c = tuple(abc)
            self.lattice[0] = a
            self.lattice[1] = b
            self.lattice[2] = c
        else:
            self.lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        # convert lattice to an np array
        self.lattice = np.array(self.lattice)
        self._calculate_reciprocal_lattice()

    def set_lattice(self, lattice, units="Angstrom"):
        """Set the lattice parameters with a 3x3 matrix.

        The lattice is stored as a numpy array
        The reciprical vectors are calculated.

        Parameters
        ----------
        lattice : 3x3 list
            The lattice vectors stored a ( a b c )
            where a, b, c are 3 vectors in xyz space.
        units : str
            A string indicating the units used in lattice.
            Internally the units are Angstrom

        Returns
        -------
        None

        """
        lattice = convert_length_units(lattice,units,self.units)
        self.lattice = np.array(lattice)
        self._calculate_reciprocal_lattice()
        return

    def write_cif(self, filename=None, filedescriptor=sys.stdout, description=None):
        """Write the crystallographic information file (CIF) representation of a structure.

        Parameters
        ----------
        filename : str, optional
            The name of the file to be written. If not provided, the CIF data is printed to stdout.
        filedescriptor : filedescriptor, optional
            The file descriptor of the file to be written. If not provided, the CIF data is printed to stdout.
        description : str, optional
            A description of the cif file

        Returns
        -------
        None

        Notes
        -----
        This method prints the CIF representation of a structure, including the space group,
        cell dimensions, angles, volume, and the fractional coordinates of the atoms.
        Any lengths are converted to Angstrom, the volumes is give as Angstrom^3
        If a filename is provided, the CIF data will be written to that file. Otherwise,
        it will be printed to standard output.
        Only one of filename or filedescriptor can be specified

        Examples
        --------
        To print the CIF representation to standard output, simply call the method without arguments:

        >>> cell.write_cif()

        To write the CIF representation to a file named "example.cif":

        >>> cell.write_cif('example.cif')

        """        
        abc = [self.a, self.b, self.c]
        volume = self.getVolume("Angstrom")
        spg_symbol, spg_number = self.find_symmetry()
        # Open the filename if it is given
        with open(filename,"w") if filename else nullcontext(filedescriptor) as file_:
            if description is not None:
                print("data_"+description,file=file_)
            else:
                print("data_",         file=file_)
            print(f"_symmetry_space_group_name_H-M '{spg_symbol}'",file=file_)
            print(f"_symmetry_Int_Tables_number      {spg_number}  ",file=file_)
            print(f"_cell_length_a      {abc[0]:12.6f}",     file=file_)
            print(f"_cell_length_b      {abc[1]:12.6f}",     file=file_)
            print(f"_cell_length_c      {abc[2]:12.6f}",     file=file_)
            print(f"_cell_angle_alpha   {self.alpha:12.6f}", file=file_)
            print(f"_cell_angle_beta    {self.beta:12.6f}",  file=file_)
            print(f"_cell_angle_gamma   {self.gamma:12.6f}", file=file_)
            print(f"_cell_volume        {volume:12.6f}",     file=file_)
            print("loop_",                                           file=file_)
            print("_atom_site_label",                                file=file_)
            print("_atom_site_type_symbol",                          file=file_)
            print("_atom_site_fract_x",                              file=file_)
            print("_atom_site_fract_y",                              file=file_)
            print("_atom_site_fract_z",                              file=file_)
            for i,(frac,el) in enumerate(zip(self.fractional_coordinates,self.element_names)):
                symbol = el+str(i+1)
                print( f"{symbol} {el} {frac[0]:12.6f} {frac[1]:12.6f} {frac[2]:12.6f}", file=file_)
            print(" ",                                               file=file_)
            print("#END",                                            file=file_)
            print(" ",                                               file=file_)
        return

    def getBoundingBox(self, originXYZ = None, originABC = None, units="Angstrom"):
        """Generate the corners and edges of a bounding box.

        This method calculates the corners and edges of a bounding box based on predefined coordinates. These coordinates are transformed using a conversion method before being paired into edges.

        Parameters
        ----------
        originABC : a list of 3 floats
            Defines the centre about which the box is centred
            The basis for the coordinates are abc
        originXYZ : a list of 3 floats
            Defines the centre about which the box is centred
            The basis for the coordinates are xyz
        units : str
            An optional unit of length required for output, default is Angstrom

        Returns
        -------
        tuple of list
            A tuple containing two elements:
            - The first element is a list of corners' coordinates after conversion (`list` of `np.ndarray`).
            - The second element is a list of tuples, each consisting of a pair of corners representing an edge (`list` of `tuple`).

        Example
        -------
        To fetch the bounding box coordinates and edges, simply call: ::

            corners_xyz, edges = object.getBoundingBox()

        """        
        if originABC is None and originXYZ is None:
            originABC = np.zeros( 3 )
            originXYZ = np.zeros( 3 )
        elif originABC is not None:
            originXYZ = self.convert_abc_to_xyz(originABC)
        elif originXYZ is not None:
            originABC = self.convert_xyz_to_abc(originXYZ)
        #
        # Return a box with 12 edges which represent the unit cell in cartesian space
        #  7---6
        # /|  /|
        # 3---2|
        # |4--|5
        # 0---1/
        # Corners in abc space
        corners = []
        corners.append( np.array( [0.0,0.0,0.0] ) )
        corners.append( np.array( [1.0,0.0,0.0] ) )
        corners.append( np.array( [1.0,1.0,0.0] ) )
        corners.append( np.array( [0.0,1.0,0.0] ) )
        corners.append( np.array( [0.0,0.0,1.0] ) )
        corners.append( np.array( [1.0,0.0,1.0] ) )
        corners.append( np.array( [1.0,1.0,1.0] ) )
        corners.append( np.array( [0.0,1.0,1.0] ) )
        labels = [ "o", "a", "", "b", "c", "", "", "" ]
        centre = np.array ( [0.0, 0.0, 0.0] )
        for corner in corners:
            centre += corner
        centre = centre / 8.0
        corners = [ corner - centre + originABC for corner  in corners ]
        # Now calculate the corners in cartesian space
        corners_xyz = [ self.convert_abc_to_xyz(corner) for corner in corners ]
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
        # Convert to the required output unit
        corners = convert_length_units(corners,self.units,units)
        edges   = convert_length_units(edges  ,self.units,units)
        return corners_xyz,edges,labels

    def getDensity(self, units="cm"):
        """Calculate the density of the crystal.
  
        Returns the density in g/cc.  If the mass is not known, then returns 1.0

        Parameters
        ----------
        units : str
            Optional parameter giving the output unit of length, default is cm, so volume will be cm3

        Returns
        -------
        float
            The density in gms/cc

        """
        volume = self.getVolume(units=units)
        mass = sum(self.atomic_masses) / avogadro_si
        if mass == 0:
            self.density = 1.0
        else:
            self.density = mass/volume
        return self.density

    def print(self):
        """Print the details of the given unit cell.

        This method prints formatted details of the unit cell object, including lattice parameters (a, b, c, alpha, beta, gamma), lattice vectors, element names, fractional coordinates, Cartesian coordinates, and molecular information if any molecules are defined within the unit cell.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        spg_symbol, spg_number = self.find_symmetry()
        print      ("Space group international symbol is: ",spg_symbol)
        print      ("Space group number is              : ",spg_number)
        print      ("Units for length are: ", self.units)
        print_reals("Unit Cell a,b,c ",[self.a, self.b, self.c], format="{:12.6f}")
        print_reals("Unit Cell alpha,beta,gamma",[self.alpha, self.beta, self.gamma], format="{:12.6f}")
        print_reals("lattice", self.lattice[0], format="{:12.6f}")
        print_reals("", self.lattice[1], format="{:12.6f}")
        print_reals("", self.lattice[2], format="{:12.6f}")
        print_strings("Element names",self.element_names)
        print_reals("Element masses",self.atomic_masses)
        if len(self.element_names) > 0:
            print("Fractional coords")
            for frac in self.fractional_coordinates:
                print_reals("",frac, format="{:12.6f}")
            print("Cartesian coords")
            for xyz in self.xyz_coordinates:
                print_reals("",xyz, format="{:12.6f}")
            if self.molecules:
                for molid,atoms in enumerate(self.molecules):
                    mass, cm_xyz, cm_frac = self.calculateCentreOfMass(atom_list=atoms,output="all")
                    molstring = "Molecule "+str(molid)+":"
                    print_ints("Atoms in "+molstring,atoms)
                    print_reals("Mass of "+molstring,[ mass ], format="{:12.6f}")
                    print_reals("Centre of Mass  (xyz) of "+molstring, cm_xyz, format="{:12.6f}")
                    print_reals("Centre of Mass (frac) of "+molstring, cm_frac,format="{:12.6f}")

    def calculateCentreOfMass(self,atom_list=None, output="xyz"):
        """Calculate the centre of mass for a given set of atoms.

        Parameters
        ----------
        atom_list : list, optional
            A list of atom indices for which the centre of mass should be calculated.
            If None, the centre of mass will be calculated for all atoms in the system. Default is None.
        output : {'xyz', 'mass', 'abc'}, optional
            The output requested:
            - 'xyz' returns the centre of mass coordinates in Cartesian (x, y, z) units.
            - 'mass' returns the total mass of the atoms in `atom_list`.
            - 'abc' returns the centre of mass in fractional (a, b, c) coordinates.
            If not specified, 'xyz' is used as default.

        Returns
        -------
        numpy.ndarray or float or tuple
            The centre of mass as requested by output:
            - If 'xyz', returns a numpy array with the x, y, z coordinates of the centre of mass.
            - If 'mass', returns a float representing the total mass of the specified atoms.
            - If 'abc', returns a numpy array with the a, b, c fractional coordinates of the centre of mass.
            - If the ouput are not recognized, a tuple containing the total mass, Cartesian coordinates, and fractional coordinates of the centre of mass is returned.

        See Also
        --------
        convert_xyz_to_abc : A method to convert Cartesian coordinates to fractional coordinates.

        Examples
        --------
        >>> calculateCentreOfMass(atom_list=[1, 2, 3], output='xyz')
        array([0.1, 0.2, 0.3])
        >>> calculateCentreOfMass(output='mass')
        50.2
        >>> calculateCentreOfMass(output='abc')
        array([0.4, 0.5, 0.6])

        """        
        # Calculate the centre of mass - if the atom list is given just use that
        # The centre of mass can be returned in as 'xyz' space or 'abc' space
        # if output='all' a tuple of (mass,cm_xyz,cm_abc) is returned
        if atom_list is None:
            atom_list = range(self.nions)
        mass = 0.0
        cm = np.zeros(3)
        for atom_index in atom_list:
            mass += self.atomic_masses[atom_index]
            cm = cm + self.atomic_masses[atom_index] * self.xyz_coordinates[atom_index]
        cm_xyz = cm / mass
        cm_fractional = self.convert_xyz_to_abc(cm_xyz)
        return {
            "xyz"  : cm_xyz,
            "mass" : mass,
            "abc"  : cm_fractional,
            "all"  : (mass, cm_xyz, cm_fractional ),
        }.get(output, (mass, cm_xyz, cm_fractional))

    def get_number_of_atoms(self):
        """Get the number of atoms.

        Parameters
        ----------
        None

        Returns
        -------
        int
            The number of atoms

        """        
        return self.nions

    def get_number_of_molecules(self):
        """Get the number of molecules.

        Parameters
        ----------
        None

        Returns
        -------
        int
            The number of molecules

        """        
        return [element_to_atomic_number[el] for el in self.element_names]

    def get_atomic_numbers(self):
        """Get the atomic numbers for the elements.

        Parameters
        ----------
        None

        Returns
        -------
        list of ints
            The atomic numbers

        """        
        return [element_to_atomic_number[el] for el in self.element_names]

    def get_atomic_masses(self):
        """Get the atomic masses for the elements.

        Parameters
        ----------
        None

        Returns
        -------
        masses : list
            The atomic masses 

        """        
        return self.atomic_masses

    def set_atomic_masses(self, masses):
        """Set the atomic masses for the elements.

        Parameters
        ----------
        masses : list
            The atomic masses to be assigned. This can be a list or a dictionary of atomic masses.

        Returns
        -------
        None

        """        
        self.atomic_masses = masses

    def convert_unitcell_to_abc(self):
        """Convert a unit cell to the equivalent a, b, c, alpha, beta, gamma designation.

        Parameters
        ----------
        None

        Returns
        -------
        a, b, c : float
            The lengths of the cell edges.
        alpha, beta, gamma : float
            The angles between the edges in degrees.

        """
        a = np.sqrt(np.dot(self.lattice[0], self.lattice[0]))
        b = np.sqrt(np.dot(self.lattice[1], self.lattice[1]))
        c = np.sqrt(np.dot(self.lattice[2], self.lattice[2]))
        gamma = np.arccos(np.dot(self.lattice[0], self.lattice[1]) / (a*b))
        beta = np.arccos(np.dot(self.lattice[0], self.lattice[2]) / (a*c))
        alpha = np.arccos(np.dot(self.lattice[1], self.lattice[2]) / (b*c))
        return a, b, c, np.degrees(alpha), np.degrees(beta), np.degrees(gamma)

    def convert_abc_to_unitcell(self, a, b, c, alpha, beta, gamma):
        """Convert a, b, c, alpha, beta, gamma to a unit cell.

        Parameters
        ----------
        a : type
            Unit cell `a`.
        b : type
            Unit cell `b`.
        c : type
            Unit cell `c`.
        alpha : type
            Unit cell alpha in degrees.
        beta : type
            Unit cell beta` in degrees.
        gamma : type
            Unit cell gamma` in degrees.

        Returns
        -------
        lattice : list of floats (3,3) list
            The lattice parameters

        """
        # This is castep convention, need to check it works with vasp
        alpha = np.radians(alpha)
        beta = np.radians(beta)
        gamma = np.radians(gamma)
        lattice = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        lattice[2] = [0.0, 0.0, c]
        lattice[1] = [0.0, b*np.sin(alpha), b*np.cos(alpha)]
        z = a * np.cos(beta)
        y = a * (np.cos(gamma) - np.cos(alpha)*np.cos(beta))/np.sin(alpha)
        x = a * np.sqrt(1.0 - np.cos(beta)**2 - (y/a)**2)
        lattice[0] = [x, y, z]
        self._calculate_reciprocal_lattice()
        return self.lattice

    def _calculate_reciprocal_lattice(self):
        """Calculate the reciprocal lattice.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Modifies
        --------
        Sets the reciprocal lattice vectors from the lattice vectors

        """
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = self.convert_unitcell_to_abc()
        self.reciprocal_lattice = np.linalg.inv(self.lattice)

    def getVolume(self, units="cm"):
        """Calculate the volume.

        Parameters
        ----------
        units : str
            An optional string indicating the required unit of length for the output value, default is cm

        Returns
        -------
        volume : float

        """
        lattice = convert_length_units(self.lattice,self.units,units)
        return np.abs(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))

    def convert_xyz_to_abc(self, xyz):
        """Convert xyz coordinates to abc lattice coordinates.

        Parameters
        ----------
        xyz: list of coordinates
            List of coordinates in xyz space

        Returns
        -------
        abc : list of 3 floats
            List of fractional coordinates

        """
        xyz = np.array(xyz)
        return np.dot(xyz, self.reciprocal_lattice)

    def convert_xyz_to_integer_abc(self, xyz):
        """Convert xyz coordinates to integer abc indices.

        Parameters
        ----------
        xyz : list 3 floats
            List of xyz

        Returns
        -------
        abc : list of abc coordinates

        """
        xyz = np.array(xyz)
        abc_array =  np.dot(xyz, self.reciprocal_lattice)
        #
        # Find the common denominator
        #
        return self.find_lcd(abc_array)

    def find_lcd(self, numbers):
        """Find the approximated lowest common denominator.

        Return the list as integers divided by the lowest common denominator

        Parameters
        ----------
        numbers : array of floats
            List of numbers

        Returns
        -------
        results : list of integers

        """
        smallest = 1.0E99
        for x in numbers:
            if abs(x) > 1.0E-2:
                smallest = min(abs(x), smallest)
        numbers = numbers/smallest
        result = 1
        # largest = int(np.max(np.round(np.abs(numbers))))
        largest = 10
        for lcd in range(1,largest):
            new = numbers*lcd
            rounded = np.round(new)
            diff = np.sum( np.abs(new-rounded) )
            if diff < 1.0E-8:
                result = lcd
                break
        results = numbers * result
        maximum = np.max(np.abs(results))
        if maximum > 20:
            results = 20/maximum * results 
        rounded = np.round(results)
        return [ round(r) for r in rounded ]

    def convert_hkl_to_xyz(self, hkl):
        """Convert hkl miller indices to xyz coordinates.

        Parameters
        ----------
        hkl : list 3 ints
            List of hkls

        Returns
        -------
        xyz : list of xyz coordinates

        """
        hkl = np.array(hkl)
        return np.dot(hkl, self.reciprocal_lattice.T)

    def convert_abc_to_xyz(self, abc):
        """Convert abc coordinates to xyz coordinates.

        Parameters
        ----------
        abc : list 3 floats
            List of fractional coordinates

        Returns
        -------
        xyz : list of coordinates

        """
        abc = np.array(abc)
        return np.dot(abc, self.lattice)

    def convert_hkl_to_xyz2(self, hkl):
        """Convert hkl coordinates to xyz coordinates (written as a check on convert_hkl_to_xyz).

        Parameters
        ----------
        hkl : list 3 ints
            List of hkls

        Returns
        -------
        xyz : list of xyz coordinates

        """
        # Use the hkl miller indices to calculate the normal to a plane in cartesian space
        p = []
        if hkl[0] != 0:
            p.append(np.array(self.lattice[0] / hkl[0]))
        if hkl[1] != 0:
            p.append(np.array(self.lattice[1] / hkl[1]))
        if hkl[2] != 0:
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
            sys.exit(1)
        return normal / norm

    def set_fractional_coordinates(self, coords):
        """Set the fractional coordinates and calculate the xyz coordinates.

        Parameters
        ----------
        coords : a list 3 floats
            Fractional coordinates

        Returns
        -------
        None

        Modifies
        --------
        fractional_coordinates
        xyz_coordinates
        nions

        """
        self.fractional_coordinates = np.array(coords)
        self.xyz_coordinates = self.convert_abc_to_xyz(coords)
        self.nions = len(coords)
        return

    def set_xyz_coordinates(self, coords, units="Angstrom"):
        """Set the xyz coordinates and calculate the fractional coordinates.

        Parameters
        ----------
        coords : list of 3 floats
            A list of xyz coordinates, the unit of length must agree with the lattice
        units : str
            A unit of length for the input values.  The default is Angstrom.

        Returns
        -------
        None

        Modifies
        --------
        xyz_coordinates
        fractional_coordinates
        nions

        """
        self.xyz_coordinates = convert_length_units(coords,units,self.units)
        self.fractional_coordinates = self.convert_xyz_to_abc(coords)
        self.nions = len(coords)
        return

    def set_element_names(self, element_names):
        """Set the names of elements in the object after cleaning them up.

        Parameters
        ----------
        element_names : list
            A list of strings representing element names to be cleaned and stored.

        Returns
        -------
        None

        Notes
        -----
        This function takes a list of element names, cleans each name using the 
        `cleanup_symbol` function (not defined here), and then updates the object's 
        `element_names` attribute with the cleaned names.

        """        
        self.element_names = []
        for el in element_names:
          el = cleanup_symbol(el)
          self.element_names.append(el)
        return

    def set_atomic_numbers(self, atomic_numbers):
        """Set the element names according to the atomic numbers for a molecule.

        Parameters
        ----------
        atomic_numbers
            A list containing the atomic numbers of the atoms

        Returns
        -------
        None

        Modifies
        --------
        element_names

        Notes
        -----
        This method updates the `atomic numbers` attribute of the molecule

        """        
        self.element_names = [atomic_number_to_element[el] for el in atomic_numbers]
        return

    def set_atom_labels(self, atom_labels):
        """Set the atom labels for a molecule.

        Parameters
        ----------
        atom_labels : list of strings
            A list containing the labels of atoms.

        Returns
        -------
        None

        Notes
        -----
        This method updates the `atom_labels` attribute of the molecule with the
        provided list of atom labels. Each element in the input `atom_labels` list 
        is appended to the `self.atom_labels` attribute.

        """        
        self.atom_labels = []
        for el in atom_labels:
          self.atom_labels.append(el)
        return

    def get_element_names(self):
        """Get the names of elements in the cell.

        Parameters
        ----------
        None

        Returns
        -------
        element_names : list
            A list of strings representing element names to be cleaned and stored.

        """        
        return self.element_names

    def get_atom_labels(self):
        """Get the atom labels for a molecule.

        If the labels have not been set, then a list of labels is created from the element names

        Parameters
        ----------
        None
            A list containing the labels of atoms.

        Returns
        -------
        A list of atom labels

        """        
        if len(self.atom_labels) == 0:
            for i,el in enumerate(self.element_names):
                label = el+str(i+1)
                self.atom_labels.append(label)
        return self.atom_labels

    def find_symmetry(self,symprec=1e-5,angle_tolerance=-1.0):
        """Find the space group symmetry of the unit cell.

        Parameters
        ----------
        symprec : float
            Determine the accuracy of the coordinates for symmetry determination
        angle_tolerance : float
            Determine the accuracy of the angles for symmetry determination

        Returns
        -------
        a tuple
            (international table symbol, number])

        """
        atomic_nos = self.get_atomic_numbers()
        if len(atomic_nos) <= 0:
            return "P 1", 1
        cell = ( self.lattice, self.fractional_coordinates, atomic_nos )
        spacegroup = spglib.get_spacegroup(cell, symprec=symprec,angle_tolerance=angle_tolerance)
        sp = spacegroup.split()
        symbol = sp[0]
        number = int(sp[1].replace("(","").replace(")",""))
        return symbol,number

    def calculate_molecular_contents(self, scale=1.1, tolerance=0.1, radii=None):
        """Find whole molecules in the unit cell.

        Does this by creating a supercell and exploring adjacent cells to see if there is any bonding to atoms in the adjacent cell
        A new unit cell is created containing whole molecules, the order of the atoms in the new cell is different.
        The routine returns the the number of moleculess.

        Parameters
        ----------
        scale : float, optional
            The scale factor applied to the covalent radii. Default is 1.1
        tolerance : float, optional. Default is 0.1
            The tolerance added to the scaled sum of radii to determine the maximum allowable distance between atoms i and j for them to be considered bonded.
        radii : a dictionary, optional
            A dictionary of covalent radii for the atoms, key is the element name.  If not given then the package radii are used from PDielec.Constants

        Returns
        -------
        int : the number of molecules

        Modifies
        --------
        fractional_coordinates
        bonds

        Notes
        -----
        The formula used to calculate the largest distance apart atoms i and j can be for a bond is:
        `scale * (radi + radj) + tolerance`
        The routine modifies its current instance
        The ordering of the atoms is maintained

        """
        # Check that the covalent radius has been given
        if radii is None:
            radii = covalent_radii.copy()
        # Calculate the contents of all the cells adjacent to the central cell
        adjacents = ( 0, -1, 1 )
        translations = [ (i, j, k) for i in adjacents for j in adjacents for k in adjacents ]
        fractional_supercell = []
        index_supercell = []
        for tr in translations:
            for l,a in enumerate(self.fractional_coordinates):
                i,j,k = tr
                index_supercell.append( l )
                new_position = [ (xyz1 + xyz2)  for xyz1, xyz2 in zip(a, tr) ]
                fractional_supercell.append( new_position )
        # Convert fractional supercell coordinates to xyz
        # xyz_supercell will be an np array
        xyz_supercell = np.empty_like(fractional_supercell)
        for i,abc in enumerate(fractional_supercell):
            xyz_supercell[i] = self.convert_abc_to_xyz(abc)
        # put the atoms into boxes of boxSize
        BoxAtoms = {}
        BoxNeighbours = {}
        boxSize = 0.0
        # calculate boxsize
        rmax = 0.0
        for el in self.element_names:
            rmax = max(rmax, radii[el])
        boxSize = 2.0*scale*rmax + 0.5 + tolerance
        # Put atoms into boxes and store the box info in Atom_box_id
        Atom_box_id = []
        for i,xyz in enumerate(xyz_supercell):
            a = int( math.floor(xyz[0]/boxSize) )
            b = int( math.floor(xyz[1]/boxSize) )
            c = int( math.floor(xyz[2]/boxSize) )
            abc = (a,b,c)
            Atom_box_id.append(abc)
            try:
                BoxAtoms[abc].append(i)
            except Exception:
                BoxAtoms[abc] = [i]
        # Calculate the neighbouring boxes for each occupied box
        for abc in BoxAtoms:
            a,b,c = abc
            BoxNeighbours[abc] = []
            for i in [ -1, 0, 1]:
              for j in [ -1, 0, 1]:
                for k in [ -1, 0, 1]:
                  BoxNeighbours[abc].append( (a+i,b+j,c+k) )
        # end for abc in Box1

        # Calculate the bonding the supercell
        bondedToAtom = {}
        for i,xyzi in enumerate(xyz_supercell):
            bondedToAtom[i] = []
            # Find the element name for this atom in the supercell
            ip = index_supercell[i]
            i_el = self.element_names[ip]
            # lookup all of the boxes that might hold a bonding atom
            cell = Atom_box_id[i]
            for abc in BoxNeighbours[cell]:
                # lookup all the atoms that are in that cell
                try:
                    for j in BoxAtoms[abc]:
                        if j < i:
                            # Find the element name for this atom in the supercell
                            jp = index_supercell[j]
                            j_el = self.element_names[jp]
                            dist1 = scale*( radii[i_el] + radii[j_el]) + tolerance
                            dist2 = calculate_distance(xyzi,xyz_supercell[j])
                            if dist2 < dist1:
                                bondedToAtom[i].append(j)
                                bondedToAtom[j].append(i)
                            # end if dist2 < dist1
                        # end if j < i
                    # end for j
                except KeyError:
                    pass
            # end for abc
        # end for i,a
        #
        # Now we have to find how many molecules we have in the cell
        # There must be at least one molecule in the cell and it must contain the first atom
        #
        # BelongsToMolecule[i] holds which molecule the ith atom belongs to
        belongsToMolecule = {}
        # molecules is a dictionary of molecules each entry is a list of atoms
        molecules = {}
        molID = -1
        # We stop when all the atoms in the original cell belong to a molecule
        remainingAtoms = [ atom for atom in range(self.nions) ]
        bonds = []
        while len(remainingAtoms) > 0:
            # create a new molecule from the first atom which has no molecule assigned to it
            molID += 1
            useAtom = remainingAtoms[0]
            belongsToMolecule[useAtom] = molID
            molecules[molID] = [ useAtom ]
            remainingAtoms.remove(useAtom)
            # Now using this new molecule with its only atom as a seed find any atoms connected to it
            # We need to continue searching for bonded atoms until there are no more to be found
            moreAtomsToBeFound = True
            while moreAtomsToBeFound:
                moreAtomsToBeFound = False
                for i in range(len(xyz_supercell)):
                    # has this atom been assigned a molecule yet?
                    if i in belongsToMolecule:
                        # atom i is already assigned to a molecule
                        useThisMolecule = belongsToMolecule[i]
                        # Go through all the atoms bonded to i and add to the current molecule
                        for j in bondedToAtom[i]:
                            jx = index_supercell[j]
                            # The image of j in the original cell might not be available, and j might be bonded
                            if jx in remainingAtoms and j not in belongsToMolecule:
                                # if j was not already in a molecule then we have new information
                                moreAtomsToBeFound = True
                                molecules[useThisMolecule].append(j)
                                belongsToMolecule[j] = useThisMolecule
                                remainingAtoms.remove(jx)
                            # The j'th atom could be already specified and we have a ring....
                            # We also need to make sure that we have unique bonds
                            if j in belongsToMolecule:
                                if i > j and (i,j) not in bonds:
                                    bonds.append( (i,j) )
                                elif i < j and (j,i) not in bonds:
                                    bonds.append( (j,i) )
                        # end for j
                    # end if i in
                # end for i
            # while moreAtomsToBeFound
        # until all the atoms belong to a molecule
        self.centres_of_mass = []
        self.total_mass = 0.0
        for mol_index in molecules:
            for atom_index in molecules[mol_index]:
                index = index_supercell[atom_index]
            # Calculate centre of mass
            mass = 0.0
            cm = np.zeros(3)
            for atom_index in molecules[mol_index]:
                index = index_supercell[atom_index]
                mass += self.atomic_masses[index]
                cm = cm + self.atomic_masses[index] * xyz_supercell[atom_index]
            cm_xyz = cm / mass
            cm_fractional = self.convert_xyz_to_abc(cm_xyz)
            self.centres_of_mass.append( cm_fractional )
        #
        # Create a new unit cell with the atoms shifted so that
        # whole molecules are ordered and within the cell
        #
        new_fractional_coordinates = np.empty_like(self.fractional_coordinates)
        for mol_index in molecules:
            cm = self.centres_of_mass[mol_index]
            shift = -np.floor(cm)
            for atom_index in molecules[mol_index]:
                old_index = index_supercell[atom_index]
                new_fractional_coordinates[old_index] = shift + np.array(fractional_supercell[atom_index])
        self.set_fractional_coordinates(new_fractional_coordinates)
        #
        # Replace the super cell id, with the index of the atom in the original cell
        #
        molecule_list = []
        for mol_index in molecules:
            atom_list = []
            for atom_index in molecules[mol_index]:
                atom_list.append(index_supercell[atom_index])
            molecule_list.append(atom_list)
        self.set_molecules(molecule_list)
        bond_list = []
        for i,j in bonds:
            bond_list.append( (index_supercell[i], index_supercell[j]) )
        self.set_bonds(bond_list)
        return len(molecule_list)

    def get_bonds(self):
        """Return a list of bonds for the unit cell.

        It also returns a list of the bond lengths in angstrom

        Parameters
        ----------
        None

        Returns
        -------
        (list of bonds, list of bondlengths)
            list of bonds is a list of pairs of integers denoting a bond
            list of bond lengths is a list of floats

        """
        bond_lengths = []
        for i,j in self.bonds:
            bond_length = calculate_distance(self.xyz_coordinates[i],self.xyz_coordinates[j])
            bond_length = convert_length_units(bond_length,self.units,"Angstrom")
            bond_lengths.append(bond_length)
        return self.bonds, bond_lengths

    def get_bond_angles(self):
        """Return a list of atoms that form bonded angles for the unit cell.

        It also returns a list of the angles in degrees

        Parameters
        ----------
        None

        Returns
        -------
        (list of 3 integer tuples, list of bond angles)
            list of angles is a list of 3 integers denoting involved in the angle
            list of angles is a list of floats

        """
        angle_list = []
        for i,j in self.bonds:
            for k,l in self.bonds:
                if i == k and j != l:
                    if l > j and (l, i, j) not in angle_list:
                        angle_list.append( (l, i, j) )
                    elif (j, i, l) not in angle_list:
                        angle_list.append( (j, i, l) )
                if i == l and j != k:
                    if k > j and (k, i, j) not in angle_list:
                        angle_list.append( (k, i, j) )
                    elif (j, i, k) not in angle_list:
                        angle_list.append( (j, i, k) )
        angles = []
        for atoms in angle_list:
            i,j,k = atoms
            a = self.xyz_coordinates[i]
            b = self.xyz_coordinates[j]
            c = self.xyz_coordinates[k]
            angle = calculate_angle(a,b,c)
            angles.append(angle)
        return angle_list, angles

    def get_torsions(self):
        """Return a list of atoms that form torsion angles for the unit cell.

        It also returns a list of the angles in degrees

        Parameters
        ----------
        None

        Returns
        -------
        (list of 4 integer tuples, list of torsion angles)
            list of torsions is a list of 4 integers denoting involved in the angle
            list of angles is a list of floats

        """
        angle_list, angles = self.get_bond_angles()
        # Generate a list of torsion angles
        torsion_list = []
        for i,j,k in angle_list:
            for l1,l2 in self.bonds:
                if l1 == i and l2 != j and (l2,i,j,k) not in torsion_list:
                    torsion_list.append( (l2,i,j,k) )
                if l2 == i and l1 != j and (l1,i,j,k) not in torsion_list:
                    torsion_list.append( (l1,i,j,k) )
                if l1 == k and l2 != j and (i,j,k,l2) not in torsion_list:
                    torsion_list.append( (i,j,k,l2) )
                if l2 == k and l1 != j and (i,j,k,l1) not in torsion_list:
                    torsion_list.append( (i,j,k,l1) )
        # Calculate the value of the torsion angle
        angles = []
        for i,j,k,l in torsion_list:
            a = self.xyz_coordinates[i]
            b = self.xyz_coordinates[j]
            c = self.xyz_coordinates[k]
            d = self.xyz_coordinates[l]
            angle = calculate_torsion(a,b,c,d)
            angles.append(angle)
        return torsion_list, angles

    def set_bonds(self, bonds):
        """Define a list of bonds for the unit cell.

        Some checking is performed.  If the bonds has duplicates but in a different order, then they are removed.

        Parameters
        ----------
        bonds : a list of integer tuple
           The bonding list for the cell

        Returns
        -------
        None

        """
        # Check that the bonds are unique
        new_list = []
        for i,j in bonds:
            if (j,i) in bonds:
                if i > j:
                    new_list.append( (i,j) )
            elif i>j:
                new_list.append( (i,j) )
            else:
                new_list.append( (j,i) )
        if len(new_list) > 0:
            self.bonds = new_list
        else:
            self.bonds = bonds
        return

    def set_molecules(self, molecules):
        """Define a list of molecules, each molecule is a list of atom coordinates.

        Parameters
        ----------
        molecules : list of 3 floats
            A list of molecules

        Returns
        -------
        None

        """
        self.molecules = molecules
        return

    def get_species(self):
        """Return a list of the unique species in the cell.

        The unique list needs to maintain the order seen in the original list

        Parameters
        ----------
        None

        Returns
        -------
        A list of strings

        """
        result = []
        for name in self.element_names:
            if name not in result:
                result.append(name)
        return result

    def fold_cell(self, symprec=1.0E-5):
        """Fold a cell so that any atoms outside a single cell are moved inside.

        Any molecule or bonding information is removed.  
        Duplication of atoms is avoided

        Parameters
        ----------
        symprec : float
            The symmetry precision

        Returns
        -------
        map_old_to_new : a list of integers
        map_new_to_old : a list of integers

        Modifies
        --------
        fractional_coordinates
        atom_labels
        atom_masses
        element_names

        """
        new_labels = []
        new_masses = []
        new_names = []
        new_positions = []
        positions = self.fractional_coordinates
        labels = self.get_atom_labels()
        masses  = self.get_atomic_masses()
        element_names  = self.get_element_names()
        map_new_to_old = []
        map_old_to_new = []
        new_index = 0
        for old_index,(p,l,m,e) in enumerate(zip(positions,labels,masses,element_names)):
            for i in range(3):
                if p[i] >= 1.0 or p[i] < 0.0:
                    p[i] = p[i] - math.floor(p[i])
            #
            # Only store the new atom position if it doesn't alread exist
            # This could be done better, but will do for the time being
            # For very large cells, an order N method will be needed
            #
            dist = 1.0
            for frac in new_positions:
                dist = min(dist, np.sqrt(np.sum(np.square(p-frac))))
            map_old_to_new.append(new_index)
            if dist > symprec:
                new_positions.append(p)
                new_labels.append(l)
                new_masses.append(m)
                new_names.append(e)
                map_new_to_old.append(old_index)
                new_index += 1
        self.set_atom_labels(new_labels)
        self.set_atomic_masses(new_masses)
        self.set_element_names(new_names)
        self.set_fractional_coordinates(new_positions)
        #
        # Remove molecule and bonding information
        #
        self.molecules = []
        self.bonds     = []
        return map_old_to_new, map_new_to_old

    def trim_cell(self):
        """Trim this cell so that any atoms outside a single cell are removed.

        As well as removing the redundant atoms, a list of the mappings
        from the original cell and to the old cell are generated

        This routine is used when generating a primitive (the new) cell, from a standard cell

        Parameters
        ----------
        None

        Returns
        -------
        map_old_to_new : a list of integers
        map_new_to_old : a list of integers

        Modifies
        --------
        fractional_coordinates
        atom_labels
        atom_masses
        element_names

        """
        positions = self.fractional_coordinates
        labels = self.get_atom_labels()
        masses  = self.get_atomic_masses()
        element_names  = self.get_element_names()
        keep_list = []
        new_positions = []
        new_labels = []
        new_names = []
        new_masses = []
        map_new_to_old = []
        map_old_to_new = []
        index_primitive = 0
        for index_cell,(p,l,m,n) in enumerate(zip(positions, labels, masses, element_names)):
            keep = True
            for i in range(3):
                if p[i] >= 1.0 or p[i] < 0.0:
                    keep = False
            map_old_to_new.append(index_primitive)
            if keep:
                keep_list.append(keep)
                new_positions.append(p)
                new_labels.append(l)
                new_masses.append(m)
                new_names.append(n)
                map_new_to_old.append(index_cell)
                index_primitive += 1
        # end for p,l,m,n
        #
        # Store the information in the current cell
        #
        self.set_atom_labels(new_labels)
        self.set_atomic_masses(new_masses)
        self.set_element_names(new_names)
        self.set_fractional_coordinates(new_positions)
        self.molecules = []
        self.bonds     = []
        return map_old_to_new, map_new_to_old

    def guess_primitive_transform(self, symprec=1.0E-5):
        """Guess the primitive cell on the basis of the centering.

        This routine is based on a similar routine in phonopy get_primitive_matrix_by_centring

        Parameters
        ----------
        symprec : float (optional)
            The symmetry precision used in determining the spacegroup

        Returns
        -------
        A guess at the new unit cell

        """
        spglib_cell = (self.lattice, 
                       self.fractional_coordinates, 
                       self.get_atomic_numbers())
        dataset = spglib.get_symmetry_dataset(spglib_cell, symprec=symprec)
        tmat = dataset.transformation_matrix
        centring = dataset.international[0]
        if centring == "P":
            pmat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype="double")
        elif centring == "F":
            pmat = np.array(
              [[0.0, 1.0 / 2, 1.0 / 2], [1.0 / 2, 0, 1.0 / 2], [1.0 / 2, 1.0 / 2, 0.0]]
            )
        elif centring == "I":
            pmat = np.array(
              [
                [-1.0 / 2, 1.0 / 2, 1.0 / 2],
                [1.0 / 2, -1.0 / 2, 1.0 / 2],
                [1.0 / 2, 1.0 / 2, -1.0 / 2],
              ]
            )
        elif centring == "A":
            pmat = np.array(
                [[1.0, 0.0, 0.0], [0.0, 1.0 / 2, -1.0 / 2], [0.0, 1.0 / 2, 1.0 / 2]]
            )
        elif centring == "C":
            pmat = np.array(
              [[1.0 / 2, 1.0 / 2, 0], [-1.0 / 2, 1.0 / 2, 0], [0.0, 0.0, 1.0]]
            )
        elif centring == "R":
            pmat = np.array(
              [
                [2.0 / 3, -1.0 / 3, -1.0 / 3],
                [1.0 / 3, 1.0 / 3, -2.0 / 3],
                [1.0 / 3, 1.0 / 3, 1.0 / 3],
              ]
            )
        else:
            pmat = np.eye(3)
            print("Centring is not recognised",centring)
        return np.dot(np.linalg.inv(tmat), pmat)

