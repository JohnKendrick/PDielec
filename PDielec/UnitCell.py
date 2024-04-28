#!/usr/bin/python
'''
Hold unit cell information and its associated calculated properties.

Copyright 2024 John Kendrick

This file is part of PDielec

This program is free software; you can redistribute it and/or modify
it under the terms of the MIT License

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

You should have received a copy of the MIT License
along with this program, if not see https://opensource.org/licenses/MIT
'''

import numpy as np
import math
from PDielec.Calculator  import calculate_distance
from PDielec.Plotter     import print_reals, print_ints, print_strings
from PDielec.Calculator  import cleanup_symbol
from PDielec.Constants   import avogadro_si, element_to_atomic_number
import sys

def convert_length_units(value, units_in, units_out):
    """"
    Convert between different length units

    The 'internal' unit is taken to be the Angstrom so units are relative to the Angstrom
    The unit strings are made lowercase, so case should be irrelevant

    Parameters
    ----------
    units_in : str
        The units of the input value(s). Can be one of '
    units_out : str
        The units of the output value(s). Must be one of 'a.u. au bohr ang angs angstrom nm um mm cm m'
    input_value : float or list of floats or a numpy array
        The value(s) for which the conversion is to be made.

    Returns
    -------
    scalar or numpy array
        The converted value(s) in the output units specified.

    Notes
    -----
    The input can be either a scalar value, a list or a numpy array of values. The function will return the converted value(s) in the output units specified.
    """
    # the conversion dictionary has a value that converts the key unit to angstroms
    angstroms = { 'a.u.'       : 0.5291772, 
                  'au'         : 0.5291772, 
                  'bohr'       : 0.5291772, 
                  'ang'        : 1.0000000,
                  'angs'       : 1.0000000,
                  'angstrom'   : 1.0000000,
                  'angstroms'  : 1.0000000,
                  'nm'         : 1.0E1    ,
                  'um'         : 1.0E4    ,
                  'mm'         : 1.0E7    ,
                  'cm'         : 1.0E8    ,
                  'm'          : 1.0E10   ,
                 }
    units_in  = units_in.lower()
    units_out = units_out.lower()
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
    """
    Hold unit cell information and its associated calculated properties.
    """
    def __init__(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None, units='Angstrom'):
        """
        Initialize the class instance with optional lattice parameters and calculate the reciprocal lattice.

        Parameters
        ----------
        a, b, c : float or array, optional
            Lattice vectors or cell lengths. If not specified, they default to [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], and [0.0, 0.0, 1.0], respectively.
        alpha, beta, gamma : float, optional
            Lattice angles (in degrees). These are only used if all three angles are specified, otherwise, the default lattice (orthorhombic) is used.
        units : str
            An optional unit such as 'a.u., au bohr angs angstrom Angs Angstrom or nm'  The default is Angstrom.

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
        self.units = units
        if gamma is not None:
            self.lattice = self.convert_abc_to_unitcell(a, b, c, alpha, beta, gamma)
        else:
            self.lattice[0] = a
            self.lattice[1] = b
            self.lattice[2] = c
        self._calculate_reciprocal_lattice(self.lattice)

    def write_cif(self, filename=None, file_=sys.stdout):
        """
        Write the crystallographic information file (CIF) representation of a structure.

        Parameters
        ----------
        filename : str, optional
            The name of the file to be written. If not provided, the CIF data is printed to stdout.
        file\_ : object, optional
            An output stream object. Default is sys.stdout.

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

        Examples
        --------
        To print the CIF representation to standard output, simply call the method without arguments:

        >>> self.write_cif()

        To write the CIF representation to a file named "example.cif":

        >>> cell.write_cif('example.cif')
        """        
        # Open the filename if it is given
        if filename is not None:
            file_ = open(filename,'w')
        abc = convert_length_units( [self.a, self.b, self.c],self.units,'Angstrom' )
        volume = self.getVolume('Angstrom')
        spg_symbol, spg_number = self.find_symmetry()
        if filename is not None:
            print('data_'+filename,file=file_)
        else:
            print('data_',         file=file_)
        print('_symmetry_space_group_name_H-M \'{}\''.format(spg_symbol),file=file_)
        print('_symmetry_Int_Tables_number      {}  '.format(spg_number),file=file_)
        print('_cell_length_a      {:12.6f}'.format(abc[0]),     file=file_)
        print('_cell_length_b      {:12.6f}'.format(abc[1]),     file=file_)
        print('_cell_length_c      {:12.6f}'.format(abc[2]),     file=file_)
        print('_cell_angle_alpha   {:12.6f}'.format(self.alpha), file=file_)
        print('_cell_angle_beta    {:12.6f}'.format(self.beta),  file=file_)
        print('_cell_angle_gamma   {:12.6f}'.format(self.gamma), file=file_)
        print('_cell_volume        {:12.6f}'.format(volume),     file=file_)
        print('loop_',                                           file=file_)
        print('_atom_site_label',                                file=file_)
        print('_atom_site_type_symbol',                          file=file_)
        print('_atom_site_fract_x',                              file=file_)
        print('_atom_site_fract_y',                              file=file_)
        print('_atom_site_fract_z',                              file=file_)
        for i,(frac,el) in enumerate(zip(self.fractional_coordinates,self.element_names)):
            symbol = el+str(i+1)
            print( '{} {} {:12.6f} {:12.6f} {:12.6f}'.format(symbol, el, frac[0], frac[1], frac[2]), file=file_)
        print(' ',                                               file=file_)
        print('#END',                                            file=file_)
        print(' ',                                               file=file_)
        return

    def getBoundingBox(self, units='Angstrom'):
        """
        Generate the corners and edges of a bounding box.

        This method calculates the corners and edges of a bounding box based on predefined coordinates. These coordinates are transformed using a conversion method before being paired into edges.

        Parameters
        ----------
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
        return corners_xyz,edges

    def getDensity(self, units='cm'):
        '''
        Calculate the density of the crystal.
  
        Returns the density in g/cc.  If the mass is not known, then returns 1.0

        Parameters
        ----------
        units : str
            Optional parameter giving the output unit of length, default is cm, so volume will be cm3

        Returns
        -------
        float
            The density in gms/cc
        '''
        volume = self.getVolume(units=units)
        mass = sum(self.atomic_masses) / avogadro_si
        if mass == 0:
            self.density = 1.0
        else:
            self.density = mass/volume
        return self.density

    def print(self):
        """
        Prints the details of the given unit cell.

        This method prints formatted details of the unit cell object, including lattice parameters (a, b, c, alpha, beta, gamma), lattice vectors, element names, fractional coordinates, Cartesian coordinates, and molecular information if any molecules are defined within the unit cell.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """        

        spg_symbol, spg_number = self.find_symmetry()
        print      ('Space group international symbol is: ',spg_symbol)
        print      ('Space group number is              : ',spg_number)
        print      ('Units for length are: ', self.units)
        print_reals('Unit Cell a,b,c ',[self.a, self.b, self.c], format='{:12.6f}')
        print_reals('Unit Cell alpha,beta,gamma',[self.alpha, self.beta, self.gamma], format='{:12.6f}')
        print_reals('lattice', self.lattice[0], format='{:12.6f}')
        print_reals('', self.lattice[1], format='{:12.6f}')
        print_reals('', self.lattice[2], format='{:12.6f}')
        print_strings('Element names',self.element_names)
        print_reals('Element masses',self.atomic_masses)
        if len(self.element_names) > 0:
            print_reals('Fractional coords',self.fractional_coordinates[0], format='{:12.6f}')
            for frac in self.fractional_coordinates[1:]:
                print_reals('',frac, format='{:12.6f}')
            print_reals('Cartesian coords',self.xyz_coordinates[0], format='{:12.6f}')
            for xyz in self.xyz_coordinates[1:]:
                print_reals('',xyz, format='{:12.6f}')
            if self.molecules:
                for molid,atoms in enumerate(self.molecules):
                    mass, cm_xyz, cm_frac = self.calculateCentreOfMass(atom_list=atoms,output='all')
                    molstring = 'Molecule '+str(molid)+':'
                    print_ints('Atoms in '+molstring,atoms)
                    print_reals('Mass of '+molstring,[ mass ], format='{:12.6f}')
                    print_reals('Centre of Mass  (xyz) of '+molstring, cm_xyz, format='{:12.6f}')
                    print_reals('Centre of Mass (frac) of '+molstring, cm_frac,format='{:12.6f}')

    def calculateCentreOfMass(self,atom_list=None, output='xyz'):
        """
        Calculate the centre of mass for a given set of atoms.

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
        if atom_list == None:
            atom_list = range(self.nions)
        mass = 0.0
        cm = np.zeros(3)
        for atom_index in atom_list:
            mass += self.atomic_masses[atom_index]
            cm = cm + self.atomic_masses[atom_index] * self.xyz_coordinates[atom_index]
        cm_xyz = cm / mass
        if output == 'xyz':
            return cm_xyz
        elif output == 'mass':
            return mass
        elif output == 'abc':
            cm_fractional = self.convert_xyz_to_abc(cm_xyz)
            return cm_fractional
        else:
            cm_fractional = self.convert_xyz_to_abc(cm_xyz)
            return mass, cm_xyz, cm_fractional

    def get_atomic_numbers(self):
        """
        Get the atomic numbers for the elements.

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
        """
        Get the atomic masses for the elements.

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
        """
        Set the atomic masses for the elements.

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
        """
        Convert a unit cell to the equivalent a, b, c, alpha, beta, gamma designation.

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

    def convert_abc_to_unitcell(self, a, b, c, alpha_degs, beta_degs, gamma_degs):
        """
        Convert a, b, c, alpha, beta, gamma to a unit cell

        Parameters
        ----------
        a : type
            Description of parameter `a`.
        b : type
            Description of parameter `b`.
        c : type
            Description of parameter `c`.
        alpha : type
            Description of parameter `alpha`.
        beta : type
            Description of parameter `beta`.
        gamma : type
            Description of parameter `gamma`.

        Returns
        -------
        lattice : list of floats (3,3) list
            The lattice parameters
        """
        # This is castep convention, need to check it works with vasp
        alpha = np.radians(alpha_degs)
        beta = np.radians(beta_degs)
        gamma = np.radians(gamma_degs)
        lattice = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        lattice[2] = [0.0, 0.0, c]
        lattice[1] = [0.0, b*np.sin(alpha), b*np.cos(alpha)]
        z = a * np.cos(beta)
        y = a * (np.cos(gamma) - np.cos(alpha)*np.cos(beta))/np.sin(alpha)
        x = a * np.sqrt(1.0 - np.cos(beta)**2 - (y/a)**2)
        lattice[0] = [x, y, z]
        self._calculate_reciprocal_lattice(lattice)
        return self.lattice

    def _calculate_reciprocal_lattice(self, lattice):
        """
        Calculate the reciprocal lattice

        Parameters
        ----------
        3x3 list of floats
            Lattice as [a, b, c] where a/b/c are [x,y,z]

        Returns
        -------
        None

        Returns
        -------
        sets the reciprocal lattice vectors
        """
        self.lattice = np.array(lattice)
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = self.convert_unitcell_to_abc()
        self.inverse_lattice = np.linalg.inv(self.lattice)
        self.reciprocal_lattice = self.inverse_lattice

    def getVolume(self, units='cm'):
        """
        Calculate the volume

        Parameters
        ----------
        units : str
            An optional string indicating the required unit of length for the output value, default is cm

        Returns
        -------
        volume : float
        """
        lattice = convert_length_units(self.lattice,self.units,units)
        volume = np.abs(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))
        return volume

    def convert_xyz_to_abc(self, xyz):
        """
        Convert xyz coordinates to abc lattice coordinates.

        Parameters
        ----------
        xyz: list of coordinates

        Returns
        -------
        abc : list of fractional coordinates
        """
        xyz = np.array(xyz)
        abc = np.dot(xyz, self.reciprocal_lattice)
        return abc

    def convert_hkl_to_xyz(self, hkl):
        """
        Convert hkl miller indices to xyz coordinates.

        Parameters
        ----------
        hkl : list of hkl coordinates

        Returns
        -------
        xyz : list of xyz coordinates
        """
        hkl = np.array(hkl)
        xyz = np.dot(hkl, self.reciprocal_lattice.T)
        return xyz

    def convert_abc_to_xyz(self, abc):
        """
        Convert abc coordinates to xyz coordinates.

        Parameters
        ----------
        abc : list of fractional coordinates

        Returns
        -------
        xyz : list of coordinates

        """
        abc = np.array(abc)
        xyz = np.dot(abc, self.lattice)
        return xyz

    def convert_hkl_to_xyz2(self, hkl):
        """
        Convert hkl coordinates to xyz coordinates (written as a check on convert_hkl_to_xyz)

        Parameters
        ----------
        hkl : list of hkl miller indices

        Returns
        -------
        xyz : list of xyz coordinates

        """
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
        """
        Set the fractional coordinates and calculate the xyz coordinates.

        Parameters
        ----------
        coords : a list of fractional coordinates

        Returns
        -------
        None
        """
        self.fractional_coordinates = coords
        self.xyz_coordinates = self.convert_abc_to_xyz(coords)
        self.nions = len(coords)
        return

    def set_xyz_coordinates(self, coords, units='Angstrom'):
        """
        Set the xyz coordinates and calculate the fractional coordinates.

        Parameters
        ----------
        coords : a list of xyz coordinates
            A list of xyz coordinates, the unit of length must agree with the lattice
        units : str
            A unit of length for the input values.  The default is Angstrom.

        Returns
        -------
        None
        """
        self.xyz_coordinates = convert_length_units(coords,units,self.units)
        self.fractional_coordinates = self.convert_xyz_to_abc(coords)
        self.nions = len(coords)
        return

    def set_element_names(self, element_names):
        """
        Set the names of elements in the object after cleaning them up.

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

    def set_atom_labels(self, atom_labels):
        """
        Set the atom labels for a molecule.

        Parameters
        ----------
        atom_labels : list
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

    def find_symmetry(self,symprec=1e-5,angle_tolerance=-1.0):
        """
        Find the space group symmetry of the unit cell.

        Parameters
        ----------
        None

        Returns
        -------
        a tuple
            (international table symbol, number])
        """
        atomic_nos = self.get_atomic_numbers()
        if len(atomic_nos) <= 0:
            return 'P 1', 1
        from spglib import get_spacegroup
        cell = ( self.lattice, self.fractional_coordinates, atomic_nos )
        spacegroup = get_spacegroup(cell, symprec=symprec,angle_tolerance=angle_tolerance)
        sp = spacegroup.split()
        symbol = sp[0]
        number = int(sp[1].replace('(','').replace(')',''))
        return symbol,number

    def calculate_molecular_contents(self, scale, toler, covalent_radii):
        """
        Finds whole molecules in the unit cell.

        Does this by creating a supercell and exploring adjacent cells to see if there is any bonding to atoms in the adjacent cell
        A new unit cell is created containing whole molecules, the order of the atoms in the new cell is different.
        The routine returns the new unit cell, the number of molecules and the old ordering of the atoms.

        Parameters
        ----------
        scale : float
            The scale factor applied to the covalent radii.
        toler : float
            The tolerance added to the scaled sum of radii to determine the maximum allowable distance between atoms i and j for them to be considered bonded.
        covalent_radii : list of float
            A list of covalent radii for the atoms

        Returns
        -------
        new_unit_cell, len(new_molecules), old_order

        Notes
        -----
        The formula used to calculate the largest distance apart atoms i and j can be for a bond is:
        `scale * (radi + radj) + toler`
        """
        #print('jk06 ')
        #self.print()
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
                #jk print('New positions ',new_position,l,i,j,k)
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
            #jk print("Element name",el)
            if covalent_radii[el] > rmax:
                rmax = covalent_radii[el]
        boxSize = 2.0*scale*rmax + 0.5 + toler
        #jk print('rmax = ',rmax)
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
            except:
                BoxAtoms[abc] = [i]
        # Calculate the neighbouring boxes for each occupied box
        for abc in BoxAtoms:
            #jk print('Box ',abc, BoxAtoms[abc])
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
            #jk print('Calculating bonding to ',i,xyzi)
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
                            dist1 = scale*( covalent_radii[i_el] + covalent_radii[j_el]) + toler
                            dist2 = calculate_distance(xyzi,xyz_supercell[j])
                            if dist2 < dist1:
                                bondedToAtom[i].append(j)
                                bondedToAtom[j].append(i)
                                #jk print('new bond', i, j, i_el, j_el, xyzi, xyz_supercell[j], dist1, dist2)
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
        # molecules is a dictionary of molecules each entry is a is a list of atoms
        molecules = {}
        molID = -1
        # We stop when all the atoms in the original cell belong to a molecule
        remainingAtoms = [ atom for atom in range(self.nions) ]
        bonds = []
        while len(remainingAtoms) > 0:
            #jk print("Remaining atoms")
            #jk print(remainingAtoms)
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
                        #jk print("Using this molecule", useThisMolecule)
                        # Go through all the atoms bonded to i and add to the current molecule
                        for j in bondedToAtom[i]:
                            jx = index_supercell[j]
                            #jk print("atom j / jx is bonded to atom i",j,jx,i)
                            # The image of j in the original cell might not be available, and j might be bonded
                            if jx in remainingAtoms and not j in belongsToMolecule:
                                # if j was not already in a molecule then we have new information
                                moreAtomsToBeFound = True
                                molecules[useThisMolecule].append(j)
                                belongsToMolecule[j] = useThisMolecule
                                #jk print("Removing atom index(j) from remaining atoms",index_supercell[j])
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
        #jk print('Number of molecules', molID+1)
        self.centres_of_mass = []
        self.total_mass = 0.0
        for mol_index in molecules:
            #jk print('Molecule ',mol_index)
            #jk print('Atoms ',molecules[mol_index])
            for atom_index in molecules[mol_index]:
                index = index_supercell[atom_index]
                #jk print('New atom index, old index', atom_index, index, self.element_names[index])
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
            #jk print('Mass', mass)
            #jk print('Centre of mass', cm_fractional)
        # Create a new unit cell with the atoms shifted so that whole molecules are ordered and within the cell
        new_molecules = []
        new_fractional = np.empty_like(self.fractional_coordinates)
        new_element_names = []
        new_index = 0
        new_masses = []
        old_order = []
        for mol_index in molecules:
            new_atom_index = []
            cm = self.centres_of_mass[mol_index]
            shift = np.array( [ 0.0, 0.0, 0.0] )
            if cm[0] < 0.0:
                shift += [ 1.0, 0.0, 0.0 ]
            elif cm[0] > 1.0:
                shift += [-1.0, 0.0, 0.0 ]
            if cm[1] < 0.0:
                shift += [ 0.0, 1.0, 0.0 ]
            elif cm[1] > 1.0:
                shift += [ 0.0,-1.0, 0.0 ]
            if cm[2] < 0.0:
                shift += [ 0.0, 0.0, 1.0 ]
            elif cm[2] > 1.0:
                shift += [ 0.0, 0.0,-1.0 ]
            for atom_index in molecules[mol_index]:
                old_index = index_supercell[atom_index]
                old_order.append(old_index)
                new_fractional[new_index] = shift + np.array(fractional_supercell[atom_index])
                new_element_names.append(self.element_names[old_index])
                new_masses.append(self.atomic_masses[old_index])
                new_atom_index.append(new_index)
                new_index += 1
            new_molecules.append(new_atom_index)
        # as well as being able to go from the new order and look up the old order
        # we need to be able to take the old order and look up what the new order is
        invert_old_order = np.zeros_like(old_order)
        for i,j in enumerate(old_order):
            invert_old_order[j] = i
        new_bonds = []
        for bond in bonds:
            i,j = bond
            ix = invert_old_order[index_supercell[i]]
            jx = invert_old_order[index_supercell[j]]
            new_bonds.append( (ix,jx) )
        new_unit_cell = UnitCell( self.a, self.b, self.c, self.alpha, self.beta, self.gamma )
        new_unit_cell.set_fractional_coordinates(new_fractional.tolist())
        new_unit_cell.set_element_names(new_element_names)
        new_unit_cell.set_atomic_masses(new_masses)
        new_unit_cell.set_molecules(new_molecules)
        new_unit_cell.set_bonds(new_bonds)
        return new_unit_cell, len(new_molecules), old_order

    def set_bonds(self, bonds):
        '''
        Define a list of bonds for the unit cell

        Parameters
        ----------
        bonds : a list of integer tuple
           The bonding list for the cell

        Returns
        -------
        None
        '''
        self.bonds = bonds
        return

    def set_molecules(self, molecules):
        '''
        Define a list of molecules, each molecule is a list of atom coordinates.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        self.molecules = molecules
        return

