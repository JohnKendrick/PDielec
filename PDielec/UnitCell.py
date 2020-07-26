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
import numpy as np
import math
from PDielec.Calculator  import calculate_distance
from PDielec.Plotter     import print_reals, print_ints, print_strings
from PDielec.Calculator  import cleanup_symbol
import sys


class UnitCell:
    """Hold unit cell information and its associated calculated properties"""
    def __init__(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None):
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
        if gamma is not None:
            self.lattice = self.convert_abc_to_unitcell(a, b, c, alpha, beta, gamma)
        else:
            self.lattice[0] = a
            self.lattice[1] = b
            self.lattice[2] = c
        self._calculate_reciprocal_lattice(self.lattice)

    def write_cif(self, filename='', file_=sys.stdout):
        print('data_'+filename,                                  file=file_)
        print('_space_group_IT_number 1',                        file=file_)
        print('_symmetry_space_group_name_H-M \'P 1\'',          file=file_)
        print('_symmetry_Int_Tables_number        1  ',          file=file_)
        print('_cell_length_a      {:12.6f}'.format(self.a),     file=file_)
        print('_cell_length_b      {:12.6f}'.format(self.b),     file=file_)
        print('_cell_length_c      {:12.6f}'.format(self.c),     file=file_)
        print('_cell_angle_alpha   {:12.6f}'.format(self.alpha), file=file_)
        print('_cell_angle_beta    {:12.6f}'.format(self.beta),  file=file_)
        print('_cell_angle_gamma   {:12.6f}'.format(self.gamma), file=file_)
        print('_cell_volume        {:12.6f}'.format(self.volume),file=file_)
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

    def getBoundingBox(self):
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
        return corners_xyz,edges


    def print_info(self):
        print_reals('Unit Cell a,b,c ',[self.a, self.b, self.c], format='{:12.6f}')
        print_reals('Unit Cell alpha,beta,gamma',[self.alpha, self.beta, self.gamma], format='{:12.6f}')
        print_reals('lattice', self.lattice[0], format='{:12.6f}')
        print_reals('', self.lattice[1], format='{:12.6f}')
        print_reals('', self.lattice[2], format='{:12.6f}')
        print_strings('Element names',self.element_names)
        print_reals('Fractional coords',self.fractional_coordinates[0], format='{:12.6f}')
        for frac in self.fractional_coordinates[1:]:
            print_reals('',frac, format='{:12.6f}')
        print_reals('Cartesian coords',self.xyz_coordinates[0], format='{:12.6f}')
        for xyz in self.xyz_coordinates[1:]:
            print_reals('',xyz, format='{:12.6f}')
        if self.molecules:
            for molid,atoms in enumerate(self.molecules):
                mass, cm_xyz, cm_frac = self.calculateCentreOfMass(atom_list=atoms,units='all')
                molstring = 'Molecule '+str(molid)+':'
                print_ints('Atoms in '+molstring,atoms)
                print_reals('Mass of '+molstring,[ mass ], format='{:12.6f}')
                print_reals('Centre of Mass  (xyz) of '+molstring, cm_xyz, format='{:12.6f}')
                print_reals('Centre of Mass (frac) of '+molstring, cm_frac,format='{:12.6f}')

    def calculateCentreOfMass(self,atom_list=None, units='xyz'):
        # Calculate the centre of mass - if the atom list is given just use that
        # The centre of mass can be returned in units of 'xyz' space or 'abc' space
        # if units='all' a tuple of (mass,cm_xyz,cm_abc) is returned
        if atom_list == None:
            atom_list = range(self.nions)
        mass = 0.0
        cm = np.zeros(3)
        for atom_index in atom_list:
            mass += self.atomic_masses[atom_index]
            cm = cm + self.atomic_masses[atom_index] * self.xyz_coordinates[atom_index]
        cm_xyz = cm / mass
        if units == 'xyz':
            return cm_xyz
        elif units == 'mass':
            return mass
        elif units == 'abc':
            cm_fractional = self.convert_xyz_to_abc(cm_xyz)
            return cm_fractional
        else:
            cm_fractional = self.convert_xyz_to_abc(cm_xyz)
            return mass, cm_xyz, cm_fractional

    def set_atomic_masses(self, masses):
        self.atomic_masses = masses

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
        lattice[1] = [0.0, b*np.sin(alpha), b*np.cos(alpha)]
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
        self.xyz_coordinates = self.convert_abc_to_xyz(coords)
        self.nions = len(coords)
        return

    def set_xyz_coordinates(self, coords):
        """Set the xyz coordinates and calculate the fractional coordinates"""
        self.xyz_coordinates = coords
        self.fractional_coordinates = self.convert_xyz_to_abc(coords)
        self.nions = len(coords)
        return

    def set_element_names(self, element_names):
        self.element_names = []
        for el in element_names:
          el = cleanup_symbol(el)
          self.element_names.append(el)
        return

    def set_atom_labels(self, atom_labels):
        self.atom_labels = []
        for el in atom_labels:
          self.atom_labels.append(el)
        return

    def find_symmetry(self):
        """Find the space group symmetry of the unit cell"""
        from spglib import get_spacegroup
        cell = ( self.lattice, self.fractional_coordinates, self.atomic_nos )
        self.spacegroup = get_spacegroup(cell, symmprec=1e-5)
        print("Symmetry space group is", self.spacegroup)

    def calculate_molecular_contents(self, scale, toler, covalent_radii):
        """Find whole molecules in the unit cell
           scale and toler are used to manipulate the covalent radii, thus the largest distance
           apart atoms i and j can be for a bond is
           scale*( radi + radj) + toler """
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
        '''Define a list of bonds for the unit cell'''
        self.bonds = bonds
        return

    def set_molecules(self, molecules):
        '''Define a list of molecules, each molecule is a list of atom coordinates'''
        self.molecules = molecules
        return

