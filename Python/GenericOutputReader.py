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
"""Generic reader of output files.  Actual reader should inherit from this class"""
from __future__ import print_function
import math
import os
import sys
import string
import numpy as np
from Python.Constants import wavenumber, avogadro_si, atomic_number_to_element
from Python.Plotter import print3x3, print_reals, print_strings, print_ints

class GenericOutputReader:
    """Generic reader of output files.  Actual reader should inherit from this class"""

    def __init__(self, filenames):
        self._outputfiles               = filenames
        self.names                      = [os.path.abspath(f) for f in filenames]
        self.debug                      = False
        self.type                       = 'Unkown'
        self.ncells                     = 0
        self.unit_cells                 = []
        self.nsteps                     = 0
        self.electrons                  = 0
        self.spin                       = 0
        self.nbands                     = 0
        self.volume                     = 0
        self.nions                      = 0
        self.nspecies                   = 0
        self.species                    = []
        self.final_free_energy          = 0.0
        self.final_energy_without_entropy = 0.0
        self.kpoints                    = 1
        self.kpoint_grid                = [ 1, 1, 1 ]
        self.energy_cutoff              = 0.0
        self.born_charges               = []
        self.manage                     = {}
        self.iterations                 = {}
        self.file_descriptor            = ''
        self.pressure                   = 0
        self.magnetization              = 0.0
        # this in epsilon infinity
        self.zerof_optical_dielectric   = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric    = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        self.elastic_constants          = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
        self.frequencies                = []
        self.mass_weighted_normal_modes = []
        self.ions_per_type              = []
        self.atom_type_list             = []
        self.masses                     = []
        self.masses_per_type            = []
        self.eckart                     = False
        self.neutral                    = False
        self.hessian_symmetrisation     = "symm"
        self.open_filename              = ""
        self.open_directory             = ""
        self._old_masses                = []
        self._masses_have_changed       = False
        return

    def read_output(self):
        """Interface to the private read output files methods"""
        self._read_output_files()
        return

    def cleanup_symbol(self,s):
        """Return a true element from the symbol"""
        s = s.capitalize()
        s = s.replace('_','')
        for i in string.digits:
            s = s.replace(i,'')
        return s

    def change_masses(self, definition, isotope_masses, average_masses, mass_dictionary):
        self._old_masses = self.masses
        old_masses_per_type = self.masses_per_type
        self._masses_have_changed = False
        self.masses = []
        self.masses_per_type = []
        if definition == "average":
            self._masses_have_changed = True
            for symbol in self.species:
                # the element name may be appended with a digit or an underscore
                element = self.cleanup_symbol(symbol)
                mass = average_masses[element]
                if element in mass_dictionary:
                    mass = mass_dictionary[element]
                self.masses_per_type.append(mass)
             # end for symbol
        elif definition == "isotope":
            self._masses_have_changed = True
            for symbol in self.species:
                # the element name may be appended with a digit or an underscore
                element = self.cleanup_symbol(symbol)
                abmax = 0.0
                for iso in isotope_masses[element]:
                    weight = iso[1]
                    abundance = iso[2]
                    if abundance > abmax:
                        abmax = abundance
                        mass = weight
                # end of for iso 
                if element in mass_dictionary:
                    mass = mass_dictionary[element]
                self.masses_per_type.append(mass)
            # end of for symbol
        elif definition == "program":
            for symbol,mass in zip(self.species,old_masses_per_type):
                # the element name may be appended with a digit or an underscore
                element = self.cleanup_symbol(symbol)
                if element in mass_dictionary:
                    mass = mass_dictionary[element]
                    self._masses_have_changed = True
                self.masses_per_type.append(mass)
            # end of for symbol
        else:
            print('Error, mass definition not recognised',definition)
            exit()
        self.masses = [ self.masses_per_type[atype] for atype in self.atom_type_list ]
        return

    def print_info(self):
        """Print information about the system"""
        # Generic printing of information
        print("")
        print("Summary of information contained in the QM/MM Reader")
        print("")
        print("Number of atoms: {:5d}".format(self.nions))
        print("")
        print("Number of species: {:5d}".format(self.nspecies))
        print_strings("Species:", self.species)
        print_ints("Number of atoms for each species:", self.ions_per_type)
        print_reals("Mass of each species:", self.masses_per_type,format="{:10.6f}")
        print_ints("Atom type list:", self.atom_type_list)
        print("")
        print("Number of kpoints: {:5d}".format(self.kpoints))
        print("")
        print("Kpoint grid      : {:5d} {:5d} {:5d}".format(self.kpoint_grid[0], self.kpoint_grid[1], self.kpoint_grid[2]))
        print("")
        print("Energy cutoff (eV): {:f}".format(self.energy_cutoff))
        print("")
        print_reals("Frequencies (cm-1):", self.frequencies)
        print_reals("Masses (amu):", self.masses,format="{:10.6f}")
        for i, charges in enumerate(self.born_charges):
            title = "Born Charges for Atom {:d}".format(i)
            print3x3(title, charges)
        print3x3("Epsilon inf: ", self.zerof_optical_dielectric)
        print3x3("Unit cell: ", self.unit_cells[-1].lattice)
        print(" ")
        print("Volume of cell: {:f}".format(self.volume))
        mtotal = 0.0
        for m in self.masses:
            mtotal = mtotal + m
        print("Total mass is: {:f} g/mol".format(mtotal))
        print("Density is: {:f} g/cc".format(mtotal/(avogadro_si * self.volume * 1.0e-24)))
        print(" ")
        return

    def _read_output_files(self):
        """Read the file names"""
        # Define the search keys to be looked for in the files
        print("Error _read_output_files must be defined by the actual file reader")
        return

    def _read_output_file(self, name):
        """Read through the files for key words.  The keywords are established in _read_output_files"""
        # Check to see if the file exists....
        if not os.path.isfile(name):
            print("Warning file is not present: ", name, file=sys.stderr)
            return
        # Open file and store file name and directory
        self.file_descriptor = open(name, 'r')
        self.open_filename = name
        self.open_directory = os.path.dirname(name)
        if self.open_directory == "":
            self.open_directory = "."
        # Loop through the contents of the file a line at a time and parse the contents
        line = self.file_descriptor.readline()
        while line != '':
            for k in self.manage:
                if self.manage[k][0].match(line):
                    method   = self.manage[k][1]
                    if self.debug:
                        print('_read_output_file: Match found {}'.format(k))
                    method(line)
                    break
                # end if
            # end for
            line = self.file_descriptor.readline()
        # end while
        self.file_descriptor.close()
        return

    def _symmetric_orthogonalisation(self, A):
        """Private routine to perform symetric orthogonalisation"""
        # The matrix A is only approximately orthogonal
        n = np.size(A, 0)
        unity = np.eye(n)
        Ak = A
        for k in range(3):
            Bk = np.dot(Ak, Ak.T)
            error = unity - Bk
            Ck = np.linalg.inv(unity + Bk)
            Kk = np.dot(error, Ck)
            Ak = np.dot((unity + Kk), Ak)
            error  = np.sum(np.abs(error))
            if self.debug:
                print("Orthogonalisation iteration: ", error)
        # end for k
        return Ak

    def calculate_mass_weighted_normal_modes(self):
        """Calculate the mass weighted normal modes"""
        #
        # Reconstruct the massweighted hessian
        # If necessary diagonalise it and extract the frequencies and normal modes
        # Any projection can be performed here
        #
        # First step is to reconstruct the dynamical matrix (D) from the frequencies and the eigenvectors
        # f^2 = UT . D . U
        # and U is a hermitian matrix so U-1 = UT
        # D = (UT)-1 f^2 U-1 = U f UT
        # Construct UT from the normal modes
        n = np.size(self.mass_weighted_normal_modes, 0)
        m = np.size(self.mass_weighted_normal_modes, 1)*3
        nmodes = 3*self.nions
        UT = np.zeros((n, m))
        for imode, mode in enumerate(self.mass_weighted_normal_modes):
            n = 0
            for atom in mode:
                # in python the first index is the row of the matrix, the second is the column
                UT[imode, n+0] = atom[0]
                UT[imode, n+1] = atom[1]
                UT[imode, n+2] = atom[2]
                n = n + 3
            # end for atom
        # end for imode
        # convert the frequencies^2 to a real diagonal array
        # Warning we have to make sure the sign is correct here
        # The convention is that if the frequency is negative
        # then it is really imaginary, so the square of the frequency
        # will be negative too.
        frequencies_a = np.array(self.frequencies)
        f2 = np.diag(np.sign(frequencies_a)*np.real(frequencies_a*frequencies_a))
        # The back transformation uses approximately orthogonal (unitary) matrices because of rounding issues on reading vectors
        # So before that lets orthogonalise them
        UT = self._symmetric_orthogonalisation(UT)
        hessian = np.dot(np.dot(UT.T, f2), UT)
        # Make sure the dynamical matrix is real
        hessian = np.real(hessian)
        # If the masses have been changed then alter the mass weighted hessian here
        if self._masses_have_changed: 
            hessian = self._modify_mass_weighting(hessian,self._old_masses,self.masses)
            self._masses_have_changed = False
        # Project out the translational modes if requested
        if self.eckart:
            hessian = self.project(hessian)
        # Find its eigenvalues and eigen vectors
        eig_val, eig_vec = np.linalg.eigh(hessian)
        self.mass_weighted_normal_modes = []
        # Store the new frequencies, using the negative convention for imaginary modes
        for i in range(nmodes):
            if eig_val[i] < 0:
                frequencies_a[i] = -math.sqrt(-eig_val[i])
            else:
                frequencies_a[i] = math.sqrt(eig_val[i])
        self.frequencies = frequencies_a.tolist()
        # Store the mass weighted normal modes
        for i in range(nmodes):
            mode = []
            n = 0
            for j in range(self.nions):
                modea = [eig_vec[n][i], eig_vec[n+1][i], eig_vec[n+2][i]]
                n = n + 3
                mode.append(modea)
            self.mass_weighted_normal_modes.append(mode)
        # end for i
        return self.mass_weighted_normal_modes

    def project(self, hessian):
        """Apply projection operators to remove translation
        Take the given matrix (np.array)
        Project out the translational modes"""
        #
        nmodes = self.nions*3
        unit = np.eye(nmodes)
        p1 = np.zeros(nmodes)
        p2 = np.zeros(nmodes)
        p3 = np.zeros(nmodes)
        for i in range(self.nions):
            p1[i*3+0] = math.sqrt(self.masses[i])
            p2[i*3+1] = math.sqrt(self.masses[i])
            p3[i*3+2] = math.sqrt(self.masses[i])
        # end for i
        # Normalise
        p1 = p1 / math.sqrt(np.dot(p1, p1))
        p2 = p2 / math.sqrt(np.dot(p2, p2))
        p3 = p3 / math.sqrt(np.dot(p3, p3))
        # Form the projection operators
        P1 = unit - np.outer(p1, p1)
        P2 = unit - np.outer(p2, p2)
        P3 = unit - np.outer(p3, p3)
        # Now project out
        hessian = np.dot(np.dot(P1.T, hessian), P1)
        hessian = np.dot(np.dot(P2.T, hessian), P2)
        hessian = np.dot(np.dot(P3.T, hessian), P3)
        return hessian

    def _read_till_phrase(self,phrase):
        """Read lines from the current file until a match with phrase is found
           Once a match is found return the matching line"""
        line = ""
        while not phrase.match(line):
            line = self.file_descriptor.readline()
        return line

    def _dynamical_matrix(self, hessian):
        """Process the dynamical matrix
        Hessian is a nxn matrix of the mass weighted force constants
        The hessian is symmetrised
        Translational modes are projected out
        The hessian is diagonalised
        Finally the frequencies and normal modes are stored"""
        #
        nmodes = self.nions*3
        # symmetrise
        uplo = 'U'
        if self.hessian_symmetrisation == "symm":
            hessian = 0.5 * (hessian + hessian.T)
            uplo = 'L'
        else:
            uplo = 'U'
        # Project out the translational modes if requested
        if self.eckart:
            hessian = self.project(hessian)
        #
        # If the masses have been changed then alter the mass weighted hessian here
        if self._masses_have_changed: 
            hessian = self._modify_mass_weighting(hessian,self._old_masses,self.masses)
            self._masses_have_changed = False
        # diagonalise
        eig_val, eig_vec = np.linalg.eigh(hessian, UPLO=uplo)
        #
        # If eig_val has negative values then we store the negative frequency
        # convert to cm-1
        frequencies_a = np.zeros(nmodes)
        for i, eig in enumerate(eig_val):
            if eig < 0:
                frequencies_a[i] = -math.sqrt(-eig) / wavenumber
            else:
                frequencies_a[i] = math.sqrt(eig) / wavenumber
            # end if
        # end for
        self.mass_weighted_normal_modes = []
        self.frequencies = frequencies_a.tolist()
        # Store the mass weighted normal modes
        for i in range(nmodes):
            mode = []
            n = 0
            for j in range(self.nions):
                modea = [eig_vec[n][i], eig_vec[n+1][i], eig_vec[n+2][i]]
                n = n + 3
                mode.append(modea)
            self.mass_weighted_normal_modes.append(mode)
        # end for i
        return

    def _born_charge_sum_rule(self):
        """Apply a simple charge sum rule to all the elements of the born matrices"""
        total = np.zeros((3, 3))
        born_charges = np.array(self.born_charges)
        new_born_charges = np.zeros_like(self.born_charges)
        total = np.sum(born_charges) / self.nions
        new_born_charges = born_charges - total
        self.born_charges = new_born_charges.tolist()
        return

    def _modify_mass_weighting(self,hessian,old,new):
        #
        # Remove the mass weighting imposed by the QM/MM program - defined in old
        # Replace by PDielec mass weighting - defined in new
        #
        ipos = -1
        for i in range(self.nions):
            for ix in range(3):
                ipos += 1
                jpos = -1
                for j in range(self.nions):
                    for jx in range(3):
                        jpos += 1
                        hessian[ipos,jpos] = hessian[ipos,jpos] * math.sqrt( old[i]*old[j] / ( new[i]*new[j] ))
                    # end for jx
                # end for j
            # end for ix
        # end for i
        return hessian
