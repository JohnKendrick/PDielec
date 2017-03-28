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
import numpy as np
from Python.Constants import wavenumber, avogadro_si
from Python.Plotter import print3x3, print_reals


class GenericOutputReader:
    """Generic reader of output files.  Actual reader should inherit from this class"""

    def __init__(self, filenames):
        self._outputfiles               = filenames
        self.names                      = [os.path.abspath(f) for f in filenames]
        self.debug                      = False
        self.type                       = 'Unkown'
        self.ncells                     = None
        self.nsteps                     = None
        self.formula                    = None
        self.nelect                     = None
        self.volume                     = None
        self.nions                      = None
        self.nspecies                   = None
        self.final_free_energy          = None
        self.kpoints                    = None
        self.kpoint_grid                = [ None, None, 0 ]
        self.energy_cutoff              = None
        self.unit_cells                  = []
        self.born_charges               = []
        self.species                    = []
        self.manage                     = {}
        self.file_descriptor            = None
        # this in epsilon infinity
        self.zerof_optical_dielectric   = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric    = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        self.elastic_constants          = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
        self.frequencies                = []
        self.mass_weighted_normal_modes = []
        self.masses                     = []
        self.ions                       = []
        self.eckart                     = False
        self.neutral                    = False
        self.hessian_symmetrisation    = "symm"
        return

    def read_output(self):
        """Interface to the private read output files methods"""
        self._read_output_files()
        return

    def print_info(self):
        """Print information about the system"""
        # Generic printing of information
        print("Number of atoms: {:5d}".format(self.nions))
        print("Number of species: {:5d}".format(self.nspecies))
        print("Number of kpoints: {:5d}".format(self.kpoints))
        print("Kpoint grid      : {:5d} {:5d} {:5d}".format(self.kpoint_grid[0], self.kpoint_grid[1], self.kpoint_grid[2]))
        print("Energy cutoff (eV): {:f}".format(self.energy_cutoff))
        print_reals("Frequencies (cm-1):", self.frequencies)
        print_reals("Masses (amu):", self.masses)
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
        self.file_descriptor = open(name, 'r')
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
        # The back transformation uses approximately orthogonal (unitary) matrices
        # So before that lets orthogonalise them
        UT = self._symmetric_orthogonalisation(UT)
        hessian = np.dot(np.dot(UT.T, f2), UT)
        # Make sure the dynamical matrix is real
        hessian = np.real(hessian)
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
