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

"""Read contents of a directory containing Crystal input and output files"""

import re
import math
import os
import numpy as np
from PDielec.Constants import amu, hartree2ev
from PDielec.UnitCell import UnitCell
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.IO import pdielec_io


class CrystalOutputReader(GenericOutputReader):
    """Read contents of a directory containing Crystal input and output files"""

    def __init__(self, filenames):
        GenericOutputReader.__init__(self, filenames)
        self.type                    = 'Crystal output'
        self.hessian_symmetrisation  = 'crystal'
        self._fractional_coordinates = []
        return

    def _read_output_files(self):
        """Read the Crystal files in the directory"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['masses']   = (re.compile(' ATOMS ISOTOPIC MASS'), self._read_masses)
        self.manage['lattice']  = (re.compile(' DIRECT LATTICE VECTORS CARTESIAN COMPONENTS'), self._read_lattice_vectors)
        self.manage['fractional']  = (re.compile(' ATOMS IN THE ASYMMETRIC UNIT'), self._read_fractional_coordinates)
        self.manage['bornCharges']  = (re.compile(' ATOMIC BORN CHARGE TENSOR'), self._read_born_charges)
        self.manage['eigenvectors']  = (re.compile(' NORMAL MODES NORMALIZ'), self._read_eigenvectors)
        self.manage['staticIonic']  = (re.compile(' SUM TENSOR OF THE VIBRATIONAL CONTRIBUTIONS TO '), self._read_ionic_dielectric)
        self.manage['epsilon']  = (re.compile(' SUSCEPTIBILITY '), self._read_epsilon)
        self.manage['kpoints']  = (re.compile(' SHRINK\. FACT\.\('), self._read_kpoints)
        self.manage['electrons']  = (re.compile(' N\. OF ELECTRONS'), self._read_electrons)
        self.manage['energy']  = (re.compile(' TOTAL ENERGY\(DFT\)'), self._read_energy)
        self.manage['energy2']  = (re.compile(' TOTAL ENERGY + DISP'), self._read_energy2)
        self.manage['energy3']  = (re.compile(' *CENTRAL POINT'), self._read_energy3)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_energy(self, line):
        self.final_free_energy = hartree2ev*float(line.split()[3])
        self.final_energy_without_entropy = hartree2ev*float(line.split()[3])

    def _read_energy2(self, line):
        self.final_free_energy = hartree2ev*float(line.split()[5])
        self.final_energy_without_entropy = hartree2ev*float(line.split()[5])

    def _read_energy3(self, line):
        self.final_free_energy = hartree2ev*float(line.split()[2])
        self.final_energy_without_entropy = hartree2ev*float(line.split()[2])

    def _read_electrons(self, line):
        self.electrons = int(line.split()[5])

    def _read_kpoints(self, line):
        self.kpoint_grid = [ int(line.split()[2]), int(line.split()[3]), int(line.split()[4]) ]
        self.kpoints = int(line.split()[12])

    def _read_epsilon(self, line):
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        optical_dielectric = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        linea = self.file_descriptor.readline().split()
        while len(linea) > 0:
            component = linea[0]
            # Problems with Crystal17 file format
            if component == "0.000":
                component = linea[1]
                epsilon = float(linea[4])
            else:
                epsilon = float(linea[3])
            if component == "XX": optical_dielectric[0][0] = epsilon
            if component == "YY": optical_dielectric[1][1] = epsilon
            if component == "ZZ": optical_dielectric[2][2] = epsilon
            if component == "XY": optical_dielectric[0][1] = epsilon
            if component == "XZ": optical_dielectric[0][2] = epsilon
            if component == "YZ": optical_dielectric[1][2] = epsilon
            # read the next line
            linea = self.file_descriptor.readline().split()
        # symmetrize
        optical_dielectric[1][0] = optical_dielectric[0][1]
        optical_dielectric[2][0] = optical_dielectric[0][2]
        optical_dielectric[2][1] = optical_dielectric[1][2]
        self.zerof_optical_dielectric = optical_dielectric
        return

    def _read_masses(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        n = 2
        self.masses_per_type = [ 0 for i in range(self.nspecies) ]
        for i in range(self.nions):
            if n > 11:
                line = self.file_descriptor.readline()
                n = 2
            mass = float(line.split()[n])
            self.masses.append(mass)
            self.masses_per_type[self.atom_type_list[i]] = mass
            n = n + 3
        return

    def _read_eigenvectors(self, line):
        if os.path.isfile(self.open_directory+"/HESSFREQ.DAT"):
            # print("Reading Hessian from HESSFREQ.DAT",file=sys.stderr)
            self._read_hessfreq_dat(self.open_directory+"/HESSFREQ.DAT")
        else:
            # print("Reading Normal Modes from output file",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            # print("The numerical precision of the input is limited",file=sys.stderr)
            # print("If possible proved a HESSFREQ.DAT from the calculation",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            self._read_output_eigenvectors(line)
        return

    def _read_hessfreq_dat(self,filename):
        fd2 = pdielec_io(filename, 'r')
        nmodes = self.nions*3
        # Create a mass weighting vector
        n = 0
        massweight = np.zeros(nmodes)
        for a in range(self.nions):
            for j in range(3):
                massweight[n] = 1.0 / math.sqrt(self.masses[a]*amu)
                n = n + 1
        # We read the hessian and store the mass weighted matrix
        hessian = np.zeros((nmodes, nmodes))
        line = fd2.readline()
        pos = 0
        for i in range(nmodes):
            for j in range(nmodes):
                hessian[i, j] = massweight[i]*massweight[j]*float(line.split()[pos])
                pos = pos + 1
                if pos >= 4:
                    line = fd2.readline()
                    pos = 0
                # end if pos
            # end for j
        # end for i
        fd2.close()
        # symmetrise, project, diagonalise and store the frequencies and normal modes
        self._dynamical_matrix(hessian)
        return

    def _read_output_eigenvectors(self, line):
        self.file_descriptor.readline()
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        n = 6
        nmodes = 3 * self.nions
        self.frequencies = []
        # Loop over all the modes
        for j in range(nmodes):
            # They come in columns of 6 so after 6 modes we read the frequencies and initialse the 6 new modes
            n = n + 1
            if n >= 6:
                n = 0
                linef = self.file_descriptor.readline().split()[1:]
                linef = [float(f) for f in linef]
                self.file_descriptor.readline()
                atoms = []
                for f in linef:
                    self.frequencies.append(f)
                    atoms.append([])
                # Read through the XYZ components for each atom and store the mode in atoms
                for i in range(self.nions):
                    mass = np.sqrt(self.masses[i])
                    linex = self.file_descriptor.readline().split()[4:]
                    liney = self.file_descriptor.readline().split()[1:]
                    linez = self.file_descriptor.readline().split()[1:]
                    for a, x, y, z in zip(atoms, linex, liney, linez):
                        x = float(x) * mass
                        y = float(y) * mass
                        z = float(z) * mass
                        a.append([x, y, z])
                for a in atoms:
                    self.mass_weighted_normal_modes.append(a)
                self.file_descriptor.readline()
            # end of if n >= 6
        for i, mode in enumerate(self.mass_weighted_normal_modes):
            marray = np.array(mode)
            sumxyz = 0.0
            for atoms in marray:
                for xyz in atoms:
                    sumxyz = sumxyz + xyz*xyz
            marray = marray / np.sqrt(sumxyz)
            self.mass_weighted_normal_modes[i] = marray.tolist()
        return

    def _read_born_charges(self, line):
        """Read the born charges from the outputfile file.
           Each row of the output refers to a given field direction
           Each column in the row refers the atomic displacement
           so the output is arranged [[a1x a1y a1z]
                                      [a2x a2y a2z]
                                      [a3x a3y a3z]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        # Extract directory containing the output file
        if os.path.isfile(self.open_directory+"/BORN.DAT"):
            # print("Reading Born charge tensor from BORN.DAT",file=sys.stderr)
            self._read_born_charges_from_born_dat(self.open_directory+"/BORN.DAT")
        else:
            # print("Reading Born Charge Tensor from output file",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            # print("The numerical precision of the input is limited",file=sys.stderr)
            # print("If possible proved a BORN.DAT from the calculation",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            self._read_born_charges_from_output(line)
        return

    def _read_born_charges_from_born_dat(self,filename):
        fd2 = pdielec_io(filename, 'r')
        self.born_charges = []
        for i in range(self.nions):
            b = []
            line = fd2.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = fd2.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = fd2.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            b_np = np.array(b)
            bt_np = b_np.T
            b = bt_np.tolist()
            self.born_charges.append(b)
        fd2.close()
        return

    def _read_born_charges_from_output(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.born_charges = []
        for i in range(self.nions):
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            b = []
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            self.born_charges.append(b)
            line = self.file_descriptor.readline()
        return

    def _read_ionic_dielectric(self, line):
        # Read the ionic contribution to the static dielectric
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric = []
        self.zerof_static_dielectric.append([float(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric.append([float(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric.append([float(f) for f in line.split()[0:3]])
        # a = np.array(ionic_dielectric)
        # self.zerof_static_dielectric = a.tolist()
        return

    def _read_lattice_vectors(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        avector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        bvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        cvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        self.unit_cells.append(UnitCell(avector, bvector, cvector))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].volume
        # The fractional coordinates are specified before the lattice vectors
        self.unit_cells[-1].set_fractional_coordinates(self._fractional_coordinates)
        self.unit_cells[-1].set_element_names(self.species_list)
        return

    def _read_fractional_coordinates(self, line):
        self.nions = int(line.split()[12])
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self._fractional_coordinates = []
        self.species = []
        self.species_list = []
        self.atom_type_list = []
        for i in range(self.nions):
            line = self.file_descriptor.readline()
            species = line.split()[3].capitalize()
            self.species_list.append(species)
            if species not in self.species:
                self.species.append(species)
            species_index = self.species.index(species)
            self.atom_type_list.append(species_index)
            self._fractional_coordinates.append([float(line.split()[4]), float(line.split()[5]), float(line.split()[6])])
        self.nspecies = len(self.species)
        self.ions_per_type = [ 0 for i in range(self.nspecies) ]
        for species_index in self.atom_type_list:
            self.ions_per_type[species_index] += 1
        return
