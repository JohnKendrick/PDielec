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
"""Read the contents of a QE output file containing QE dynmical matrix"""
import re
import math
import numpy as np
from Python.Constants import amu, angs2bohr
from Python.UnitCell import UnitCell
from Python.GenericOutputReader import GenericOutputReader


class QEOutputReader(GenericOutputReader):
    """Read the contents of a QE output file containing QE dynmical matrix"""

    def __init__(self, filenames):
        GenericOutputReader.__init__(self, filenames)
        self.type                    = 'QE output'
        self.alat                    = None
        self.mass_species            = None
        self.magnetization           = None
        self.pressure                = None
        return

    def _read_output_files(self):
        """Read the QE file names"""
        # Define the search keys to be looked for in the files
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['header']   = (re.compile('Dynamical matrix file'), self._read_header)
        self.manage['lattice']  = (re.compile('Basis vectors'), self._read_lattice_vectors)
        self.manage['lattice2']  = (re.compile('cubic'), self._read_lattice_vectors)
        self.manage['dynamical']  = (re.compile(' *Dynamical  Matrix in c'), self._read_dynamical)
        self.manage['epsilon']  = (re.compile(' *Dielectric Tensor:'), self._read_epsilon)
        self.manage['charges']  = (re.compile(' *Effective Charges E-U:'), self._read_born_charges)
        self.manage['energy_cutoff']  = (re.compile(' *kinetic-energy cutoff'), self._read_energy_cutoff)
        self.manage['kpoints']  = (re.compile(' *number of k points'), self._read_kpoints)
        self.manage['kpoint_grid']  = (re.compile('K_POINTS automatic'), self._read_kpoint_grid)
        self.manage['electrons']  = (re.compile('^ *number of electrons'), self._read_electrons)
        self.manage['energy']  = (re.compile('^ *total energy  *='), self._read_energy)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_electrons(self, line):
        self.electrons = float(line.split()[4])

    def _read_energy(self, line):
        self.final_energy_without_entropy = float(line.split()[3]) * 27.21 / 2.0
        self.final_free_energy = float(line.split()[3]) * 27.21 / 2.0

    def _read_energy_cutoff(self, line):
        self.energy_cutoff = float(line.split()[3]) * 27.21 / 2.0

    def _read_kpoints(self, line):
        self.kpoints = int(line.split()[4])

    def _read_kpoint_grid(self, line):
        line = self.file_descriptor.readline()
        self.kpoint_grid = [ float(f) for f in line.split()[0:3] ]

    def _read_header(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.nspecies = int(line.split()[0])
        self.nions    = int(line.split()[1])
        self.alat     = float(line.split()[3])

    def _read_epsilon(self, line):
        self.file_descriptor.readline()
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric = []
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        return

    def _read_masses(self):
        self.mass_species = []
        self.species      = []
        for i in range(self.nspecies):
            linea = self.file_descriptor.readline().replace('\'', '').split()
            self.species.append(linea[1])
            # The factor of two is because au in pwscf are half mass of electron
            self.mass_species.append(float(linea[2])*2/amu)
        self._read_cartesian_coordinates()
        return

    def _read_dynamical(self, line):
        nmodes = self.nions*3
        hessian = np.zeros((nmodes, nmodes))
        self.file_descriptor.readline()
        linea = self.file_descriptor.readline().split()
        # We only want to read the hessian at gamma
        q = [float(q) for q in linea[3:6]]
        qsum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2]
        if qsum > 0.0001:
            return
        # We read the hessian and store the mass weighted matrix
        linea = self.file_descriptor.readline().split()
        for a in range(self.nions):
            for b in range(self.nions):
                self.file_descriptor.readline()
                for ixyz in range(3):
                    ipos = a*3 + ixyz
                    linea = self.file_descriptor.readline().split()
                    for jxyz in range(3):
                        jpos = b*3 + jxyz
                        # factor of 0.5 'cos of au units in pwscf
                        hessian[ipos, jpos] = 0.5*float(linea[2*jxyz])/(amu*math.sqrt(self.masses[a]*self.masses[b]))
                    # end for jxyz
                # end for ixyz
            # end for b
        # end for a
        self._dynamical_matrix(hessian)

    def _read_born_charges(self, line):
        self.born_charges = []
        line = self.file_descriptor.readline()
        for i in range(self.nions):
            b = []
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            self.born_charges.append(b)
        if self.neutral:
            self._born_charge_sum_rule()
        return

    def _read_lattice_vectors(self, line):
        linea = self.file_descriptor.readline().split()
        avector = [float(f)*self.alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        bvector = [float(f)*self.alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        cvector = [float(f)*self.alat/angs2bohr for f in linea[0:3]]
        self.unit_cells.append(UnitCell(avector, bvector, cvector))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].volume
        self._read_masses()
        return

    def _read_cartesian_coordinates(self):
        self.ions = []
        self.masses = []
        for i in range(self.nions):
            linea = self.file_descriptor.readline().split()
            species = int(linea[1])
            self.ions.append([float(linea[2]), float(linea[3]), float(linea[4])])
            self.masses.append(self.mass_species[species-1])
        return
