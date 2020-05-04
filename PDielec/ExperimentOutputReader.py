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
"""Read the contents of a directory containg Experiment input and output files"""
import re
import os
import numpy as np
from PDielec.UnitCell import UnitCell
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.Constants           import atomic_number_to_element
from PDielec.Calculator          import initialise_diagonal_tensor


class ExperimentOutputReader(GenericOutputReader):
    """Read the contents of a directory containg Experiment input and output files"""

    def __init__(self, names):
        GenericOutputReader.__init__(self, names)
        self.type                    = 'Experimental output'
        self._ion_type_index = {}
        return

    def _read_output_files(self):
        """Read the Experiment files in the directory"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['lattice']      = (re.compile('lattice'), self._read_lattice_vectors)
        self.manage['species']      = (re.compile('species'), self._read_species)
        self.manage['fractional']   = (re.compile('unitcell'), self._read_fractional_coordinates)
        self.manage['static']       = (re.compile('static'), self._read_static_dielectric)
        self.manage['frequencies']  = (re.compile('frequencies'), self._read_frequencies)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_frequencies(self, line):
        nfreq = int(line.split()[1])
        self.frequencies = []
        self.oscillator_strengths = []
        for i in range(nfreq):
            line = self.file_descriptor.readline()
            self.frequencies.append(float(line.split()[0]))
            strength = float(line.split()[1])
            self.oscillator_strengths.append(initialise_diagonal_tensor( [strength, strength, strength] ) )
        return

    def _read_species(self, line):
        nspecies = int(line.split()[1])
        self.species = []
        self.masses_per_type = []
        for i in range(nspecies):
            line = self.file_descriptor.readline()
            species = line.split()[0]
            self.species.append(species)
            self.masses_per_type.append(float(line.split()[1]))
            self._ion_type_index[species] = self.nspecies
            self.nspecies = len(self.species)
        return

    def _read_lattice_vectors(self, line):
        scalar = float(line.split()[1])
        line = self.file_descriptor.readline()
        avector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        line = self.file_descriptor.readline()
        bvector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        line = self.file_descriptor.readline()
        cvector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        cell = UnitCell(avector, bvector, cvector)
        self.unit_cells.append(cell)
        self.ncells = len(self.unit_cells)
        self.volume = cell.volume
        return

    def _read_fractional_coordinates(self, line):
        n = 0
        ions = []
        self.nions = int(line.split()[1])
        self.ions_per_type = [ 0 for s in self.species ]
        self.masses = []
        species_list = []
        for n in range(self.nions):
            line = self.file_descriptor.readline()
            species = line.split()[0]
            index = self._ion_type_index[species]
            self.atom_type_list.append(index)
            species_list.append(species)
            ions.append([float(f) for f in line.split()[1:4]])
            self.ions_per_type[index] += 1
            self.masses.append(self.masses_per_type[index])
        self.unit_cells[-1].set_fractional_coordinates(ions)
        self.unit_cells[-1].set_element_names(species_list)
        return

    def _read_static_dielectric(self, line):
        # the is epsilon infinity
        od = []
        line = self.file_descriptor.readline()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        od.append([complex(f) for f in line.split()[0:3]])
        # If we have complex input return a complex list, otherwise return a real list
        odc = np.array(od)
        odi = np.absolute(np.imag(odc))
        sumi = np.sum(odi)
        if sumi < 1.0e-12:
            odc = np.real(odc)
        self.zerof_optical_dielectric = odc.tolist()
        return

    def calculate_mass_weighted_normal_modes(self):
        print("Calculate mass weighted normal modes is not needed")
        return []
