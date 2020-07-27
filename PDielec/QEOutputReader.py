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
from PDielec.Constants import amu, angs2bohr, hartree2ev
from PDielec.UnitCell import UnitCell
from PDielec.GenericOutputReader import GenericOutputReader


class QEOutputReader(GenericOutputReader):
    """Read the contents of a QE output file containing QE dynmical matrix"""

    def __init__(self, filenames):
        GenericOutputReader.__init__(self, filenames)
        self.type                    = 'QE output'
        self._alat                    = None
        return

    def _read_output_files(self):
        """Read the QE file names"""
        # Define the search keys to be looked for in the files
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['header']   = (re.compile('Dynamical matrix file'), self._read_header)
        self.manage['lattice']  = (re.compile('Basis vectors'), self._read_lattice_vectors)
        self.manage['lattice2']  = (re.compile('cubic'), self._read_lattice_vectors)
        self.manage['lattice3']  = (re.compile('CELL_PARAMETERS'), self._read_cell_parameters)
        self.manage['positions']  = (re.compile('ATOMIC_POSITIONS'), self._read_fractional_coordinates)
        self.manage['dynamical']  = (re.compile(' *Dynamical  Matrix in c'), self._read_dynamical)
        self.manage['epsilon']  = (re.compile(' *Dielectric Tensor:'), self._read_epsilon)
        self.manage['charges']  = (re.compile(' *Effective Charges E-U:'), self._read_born_charges)
        self.manage['energy_cutoff']  = (re.compile(' *kinetic-energy cutoff'), self._read_energy_cutoff)
        self.manage['kpoints']  = (re.compile(' *number of k points'), self._read_kpoints)
        self.manage['kpoint_grid']  = (re.compile('K_POINTS automatic'), self._read_kpoint_grid)
        self.manage['electrons']  = (re.compile('^ *number of electrons'), self._read_electrons)
        self.manage['energy']  = (re.compile('^ *total energy  *='), self._read_energy)
        self.manage['alat']  = (re.compile('^ *lattice parameter'), self._read_alat)
        self.manage['alat']  = (re.compile('^ *A = '), self._read_alat2)
        self.manage['celldm1']  = (re.compile('^ *celldm.1. ='), self._read_celldm1)
        self.manage['pressure']  = (re.compile('^ *total *stress *.Ry'), self._read_pressure)
        self.manage['nions']  = (re.compile('^ *number of atoms/cell'), self._read_nions)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_nions(self, line):
        #
        # Read from log file
        #
        self.nions = int(line.split()[4])
        if self.debug:
            print('_read_nions: nions={}'.format(self.nions) )
        return

    def _read_pressure(self, line):
        self.pressure = float(line.split()[5])/10.0
        if self.debug:
            print('_read_pressure: pressure={}'.format(self.pressure) )
        return

    def _read_celldm1(self, line):
        #
        # Read from log file
        #
        string = line.split()[2]
        string = string.replace(',','')
        t = float(string)
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        # There are rounding errors when reading from the log file
        # So only read if there is no alternative
        if self._alat is None:
            self._alat = t
        if self.debug:
            print('_read_celldm1: _alat={}'.format(self._alat) )
        return

    def _read_alat(self, line):
        #
        # Read from log file
        #
        t = float(line.split()[4])
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        # There are rounding errors when reading from the log file
        # So only read if there is no alternative
        if self._alat is None:
            self._alat = t
        if self.debug:
            print('_read_alat: _alat={}'.format(self._alat) )
        return

    def _read_alat2(self, line):
        #
        # Read from log file
        #
        t = float(line.split()[2])
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        # There are rounding errors when reading from the log file
        # So only read if there is no alternative
        if self._alat is None:
            self._alat = t
        if self.debug:
            print('_read_alat2: _alat={}'.format(self._alat) )
        return

    def _read_electrons(self, line):
        #
        # Read from log file
        #
        self.electrons = float(line.split()[4])
        if self.debug:
            print('_read_electrons: electrons={}'.format(self.electrons) )
        return

    def _read_energy(self, line):
        #
        # Read from log file
        #
        self.final_energy_without_entropy = float(line.split()[3]) * hartree2ev / 2.0
        self.final_free_energy = float(line.split()[3]) * hartree2ev / 2.0
        if self.debug:
            print('_read_energy: energy={}'.format(self.final_free_energy) )
        return

    def _read_energy_cutoff(self, line):
        #
        # Read from log file
        #
        self.energy_cutoff = float(line.split()[3]) * hartree2ev / 2.0
        if self.debug:
            print('_read_energy_cutoff: energy_cutoff={}'.format(self.energy_cutoff) )
        return

    def _read_kpoints(self, line):
        #
        # Read from log file
        #
        self.kpoints = int(line.split()[4])
        if self.debug:
            print('_read_kpoints: kpoints={}'.format(self.kpoints) )
        return

    def _read_kpoint_grid(self, line):
        #
        # Read from log file
        #
        line = self.file_descriptor.readline()
        self.kpoint_grid = [ float(f) for f in line.split()[0:3] ]
        if self.debug:
            print('_read_kpoints_grid kpoint_grid={}'.format(self.kpoint_grid) )
        return

    def _read_header(self, line):
        #
        # Read from dynG
        #
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.nspecies = int(line.split()[0])
        self.nions    = int(line.split()[1])
        t              = float(line.split()[3])
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        self._alat = t
        if self.debug:
            print('_read_header alat={}'.format(self._alat) )
        return

    def _read_epsilon(self, line):
        #
        # Read from dynG
        #
        self.file_descriptor.readline()
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric = []
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        if self.debug:
            print('_read_espilon zerof_optical_dielectric={}'.format(self.zerof_optical_dielectric) )
        return

    def _read_masses(self):
        #
        # Read from dynG
        #
        self.masses_per_type = []
        self.species      = []
        for i in range(self.nspecies):
            linea = self.file_descriptor.readline().replace('\'', '').split()
            self.species.append(linea[1].capitalize())
            # The factor of two is because au in pwscf are half mass of electron
            self.masses_per_type.append(float(linea[2])*2/amu)
        if self.debug:
            print('_read_masses masses={}'.format(self.masses_per_type) )
        line = ''
        self._read_dyng_coordinates(line)
        return

    def _read_dynamical(self, line):
        #
        # Read from dynG
        #
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
        if self.debug:
            print('_read_dynamical')
        return

    def _read_born_charges(self, line):
        #
        # Read from dynG
        #
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
        if self.debug:
            print('_read_born_charges')
        return

    def _read_cell_parameters(self, line):
        #
        # Read from the log file
        #
        line = line.replace( '(', ' ')
        line = line.replace( ')', ' ')
        linea = line.split()
        # Overwrite alat here as it is more accurate than anywhere else in the log file
        if len(linea) > 2:
            self._alat = float(linea[2])
        linea = self.file_descriptor.readline().split()
        avector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        bvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        cvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        self.unit_cells.append(UnitCell(avector, bvector, cvector))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].volume
        if self.debug:
            print('_read_cell_parameters: volume={}'.format(self.volume))
        return

    def _read_lattice_vectors(self, line):
        #
        # Read from dynG
        #
        linea = self.file_descriptor.readline().split()
        avector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        bvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        cvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        self.unit_cells.append(UnitCell(avector, bvector, cvector))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].volume
        if self.debug:
            print('_read_lattices_vectors: volume={}'.format(self.volume))
        self._read_masses()
        return

    def _read_fractional_coordinates(self,line):
        #
        # Read from log file
        #
        if self.nions <= 0:
            return
        species_list = []
        fractional_coordinates = []
        for i in range(self.nions):
            linea = self.file_descriptor.readline().split()
            species_list.append( linea[0] )
            fractional_coordinates.append([float(linea[1]), float(linea[2]), float(linea[3])])
        self.unit_cells[-1].set_fractional_coordinates(fractional_coordinates)
        self.unit_cells[-1].set_element_names(species_list)
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].volume
        if self.debug:
            print('_read_fractional_coordinates: volume={}'.format(self.volume))
        return

    def _read_dyng_coordinates(self,line):
        #
        # Read from dynG
        #
        self.masses = []
        self.atom_type_list = []
        self.ions_per_type = [ 0 for i in range(self.nspecies) ]
        species_list = []
        xyz_coordinates = []
        # It took a long time to work out that alat is in bohr
        const = self._alat/angs2bohr
        for i in range(self.nions):
            linea = self.file_descriptor.readline().split()
            species_index = int(linea[1])
            xyz_coordinates.append([const*float(linea[2]), const*float(linea[3]), const*float(linea[4])])
            self.masses.append(self.masses_per_type[species_index-1])
            self.atom_type_list.append(species_index-1)
            self.ions_per_type[species_index-1] += 1
            species_list.append(self.species[species_index-1])
        self.unit_cells[-1].set_xyz_coordinates(xyz_coordinates)
        self.unit_cells[-1].set_element_names(species_list)
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].volume
        if self.debug:
            print('_read_dyng_coordinates: volume={}'.format(self.volume))
        return
