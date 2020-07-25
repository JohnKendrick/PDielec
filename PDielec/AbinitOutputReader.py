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
"""Read the contents of a directory containing Abinit input and output files
   Inherit the following from the GenericOutputReader
   __init__
   print_info
   _read_output_file
"""

import re
import math
import numpy as np
from PDielec.Constants import amu, angs2bohr, atomic_number_to_element, hartree2ev
from PDielec.UnitCell import UnitCell
from PDielec.GenericOutputReader import GenericOutputReader


class AbinitOutputReader(GenericOutputReader):
    """Read the contents of a directory containing Abinit input and output files
       Inherit the following from the GenericOutputReader
       __init__
       print_info
       _read_output_file
"""

    def __init__(self, filenames):
        GenericOutputReader.__init__(self, filenames)
        self.type = 'Abinit output files'
        self._acell = None
        self._charges = None
        return

    def _read_output_files(self):
        """Read the Abinit file names"""
        # Define the search keys to be looked for in the files
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['dynamical']  = (re.compile('  Dynamical matrix,'), self._read_dynamical)
        self.manage['bornCharges']  = (re.compile('  Effective charges,'), self._read_born_charges)
        self.manage['epsilon']  = (re.compile('  Dielectric tensor,'), self._read_epsilon)
        self.manage['masses']   = (re.compile('              amu '), self._read_masses)
        self.manage['nions']    = (re.compile('            natom '), self._read_natom)
        self.manage['lattice']  = (re.compile('            rprim '), self._read_lattice_vectors)
        self.manage['xred']  = (re.compile('             xred '), self._read_xred)
        self.manage['typat']    = (re.compile('            typat '), self._read_typat)
        self.manage['ntypat']   = (re.compile('           ntypat '), self._read_ntypat)
        self.manage['acell']    = (re.compile('            acell '), self._read_acell)
        self.manage['nkpt']    = (re.compile('             nkpt '), self._read_kpoints)
        self.manage['band']    = (re.compile('            nband '), self._read_band)
        self.manage['band1']    = (re.compile('            nband1 '), self._read_band)
        self.manage['occupancy']    = (re.compile('              occ '), self._read_occupancy)
        self.manage['occupancy1']    = (re.compile('              occ1 '), self._read_occupancy)
        self.manage['ecut']    = (re.compile('^ *ecut '), self._read_energy_cutoff)
        self.manage['kptrlatt']    = (re.compile('         kptrlatt '), self._read_kpoint_grid)
        self.manage['electrons']    = (re.compile('  fully or partial'), self._read_electrons)
        self.manage['pressure']    = (re.compile('-Cartesian.*GPa'), self._read_pressure)
        self.manage['znucl']    = (re.compile('^  *znucl '), self._read_znucl)
        self.manage['totalenergy']    = (re.compile('^  *Total energy '), self._read_total_energy)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_total_energy(self, line):
        self.final_energy_without_entropy = float(line.split()[4]) * hartree2ev
        self.final_free_energy = float(line.split()[4]) * hartree2ev
        return

    def _read_znucl(self, line):
        self.species = []
        for z in line.split()[1:]:
            iz = int(float(z)+0.001)
            self.species.append(atomic_number_to_element[iz].capitalize())
        self.nspecies = len(self.species)
        return

    def _read_band(self, line):
        self.nbands = int(line.split()[1])
        return

    def _read_occupancy(self, line):
        occs = []
        occupancies = line.split()[1:]
        while len(occs) < self.nbands:
            occs+= [ float(f) for f in occupancies ]
            occupancies = self.file_descriptor.readline().split()
        sum = 0.0
        for f in occs:
            sum += f
        self.electrons = int(sum + 0.0001)
        return

    def _read_pressure(self, line):
        self.pressure = float(line.split()[7])
        return

    def _read_electrons(self, line):
        self.electrons = float(line.split()[6])
        return

    def _read_kpoint_grid(self, line):
        self.kpoint_grid = [ int(line.split()[1]), int(line.split()[5]), int(line.split()[9]) ]
        return

    def _read_kpoints(self, line):
        self.kpoints = int(line.split()[1])
        return

    def _read_energy_cutoff(self, line):
        self.energy_cutoff = hartree2ev * float(line.split()[1])
        return

    def _read_acell(self, line):
        self._acell = [float(f)/angs2bohr for f in line.split()[1:4]]
        return

    def _read_ntypat(self, line):
        self.nspecies = int(line.split()[1])
        return

    def _read_typat(self, line):
        # typat occurs last in the list of data items we need from the output file
        self.atom_type_list = [int(i)-1 for i in line.split()[1:]]
        self.masses = [None for i in range(self.nions)]
        self.ions_per_type = [0 for i in range(self.nspecies)]
        for i, a in enumerate(self.atom_type_list):
            self.ions_per_type[a-1] += 1
            self.masses[i] = self.masses_per_type[a]
        if self.species:
            species_list = [ self.species[i] for i in self.atom_type_list ]
            self.unit_cells[-1].set_element_names(species_list)
        return

    def _read_epsilon(self, line):
        for i in range(3):
            linea = self.file_descriptor.readline().split()
        nlines = 9
        for i in range(nlines):
            linea = self.file_descriptor.readline().split()
            if not linea:
                linea = self.file_descriptor.readline().split()
            j = int(linea[0])
            k = int(linea[2])
            self.zerof_optical_dielectric[j-1][k-1] = float(linea[4])
        return

    def _read_natom(self, line):
        self.nions = int(line.split()[1])
        # We can only create this once we know the number of ions
        self._charges = np.zeros((self.nions, 3, 3))
        return

    def _read_masses(self, line):
        self.masses_per_type = [float(f) for f in line.split()[1:]]
        return

    def _read_dynamical(self, line):
        # Read the dynamical matrix
        nmodes = self.nions*3
        hessian = np.zeros((nmodes, nmodes))
        for i in range(4):
            self.file_descriptor.readline()
        nlines = nmodes*nmodes
        for i in range(nlines):
            linea = self.file_descriptor.readline().split()
            if not linea:
                linea = self.file_descriptor.readline().split()
            diri  = int(linea[0])
            atomi = int(linea[1])
            dirj  = int(linea[2])
            atomj = int(linea[3])
            ipos  = (atomi - 1)*3 + diri - 1
            jpos  = (atomj - 1)*3 + dirj - 1
            # store the massweighted matrix
            hessian[ipos][jpos] = float(linea[4])/(amu*math.sqrt(self.masses[atomi-1]*self.masses[atomj-1]))
        # symmetrise, project diagonalise and store frequencies and normal modes
        self._dynamical_matrix(hessian)
        return

    def _read_born_charges(self, line):
        """Read the born charges from the outputfile file.
           Each row of the output refers to a given field direction
           Each column in the row refers the atomic displacement
           so the output is arranged [[a1x a1y a1z]
                                      [a2x a2y a2z]
                                      [a3x a3y a3z]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        for i in range(5):
            self.file_descriptor.readline()
        #  The charges are calculated in two ways, we take the mean of the phonon and the field
        nlines = 9*self.nions
        for i in range(nlines):
            linea = self.file_descriptor.readline().split()
            if not linea:
                linea = self.file_descriptor.readline().split()
            if int(linea[3]) > self.nions:
                ifield = int(linea[2])
                ixyz   = int(linea[0])
                iatom  = int(linea[1])
            else:
                ifield = int(linea[0])
                ixyz   = int(linea[2])
                iatom  = int(linea[3])
            self._charges[iatom-1][ifield-1][ixyz-1] += 0.5*float(linea[4])
        # Convert the charges
        self.born_charges = []
        for i in range(self.nions):
            atom = []
            for ifield in range(3):
                b = self._charges[i][ifield][:].tolist()
                atom.append(b)
            self.born_charges.append(atom)
        return

    def _read_xred(self, line):
        linea = line.split()[1:]
        fractional = []
        fractional.append( [ float(xyz) for xyz in linea ] )
        for i in range(self.nions-1):
            linea = self.file_descriptor.readline().split()
            fractional.append( [ float(xyz) for xyz in linea ] )
        # end for i
        self.unit_cells[-1].set_fractional_coordinates(fractional)
        if self.species:
            species_list = [ self.species[i] for i in self.atom_type_list ]
            self.unit_cells[-1].set_element_names(species_list)
    # end def


    def _read_lattice_vectors(self, line):
        linea = line.split()
        avector = [float(linea[1]), float(linea[2]), float(linea[3])]
        linea = self.file_descriptor.readline().split()
        bvector = [float(linea[0]), float(linea[1]), float(linea[2])]
        linea = self.file_descriptor.readline().split()
        cvector = [float(linea[0]), float(linea[1]), float(linea[2])]
        avector = [f * self._acell[0] for f in avector]
        bvector = [f * self._acell[1] for f in bvector]
        cvector = [f * self._acell[2] for f in cvector]
        self.unit_cells.append(UnitCell(avector, bvector, cvector))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].volume
        return
