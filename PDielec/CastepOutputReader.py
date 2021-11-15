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
"""Read the contents of a Castep output file
   Inherit the following from the GenericOutputReader
   __init__
   print_info
   _read_output_file
   _DyanmicalMatrix
"""
import re
import os
import numpy as np
from PDielec.UnitCell import UnitCell
from PDielec.GenericOutputReader import GenericOutputReader


class CastepOutputReader(GenericOutputReader):
    """Read the contents of a Castep output file
       Inherit the following from the GenericOutputReader
       __init__
       print_info
       _read_output_file
       _DyanmicalMatrix
       _read_till_phrase
"""

    def __init__(self, filenames):
        GenericOutputReader.__init__(self, filenames)
        if filenames[0].find(".castep"):
            seedname, ext = os.path.splitext(filenames[0])
        elif filenames[0].find(".phonon"):
            seedname, ext = os.path.splitext(filenames[0])
        self._castepfile             = seedname+".castep"
        self._phononfile             = seedname+".phonon"
        self.names                   = [self._castepfile, self._phononfile]
        self._outputfiles            = [self._castepfile, self._phononfile]
        self.type                    = 'Castep output'
        # Specific Castep Reader Variables
        self._pspots                  = {}
        self._ediff                   = 0.0
        self._epsilon                 = None
        self._nbranches              = 0
        self._pulay                  = None
        self._ion_type_index         = {}
        self._ion_index_type         = {}
        self._intensities             = None
        self.de_ion                  = []
        self.fmax                    = []
        self.dr_max                  = []
        self.smax                    = []
        return

    def _read_output_files(self):
        """For the .castep file"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['spin']          = (re.compile(' *net spin   of'), self._read_spin)
        self.manage['nelect']        = (re.compile(' *number of  electrons'), self._read_nelect)
        self.manage['cellcontents']  = (re.compile(' *Unit Cell'), self._read_cellcontents)
        self.manage['pspots']        = (re.compile(' *Files used for pseudopotentials:'), self._read_pspot)
        self.manage['masses']        = (re.compile(' *Mass of species in AMU'), self._read_masses)
        self.manage['kpoints']       = (re.compile(' *Number of kpoints used'), self._read_kpoints)
        self.manage['kpoint_grid']   = (re.compile(' *MP grid size for SCF'), self._read_kpoint_grid)
        self.manage['energies']   = (re.compile(' *Final energy, E'), self._read_energies)
        self.manage['energies2']   = (re.compile('Final energy ='), self._read_energies2)
        self.manage['energies3']        = (re.compile('Dispersion corrected final energy'), self._read_energies3)
        self.manage['energy_cutoff'] = (re.compile(' *plane wave basis set cut'), self._read_energy_cutoff)
        self.manage['convergence'] = (re.compile('.*finished iteration'), self._read_convergence)
        self.manage['nbands']        = (re.compile(' *number of bands'), self._read_nbands)
        self.manage['pressure']      = (re.compile(' *\* *Pressure: '), self._read_external_pressure)
        self.manage['opticalDielectric']  = (re.compile(' *Optical Permittivity'), self._read_dielectric)
        self.manage['bornCharges']    = (re.compile(' *Born Effective Charges'), self._read_born_charges)
        #  For the .phonon file
        self.manage['frequency']      = (re.compile('     q-pt=    1    0.000000  0.000000  0.000000      1.0000000000 *$'), self._read_frequencies)
        self.manage['nbranches']      = (re.compile(' Number of branches'), self._read_nbranches)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_nbranches(self, line):
        # phonon file being read
        self._nbranches = int(line.split()[3])

    def _read_frequencies(self, line):
        # phonon file being read
        self.mass_weighted_normal_modes = []
        # maybe should use _nbranches
        frequencies = []
        intensities = []
        normal_modes = []
        for imode in range(self._nbranches):
            line = self.file_descriptor.readline()
            frequencies.append(float(line.split()[1]))
            intensities.append(float(line.split()[2]))

        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        for imode in range(self._nbranches):
            a = []
            for ion in range(self.nions):
                line = self.file_descriptor.readline()
                a.append([float(line.split()[2]), float(line.split()[4]), float(line.split()[6])])
            # end for ion
            normal_modes.append(a)
        # end for imode
        self.frequencies = []
        self._intensities = []
        self.mass_weighted_normal_modes = []
        # now reads all frequencies imaginary or not
        # imaginary frequencies are indicated by real negative values
        for i, freq in enumerate(frequencies):
            self.frequencies.append(freq)
            self._intensities.append(intensities[i])
            self.mass_weighted_normal_modes.append(normal_modes[i])
            # end of if freq
        # end of for freq
        return

    def _read_kpoint_grid(self, line):
        self.kpoint_grid = [ int(line.split()[7]), int(line.split()[8]), int(line.split()[9]) ]
        return

    def _read_kpoints(self, line):
        self.kpoints = int(line.split()[5])
        return

    def _read_nbands(self, line):
        self.nbands = int(line.split()[4])
        return

    def _read_pseudoatom(self, line):
        species = line.split()[5].capitalize()
        # These are two private dictionary to map the species name to a type index
        self.species.append(species)
        self._ion_type_index[species] = self.nspecies
        self._ion_index_type[self.nspecies] = species
        self.nspecies += 1
        return

    def _read_cellcontents(self, line):
        line = self._read_till_phrase(re.compile(' *Real Lattice'))
        line = self.file_descriptor.readline()
        avector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        bvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        cvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        self.unit_cells.append(UnitCell(avector, bvector, cvector))
        self.ncells = len(self.unit_cells)
        line = self._read_till_phrase(re.compile(' *Lattice*'))
        line = self._read_till_phrase(re.compile(' *Current cell volume'))
        self.volume = float(line.split()[4])
        self.volumes.append(float(line.split()[4]))
        if self.nions == 0:
            line = self._read_till_phrase(re.compile(' *Total number of ions in cell'))
            if len(self.unit_cells) == 1:
                self.nions = int(line.split()[7])
                line = self.file_descriptor.readline()
                self.nspecies = int(line.split()[7])
                line = self.file_descriptor.readline()
                line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
        # end if self.nions
        else:
            # Not all Unit Cell output is the same in Castep 17
            line = self._read_till_phrase(re.compile(' *xxxxxxxxxxxxxxxxxxxxxxxxxxxx'))
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
        # end of else
        fractional_coordinates = []
        # ions_per_type is a dictionary here and is local
        ions_per_type = {}
        self.atom_type_list = []
        self.species = []
        species_list = []
        for i in range(self.nions):
            line = self.file_descriptor.readline()
            atom_type = line.split()[1].capitalize()
            species_list.append(atom_type)
            if atom_type not in ions_per_type:
                self.species.append(atom_type)
                ions_per_type[atom_type] = 0
            self.atom_type_list.append(self.species.index(atom_type))
            ions_per_type[atom_type] += 1
            atom_frac = [float(f) for f in line.split()[3:6]]
            fractional_coordinates.append(atom_frac)
        # At some point we should store the fractional coordinates
        for species in self.species:
            n = ions_per_type[species]
            # self.ions_per_type is a list
            self.ions_per_type.append(n)
        self.unit_cells[-1].set_fractional_coordinates(fractional_coordinates)
        self.unit_cells[-1].set_element_names(species_list)
        return

    def _read_masses(self, line):
        self.masses = []
        self.masses_per_type = []
        self.species = []
        self.nspecies = 0
        line = self.file_descriptor.readline()
        while len(line.split()) != 0:
            species = line.split()[0].capitalize()
            mass = float(line.split()[1])
            nions = self.ions_per_type[self.nspecies]
            self.species.append(species)
            # These are two private dictionary to map the species name to a type index
            # self._ion_type_index[species] = self.nspecies
            self._ion_index_type[self.nspecies] = species
            self.masses_per_type.append(mass)
            self.nspecies += 1
            for j in range(nions):
                self.masses.append(mass)
            line = self.file_descriptor.readline()
        # end while loop
        return

    def _read_born_charges(self, line):
        """Read the born charges from the castep file.
           Each column of the output refers to a given field direction
           Each row in the row refers the atomic displacement
           So the numbers in the output are arrange  [[a1x a2x a3x]
                                                      [a1y a2y a3y]
                                                      [a1z a2z a3z]]
           The output tensor needs them arranged     [[a1x a1y a1z]
                                                      [a2x a2y a2z]
                                                      [a3x a3y a3z]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        line = self.file_descriptor.readline()
        self.born_charges = []
        for i in range(self.nions):
            b = []
            line = self.file_descriptor.readline()
            b.append([float(f) for f in line.split()[2:5]])
            line = self.file_descriptor.readline()
            b.append([float(f) for f in line.split()[0:3]])
            line = self.file_descriptor.readline()
            b.append([float(f) for f in line.split()[0:3]])
            B = np.array(b)
            C = B.T
            self.born_charges.append(C.tolist())
        return

    def _read_dielectric(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        # this is epsilon infinity
        self.zerof_optical_dielectric = []
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric = []
        self.zerof_optical_dielectric.append([float(f) for f in line.split()[0:3]])
        self.zerof_static_dielectric.append([float(f) for f in line.split()[3:6]])
        line = self.file_descriptor.readline()
        self.zerof_optical_dielectric.append([float(f) for f in line.split()[0:3]])
        self.zerof_static_dielectric.append([float(f) for f in line.split()[3:6]])
        line = self.file_descriptor.readline()
        self.zerof_optical_dielectric.append([float(f) for f in line.split()[0:3]])
        self.zerof_static_dielectric.append([float(f) for f in line.split()[3:6]])
        return

    def _read_skip4(self, line):
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        return

    def _read_external_pressure(self, line):
        self.pressure = float(line.split()[2])
        return

    def _read_pspot(self, line):
        for i in range(self.nspecies):
            line = self.file_descriptor.readline()
            self._pspots[line.split()[0]] = line.split()[1]
        return

    def _read_spin(self, line):
        self.spin = int(float(line.split()[5]))
        return

    def _read_energy_cutoff(self, line):
        self.energy_cutoff = float(line.split()[6])
        return

    def _read_ediff(self, line):
        self._ediff = float(line.split()[2])
        return

    def _read_nelect(self, line):
        self.electrons = int(float(line.split()[4]))
        return

    def _read_epsilon(self, line):
        var = line.split()[2].lower()
        self.epsilon = False
        if var == '.true.' or var == 't' or var == 'true':
            self.epsilon = True
        return

    def _read_energies(self, line):
        self.energiesDFT.append(float(line.split()[4]))
        self.final_energy_without_entropy = float(line.split()[4])
        line = self.file_descriptor.readline()
        self.final_free_energy = float(line.split()[5])
        self.geomsteps = int(len(self.energiesDFT))
        return

    def _read_energies2(self, line):
        self.energiesDFT.append(float(line.split()[3]))
        self.final_energy_without_entropy = float(line.split()[3])
        self.final_free_energy = float(line.split()[3])
        self.geomsteps = int(len(self.energiesDFT))
        return

        
    def _read_energies3(self, line):
        self.energiesDFT_disp.append(float(line.split()[5]))
        self.final_energy_without_entropy = float(line.split()[5])
        self.final_free_energy = float(line.split()[5])
        self.geomsteps = int(len(self.energiesDFT_disp))
        return
        
    def _read_convergence(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.de_ion.append(float(line.split()[3]))
        line = self.file_descriptor.readline()
        self.fmax.append(float(line.split()[3]))
        line = self.file_descriptor.readline()
        self.dr_max.append(float(line.split()[3]))
        line = self.file_descriptor.readline()
        self.smax.append(float(line.split()[3]))
        line = self.file_descriptor.readline()
        return    
        
            

