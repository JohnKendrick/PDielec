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
import numpy as np
from PDielec.UnitCell import UnitCell
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.Calculator          import initialise_diagonal_tensor
from PDielec.DielectricFunction  import DielectricFunction


class ExperimentOutputReader(GenericOutputReader):
    """Read the contents of a directory containg Experiment input and output files"""

    def __init__(self, names):
        GenericOutputReader.__init__(self, names)
        self.type                    = 'Experimental output'
        self._ion_type_index = {}
        self.CrystalPermittivity = None
        self.zerof_optical_dielectric = None
        self.oscillator_strengths = None
        self.frequencies = None
        return

    def _read_line(self):
        line = self.file_descriptor.readline()
        while not line.strip() or line[0] == '#':
            line = self.file_descriptor.readline()
        return line

    def _read_output_files(self):
        """Read the Experiment files in the directory"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['lattice']       = (re.compile('lattice'),          self._read_lattice_vectors)
        self.manage['CPK_LATTICE']   = (re.compile('&CELL'),            self._read_cpk_lattice_vectors)
        self.manage['cpk_lattice']   = (re.compile('&cell'),            self._read_cpk_lattice_vectors)
        self.manage['species']       = (re.compile('species'),          self._read_species)
        self.manage['fractional']    = (re.compile('unitcell'),         self._read_fractional_coordinates)
        self.manage['CPK_CARTESIANS']= (re.compile('&COORD'),           self._read_cpk_coords)
        self.manage['cpk_cartesians']= (re.compile('&coord'),           self._read_cpk_coords)
        self.manage['static']        = (re.compile('static'),           self._read_static_dielectric)
        self.manage['epsinf']        = (re.compile('epsinf'),           self._read_static_dielectric)
        self.manage['fpsq']          = (re.compile('fpsq'),             self._read_fpsq_model)
        self.manage['drude-lorentz'] = (re.compile('drude-lorentz'),    self._read_drude_lorentz_model)
        self.manage['constant']      = (re.compile('constant'),         self._read_constant_model)
        # The order of these interpolates is important !
        self.manage['interpolate3']  = (re.compile('interpolate_3'),    self._read_interpolate3_model)
        self.manage['interpolate6']  = (re.compile('interpolate_6'),    self._read_interpolate6_model)
        self.manage['interpolate']   = (re.compile('interpolate'),      self._read_interpolate1_model)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_constant_model(self, line):
        """
        Read in a full constant dielectric tensor 3 values on each line
        """
        od = []
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        # If we have complex input return a complex list, otherwise return a real list
        odc = np.array(od,dtype=complex)
        odi = np.absolute(np.imag(odc))
        sumi = np.sum(odi)
        if sumi < 1.0e-12:
            odc = np.real(odc)
        parameters = odc.tolist()
        self.CrystalPermittivity = DielectricFunction('constant',parameters=parameters)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_interpolate1_model(self, line):
        """
        Read in a tabulated permittivity and use it for interpolation -
        """
        full_eps = []
        contributions = []
        line = self._read_line()
        line = line.lower()
        line = line.replace(',',' ')
        split_line = line.split()
        while split_line[0] != '&end' and split_line[0] != 'end':
            omega  = float(split_line[0])
            epsrxx = float(split_line[1])
            epsixx = float(split_line[2])
            full_eps.append( (omega, epsrxx, epsixx, epsrxx, epsixx, epsrxx, epsixx) )
            line = self._read_line()
            line = line.lower()
            line = line.replace(',',' ')
            split_line = line.split()
        # end for i
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction('interpolate_3', parameters=full_eps)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_interpolate3_model(self, line):
        """
        Read in a tabulated permittivity and use it for interpolation -
        """
        full_eps = []
        contributions = []
        line = self._read_line()
        line = line.lower()
        line = line.replace(',',' ')
        split_line = line.split()
        while split_line[0] != '&end' and split_line[0] != 'end':
            omega  = float(split_line[0])
            epsrxx = float(split_line[1])
            epsixx = float(split_line[2])
            epsryy = float(split_line[3])
            epsiyy = float(split_line[4])
            epsrzz = float(split_line[5])
            epsizz = float(split_line[6])
            full_eps.append( (omega, epsrxx, epsixx, epsryy, epsiyy, epsrzz, epsizz) )
            line = self._read_line()
            line = line.lower()
            line = line.replace(',',' ')
            split_line = line.split()
        # end for i
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction('interpolate_3', parameters=full_eps)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_interpolate6_model(self, line):
        """
        Read in a tabulated permittivity and use it for interpolation -
        """
        full_eps = []
        contributions = []
        line = self._read_line()
        line = line.lower()
        line = line.replace(',',' ')
        split_line = line.split()
        while split_line[0] != '&end' and split_line[0] != 'end':
            omega  = float(split_line[0])
            epsrxx = float(split_line[1])
            epsixx = float(split_line[2])
            epsryy = float(split_line[3])
            epsiyy = float(split_line[4])
            epsrzz = float(split_line[5])
            epsizz = float(split_line[6])
            epsrxy = float(split_line[7])
            epsixy = float(split_line[8])
            epsrxz = float(split_line[9])
            epsixz = float(split_line[10])
            epsryz = float(split_line[11])
            epsiyz = float(split_line[12])
            full_eps.append( (omega, epsrxx, epsixx, epsryy, epsiyy, epsrzz, epsizz, epsrxy, epsixy, epsrxz, epsixz, epsryz, epsiyz) )
            line = self._read_line()
            line = line.lower()
            line = line.replace(',',' ')
            split_line = line.split()
        # end for i
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction('interpolate_6', parameters=full_eps)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_drude_lorentz_model(self, line):
        """
        Read in the drude_lorentz model parameters
        A simple example 2 oscillator model is given below
        All units are in cm-1
        drude-lorentz
        xx  1
        413.7 1050.0  22.2
        yy  1
        413.7 1050.0  22.2
        zz  1
        413.7 1050.0  22.2
        """
        diag_eps = []
        for diag in range(0,3):
            line = self._read_line().split()
            element = line[0]
            n = int(line[1])
            contributions = []
            for i in range(n):
                line = self._read_line().split()
                omega    = float(line[0])
                strength = float(line[1])
                gamma    = float(line[2])
                contributions.append( (omega, strength, gamma) )
            # end for i
            diag_eps.append(contributions)
        # end for diag
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction('drude-lorentz', parameters=diag_eps)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_fpsq_model(self, line):
        """
        Read in the fpsq model parameters
        There is a seperator between each diagonal contribution to the permittivity
        For alpha quartz the input looks like this;
        epsxx 4 # 1.5K Data # A2 Phonon modes
          360.7     1.5        384.8    1.5
          497.9     3.1        553.6    2.8
          773.7     5.4        789.9    6.3
         1073.0     6.2       1238.7   12.4
        epsyy 4
          360.7     1.5        384.8    1.5
          497.9     3.1        553.6    2.8
          773.7     5.4        789.9    6.3
         1073.0     6.2       1238.7   12.4
        epszz 6 # E Phonon modes
          391.5     1.5        403.0    1.5
          454.0     2.6        510.5    1.5
          695.9     4.9        698.4    4.0
          797.2     4.8        810.0    4.3
         1063.7     6.1       1230.7    8.2
         1157.2     6.2       1154.9    6.1
         """
        diag_eps = []
        for diag in range(0,3):
            line = self._read_line().split()
            element = line[0]
            n = int(line[1])
            contributions = []
            for i in range(n):
                line = self._read_line().split()
                omega_to = float(line[0])
                gamma_to = float(line[1])
                omega_lo = float(line[2])
                gamma_lo = float(line[3])
                contributions.append( (omega_to, gamma_to, omega_lo, gamma_lo) )
            # end for i
            diag_eps.append(contributions)
        # end for diag
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction('fpsq', parameters=diag_eps)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_frequencies(self, line):
        nfreq = int(line.split()[1])
        self.frequencies = []
        self.oscillator_strengths = []
        for i in range(nfreq):
            line = self._read_line()
            self.frequencies.append(float(line.split()[0]))
            strength = float(line.split()[1])
            self.oscillator_strengths.append(initialise_diagonal_tensor( [strength, strength, strength] ) )
        return

    def _read_species(self, line):
        nspecies = int(line.split()[1])
        self.species = []
        self.masses_per_type = []
        for i in range(nspecies):
            line = self._read_line()
            species = line.split()[0]
            self.species.append(species)
            self.masses_per_type.append(float(line.split()[1]))
            self._ion_type_index[species] = self.nspecies
            self.nspecies = len(self.species)
        return

    def _read_cpk_lattice_vectors(self, line):
        """Process a cpk input file with cell information"""
        line = line.lower()
        split_line = line.split()
        alpha = 90.0
        beta = 90.0
        gamma = 90.0
        while split_line[0] != '&end' and split_line[0] != 'end':
            if  split_line[0] == 'abc':
                a = float(split_line[1])
                b = float(split_line[2])
                c = float(split_line[3])
            elif split_line[0] == 'alpha_beta_gamma':
                alpha = float(split_line[1])
                beta = float(split_line[2])
                gamma = float(split_line[3])
            # end if split-line
            line = self._read_line()
            line = line.lower()
            split_line = line.split()
        # end while
        cell = UnitCell(a,b,c,alpha,beta,gamma)
        self.unit_cells.append(cell)
        self.ncells = len(self.unit_cells)
        self.volume = cell.volume
        if self.CrystalPermittivity:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_lattice_vectors(self, line):
        line = self._read_line()
        scalar = float(line.split()[0])
        line = self._read_line()
        avector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        line = self._read_line()
        bvector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        line = self._read_line()
        cvector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        cell = UnitCell(avector, bvector, cvector)
        self.unit_cells.append(cell)
        self.ncells = len(self.unit_cells)
        self.volume = cell.volume
        if self.CrystalPermittivity:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_cpk_coords(self, line):
        line = self._read_line()
        line = line.lower()
        line = line.replace(',',' ')
        split_line = line.split()
        species_list = []
        ions = []
        self.nspecies = 0
        while split_line[0] != '&end' and split_line[0] != 'end':
            self.nions += 1
            species = line.split()[0]
            if species not in self.species:
                self.species.append(species)
                self.masses_per_type.append(1.0)
                self._ion_type_index[species] = self.nspecies
                self.ions_per_type.append(1)
                self.nspecies = len(self.species)
            index = self._ion_type_index[species]
            self.atom_type_list.append(index)
            species_list.append(species)
            ions.append([float(f) for f in line.split()[1:4]])
            self.ions_per_type[index] += 1
            self.masses.append(self.masses_per_type[index])
            line = self._read_line()
            line = line.lower()
            line = line.replace(',',' ')
            split_line = line.split()
        # end while
        self.unit_cells[-1].set_xyz_coordinates(ions)
        self.unit_cells[-1].set_element_names(species_list)
        if self.oscillator_strengths == None:
            self.oscillator_strengths = np.zeros( (3*self.nions,3,3) )
        if self.frequencies == None:
            self.frequencies = np.zeros( (3*self.nions) )

    def _read_fractional_coordinates(self, line):
        n = 0
        ions = []
        self.nions = int(line.split()[1])
        self.ions_per_type = [ 0 for s in self.species ]
        self.masses = []
        species_list = []
        for n in range(self.nions):
            line = self._read_line()
            species = line.split()[0]
            index = self._ion_type_index[species]
            self.atom_type_list.append(index)
            species_list.append(species)
            ions.append([float(f) for f in line.split()[1:4]])
            self.ions_per_type[index] += 1
            self.masses.append(self.masses_per_type[index])
        self.unit_cells[-1].set_fractional_coordinates(ions)
        self.unit_cells[-1].set_element_names(species_list)
        if self.oscillator_strengths == None:
            self.oscillator_strengths = np.zeros( (3*self.nions,3,3) )
        if self.frequencies == None:
            self.frequencies = np.zeros( (3*self.nions) )
        return

    def _read_static_dielectric(self, line):
        # the is epsilon infinity
        od = []
        line = self._read_line()
        od.append([np.complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([np.complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([np.complex(f) for f in line.split()[0:3]])
        # If we have complex input return a complex list, otherwise return a real list
        odc = np.array(od,dtype=complex)
        odi = np.absolute(np.imag(odc))
        sumi = np.sum(odi)
        if sumi < 1.0e-12:
            odc = np.real(odc)
        self.zerof_optical_dielectric = odc.tolist()
        if self.CrystalPermittivity:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        return

    def calculate_mass_weighted_normal_modes(self):
        self.mass_weighted_normal_modes = np.zeros( (3*self.nions,self.nions,3) )
        return self.mass_weighted_normal_modes
