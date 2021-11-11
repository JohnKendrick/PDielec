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
"""Read the contents of a directory containg VASP input and output files"""
import re
import numpy as np
from PDielec.UnitCell import UnitCell
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.Constants           import atomic_number_to_element


def myfloat(string):
    '''A replacement for float() which will return a large number if it finds a * in the string '''
    if '*' in string:
        return 9999.999
    else:
        return float(string)


class VaspOutputReader(GenericOutputReader):
    """Read the contents of a directory containg VASP input and output files"""

    def __init__(self, names):
        GenericOutputReader.__init__(self, names)
        self.type                    = 'Vasp output'
        self._pspots                  = {}
        self._ediff                   = 0.0
        self._pulay                  = None
        self._ibrion                  = 0
        self._potim                   = 0.0
        return

    def _read_output_files(self):
        """Read the vasp files in the directory"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['ionspertype']  = (re.compile('   ions per type ='), self._read_ionspertype)
        # self.manage['masses_skip']  = (re.compile('  Mass of Ions in am'), self._read_skip4)
        self.manage['pspots']       = (re.compile(' POTCAR:'), self._read_pspot)
        self.manage['arrays']       = (re.compile(' Dimension of arrays:'), self._read_array_dimensions)
        #self.manage['masses']       = (re.compile('   POMASS ='), self._read_masses)
        self.manage['spin']         = (re.compile('   ISPIN  = '), self._read_spin)
        self.manage['encut']        = (re.compile('   ENCUT  = '), self._read_encut)
        self.manage['ediff']        = (re.compile('   EDIFF  = '), self._read_ediff)
        self.manage['ibrion']       = (re.compile('   IBRION = '), self._read_ibrion)
        self.manage['potim']        = (re.compile('   POTIM  = '), self._read_potim)
        self.manage['nelect']       = (re.compile('   NELECT = '), self._read_nelect)
        self.manage['lattice']      = (re.compile('  volume of cell :'), self._read_lattice_vectors)
        self.manage['fractional']   = (re.compile(' position of ions in fractional coordinates'), self._read_fractional_coordinates)
        self.manage['forces']       = (re.compile(' POSITION  *TOTAL-FORCE'), self._read_forces)
        self.manage['energy']       = (re.compile('  FREE ENERGIE OF THE ION'), self._read_energy)
        self.manage['magnet']       = (re.compile(' number of electron '), self._read_magnet)
        self.manage['pressure']     = (re.compile('  external pressure ='), self._read_external_pressure)
        self.manage['skip1']        = (re.compile(' old parameters found'), self._read_skip4)
        self.manage['staticDielectric']  = (re.compile(' MACROSCOPIC STATIC DIELECTRIC TENSOR .including'), self._read_static_dielectric)
        self.manage['staticIonic']  = (re.compile(' MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONT'), self._read_ionic_dielectric)
        self.manage['bornCharges']  = (re.compile(' BORN EFFECTIVE CHARGES'), self._read_born_charges)
        self.manage['eigenvectors'] = (re.compile(' Eigenvectors and eige'), self._read_eigenvectors)
        self.manage['eigenskip']    = (re.compile(' Eigenvectors after division'), self._read_skip4)
        self.manage['elastic']      = (re.compile(' TOTAL ELASTIC MODULI'), self._read_elastic_constants)
        self.manage['kpoint']      = (re.compile('^Gamma'), self._read_kpoint_grid)
        self.manage['species']      = (re.compile('^ *Atomic configuration'), self._read_species)
        self.manage['newmasses']  = (re.compile('  Mass of Ions in am'), self._read_newmasses)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_forces(self, line):
        line = self.file_descriptor.readline()
        maxf = 0.0
        rmsf = 0.0
        for i in range(self.nions):
            line = self.file_descriptor.readline()
            forces = [ float(f) for f in line.split()[3:6] ]
            for f in forces:
                rmsf += f*f
                if abs(f) > maxf:
                    maxf = abs(f)
                # end if
            # end for f
        #end for i
        if not "max_force" in self.iterations:
            self.iterations["max_force"] = []
        self.iterations["max_force"].append(maxf)
        if not "rms_force" in self.iterations:
            self.iterations["rms_force"] = []
        self.iterations["rms_force"].append(rmsf)
        return

    def _read_species(self, line):
        line = self.file_descriptor.readline()
        nlines = int(line.split()[0])
        line = self.file_descriptor.readline()
        zcharge = 0.0
        for i in range(nlines):
            line = self.file_descriptor.readline()
            zcharge = zcharge + float(line.split()[4])
        self.species.append(atomic_number_to_element[int(zcharge+0.001)])
        self.nspecies = len(self.species)
        return

    def _read_kpoint_grid(self, line):
        line = self.file_descriptor.readline()
        self.kpoint_grid = [int(f) for f in line.split()[0:3] ]
        return

    def _read_ionspertype(self, line):
        self.ions_per_type = [int(i) for i in line.split()[4:]]
        self.nspecies = len(self.ions_per_type)
        return

    def _read_newmasses(self, line):
        self.masses_per_type = []
        line = self.file_descriptor.readline()
        mass_string = line[12:]
        start = 0
        increment = 6
        for i in range(self.nspecies):
            mass = mass_string[start:start+increment]
            self.masses_per_type.append(float(mass))
            start = start + increment
        self.masses = []
        self.atom_type_list = []
        self.species_list = []
        for k, mass in enumerate(self.masses_per_type):
            n = self.ions_per_type[k]
            for i in range(n):
                self.atom_type_list.append(k)
                self.masses.append(mass)
                self.species_list.append(self.species[k])
            # end loop over i
        # end look over current know types
        return

    def _read_masses(self, line):
        mass_string = line.split()[2]
        mass_string = mass_string.replace(";", "")
        self.masses_per_type.append(float(mass_string))
        return

    def _read_eigenvectors(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        n = 3 * self.nions
        for i in range(n):
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            imaginary = (line.split()[1] == "f/i=")
            if imaginary:
                # represent imaginary by negative real
                freq = -float(line.split()[6])
            else:
                freq = float(line.split()[7])
            # end if
            self.frequencies.append(freq)
            line = self.file_descriptor.readline()
            a = []
            for j in range(self.nions):
                line = self.file_descriptor.readline()
                a.append([float(line.split()[3]), float(line.split()[4]), float(line.split()[5])])
            # end for j
            self.mass_weighted_normal_modes.append(a)
        # end of for i in range(n)
        return

    def _read_born_charges(self, line):
        """Read the born charges from the OUTCAR file.
           Each row of the output refers to a given field direction
           Each column in the row refers the atomic displacement
           so the output is arranged [[a1x a1y a1z]
                                      [ a2x a2y a2z]
                                      [ a3x a3y a3z]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.born_charges = []
        for i in range(self.nions):
            line = self.file_descriptor.readline()
            b = []
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            self.born_charges.append(b)
        return

    def _read_elastic_constants(self, line):
        # Read the total elastic constants
        elastic_constants = []
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        econs = np.array(elastic_constants)
        # convert from kBar to GPa
        econs = econs / 10.0
        self.elastic_constants = econs.tolist()
        return

    def _read_ionic_dielectric(self, line):
        # Read the ionic contribution to the static dielectric and use it to computet
        # the full static dielectric constant
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        # the is zero frequency ionic contribution to the static dielectric
        ionic_dielectric = []
        if '*' in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        ionic_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if '*' in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        ionic_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if '*' in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        ionic_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        array1 = np.array(self.zerof_optical_dielectric)
        array2 = np.array(ionic_dielectric)
        array3 = array1 + array2
        self.zerof_static_dielectric = array3.tolist()
        return

    def _read_static_dielectric(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        # the is epsilon infinity
        self.zerof_optical_dielectric = []
        if '*' in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        self.zerof_optical_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if '*' in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        self.zerof_optical_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if '*' in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        self.zerof_optical_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        return

    def _read_skip4(self, line):
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        return

    def _read_external_pressure(self, line):
        # Vasp writes out kbar so convert to GPa
        self.pressure = float(line.split()[3])/10.0
        self.pressures.append(float(line.split()[3])/10.0)
        self._pulay = float(line.split()[8])/10.0
        return

    def _read_pspot(self, line):
        self._pspots[line.split()[2]] = line.split()[1]
        return

    def _read_array_dimensions(self, line):
        line = self.file_descriptor.readline()
        self.kpoints = int(line.split()[3])
        self.nbands = int(line.split()[14])
        line = self.file_descriptor.readline()
        self.nions = int(line.split()[11])
        return

    def _read_lattice_vectors(self, line):
        self.volume = float(line.split()[4])
        self.volumes.append(float(line.split()[4]))
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        avector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        bvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        cvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        cell = UnitCell(avector, bvector, cvector)
        if self.ncells > 0:
            cell.set_fractional_coordinates(self.unit_cells[-1].fractional_coordinates)
            cell.set_element_names(self.species_list)
        self.unit_cells.append(cell)
        self.ncells = len(self.unit_cells)
        return

    def _read_fractional_coordinates(self, line):
        n = 0
        ions = []
        for n in range(self.nions):
            line = self.file_descriptor.readline()
            ions.append([float(f) for f in line.split()[0:3]])
        self.unit_cells[-1].set_fractional_coordinates(ions)
        self.unit_cells[-1].set_element_names(self.species_list)
        return

    def _read_spin(self, line):
        self.spin = int(line.split()[2])
        return

    def _read_encut(self, line):
        self.energy_cutoff = float(line.split()[2])
        return

    def _read_ediff(self, line):
        self._ediff = float(line.split()[2])
        return

    def _read_ibrion(self, line):
        self._ibrion = int(line.split()[2])
        return

    def _read_potim(self, line):
        self._potim = float(line.split()[2])
        return

    def _read_nelect(self, line):
        self.electrons = int(float(line.split()[2]))
        return

    def _read_magnet(self, line):
        self.electrons = float(line.split()[3])
        if len(line.split()) > 5:
            self.magnetization = float(line.split()[5])
        else:
            self.magnetization = 0.0

    def _read_energy(self, line):
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.final_free_energy = float(line.split()[4])
        self.final_free_energies.append(float(line.split()[4]))
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.final_energy_without_entropy = float(line.split()[3])
        self.final_energies_without_entropy.append(float(line.split()[3]))
        return
