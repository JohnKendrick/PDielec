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
"""Read the VASP OUTCAR and other files"""
import string
import re
import numpy as np
import math 
import os, sys
from Python.Constants import *
from Python.UnitCell import *
from Python.GenericOutputReader import *
    
class VaspOutputReader(GenericOutputReader):
    """Read the contents of a directory containg VASP input and output files"""

    def __init__(self,names):
        GenericOutputReader.__init__(self,names)
        if names[0] == "OUTCAR":
            directory = "."
        else:
            directory = names[0]
        self._outputfiles            = [ os.path.join(directory,'OUTCAR') ]
        self.name                    = os.path.abspath(directory)
        self.type                    = 'Vasp output'
        self.pspots                  = {}
        self.spin                    = 0
        self.encut                   = 0.0
        self.ediff                   = 0.0
        self.nelect                  = 0
        self.nkpts                   = 0
        self.nbands                  = 0
        self.electrons               = 0.0
        self.magnetization           = 0.0
        self.final_energy_without_entropy = 0
        self.final_free_energy       = 0
        self.pressure                = None
        self._pulay                   = None
        self.elastic_constants       = [ [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0] ]
        self.masses_list         = []
        self.ions_per_type             = []
        self.ibrion                 = 0
        self.potim                  = 0.0
        self.epsilon                = None
        return

    def _ReadOutputFiles(self):
        """Read the vasp files in the directory"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['ionspertype']  = (re.compile('   ions per type ='),self._read_ionspertype)
        self.manage['masses_skip']   = (re.compile('  Mass of Ions in am'),self._read_skip4)
        self.manage['pspots']   = (re.compile(' POTCAR:'),self._read_pspot)
        self.manage['arrays']   = (re.compile(' Dimension of arrays:'),self._read_array_dimensions)
        self.manage['masses']   = (re.compile('   POMASS ='),self._read_masses)
        self.manage['spin']     = (re.compile('   ISPIN  = '),self._read_spin)
        self.manage['encut']    = (re.compile('   ENCUT  = '),self._read_encut)
        self.manage['ediff']    = (re.compile('   EDIFF  = '),self._read_ediff)
        self.manage['ibrion']   = (re.compile('   IBRION = '),self._read_ibrion)
        self.manage['potim']    = (re.compile('   POTIM  = '),self._read_potim)
        self.manage['nelect']   = (re.compile('   NELECT = '),self._read_nelect)
        self.manage['epsilon']  = (re.compile('   LEPSILON= '),self._read_epsilon)
        self.manage['lattice']  = (re.compile('  volume of cell :'),self._read_lattice_vectors)
        self.manage['fractional']  = (re.compile(' position of ions in fractional coordinates'),self._read_fractional_coordinates)
        self.manage['energy']  = (re.compile('  FREE ENERGIE OF THE ION'),self._read_energy)
        self.manage['magnet']  = (re.compile(' number of electron '),self._read_magnet)
        self.manage['pressure']  = (re.compile('  external pressure ='),self._read_external_pressure)
        self.manage['skip1']  = (re.compile(' old parameters found'),self._read_skip4)
        self.manage['staticDielectric']  = (re.compile(' MACROSCOPIC STATIC DIELECTRIC TENSOR .including'),self._read_static_dielectric)
        self.manage['staticIonic']  = (re.compile(' MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONT'),self._read_ionic_dielectric)
        self.manage['bornCharges']  = (re.compile(' BORN EFFECTIVE CHARGES'),self._read_born_charges)
        self.manage['eigenvectors']  = (re.compile(' Eigenvectors and eige'),self._read_eigenvectors)
        self.manage['eigenskip']  = (re.compile(' Eigenvectors after division'),self._read_skip4)
        self.manage['elastic']  = (re.compile(' TOTAL ELASTIC MODULI'),self._read_elastic_constants)
        for f in self._outputfiles:
            self._ReadOutputFile(f)
        return


    def _read_ionspertype(self,line):
        self.ions_per_type=[ int(i) for i in line.split()[4:] ]
        self.nspecies=len(self.ions_per_type)
        self.masses = []
        for k,mass in enumerate(self.masses_list):
            n = self.ions_per_type[k]
            for i in range(0,n):
                self.masses.append(mass)
            # end loop over i
        # end look over current know types
        return

    def _read_masses(self,line):
        mass_string = line.split()[2]
        mass_string = mass_string.replace(";","")
        self.masses_list.append ( float(mass_string) ) 
        return

    def _read_eigenvectors(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        n = 3 * self.nions
        for i in range(n) :
          line = self.fd.readline()
          line = self.fd.readline()
          imaginary = ( line.split()[1] == "f/i=" )
          if imaginary :
            # represent imaginary by negative real 
            freq = -float(line.split()[6])
          else :
            freq = float(line.split()[7])
          self.frequencies.append(freq)
          line = self.fd.readline()
          a = []
          for j in range(self.nions) :
            line = self.fd.readline()
            a.append([ float(line.split()[3]), float(line.split()[4]), float(line.split()[5])] )
            # end if
          # end for j 
          self.mass_weighted_normal_modes.append(a)
        # end of for i in range(n)
        return

    def _read_born_charges(self,line):
        """Read the born charges from the OUTCAR file.  
           Each row of the output refers to a given field direction
           Each column in the row refers the atomic displacement 
           so the output is arranged [ [ a1x a1y a1z ] 
                                       [ a2x a2y a2z ] 
                                       [ a3x a3y a3z ]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        line = self.fd.readline()
        line = self.fd.readline()
        self.born_charges = []
        for i in range(self.nions) :
          line = self.fd.readline()
          b = []
          b.append( [ float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ] )
          line = self.fd.readline()
          self.born_charges.append(b)
        if self.neutral:
            self._BornChargeSumRule()
        return

    def _read_elastic_constants(self,line):
        # Read the total elastic constants
        elastic_constants = []
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        elastic_constants.append( [ float(f) for f in line.split()[1:7] ] )
        line = self.fd.readline()
        elastic_constants.append( [ float(f) for f in line.split()[1:7] ] )
        line = self.fd.readline()
        elastic_constants.append( [ float(f) for f in line.split()[1:7] ] )
        line = self.fd.readline()
        elastic_constants.append( [ float(f) for f in line.split()[1:7] ] )
        line = self.fd.readline()
        elastic_constants.append( [ float(f) for f in line.split()[1:7] ] )
        line = self.fd.readline()
        elastic_constants.append( [ float(f) for f in line.split()[1:7] ] )
        a1 = np.array(elastic_constants)
        # convert from kBar to GPa
        a1 = a1 / 10.0
        self.elastic_constants = a1.tolist()
        return

    def _read_ionic_dielectric(self,line):
        # Read the ionic contribution to the static dielectric and use it to computet
        # the full static dielectric constant
        line = self.fd.readline()
        line = self.fd.readline()
        # the is zero frequency ionic contribution to the static dielectric
        ionic_dielectric = []
        ionic_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        line = self.fd.readline()
        ionic_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        line = self.fd.readline()
        ionic_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        a1 = np.array(self.zerof_optical_dielectric)
        a2 = np.array(ionic_dielectric)
        a3 = a1 + a2
        self.zerof_static_dielectric = a3.tolist()
        return

    def _read_static_dielectric(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        # the is epsilon infinity
        self.zerof_optical_dielectric = []
        self.zerof_optical_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        line = self.fd.readline()
        self.zerof_optical_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        line = self.fd.readline()
        self.zerof_optical_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        return

    def _read_skip4(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        return


    def _read_external_pressure(self,line):
        self.pressure = float(line.split()[3])
        self._pulay = float(line.split()[8])
        return

    def _read_pspot(self,line):
        self.pspots[line.split()[2]] = line.split()[1]
        return

    def _read_array_dimensions(self,line):
        line = self.fd.readline()
        self.nkpts = int(line.split()[3])
        self.nbands = int(line.split()[14])
        line = self.fd.readline()
        self.nions = int(line.split()[11])
        return

    def _read_lattice_vectors(self,line):
        self.volume = float(line.split()[4])
        line = self.fd.readline()
        line = self.fd.readline()
        aVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        line = self.fd.readline()
        bVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        line = self.fd.readline()
        cVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        cell = UnitCell(aVector, bVector, cVector)
        if self.ncells > 0:
            cell.fractionalCoordinates(self.unitCells[-1].fractional_coordinates)
        self.unitCells.append(cell)
        self.ncells = len(self.unitCells)
        return

    def _read_fractional_coordinates(self,line):
        n = 0
        ions = []
        for n in range(self.nions) :
            line = self.fd.readline()
            ions.append( [ float(f) for f in line.split()[0:3] ] )
        self.unitCells[-1].fractionalCoordinates(ions)
        return

    def _read_spin(self,line):
        self.spin = int(line.split()[2]) 
        return
    def _read_encut(self,line):
        self.encut = float(line.split()[2]) 
        return
    def _read_ediff(self,line):
        self.ediff = float(line.split()[2]) 
        return
    def _read_ibrion(self,line):
        self.ibrion = int(line.split()[2]) 
        return
    def _read_potim(self,line):
        self.potim = float(line.split()[2]) 
        return
    def _read_nelect(self,line):
        self.nelect = int(float(line.split()[2]) )
        return
    def _read_epsilon(self,line):
        var = line.split()[2].lower()
        self.epsilon = False
        if var == '.true.' or var == 't' or var == 'true' :
          self.epsilon = True
        return

    def _read_magnet(self,line):
        self.electrons     = float(line.split()[3])
        if len( line.split() ) > 5:
            self.magnetization = float(line.split()[5])
        else:
            self.magnetization = 0.0

    def _read_energy(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        self.final_free_energy = float(line.split()[4])
        line = self.fd.readline()
        line = self.fd.readline()
        self.final_energy_without_entropy = float(line.split()[3])
        return

