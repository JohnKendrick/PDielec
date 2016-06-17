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
"""Read the AbInit output files"""
import string
import re
import numpy as np
import math 
import os, sys
from Python.Constants import *
from Python.UnitCell import *
from Python.GenericOutputReader import *
    
class QEOutputReader(GenericOutputReader):
    """Read the contents of a dmat file containing QE dynmical matrix"""

    def __init__(self,filenames):
        GenericOutputReader.__init__(self,filenames)
        self.type                    = 'QE output'
        return

    def _ReadOutputFiles(self):
        """Read the AbInit file names"""
        # Define the search keys to be looked for in the files
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['header']   = (re.compile('Dynamical matrix file'),self._read_header)
        self.manage['lattice']  = (re.compile('Basis vectors'),self._read_lattice_vectors)
        self.manage['lattice2']  = (re.compile('cubic'),self._read_lattice_vectors)
        self.manage['dynamical']  = (re.compile(' *Dynamical  Matrix in c'),self._read_dynamical)
        self.manage['epsilon']  = (re.compile(' *Dielectric Tensor:'),self._read_epsilon)
        self.manage['charges']  = (re.compile(' *Effective Charges E-U:'),self._read_born_charges)
        for f in self._outputfiles:
            self._ReadOutputFile(f)
        return

    def _read_header(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        self.nspecies = int(line.split()[0])
        self.nions    = int(line.split()[1])
        self.alat     = float(line.split()[3])
 
    def _read_epsilon(self,line):
        line = self.fd.readline()
        linea = self.fd.readline().split()
        self.zerof_optical_dielectric = []
        self.zerof_optical_dielectric.append( [float(f) for f in linea[0:3] ])
        linea = self.fd.readline().split()
        self.zerof_optical_dielectric.append( [float(f) for f in linea[0:3] ])
        linea = self.fd.readline().split()
        self.zerof_optical_dielectric.append( [float(f) for f in linea[0:3] ])
        return

    def _read_masses(self):
        self.mass_species = []
        self.species      = []
        for i in range(self.nspecies):
            linea = self.fd.readline().replace('\'','').split()
            self.species.append(linea[1])
            # The factor of two is because au in pwscf are half mass of electron
            self.mass_species.append(float(linea[2])*2/amu)
        self._read_cartesian_coordinates()
        return

    def _read_dynamical(self,line):
        nmodes = self.nions*3
        hessian = np.zeros( (nmodes,nmodes) )
        line = self.fd.readline()
        linea = self.fd.readline().split()
        # We only want to read the hessian at gamma
        q = [ float(q) for q in linea[3:6] ]
        qsum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2]
        if qsum > 0.0001:
            return
        # We read the hessian and store the mass weighted matrix
        linea = self.fd.readline().split()
        for a in range(self.nions) :
            for b in range(self.nions) :
               line = self.fd.readline()
               for ixyz in range(3):
                   ipos = a*3 + ixyz
                   linea = self.fd.readline().split()
                   for jxyz in range(3):
                       jpos = b*3 + jxyz
                       # factor of 0.5 'cos of au units in pwscf
                       hessian[ipos,jpos] = 0.5*float(linea[2*jxyz])/(amu*math.sqrt(self.masses[a]*self.masses[b]))
                   # end for jxyz
               # end for ixyz
            # end for b
        # end for a
        self._DynamicalMatrix(hessian)

    def _read_born_charges(self,line):
        self.born_charges = []
        line = self.fd.readline()
        for i in range(self.nions) :
          b = []
          line = self.fd.readline()
          line = self.fd.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          self.born_charges.append(b)
        if self.neutral:
            self._BornChargeSumRule()
        return


    def _read_lattice_vectors(self,line):
        linea = self.fd.readline().split()
        aVector = [ float(f)*self.alat/angs2bohr for f in linea[0:3] ]
        linea = self.fd.readline().split()
        bVector = [ float(f)*self.alat/angs2bohr for f in linea[0:3] ]
        linea = self.fd.readline().split()
        cVector = [ float(f)*self.alat/angs2bohr for f in linea[0:3] ]
        self.unitCells.append(UnitCell(aVector, bVector, cVector))
        self.ncells = len(self.unitCells)
        self.volume = self.unitCells[-1].volume
        self._read_masses()
        return

    def _read_cartesian_coordinates(self):
        self.ions = []
        self.masses = []
        for i in range(self.nions) :
            linea = self.fd.readline().split()
            species = int(linea[1])
            self.ions.append( [ float(linea[2]), float(linea[3]), float(linea[4]) ] )
            self.masses.append(self.mass_species[species-1])
        print "MASSES:",self.masses
        return

