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
        self._ReadOutputFiles() 

    def _ReadOutputFiles(self):
        """Read the AbInit file names"""
        # Define the search keys to be looked for in the files
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['header']   = (re.compile('Dynamical matrix file'),self._read_header)
        self.manage['lattice']  = (re.compile('Basis vectors'),self._read_lattice_vectors)
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

    def _read_masses(self,line):
        self.mass_species = []
        self.species      = []
        for i in range(self.nspecies):
            linea = self.fd.readline().replace('\'').split()
            self.species.append(linea[1])
            self.mass_species.append(float(linea[2])/amu)
        self._read_cartesian_coordinates()
        return

    def _read_eigenvectors(self,line):
        if os.path.isfile("HESSFREQ.DAT"):
            print "Reading Hessian from HESSFREQ.DAT"
            self._read_hessfreq_dat();
        else:
            print "Reading Normal Modes from output file"
            print "WARNING! WARNING! WARNING! WARNING! WARNING!"
            print "The numerical precision of the input is limited"
            print "WARNING! WARNING! WARNING! WARNING! WARNING!"
            self._read_output_eigenvectors(line)
        return

    def _read_dynamical(self):
        nmodes = self.nions*3
        hessian = np.zeros( (nmodes,nmodes)
        line = self.fd.readline()
        linea = self.fd.readline().split()
        # We only want to read the hessian at gamma
        q = [ q for q in linea[3:6] ]
        qsum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2]
        if qsum > 0.0001:
            return
        # We read the hessian and store the mass weighted matrix
        linea = self.fd.readline().split()
        for a in range(self.nions) :
            for b in range(self.nions) :
               linea = self.fd.readline().split()
               for ixyz in range(3):
                   ipos = (a-1)*3 + ixyz
                   for jxyz in range(3):
                       ipos = (b-1)*3 + jxyz
                       linea = self.fd.readline().split()
                       hessian[ipos,jpos] = self.masses[a]*self.masses[b]*float(linea[2*jxyz])
                   # end for jxyz
               # end for ixyz
            # end for b
        # end for a
        self._DynamicalMatrix(hessian)

    def _read_born_charges(self):
        self.born_charges = []
        for i in range(self.nions) :
          b = []
          line = self.fd.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          self.born_charges.append(b)
        return


    def _read_lattice_vectors(self,line):
        line = self.fd.readline()
        aVector = [ float(f)*alat for f in linea[0:3] ]
        linea = self.fd.readline().split()
        bVector = [ float(f)*alat for f in linea[0:3] ]
        linea = self.fd.readline().split()
        cVector = [ float(f)*alat for f in linea[0:3] ]
        self.unitCells.append(UnitCell(aVector, bVector, cVector))
        self.ncells = len(self.unitCells)
        self.volume = self.unitCells[-1].volume
        self._read_masses()
        return

    def _read_cartesian_coordinates(self,line):
        self.ions = []
        self.masses = []
        for i in range(self.nions) :
            linea = self.fd.readline().split()
            species = int(linea[1])
            self.ions.append( [ float(linea[2]), float(linea[3]), float(linea[4]) ] )
            self.masses.append(mass_species[species-1])
        return

