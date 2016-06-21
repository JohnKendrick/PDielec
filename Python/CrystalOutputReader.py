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
"""Read the Crystal output files"""
import string
import re
import numpy as np
import math 
import os, sys
from Python.Constants import *
from Python.UnitCell import *
from Python.GenericOutputReader import *
    
class CrystalOutputReader(GenericOutputReader):
    """Read contents of a directory containing Crystal input and output files"""

    def __init__(self,filenames):
        GenericOutputReader.__init__(self, filenames)
        self.type                    = 'Crystal output'
        self.species                 = []
        return

    def _ReadOutputFiles(self):
        """Read the Crystal files in the directory"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['masses']   = (re.compile(' ATOMS ISOTOPIC MASS'),self._read_masses)
        self.manage['lattice']  = (re.compile(' DIRECT LATTICE VECTORS CARTESIAN COMPONENTS'),self._read_lattice_vectors)
        self.manage['fractional']  = (re.compile(' ATOMS IN THE ASYMMETRIC UNIT'),self._read_fractional_coordinates)
        self.manage['bornCharges']  = (re.compile(' ATOMIC BORN CHARGE TENSOR'),self._read_born_charges)
        self.manage['eigenvectors']  = (re.compile(' NORMAL MODES NORMALIZ'),self._read_eigenvectors)
        self.manage['staticIonic']  = (re.compile(' SUM TENSOR OF THE VIBRATIONAL CONTRIBUTIONS TO '),self._read_ionic_dielectric)
        self.manage['noeckart']  = (re.compile('.* REMOVING ECKART CONDITIONS'),self._read_eckart)
        self.manage['epsilon']  = (re.compile(' SUSCEPTIBILITY '),self._read_epsilon)
        for f in self._outputfiles:
            self._ReadOutputFile(f)
        return

    def _read_epsilon(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        optical_dielectric= [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
        optical_dielectric[0][0] = float(self.fd.readline().split()[3])
        optical_dielectric[0][1] = float(self.fd.readline().split()[3])
        optical_dielectric[0][2] = float(self.fd.readline().split()[3])
        optical_dielectric[1][1] = float(self.fd.readline().split()[3])
        optical_dielectric[1][2] = float(self.fd.readline().split()[3])
        optical_dielectric[2][2] = float(self.fd.readline().split()[3])
        optical_dielectric[1][0] = optical_dielectric[0][1]
        optical_dielectric[2][0] = optical_dielectric[0][2] 
        optical_dielectric[2][1] = optical_dielectric[1][2] 
        self.zerof_optical_dielectric = optical_dielectric
        return

    def _read_eckart(self,line):
        self.eckart = False
        return

    def _read_masses(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        n = 2
        for i in range(self.nions):
          if n > 11 :
            line = self.fd.readline()
            n = 2 
          self.masses.append ( float(line.split()[n]) ) 
          n = n + 3
        return

    def _read_eigenvectors(self,line):
        if os.path.isfile("HESSFREQ.DAT"):
            print("Reading Hessian from HESSFREQ.DAT")
            self._read_hessfreq_dat();
        else:
            print("Reading Normal Modes from output file")
            print("WARNING! WARNING! WARNING! WARNING! WARNING!")
            print("The numerical precision of the input is limited")
            print("If possible proved a HESSFREQ.DAT from the calculation")
            print("WARNING! WARNING! WARNING! WARNING! WARNING!")
            self._read_output_eigenvectors(line)
        return

    def _read_hessfreq_dat(self):
        fd2 = open("HESSFREQ.DAT",'r')
        nmodes = self.nions*3
        # Create a mass weighting vector
        n = 0
        massweight = np.zeros( nmodes )
        frequencies = np.zeros( nmodes )
        for a in range(self.nions):
            for j in range(3):
                massweight[n] = 1.0 / math.sqrt(self.masses[a]*amu)
                n = n + 1
        # We read the hessian and store the mass weighted matrix
        hessian = np.zeros( (nmodes,nmodes) )
        line = fd2.readline()
        pos = 0 
        for i in range(nmodes) :
            for j in range(nmodes) :
              hessian[i,j] = massweight[i]*massweight[j]*float(line.split()[pos])
              pos = pos + 1
              if pos >= 4:
                  line = fd2.readline()
                  pos = 0
              # end if pos
            # end for j
        # end for i
        fd2.close()
        # symmetrise, project, diagonalise and store the frequencies and normal modes
        self._DynamicalMatrix(hessian)
        return

    def _read_output_eigenvectors(self,line):
        line = self.fd.readline()
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        n = 6
        nmodes = 3 * self.nions
        self.frequencies = []
        # Loop over all the modes
        for j in range(nmodes):
          # They come in columns of 6 so after 6 modes we read the frequencies and initialse the 6 new modes
          n = n + 1
          if n >= 6:
            n = 0
            linef = self.fd.readline().split()[1:]
            linef = [ float(f) for f in linef ]
            self.fd.readline()
            atoms = []
            for f in linef:
              self.frequencies.append(f)
              atoms.append([])
            # Read through the XYZ components for each atom and store the mode in atoms
            for i in range(self.nions):
              mass = np.sqrt( self.masses[i] )
              linex = self.fd.readline().split()[4:]
              liney = self.fd.readline().split()[1:]
              linez = self.fd.readline().split()[1:]
              for a,x,y,z in zip(atoms,linex,liney,linez):
                x = float(x) * mass
                y = float(y) * mass
                z = float(z) * mass
                a.append([x,y,z])
            for a in atoms:
              self.mass_weighted_normal_modes.append(a)
            line = self.fd.readline()
        
        for i,mode in enumerate(self.mass_weighted_normal_modes):
          marray = np.array( mode )
          sum = 0.0
          for atoms in marray:
            for xyz in atoms:
              sum = sum + xyz*xyz
          marray = marray / np.sqrt(sum)
          self.mass_weighted_normal_modes[i] = marray.tolist()
            
        return

    def _read_born_charges(self,line):
        """Read the born charges from the outputfile file.  
           Each row of the output refers to a given field direction
           Each column in the row refers the atomic displacement 
           so the output is arranged [ [ a1x a1y a1z ] 
                                       [ a2x a2y a2z ] 
                                       [ a3x a3y a3z ]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        if os.path.isfile("BORN.DAT"):
            print("Reading Born charge tensor from BORN.DAT")
            self._read_born_charges_from_born_dat();
        else:
            print("Reading Born Charge Tensor from output file")
            print("WARNING! WARNING! WARNING! WARNING! WARNING!")
            print("The numerical precision of the input is limited")
            print("If possible proved a BORN.DAT from the calculation")
            print("WARNING! WARNING! WARNING! WARNING! WARNING!")
            self._read_born_charges_from_output(line)
        if self.neutral:
            self._BornChargeSumRule()
        return

    def _read_born_charges_from_born_dat(self):
        fd2 = open("BORN.DAT",'r')
        self.born_charges = []
        for i in range(self.nions) :
          b = []
          line = fd2.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          line = fd2.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          line = fd2.readline()
          b.append( [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ] )
          self.born_charges.append(b)
        fd2.close()
        return


    def _read_born_charges_from_output(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        self.born_charges = []
        for i in range(self.nions):
          line = self.fd.readline()
          line = self.fd.readline()
          line = self.fd.readline()
          b = []
          b.append( [ float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ] )
          line = self.fd.readline()
          b.append( [ float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ] )
          line = self.fd.readline()
          self.born_charges.append(b)
          line = self.fd.readline()
        return

    def _read_ionic_dielectric(self,line):
        # Read the ionic contribution to the static dielectric 
        line = self.fd.readline()
        line = self.fd.readline()
        ionic_dielectric = []
        ionic_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        line = self.fd.readline()
        ionic_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        line = self.fd.readline()
        ionic_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        a2 = np.array(ionic_dielectric)
        self.zerof_static_dielectric = a2.tolist()
        return

    def _read_lattice_vectors(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        aVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        line = self.fd.readline()
        bVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        line = self.fd.readline()
        cVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        self.unitCells.append(UnitCell(aVector, bVector, cVector))
        self.ncells = len(self.unitCells)
        self.volume = self.unitCells[-1].volume
        # The fractional coordinates are specified before the lattice vectors
        self.unitCells[-1].fractionalCoordinates(self.ions)
        return

    def _read_fractional_coordinates(self,line):
        self.nions = int(line.split()[12])
        line = self.fd.readline()
        line = self.fd.readline()
        self.ions = []
        for i in range(self.nions) :
            line = self.fd.readline()
            self.species.append( line.split()[3] )
            self.ions.append( [ float(line.split()[4]), float(line.split()[5]), float(line.split()[6]) ] )
        return

