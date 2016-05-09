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
    
class CrystalOutputReader:
    """Read the contents of a directory containg Crystal input and output files"""

    def __init__(self,filename):
        self._outputfile             = filename
        self.name                    = os.path.abspath(self._outputfile)
        self.debug                   = False
        self.title                   = None
        self.type                    = 'Crystal output'
        self.ncells                  = 0
        self.nsteps                  = 0
        self.formula                 = None
        self.nelect                  = 0
        self.volume                  = 0.0
        self.nions                   = 0
        self.nspecies                = 0
        self.final_energy_without_entropy = 0
        self.final_free_energy       = 0
        self.unitCells               = []
        self.born_charges            = []
        self.species                 = []
        # this in epsilon infinity
        self.zerof_optical_dielectric= [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric = [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
        self.elastic_constants       = [ [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0] ]
        self.frequencies             = []
        self.mass_weighted_normal_modes = []
        self.masses                  = []
        self.ions_per_type             = []
        self.epsilon                = None
        self._ReadOutputFile() 

    def printInfo (self):
        print "Number of atoms: ", self.nions
        print "Number of species: ", self.nspecies
        print "Frequencies: ", self.frequencies
        print "Masses: ", self.masses
        print "Born Charge: ", self.born_charges
        print "Epsilon inf: ", self.zerof_optical_dielectric
        print "Volume of cell: ", self.volume
        print "Unit cell: ", self.unitCells[-1].lattice
        print "Fractional Coordinates: ", self.unitCells[-1].fractional_coordinates
        mtotal = 0.0
        for m in self.masses :
           mtotal = mtotal + m
        print "Total mass is: ", mtotal, " g/mol"
        print "Density is: ", mtotal/( avogadro_si * self.volume * 1.0e-24) , " g/cc"
        return


    def _ReadOutputFile(self):
        """Read the Crystal files in the directory"""
        self.fd = open(self._outputfile,'r')
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['masses']   = (re.compile(' ATOMS ISOTOPIC MASS'),self._read_masses)
        self.manage['lattice']  = (re.compile(' DIRECT LATTICE VECTORS CARTESIAN COMPONENTS'),self._read_lattice_vectors)
        self.manage['fractional']  = (re.compile(' ATOMS IN THE ASYMMETRIC UNIT'),self._read_fractional_coordinates)
        self.manage['bornCharges']  = (re.compile(' ATOMIC BORN CHARGE TENSOR'),self._read_born_charges)
        self.manage['eigenvectors']  = (re.compile(' NORMAL MODES NORMALIZ'),self._read_eigenvectors)
        self.manage['staticIonic']  = (re.compile(' SUM TENSOR OF THE VIBRATIONAL CONTRIBUTIONS TO '),self._read_ionic_dielectric)
 
        # Loop through the contents of the file a line at a time and parse the contents
        line = self.fd.readline()
        while line != '' :
            for k in self.manage.keys():
                if self.manage[k][0].match(line): 
                    method   = self.manage[k][1]
                    if self.debug:
                        print 'Match found %s' % k
                        print self.manage[k]
                    method(line)
                    break
                #end if
            #end for
            line = self.fd.readline()
        #end while
        self.fd.close()

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
            print "Reading Hessian from HESSFREQ.DAT"
            self._read_hessfreq_dat();
        else:
            print "Reading Normal Modes from output file"
            print "WARNING! WARNING! WARNING! WARNING! WARNING!"
            print "The numerical precision of the input is limited"
            print "WARNING! WARNING! WARNING! WARNING! WARNING!"
            self._read_output_eigenvectors(line)
        return

    def _read_hessfreq_dat(self):
        fd2 = open("HESSFREQ.DAT",'r')
        nmodes = self.nions*3
        # Create a mass weighting vector
        n = 0
        massweight = np.zeros( nmodes )
        frequencies = np.zeros( nmodes )
        print "MASSES", self.masses
        print "NIONS", self.nions
        for a in range(self.nions):
            for j in range(3):
                massweight[n] = 1.0 / math.sqrt(amu*self.masses[a])
                n = n + 1
        print "Mass weighting"
        print massweight
        # We read the hessian and store the mass weighted matrix
        hessian = np.zeros( (nmodes,nmodes) )
        self.born_charges = []
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
        # Symmetrise the hessian
        hessian = 0.5 * ( hessian + hessian.T )
        # Find the eigen vectors and normal modes
        eig_val1, eig_vec1 = np.linalg.eigh(hessian)
        print "Eigenvalues before projection"
        print eig_val1
        # Look at the sum rule for the hessian, taken from equation 2 of 
        # Practical methods in abinitio  Lattice Dynamics by G J Ackland, M C Warren and S J Clark
        # J. Phys. Condens. Natter 9 (1997) 7861-7872.
        for i in range(self.nions):
            for alpha in range(3):
                ialpha = i*3 + alpha
                for beta in range(3):
                    ibeta = i*3 + beta 
                    sum = 0.0
                    for j in range(self.nions):
                        jbeta = j*3 + beta 
                        sum = sum + hessian[ialpha][jbeta]
                    #end for j
                hessian[ialpha,ibeta] = hessian[ialpha,ibeta]-sum
                #end for beta
            #end for alpha
        #end for i 
        # Find the eigen vectors and normal modes
        eig_val3, eig_vec3 = np.linalg.eigh(hessian)
        print "Eigenvalues after sum rule"
        print eig_val3
        # Project out the translational modes
        unit = np.eye( nmodes )
        p1 = np.zeros( nmodes )
        p2 = np.zeros( nmodes )
        p3 = np.zeros( nmodes )
        for i in range(self.nions):
            p1[i*3+0] = math.sqrt(self.masses[i])
            p2[i*3+1] = math.sqrt(self.masses[i])
            p3[i*3+2] = math.sqrt(self.masses[i])
        # end for i
        # Normalise
        p1 = p1 / math.sqrt(np.dot(p1,p1))
        p2 = p2 / math.sqrt(np.dot(p2,p2))
        p3 = p3 / math.sqrt(np.dot(p3,p3))
        P1 = unit - np.outer(p1,p1)
        P2 = unit - np.outer(p2,p2)
        P3 = unit - np.outer(p3,p3)
        # Now project out
        hessian = np.dot ( np.dot(P1.T, hessian), P1)
        hessian = np.dot ( np.dot(P2.T, hessian), P2)
        hessian = np.dot ( np.dot(P3.T, hessian), P3)
        # Find the eigen vectors and normal modes
        eig_val2, eig_vec2 = np.linalg.eigh(hessian)
        eig_val2 = eig_val2 
        # If eig_val has negative values then we store the negative frequency
        # convert to cm-1
        for i,eig in enumerate(eig_val1):
            if eig < 0 :
                frequencies[i] = -math.sqrt(-eig) / wavenumber
            else :
                frequencies[i] = math.sqrt( eig) / wavenumber
            # end if
        #end for
        print "Eigenvalues after projection", eig_val2
        print "Eigenvalues changes projection", eig_val1-eig_val2
        print "Frequencies", frequencies
        self.mass_weighted_normal_modes = []
        self.frequencies = frequencies.tolist()
        for i in range(nmodes):
           mode = []
           n = 0
           for j in range(self.nions):
             ma = [ eig_vec1[i][n], eig_vec1[i][n+1], eig_vec1[i][n+2] ]
             n = n + 3
             mode.append(ma)
           self.mass_weighted_normal_modes.append(mode)
        print "Frequencies", self.frequencies
        # end for i
        fd2.close()

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
            print "Reading Born Charge Tensor from BORN.DAT"
            self._read_born_charges_from_born_dat();
        else:
            print "Reading Born Charge Tensor from output file"
            print "WARNING! WARNING! WARNING! WARNING! WARNING!"
            print "The numerical precision of the input is limited"
            print "WARNING! WARNING! WARNING! WARNING! WARNING!"
            self._read_born_charges_from_output(line)
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
        return

    def _read_fractional_coordinates(self,line):
        self.nions = int(line.split()[12])
        line = self.fd.readline()
        line = self.fd.readline()
        ions = []
        for i in range(self.nions) :
            line = self.fd.readline()
            self.species.append( line.split()[3] )
            ions.append( [ float(line.split()[4]), float(line.split()[5]), float(line.split()[6]) ] )
        self.unitCells[-1].fractionalCoordinates(ions)
        return

