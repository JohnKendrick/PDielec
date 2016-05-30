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
"""Generic reader for output files"""
import string
import re
import numpy as np
import math 
import os, sys
from Python.Constants import *
from Python.UnitCell import *
    
class GenericOutputReader:
    """Generic reader of output files.  Actual reader should inherit from this class"""

    def __init__(self,filenames):
        self._outputfiles            = filenames
        self.names                   = [ os.path.abspath(f) for f in filenames ] 
        self.debug                   = False
        self.type                    = 'Unkown'
        self.ncells                  = 0
        self.nsteps                  = 0
        self.formula                 = None
        self.nelect                  = 0
        self.volume                  = 0.0
        self.nions                   = 0
        self.nspecies                = 0
        self.final_free_energy       = 0
        self.unitCells               = []
        self.born_charges            = []
        self.species                 = []
        self.manage                  = {}
        # this in epsilon infinity
        self.zerof_optical_dielectric= [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric = [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
        self.elastic_constants       = [ [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0] ]
        self.frequencies             = []
        self.mass_weighted_normal_modes = []
        self.masses                  = []
        self.ions                    = []
        self.eckart = True
        # The child will issue a call to _ReadOutputfiles

    def PrintInfo (self):
        # Generic printing of information
        print "Number of atoms: ", self.nions
        print "Number of species: ", self.nspecies
        print "Frequencies: ", self.frequencies
        print "Masses: ", self.masses
        print "Born Charge: ", self.born_charges
        print "Epsilon inf: ", self.zerof_optical_dielectric
        print "Volume of cell: ", self.volume
        print "Unit cell: ", self.unitCells[-1].lattice
        mtotal = 0.0
        for m in self.masses :
           mtotal = mtotal + m
        print "Total mass is: ", mtotal, " g/mol"
        print "Density is: ", mtotal/( avogadro_si * self.volume * 1.0e-24) , " g/cc"
        return


    def _ReadOutputFiles(self):
        """Read the file names"""
        # Define the search keys to be looked for in the files
        print "Error _ReadOutputFiles must be defined by the actual file reader"
        return

    def _ReadOutputFile(self,name):
        """Read through the files for key words.  The keywords are established in _ReadOutputFiles"""
        self.fd = open(name,'r')
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
        return

    def _DynamicalMatrix(self,hessian):
        # Process the diagonal matrix
        # Hessian is a nxn matrix of the mass weighted force constants
        # The hessian is symmetrised
        # Translational modes are projected out
        # The hessian is diagonalised
        # Finally the frequencies and normal modes are stored
        nmodes = self.nions*3
        # symmetrise
        hessian = 0.5 * (hessian + hessian.T)
        # diagonalise
        eig_val, eig_vec = np.linalg.eigh(hessian)
        print "Eigen values of hessian"
        print eig_val
        # Project out the translational modes if requested
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
        # Form the projection operators
        P1 = unit - np.outer(p1,p1)
        P2 = unit - np.outer(p2,p2)
        P3 = unit - np.outer(p3,p3)
        if self.eckart:
            # Now project out
            print "Hessian will be modified by projecting out pure translation"
            hessian = np.dot ( np.dot(P1.T, hessian), P1)
            hessian = np.dot ( np.dot(P2.T, hessian), P2)
            hessian = np.dot ( np.dot(P3.T, hessian), P3)
        else:
            print "No projection of the hessian has been performed (NOECKART)"

        # diagonalise
        eig_val, eig_vec = np.linalg.eigh(hessian)
        print "Eigen values after projection"
        print eig_val
        # If eig_val has negative values then we store the negative frequency
        # convert to cm-1
        frequencies = np.zeros( nmodes )
        for i,eig in enumerate(eig_val):
            if eig < 0 :
                frequencies[i] = -math.sqrt(-eig) / wavenumber
            else :
                frequencies[i] = math.sqrt( eig) / wavenumber
            # end if
        #end for
        #jk print "Eigenvalues"
        #jk print ( ''.join('{:15.8}'.format(e) for e in eig_val) )
        #jk print "Frequencies"
        #jk print ( ''.join('{:15.8}'.format(f) for f in frequencies) )
        self.mass_weighted_normal_modes = []
        self.frequencies = frequencies.tolist()
        # Store the mass weighted normal modes
        for i in range(nmodes):
           mode = []
           n = 0
           for j in range(self.nions):
             ma = [ eig_vec[n][i], eig_vec[n+1][i], eig_vec[n+2][i] ]
             n = n + 3
             mode.append(ma)
           self.mass_weighted_normal_modes.append(mode)
        # end for i
        return
 

