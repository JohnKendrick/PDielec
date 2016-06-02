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
"""Read the contents of a directory containing Castep output and create a csv style file of information"""
import string
import re
import numpy as np
import math 
import os, sys
from Python.UnitCell import *
from Python.Constants import *
from Python.GenericOutputReader import *

class CastepOutputReader(GenericOutputReader):
    """Read the contents of a Castep output file
       Inherit the following from the GenericOutputReader
       __init__
       printInfo
       _ReadOutputFile
       _dyanmicalMatrix
"""

    def __init__(self,filenames):
        GenericOutputReader.__init__(self, filenames)
        if filenames[0].find(".castep"):
            seedname,ext = os.path.splitext(filenames[0])
        elif filenames[0].find(".phonon"):
            seedname,ext = os.path.splitext(filenames[0])
        self._castepfile             = seedname+".castep"
        self._phononfile             = seedname+".phonon"
        self.names                   = [ self._castepfile, self._phononfile ]
        self._outputfiles            = [ self._castepfile, self._phononfile ]
        self.type                    = 'Castep output'
        # Specific Castep Reader Variables
        self.pspots                  = {}
        self.spin                    = 0
        self.encut                   = 0.0
        self.ediff                   = 0.0
        self.epsilon                 = None
        self.nkpts                   = 0
        self.nbands                  = 0
        self._nbranches              = 0
        self.species                 = []
        self.final_energy_without_entropy = 0.0
        self.pressure                = None
        self._pulay                   = None
        self.masses_dictionary       = {}
        self.ions_per_type           = []
        self._ion_type_index         = {}
        self._ion_index_type         = {}
        self._ReadOutputFiles() 

    def _ReadOutputFiles(self):
        """For the .castep file"""
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['spin']          = (re.compile(' *net spin   of'),       self._read_spin)
        self.manage['nelect']        = (re.compile(' *number of  electrons'),self._read_nelect)
        self.manage['cellcontents']  = (re.compile(' *Unit Cell'),self._read_cellcontents)
        self.manage['pspots']        = (re.compile(' *Files used for pseudopotentials:'),self._read_pspot)
        self.manage['masses']        = (re.compile(' *Mass of species in AMU'),self._read_masses)
        self.manage['kpts']          = (re.compile(' *Number of kpoints used'),self._read_kpoints)
        self.manage['finalenergy']   = (re.compile(' *Final energy, E'),self._read_energies)
        self.manage['encut']         = (re.compile(' *plane wave basis set cut'),self._read_encut)
        self.manage['nbands']        = (re.compile(' *number of bands'),self._read_nbands)
        self.manage['pressure']      = (re.compile(' *\* Pressure: '),self._read_external_pressure)
        self.manage['opticalDielectric']  = (re.compile(' *Optical Permittivity'),self._read_dielectric)
        self.manage['bornCharges']  = (re.compile(' *Born Effective Charges'),self._read_born_charges)
        """For the .phonon file"""
        self.manage['frequency']  = (re.compile('     q-pt=    1    0.000000  0.000000  0.000000      1.0000000000 *$'),self._read_frequencies)
        self.manage['nbranches']  = (re.compile(' Number of branches'),self._read_nbranches)
        for f in self._outputfiles:
            self._ReadOutputFile(f)
        return

    def _read_nbranches(self,line):
        # phonon file being read
        self._nbranches = int(line.split()[3])

    def _read_frequencies(self,line):
        # phonon file being read
        self.mass_weighted_normal_modes = []
        # maybe should use _nbranches
        frequencies = []
        intensities = []
        normal_modes = []
        for imode in range(self._nbranches):
            line = self.fd.readline()
            frequencies.append(float(line.split()[1]))
            intensities.append(float(line.split()[2]))

        line = self.fd.readline()
        line = self.fd.readline()
        for imode in range(self._nbranches):
            a = []
            for ion in range(self.nions):
                line = self.fd.readline()
                a.append([ float(line.split()[2]), float(line.split()[4]), float(line.split()[6])] )
            # end for ion
            normal_modes.append(a)
        # end for imode
        self.frequencies = []
        self.intensities = []
        self.mass_weighted_normal_modes = []
        # now reads all frequencies imaginary or not
        for i,freq in enumerate(frequencies):
             if freq < 0.0 :
                 self.frequencies.append(complex(0,-freq))
             else:
                 self.frequencies.append(freq)
             self.intensities.append(intensities[i])
             self.mass_weighted_normal_modes.append(normal_modes[i])
            # end of if freq
        # end of for freq
        return

    def _ReadCastep(self):
        """Read the .castep file"""
        self.fd = open(self._castepfile,'r')
        self.manage = {}   # Empty the dictionary matching phrases
        #jk self.manage['pseudoatom']    = (re.compile(' Pseudo atomic calculation'),  self._read_pseudoatom)
        self.manage['spin']          = (re.compile(' *net spin   of'),       self._read_spin)
        self.manage['nelect']        = (re.compile(' *number of  electrons'),self._read_nelect)
        self.manage['cellcontents']  = (re.compile(' *Unit Cell'),self._read_cellcontents)
        self.manage['pspots']        = (re.compile(' *Files used for pseudopotentials:'),self._read_pspot)
        self.manage['masses']        = (re.compile(' *Mass of species in AMU'),self._read_masses)
        self.manage['kpts']          = (re.compile(' *Number of kpoints used'),self._read_kpoints)
        self.manage['finalenergy']   = (re.compile(' *Final energy, E'),self._read_energies)
        self.manage['encut']         = (re.compile(' *plane wave basis set cut'),self._read_encut)
        self.manage['nbands']        = (re.compile(' *number of bands'),self._read_nbands)
        self.manage['pressure']      = (re.compile(' *\* Pressure: '),self._read_external_pressure)
        self.manage['opticalDielectric']  = (re.compile(' *Optical Permittivity'),self._read_dielectric)
        self.manage['bornCharges']  = (re.compile(' *Born Effective Charges'),self._read_born_charges)
        # Loop through the contents of the file a line at a time and parse the contents
        line = self.fd.readline()
        while line != '' :
            for k in self.manage.keys():
                if self.manage[k][0].match(line): 
                    method   = self.manage[k][1]
                    if self.debug:
                        print >> sys.stderr, 'Match found %s' % k
                        print >> sys.stderr, self.manage[k]
                    method(line)
                    break
                #end if
            #end for
            line = self.fd.readline()
        #end while
        self.fd.close()

    def _read_kpoints(self,line):
        self.nkpts = int(line.split()[5])
        return

    def _read_nbands(self,line):
        self.nbands = int(line.split()[4])
        return

    def _read_pseudoatom(self,line):
        species = line.split()[5]
        # These are two private dictionary to map the species name to a type index
        self.species.append(species)
        self._ion_type_index[species] = self.nspecies
        self._ion_index_type[self.nspecies] = species
        self.nspecies += 1
        return

    def _read_cellcontents(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        aVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        line = self.fd.readline()
        bVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        line = self.fd.readline()
        cVector = [ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ]
        self.unitCells.append(UnitCell(aVector, bVector, cVector))
        self.ncells = len(self.unitCells)
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        self.volume = float(line.split()[4])
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        if len(self.unitCells) == 1 :
            self.nions=int(line.split()[7])
            line = self.fd.readline()
            self.nspecies=int(line.split()[7]) 
            line = self.fd.readline()
            line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        fractional_coordinates=[]
        # ions_per_type is a dictionary here and is local
        ions_per_type = {}
        if len(self.unitCells) == 1 :
            self.species = []
            for i in range(self.nions) :
                line = self.fd.readline()
                atom_type = line.split()[1]
                if atom_type not in ions_per_type :
                    self.species.append(atom_type)
                    ions_per_type[atom_type] = 0
                ions_per_type[atom_type] += 1
                atom_frac = [ float(f) for f in (line.split()[3:6]) ]
                fractional_coordinates.append(atom_frac)
            # At some point we should store the fractional coordinates
            for species in self.species :
                n = ions_per_type[species]
                # self.ions_per_type is a list
                self.ions_per_type.append(n)
        self.unitCells[-1].fractionalCoordinates(fractional_coordinates)
        return

    def _read_masses(self,line):
        self.masses = []
        for i in range(self.nspecies) :
          line = self.fd.readline()
          species = line.split()[0]
          mass = float(line.split()[1])
          self.masses_dictionary[species] = mass
          nions = self.ions_per_type[i]
          for j in range(nions):
            self.masses.append(mass)
        return

    def _read_born_charges(self,line):
        """Read the born charges from the castep file.  
           Each column of the output refers to a given field direction
           Each row in the row refers the atomic displacement 
           So the numbers in the output are arrange  [ [ a1x a2x a3x ]  
                                                       [ a1y a2y a3y ] 
                                                       [ a1z a2z a3z ]]
           The output tensor needs them arranged     [ [ a1x a1y a1z ]
                                                       [ a2x a2y a2z ]
                                                       [ a3x a3y a3z ]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        line = self.fd.readline()
        self.born_charges = []
        for i in range(self.nions) :
          b = []
          line = self.fd.readline()
          b.append( [ float(f) for f in line.split()[2:5] ] )
          line = self.fd.readline()
          b.append( [ float(f) for f in line.split()[0:3] ] )
          line = self.fd.readline()
          b.append( [ float(f) for f in line.split()[0:3] ] )
          B = np.array(b)
          C = B.T
          self.born_charges.append(C.tolist())
        if self.neutral:
            self._BornChargeSumRule()
        return

    def _read_dielectric(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        # this is epsilon infinity
        self.zerof_optical_dielectric = []
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric = []
        self.zerof_optical_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        self.zerof_static_dielectric.append(  [ float(f) for f in line.split()[3:6] ] )
        line = self.fd.readline()
        self.zerof_optical_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        self.zerof_static_dielectric.append(  [ float(f) for f in line.split()[3:6] ] )
        line = self.fd.readline()
        self.zerof_optical_dielectric.append( [ float(f) for f in line.split()[0:3] ] )
        self.zerof_static_dielectric.append(  [ float(f) for f in line.split()[3:6] ] )
        return

    def _read_skip4(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        line = self.fd.readline()
        return


    def _read_external_pressure(self,line):
        self.pressure = float(line.split()[2])
        return

    def _read_pspot(self,line):
        for i in range(self.nspecies) :
          line = self.fd.readline()
          self.pspots[line.split()[0]] = line.split()[1]
        return

    def _read_spin(self,line):
        self.spin = int(float(line.split()[5]))
        return
    def _read_encut(self,line):
        self.encut = float(line.split()[6]) 
        return
    def _read_ediff(self,line):
        self.ediff = float(line.split()[2]) 
        return
    def _read_nelect(self,line):
        self.nelect = int(float(line.split()[4]) )
        return
    def _read_epsilon(self,line):
        var = line.split()[2].lower()
        self.epsilon = False
        if var == '.true.' or var == 't' or var == 'true' :
          self.epsilon = True
        return
    def _read_energies(self,line):
        self.final_energy_without_entropy = float(line.split()[4])
        line = self.fd.readline()
        self.final_free_energy = float(line.split()[5])
        #jk line = self.fd.readline()
        #jk line = self.fd.readline()
        #jk self.final_0K_energy = float(line.split()[5])
        return
