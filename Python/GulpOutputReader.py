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
"""Read the contents Gulp output and create a csv style file of information"""
import string
import re
import numpy as np
import math 
import os, sys
from Python.UnitCell import *
from Python.Constants import *

class GulpOutputReader:
    """Read the contents of a Gulp output file
       It is quite difficult to work out what is a suitable set of commands for gulp
       The following seems to work OK, important are 
       phon - a phonon calculation
       intens - calculate the IR intensities
       eigen - print the eigen vectors of the dynamical matrix
       cart - print the cartesian coordinates (this is the only way of see what all the atom types are)
       nosym - make sure the calculation is done in P1
       qeq molq optimise conp qok nomodcoord prop phon intens eigen cart
    """

    def __init__(self,name):
        self._gulpfile               = name
        self.name                    = os.path.abspath(self._gulpfile)
        self.debug                   = False
        self.title                   = None
        self.type                    = 'Gulp output'
        self.ncells                  = 0
        self.nsteps                  = 0
        self.formula                 = None
        self.volume                  = 0.0
        self.nions                   = 0
        self.shells                  = 0
        self.ncores                  = 0
        self.nspecies                = 0
        self.species                 = []
        self.cartesian_coordinates   = []
        self.fractional_coordinates  = []
        self.final_energy_without_entropy = 0.0
        self.final_free_energy       = 0.0
        self.pressure                = None
        self._pulay                   = None
        self.unitCells               = []
        self.born_charges            = []
        # this is epsilon infinity
        self.zerof_optical_dielectric= [ [1,0,0], [0,1,0], [0,0,1] ]
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric = []
        self.frequencies             = []
        self.mass_weighted_normal_modes = []
        self.masses                  = []
        self.mass_dictionary         = {}
        self.ions_per_type           = []
        self._ion_type_index         = {}
        self._ion_index_type         = {}
        self._ReadGulpOutput() 

    def _ReadGulpOutput(self):
        """Read the .gout file"""
        self.fd = open(self._gulpfile,'r')
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage['formula']  = (re.compile(' *Formula'),self._read_formula)
        self.manage['nions']  = (re.compile(' *Total number atoms'),self._read_total_number_of_atoms)
        self.manage['nions_irreducible']  = (re.compile(' *Number of irreducible atoms'),self._read_number_of_irreducible_atoms)
        self.manage['cellcontents']  = (re.compile(' *Final Cartesian coor'),self._read_cellcontents)
        self.manage['cellcontentsf']  = (re.compile(' *Final fractional coor'),self._read_cellcontentsf)
        self.manage['species']        = (re.compile(' *Species output for all configurations'),self._read_species)
        self.manage['lattice']          = (re.compile(' *Final Cartesian lattice'),self._read_lattice)
        self.manage['finalenergy']   = (re.compile(' *Total lattice energy'),self._read_energies)
        self.manage['pressure']      = (re.compile(' *Pressure of configuration'),self._read_external_pressure)
        self.manage['temperature']      = (re.compile(' *Temperature of configuration'),self._read_temperature)
        self.manage['staticDielectric']  = (re.compile(' *Static dielectric c'),self._read_static_dielectric)
        self.manage['opticalDielectric']  = (re.compile(' *High frequency diel'),self._read_optical_dielectric)
        self.manage['bornCharges']  = (re.compile(' *Born effective charge'),self._read_born_charges)
        self.manage['elasticConstants']  = (re.compile(' *Elastic Constant Matrix'),self._read_elastic_constants)
        self.manage['frequencies']  = (re.compile(' *Frequencies .cm-1.'),self._read_frequencies)
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

    def printInfo (self):
        print "Number of ions: ", self.nions
        print "Number of cores: ", self.ncores
        print "Number of shells: ", self.nshells
        print "Number of species: ", self.nspecies
        print "Frequencies: ", self.frequencies
        print "Masses: ", self.masses
        print "Born Charge: ", self.born_charges
        print "Epsilon inf: ", self.zerof_optical_dielectric
        print "Volume of cell: ", self.volume
        print "Unit cell: ", self.unitCells[-1].lattice
        print "Fractional coordinates: ", self.unitCells[-1].fractional_coordinates
        mtotal = 0.0
        for m in self.masses :
           mtotal = mtotal + m
        print "Total mass is: ", mtotal, " g/mol"
        print "Density is: ", mtotal/( avogadro_si * self.volume * 1.0e-24) , " g/cc"
        return

    def _read_formula(self,line):
        self.formula = line.split()[2]
        return
 
    def _read_elastic_constants(self,line):
        for skip in range(0,5):
            line = self.fd.readline()
        elastic = []
        for skip in range(0,6):
             elastic.append( [ float(f) for f in  line.split()[1:] ] )
             line = self.fd.readline()
        self.elastic_constant_tensor = np.array(elastic)
        print "Elastic constant tensor"
        print self.elastic_constant_tensor
        return
 
    def _read_frequencies(self,line):
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        line = self.fd.readline()
        nmodes = 0
        nions = self.ncores
        while nmodes < 3*nions :
            line = self.fd.readline()
            if line.split()[0] == "Note:":
                line = self.fd.readline()
                line = self.fd.readline()
            freq = [ float(f) for f in (line.split()[1:]) ]
            ncolumns = len(freq)
            nmodes = nmodes + ncolumns
            self.frequencies += freq
            line = self.fd.readline() # IR
            line = self.fd.readline() # IR x
            line = self.fd.readline() # IR y
            line = self.fd.readline() # IR z
            line = self.fd.readline() # Raman
            line = self.fd.readline() 
            columns = []
            for n in range(0,ncolumns):
                columns.append( [] )
            for i in range(0,nions):
                line = self.fd.readline() 
                modex = [ float(f) for f in  line.split()[2:] ]
                line = self.fd.readline() 
                modey = [ float(f) for f in  line.split()[2:] ]
                line = self.fd.readline() 
                modez = [ float(f) for f in  line.split()[2:] ]
                n = 0
                for mx,my,mz in zip(modex,modey,modez):
                    columns[n].append( [mx, my, mz] )
                    n += 1
                 # end for mx,my,mz (columns)
            # end loop over atoms
            #print columns
            for mode in columns:
                 self.mass_weighted_normal_modes.append(mode)
            line = self.fd.readline()
            line = self.fd.readline()
        #jk print "Normal modes"
        #jk print self.mass_weighted_normal_modes
 
    def _read_total_number_of_atoms(self,line):
        self.nions = int(line.split()[4])
        #jk print "Nions",self.nions
 
    def _read_number_of_irreducible_atoms(self,line):
        self.nions_irreducible = int(line.split()[5])
        #jk print "Nions irreducible",self.nions_irreducible
 
    def _read_lattice(self,line):
        line = self.fd.readline()
        line = self.fd.readline()
        ax = float(line.split()[0])
        ay = float(line.split()[1])
        az = float(line.split()[2])
        line = self.fd.readline()
        bx = float(line.split()[0])
        by = float(line.split()[1])
        bz = float(line.split()[2])
        line = self.fd.readline()
        cx = float(line.split()[0])
        cy = float(line.split()[1])
        cz = float(line.split()[2])
        aVector = [ ax, ay, az ]
        bVector = [ bx, by, bz ]
        cVector = [ cx, cy, cz ]
        cell = UnitCell(aVector, bVector, cVector)
        self.unitCells.append(cell)
        self.volume = cell.volume
        self.ncells = len(self.unitCells)
        print "GULP: volume = ", self.volume
        # Convert fractional coordinates to cartesians
        if len(self.cartesian_coordinates) == 0:
            if len(self.fractional_coordinates) == 0:
                print "Error no coordinates fraction or cartesian found"
                exit()
            for atom_frac in self.fractional_coordinates:
                atom_cart = cell.convertAbc2Xyz(atom_frac)
                self.cartesian_coordinates.append(atom_cart)
             # end for
        # end if
        self.unitCells[-1].fractionalCoordinates(self.fractional_coordinates)

    def _read_cellcontents(self,line):
        self.atom_types=[]
        self.masses=[]
        self.cartesian_coordinates=[]
        self.atomic_charges=[]
        self.ncores = 0
        self.nshells = 0
        for skip in range(0,5):
            line = self.fd.readline()
        for ion in range(0, self.nions):
            line = self.fd.readline()
            atom_type = line.split()[1]
            coreshell = line.split()[2]
            if coreshell == 's':
                self.nshells += 1
            else:
                self.ncores += 1
                atom_frac = [ float(f) for f in (line.split()[3:6]) ]
                q = float(line.split()[6])
                self.atom_types.append(atom_type)
                self.masses.append(self.mass_dictionary[atom_type])
                self.atomic_charges.append(q)
                self.cartesian_coordinates.append(atom_frac)
        #jk print "Number of cores: ", self.ncores
        #jk print "Masses ", self.masses
        self.nions = self.ncores
        if len(self.born_charges) == 0:
            for q in self.atomic_charges:
                a = np.zeros( (3,3) )
                a[0,0] = q
                a[1,1] = q
                a[2,2] = q
                self.born_charges.append(a)
            # end loop over charges
        # end if len()
        return

    def _read_cellcontentsf(self,line):
        self.atom_types=[]
        self.masses=[]
        self.fractional_coordinates=[]
        self.atomic_charges=[]
        self.ncores = 0
        self.nshells = 0
        for skip in range(0,5):
            line = self.fd.readline()
        for ion in range(0, self.nions):
            line = self.fd.readline()
            atom_type = line.split()[1]
            coreshell = line.split()[2]
            if coreshell == 's':
                self.nshells += 1
            else:
                self.ncores += 1
                atom_frac = [ float(f) for f in (line.split()[3:6]) ]
                self.atom_types.append(atom_type)
                self.masses.append(self.mass_dictionary[atom_type])
                self.fractional_coordinates.append(atom_frac)
        #jk print "Number of cores: ", self.ncores
        #jk print "Masses ", self.masses
        self.nions = self.ncores
        if len(self.born_charges) == 0:
            for q in self.atomic_charges:
                a = np.zeros( (3,3) )
                a[0,0] = q
                a[1,1] = q
                a[2,2] = q
                self.born_charges.append(a)
            # end loop over charges
        # end if len()
        return

    def _read_species(self,line):
        self.mass_dictionary = {}
        for skip in range(0,6):
            line = self.fd.readline()
        n = len(line.split())
        while n > 1 :
          species = line.split()[0]
          coreshell = line.split()[1]
          if coreshell == "Core":
              self.mass_dictionary[species] = (float(line.split()[3]))
          line = self.fd.readline()
          n = len(line.split())
        return


    def _read_born_charges(self,line):
        """Read the born charges from the gulp file.  
           Each column of the output refers to a given field direction
           Each row in the row refers the atomic displacement 
           So the numbers in the output are arrange  [ [ a1x a1y a1z ]  
                                                       [ a2x a2y a3z ] 
                                                       [ a3x a3y a3z ]]
           The output tensor needs them arranged     [ [ a1x a1y a1z ] [ a2x a2y a2z ] [ a3x a3y a3z ]]
           where 1,2,3 are the field directions and x, y, z are the atomic displacements"""
        #jk print "Born Charges no. of ions " ,self.nions
        self.born_charges = []
        for skip in range(0,5):
            line = self.fd.readline()
        for i in range(self.nions) :
          b = []
          b.append( [ float(f) for f in line.split()[3:6] ] )
          line = self.fd.readline()
          b.append( [ float(f) for f in line.split()[1:4] ] )
          line = self.fd.readline()
          b.append( [ float(f) for f in line.split()[1:4] ] )
          #jk B = np.array(b)
          #jk C = B.T
          #jk self.born_charges.append(B.tolist())
          self.born_charges.append(b)
          line = self.fd.readline()
          line = self.fd.readline()
        #jk print self.born_charges
        return

    def _read_optical_dielectric(self,line):
        for skip in range(0,5):
            line = self.fd.readline()
        # this is the zero frequency optical dielectric constant
        self.zerof_optical_dielectric = []
        self.zerof_optical_dielectric.append(  [ float(f) for f in line.split()[1:4] ] )
        line = self.fd.readline()
        self.zerof_optical_dielectric.append(  [ float(f) for f in line.split()[1:4] ] )
        line = self.fd.readline()
        self.zerof_optical_dielectric.append(  [ float(f) for f in line.split()[1:4] ] )
        return

    def _read_static_dielectric(self,line):
        for skip in range(0,5):
            line = self.fd.readline()
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric = []
        self.zerof_static_dielectric.append(  [ float(f) for f in line.split()[1:4] ] )
        line = self.fd.readline()
        self.zerof_static_dielectric.append(  [ float(f) for f in line.split()[1:4] ] )
        line = self.fd.readline()
        self.zerof_static_dielectric.append(  [ float(f) for f in line.split()[1:4] ] )
        return

    def _read_temperature(self,line):
        self.temperature = float(line.split()[4])
        return

    def _read_external_pressure(self,line):
        self.pressure = float(line.split()[4])
        return

    def _read_energies(self,line):
        #self.final_energy_without_entropy = unkown
        #line = self.fd.readline()
        #self.final_free_energy = float(line.split()[5])
        #jk line = self.fd.readline()
        #jk line = self.fd.readline()
        #jk self.final_0K_energy = float(line.split()[5])
        return
