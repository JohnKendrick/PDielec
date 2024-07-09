#!/usr/bin/python
#
# Copyright 2024 John Kendrick & Andrew Burnett
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
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""Module to read the contents of a Gulp output file."""

import os
import re
import sys

import numpy as np

from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell


class GulpOutputReader(GenericOutputReader):
    """Read the contents of a Gulp output file.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    This function provides a way to read and interpret the contents of a Gulp output file, 
    advising on a suitable set of command parameters that can be used for various calculations 
    such as phonon calculations, IR intensity calculations, printing eigen vectors of the dynamical 
    matrix, and printing cartesian coordinates. It emphasizes the importance of certain parameters 
    to ensure accurate and comprehensive results.

    Parameters
    ----------
    names : list
        A list containing the name(s) of the Gulp output file(s).

    Attributes
    ----------
    _gulpfile : str
        The name of the first Gulp output file in the `names` list.
    name : str
        The absolute path of the `_gulpfile`.
    type : str
        A string indicating the type of the reader object, set to 'Gulp output'.
    shells : int
        Number of shells, initialized to 0.
    ncores : int
        Number of cores, initialized to 0.
    _cartesian_coordinates : list
        List to store cartesian coordinates, initially empty.
    _fractional_coordinates : list
        List to store fractional coordinates, initially empty.
    atomic_charges : list
        List to store atomic charges, initially empty.
    _mass_dictionary : dict
        Dictionary to store mass information, initially empty.
    temperature : NoneType or float
        Temperature, initialized as None.
    elastic_constant_tensor : NoneType or ndarray
        Elastic constant tensor, initialized as None.
    nshells : NoneType or int
        Number of shells, initialized as None (may be updated later).
    nions_irreducible : NoneType or int
        Number of irreducible ions, initialized as None.
    _atom_types : NoneType or list
        List of atom types, initialized as None.

    Notes
    -----
    It is quite difficult to work out what is a suitable set of commands for Gulp. The following seems 
    to work OK, with emphasis on these important parameters:
    - phon: a phonon calculation.
    - intens: calculate the IR intensities.
    - eigen: print the eigen vectors of the dynamical matrix.
    - cart: print the cartesian coordinates (this is the only way of seeing what all the atom types are).
    - nosym: make sure the calculation is done in P1.
    Also includes options for various modifications and optimizations like qeq, molq, optimise, conp, qok, 
    nomodcoord, and prop.

    """

    def __init__(self, names):
        """Initialise the Gulp output reader object.

        Parameters
        ----------
        names : list
            A list containing the name(s) of the Gulp output file(s).

        Attributes
        ----------
        _gulpfile : str
            The name of the first Gulp output file in the `names` list.
        name : str
            The absolute path of the `_gulpfile`.
        type : str
            A string indicating the type of the reader object, set to 'Gulp output'.
        shells : int
            Number of shells, initialized to 0.
        ncores : int
            Number of cores, initialized to 0.
        _cartesian_coordinates : list
            List to store cartesian coordinates, initially empty.
        _fractional_coordinates : list
            List to store fractional coordinates, initially empty.
        atomic_charges : list
            List to store atomic charges, initially empty.
        _mass_dictionary : dict
            Dictionary to store mass information, initially empty.
        temperature : NoneType or float
            Temperature, initialized as None.
        elastic_constant_tensor : NoneType or ndarray
            Elastic constant tensor, initialized as None.
        nshells : NoneType or int
            Number of shells, initialized as None (may be updated later).
        nions_irreducible : NoneType or int
            Number of irreducible ions, initialized as None.
        _atom_types : NoneType or list
            List of atom types, initialized as None.

        Notes
        -----
        This constructor inherits from `GenericOutputReader` and is designed specifically for reading and parsing Gulp output files.

        """        
        GenericOutputReader.__init__(self, names)
        self._gulpfile               = names[0]
        self.name                    = os.path.abspath(self._gulpfile)
        self.type                    = "Gulp output"
        self.shells                  = 0
        self.ncores                  = 0
        self._cartesian_coordinates  = []
        self._fractional_coordinates  = []
        self.atomic_charges          = []
        self._mass_dictionary         = {}
        self.temperature             = None
        self.elastic_constant_tensor = None
        self.nshells                 = None
        self.nions_irreducible       = None
        self._atom_types              = None
        return

    def _read_output_files(self):
        """Initialise the phrases needed to search the .gout file and the routines to use when a match is found. Read the .gout file.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage["nions"]  = (re.compile(" *Total number atoms"), self._read_total_number_of_atoms)
        self.manage["nions_irreducible"]  = (re.compile(" *Number of irreducible atoms"), self._read_number_of_irreducible_atoms)
        self.manage["cellcontents"]  = (re.compile(" *Final Cartesian coor"), self._read_cellcontents)
        self.manage["cellcontentsf"]  = (re.compile(" *Final fractional coor"), self._read_cellcontentsf)
        self.manage["species"]        = (re.compile(" *Species output for all configurations"), self._read_species)
        self.manage["lattice"]          = (re.compile(" *Final Cartesian lattice"), self._read_lattice)
        self.manage["finalenergy"]   = (re.compile(" *Total lattice energy"), self._read_energies)
        self.manage["pressure"]      = (re.compile(" *Pressure of configuration"), self._read_external_pressure)
        self.manage["temperature"]      = (re.compile(" *Temperature of configuration"), self._read_temperature)
        self.manage["staticDielectric"]  = (re.compile(" *Static dielectric c"), self._read_static_dielectric)
        self.manage["opticalDielectric"]  = (re.compile(" *High frequency diel"), self._read_optical_dielectric)
        self.manage["bornCharges"]  = (re.compile(" *Born effective charge"), self._read_born_charges)
        self.manage["elasticConstants"]  = (re.compile(" *Elastic Constant Matrix"), self._read_elastic_constants)
        self.manage["frequencies"]  = (re.compile(" *Frequencies .cm-1."), self._read_frequencies)
        self._read_output_file(self._gulpfile)

    def _read_elastic_constants(self, line):
        """Read the elastic constants.

        Parameters
        ----------
        line : str (not used)
            The initial line

        Returns
        -------
        None

        """
        for _skip in range(0, 5):
            line = self.file_descriptor.readline()
        elastic = []
        for _skip in range(0, 6):
            elastic.append([float(f) for f in line.split()[1:]])
            line = self.file_descriptor.readline()
        self.elastic_constant_tensor = np.array(elastic)
        if self.debug:
            print("Elastic constant tensor")
            print(self.elastic_constant_tensor)
        return

    def _read_frequencies(self, line):
        """Read the frequencies.

        Parameters
        ----------
        line : str (not used)
            The initial line

        """
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        line = self.file_descriptor.readline()
        nmodes = 0
        nions = self.ncores
        while nmodes < 3*nions:
            line = self.file_descriptor.readline()
            if line.split()[0] == "Note:":
                line = self.file_descriptor.readline()
                line = self.file_descriptor.readline()
            freq = [float(f) for f in line.split()[1:]]
            ncolumns = len(freq)
            nmodes = nmodes + ncolumns
            self.frequencies += freq
            line = self.file_descriptor.readline()  # IR
            line = self.file_descriptor.readline()  # IR x
            line = self.file_descriptor.readline()  # IR y
            line = self.file_descriptor.readline()  # IR z
            line = self.file_descriptor.readline()  # Raman
            line = self.file_descriptor.readline()
            columns = []
            for _n in range(0, ncolumns):
                columns.append([])
            for _i in range(0, nions):
                line = self.file_descriptor.readline()
                modex = [float(f) for f in line.split()[2:]]
                line = self.file_descriptor.readline()
                modey = [float(f) for f in line.split()[2:]]
                line = self.file_descriptor.readline()
                modez = [float(f) for f in line.split()[2:]]
                for n, (x, y, z) in enumerate(zip(modex, modey, modez)):
                    columns[n].append([x, y, z])
                # end for x, y, z (columns)
            # end loop over atoms
            for mode in columns:
                self.mass_weighted_normal_modes.append(mode)
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()

    def _read_total_number_of_atoms(self, line):
        """Read the number of atoms.

        Parameters
        ----------
        line : str 
            The initial line

        Returns
        -------
        None

        """
        self.nions = int(line.split()[4])

    def _read_number_of_irreducible_atoms(self, line):
        """Read the number of irreducible atoms.

        Parameters
        ----------
        line : str 
            The initial line

        Returns
        -------
        None

        """
        self.nions_irreducible = int(line.split()[5])

    def _read_lattice(self, line):
        """Read the lattice vectors.

        Parameters
        ----------
        line : str (not used)
            The initial line

        Returns
        -------
        None

        """
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        x = float(line.split()[0])
        y = float(line.split()[1])
        z = float(line.split()[2])
        avector = [x, y, z]
        line = self.file_descriptor.readline()
        x = float(line.split()[0])
        y = float(line.split()[1])
        z = float(line.split()[2])
        bvector = [x, y, z]
        line = self.file_descriptor.readline()
        x = float(line.split()[0])
        y = float(line.split()[1])
        z = float(line.split()[2])
        cvector = [x, y, z]
        cell = UnitCell(avector, bvector, cvector, "Angstrom")
        # Convert fractional coordinates to cartesians
        if len(self._cartesian_coordinates) == 0:
            if len(self._fractional_coordinates) == 0:
                print("Error no coordinates fraction or cartesian found")
                sys.exit()
            for atom_frac in self._fractional_coordinates:
                atom_cart = cell.convert_abc_to_xyz(atom_frac)
                self.cartesian_coordinates.append(atom_cart)
            # end for
        # end if
        cell.set_fractional_coordinates(self._fractional_coordinates)
        cell.set_element_names(self._atom_types)
        self.unit_cells.append(cell)
        self.volume = cell.getVolume("Angstrom")
        self.ncells = len(self.unit_cells)

    def _read_cellcontents(self, line):
        """Read the cell contents in xyz space.

        Parameters
        ----------
        line : str (not used)
            The initial line

        Returns
        -------
        None

        """
        self._atom_types = []
        self.masses = []
        self._cartesian_coordinates = []
        self.atomic_charges = []
        self.ncores = 0
        self.nshells = 0
        self.atom_type_list = []
        self.ions_per_type = [ 0 for i in range(self.nspecies) ]
        for _skip in range(0, 5):
            line = self.file_descriptor.readline()
        for _ion in range(0, self.nions):
            line = self.file_descriptor.readline()
            atom_type = line.split()[1].capitalize()
            coreshell = line.split()[2]
            if coreshell == "s":
                self.nshells += 1
            else:
                self.ncores += 1
                atom_frac = [float(f) for f in line.split()[3:6]]
                q = float(line.split()[6])
                self._atom_types.append(atom_type)
                self.masses.append(self._mass_dictionary[atom_type])
                species_index=self.species.index(atom_type)
                self.atom_type_list.append(species_index)
                self.ions_per_type[species_index] += 1
                self.atomic_charges.append(q)
                self._cartesian_coordinates.append(atom_frac)
        self.nions = self.ncores
        if len(self.born_charges) == 0:
            for q in self.atomic_charges:
                a = np.zeros((3, 3))
                a[0, 0] = q
                a[1, 1] = q
                a[2, 2] = q
                self.born_charges.append(a)
            # end loop over charges
        # end if len()
        return

    def _read_cellcontentsf(self, line):
        """Read the cell contents in fractional space.

        Parameters
        ----------
        line : str (not used)
            The initial line

        Returns
        -------
        None

        """
        self._atom_types = []
        self.masses = []
        self._fractional_coordinates = []
        self.atomic_charges = []
        self.ncores = 0
        self.nshells = 0
        self.atom_type_list = []
        for _skip in range(0, 5):
            line = self.file_descriptor.readline()
        for _ion in range(0, self.nions):
            line = self.file_descriptor.readline()
            atom_type = line.split()[1].capitalize()
            coreshell = line.split()[2]
            if coreshell == "s":
                self.nshells += 1
            else:
                self.ncores += 1
                atom_frac = [float(f) for f in line.split()[3:6]]
                self._atom_types.append(atom_type)
                self.masses.append(self._mass_dictionary[atom_type])
                self.atom_type_list.append(self.species.index(atom_type))
                self._fractional_coordinates.append(atom_frac)
        self.nions = self.ncores
        if len(self.born_charges) == 0:
            for q in self.atomic_charges:
                a = np.zeros((3, 3))
                a[0, 0] = q
                a[1, 1] = q
                a[2, 2] = q
                self.born_charges.append(a)
            # end loop over charges
        # end if len()
        return

    def _read_species(self, line):
        """Read species information.

        Parameters
        ----------
        line : str (not used)
            The initial line

        Returns
        -------
        None

        """
        self.species = []
        self.mass_per_type = []
        self._mass_dictionary = {}
        for _skip in range(0, 6):
            line = self.file_descriptor.readline()
        n = len(line.split())
        while n > 1:
            species = line.split()[0].capitalize()
            coreshell = line.split()[1]
            if coreshell == "Core":
                self._mass_dictionary[species] = (float(line.split()[3]))
                if species not in self.species:
                    self.species.append(species)
                    self.masses_per_type.append(self._mass_dictionary[species])
            line = self.file_descriptor.readline()
            n = len(line.split())
        self.nspecies = len(self.species)
        return

    def _read_born_charges(self, line):
        """Read the born charges from the gulp file.

        This function reads the born charges as arranged in the gulp file. Each column 
        in the output fetched refers to a given field direction, and each row refers 
        to the atomic displacement. The function rearranges the output to match the 
        required tensor format where 1, 2, and 3 are the field directions, and x, y, 
        and z are the atomic displacements.

        Parameters
        ----------
        line : str
            Ignored

        Returns
        -------
        numpy.ndarray
            The born charges arranged in the tensor format as specified. The format
            of the tensor is a square matrix where each row is [ax, ay, az] for atomic
            displacements under a particular field direction.

        Examples
        --------
        Consider the output format [[a1x a1y a1z]
                                    [a2x a2y a2z]
                                    [a3x a3y a3z]]

        The function rearranges it to the required format:
        [[a1x a1y a1z] [a2x a2y a2z] [a3x a3y a3z]]

        Notes
        -----
        The original output from the gulp file is read line by line, with each line representing
        the atomic displacement in each of the field directions. The reformatted output aligns
        these displacements into a more conventional tensor format for further processing or
        analysis.

        """
        self.born_charges = []
        for _skip in range(0, 5):
            line = self.file_descriptor.readline()
        for _i in range(self.nions):
            b = []
            b.append([float(f) for f in line.split()[3:6]])
            line = self.file_descriptor.readline()
            b.append([float(f) for f in line.split()[1:4]])
            line = self.file_descriptor.readline()
            b.append([float(f) for f in line.split()[1:4]])
            # jk B = np.array(b)
            # jk C = B.T
            # jk self.born_charges.append(B.tolist())
            self.born_charges.append(b)
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
        return

    def _read_optical_dielectric(self, line):
        """Read optical dielectric constant.

        Parameters
        ----------
        line : str (not used)
            The initial line

        """
        for _skip in range(0, 5):
            line = self.file_descriptor.readline()
        # this is the zero frequency optical dielectric constant
        self.zerof_optical_dielectric = []
        self.zerof_optical_dielectric.append([float(f) for f in line.split()[1:4]])
        line = self.file_descriptor.readline()
        self.zerof_optical_dielectric.append([float(f) for f in line.split()[1:4]])
        line = self.file_descriptor.readline()
        self.zerof_optical_dielectric.append([float(f) for f in line.split()[1:4]])
        return

    def _read_static_dielectric(self, line):
        """Read static dielectric constant.

        Parameters
        ----------
        line : str (not used)
            The initial line

        Returns
        -------
        None

        """
        for _skip in range(0, 5):
            line = self.file_descriptor.readline()
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric = []
        self.zerof_static_dielectric.append([float(f) for f in line.split()[1:4]])
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric.append([float(f) for f in line.split()[1:4]])
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric.append([float(f) for f in line.split()[1:4]])
        return

    def _read_temperature(self, line):
        """Read temperature.

        Parameters
        ----------
        line : str
            The initial line

        Returns
        -------
        None

        """
        self.temperature = float(line.split()[4])
        return

    def _read_external_pressure(self, line):
        """Read external pressure.

        Parameters
        ----------
        line : str
            The initial line

        Returns
        -------
        None

        """
        self.pressure = float(line.split()[4])
        return

    def _read_energies(self, line):
        """Read energies.

        Parameters
        ----------
        line : str
            The initial line

        Returns
        -------
        None

        """
        #line = self.file_descriptor.readline()
        self.final_free_energy = float(line.split()[4])
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        #self.final_0K_energy = float(line.split()[5])
        self.final_energy_without_entropy = self.final_free_energy
        return
