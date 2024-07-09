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
"""Read the contents of a directory containing Experiment input and output files."""
import re

import numpy as np

from PDielec import DielectricFunction
from PDielec.Calculator import initialise_diagonal_tensor
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell


class ExperimentOutputReader(GenericOutputReader):
    """Read the contents of a directory containing Experiment input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    names : Various types
        Initial parameter used for construction, typically a string or list of identifiers that this instance will process or hold.

    Attributes
    ----------
    type : str
        A string indicating the type of output, set to 'Experimental output'.
    _ion_type_index : dict
        A dictionary to index ion types, initialized as empty.
    CrystalPermittivity : NoneType or Various types
        Initially None, intended to hold information about crystal permittivity.
    zerof_optical_dielectric : NoneType or Various types
        Initially None, intended to hold zero-frequency optical dielectric information.
    oscillator_strengths : NoneType or Various types
        Initially None, intended for storing oscillator strengths data.
    frequencies : NoneType or Various types
        Initially None, intended for storing frequencies data.

    """

    def __init__(self, names):
        """Initialize an instance of the class.

        Parameters
        ----------
        names : Various types
            Initial parameter used for construction, typically a string or list of identifiers that this instance will process or hold.

        Attributes
        ----------
        type : str
            A string indicating the type of output, set to 'Experimental output'.
        _ion_type_index : dict
            A dictionary to index ion types, initialized as empty.
        CrystalPermittivity : NoneType or Various types
            Initially None, intended to hold information about crystal permittivity.
        zerof_optical_dielectric : NoneType or Various types
            Initially None, intended to hold zero-frequency optical dielectric information.
        oscillator_strengths : NoneType or Various types
            Initially None, intended for storing oscillator strengths data.
        frequencies : NoneType or Various types
            Initially None, intended for storing frequencies data.

        """        
        GenericOutputReader.__init__(self, names)
        self.type                    = "Experimental output"
        self._ion_type_index = {}
        self.CrystalPermittivity = None
        self.zerof_optical_dielectric = None
        self.oscillator_strengths = None
        self.frequencies = None

    def _read_line(self):
        """Read and return the next non-comment, non-empty line from the file.

        This method reads lines from a file, skipping over any lines that are
        either empty, contain only whitespace, or start with the '#' character
        (indicating a comment).

        Parameters
        ----------
        None

        Returns
        -------
        str
            The next non-empty, non-comment line from the file. If the end of the
            file is reached, the returned string will be empty.

        Notes
        -----
        - Lines consisting solely of whitespace characters (spaces, tabs, newlines) are considered empty.

        """        
        line = self.file_descriptor.readline()
        while not line.strip() or line[0] == "#":
            line = self.file_descriptor.readline()
        return line

    def _read_output_files(self):
        """Define the search strings and read the Experiment files in the directory.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage["lattice"]       = (re.compile("lattice"),          self._read_lattice_vectors)
        self.manage["CPK_LATTICE"]   = (re.compile("&CELL"),            self._read_cpk_lattice_vectors)
        self.manage["cpk_lattice"]   = (re.compile("&cell"),            self._read_cpk_lattice_vectors)
        self.manage["species"]       = (re.compile("species"),          self._read_species)
        self.manage["fractional"]    = (re.compile("unitcell"),         self._read_fractional_coordinates)
        self.manage["CPK_CARTESIANS"]= (re.compile("&COORD"),           self._read_cpk_coords)
        self.manage["cpk_cartesians"]= (re.compile("&coord"),           self._read_cpk_coords)
        self.manage["static"]        = (re.compile("static"),           self._read_static_dielectric)
        self.manage["epsinf"]        = (re.compile("epsinf"),           self._read_static_dielectric)
        self.manage["fpsq"]          = (re.compile("fpsq"),             self._read_fpsq_model)
        self.manage["drude-lorentz"] = (re.compile("drude-lorentz"),    self._read_drude_lorentz_model)
        self.manage["constant"]      = (re.compile("constant"),         self._read_constant_model)
        # The order of these interpolates is important !
        self.manage["interpolate3"]  = (re.compile("interpolate_3"),    self._read_interpolate3_model)
        self.manage["interpolate6"]  = (re.compile("interpolate_6"),    self._read_interpolate6_model)
        self.manage["interpolate"]   = (re.compile("interpolate"),      self._read_interpolate1_model)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_constant_model(self, line):
        """Read in a full constant dielectric tensor 3 values on each line.

        Parameters
        ----------
        line : str
            Line to be read (not used)

        Returns
        -------
        None

        """
        od = []
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        # If we have complex input return a complex list, otherwise return a real list
        odc = np.array(od,dtype=complex)
        odi = np.absolute(np.imag(odc))
        sumi = np.sum(odi)
        if sumi < 1.0e-12:
            odc = np.real(odc)
        self.CrystalPermittivity = DielectricFunction.ConstantScalar(odc)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_interpolate1_model(self, line):
        """Read in a tabulated permittivity and use it for interpolation.

        Only a single permittivity is given (real & imaginary) which is taken to be isotropic

        Parameters
        ----------
        line : str
            Line to be read (not used)

        Returns
        -------
        None

        """
        line = self._read_line()
        line = line.lower()
        line = line.replace(","," ")
        split_line = line.split()
        omegas = []
        eps    = []
        while split_line[0] != "&end" and split_line[0] != "end":
            omega  = float(split_line[0])
            epsr   = float(split_line[1])
            epsi   = float(split_line[2])
            omegas.append(omega)
            eps.append(complex(epsr,epsi))
            line = self._read_line()
            line = line.lower()
            line = line.replace(","," ")
            split_line = line.split()
        # end for i
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction.Tabulate1(omegas,eps)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_interpolate3_model(self, line):
        """Read in a tabulated permittivity and use it for interpolation.

        Three permittivities are given from the diagonal components of the tensor

        Parameters
        ----------
        line : str
            Line to be read (not used)

        Returns
        -------
        None

        """
        line = self._read_line()
        line = line.lower()
        line = line.replace(","," ")
        split_line = line.split()
        omegas = []
        epsxx = []
        epsyy = []
        epszz = []
        while split_line[0] != "&end" and split_line[0] != "end":
            omega  = float(split_line[0])
            epsrxx = float(split_line[1])
            epsixx = float(split_line[2])
            epsryy = float(split_line[3])
            epsiyy = float(split_line[4])
            epsrzz = float(split_line[5])
            epsizz = float(split_line[6])
            omegas.append(omega)
            epsxx.append(complex(epsrxx,epsixx))
            epsyy.append(complex(epsryy,epsiyy))
            epszz.append(complex(epsrzz,epsizz))
            line = self._read_line()
            line = line.lower()
            line = line.replace(","," ")
            split_line = line.split()
        # end for i
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction.Tabulate3(omegas,epsxx,epsyy,epszz)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_interpolate6_model(self, line):
        """Read in a tabulated permittivity and use it for interpolation.

        This is a full 6 parameter permittivity matrix

        Parameters
        ----------
        line : str
            Line to be read

        Returns
        -------
        None

        """
        line = self._read_line()
        line = line.lower()
        line = line.replace(","," ")
        split_line = line.split()
        omegas = []
        epsxx = []
        epsyy = []
        epszz = []
        epsxy = []
        epsxz = []
        epsyz = []
        while split_line[0] != "&end" and split_line[0] != "end":
            omega  = float(split_line[0])
            epsrxx = float(split_line[1])
            epsixx = float(split_line[2])
            epsryy = float(split_line[3])
            epsiyy = float(split_line[4])
            epsrzz = float(split_line[5])
            epsizz = float(split_line[6])
            epsrxy = float(split_line[7])
            epsixy = float(split_line[8])
            epsrxz = float(split_line[9])
            epsixz = float(split_line[10])
            epsryz = float(split_line[11])
            epsiyz = float(split_line[12])
            omegas.append(omega)
            epsxx.append(complex(epsrxx,epsixx))
            epsyy.append(complex(epsryy,epsiyy))
            epszz.append(complex(epsrzz,epsizz))
            epsxy.append(complex(epsrxy,epsixy))
            epsxz.append(complex(epsrxz,epsixz))
            epsyz.append(complex(epsryz,epsiyz))
            line = self._read_line()
            line = line.lower()
            line = line.replace(","," ")
            split_line = line.split()
        # end for i
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction.Tablulate6(omegas,epsxx,epsyy,epszz,epsxy,epsxz,epsyz)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_drude_lorentz_model(self, line):
        """Read in the drude_lorentz model parameters.

        Parameters
        ----------
        line : str
            line to be read

        Returns
        -------
        None
 
        Notes
        -----
        A simple example of a 2 oscillator model is given below. All units are in cm-1.

        ```
        drude-lorentz
        xx  1
        413.7 1050.0  22.2
        yy  1
        413.7 1050.0  22.2
        zz  1
        413.7 1050.0  22.2
        ```

        """
        omegas_all = []
        strengths_all = []
        gammas_all = []
        for _diag in range(0,3):
            line = self._read_line().split()
            n = int(line[1])
            omegas = []
            strengths = []
            gammas = []
            for _i in range(n):
                line = self._read_line().split()
                omegas.append(float(line[0]))
                strengths.append(float(line[1]))
                gammas.append(float(line[2]))
            # end for i
            omegas_all.append(omegas)
            strengths_all.append(strengths)
            gammas_all.append(gammas)
        # end for diag
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction.DrudeLorentz(omegas_all,strengths_all,gammas_all)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_fpsq_model(self, line):
        """Read in the fpsq model parameters.

        There is a separator between each diagonal contribution to the permittivity.

        Parameters
        ----------
        line : str
            Line to be processed

        Returns
        -------
        None

        Notes
        -----
        The input format is expected to follow the example below for alpha quartz: ::

            epsxx 4 # 1.5K Data # A2 Phonon modes
              360.7     1.5        384.8    1.5
              497.9     3.1        553.6    2.8
              773.7     5.4        789.9    6.3
             1073.0     6.2       1238.7   12.4
            epsyy 4
              360.7     1.5        384.8    1.5
              497.9     3.1        553.6    2.8
              773.7     5.4        789.9    6.3
             1073.0     6.2       1238.7   12.4
            epszz 6 # E Phonon modes
              391.5     1.5        403.0    1.5
              454.0     2.6        510.5    1.5
              695.9     4.9        698.4    4.0
              797.2     4.8        810.0    4.3
             1063.7     6.1       1230.7    8.2
             1157.2     6.2       1154.9    6.1

        Each section starts with `epsxx`, `epsyy`, or `epszz` followed by the number of modes. 
        After that, data lines include mode frequencies and corresponding contributions.

        Returns
        -------
        None

        """
        omega_tos_all = []
        gamma_tos_all = []
        omega_los_all = []
        gamma_los_all = []
        for _diag in range(0,3):
            line = self._read_line().split()
            n = int(line[1])
            omega_tos = []
            gamma_tos = []
            omega_los = []
            gamma_los = []
            for _i in range(n):
                line = self._read_line().split()
                omega_tos.append(float(line[0]))
                gamma_tos.append(float(line[1]))
                omega_los.append(float(line[2]))
                gamma_los.append(float(line[3]))
            omega_tos_all.append(omega_tos)
            gamma_tos_all.append(gamma_tos)
            omega_los_all.append(omega_los)
            gamma_los_all.append(gamma_los)
            # end for i
        # end for diag
        # Create a dielectric function for use in calculations
        self.CrystalPermittivity = DielectricFunction.FPSQ(omega_tos_all,gamma_tos_all,omega_los_all,gamma_los_all)
        if self.zerof_optical_dielectric:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_frequencies(self, line):
        """Read frequencies and oscillator strengths from a line and append them to the object's attributes.

        Parameters
        ----------
        line : str
            A string containing the frequency count, followed by lines of frequencies and strengths.

        Returns
        -------
        None.

        Notes
        -----
        This method reads the 'nfreq' (number of frequencies) from the initial given line. Then, for each frequency, it reads
        the corresponding line to extract the frequency and its oscillator strength. These values are stored in the object's
        'frequencies' and 'oscillator_strengths' lists, respectively. Each oscillator strength is initialized as a diagonal
        tensor with identical values across the diagonal.

        The method assumes 'line' is in a specific format, where the first entry in a line after splitting is the frequency
        and the second entry is the oscillator strength. It also uses '_read_line' method, which is not defined in the snippet,
        to read each subsequent line for the frequencies and their oscillator strengths. Furthermore, 'initialise_diagonal_tensor'
        is used to convert the scalar strength values into tensor form, which is not detailed here.

        """        
        nfreq = int(line.split()[1])
        self.frequencies = []
        self.oscillator_strengths = []
        for _i in range(nfreq):
            line = self._read_line()
            self.frequencies.append(float(line.split()[0]))
            strength = float(line.split()[1])
            self.oscillator_strengths.append(initialise_diagonal_tensor( [strength, strength, strength] ) )
        return

    def _read_species(self, line):
        """Read species information from a line.

        Parameters
        ----------
        line : str
            A string containing the number of species followed by species' information in subsequent lines.

        Returns
        -------
        None

        Notes
        -----
        This function updates the object's `species`, `masses_per_type`, and `ion_type_index` attributes with the information read from the input line(s). Initially, it reads the total number of species from the first line. Then, for each species, it reads its name and mass, appends these to the `species` and `masses_per_type` lists, respectively, and updates the `ion_type_index` dictionary to map the species name to its type index. The `nspecies` attribute is updated with the current number of species after each iteration.

        """        
        nspecies = int(line.split()[1])
        self.species = []
        self.masses_per_type = []
        for _i in range(nspecies):
            line = self._read_line()
            species = line.split()[0]
            self.species.append(species)
            self.masses_per_type.append(float(line.split()[1]))
            self._ion_type_index[species] = self.nspecies
            self.nspecies = len(self.species)
        return

    def _read_cpk_lattice_vectors(self, line):
        """Process a cpk input file with cell information.

        The lattice vectors are assumed to be in Angstrom

        Parameters
        ----------
        line : str
            Line to be processed (not used)

        Returns
        -------
        None

        """
        line = line.lower()
        split_line = line.split()
        alpha = 90.0
        beta = 90.0
        gamma = 90.0
        while split_line[0] != "&end" and split_line[0] != "end":
            if  split_line[0] == "abc":
                a = float(split_line[1])
                b = float(split_line[2])
                c = float(split_line[3])
            elif split_line[0] == "alpha_beta_gamma":
                alpha = float(split_line[1])
                beta = float(split_line[2])
                gamma = float(split_line[3])
            # end if split-line
            line = self._read_line()
            line = line.lower()
            split_line = line.split()
        # end while
        cell = UnitCell(a,b,c,alpha,beta,gamma,units="Angstrom")
        self.unit_cells.append(cell)
        self.ncells = len(self.unit_cells)
        self.volume = cell.getVolume(units="Angstrom")
        if self.CrystalPermittivity:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_lattice_vectors(self, line):
        """Read and process lattice vectors from an input line, then updates internal structures with the new unit cell.

        This method is intended to parse a series of lines that describe lattice vectors, scale them appropriately, and update the object's list of unit cells as well as other relevant properties like volume and, optionally, crystal permittivity.
        The lattice vectors are assumed to be in Angstrom

        Parameters
        ----------
        line : str
            The first line from which to begin reading the lattice vector information.

        Returns
        -------
        None

        See Also
        --------
        UnitCell : A class or function used to represent and manipulate unit cell information.

        """        
        line = self._read_line()
        scalar = float(line.split()[0])
        line = self._read_line()
        avector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        line = self._read_line()
        bvector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        line = self._read_line()
        cvector = [scalar*float(line.split()[0]), scalar*float(line.split()[1]), scalar*float(line.split()[2])]
        cell = UnitCell(avector, bvector, cvector, units="Angstrom")
        self.unit_cells.append(cell)
        self.ncells = len(self.unit_cells)
        self.volume = cell.getVolume("Angstrom")
        if self.CrystalPermittivity:
            self.CrystalPermittivity.setVolume(self.volume)
        return

    def _read_cpk_coords(self, line):
        """Read CPK coordinates from a given line and update class attributes accordingly.

        The coordinates are assumed to be in Angstrom

        Parameters
        ----------
        line : str
            The initial line input to start reading CPK coordinates from.

        Returns
        -------
        None

        Notes
        -----
        - This method directly modifies several attributes of the class instance it 
          belongs to, including lists of species, ions, and masses, as well as data 
          structures related to unit cells.

        - This method utilizes `self.nspecies`, `self.nions`, `self.species`, `self.masses_per_type`,
          `self._ion_type_index`, `self.ions_per_type`, `self.atom_type_list`, `self.masses`, 
          `self.unit_cells`, `self.oscillator_strengths`, and `self.frequencies` as 
          relevant attributes that get updated.

        - Initial coordinates and species names for each unit cell are eventually encapsulated 
          in `self.unit_cells` with their respective setter methods.

        - The method is designed to handle CPK format data incrementally, assuming that the 
          beginning of a new entry is marked by a line that doesn't start with '&end' or 'end'.

        Raises
        ------
        - This function implicitly assumes correct formatting and content of the input data. 
          If the data is not correctly formatted or if `_read_line` fails to read further 
          lines as expected, unexpected behavior or errors could occur.

        See Also
        --------
        _read_line : Method used to read the next line for continued parsing.

        """        
        line = self._read_line()
        line = line.lower()
        line = line.replace(","," ")
        split_line = line.split()
        species_list = []
        ions = []
        self.nspecies = 0
        while split_line[0] != "&end" and split_line[0] != "end":
            self.nions += 1
            species = line.split()[0]
            if species not in self.species:
                self.species.append(species)
                self.masses_per_type.append(1.0)
                self._ion_type_index[species] = self.nspecies
                self.ions_per_type.append(1)
                self.nspecies = len(self.species)
            index = self._ion_type_index[species]
            self.atom_type_list.append(index)
            species_list.append(species)
            ions.append([float(f) for f in line.split()[1:4]])
            self.ions_per_type[index] += 1
            self.masses.append(self.masses_per_type[index])
            line = self._read_line()
            line = line.lower()
            line = line.replace(","," ")
            split_line = line.split()
        # end while
        self.unit_cells[-1].set_xyz_coordinates(ions,units="Angstrom")
        self.unit_cells[-1].set_element_names(species_list)
        if self.oscillator_strengths is None:
            self.oscillator_strengths = np.zeros( (3*self.nions,3,3) )
        if self.frequencies is None:
            self.frequencies = np.zeros( 3*self.nions )

    def _read_fractional_coordinates(self, line):
        """Read and process fractional coordinates from a string line.

        This method reads a string representing a line containing the number of ions followed by the fractional coordinates and species of each ion. It updates the instance variables associated with the ions' type, position, mass, and further initializes oscillator strengths and frequencies arrays if they haven't been initialized yet.

        Parameters
        ----------
        line : str
            The line from which the number of ions and their fractional coordinates will be extracted. The line format is expected to be 'num_ions species x y z' for each ion, where 'num_ions' is the total number of ions, 'species' is the species identifier, and 'x y z' are the fractional coordinates.

        Returns
        -------
        None

        Notes
        -----
        This method updates several instance variables:
        - `nions`: The total number of ions as determined from `line`.
        - `ions_per_type`: A list that tracks the number of ions per species.
        - `masses`: A list of masses for each ion, ordered accordingly.
        - `atom_type_list`: A list mapping each ion to its species index.
        - `unit_cells[-1]`: The coordinates and element names for ions in the most recent unit cell are set based on the processed line.
        - Oscillator strengths `oscillator_strengths` and frequencies `frequencies` arrays are initialized as zero arrays if they haven't been already.

        Exceptions related to file reading or value conversion within the method are implicitly assumed to be handled outside of its scope.

        """        
        ions = []
        self.nions = int(line.split()[1])
        self.ions_per_type = [ 0 for s in self.species ]
        self.masses = []
        species_list = []
        for _n in range(self.nions):
            line = self._read_line()
            species = line.split()[0]
            index = self._ion_type_index[species]
            self.atom_type_list.append(index)
            species_list.append(species)
            ions.append([float(f) for f in line.split()[1:4]])
            self.ions_per_type[index] += 1
            self.masses.append(self.masses_per_type[index])
        self.unit_cells[-1].set_fractional_coordinates(ions)
        self.unit_cells[-1].set_element_names(species_list)
        if self.oscillator_strengths is None:
            self.oscillator_strengths = np.zeros( (3*self.nions,3,3) )
        if self.frequencies is None:
            self.frequencies = np.zeros( 3*self.nions )
        return

    def _read_static_dielectric(self, line):
        """Read and process static dielectric data.

        Parameters
        ----------
        line : str
            Not used but overwritten

        Returns
        -------
        None

        Notes
        -----
        This method is designed to be used within a class that handles reading and processing
        simulation output data. It assumes access to a `_read_line` method which reads individual
        lines from the output, and potentially a `CrystalPermittivity` object for storing results.
        This method processes three lines of dielectric data, converts them into complex numbers,
        and checks if the imaginary part of the overall dielectric constant is negligible. If it is,
        the dielectric constant is treated as real. The result is stored in the class instance and,
        if applicable, updates the `CrystalPermittivity` object's epsilon infinity value.

        Raises
        ------
        ValueError
            If the input lines cannot be converted into complex numbers.

        """        
        # the is epsilon infinity
        od = []
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        line = self._read_line()
        od.append([complex(f) for f in line.split()[0:3]])
        # If we have complex input return a complex list, otherwise return a real list
        odc = np.array(od,dtype=complex)
        odi = np.absolute(np.imag(odc))
        sumi = np.sum(odi)
        if sumi < 1.0e-12:
            odc = np.real(odc)
        self.zerof_optical_dielectric = odc.tolist()
        if self.CrystalPermittivity:
            self.CrystalPermittivity.setEpsilonInfinity(self.zerof_optical_dielectric)
        return

    def calculate_mass_weighted_normal_modes(self):
        """Calculate the mass weight normal modes.

        Parameters
        ----------
        None

        Returns
        -------
        np.array
            Returns an array for the mass weighted normal modes which in this case is zero
            The array has a shape 3*nions, nions, 3

        """        
        self.mass_weighted_normal_modes = np.zeros( (3*self.nions,self.nions,3) )
        return self.mass_weighted_normal_modes
