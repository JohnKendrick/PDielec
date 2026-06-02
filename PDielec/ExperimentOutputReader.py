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
import logging
import math
import re

import numpy as np

from PDielec import DielectricFunction
from PDielec.Calculator import initialise_diagonal_tensor, normal_modes, oscillator_strengths
from PDielec.Constants import amu, wavenumber
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell

logger = logging.getLogger(__name__)



class ExperimentOutputReader(GenericOutputReader):
    """Read the contents of a directory containing Experiment input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    names : Various
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
        names : Various
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
        self.manage["raman_tensors"] = (re.compile("raman_tensors"),    self._read_raman_tensors)
        self.manage["normal_modes"]  = (re.compile("normal_modes"),     self._read_normal_modes)
        self.manage["born_charges"]  = (re.compile("born_charges"),     self._read_born_charges)
        self.manage["frequencies"]   = (re.compile("frequencies"),      self._read_frequencies)
        for f in self._outputfiles:
            self._read_output_file(f)
        self._calculate_oscillator_strengths_from_born_charges()
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
            self.CrystalPermittivity.set_epsilon_infinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.set_volume(self.volume)
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
            self.CrystalPermittivity.set_epsilon_infinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.set_volume(self.volume)
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
            self.CrystalPermittivity.set_epsilon_infinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.set_volume(self.volume)
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
            self.CrystalPermittivity.set_epsilon_infinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.set_volume(self.volume)
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
            self.CrystalPermittivity.set_epsilon_infinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.set_volume(self.volume)
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
            self.CrystalPermittivity.set_epsilon_infinity(self.zerof_optical_dielectric)
        if self.volume:
            self.CrystalPermittivity.set_volume(self.volume)
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
        This method reads the 'nfreq' (number of frequencies) from the initial given line. Then, for each frequency, it
        reads the corresponding line to extract the frequency and its oscillator strength. These values are stored in
        the object's 'frequencies' and 'oscillator_strengths' lists, respectively. Each oscillator strength is
        initialized as a diagonal tensor with identical values across the diagonal.

        The method assumes 'line' is in a specific format, where the first entry in a line after splitting is the
        frequency and the second entry is the oscillator strength. It also uses '_read_line' method, which is not
        defined in the snippet, to read each subsequent line for the frequencies and their oscillator strengths.
        Furthermore, 'initialise_diagonal_tensor' is used to convert the scalar strength values into tensor form, which
        is not detailed here.

        """        
        nfreq = int(line.split()[1])
        self.frequencies = []
        self.oscillator_strengths = []
        for _i in range(nfreq):
            line = self._read_line()
            parts = line.split("#", maxsplit=1)[0].split()
            self.frequencies.append(float(parts[0]))
            if len(parts) >= 4:
                sxx, syy, szz = float(parts[1]), float(parts[2]), float(parts[3])
            elif len(parts) == 2:
                sxx = syy = szz = float(parts[1])
            else:
                sxx = syy = szz = 0.0
            self.oscillator_strengths.append(initialise_diagonal_tensor( [sxx, syy, szz] ) )
        return

    def _read_born_charges(self, line):
        """Read Born effective charge tensors.

        Parameters
        ----------
        line : str
            The trigger line, expected format: ``born_charges N`` where
            ``N`` is the number of atom tensors to read.

        Notes
        -----
        Each atom tensor is read from three non-comment lines containing
        the 3x3 Born effective charge matrix in electron units.  An optional
        atom label line may precede each tensor; labels are ignored.

        Example::

            born_charges 2
            # Cd1
              2.172 -0.052 -0.077
              0.064  2.024 -0.004
              0.059  0.065  2.244
            # N1
             -2.713  0.415  0.133
              0.402 -3.058  0.114
              0.218  0.069 -2.376

        """
        n = int(line.split()[1])
        charges = []
        while len(charges) < n:
            first = self._read_line()
            parts = first.split()
            try:
                row0 = [float(x) for x in parts[:3]]
            except ValueError:
                row0 = [float(x) for x in self._read_line().split()[:3]]
            row1 = [float(x) for x in self._read_line().split()[:3]]
            row2 = [float(x) for x in self._read_line().split()[:3]]
            charges.append([row0, row1, row2])
        self.born_charges = charges
        return

    def _calculate_oscillator_strengths_from_born_charges(self):
        """Calculate oscillator strengths when Born charges and normal modes are present."""
        if len(self.born_charges) != self.nions:
            return
        if not isinstance(self.mass_weighted_normal_modes, np.ndarray):
            return
        if self.mass_weighted_normal_modes.shape != (3*self.nions, self.nions, 3):
            return
        if not np.any(self.mass_weighted_normal_modes):
            return
        masses = np.array(self.masses) * amu
        xyz_normal_modes = normal_modes(masses, self.mass_weighted_normal_modes)
        self.oscillator_strengths = oscillator_strengths(xyz_normal_modes, np.array(self.born_charges))
        self._calculate_hessian_from_frequencies_and_normal_modes()
        return

    def _calculate_hessian_from_frequencies_and_normal_modes(self):
        """Reconstruct a mass-weighted Hessian from read frequencies and normal modes."""
        if len(self.frequencies) != 3*self.nions:
            return
        if self.mass_weighted_normal_modes.shape != (3*self.nions, self.nions, 3):
            return
        ncoords = 3*self.nions
        ut = np.zeros((ncoords, ncoords))
        for imode, mode in enumerate(self.mass_weighted_normal_modes):
            ut[imode, :] = np.array(mode, dtype=float).reshape(ncoords)
        frequencies_au = np.array(self.frequencies, dtype=float) * wavenumber
        f2 = np.diag(np.sign(frequencies_au) * np.real(frequencies_au * frequencies_au))
        self.hessian = np.real(ut.T @ f2 @ ut)
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
        This function updates the object's `species`, `masses_per_type`, and `ion_type_index` attributes with the
        information read from the input line(s). Initially, it reads the total number of species from the first line.
        Then, for each species, it reads its name and mass, appends these to the `species` and `masses_per_type` lists,
        respectively, and updates the `ion_type_index` dictionary to map the species name to its type index. The
        `nspecies` attribute is updated with the current number of species after each iteration.

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
        self.volume = cell.get_volume(units="Angstrom")
        if self.CrystalPermittivity:
            self.CrystalPermittivity.set_volume(self.volume)
        return

    def _read_lattice_vectors(self, line):
        """Read and process lattice vectors from an input line, then updates internal structures with the new unit cell.

        This method is intended to parse a series of lines that describe lattice vectors, scale them appropriately, and
        update the object's list of unit cells as well as other relevant properties like volume and, optionally, crystal
        permittivity. The lattice vectors are assumed to be in Angstrom

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
        self.volume = cell.get_volume("Angstrom")
        if self.CrystalPermittivity:
            self.CrystalPermittivity.set_volume(self.volume)
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

        This method reads a string representing a line containing the number of ions followed by the fractional
        coordinates and species of each ion. It updates the instance variables associated with the ions' type, position,
        mass, and further initializes oscillator strengths and frequencies arrays if they haven't been initialized yet.

        Parameters
        ----------
        line : str
            The line from which the number of ions and their fractional coordinates will be extracted. The line format is expected to be 'num_ions species x y z' for each ion, where 'num_ions' is the total number of ions, 'species' is the species identifier, and 'x y z' are the fractional coordinates.

        Returns
        -------
        None

        Notes
        -----
        This method updates several instance variables: - `nions`: The total number of ions as determined from `line`. -
        `ions_per_type`: A list that tracks the number of ions per species. - `masses`: A list of masses for each ion,
        ordered accordingly. - `atom_type_list`: A list mapping each ion to its species index. - `unit_cells[-1]`: The
        coordinates and element names for ions in the most recent unit cell are set based on the processed line. -
        Oscillator strengths `oscillator_strengths` and frequencies `frequencies` arrays are initialized as zero arrays
        if they haven't been already.


        Exceptions related to file reading or value conversion within the method are implicitly assumed to be handled
        outside of its scope.

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
            self.CrystalPermittivity.set_epsilon_infinity(self.zerof_optical_dielectric)
        return

    def _read_raman_tensors(self, line):
        """Read Raman tensors for each normal mode.

        Parameters
        ----------
        line : str
            The trigger line, expected format: ``raman_tensors N [units]``
            where *N* is the total number of modes and the optional *units*
            keyword is ``castep`` (default), ``epsilon``, or ``skelton``.

        Returns
        -------
        None

        Notes
        -----
        Each mode occupies four consecutive non-comment lines:

        1. Frequency in cm-1.
        2. First row of the 3×3 Raman tensor.
        3. Second row.
        4. Third row.

        Modes must be listed in ascending-frequency order (acoustic / silent
        modes first, with zero tensors), consistent with the ordering used
        throughout PDielec.

        Two unit conventions are supported via the optional keyword on the
        header line:

        ``castep`` (default)
            Elements are CASTEP-style polarizability-volume tensors divided by
            ``sqrt(Vcell)`` and are multiplied by ``4π`` on read, matching
            :class:`PDielec.CastepOutputReader.CastepOutputReader`.
        ``epsilon``
            Elements are already in PDielec's ``R_epsilon = sqrt(V)dε/dQ``
            convention, with units ``(Å/amu)^{0.5}``, and are stored directly.
        ``skelton``
            Elements are Skelton-style polarizability-volume derivatives
            ``R_alpha = dα_vol/dQ`` in ``Å²·amu^{-0.5}``.  They are converted
            to ``R_epsilon`` by multiplying by ``4π/sqrt(Vcell)``.

        Frequencies are taken from a preceding ``frequencies`` block; this
        section contains only the tensor data.

        Example input (ZnO, 12 modes, CASTEP units)::

            raman_tensors  12
            # 3x3 tensor rows for each mode, in ascending-frequency order
              0.000   0.000   0.000
              0.000   0.000   0.000
              0.000   0.000   0.000
              0.0487  -0.0101   0.0041
             -0.0101   0.0487   0.0041
              0.0041   0.0041   0.1208
            ...

        """
        parts = line.split()
        n = int(parts[1])
        units = parts[2].lower() if len(parts) > 2 else "castep"
        tensors = []
        for _i in range(n):
            row0 = [float(x) for x in self._read_line().split()[:3]]
            row1 = [float(x) for x in self._read_line().split()[:3]]
            row2 = [float(x) for x in self._read_line().split()[:3]]
            tensors.append(np.array([row0, row1, row2], dtype=float))
        if units == "castep":
            factor = 4.0 * math.pi
            tensors = [t * factor for t in tensors]
        elif units == "skelton":
            factor = 4.0 * math.pi / math.sqrt(self.volume)
            tensors = [t * factor for t in tensors]
        elif units != "epsilon":
            logger.warning(f"Unknown experimental Raman tensor units '{units}', assuming R_epsilon")
        self.raman_tensors = tensors
        return

    def _read_normal_modes(self, line):
        """Read mass-weighted normal-mode eigenvectors for each mode.

        Parameters
        ----------
        line : str
            The trigger line, expected format: ``normal_modes N``
            where *N* is the total number of modes.

        Returns
        -------
        None

        Notes
        -----
        Frequencies are taken from a preceding ``frequencies`` block.
        Each mode occupies *nions* consecutive non-comment lines, each
        containing three floating-point values (the x, y, z components of
        the mass-weighted displacement for that atom).

        The eigenvectors should be normalised to unit length
        (``Σ |u|² = 1``) and ordered to match the ``frequencies`` block —
        acoustic / imaginary modes first, then optical in ascending
        frequency order.

        Example input (ZnO, 12 modes, 4 atoms)::

            normal_modes  12
            # mode 0 (acoustic)
             0.19709682  -0.22840953   0.00000001
             0.19709682  -0.22840953   0.00000001
             0.41779946  -0.48417512   0.00000003
             0.41779946  -0.48417512   0.00000003
            # mode 1 ...

        """
        n = int(line.split()[1])
        modes = []
        for _i in range(n):
            mode = []
            for _atom in range(self.nions):
                mode.append([float(x) for x in self._read_line().split()[:3]])
            modes.append(mode)
        self.mass_weighted_normal_modes = np.array(modes, dtype=float)
        return

    def calculate_mass_weighted_normal_modes(self):
        """Return mass-weighted normal modes, preserving any read from file.

        Parameters
        ----------
        None

        Returns
        -------
        np.array
            Mass-weighted normal modes with shape (3*nions, nions, 3).
            If a ``normal_modes`` block was present in the input file the
            stored values are returned unchanged.  Otherwise a zero array
            is returned (no phonon data available).

        """
        if (isinstance(self.mass_weighted_normal_modes, np.ndarray)
                and self.mass_weighted_normal_modes.shape == (3*self.nions, self.nions, 3)
                and np.any(self.mass_weighted_normal_modes)):
            return self.mass_weighted_normal_modes
        self.mass_weighted_normal_modes = np.zeros( (3*self.nions,self.nions,3) )
        return self.mass_weighted_normal_modes
