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
"""AbinitOutputReader module."""

import math
import re

import numpy as np

from PDielec.Constants import amu, angs2bohr, atomic_number_to_element, hartree2ev
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell


class AbinitOutputReader(GenericOutputReader):
    """Read the contents of a directory containing Abinit input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    filenames : list or str
        A list of filenames or a single filename string corresponding to the Abinit output files to be processed.

    Attributes
    ----------
    type : str
        A string descriptor for the type of files being dealt with, set to 'Abinit output files'.
    _acell : NoneType
        Initially set to None, this attribute is intended to store cell dimensions after being processed.
    _charges : NoneType
        Initially set to None, this attribute is intended to store charge information after being processed.

    Notes
    -----
    This method calls the `__init__` method of its superclass, `GenericOutputReader`, passing the `filenames` to it for any necessary higher-level initialization tasks. It then sets the `type` attribute to describe the types of files it will process and initializes `_acell` and `_charges` attributes to `None`, indicating that they are yet to be defined or processed.

    """

    def __init__(self, filenames):
        """Initialize the AbinitOutputReader with given filenames.

        This method initializes an instance of the class, setting up preliminary attributes and conditions for further operations on Abinit output files.

        Parameters
        ----------
        filenames : list or str
            A list of filenames or a single filename string corresponding to the Abinit output files to be processed.

        """        
        GenericOutputReader.__init__(self, filenames)
        self.type = "Abinit output files"
        self._acell = None
        self._charges = None
        return

    def _read_output_files(self):
        """Define the strings to be used for searching and Read the Abinit file names.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # Define the search keys to be looked for in the files
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage["dynamical"]  = (re.compile("  Dynamical matrix,"), self._read_dynamical)
        self.manage["bornCharges"]  = (re.compile("  Effective charges,"), self._read_born_charges)
        self.manage["epsilon"]  = (re.compile("  Dielectric tensor,"), self._read_epsilon)
        self.manage["masses"]   = (re.compile("              amu "), self._read_masses)
        self.manage["nions"]    = (re.compile("            natom "), self._read_natom)
        self.manage["rprim"]  = (re.compile("            rprim "), self._read_lattice_vectors)
        self.manage["xred"]  = (re.compile("             xred "), self._read_xred)
        self.manage["typat"]    = (re.compile("            typat "), self._read_typat)
        self.manage["ntypat"]   = (re.compile("           ntypat "), self._read_ntypat)
        self.manage["acell"]    = (re.compile("            acell "), self._read_acell)
        self.manage["nkpt"]    = (re.compile("             nkpt "), self._read_kpoints)
        self.manage["band"]    = (re.compile("            nband "), self._read_band)
        self.manage["band1"]    = (re.compile("            nband1 "), self._read_band)
        self.manage["occupancy"]    = (re.compile("              occ "), self._read_occupancy)
        self.manage["occupancy1"]    = (re.compile("              occ1 "), self._read_occupancy)
        self.manage["ecut"]    = (re.compile("^ *ecut "), self._read_energy_cutoff)
        self.manage["kptrlatt"]    = (re.compile("         kptrlatt "), self._read_kpoint_grid)
        self.manage["electrons"]    = (re.compile("  fully or partial"), self._read_electrons)
        self.manage["pressure"]    = (re.compile("-Cartesian.*GPa"), self._read_pressure)
        self.manage["znucl"]    = (re.compile("^  *znucl "), self._read_znucl)
        self.manage["totalenergy"]    = (re.compile("^  *Total energy "), self._read_total_energy)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_total_energy(self, line):
        """Read and store total energy values from a line.

        Parameters
        ----------
        line : str
            A string containing energy information, expected to have a numeric value at the 5th position (index 4) which represents energy in Hartrees.

        Returns
        -------
        None

        """        
        self.final_energy_without_entropy = float(line.split()[4]) * hartree2ev
        self.final_free_energy = float(line.split()[4]) * hartree2ev
        return

    def _read_znucl(self, line):
        """Parse and store atomic species information from a line input.

        This method extracts atomic species from a given line, converts them to 
        their element representation using their atomic numbers, and stores them 
        in a list. The total number of species processed is also updated.

        Parameters
        ----------
        line : str
            A string containing atomic numbers of elements, typically read from a file.

        Returns
        -------
        None

        Notes
        -----
        - Atomic numbers in the line are expected to be separated by space.
        - The atomic numbers are converted to integers with a small increment to 
          account for any floating-point representation issues.
        - This method updates the instance attributes `species` with the element 
          names capitalized, and `nspecies` with the count of unique species.

        See Also
        --------
        atomic_number_to_element : Dictionary mapping atomic numbers to element symbols.

        """        
        self.species = []
        for z in line.split()[1:]:
            iz = int(float(z)+0.001)
            self.species.append(atomic_number_to_element[iz].capitalize())
        self.nspecies = len(self.species)
        return

    def _read_band(self, line):
        """Read the number of bands from a line and update the band count.

        This private method parses a given line to extract the number of bands
        and updates the `nbands` attribute of the object accordingly.

        Parameters
        ----------
        line : str
            The line containing the band information, expected to have the number
            of bands as its second element when split by whitespace.

        Returns
        -------
        None

        """        
        self.nbands = int(line.split()[1])
        return

    def _read_occupancy(self, line):
        """Read and calculate the total occupancy from a line or lines in a file.

        This method reads occupancies from a given line, and if necessary, continues reading the next lines from the file until it gathers occupancies equal to the number of bands (`self.nbands`). It then calculates the sum of these occupancies and updates the total number of electrons.

        Parameters
        ----------
        line : str
            A string containing the occupancies for one or multiple bands. Each occupancy is expected to be separated by a space. The first item, if present, is ignored as it is not an occupancy value.

        Returns
        -------
        None

        Attributes Modified
        -------------------
        self.electrons : int
            Sets the `self.electrons` attribute to the integer value closest to the total sum of the occupancies read. The total is slightly adjusted before conversion to ensure correct rounding.

        Raises
        ------
        ValueError
            If the input line does not contain the expected format or types for occupancies.

        Notes
        -----
        - This method modifies the state of the object by setting `self.electrons` based on the read occupancies.
        - It implicitly relies on `self.nbands` to know how many occupancy values need to be read.

        """        
        occs = []
        occupancies = line.split()[1:]
        while len(occs) < self.nbands:
            occs+= [ float(f) for f in occupancies ]
            occupancies = self.file_descriptor.readline().split()
        sum = 0.0
        for f in occs:
            sum += f
        self.electrons = int(sum + 0.0001)
        return

    def _read_pressure(self, line):
        """Read and update the pressure attribute from a line of text.

        Parameters
        ----------
        line : str
            A line of text containing pressure data. The pressure value is expected to be the eighth element (index 7) in a space-separated list.

        Returns
        -------
        None

        Notes
        -----
        The function expects `line` to be a well-formatted string containing at least 8 space-separated values. The eighth value is converted to a float and stored in the `pressure` attribute of the class instance. This function does not return any value.

        """        
        self.pressure = float(line.split()[7])
        return

    def _read_electrons(self, line):
        """Read and set the number of electrons from a given line.

        This method sets the number of electrons for the object based on the string input. It expects the string to be formatted in a specific way where the seventh item (index 6, as indexing is zero-based) when split by whitespace, is convertible to float and represents the number of electrons.

        Parameters
        ----------
        line : str
            The input string containing electron information. It's expected that the seventh item in this string, when split by whitespace, can be converted to a float.

        Returns
        -------
        None

        Notes
        -----
        This method modifies the object's state by setting the 'electrons' attribute.

        """        
        self.electrons = float(line.split()[6])
        return

    def _read_kpoint_grid(self, line):
        """Parse and set the k-point grid from a given line of text.

        Parameters
        ----------
        line : str
            A string containing k-point grid information, expected to have integers
            at specific positions (1st, 5th, and 9th words).

        Returns
        -------
        None

        Notes
        -----
        This function is intended to be used internally within its class.
        It reads a line of text, extracts integers from specific positions,
        and sets them as the k-point grid dimensions for the instance.

        """        
        self.kpoint_grid = [ int(line.split()[1]), int(line.split()[5]), int(line.split()[9]) ]
        return

    def _read_kpoints(self, line):
        """Read and store the k-points from a line of text.

        This method processes a line of text to extract the number of k-points
        and stores it in the instance variable `kpoints`.

        Parameters
        ----------
        line : str
            A line of text that includes the number of k-points, expected to be
            in a format where the number of k-points follows a keyword, separated
            by spaces (e.g., "KPOINTS 100").

        Returns
        -------
        None

        Notes
        -----
        The method does not return any value but updates the instance variable
        `kpoints` with the integer value of the k-points extracted from the input line.

        """        
        self.kpoints = int(line.split()[1])
        return

    def _read_energy_cutoff(self, line):
        """Process and update the energy cutoff value from a given line.

        This function reads a line containing energy information, converts the energy value from Hartree to electron volts (eV), and updates the object's energy_cutoff attribute accordingly.

        Parameters
        ----------
        line : str
            A string containing the energy cutoff information. It is expected that the energy value is in Hartree units and is the second token in the string when split by whitespace.

        Returns
        -------
        None

        """        
        self.energy_cutoff = hartree2ev * float(line.split()[1])
        return

    def _read_acell(self, line):
        """Read and convert cell parameters from bohr to Ångström.

        Parameters
        ----------
        line : str
            A string containing the cell parameters in Ångström units separated by spaces. This string should start with an identifier followed by three numerical values corresponding to the cell parameters.

        Returns
        -------
        None
            This function does not return a value but updates the `_acell` attribute of the object.

        """        
        self._acell = [float(f)/angs2bohr for f in line.split()[1:4]]
        # Create a unit cell with diagonal elements of acell
        a = [ self._acell[0], 0.0, 0.0]
        b = [ 0.0, self._acell[1], 0.0]
        c = [ 0.0, 0.0, self._acell[2]]
        # Append the cell to the list of cells
        self.unit_cells.append(UnitCell(a,b,c,units="Angstrom"))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        return

    def _read_ntypat(self, line):
        """Read and set the number of species from a given line.

        Parameters
        ----------
        line : str
            A string containing the information of the number of species, 
            where the relevant data is expected to be in the second token 
            when splitting the string by whitespace.

        Returns
        -------
        None

        """        
        self.nspecies = int(line.split()[1])
        return

    def _read_typat(self, line):
        # typat occurs last in the list of data items we need from the output file
        """Process and read atom type patterns from a line of text.

        Parameters
        ----------
        line : str
            The line of text containing atom type patterns to be processed.

        Returns
        -------
        None

        Attributes Set
        --------------
        atom_type_list : list of int
            This list contains the atom type indices from the processed line, adjusted by -1.
        masses : list of optional[float]
            A list initialized with None values, set to the masses corresponding to the atom types
            read from the input line. The length of the list is equal to the number of ions (`nions`).
        ions_per_type : list of int
            Initialized as a list of zeros with length equal to the number of species (`nspecies`). Each
            entry represents the number of ions for each species, updated based on the atom types read.
        species_list : list of str, optional
            A list of species names corresponding to the atom types if the species attribute is set. This
            list is used to assign element names to the unit cells.

        """        
        self.atom_type_list = [int(i)-1 for i in line.split()[1:]]
        self.masses = [None for i in range(self.nions)]
        self.ions_per_type = [0 for i in range(self.nspecies)]
        for i, a in enumerate(self.atom_type_list):
            self.ions_per_type[a-1] += 1
            self.masses[i] = self.masses_per_type[a]
        if self.species:
            species_list = [ self.species[i] for i in self.atom_type_list ]
            self.unit_cells[-1].set_element_names(species_list)
        return

    def _read_epsilon(self, line):
        """Read zero frequency optical permittivity.

        Parse epsilon related data, and populates the `zerof_optical_dielectric` array
        based on the read data. 

        Parameters
        ----------
        line : str
            The current line being processed. This parameter is not directly used
            in the function but typically represents the context in which this method is called.

        Returns
        -------
        None

        This method does not return any value.

        See Also
        --------
        file_descriptor : An attribute that should be an open file object from which
                          the data is read.

        zerof_optical_dielectric : An attribute, typically a 2D list or an array, that this
                                   method updates based on the data read from the file.

        """        
        for _i in range(3):
            linea = self.file_descriptor.readline().split()
        nlines = 9
        for _i in range(nlines):
            linea = self.file_descriptor.readline().split()
            if not linea:
                linea = self.file_descriptor.readline().split()
            j = int(linea[0])
            k = int(linea[2])
            self.zerof_optical_dielectric[j-1][k-1] = float(linea[4])
        return

    def _read_natom(self, line):
        """Read and process the number of atoms from a given line.

        This method reads the number of atoms from a given line of text, sets the object's 'nions' attribute to this number, and initializes the '_charges' attribute to a zero array with shape `(nions, 3, 3)`. This method does not return any value.

        Parameters
        ----------
        line : str
            A line of text containing the number of atoms.

        Returns
        -------
        None

        See Also
        --------
        numpy.zeros : Used to initialize the '_charges' attribute.

        """        
        self.nions = int(line.split()[1])
        # We can only create this once we know the number of ions
        self._charges = np.zeros((self.nions, 3, 3))
        return

    def _read_masses(self, line):
        """Parse and save masses from a given line of text.

        Extracts and converts to float each mass value from a string,
        excluding the first token, and saves them into `masses_per_type`.

        Parameters
        ----------
        line : str
            The input string containing mass values starting from the second token.

        Returns
        -------
        None

        """        
        self.masses_per_type = [float(f) for f in line.split()[1:]]
        finished = False
        while not finished:
            linea = self.file_descriptor.readline().split()
            for f in linea:
                try:
                    self.masses_per_type.append(float(f))
                except ValueError:
                    finished = True
                    break
        return

    def _read_dynamical(self, line):
        # Read the dynamical matrix
        """Read and populate the Hessian matrix from a file.

        This method takes a line from a file, reads the dynamical matrix data from subsequent lines, and computes the Hessian matrix from these data. 

        Parameters
        ----------
        line : str
            The current line being read from the file. This parameter is not directly used in the function but indicates the starting point in the file from where to read the dynamical matrix data.

        Returns
        -------
        None

        Notes
        -----
        This method assumes the file format is structured in a specific way, where each line after the first four lines (which are skipped) contains the indices and value needed to populate the Hessian matrix.
        The Hessian matrix is used in molecular dynamics and vibrational analysis to describe the second-order interaction between atoms.

        See Also
        --------
        _dynamical_matrix : Method called with the computed Hessian matrix.

        """        
        nmodes = self.nions*3
        hessian = np.zeros((nmodes, nmodes))
        for _i in range(4):
            self.file_descriptor.readline()
        nlines = nmodes*nmodes
        for _i in range(nlines):
            linea = self.file_descriptor.readline().split()
            if not linea:
                linea = self.file_descriptor.readline().split()
            diri  = int(linea[0])
            atomi = int(linea[1])
            dirj  = int(linea[2])
            atomj = int(linea[3])
            ipos  = (atomi - 1)*3 + diri - 1
            jpos  = (atomj - 1)*3 + dirj - 1
            # store the massweighted matrix
            hessian[ipos][jpos] = float(linea[4])/(amu*math.sqrt(self.masses[atomi-1]*self.masses[atomj-1]))
        # symmetrise, project diagonalise and store frequencies and normal modes
        self._dynamical_matrix(hessian)
        return

    def _read_born_charges(self, line):
        """Read the born charges from the output file.

        Each row of the output refers to a given field direction.
        Each column in the row refers to the atomic displacement,
        so the output is arranged as follows: ::

            [[a1x, a1y, a1z],
             [a2x, a2y, a2z],
             [a3x, a3y, a3z]]

        where 1, 2, 3 are the field directions and x, y, z are the atomic displacements.

        Parameters
        ----------
        line : str
            Not used

        Returns
        -------
        numpy.ndarray
            A 2D array where each row represents a field direction and each column represents an atomic displacement in the x, y, and z directions.

        """
        for _i in range(5):
            self.file_descriptor.readline()
        #  The charges are calculated in two ways, we take the mean of the phonon and the field
        nlines = 9*self.nions
        for _i in range(nlines):
            linea = self.file_descriptor.readline().split()
            if not linea:
                linea = self.file_descriptor.readline().split()
            if int(linea[3]) > self.nions:
                ifield = int(linea[2])
                ixyz   = int(linea[0])
                iatom  = int(linea[1])
            else:
                ifield = int(linea[0])
                ixyz   = int(linea[2])
                iatom  = int(linea[3])
            self._charges[iatom-1][ifield-1][ixyz-1] += 0.5*float(linea[4])
        # Convert the charges
        self.born_charges = []
        for i in range(self.nions):
            atom = []
            for ifield in range(3):
                b = self._charges[i][ifield][:].tolist()
                atom.append(b)
            self.born_charges.append(atom)
        return

    def _read_xred(self, line):
        """Read the fractional coordinates from a file and update the last unit cell's fractional coordinates and element names.

        Parameters
        ----------
        line : str
            The initial line containing part of the fractional coordinates.

        Returns
        -------
        None

        Notes
        -----
        - This method assumes that the current context includes attributes for 'nions' (the number of ions),
          'file_descriptor' (a file object for reading), 'unit_cells' (a list of unit cell objects),
          and 'species' (a list of element species names) along with 'atom_type_list' (a list of indices for species).
        - It is further assumed that the 'unit_cells' list has at least one unit cell object that supports
          'set_fractional_coordinates' and 'set_element_names' methods for updating its state.
        - The method reads a line, splits it to extract fractional coordinates, and then reads additional lines
          corresponding to the remaining ions ('nions' - 1).
        - The extracted fractional coordinates are stored in a nested list, which is passed to the last unit cell's
          'set_fractional_coordinates' method.
        - If the 'species' attribute is provided along with a coherent 'atom_type_list', the method will also set element names for the 
          last unit cell using these mappings.

        No direct returns but updates the state of the last unit cell in the 'unit_cells' list by setting its fractional coordinates and optionally, element names.

        """        
        linea = line.split()[1:]
        fractional = []
        fractional.append( [ float(xyz) for xyz in linea ] )
        for _i in range(self.nions-1):
            linea = self.file_descriptor.readline().split()
            fractional.append( [ float(xyz) for xyz in linea ] )
        # end for i
        self.unit_cells[-1].set_fractional_coordinates(fractional)
        if self.species:
            species_list = [ self.species[i] for i in self.atom_type_list ]
            self.unit_cells[-1].set_element_names(species_list)
    # end def


    def _read_lattice_vectors(self, line):
        """Read lattice vectors from a line of a file and update the object state.

        Reads three lines starting from the given line from the file associated with the object's file descriptor, interprets them as lattice vectors, scales them according to the object's _acell attribute, and updates the object by appending a new `UnitCell` to its `unit_cells` list and updating the `ncells` and `volume` attributes.

        Parameters
        ----------
        line : str
            The line from which to start reading lattice vectors. This line should contain the first lattice vector.

        Returns
        -------
        None

        Notes
        -----
        This method modifies the internal state of the object by adding a new `UnitCell` to `unit_cells`, updating `ncells` to the new total number of unit cells, and updating `volume` to the volume of the most recently added `UnitCell`.

        This method assumes that `line` is formatted correctly, with lattice vector components being float-convertible strings separated by spaces. It also assumes that `self.file_descriptor` is an open file object positioned appropriately to read the subsequent lines for the b and c vectors.

        This method requires that `self._acell` is already set to an iterable of three scaling factors corresponding to each of the a, b, and c lattice vectors.

        See Also
        --------
        UnitCell : The class used to represent a unit cell, which takes vectors a, b, and c as parameters to its constructor.

        """        
        linea = line.split()
        avector = [float(linea[1]), float(linea[2]), float(linea[3])]
        linea = self.file_descriptor.readline().split()
        bvector = [float(linea[0]), float(linea[1]), float(linea[2])]
        linea = self.file_descriptor.readline().split()
        cvector = [float(linea[0]), float(linea[1]), float(linea[2])]
        avector = [f * self._acell[0] for f in avector]
        bvector = [f * self._acell[1] for f in bvector]
        cvector = [f * self._acell[2] for f in cvector]
        # Overwrite the unit cell created by the acell directive
        self.unit_cells[-1] = UnitCell(avector, bvector, cvector,units="Angstrom")
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        return
