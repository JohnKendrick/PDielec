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
"""A module to read contents of a directory containing Crystal input and output files."""
import math
import os
import re

import numpy as np

from PDielec.Constants import amu, hartree2ev
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.IO import pdielec_io
from PDielec.UnitCell import UnitCell


class CrystalOutputReader(GenericOutputReader):
    """Read contents of a directory containing Crystal input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    filenames : iterable
        An iterable object (e.g., list or tuple) containing the file names to be processed.

    Notes
    -----
    This function calls the initializer of the parent class 'GenericOutputReader' with the given filenames. It initializes several class attributes including the type of output (set to 'Crystal output'), the method for hessian symmetrisation (set to 'crystal'), and an empty list for fractional coordinates.

    """

    def __init__(self, filenames):
        """Initialize an instance of the class which handles crystal output.

        Parameters
        ----------
        filenames : iterable
            An iterable object (e.g., list or tuple) containing the file names to be processed.

        Notes
        -----
        This function calls the initializer of the parent class 'GenericOutputReader' with the given filenames. It initializes several class attributes including the type of output (set to 'Crystal output'), the method for hessian symmetrisation (set to 'crystal'), and an empty list for fractional coordinates.

        """        
        GenericOutputReader.__init__(self, filenames)
        self.type                    = "Crystal output"
        self.hessian_symmetrisation  = "crystal"
        self._fractional_coordinates = []
        return

    def _read_output_files(self):
        """Read the Crystal files in the directory.

        Define the strings used to search for the required information and then process the file

        Returns
        -------
        None

        """
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage["masses"]   = (re.compile(" ATOMS ISOTOPIC MASS"), self._read_masses)
        self.manage["lattice"]  = (re.compile(" DIRECT LATTICE VECTORS CARTESIAN COMPONENTS"), self._read_lattice_vectors)
        self.manage["fractional"]  = (re.compile(" ATOMS IN THE ASYMMETRIC UNIT"), self._read_fractional_coordinates)
        self.manage["bornCharges"]  = (re.compile(" ATOMIC BORN CHARGE TENSOR"), self._read_born_charges)
        self.manage["eigenvectors"]  = (re.compile(" NORMAL MODES NORMALIZ"), self._read_eigenvectors)
        self.manage["staticIonic"]  = (re.compile(" SUM TENSOR OF THE VIBRATIONAL CONTRIBUTIONS TO "), self._read_ionic_dielectric)
        self.manage["epsilon"]  = (re.compile(" SUSCEPTIBILITY "), self._read_epsilon)
        self.manage["kpoints"]  = (re.compile(" SHRINK\\. FACT\\.\\("), self._read_kpoints)
        self.manage["electrons"]  = (re.compile(" N\\. OF ELECTRONS"), self._read_electrons)
        self.manage["energy"]  = (re.compile(" TOTAL ENERGY\\(DFT\\)"), self._read_energy)
        self.manage["energy2"]  = (re.compile(" TOTAL ENERGY + DISP"), self._read_energy2)
        self.manage["energy3"]  = (re.compile(" *CENTRAL POINT"), self._read_energy3)
        for f in self._outputfiles:
            self._read_output_file(f)
        return

    def _read_energy(self, line):
        """Read and set the free energy and energy without entropy from a line.

        Parameters
        ----------
        line : str
            A string containing the energy information, expected to be in a format
            where the fourth element (index 3) after splitting contains the energy value
            in Hartrees.

        Returns
        -------
        None

        Notes
        -----
        This method modifies the instance attributes `final_free_energy` and
        `final_energy_without_entropy` by converting the energy value from Hartrees to
        electronvolts (eV) using a conversion factor stored in `hartree2ev`.

        There are no return values. Modifications are made directly to the instance attributes.

        """        
        self.final_free_energy = hartree2ev*float(line.split()[3])
        self.final_energy_without_entropy = hartree2ev*float(line.split()[3])

    def _read_energy2(self, line):
        """Read and set the final free energy and final energy without entropy from a line of text.

        This method reads a line of text that contains energy values and updates the
        instance's final free energy and final energy without entropy attributes.
        The energies are assumed to be in Hartree units and are converted to electronvolts (eV).

        Parameters
        ----------
        line : str
            A string containing the desired energy values. The sixth value in the line
            (when split by whitespace) is assumed to be the energy in Hartree units.

        Returns
        -------
        None

        """        
        self.final_free_energy = hartree2ev*float(line.split()[5])
        self.final_energy_without_entropy = hartree2ev*float(line.split()[5])

    def _read_energy3(self, line):
        """Parse energy information from a given line of text and update object properties.

        Parameters
        ----------
        line : str
            A line of text containing energy information where the 3rd element (index 2)
            after splitting by whitespace contains the energy value in Hartree units.

        Returns
        -------
        None

        """        
        self.final_free_energy = hartree2ev*float(line.split()[2])
        self.final_energy_without_entropy = hartree2ev*float(line.split()[2])

    def _read_electrons(self, line):
        """Parse the number of electrons from a string and assign it to the instance variable.

        Parameters
        ----------
        line : str
            The number of electrons is the sixth element in a space-separated list.

        Returns
        -------
        None

        """        
        self.electrons = int(line.split()[5])

    def _read_kpoints(self, line):
        """Parse a line of text to extract k-point grid dimensions and the total number of k-points.

        The function sets the object's kpoint_grid attribute to a list of three integers representing the grid dimensions, and the kpoints attribute to an integer representing the total number of k-points.

        Parameters
        ----------
        line : str
            A string from which k-point grid dimensions and the total number of k-points are extracted. The line is expected to have specific data at predefined positions: the k-point grid dimensions at indices 2, 3, and 4, and the total number of k-points at index 12, according to the split by whitespace.

        Returns
        -------
        None

        """        
        self.kpoint_grid = [ int(line.split()[2]), int(line.split()[3]), int(line.split()[4]) ]
        self.kpoints = int(line.split()[12])

    def _read_epsilon(self, line):
        """Read epsilon values from a formatted file and update the optical dielectric matrix.

        This internal method reads lines from a previously opened file (using the file_descriptor attribute of the object), interpreting data that specifies the components of the optical dielectric tensor. The method updates the object's `zerof_optical_dielectric` attribute with a 3x3 matrix representing the read optical dielectric values. Optical dielectric components are identified by their 'XX', 'YY', 'ZZ', 'XY', 'XZ', and 'YZ' labels in the file. The method expects the values in a specific format and parses them accordingly.


        Parameters
        ----------
        line : str
            The initial line passed to the method. Currently, this parameter is not used within the method.

        Returns
        -------
        None

        Notes
        -----
        - The method does not return any value but instead updates the `zerof_optical_dielectric` attribute of the object.
        - This method assumes that the file being read is opened and managed elsewhere in the class and that the `file_descriptor` attribute contains a valid file object.
        - The method reads lines sequentially from the `file_descriptor` and expects the file to be formatted in a specific way, with optical dielectric components starting from a certain line.
        - Each line read corresponds to a different component of the optical dielectric tensor, and the method parses these lines to identify the component and its value.
        - After reading and storing the components' values, the method also updates the symmetric components of the tensor (i.e., sets optical_dielectric[1][0], optical_dielectric[2][0], and optical_dielectric[2][1] to their corresponding symmetric values).

        """        
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        optical_dielectric = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        linea = self.file_descriptor.readline().split()
        while len(linea) > 0:
            component = linea[0]
            # Problems with Crystal17 file format
            if component == "0.000":
                component = linea[1]
                epsilon = float(linea[4])
            else:
                epsilon = float(linea[3])
            if component == "XX" : 
                optical_dielectric[0][0] = epsilon
            elif component == "YY" : 
                optical_dielectric[1][1] = epsilon
            elif component == "ZZ" : 
                optical_dielectric[2][2] = epsilon
            elif component == "XY" : 
                optical_dielectric[0][1] = epsilon
            elif component == "XZ" : 
                optical_dielectric[0][2] = epsilon
            elif component == "YZ" : 
                optical_dielectric[1][2] = epsilon
            # read the next line
            linea = self.file_descriptor.readline().split()
        # symmetrize
        optical_dielectric[1][0] = optical_dielectric[0][1]
        optical_dielectric[2][0] = optical_dielectric[0][2]
        optical_dielectric[2][1] = optical_dielectric[1][2]
        self.zerof_optical_dielectric = optical_dielectric
        return

    def _read_masses(self, line):
        """Read and parse the masses from a file descriptor line by line.

        This private method reads mass information from consecutive lines of a
        file, updating a list of masses and a list of masses per species type.

        Notes
        -----
        The function uses a file descriptor (self.file_descriptor) that should
        be opened before calling this function. It reads the masses from the
        file, expecting the mass values to be placed at specific positions in
        the lines (every third value starting from the third value on a new
        line after skipping one line). The method updates the instance
        variables `masses` (a list of all masses) and `masses_per_type` (a list
        of masses indexed by species type).

        The function assumes that the lines are formatted in a specific way and
        that `self.nions`, `self.nspecies`, and `self.atom_type_list` are set
        correctly before calling this method.

        Parameters
        ----------
        line : str
            The current line from which to start reading. Note that this
            parameter is not used in the function body and the actual reading
            starts from the next line of the file descriptor.

        Returns
        -------
        None

        """        
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        n = 2
        self.masses_per_type = [ 0 for i in range(self.nspecies) ]
        for i in range(self.nions):
            if n > 11:
                line = self.file_descriptor.readline()
                n = 2
            mass = float(line.split()[n])
            self.masses.append(mass)
            self.masses_per_type[self.atom_type_list[i]] = mass
            n = n + 3
        return

    def _read_eigenvectors(self, line):
        """Read eigenvectors from appropriate file.

        This function attempts to read eigenvectors for a given molecular system. It first checks for
        the file 'HESSFREQ.DAT' in the specified `open_directory`. If the file exists, it reads the
        eigenvectors from this file using the method `_read_hessfreq_dat`. If the file does not exist,
        it reads the eigenvectors from output data indicated by 'line' using the method
        `_read_output_eigenvectors`.

        Parameters
        ----------
        line : str
            The line or identifier used to read the eigenvectors from output data, if 'HESSFREQ.DAT'
            is not available.

        Returns
        -------
        None

        """        
        if os.path.isfile(self.open_directory+"/HESSFREQ.DAT"):
            # print("Reading Hessian from HESSFREQ.DAT",file=sys.stderr)
            self._read_hessfreq_dat(self.open_directory+"/HESSFREQ.DAT")
        else:
            # print("Reading Normal Modes from output file",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            # print("The numerical precision of the input is limited",file=sys.stderr)
            # print("If possible provide a HESSFREQ.DAT from the calculation",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            self._read_output_eigenvectors(line)
        return

    def _read_hessfreq_dat(self,filename):
        """Read a Hessian frequency data from a given file and update the dynamical matrix.

        Parameters
        ----------
        filename : str
            The path to the file containing Hessian frequency data.

        Notes
        -----
        This function reads a Hessian frequency data file, processes the data by mass weighting the Hessian matrix, and then updates the dynamical matrix of the system. 

        See Also
        --------
        _dynamical_matrix : Updates the dynamical matrix of the system using the processed Hessian matrix.

        """        
        fd2 = pdielec_io(filename, "r")
        nmodes = self.nions*3
        # Create a mass weighting vector
        n = 0
        massweight = np.zeros(nmodes)
        for a in range(self.nions):
            for _j in range(3):
                massweight[n] = 1.0 / math.sqrt(self.masses[a]*amu)
                n = n + 1
        # We read the hessian and store the mass weighted matrix
        hessian = np.zeros((nmodes, nmodes))
        line = fd2.readline()
        pos = 0
        for i in range(nmodes):
            for j in range(nmodes):
                hessian[i, j] = massweight[i]*massweight[j]*float(line.split()[pos])
                pos = pos + 1
                if pos >= 4:
                    line = fd2.readline()
                    pos = 0
                # end if pos
            # end for j
        # end for i
        fd2.close()
        # symmetrise, project, diagonalise and store the frequencies and normal modes
        self._dynamical_matrix(hessian)
        return

    def _read_output_eigenvectors(self, line):
        """Read and process output eigenvectors from a file descriptor.

        This method reads eigenvector data for each frequency from a file, processes this data by scaling it with the square root of the atom's masses to generate mass-weighted normal modes, and normalizes these modes.

        Parameters
        ----------
        line : str
            The current line in the file from which to start reading eigenvectors. This parameter is not explicitly used in the processing inside the function but hints at expecting a certain positioning within the file content.

        Returns
        -------
        None

        Notes
        -----
        - This method directly modifies the `frequencies` and `mass_weighted_normal_modes` attributes of the class instance.
        - Atoms' masses are used to compute mass-weighted normal modes
        - Normalization of each mass-weighted mode is performed at the end of the method.

        """        
        self.file_descriptor.readline()
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        n = 6
        nmodes = 3 * self.nions
        self.frequencies = []
        # Loop over all the modes
        for _j in range(nmodes):
            # They come in columns of 6 so after 6 modes we read the frequencies and initialse the 6 new modes
            n = n + 1
            if n >= 6:
                n = 0
                linef = self.file_descriptor.readline().split()[1:]
                linef = [float(f) for f in linef]
                self.file_descriptor.readline()
                atoms = []
                for f in linef:
                    self.frequencies.append(f)
                    atoms.append([])
                # Read through the XYZ components for each atom and store the mode in atoms
                for i in range(self.nions):
                    mass = np.sqrt(self.masses[i])
                    linex = self.file_descriptor.readline().split()[4:]
                    liney = self.file_descriptor.readline().split()[1:]
                    linez = self.file_descriptor.readline().split()[1:]
                    for a, x, y, z in zip(atoms, linex, liney, linez):
                        x = float(x) * mass
                        y = float(y) * mass
                        z = float(z) * mass
                        a.append([x, y, z])
                for a in atoms:
                    self.mass_weighted_normal_modes.append(a)
                self.file_descriptor.readline()
            # end of if n >= 6
        for i, mode in enumerate(self.mass_weighted_normal_modes):
            marray = np.array(mode)
            sumxyz = 0.0
            for atoms in marray:
                for xyz in atoms:
                    sumxyz = sumxyz + xyz*xyz
            marray = marray / np.sqrt(sumxyz)
            self.mass_weighted_normal_modes[i] = marray.tolist()
        return

    def _read_born_charges(self, line):
        """Read the born charges from the output file.

        Parameters
        ----------
        line : str
            The current line

        Returns
        -------
        numpy.ndarray
            A 3x3 array where each row of the output refers to a given field direction and each column in the row refers to the atomic displacement. The output is arranged as: ::
    
                [[a1x a1y a1z]
                 [a2x a2y a2z]
                 [a3x a3y a3z]]
    
            where 1, 2, 3 represent the field directions and x, y, z are the atomic displacements.

        """
        # Extract directory containing the output file
        if os.path.isfile(self.open_directory+"/BORN.DAT"):
            # print("Reading Born charge tensor from BORN.DAT",file=sys.stderr)
            self._read_born_charges_from_born_dat(self.open_directory+"/BORN.DAT")
        else:
            # print("Reading Born Charge Tensor from output file",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            # print("The numerical precision of the input is limited",file=sys.stderr)
            # print("If possible provide a BORN.DAT from the calculation",file=sys.stderr)
            # print("WARNING! WARNING! WARNING! WARNING! WARNING!",file=sys.stderr)
            self._read_born_charges_from_output(line)
        return

    def _read_born_charges_from_born_dat(self,filename):
        """Read born charges from a .born.dat file and store them in the instance variable.

        This method reads born charges for each ion from a specified file, which contains 
        data in three lines per ion. It transposes the read charge matrix and stores 
        the result in the `born_charges` instance variable.

        Parameters
        ----------
        filename : str
            The path to the .born.dat file to be read.

        Returns
        -------
        None

        """        
        fd2 = pdielec_io(filename, "r")
        self.born_charges = []
        for _i in range(self.nions):
            b = []
            line = fd2.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = fd2.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = fd2.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            b_np = np.array(b)
            bt_np = b_np.T
            b = bt_np.tolist()
            self.born_charges.append(b)
        fd2.close()
        return

    def _read_born_charges_from_output(self, line):
        """Read Born charges from output file.

        This method is responsible for parsing and reading Born effective charges from a given output file, beginning from the current line position of the file descriptor. It stores the read charges into the `born_charges` list attribute.

        Parameters
        ----------
        line : str
            The current line (not used)

        Returns
        -------
        None

        Note
        ----
        This method directly modifies the `born_charges` attribute of the class instance by appending the read Born charges. It assumes that the file is correctly positioned before the starts of the Born charges and that it's formatted in a specific way expected by the logic of this method.

        """        
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.born_charges = []
        for _i in range(self.nions):
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            b = []
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            self.born_charges.append(b)
            line = self.file_descriptor.readline()
        return

    def _read_ionic_dielectric(self, line):
        """Read ionic dielectric data from the file and store in an instance variable.

        This method reads the next four lines from the file associated with the
        instance's file descriptor. The second, third, and fourth lines are expected 
        to contain ionic dielectric data, each consisting of three floating-point 
        numbers. These numbers are stored in a list of lists in the instance variable 
        `zerof_static_dielectric`.

        Parameters
        ----------
        line : str
            This parameter is not used. The function reads lines directly from the
            instance's file descriptor.

        Returns
        -------
        None

        Notes
        -----
        The method directly modifies the instance variable `zerof_static_dielectric`.
        It is assumed that the file has been adequately positioned before this method
        is called, as it starts reading from the current position of the file
        descriptor's pointer.

        """        
        # Read the ionic contribution to the static dielectric
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric = []
        self.zerof_static_dielectric.append([float(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric.append([float(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        self.zerof_static_dielectric.append([float(f) for f in line.split()[0:3]])
        # a = np.array(ionic_dielectric)
        # self.zerof_static_dielectric = a.tolist()
        return

    def _read_lattice_vectors(self, line):
        """Read lattice vectors from a file and update the corresponding properties.

        This method reads three consecutive lattice vectors (a, b, c) from the current position in an open file, creates a new `UnitCell` object with these vectors, and updates the unit cell collection of the object. It also updates the number of cells, the volume of the last cell, and sets the fractional coordinates and element names for the last unit cell based on existing class attributes.

        Modifications to the instance attributes include:
        - Appending a new `UnitCell` object to the `unit_cells` list.
        - Updating the `ncells` attribute to match the length of the `unit_cells` list.
        - Setting the volume of the last unit cell in `unit_cells`.
        - Updating the fractional coordinates and element names of the last unit cell.

        Parameters
        ----------
        line : str
            The current line from the file where the method starts reading. This parameter is actually not used as the method immediately reads new lines from the file, implying a design choice where the `line` parameter could be omitted or revised.

        """        
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        avector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        bvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        cvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        self.unit_cells.append(UnitCell(avector, bvector, cvector,units="Angstrom"))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        # The fractional coordinates are specified before the lattice vectors
        self.unit_cells[-1].set_fractional_coordinates(self._fractional_coordinates)
        self.unit_cells[-1].set_element_names(self.species_list)
        return

    def _read_fractional_coordinates(self, line):
        """Read and store fractional coordinates from a file line.

        This method processes a line in the file to extract the number of ions and 
        then reads the successive lines to acquire the fractional coordinates, species, 
        and other related atomic information, storing them into the class attributes.

        Parameters
        ----------
        line : str
            The initial line to be processed for extracting the total number of ions.

        Returns
        -------
        None

        Attributes Updated
        ------------------
        nions : int
            The total number of ions.
        _fractional_coordinates : list of list of float
            The fractional coordinates of each ion.
        species : list of str
            The unique list of species types.
        species_list : list of str
            List of species corresponding to each ion.
        atom_type_list : list of int
            Index list specifying the species type for each ion.
        nspecies : int
            The total number of unique species.
        ions_per_type : list of int
            A count of ions per species type.

        """        
        self.nions = int(line.split()[12])
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self._fractional_coordinates = []
        self.species = []
        self.species_list = []
        self.atom_type_list = []
        for _i in range(self.nions):
            line = self.file_descriptor.readline()
            species = line.split()[3].capitalize()
            self.species_list.append(species)
            if species not in self.species:
                self.species.append(species)
            species_index = self.species.index(species)
            self.atom_type_list.append(species_index)
            self._fractional_coordinates.append([float(line.split()[4]), float(line.split()[5]), float(line.split()[6])])
        self.nspecies = len(self.species)
        self.ions_per_type = [ 0 for i in range(self.nspecies) ]
        for species_index in self.atom_type_list:
            self.ions_per_type[species_index] += 1
        return
