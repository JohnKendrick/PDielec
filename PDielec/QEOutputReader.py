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
"""QEOutputReader: Read the contents of a QE output file containing QE dynamical matrix."""
import math
import re
import xml.etree.ElementTree as ET

import numpy as np

from PDielec.Constants import amu, angs2bohr, au2GPa, hartree2ev
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell


class QEOutputReader(GenericOutputReader):
    """Read the contents of a QE output file containing QE dynamical matrix.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    filenames : str or list
        The filename or a list of filenames for QE output files.

    Notes
    -----
    This constructor initializes a ``GenericOutputReader`` with the given filenames,
    sets the type attribute to 'QE output', and initializes `_alat` attribute to None.

    """

    def __init__(self, filenames):
        """Initialize a QE output reader.

        Parameters
        ----------
        filenames : str or list
            The filename or a list of filenames for QE output files.

        Notes
        -----
        This constructor initializes a ``GenericOutputReader`` with the given filenames,
        sets the type attribute to 'QE output', and initializes `_alat` attribute to None.

        """
        GenericOutputReader.__init__(self, filenames)
        self.type                    = "QE output"
        self._alat                   = None
        self._alat_from_xml          = False
        return

    def _read_output_files(self):
        """Process the QE files.

        Notes
        -----
        Defines the search patterns for QE output files

        """
        # Define the search keys to be looked for in the files
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage["header"]   = (re.compile("Dynamical matrix file"), self._read_dyng_file)
        self.manage["lattice"]  = (re.compile("Basis vectors"), self._read_dyng_basis_vectors)
        self.manage["lattice2"]  = (re.compile("cubic"), self._read_dyng_basis_vectors)
        self.manage["lattice3"]  = (re.compile("CELL_PARAMETERS"), self._read_cell_parameters)
        self.manage["positions"]  = (re.compile("ATOMIC_POSITIONS"), self._read_fractional_coordinates)
        self.manage["dynamical"]  = (re.compile(" *Dynamical  Matrix in c"), self._read_dyng_dynamical)
        self.manage["epsilon"]  = (re.compile(" *Dielectric Tensor:"), self._read_dyng_epsilon)
        self.manage["charges"]  = (re.compile(" *Effective Charges E-U:"), self._read_born_charges)
        self.manage["energy_cutoff"]  = (re.compile(" *kinetic-energy cutoff"), self._read_energy_cutoff)
        self.manage["kpoints"]  = (re.compile(" *number of k points"), self._read_kpoints)
        self.manage["kpoint_grid"]  = (re.compile("K_POINTS automatic"), self._read_kpoint_grid)
        self.manage["electrons"]  = (re.compile("^ *number of electrons"), self._read_electrons)
        self.manage["energy"]  = (re.compile("^ *total energy  *="), self._read_energy)
        self.manage["alat"]  = (re.compile("^ *lattice parameter"), self._read_alat)
        self.manage["alat"]  = (re.compile("^ *A = "), self._read_alat2)
        self.manage["celldm1"]  = (re.compile("^ *celldm.1. ="), self._read_celldm1)
        self.manage["pressure"]  = (re.compile("^ *total *stress *.Ry"), self._read_pressure)
        self.manage["nions"]  = (re.compile("^ *number of atoms/cell"), self._read_nions)
        for f in self._outputfiles:
            if f.lower().endswith(".xml"):
                self._read_xml(f)
            else:
                self._read_output_file(f)
        return

    def _read_nions(self, line):
        #
        # Read from log file
        #
        """Read the number of ions from a line and update the object's nions attribute.

        The number of ions is expected to be the fifth element in the line when split by space.

        Parameters
        ----------
        line : str
            The line from which to extract the number of ions.

        Returns
        -------
        None

        Notes
        -----
        If 'debug' attribute of the object is True, this method will print the number of ions found.

        """
        self.nions = int(line.split()[4])
        if self.debug:
            print(f"_read_nions: nions={self.nions}" )
        return

    def _read_pressure(self, line):
        """Parse a pressure reading from a given line and update the object's pressure attribute.

        Parameters
        ----------
        line : str
            A string containing the pressure reading among other space-separated values. 
            The pressure value is expected to be at the 6th position in this string 
            (index 5 when zero-indexed) and will be divided by 10.0 to convert from kbar to GPa

        Returns
        -------
        None

        Notes
        -----
        This function updates the pressure attribute of the object it is called on and optionally prints the updated pressure value if the object's `debug` attribute is set to True.

        """
        self.pressure = float(line.split()[5])/10.0
        if self.debug:
            print(f"_read_pressure: pressure={self.pressure}" )
        return

    def _read_celldm1(self, line):
        #
        # Read from log file
        #
        """Process and update the class attribute _alat based on a line from a file.

        This function is intended to read and process a line that contains a unit of length. It updates the class attribute `_alat` after converting the value properly and checks for conformity within a tolerance. If `debug` is set to True in the class, it also prints the updated `_alat` value.

        Parameters
        ----------
        line : str
            The line read from a file which contains the value to be processed. It's expected that the relevant value is the third element of the line when split by whitespace, and it may contain a comma which should be ignored.

        Returns
        -------
        None

        Notes
        -----
        - Presumes the existence of a class attribute `_alat` which is updated by this function.
        - Assumes `angs2bohr` is a predefined constant available within the scope of this method for conversion or comparison purposes.
        - Uses a tolerance of 0.0001 for deciding whether the processed value matches `angs2bohr`, in which case it precisely sets `_alat` to `angs2bohr`.
        - The `debug` attribute of the class controls whether an update operation is printed or not.

        """
        string = line.split()[2]
        string = string.replace(",","")
        t = float(string)
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        # There are rounding errors when reading from the log file
        # So only read if there is no alternative
        if not self._alat_from_xml and self._alat is None:
            self._alat = t
        if self.debug:
            print(f"_read_celldm1: _alat={self._alat}" )
        return

    def _read_alat(self, line):
        #
        # Read from log file
        #
        """Read and set the lattice parameter from a given line of text.

        This method processes a line of text, aiming to extract and set the lattice parameter
        for the current object. The value is directly sought from a specified position in the line,
        then adjusted if necessary to match a predefined standard unit conversion threshold. 
        This method does not return any value.

        Parameters
        ----------
        line : str
            The line of text from which the lattice parameter is to be read.

        Returns
        -------
        None

        Notes
        -----
        The function checks if the extracted lattice parameter is within a certain threshold (`0.0001`)
        of a predefined conversion factor (`angs2bohr`). If so, it adjusts the lattice parameter to exactly
        match this conversion factor. This functionality assumes that `angs2bohr` is defined elsewhere
        and is accessible within the scope of this method. Also assumes that `self._alat` has been 
        initialized prior, possibly as `None`, and `self.debug` is a boolean indicating whether debug
        prints should be activated.

        """
        t = float(line.split()[4])
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        # There are rounding errors when reading from the log file
        # So only read if there is no alternative
        if not self._alat_from_xml and self._alat is None:
            self._alat = t
        if self.debug:
            print(f"_read_alat: _alat={self._alat}" )
        return

    def _read_alat2(self, line):
        #
        # Read from log file
        #
        """Read and set the lattice parameter from a given line.

        Parameters
        ----------
        line : str
            The line from which the lattice parameter is read. It is assumed that the lattice parameter can be found at the third whitespace-separated value in the line.

        Returns
        -------
        None

        Attributes Modified
        -------------------
        _alat
            The lattice parameter is set or adjusted within this instance based on the value read from `line`. It is adjusted to `angs2bohr` if the absolute difference between the read value and `angs2bohr` is less than 0.0001.

        Notes
        -----
        This method does not return any value. It updates the `_alat` attribute of the class instance with the value read from the input line, after applying necessary conversions and checks.

        """
        t = float(line.split()[2])
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        # There are rounding errors when reading from the log file
        # So only read if there is no alternative
        if not self._alat_from_xml and self._alat is None:
            self._alat = t
        if self.debug:
            print(f"_read_alat2: _alat={self._alat}" )
        return

    def _read_electrons(self, line):
        #
        # Read from log file
        #
        """Read and store the number of electrons from a line.

        This private method extracts and stores the number of electrons 
        from a given line, assuming the number of electrons is specified 
        in the 5th position of the line when split by whitespace. 
        If the debug attribute is set to True, it prints out the value 
        of electrons read.

        Parameters
        ----------
        line : str
            The line from which to read the number of electrons.

        Returns
        -------
        None

        Notes
        -----
        This method modifies the state of the object by setting the 
        'electrons' attribute based on the input line.

        """
        self.electrons = float(line.split()[4])
        if self.debug:
            print(f"_read_electrons: electrons={self.electrons}" )
        return

    def _read_energy(self, line):
        #
        # Read from log file
        #
        """Parse the energy values from a given line and set them to the instance variables.

        This function extracts energy-related values from a provided line of text, then calculates and assigns these values (in electron volts) to instance attributes for final energy without entropy and final free energy, while considering conversion factors. If debugging is enabled, it prints the final free energy.

        Parameters
        ----------
        line : str
            A string containing energy-related information, expected to follow a pre-determined format where the relevant energy value is the fourth element in a whitespace-separated list.

        Returns
        -------
        None

        See Also
        --------
        hartree2ev : Constant used for conversion from Hartree to electron volts.

        Notes
        -----
        - Assumes that 'hartree2ev' is a predefined conversion constant available in the scope where this function is defined.
        - Utilizes the fourth element in the whitespace-separated list from the input string, after splitting, for the relevant energy-related value.
        - The energy value is converted to electron volts, halved, and then assigned to both 'final_energy_without_entropy' and 'final_free_energy' attributes of the object.
        - If the 'debug' attribute of the object is True, the final free energy is printed.

        """
        self.final_energy_without_entropy = float(line.split()[3]) * hartree2ev / 2.0
        self.final_free_energy = float(line.split()[3]) * hartree2ev / 2.0
        if self.debug:
            print(f"_read_energy: energy={self.final_free_energy}" )
        return

    def _read_energy_cutoff(self, line):
        #
        # Read from log file
        #
        """Parse and set the energy cutoff from a line of text.

        Parameters
        ----------
        line : str
            The line of text containing the energy cutoff information. The expected
            format is such that the relevant numerical value is the fourth element (index 3)
            of the whitespace-separated entries.

        Returns
        -------
        None

        Notes
        -----
        This method directly modifies the `energy_cutoff` attribute of the class instance
        by setting it to the parsed value converted from Hartree to electron volts (eV) and
        divided by 2. The conversion factor `hartree2ev` is expected to be defined elsewhere
        in the code.

        This method optionally prints a debug message with the new energy cutoff value if
        the instance's `debug` attribute is `True`.

        """
        self.energy_cutoff = float(line.split()[3]) * hartree2ev / 2.0
        if self.debug:
            print(f"_read_energy_cutoff: energy_cutoff={self.energy_cutoff}" )
        return

    def _read_kpoints(self, line):
        #
        # Read from log file
        #
        """Read and store the number of k-points from a given line.

        This is an internal method that extracts the number of k-points from a provided line (assumed to contain this information in a predefined format) and stores it in the class instance. Optionally, it prints debug information if the debug mode is activated.

        Parameters
        ----------
        line : str
            The line from which the number of k-points is extracted. The line should contain the number of k-points at the 5th position (index 4 when split by whitespace).

        Returns
        -------
        None

        Notes
        -----
        This method modifies the state of the object by setting the `kpoints` attribute.

        """
        self.kpoints = int(line.split()[4])
        if self.debug:
            print(f"_read_kpoints: kpoints={self.kpoints}" )
        return

    def _read_kpoint_grid(self, line):
        #
        # Read from log file
        #
        """Read and set the k-point grid from a line in a file.

        This method reads a line from the file associated with the file_descriptor attribute, parses the first three elements of the line as floating-point numbers, and sets them as the k-point grid. If the debug attribute is set to True, it prints the k-point grid.

        Parameters
        ----------
        line : str
            The line from which to read the k-point grid. Note that while this parameter is named 'line', it is not used directly in the function as the function reads the next line from the file itself.

        Returns
        -------
        None

        Notes
        -----
        This function modifies the state of the object by setting the `kpoint_grid` attribute to the parsed k-point grid values.

        """
        line = self.file_descriptor.readline()
        self.kpoint_grid = [ int(f) for f in line.split()[0:3] ]
        if self.debug:
            print(f"_read_kpoint_grid kpoint_grid={self.kpoint_grid}" )
        return

    def _read_dyng_header(self):
        """Read and process the header line from the dynG file.

        This method updates the object's state by reading the header line from the set file descriptor and parsing it to extract the number of species (nspecies), the number of ions (nions), and the lattice parameter (alat) after converting it to Bohr units if it is within a tolerance range of the conversion factor.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        The function assumes that the file has been opened and that the file descriptor is stored in `self.file_descriptor`. The header line is expected to contain at least four space-separated values, with the first, second, and fourth values being integers representing the number of species, the number of ions, and a float representing the lattice parameter, respectively. The `angs2bohr` constant is assumed to be defined elsewhere in the code and is used as a threshold for converting the lattice parameter to Bohr units. If the `debug` attribute of the object is set to True, it prints the lattice parameter after conversion or identification.

        This method modifies the state of the object by setting `nspecies`, `nions`, and `_alat` based on the contents of the header line.

        """
        # Skip a couple of lines
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.nspecies = int(line.split()[0])
        self.nions    = int(line.split()[1])
        t             = float(line.split()[3])
        if abs(t - angs2bohr) < 1.0e-4:
            t = angs2bohr
        # Only overwrite alat if it hasn't come from the xml file
        if not self._alat_from_xml:
             self._alat = t
        if self.debug:
            print(f"_read_dyng_header alat={self._alat}" )
        return

    def _read_dyng_epsilon(self, line):
        """Read and process zero frequency optical permittivity data from the dynG file.

        This method reads the next four lines from a file, processes them to extract optical dielectric constant data, and stores this data in the instance's `zerof_optical_dielectric` attribute. Specifically, it skips the first line, then reads three lines of numeric values, where each line contains at least three numeric values which are converted to floats and stored in a 2D list assigned to `zerof_optical_dielectric`.

        Parameters
        ----------
        line : str
            The last line read

        Returns
        -------
        None

        Notes
        -----
        - This method modifies the instance attribute `zerof_optical_dielectric` directly.
        - If the instance attribute `debug` is set to True, this method prints the `zerof_optical_dielectric` attribute after it is updated.

        """
        self.file_descriptor.readline()
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric = []
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        linea = self.file_descriptor.readline().split()
        self.zerof_optical_dielectric.append([float(f) for f in linea[0:3]])
        if self.debug:
            print(f"_read_espilon zerof_optical_dielectric={self.zerof_optical_dielectric}" )
        return

    def _read_dyng_masses(self):
        """Read the masses of species from an dynG file and populate mass and species lists.

        This method parses lines from an already opened file object (`file_descriptor`), extracting species names and their corresponding masses, which are then adjusted by multiplying by 2 and dividing by the atomic mass unit (amu). The extracted species names are capitalized and stored, along with the calculated masses, in their respective class attributes. If debugging is enabled, it prints the loaded masses. It also calls `_read_dyng_coordinates` method before finishing, passing an empty string as an argument.

        Parameters
        ----------
        None 

        Returns
        -------
        None

        Attributes Modified
        -------------------
        - `masses_per_type` : list
            A list that this method populates with the calculated masses for each species.
        - `species` : list
            A list that this method populates with the capitalized species names.

        Notes
        -----
        - The atomic mass unit (`amu`) is defined elsewhere
        - It's assumed that `_read_dyng_coordinates` is a method within the same class that takes a single string argument.

        """
        self.masses_per_type = []
        self.species      = []
        for _i in range(self.nspecies):
            linea = self.file_descriptor.readline().replace("'", "").split()
            self.species.append(linea[1].capitalize())
            # The factor of two is because au in pwscf are half mass of electron
            self.masses_per_type.append(float(linea[2])*2/amu)
        if self.debug:
            print(f"_read_dyng_masses masses={self.masses_per_type}" )
            print(f"_read_dyng_masses species={self.species}" )
        return

    def _read_dyng_dynamical(self,line):
        """Read and process the dynamical matrix from the dynG file.

        This method reads a specific format from the current position of the file descriptor,
        computes the Hessian matrix based on the read values, and updates the dynamical matrix.
        The routine only reads the gamma point matrix

        Parameters
        ----------
        line : str
            The last line read

        Returns
        -------
        None

        """
        nmodes = self.nions*3
        hessian = np.zeros((nmodes, nmodes))
        self.file_descriptor.readline()
        linea = self.file_descriptor.readline().split()
        # We only want to read the hessian at gamma
        q = [float(q) for q in linea[3:6]]
        qsum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2]
        if qsum > 0.0001:
            return
        # We read the hessian and store the mass weighted matrix
        linea = self.file_descriptor.readline().split()
        for a in range(self.nions):
            for b in range(self.nions):
                self.file_descriptor.readline()
                for ixyz in range(3):
                    ipos = a*3 + ixyz
                    linea = self.file_descriptor.readline().split()
                    for jxyz in range(3):
                        jpos = b*3 + jxyz
                        # factor of 0.5 'cos of au units in pwscf
                        hessian[ipos, jpos] = 0.5*float(linea[2*jxyz])/(amu*math.sqrt(self.masses[a]*self.masses[b]))
                    # end for jxyz
                # end for ixyz
            # end for b
        # end for a
        self._dynamical_matrix(hessian)
        if self.debug:
            print("_read_dyng_dynamical")
        return

    def _read_born_charges(self, line):
        """Read and parse Born charges from a file.

        This internal method reads the Born charges of `nions` ions from the file specified by `file_descriptor` and updates the `born_charges` attribute of the class instance.

        Parameters
        ----------
        line : str
            Dummy argument which is overwritten
        self.nions : int
            The number of ions to read Born charges for.

        Returns
        -------
        None

        Notes
        -----
        This method reads three lines for each ion, each line corresponding to one of the three Cartesian coordinates (x, y, z). It expects each line to contain three space-separated values, which represent the components of the Born charge tensor for that direction.

        This method directly modifies the `born_charges` attribute of the class instance, appending a 3x3 list (representing the Born charge tensor) for each ion.

        If `self.debug` is True, it prints a debug message upon completion.

        """
        self.born_charges = []
        line = self.file_descriptor.readline()
        for _i in range(self.nions):
            b = []
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            self.born_charges.append(b)
        if self.debug:
            print("_read_born_charges")
        return

    def _read_cell_parameters(self, line):
        """Read and parse cell parameters from a given line in the input file.

        This method parses the cell parameters from a specific line in the input file, updating the internal
        state with the scaled lattice vectors and their corresponding unit cell added to the list of unit cells.
        Additionally, this method updates the total number of cells and the volume of the last unit cell read.

        Parameters
        ----------
        line : str
            A string representing a line from the input file that contains cell parameters.

        Returns
        -------
        None

        Notes
        -----
        The function expects the `line` input to contain the lattice parameter as its third value after removing parentheses.
        Each subsequent line read from the input file (after the parameter line) is expected to represent the lattice vectors,
        scaled by the lattice parameter and converted from Angstroms to Bohr radii.
        The method directly modifies the object's state by updating the `unit_cells` list, `ncells`, and `volume` attributes.

        """
        line = line.replace( "(", " ")
        line = line.replace( ")", " ")
        linea = line.split()
        # Overwrite alat here as it is more accurate than anywhere else in the log file
        if not self._alat_from_xml and len(linea) > 2:
            self._alat = float(linea[2])
        linea = self.file_descriptor.readline().split()
        avector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        bvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        cvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        self.unit_cells.append(UnitCell(avector, bvector, cvector,units="Angstrom"))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        if self.debug:
            print(f"_read_cell_parameters: volume={self.volume}")
        return

    def _read_dyng_basis_vectors(self):
        """Read and process lattice vectors from a dynG file, then updates unit cell information.

        This private method reads the next three lines from an open file, each representing the a, b, and c vectors of a unit cell, respectively. Each vector is scaled according to the attribute `_alat` and a conversion factor from Angstroms to Bohrs. It then updates the list of unit cells with a new `UnitCell` object created from these vectors, recalculates the total number of cells, updates the volume attribute to the volume of the latest unit cell, and, if debug mode is active, prints the volume. Finally, it calls another method to read masses.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Attributes Updated
        ------------------
        unit_cells : list
            A list of `UnitCell` objects representing the unit cells read from the file. Updated by appending a new `UnitCell`.
        ncells : int
            The total number of unit cells, updated after reading new lattice vectors.
        volume : float
            The volume of the last read unit cell, updated after appending a new `UnitCell`.

        """
        linea = self.file_descriptor.readline().split()
        if linea[0] == "Basis" or linea[0] == "cubic":
            linea = self.file_descriptor.readline().split()
        avector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        bvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        linea = self.file_descriptor.readline().split()
        cvector = [float(f)*self._alat/angs2bohr for f in linea[0:3]]
        self.unit_cells.append(UnitCell(avector, bvector, cvector,units="Angstrom"))
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        if self.debug:
            print(f"_read_lattic_vectors: volume={self.volume}")
        return

    def _read_fractional_coordinates(self,line):
        """Read the log file and sets the fractional coordinates and element names for the last unit cell in the object.

        This function iterates over a specified number of lines (defined by `self.nions`) from the object's file descriptor, parsing each line into species names and their fractional coordinates. These species names and coordinates are then used to set the corresponding properties for the last unit cell managed by the object. Additional side-effects include updating the total number of unit cells (`self.ncells`) and the volume of the last unit cell (`self.volume`), with the option to print the volume if debugging is enabled.

        Parameters
        ----------
        line : str
            The initial line from which to start reading the fractional coordinates. The actual parameter is not used within the function but symbolically represents the starting point for reading. 

        Returns
        -------
        None

        Notes
        -----
        - The function directly modifies the object's state by setting properties of the last unit cell (`self.unit_cells[-1]`) and updating `self.ncells` and `self.volume`.
        - Only operates if `self.nions` is greater than 0.
        - The function is intended to be called in a context where `self.file_descriptor` is already open and positioned correctly.

        """
        if self.nions <= 0 or len(self.unit_cells) <=0 :
            if self.debug:
                print(f"_read_fractional_coordinates: nions={self.nions}")
                print(f"_read_fractional_coordinates: number of cells={len(self.unit_cells)}")
                print("_read_fractional_coordinates: Unable to store fractional coordinates aborting")
            return
        species_list = []
        fractional_coordinates = []
        for _i in range(self.nions):
            linea = self.file_descriptor.readline().split()
            species_list.append( linea[0] )
            fractional_coordinates.append([float(linea[1]), float(linea[2]), float(linea[3])])
        self.unit_cells[-1].set_fractional_coordinates(fractional_coordinates)
        self.unit_cells[-1].set_element_names(species_list)
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        if self.debug:
            print(f"_read_fractional_coordinates: volume={self.volume}")
        return

    def _read_dyng_coordinates(self):
        """Read dynamic coordinates from a dynG file and update related properties.

        This method reads atomic positions from a line, updating the mass, atomic type list, ions per type, and species list. It scales the coordinates according to the conversion from Angstroms to Bohr units and updates the corresponding `unit_cells` properties.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Attributes Updated
        ------------------
        self.masses : list
            Updated list of masses based on the atomic positions read.
        self.atom_type_list : list
            Updated list of atom types corresponding to each atom based on the read positions.
        self.ions_per_type : list
            Updated list with the count of ions per type.
        self.unit_cells : list
            List of unit cells, where the last one is updated with new xyz coordinates and element names based on read data.
        self.ncells : int
            Updated count of unit cells.
        self.volume : float
            Updated volume of the last unit cell in the list.
    
        See Also
        --------
        `set_xyz_coordinates` : Method used to update xyz coordinates of the last unit cell.
        `set_element_names` : Method used to update element names of the last unit cell.

        """
        self.masses = []
        self.atom_type_list = []
        self.ions_per_type = [ 0 for i in range(self.nspecies) ]
        species_list = []
        xyz_coordinates = []
        # It took a long time to work out that alat is in bohr
        const = self._alat/angs2bohr
        for _i in range(self.nions):
            linea = self.file_descriptor.readline().split()
            species_index = int(linea[1])
            xyz_coordinates.append([const*float(linea[2]), const*float(linea[3]), const*float(linea[4])])
            self.masses.append(self.masses_per_type[species_index-1])
            self.atom_type_list.append(species_index-1)
            self.ions_per_type[species_index-1] += 1
            species_list.append(self.species[species_index-1])
        self.unit_cells[-1].set_xyz_coordinates(xyz_coordinates,units="Angstrom")
        self.unit_cells[-1].set_element_names(species_list)
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        if self.debug:
            print(f"_read_dyng_coordinates: volume={self.volume}")
        return

    def _read_xml(self, filename):
        """Process the QE xml files.

        Parameters
        ----------
        filename : str
            The file name of the xml file to be processed

        Notes
        -----
        The routine only handles the pwscf output file pwscf.xml
        This file is created when a single point calculation is performed prior to the born calculation
        The routine should also be able to parse an optimisation run.
        In this case there is a 'steps' entry in the xml

        """
        # For each file generate an xml tree to process.
        tree = ET.parse(filename)
        root = tree.getroot()
        if self.debug:
            print("QEOutputReader _read_xml: filename = ",filename)
            print("QEOutputReader _read_xml: tree = ",tree)
            print("QEOutputReader _read_xml: root = ",root)
        inputxml = root.find("input")
        outputxml = root.find("output")
        stepsxml = None
        if outputxml is not None:
            stepsxml = outputxml.find("steps")
        # Handle the case of an optimisation or a single point calculation
        use_this_xml = stepsxml if stepsxml is not None else outputxml
        tensorsxml = root.find("EF_TENSORS")
        if self.debug:
            print("QEOutputReader _read_xml: inputxml = ",inputxml)
            print("QEOutputReader _read_xml: outputxml = ",outputxml)
            print("QEOutputReader _read_xml: stepsxml = ",stepsxml)
            print("QEOutputReader _read_xml: use_this_xml = ",use_this_xml)
            print("QEOutputReader _read_xml: tensorsxml = ",tensorsxml)
        # Handle input
        if inputxml is not None:
            self._basis(inputxml.findall("basis")[-1])
            self._atomic_species(inputxml.findall("atomic_species")[-1])
            self._kpoints(inputxml.findall("k_points_IBZ")[-1])
        # Handle output 
        if outputxml is not None:
            self._pressure(outputxml.findall("stress")[-1])
            self._band_structure(outputxml.findall("band_structure")[-1])
        # Handle the case of optimisation
        if use_this_xml is not None:
            self._atomic_structure(use_this_xml.findall("atomic_structure"))
            self._total_energy(use_this_xml.findall("total_energy"))
        # Handle ef tensors
        if tensorsxml is not None:
            self._dielectric_constant(tensorsxml.find("DIELECTRIC_CONSTANT"))
            self._effective_charges(tensorsxml.find("EFFECTIVE_CHARGES_EU"))
        return

    def _effective_charges(self,effective_charges_xml):
        """Process the effective_charges element(s).

        For a given atom the ordering of the elements in the tensors.xml file appears to be:
        d2E/dFxdX d2E/dFydX d2E/dFzdX        0 1 2
        d2E/dFxdY d2E/dFydY d2E/dFzdY        3 4 5
        d2E/dFxdZ d2E/dFydZ d2E/dFzdZ        6 7 8

        The order used with PDielec is:
        d2E/dFxdX d2E/dFxdY d2E/dFxdZ        0 3 6
        d2E/dFydX d2E/dFydY d2E/dFydZ        1 4 7
        d2E/dFzdX d2E/dFzdY d2E/dFzdZ        2 5 8

        Parameters
        ----------
        effective_charges_xml : an xml element
            An effective_charges xml element

        Set
        ---
        self.born_charges
            The born_charges 

        """
        self.born_charges = []
        charges_string = effective_charges_xml.text.strip()
        charges = [ float(f) for f in charges_string.split() ]
        index = 0
        for _i in range(self.nions):
            b = []
            b.append([ charges[index+0], charges[index+3], charges[index+6] ])
            b.append([ charges[index+1], charges[index+4], charges[index+7] ])
            b.append([ charges[index+2], charges[index+5], charges[index+8] ])
            index = index + 9
            self.born_charges.append(b)
        if self.debug:
            print(f"_basis: born_charges={self.born_charges}")
        return

    def _dielectric_constant(self,dielectric_constant_xml):
        """Process the dielectric_constant element(s).

        Parameters
        ----------
        dielectric_constant_xml : an xml element
            A dielectric_constant xml element

        Set
        ---
        self.zerof_optical_dielectric
            The zero frequency optical dielectric

        """
        eps_string = dielectric_constant_xml.text.strip()
        eps = [ float(f) for f in eps_string.split() ]
        eps = np.array(eps).reshape(3,3)
        self.zerof_optical_dielectric = eps.tolist()
        if self.debug:
            print(f"_basis: zerof_optical_dielectric={self.zerof_optical_dielectric}")
        return

    def _basis(self,basis_xml):
        """Process the basis element(s).

        Store the energy in eV

        Parameters
        ----------
        basis_xml : a basis xml element
            A basis xml element

        Set
        ---
        self.energy_cutoff
            The energy cutoff in hartree
        self.rho_cutoff
            The charge density cutoff in hartree

        """
        ecutwfc_xml = basis_xml.find("ecutwfc")
        ecutrho_xml = basis_xml.find("ecutrho")
        self.energy_cutoff = float( ecutwfc_xml.text )
        self.rho_cutoff    = float( ecutrho_xml.text )
        if self.debug:
            print(f"_basis: energy_cutoff={self.energy_cutoff}")
            print(f"_basis: rho_cutoff={self.rho_cutoff}")
        return

    def _total_energy(self,total_energy_xml):
        """Process the total_energy element(s).

        Store the energy in eV

        Parameters
        ----------
        total_energy_xml : a list of xml element
            A list of total_energy elements

        Set
        ---
        self.energies :  a list of floats
            The final energies in eV of each step in an optimisation
        self.final_energy_without_entropy : float
            The final energy in eV
        self.final_free_energy : float
            The final energy in eV

        """
        self.energies = []
        for step_xml in total_energy_xml:
            etot_xml = step_xml.find("etot")
            self.final_energy_without_entropy = hartree2ev * float( etot_xml.text )
            self.final_free_energy = self.final_energy_without_entropy
            self.final_free_energies.append(self.final_free_energy)
            self.final_energies_without_entropy.append(self.final_free_energy)
        if self.debug:
            print(f"_total_energy: final_energy_without_entropy={self.final_energy_without_entropy}")
            print(f"_total_energy: final_free_energy={self.final_free_energy}")
            print(f"_total_energy: energies={self.energies}")
        return

    def _band_structure(self,band_structure_xml):
        """Process the band_structure element.

        Parameters
        ----------
        band_structure_xml : an xml element
            A band structure element

        Set
        ---
        self.electrons : int
            The number of electrons

        """
        nelec_xml = band_structure_xml.find("nelec")
        self.electrons = int(float(nelec_xml.text)+0.00000001)
        if self.debug:
            print(f"_band_structure: electrons={self.electrons}")
        return

    def _atomic_structure(self,atomic_structure_xml):
        """Process the atomic_structure element(s).

        Parameters
        ----------
        atomic_structure_xml : a list of xml elements
            A list of atomic structure elements

        Set
        ---
        self.nions : int
            The number of ions
        self.unit_cells : unitCell
            A list of the unit cells
        self.volume :  float
            A volume of the last unit cell
        self.ncell :  int
            A number of unit cells

        """
        for structure_xml in atomic_structure_xml:
            self.nions = int(structure_xml.attrib["nat"])
            self._alat  = float(structure_xml.attrib["alat"])
            self._alat_from_xml = True
            species_list = []
            xyz_coordinates = []
            # First find the cell dimensions and volume
            cell_xml = structure_xml.find("cell")
            a1 = [ float(f)/angs2bohr for f in cell_xml[0].text.split() ]
            a2 = [ float(f)/angs2bohr for f in cell_xml[1].text.split() ]
            a3 = [ float(f)/angs2bohr for f in cell_xml[2].text.split() ]
            self.unit_cells.append(UnitCell(a1, a2, a3, units="Angstrom"))
            # Now find the fractional coordinates
            positions_xml = structure_xml.find("atomic_positions")
            self.atom_type_list = []
            self.ions_per_type = [ 0 for i in range(self.nspecies) ]
            for atom_xml in positions_xml.findall("atom"):
                species = atom_xml.attrib["name"]
                species_index = self.species.index(species)
                species_list.append( species )
                coords = [ float(f)/angs2bohr for f in atom_xml.text.split() ]
                xyz_coordinates.append(coords)
                self.atom_type_list.append(species_index)
                self.ions_per_type[species_index] += 1
            self.unit_cells[-1].set_xyz_coordinates(xyz_coordinates)
            self.unit_cells[-1].set_element_names(species_list)
            # Book keeping
            self.ncells = len(self.unit_cells)
            self.volume = self.unit_cells[-1].getVolume("Angstrom")
            if self.debug:
                print(f"_atomic_structure: volume={self.volume}")
        return
            
    def _atomic_species(self,atomic_species_xml):
        """Process the atomic_species element.

        Parameters
        ----------
        atomic_species_xml : an xml element
            The atomic species element

        Set
        ---
        self.nspecies
            The number of species
        self.species
            A list of the names of the species
        self.masses_per_type
            A list of the masses of the species

        """
        self.nspecies = int(atomic_species_xml.attrib["ntyp"])
        self.masses_per_type = []
        self.species      = []
        for child in atomic_species_xml.findall("species"):
            self.species.append( child.attrib["name"] )
            self.masses_per_type.append ( float(child.find("mass").text) )
        if self.debug:
            print(f"_atomic_species: nspecies={self.nspecies}")
            print(f"_atomic_species: species={self.species}")
            print(f"_atomic_species: masses_per_type={self.masses_per_type}")
        return
            

    def _kpoints(self,kpoints_xml):
        """Process the k_points element.

        Parameters
        ----------
        kpoints_xml : an xml element
            The k-point definition using Monkhorst-Pack

        Set
        ---
        self.kpoints
            A list of the number of kpoints to be used in a Monkhorst-Pack grid

        """
        kpt_dict = kpoints_xml[0].attrib
        self.kpoint_grid[0] = int(kpt_dict["nk1"])
        self.kpoint_grid[1] = int(kpt_dict["nk2"])
        self.kpoint_grid[2] = int(kpt_dict["nk3"])
        if self.debug:
            print(f"_read_kpoints: kpoint_grid={self.kpoint_grid}" )
        return

    def _pressure(self,pressure_xml):
        """Process the pressure xml element.

        Parameters
        ----------
        pressure_xml : an xml element
            The pressure is calculated for the average of the trace of the stress tensor

        Set
        ---
        self.pressure
            The pressure in GPa
            The xml file holds the pressure in Hartree/Bohr3
            Conversion to SI units is done by autoGPa from Constants.py

        """
        stress = [ float(f) for f in pressure_xml.text.split()]
        # Take the trace of the (flat) matrix and convert from atomic units to GPa
        self.pressure = au2GPa * (stress[0] + stress[4] + stress[8]) / 3.0 
        if self.debug:
            print(f"_pressure: pressure={self.pressure}" )
        return

    def _read_dyng_file(self,line):
        """Read a dynG file.

        This method the dynG output file from a qe ph.x calculation of IR intensities

        Parameters
        ----------
        line : str
            The first line of input in the dyng file

        Returns
        -------
        None

        Notes
        -----
        Calls
        - `_read_dyng_header`
        - `_read_dyng_basis_vectors`
        - `_read_dyng_masses`
        - `_read_dyng_coordinates`

        """
        if self.debug:
            print("_read_dyng:" )
        # Read the header
        self._read_dyng_header()
        # Read the basis vectors
        self._read_dyng_basis_vectors()
        # Read the species masses
        self._read_dyng_masses()
        # Read the fractional coordinates
        self._read_dyng_coordinates()
