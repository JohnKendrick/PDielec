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
"""VASP output reader."""

import re
import xml.etree.ElementTree as ET

import numpy as np

from PDielec.Calculator import calculate_normal_modes_and_frequencies
from PDielec.Constants import atomic_number_to_element, hertz
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell


def read_xml_element(ele):
    """Read an xml element and return a value.

    Parameters
    ----------
    ele : an xml element
        An xml element

    Returns
    -------
    name : str
        The name of the values returned
    values : float, int or bool (or a list of them)
        The return type depends on the contents of the element

    """
    if ele.tag == "separator":
        return None, None
    if ele.tag == "parameters":
        return None, None
    name = ele.attrib["name"]
    values = []
    if "type" in ele.attrib and ele.attrib["type"] == "int":
        values = [ int(f) for f in ele.text.strip().split() ]
    elif "type" in ele.attrib and ele.attrib["type"] == "logical":
        values = [ f == "T" for f in ele.text.strip().split() ]
    elif "type" in ele.attrib and ele.attrib["type"] == "string":
        values = [ None ] if ele.text is None else [ ele.text.strip() ]
    else:
        values = [ float(f) for f in ele.text.strip().split() ]
    if ele.tag == "i":
        values = values[0]
    return name,values


def myfloat(string):
    """Replace float(), return a large number if it finds a * in the string.

    Parameters
    ----------
    string : str
        The string to be analysed as a float

    Returns
    -------
    float
        A standard float conversion of the input or a large number if '*' is found in the input.

    """
    if "*" in string:
        return 9999.999
    return float(string)

class VaspOutputReader(GenericOutputReader):
    """Read the contents of a directory containing VASP input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    names : list of str
        A list of file names to be used in processing the DFT output files

    Attributes
    ----------
    type : str
        A string attribute that specifies the type of output, set to 'Vasp output'.
    _pspots : dict
        A private dictionary for the pseudo potentials
    _ediff : float
        A private float for the change in energy on each iteration
    _pulay : float
        The Pulay pressure
    _ibrion : int
        The value of ibrion
    _potim : float
        The value of potim

    """

    def __init__(self, names):
        """Initialize a new instance of the class.

        Parameters
        ----------
        names : list of str
            A list of file names to be used in processing the DFT output files

        Attributes
        ----------
        type : str
            A string attribute that specifies the type of output, set to 'Vasp output'.
        _pspots : dict
            A private dictionary for the pseudo potentials
        _ediff : float
            A private float for the change in energy on each iteration
        _pulay : float
            The Pulay pressure
        _ibrion : int
            The value of ibrion
        _potim : float
            The value of potim

        """        
        GenericOutputReader.__init__(self, names)
        self.type                    = "Vasp output"
        self._pspots                  = {}
        self._ediff                   = 0.0
        self._pulay                  = None
        self._ibrion                  = 0
        self._potim                   = 0.0
        return

    def _read_output_files(self):
        """Read the vasp files in the directory.

        Define the search strings neeeded and the associated routines will read the required values
        Process the files accordingly

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.manage = {}   # Empty the dictionary matching phrases
        self.manage["ionspertype"]  = (re.compile("   ions per type ="), self._read_ionspertype)
        # self.manage['masses_skip']  = (re.compile('  Mass of Ions in am'), self._read_skip4)
        self.manage["pspots"]       = (re.compile(" POTCAR:"), self._read_pspot)
        self.manage["arrays"]       = (re.compile(" Dimension of arrays:"), self._read_array_dimensions)
        #self.manage['masses']       = (re.compile('   POMASS ='), self._read_masses)
        self.manage["spin"]         = (re.compile("   ISPIN  = "), self._read_spin)
        self.manage["encut"]        = (re.compile("   ENCUT  = "), self._read_encut)
        self.manage["ediff"]        = (re.compile("   EDIFF  = "), self._read_ediff)
        self.manage["ibrion"]       = (re.compile("   IBRION = "), self._read_ibrion)
        self.manage["potim"]        = (re.compile("   POTIM  = "), self._read_potim)
        self.manage["nelect"]       = (re.compile("   NELECT = "), self._read_nelect)
        self.manage["lattice"]      = (re.compile("  volume of cell :"), self._read_lattice_vectors)
        self.manage["fractional"]   = (re.compile(" position of ions in fractional coordinates"), self._read_fractional_coordinates)
        self.manage["forces"]       = (re.compile(" POSITION  *TOTAL-FORCE"), self._read_forces)
        self.manage["energy"]       = (re.compile("  FREE ENERGIE OF THE ION"), self._read_energy)
        self.manage["magnet"]       = (re.compile(" number of electron "), self._read_magnet)
        self.manage["pressure"]     = (re.compile("  external pressure ="), self._read_external_pressure)
        self.manage["skip1"]        = (re.compile(" old parameters found"), self._read_skip4)
        self.manage["staticDielectric"]  = (re.compile(" MACROSCOPIC STATIC DIELECTRIC TENSOR .including"), self._read_static_dielectric)
        self.manage["staticIonic"]  = (re.compile(" MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONT"), self._read_ionic_dielectric)
        self.manage["bornCharges"]  = (re.compile(" BORN EFFECTIVE CHARGES"), self._read_born_charges)
        self.manage["eigenvectors"] = (re.compile(" Eigenvectors and eige"), self._read_eigenvectors)
        self.manage["eigenskip"]    = (re.compile(" Eigenvectors after division"), self._read_skip4)
        self.manage["elastic"]      = (re.compile(" TOTAL ELASTIC MODULI"), self._read_elastic_constants)
        self.manage["kpoint"]      = (re.compile("^Gamma"), self._read_kpoint_grid)
        self.manage["species"]      = (re.compile("^ *Atomic configuration"), self._read_species)
        self.manage["newmasses"]  = (re.compile("  Mass of Ions in am"), self._read_newmasses)
        for f in self._outputfiles:
            if f.lower().endswith(".xml"):
                self._read_xml(f)
            else:
                self._read_output_file(f)
        return

    def _read_forces(self, line):
        """Read force data from a file and update the iterations datasets with maximum and RMS forces.

        After execution, `self.iterations` will be updated with the new `max_force` and `rms_force` values for the current iteration.

        Parameters
        ----------
        line : str
            The current line of the file to start reading from.

        Returns
        -------
        None

        """        
        line = self.file_descriptor.readline()
        maxf = 0.0
        rmsf = 0.0
        for _i in range(self.nions):
            line = self.file_descriptor.readline()
            forces = [ float(f) for f in line.split()[3:6] ]
            for f in forces:
                rmsf += f*f
                maxf = max(maxf, abs(f))
                # end if
            # end for f
        #end for i
        if "max_force" not in self.iterations:
            self.iterations["max_force"] = []
        self.iterations["max_force"].append(maxf)
        if "rms_force" not in self.iterations:
            self.iterations["rms_force"] = []
        self.iterations["rms_force"].append(rmsf)
        return

    def _read_species(self, line):
        """Read species from a file and updates species information.

        This method reads the number of lines for a species, computes the total charge from these lines, and appends the corresponding element to the species list. It utilizes the atomic_number_to_element dictionary to map the computed total charge to an element symbol.

        Parameters
        ----------
        line : str
            The initial line to start processing from, typically provides the number of lines for a species.

        Returns
        -------
        None

        See Also
        --------
        atomic_number_to_element : A dictionary mapping total charges to element symbols.

        """        
        line = self.file_descriptor.readline()
        nlines = int(line.split()[0])
        line = self.file_descriptor.readline()
        zcharge = 0.0
        for _i in range(nlines):
            line = self.file_descriptor.readline()
            zcharge = zcharge + float(line.split()[4])
        self.species.append(atomic_number_to_element[int(zcharge+0.001)])
        self.nspecies = len(self.species)
        return

    def _read_kpoint_grid(self, line):
        """Read and parse the k-point grid dimensions from a file line.

        This method assumes that the current line from which `line` is read contains three integers that represent the dimensions of the k-point grid. These integers are read, converted to integers, and stored in the instance variable `kpoint_grid`.

        Parameters
        ----------
        line : str
            The current line of the file to read from. This parameter is not used, and the actual line reading happens from `self.file_descriptor`.

        Returns
        -------
        None

        Notes
        -----
        - This method modifies the state of the object by setting the `kpoint_grid` attribute.
        - The method reads the next line from `self.file_descriptor`, not from the `line` parameter.

        """        
        line = self.file_descriptor.readline()
        self.kpoint_grid = [int(f) for f in line.split()[0:3] ]
        return

    def _read_ionspertype(self, line):
        """Read and store the number of ions per species from a given line.

        This method processes a string that contains information about the number of ions for each type/species. It updates the attributes `ions_per_type` and `nspecies` of the class.

        Parameters
        ----------
        line : str
            The line from which the ions per type information will be extracted. The relevant data is expected to start from the fifth element when splitting the line by spaces.

        Returns
        -------
        None

        """        
        self.ions_per_type = [int(i) for i in line.split()[4:]]
        self.nspecies = len(self.ions_per_type)
        return

    def _read_newmasses(self, line):
        """Read and process the "newmasses" from a line in the file descriptor.

        This internal method parses a single line from the file descriptor to read the masses associated with different atom types in a molecular system. Based on the number of species (`nspecies`) and ions per type, it updates `masses_per_type`, `masses`, `atom_type_list`, and `species_list` lists with the parsed and computed information.

        Parameters
        ----------
        line : str
            Line is a string (not used)

        Returns
        -------
        None
            The string line from which the new masses information is to be read. Note that this parameter is initially passed but not used directly, as the function immediately reads a new line from `self.file_descriptor`.

        """        
        self.masses_per_type = []
        line = self.file_descriptor.readline()
        mass_string = line[12:]
        start = 0
        increment = 6
        for _i in range(self.nspecies):
            mass = mass_string[start:start+increment]
            self.masses_per_type.append(float(mass))
            start = start + increment
        self.masses = []
        self.atom_type_list = []
        self.species_list = []
        for k, mass in enumerate(self.masses_per_type):
            n = self.ions_per_type[k]
            for _i in range(n):
                self.atom_type_list.append(k)
                self.masses.append(mass)
                self.species_list.append(self.species[k])
            # end loop over i
        # end look over current know types
        return

    def _read_masses(self, line):
        """Parse and store the mass value from a line of text.

        This method extracts the third element from a whitespace-separated list (assumed
        to be a mass), removes semicolons, converts it to a float, and appends it to the 
        instance's list of masses.

        Parameters
        ----------
        line : str
            A string from which the mass value is to be extracted. It is expected to
            be whitespace-separated with the mass as the third element, potentially
            followed by a semicolon.

        Returns
        -------
        None

        """        
        mass_string = line.split()[2]
        mass_string = mass_string.replace(";", "")
        self.masses_per_type.append(float(mass_string))
        return

    def _read_eigenvectors(self, line):
        """Read and process eigenvectors from a file.

        This method extracts frequencies and mass-weighted normal modes from a read file, storing them in the class's frequencies and mass_weighted_normal_modes attributes, respectively. It assumes the file follows a specific format where frequencies and their corresponding eigenvectors are listed, considering imaginary frequencies if present.

        Parameters
        ----------
        line : str
            This method operates on the instance's file_descriptor and nions attributes, assuming they are correctly initialized.

        Returns
        -------
        None

        """        
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.frequencies = []
        self.mass_weighted_normal_modes = []
        n = 3 * self.nions
        for _i in range(n):
            line = self.file_descriptor.readline()
            line = self.file_descriptor.readline()
            imaginary = (line.split()[1] == "f/i=")
            freq = -float(line.split()[6]) if imaginary else float(line.split()[7])
            # end if
            self.frequencies.append(freq)
            line = self.file_descriptor.readline()
            a = []
            for _j in range(self.nions):
                line = self.file_descriptor.readline()
                a.append([float(line.split()[3]), float(line.split()[4]), float(line.split()[5])])
            # end for j
            self.mass_weighted_normal_modes.append(a)
        # end of for i in range(n)
        return

    def _read_born_charges(self, line):
        """Read the born charges from the OUTCAR file.

        Each row of the output refers to a given field direction and each column in the row refers to the atomic displacement. The output is arranged as: ::

            [[a1x a1y a1z]
             [a2x a2y a2z]
             [a3x a3y a3z]]

        where 1, 2, 3 are the field directions and x, y, z are the atomic displacements.

        Parameters
        ----------
        line : str 
            Not used and overridden

        Returns
        -------
        numpy.ndarray
            A numpy array containing the born charges with its shape determined by the field directions and atomic displacements as described above.

        """
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.born_charges = []
        for _i in range(self.nions):
            line = self.file_descriptor.readline()
            b = []
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            b.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            line = self.file_descriptor.readline()
            self.born_charges.append(b)
        return

    def _read_elastic_constants(self, line):
        """Read and process elastic constants from a file.

        This method reads lines from an opened file pointed by the instance's file
        descriptor, extracts numeric values assumed to be elastic constants, and
        stores a processed version of these constants in the instance attribute
        `elastic_constants`.

        Parameters
        ----------
        line : str
            The initial line from where to start reading the elastic constants in
            the file. This parameter is somewhat misleading in this context, as the
            function does not actually use this input variable but rather immediately
            overrides it with `self.file_descriptor.readline()`.

        Returns
        -------
        None

        Notes
        -----
        Each line of interest contains numeric values, starting from the second
        value up to the seventh one, inclusive. These values are converted to floats,
        scaled down by a factor of 10 (to convert to GPa), assembled into a 2D array (list of lists),
        converted into a NumPy array, scaled, and then stored as a list in the
        `self.elastic_constants` attribute.

        """        
        # Read the total elastic constants
        elastic_constants = []
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        line = self.file_descriptor.readline()
        elastic_constants.append([myfloat(f) for f in line.split()[1:7]])
        econs = np.array(elastic_constants)
        # convert from kBar to GPa
        econs = econs / 10.0
        self.elastic_constants = econs.tolist()
        return

    def _read_ionic_dielectric(self, line):
        """Read ionic dielectric constants from a file and update the zero-frequency static dielectric constant.

        This method assumes the presence of two successive lines in the file corresponding to two dielectric constant values. If any of these lines are marked with an asterisk (*) as missing or do not contain at least three numeric values, default values are used instead. The method calculates the sum of the existing zero-frequency optical dielectric constants and the newly read ionic dielectric constants and updates the corresponding class attribute with this sum.

        Parameters
        ----------
        line : str
            The initial line from which to start reading the ionic dielectric constants. This parameter is not directly used as the method starts by reading a new line from the file, but it's implicitly assumed to represent the current position in the file.

        Returns
        -------
        None

        Notes
        -----
        - Requires that `self.file_descriptor` is an open file object positioned at the correct line where dielectric information starts.
        - The method directly modifies `self.zerof_static_dielectric` based on read values and calculated sums.
        - This method performs no direct output but modifies the class state.

        """        
        # Read the ionic contribution to the static dielectric and use it to computet
        # the full static dielectric constant
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        # the is zero frequency ionic contribution to the static dielectric
        ionic_dielectric = []
        if "*" in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        ionic_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if "*" in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        ionic_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if "*" in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        ionic_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        array1 = np.array(self.zerof_optical_dielectric)
        array2 = np.array(ionic_dielectric)
        array3 = array1 + array2
        self.zerof_static_dielectric = array3.tolist()
        return

    def _read_static_dielectric(self, line):
        """Read and parse static dielectric data from a file.

        This method is designed to read lines from a file, which the instance's file_descriptor points to,
        and parse the static dielectric data. It expects the data to be in specific line positions within the file.
        If a line does not contain data in the expected format or contains an asterisk ('*') indicating missing
        or invalid data, a default value is used instead. The parsed data is stored in the instance's
        zerof_optical_dielectric attribute as a list of lists, where each sublist contains three floating-point numbers.

        Parameters
        ----------
        line : str
            The initial line from where to begin parsing. (This parameter is currently not used in the implementation
            and the method reads lines directly using the object's file_descriptor.)

        Returns
        -------
        None

        """        
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        # the is epsilon infinity
        self.zerof_optical_dielectric = []
        if "*" in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        self.zerof_optical_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if "*" in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        self.zerof_optical_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        line = self.file_descriptor.readline()
        if "*" in line or len(line.split()) < 3:
            line = "99999.999 99999.999 99999.999"
        self.zerof_optical_dielectric.append([myfloat(f) for f in line.split()[0:3]])
        return

    def _read_skip4(self, line):
        """Read and skip the next 4 lines from the current position.

        This method reads and discards the next four lines from the file
        associated with the `file_descriptor` attribute of the calling object.

        Parameters
        ----------
        line : not used
            This parameter is not used but kept for function signature purposes.

        Returns
        -------
        None

        """        
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        self.file_descriptor.readline()
        return

    def _read_external_pressure(self, line):
        """Read and parse external pressure from a given line of text, updating instance attributes.

        Parameters
        ----------
        line : str
            A line of text expected to contain pressure data at specific positions.

        Returns
        -------
        None

        Notes
        -----
        This method updates the instance's `pressure` and `_pulay` attributes based on values extracted from the input string. The `pressure` value is obtained from the fourth item in the split string (`line.split()[3]`), and the `_pulay` value is from the ninth item (`line.split()[8]`). Both values are divided by 10.0 to convert from kbar to GPa.

        """        
        # Vasp writes out kbar so convert to GPa
        self.pressure = float(line.split()[3])/10.0
        self.pressures.append(float(line.split()[3])/10.0)
        self._pulay = float(line.split()[8])/10.0
        return

    def _read_pspot(self, line):
        """Read and store pseudopotential information from a line.

        This method parses a line to extract pseudopotential information and stores it in the `_pspots` dictionary attribute of the object. The pseudopotential label is used as the key, and its corresponding value as the value in the `_pspots` dictionary.

        Parameters
        ----------
        line : str
            A string containing the pseudopotential information to be parsed and stored. Expected format: "any_text value label any_other_text", where "value" and "label" are the parts of interest for storage.

        Returns
        -------
        None

        """        
        self._pspots[line.split()[2]] = line.split()[1]
        return

    def _read_array_dimensions(self, line):
        """Read and set the array dimensions from a line in a file.

        This method reads lines from a file pointed to by `self.file_descriptor`,
        extracts dimensions for kpoints, nbands, and nions from these lines, and
        sets the corresponding attributes of the class instance.

        Parameters
        ----------
        line : str
            The initial line read from the file, not used after modification.

        Returns
        -------
        None

        Notes
        -----
        Assumes specific formatting of the lines being read, where the required
        numbers are located at fixed positions in the lines:
        - `kpoints` is the 4th value in the first line read within the function.
        - `nbands` is the 15th value in the same line.
        - `nions` is the 12th value in the next line read.

        """        
        line = self.file_descriptor.readline()
        self.kpoints = int(line.split()[3])
        self.nbands = int(line.split()[14])
        line = self.file_descriptor.readline()
        self.nions = int(line.split()[11])
        return

    def _read_lattice_vectors(self, line):
        """Read and process lattice vectors from a line of a file.

        This method reads the volume and lattice vectors from specified lines of an input file,
        constructs a new `UnitCell` object from the read vectors, and appends it to the unit cells list.
        If there is at least one `UnitCell` already present, it copies fractional coordinates and element
        names from the last `UnitCell` in the list to the new one.

        Parameters
        ----------
        line : str
            The line from the input file containing the volume of the UnitCell.

        Returns
        -------
        None

        Attributes Modified
        -------------------
        volume : float
            The volume of the UnitCell, updated by parsing the input line.
        volumes : list of float
            List of volumes, appended with the current UnitCell's volume.
        unit_cells : list of UnitCell
            List of `UnitCell` objects, appended with the newly constructed `UnitCell`.

        """        
        self.volume = float(line.split()[4])
        self.volumes.append(float(line.split()[4]))
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        avector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        bvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        line = self.file_descriptor.readline()
        cvector = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
        cell = UnitCell(avector, bvector, cvector,units="Angstrom")
        if self.ncells > 0:
            cell.set_fractional_coordinates(self.unit_cells[-1].fractional_coordinates)
            cell.set_element_names(self.species_list)
        self.unit_cells.append(cell)
        self.ncells = len(self.unit_cells)
        return

    def _read_fractional_coordinates(self, line):
        """Read fractional coordinates from a file and set them for the current unit cell.

        This function reads `nions` lines from the `file_descriptor` attribute, extracts
        the first three space-separated values as fractional coordinates, and appends them
        to the `ions` list. It then sets the fractional coordinates and element names for
        the last unit cell in the `unit_cells` list.

        Parameters
        ----------
        line : str
            The line from which to start reading the fractional coordinates. This parameter
            is passed but not used directly in the function; instead, subsequent lines are
            read from the file.

        Returns
        -------
        None

        """        
        ions = []
        for _n in range(self.nions):
            line = self.file_descriptor.readline()
            ions.append([float(f) for f in line.split()[0:3]])
        self.unit_cells[-1].set_fractional_coordinates(ions)
        self.unit_cells[-1].set_element_names(self.species_list)
        return

    def _read_spin(self, line):
        """Parse and set the spin from a given line.

        Parameters
        ----------
        line : str
            The line of text from which the spin value will be extracted.

        Returns
        -------
        None

        Notes
        -----
        This method expects the spin value to be at the third position in the line. The method doesn't return any value but sets the object's `spin` attribute based on the extracted integer value.

        """        
        self.spin = int(line.split()[2])
        return

    def _read_encut(self, line):
        """Parse and set the energy cutoff value from a given line.

        Parameters
        ----------
        line : str
            A string containing the energy cutoff value. The energy cutoff value is the third element in the string when split by whitespace.

        Returns
        -------
        None

        Notes
        -----
        This method does not return any value. It sets the value of `energy_cutoff` attribute of the instance based on the input line.

        """        
        self.energy_cutoff = float(line.split()[2])
        return

    def _read_ediff(self, line):
        """Read and set the energy difference value from a given line string.

        Parameters
        ----------
        line : str
            A string containing the energy difference information, where the
            energy difference value is expected to be at the third position
            (index 2) when split by whitespace.

        Returns
        -------
        None

        """        
        self._ediff = float(line.split()[2])
        return

    def _read_ibrion(self, line):
        """Read and set the IBRION parameter from a given line.

        This method extracts the IBRION value from a formatted line and sets it to the object's _ibrion attribute.

        Parameters
        ----------
        line : str
            String to be processed

        Returns
        -------
        None
            The line of text containing the IBRION value, expected to be at the third position (index 2) when split by whitespace.

        """        
        self._ibrion = int(line.split()[2])
        return

    def _read_potim(self, line):
        """Parse and set the POTIM parameter from a line of text.

        This method extracts the POTIM parameter value from a provided string
        and sets it as a floating-point number within the object.

        Parameters
        ----------
        line : str
            A string typically representing a line from a file containing
            the POTIM parameter. The expected format is: key-value pairs
            separated by spaces with the POTIM value being the third element
            in the sequence (index 2 when zero-indexed).

        Returns
        -------
        None

        """        
        self._potim = float(line.split()[2])
        return

    def _read_nelect(self, line):
        """Read and set the number of electrons from a given line.

        Parameters
        ----------
        line : str
            A string representation possibly containing information about
            the number of electrons.

        Notes
        -----
        This method updates the object's 'electrons' attribute with the number of electrons
        extracted from the given line. It assumes that the number of electrons is indicated
        at the third position (index 2) in the list produced by splitting `line` on whitespace.
        It also assumes that this number can be represented as a float, but converts it to
        an integer before setting the 'electrons' attribute.

        Returns
        -------
        None

        """        
        self.electrons = int(float(line.split()[2]))
        return

    def _read_magnet(self, line):
        """Read magnetization properties from a line and update the object attributes.

        Parameters
        ----------
        line : str
            A string expected to contain information about electrons and,
            optionally, magnetization. The expected format is whitespace-separated
            values where the fourth value is the number of electrons and the sixth

        Returns
        -------
        None
            value, if present, represents the magnetization.

        Notes
        -----
        Modifies the object's `electrons` attribute by setting it to the value parsed
        from the line. If magnetization information is present in the line, it also
        updates the `magnetization` attribute with that value; otherwise, it sets
        `magnetization` to 0.0.

        """        
        self.electrons = float(line.split()[3])
        if len(line.split()) > 5:
            self.magnetization = float(line.split()[5])
        else:
            self.magnetization = 0.0

    def _read_energy(self, line):
        """Read the energy-related values from a file and update the instance attributes with these values.

        Parameters
        ----------
        line : str
            The current line in the file from which to start reading. This parameter is not used in the function as the function reads new lines itself.

        Returns
        -------
        None

        Notes
        -----
        This method directly modifies the instance attributes `final_free_energy`, `final_free_energies`, `final_energy_without_entropy`, and `final_energies_without_entropy` by reading values from the file associated with `file_descriptor`. The relevant energy values can be extracted from the 5th and 4th positions (zero-indexed) of the line text split by spaces, on specific lines read in sequence.

        """        
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.final_free_energy = float(line.split()[4])
        self.final_free_energies.append(float(line.split()[4]))
        line = self.file_descriptor.readline()
        line = self.file_descriptor.readline()
        self.final_energy_without_entropy = float(line.split()[3])
        self.final_energies_without_entropy.append(float(line.split()[3]))
        return

    def _read_xml(self, filename):
        """Process the VASP xml files.

        Parameters
        ----------
        filename : str
            The file name of the xml file to be processed

        Notes
        -----
        The routine only handles the Vasp output file vasprun.xml

        """
        # For each file generate an xml tree to process.
        tree = ET.parse(filename)
        root = tree.getroot()
        if self.debug:
            print(f"_read_xml filename = {filename}",flush=True)
            print(f"_read_xml tree = {tree}",flush=True)
            print(f"_read_xml root = {root}",flush=True)
        self._handle_kpoints_xml(root.findall("kpoints"))
        self._handle_parameters_xml(root.findall("parameters"))
        self._handle_atominfo_xml(root.findall("atominfo"))
        self._handle_structures_xml(root.findall("structure"))
        self._handle_calculations_xml(root.findall("calculation"))

    def _handle_calculations_xml(self,calculations_xml):
        """Process the calculations xml child.

        Parameters
        ----------
        calculations_xml : a list of xml element
            A list of calculations_xml elements

        Set
        ---
        self.nions : int
            The number of atoms
        self.nspecies : int
            The number of atoms
        self.atom_type_list : list of ints
            A list of the atom types in the order they will be given in the coordinate list
        self.zerof_static_dielectric : 3x3 list of floats
            The ionic contribution to the static dielectric 
        self.zerof_optical_dielectric : 3x3 list of floats
            The zero frequency optical dielectric dielectric 
        self.final_free_energies : a list of floats
            A list of the free energies at each geometry optimisation
        self.final_energies_without_entropy : a list of floats
            A list of the energies without entropy at each geometry optimisation
        self.borncharges : a list of floats nionsx3x3
            A list of the energies without entropy at each geometry optimisation
        self.mass_weighted_normal_modes : ndarray floats (3*nions x n*nions)
            The mass weighted normal modes
        self.frequencies : ndarray floats (3*nions)
            The frequencies in cm-1

        """
        dipole = None
        hessian = None
        dielectric_tensor = None
        forces_list = []
        stress_list = []
        unit_cell_list = []
        dipole_list = []
        # Loop around the calculations
        for calcxml in calculations_xml:
            lastscf = calcxml.findall("scstep/energy")[-1]
            dipoles = calcxml.findall("scstep/dipoles")
            lastdipole = dipoles[-1] if dipoles else None
            # Energies from the last step in this calculation
            # Correction because Vasp does not include the DFT energy
            # properly in the calculate/energy entry
            de = (float(lastscf.find('i[@name="e_0_energy"]').text) -
                  float(lastscf.find('i[@name="e_fr_energy"]').text))
            free_energy = float(calcxml.find('energy/i[@name="e_fr_energy"]').text)
            energy = free_energy + de
            self.final_free_energies.append(energy)
            energy = float(calcxml.find('energy/i[@name="e_wo_entrp"]').text)
            self.final_energies_without_entropy.append(energy)
            # Unit cell
            unit_cell = self._get_unit_cell_from_xml(calcxml.find("structure"))
            unit_cell_list.append(unit_cell)
            # Forces
            forces = None
            xml = calcxml.find('varray[@name="forces"]')
            if xml is not None:
                forces = []
                for vxml in xml:
                    forces.append( [float(f) for f in vxml.text.split()] )
            forces_list.append(forces)
            # stresses 
            stress = None
            xml = calcxml.find('varray[@name="stress"]')
            if xml is not None:
                stress = []
                for vxml in xml:
                    stress.append( [float(f)*0.1 for f in vxml.text.split()])
            stress_list.append(stress)
            # dipoles
            if lastdipole is not None:
                xml = lastdipole.find('v[@name="dipole"]')
                if xml is not None:
                    dipole = [ float(f) for f in xml.text.split() ] 
            xml = calcxml.find('dipole/v[@name="dipole"]')
            if xml is not None:
                dipole = [float(f) for f in xml.text.split()]
            dipole_list.append(dipole)
            # DFPT properties (There should only be one of these in a collection of calculations)
            # epsilon ion
            xml = calcxml.find('varray[@name="epsilon_ion"]')
            if xml is not None:
                zerof_static_dielectric = []
                for vxml in xml:
                    zerof_static_dielectric.append( [ float(f) for f in vxml.text.split()] )
                # print('zerof_static_dielectric',zerof_static_dielectric,flush=True)

            # epsilon optical
            xml = calcxml.find('varray[@name="epsilon"]')
            if xml is not None:
                self.zerof_optical_dielectric = []
                for vxml in xml:
                    self.zerof_optical_dielectric.append( [ float(f) for f in vxml.text.split()] )
                # print('self.zerof_optical_dielectric',self.zerof_optical_dielectric,flush=True)

            # hessian tensor
            xml = calcxml.find('dynmat/i[@name="unit"]')
            if xml is not None:
                unit=xml
            else:
                unit=None
            xml = calcxml.find('dynmat/varray[@name="hessian"]')
            if xml is not None:
                hessian = []
                for v in xml.iter():
                    if v.tag == "v":
                        hessian.append( [ float(f) for f in v.text.split() ] )
            # dielectric tensor (epsilon_scf)
            xml = calcxml.find('varray[@name="epsilon_scf"]')
            if xml is not None:
                dielectric_tensor = []
                for vxml in xml:
                    dielectric_tensor.append( [ float(f) for f in vxml.text.split()] )
            # dielectric tensor (dielectric_dft)
            xml = calcxml.find('varray[@name="dielectric_dft"]')
            if xml is not None:
                dielectric_tensor = []
                for vxml in xml:
                    dielectric_tensor.append( [ float(f) for f in vxml.text.split()] )
                # print('dielectric_tensor',dielectric_tensor,flush=True)
            xml = calcxml.find('array[@name="born_charges"]')
            if xml is not None:
                # Born effective charges
                #        The order used with PDielec is:
                #        d2E/dFxdX d2E/dFxdY d2E/dFxdZ        0 3 6
                #        d2E/dFydX d2E/dFydY d2E/dFydZ        1 4 7
                #        d2E/dFzdX d2E/dFzdY d2E/dFzdZ        2 5 8
                self.born_charges = []
                for element_xml in xml[1:]:  # 1. element = dimension
                    charges = []
                    for vxml in element_xml:
                        charges.append( [ float(f) for f in vxml.text.split() ] )
                    self.born_charges.append(charges)
               #  print('born_charges',self.born_charges,flush=True)
        #
        # Convert Hessian units, first to hertz, then to atomic units
        # NB. this won't work with VASP 6
        #
        evtoj = 1.60217733E-19    # Taken from VASP5.4
        amtokg = 1.6605402E-27    # Taken from VASP5.4
        if unit is not None:
            # Convert to hertz
            conversion = -1.0E24
        else:
            # convert to hertz        # Taken from VASP5.4
            conversion = -evtoj /amtokg / (4*np.pi *np.pi) / 1.0E-10**2 
        # convert to hertz to au
        conversion *= hertz**2
        if hessian is not None:
            hessian = conversion * np.array(hessian)
            self.mass_weighted_normal_modes, self.frequencies = calculate_normal_modes_and_frequencies(hessian)
        return

    def _handle_structures_xml(self,structures_xml):
        """Process the structures xml child.

        Parameters
        ----------
        structures_xml : a list of xml element
            A list of structurexml elements

        Set
        ---
        self.ncells : int
            The number of unit cells read so far
        self.volume : float
            The volume of the last unit cell read
        self.unit_cells : list of unit_cells
            A list of the unit cells

        """
        # Create a dictionary of structures
        structures = dict()
        for elem in structures_xml:
            name = elem.attrib.get("name")
            structures[name] = elem
        # Process the initial position
        self.unit_cells.append(self._get_unit_cell_from_xml(structures["initialpos"]))
        # Book keeping
        self.ncells = len(self.unit_cells)
        self.volume = self.unit_cells[-1].getVolume("Angstrom")
        if self.debug:
            print(f"_handle_structure_xml: volume={self.volume}",flush=True)
        return

    def _get_unit_cell_from_xml(self,structure_xml):
        """Read the xml and return a unit cell with the dimensions and positions.

        Parameters
        ----------
        structure_xml : a structure xml element
            A structure xml element

        Returns
        -------
        cell : a unit_cell

        """
        cell = []
        for v in structure_xml.find("crystal/varray[@name='basis']"):
            cell.append([ float(f) for f in v.text.split()])
        positions = []
        for v in structure_xml.find("varray[@name='positions']"):
            positions.append([ float(f) for f in v.text.split()])
        unit_cell = (UnitCell(cell[0], cell[1], cell[2], units="Angstrom"))
        unit_cell.set_fractional_coordinates(positions)
        unit_cell.set_element_names(self.species_list)
        return unit_cell

    def _handle_atominfo_xml(self,atominfos_xml):
        """Process the atminfos xml child.

        Parameters
        ----------
        atominfos_xml : a list of xml element
            A list of atminfo xml elements

        Set
        ---
        self.nions : int
            The number of atoms
        self.nspecies : int
            The number of atoms
        self.atom_type_list : list of ints
            A list of the atom types in the order they will be given in the coordinate list
        self.species : a list strings
            The list of unique species in the calculation. ie the list of atom types
        self.species_list : a list strings
            The list of species for all atoms
        self.ions_per_tpye : a list ints
            The list of atom types, could be used to lookup self.species
        self.masses_per_type : a list floats
            The list of masses for each atom type
        self.masses : a list floats
            The list of masses for each atom in the unit cell

        """
        self.nions = int( atominfos_xml[-1].find("atoms").text )
        self.nspecies = int( atominfos_xml[-1].find("types").text )
        self.atom_type_list = []
        self.ions_per_type = [ 0 for i in range(self.nspecies) ]
        self.species_list = []
        self.species = []
        for entry in atominfos_xml[-1].find("array/[@name='atoms']/set"):
            species = entry[0].text.strip()
            index = int(entry[1].text.strip())
            self.species_list.append(species)
            if species not in self.species:
                self.species.append(species)
            self.atom_type_list.append(index-1)
            self.ions_per_type[index-1] += 1
        self.masses_per_type = []
        for entry in atominfos_xml[-1].findall("array/[@name='atomtypes']/set/rc"):
            self.masses_per_type.append(float( entry.findall("c")[2].text.strip()))
        self.masses = []
        for atom_type in self.atom_type_list:
            self.masses.append(self.masses_per_type[atom_type])
        if self.debug:
            print("_handle_atomicinfo_xml: nspecies", self.nspecies,flush=True)
            print("_handle_atomicinfo_xml: nions", self.nions,flush=True)
            print("_handle_atomicinfo_xml: ions_per_type", self.ions_per_type,flush=True)
            print("_handle_atomicinfo_xml: atom_type_list", self.atom_type_list,flush=True)
            print("_handle_atomicinfo_xml: species", self.species,flush=True)
            print("_handle_atomicinfo_xml: species_list", self.species_list,flush=True)
            print("_handle_atomicinfo_xml: masses", self.masses,flush=True)
            print("_handle_atomicinfo_xml: masses_per_type", self.masses_per_type,flush=True)
        return

    def _handle_parameters_xml(self,parameters_xml):
        """Process the parameters xml child.

        Parameters
        ----------
        parameters_xml : a list of xml element
            A list of parameters xml elements

        Set
        ---
        self.electrons : int
            The number of electrons
        self.spin : int
            The spin of the system

        """
        self.electrons = 0
        parameters_xml=parameters_xml[-1]
        parameters = dict()
        for ele in parameters_xml.iter():
            name,value = read_xml_element(ele)
            if name is not None:
                parameters[name] = value
        self.electrons = int(parameters["NELECT"]+0.0001)
        self.spin = parameters["ISPIN"]
        # self.masses_per_type = parameters["POMASS"] 
        if self.debug:
            print("_handle_parameters_xml: nelect", self.electrons,flush=True)
            print("_handle_parameters_xml: spin", self.spin,flush=True)
        return

    def _handle_kpoints_xml(self,kpoints_xml):
        """Process the kpoints xml child.

        Parameters
        ----------
        kpoints_xml : a list of xml element
            A list of kpoint xml elements

        Set
        ---
        self.kpoint_grid : a list of 3 integers
            The Monkorst Pack grid
        self.kpoints : int
            The number of kpoints

        """
        kpoints_xml=kpoints_xml[-1]
        self.kpoint_grid = [1, 1, 1]
        xml = kpoints_xml.find("generation").find("v/[@name='divisions']")
        self.kpoint_grid = [ int(f) for f in xml.text.split() ]
        if self.debug:
            print("_handle_kpoints_xml: kpoint_grid", self.kpoint_grid ,flush=True)
        # Number of kpoints
        xml = kpoints_xml.find("varray/[@name='kpointlist']")
        vs = xml.findall("v")
        self.kpoints = len(vs)
        if self.debug:
            print("_handle_kpoints_xml: kpoints", self.kpoints,flush=True)
        return
