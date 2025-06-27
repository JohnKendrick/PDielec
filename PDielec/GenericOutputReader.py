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
"""Generic reader of output files. An actual reader should inherit from this class."""
import math
import os
import sys

import numpy as np

from PDielec.Calculator import calculate_normal_modes_and_frequencies, cleanup_symbol
from PDielec.Constants import amu, avogadro_si, wavenumber
from PDielec.IO import pdielec_io
from PDielec.Plotter import print3x3, print_ints, print_reals, print_strings


class GenericOutputReader:
    """Generic reader of output files. Actual reader should inherit from this class.

    Examples of readers with this base class:

    - :class:`~PDielec.AbinitOutputReader`
    - :class:`~PDielec.CastepOutputReader`
    - :class:`~PDielec.CrystalOutputReader`
    - :class:`~PDielec.ExperimentOutputReader`
    - :class:`~PDielec.GulpOutputReader`
    - :class:`~PDielec.PhonopyOutputReader`
    - :class:`~PDielec.QEOutputReader`
    - :class:`~PDielec.VaspOutputReader`

    Parameters
    ----------
    filenames : list
        A list containing file names (strings) to be associated with the instance. These
        are converted into the absolute paths and stored.

    Attributes
    ----------
    _outputfiles : list
        Stores the initial filenames provided during instantiation.

    names : list
        The absolute paths of the files specified in `filenames`. Computed by resolving
        each filename to its absolute path.

    debug : bool
        Indicates if debug mode is on or off. Defaults to False.

    type : str
        The type of the object, initially set as 'Unkown'.

    ncells, nsteps, electrons, spin, nbands, volume, nions, nspecies, geomsteps : int
        Attributes initialized as 0, representing various numerical properties.

    unit_cells, volumes, species, energiesDFT, energiesDFT_disp, final_free_energies, final_energies_without_entropy, volumes, pressures, born_charges, frequencies, mass_weighted_normal_modes, ions_per_type, atom_type_list, masses, masses_per_type, elastic_constants, zerof_optical_dielectric, zerof_static_dielectric : list
        Various lists to store computational results related to the object. Initialized as empty lists.
    
    final_free_energy, final_energy_without_entropy, magnetization, energy_cutoff, pressure : float
        Floating point attributes initialized to represent different scalar quantities. Defaults to 0 or 0.0.

    kpoints : int
        The number of kpoints, initialized to 1.

    kpoint_grid : list
        The kpoint grid dimensions, initialized as [1, 1, 1].

    manage : set
        The set contains phrases which are matched when looking through a file

    iterations : dict
        Holds the values of 'max force' and 'rms force' at each iteration

    file_descriptor : str
        A string to hold file descriptor, initialized as an empty string.

    eckart : bool
        A boolean flag to indicate the use of Eckart conditions, defaults to False.

    hessian_symmetrisation : str
        The method for symmetrization of Hessian, defaults to 'symm'.

    open_filename, open_directory : str
        Strings for tracking open files or directories, initialized as empty strings.

    _old_masses : list
        A list to keep track of old mass values, initialized as empty.

    nomass_hessian : type or None
        Placeholder for hessian without mass contribution, initially None.

    nomass_hessian_has_been_set : bool
        Indicates if the `nomass_hessian` has been provided, defaults to False.

    original_born_charges, CrystalPermittivity, oscillator_strengths, edited_masses : Various types
        Attributes initialized as `None` to represent unassigned computational properties.

    original_born_charges_are_being_used : bool
        Signifies whether original Born charges are in use, defaults to True.

    """

    def __init__(self, filenames):
        """Initialise an instance of the reader.

        Parameters
        ----------
        filenames : list
            A list containing file names (strings) to be associated with the instance. These
            are converted into the absolute paths and stored.

        """        
        # Test to see if all the files needed are there
        if not all( [ os.path.isfile(f) for f in filenames ] ):
            print("Error in calling the generic output reader some files are not present")
            print(filenames)
            return
        self._outputfiles               = filenames
        self.names                      = [os.path.abspath(f) for f in filenames]
        self.debug                      = False
        self.type                       = "Unkown"
        self.ncells                     = 0
        self.unit_cells                 = []
        self.nsteps                     = 0
        self.electrons                  = 0
        self.spin                       = 0
        self.nbands                     = 0
        self.volume                     = 0
        self.volumes                    = []
        self.nions                      = 0
        self.nspecies                   = 0
        self.geomsteps                  = 0
        self.species                    = []       # a list of the unique species
        self.species_list               = []       # a list of the species for all atoms
        self.energiesDFT                = []
        self.energiesDFT_disp           = []
        self.final_free_energy          = 0.0
        self.final_free_energies        = []
        self.final_energy_without_entropy = 0.0
        self.final_energies_without_entropy = []
        self.kpoints                    = 1
        self.kpoint_grid                = [ 1, 1, 1 ]
        self.energy_cutoff              = 0.0
        self.born_charges               = []
        self.manage                     = {}
        self.iterations                 = {}
        self.file_descriptor            = ""
        self.pressure                   = 0
        self.pressures                  = []
        self.magnetization              = 0.0
        # this in epsilon infinity
        self.zerof_optical_dielectric   = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        # this is the zero frequency static dielectric constant
        self.zerof_static_dielectric    = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        self.elastic_constants          = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
        self.frequencies                = []
        self.mass_weighted_normal_modes = []
        self.ions_per_type              = []
        self.atom_type_list             = []
        self.masses                     = []
        self.masses_per_type            = []
        self.program_mass_dictionary    = {}
        self.eckart                     = False
        self.hessian_symmetrisation     = "symm"
        self.open_filename              = ""
        self.open_directory             = ""
        self._old_masses                = []
        self.nomass_hessian             = None
        self.nomass_hessian_has_been_set= False
        self.original_born_charges      = None
        self.original_born_charges_are_being_used = True
        self.CrystalPermittivity        = None
        self.oscillator_strengths       = None
        self.edited_masses              = None
        self.primitive_transformation   = None
        return

    def read_output(self):
        """Interface to the private read output files methods.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self._read_output_files()
        return

    def reset_masses(self):
        """Reset the mass values to the program's defaults.

        This function resets the mass dictionary of an object to the values specified
        by the object's `program_mass_dictionary`. Typically used when needing to revert
        any runtime changes to the object's mass values to their original state.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        If `debug` attribute of the object is True, it prints a message indicating 
        that the masses are being reset. This function does not return any value
        but modifies the object's state by changing its mass dictionary to match
        the `program_mass_dictionary`.

        Examples
        --------
        Assuming an object `obj` of a class that has `reset_masses` method, `debug` attribute,
        and `program_mass_dictionary`:

            >>> obj.debug = True
            >>> obj.reset_masses()
            Re setting mass dictionary to the program values
    
        Make sure that `program_mass_dictionary` is set properly before calling this
        method to avoid setting the masses to an unintended state.

        """        
        mass_dictionary = {}
        if self.debug:
            print("Re setting mass dictionary to the program values")
        if self.program_mass_dictionary:
            self.change_masses(self.program_mass_dictionary,mass_dictionary)

    def getSpecies(self):
        """Return a list of cleaned species symbols.

        This method applies a specified cleanup function to each element of 
        the `species` attribute and returns the resulting list.

        Parameters
        ----------
        None

        Returns
        -------
        list
            A list of cleaned species symbols.

        """        
        return [ cleanup_symbol(el) for el in self.species ]

    def mass_dictionary(self):
        """Generate a dictionary mapping chemical species symbols to their masses.

        Parameters
        ----------
        None

        Returns
        -------
        dictionary : dict
            A dictionary where keys are cleaned up species symbols and values are their corresponding masses.

        Notes
        -----
        The function uses a helper function `cleanup_symbol` to clean up the species symbols before using them as keys in the dictionary. If the instance's `debug` attribute is set to True, the function will print the newly created mass dictionary.

        Examples
        --------
        Assuming an object `molecule` with species `['H2', 'O']`, masses_per_type `[2.016, 15.999]`, and a debug attribute set to True, calling `molecule.mass_dictionary()` would print `{'H2': 2.016, 'O': 15.999}` and return this dictionary.

        """        
        dictionary = {}
        for symbol,mass in zip(self.species,self.masses_per_type):
            # the element name may be appended with a digit or an underscore
            element = cleanup_symbol(symbol)
            dictionary[element] = mass
        if self.debug:
            print("new mass_dictionary", dictionary)
        return dictionary

    def set_edited_masses(self,masses):
        """Set the edited masses attribute if the length matches the original masses attribute.

        Parameters
        ----------
        masses : list
            A list of mass values to be set as the edited masses.

        Returns
        -------
        None

        Notes
        -----
        The method updates the `edited_masses` attribute only if the length of the input list matches the length of the `masses` attribute of the object. If they do not match, an error message is printed, and `edited_masses` is set to `None`.

        Returns
        -------
        None

        """        
        if len(masses) == len(self.masses):
            self.edited_masses = masses
        else:
            print("Error unable to edited_masses")
            self.edited_masses = None
        return

    def change_masses(self, new_masses, mass_dictionary):
        """Change the masses of elements in a output reader.

        This function updates the masses of the elements within the reader according to the new masses provided. If the object has already had its masses edited, those edited masses will be used. Otherwise, it sets up or updates the simulation's internal mass dictionary based on the provided `new_masses` and an optional external `mass_dictionary` for overrides. If debugging is enabled, various diagnostic messages will be printed during the function's execution.

        Parameters
        ----------
        new_masses : dict
            A dictionary where keys are element symbols as strings and values are the new masses for these elements.
        mass_dictionary : dict, optional
            An additional dictionary provided for mass overrides. If any element present in `new_masses` also exists in `mass_dictionary`, the mass from `mass_dictionary` will be used.

        Returns
        -------
        None

        Notes
        -----
        - It assumes that `cleanup_symbol()` is a function used to standardize or validate the element symbols in `self.species`.
        - It alters `self.masses` and `self.masses_per_type` to reflect changes in masses as per `new_masses` and `mass_dictionary`.
        - If `self.edited_masses` is already set to a truthy value, it bypasses the update process and uses these values instead.
        - `self.species` and `self.atom_type_list` are expected to be iterable attributes of the object containing symbols for elements and types of atoms, respectively.
        - Debugging messages are conditionally printed based on the boolean attribute `self.debug`.

        """        
        if self.edited_masses:
            # This is pretty crude!  If the reader has this variable set then we
            # only use the masses stored in the edit_masses list
            # Once this variable has been set then only these masses are used
            if self.debug:
                print("Using the edited masses")
            self.masses = self.edited_masses
            return
        if not self.program_mass_dictionary:
            # We only want to do this once - remember the program masses as a dictionary
            if self.debug:
                print("Setting program mass dictionary")
            for symbol,mass in zip(self.species,self.masses_per_type):
                element = cleanup_symbol(symbol)
                self.program_mass_dictionary[element] = mass
        if self.debug:
            print("changing masses", self.program_mass_dictionary)
        self.masses = []
        self.masses_per_type = []
        for symbol in self.species:
            # the element name may be appended with a digit or an underscore
            element = cleanup_symbol(symbol)
            mass = new_masses[element]
            if element in mass_dictionary:
                mass = mass_dictionary[element]
            self.masses_per_type.append(mass)
        # end for symbol
        self.masses = [ self.masses_per_type[atype] for atype in self.atom_type_list ]
        if self.debug:
            print("new masses", self.masses)
        return

    def print(self):
        """Print information about the reader.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This function prints out detailed information about the reader it is run on.

        """
        # Generic printing of information
        print("")
        print("Summary of information contained in the QM/MM Reader")
        print("")
        print(f"Number of atoms: {self.nions:5d}")
        print("")
        print(f"Number of species: {self.nspecies:5d}")
        print_strings("Species:", self.species)
        print_ints("Number of atoms for each species:", self.ions_per_type)
        print_reals("Mass of each species:", self.masses_per_type,format="{:10.6f}")
        print_ints("Atom type list:", self.atom_type_list)
        print("")
        print(f"Number of kpoints: {self.kpoints:5d}")
        print("")
        print(f"Kpoint grid      : {self.kpoint_grid[0]:5d} {self.kpoint_grid[1]:5d} {self.kpoint_grid[2]:5d}")
        print("")
        print(f"Energy cutoff (eV): {self.energy_cutoff:f}")
        print("")
        print(f"final_free_energy(eV): {self.final_free_energy:f}")
        print("")
        print(f"geomsteps: {self.geomsteps:f}")
        print("")
        print_reals("DFT energies (eV):", self.energiesDFT,format="{:10.8f}")
        print("")
        print_reals("DFT energies including dispersion (eV):", self.energiesDFT_disp,format="{:10.8f}")
        print("")
        print_reals("Volumes:", self.volumes,format="{:10.8f}")
        print("")
        print_reals("Frequencies (cm-1):", self.frequencies)
        print_reals("Masses (amu):", self.masses,format="{:10.6f}")
        for i, charges in enumerate(self.born_charges):
            title = f"Born Charges for Atom {i:d}"
            print3x3(title, charges)
        print3x3("Epsilon inf: ", self.zerof_optical_dielectric)
        #jk print3x3("Unit cell: ", self.unit_cells[-1].lattice)
        self.unit_cells[-1].print()
        print(" ")
        print(f"Volume of cell: {self.volume:f}")
        mtotal = 0.0
        for m in self.masses:
            mtotal = mtotal + m
        print(f"Total mass is: {mtotal:f} g/mol")
        print(f"Density is: {mtotal/(avogadro_si * self.volume * 1.0e-24):f} g/cc")
        print(" ")
        return

    def get_unit_cell(self):
        """Return the last unit cell in the reader.

        The routine adds the current set of masses to the unit cell.
        It returns None if no unit cells have been read.

        Parameters
        ----------
        None

        Returns
        -------
        unitCell object
            The last unit cell read by the reader

        """
        if len(self.unit_cells) == 0:
            return None
        # Access the last unit cell in the reader
        cell = self.unit_cells[-1]
        # Add the current masses in the reader to the unit cell
        cell.set_atomic_masses(self.masses)
        # return the cell
        return cell

    def get_crystal_density(self):
        """Return the crystal density in g/cc.

        The volume is in angstrom^3, the masses are in atomic mass units, the density is in g/cc

        Parameters
        ----------
        None

        Returns
        -------
        float
            Density in g/cc

        """
        mtotal = 0.0
        for m in self.masses:
            mtotal = mtotal + m
        return mtotal/(avogadro_si * self.volume * 1.0e-24)

    def _read_output_files(self):
        """Read the through the output files.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method is over-ridden by the child class

        """
        # Define the search keys to be looked for in the files
        print("Error _read_output_files must be defined by the actual file reader")
        return

    def _read_output_file(self, name):
        """Read through the file 'name' for key words. The keywords are established in _read_output_files.

        Parameters
        ----------
        name : str
            The file name of the output file

        Returns
        -------
        None

        Notes
        -----
        The actual implementation for identifying keywords is handled by the `_read_output_files` method.

        """
        # Check to see if the file exists....
        if not os.path.isfile(name):
            print("Warning file is not present: ", name, file=sys.stderr)
            return
        # Open file and store file name and directory
        self.file_descriptor = pdielec_io(name, "r")
        self.open_filename = name
        self.open_directory = os.path.dirname(name)
        if self.open_directory == "":
            self.open_directory = "."
        # Loop through the contents of the file a line at a time and parse the contents
        line = self.file_descriptor.readline()
        while line != "":
            for k in self.manage:
                if self.manage[k][0].match(line):
                    method   = self.manage[k][1]
                    if self.debug:
                        print(f"_read_output_file({name}): Match found {k}")
                    method(line)
                    break
                # end if
            # end for
            line = self.file_descriptor.readline()
        # end while
        self.file_descriptor.close()
        return

    def _symmetric_orthogonalisation(self, A):
        """Private routine to perform symmetric orthogonalization.

        Parameters
        ----------
        A : ndarray
            matrix to be orthogonalised

        Returns
        -------
        None

        """
        # The matrix A is only approximately orthogonal
        n = np.size(A, 0)
        unity = np.eye(n)
        Ak = A
        for _k in range(3):
            Bk = np.dot(Ak, Ak.T)
            error = unity - Bk
            Ck = np.linalg.inv(unity + Bk)
            Kk = np.dot(error, Ck)
            Ak = np.dot((unity + Kk), Ak)
            error  = np.sum(np.abs(error))
            if self.debug:
                print("Orthogonalisation iteration: ", error)
        # end for k
        return Ak

    def calculate_mass_weighted_normal_modes(self):
        """Calculate the mass weighted normal modes from the hessian.
  
        The hessian itself is constructed from the frequencies and normal modes
        Any changes to the atomic masses is applied
        Eckart conditions are applied if requested

        Parameters
        ----------
        None

        Returns
        -------
        mass_weighted_normal_modes : np.array
           the mass weighted normal modes

        """
        #
        # Reconstruct the massweighted hessian
        # If necessary diagonalise it and extract the frequencies and normal modes
        # Any projection can be performed here
        #
        # First step is to reconstruct the dynamical matrix (D) from the frequencies and the eigenvectors
        # f^2 = UT . D . U
        # and U is a hermitian matrix so U-1 = UT
        # D = (UT)-1 f^2 U-1 = U f UT
        # Construct UT from the normal modes
        if self.debug:
            print("calculate mass weighted normal modes")
        n = np.size(self.mass_weighted_normal_modes, 0)
        m = np.size(self.mass_weighted_normal_modes, 1)*3
        UT = np.zeros((n, m))
        frequencies_a = np.array(self.frequencies) * wavenumber
        if self.debug:
            print("frequencies_a",frequencies_a)
        masses = np.array(self.masses)*amu
        # if the non mass-weighted hasn't been set, set it
        if not self.nomass_hessian_has_been_set:
            if self.debug:
                print("hessian was not set")
            self.nomass_hessian_has_been_set = True
            for imode, mode in enumerate(self.mass_weighted_normal_modes):
                n = 0
                for atom in mode:
                    # in python the first index is the row of the matrix, the second is the column
                    UT[imode, n+0] = atom[0]
                    UT[imode, n+1] = atom[1]
                    UT[imode, n+2] = atom[2]
                    n = n + 3
                # end for atom
            # end for imode
            # convert the frequencies^2 to a real diagonal array
            # Warning we have to make sure the sign is correct here
            # The convention is that if the frequency is negative
            # then it is really imaginary, so the square of the frequency
            # will be negative too.
            frequencies_a = np.array(self.frequencies) * wavenumber
            f2 = np.diag(np.sign(frequencies_a)*np.real(frequencies_a*frequencies_a))
            # The back transformation uses approximately orthogonal (unitary) matrices because of rounding issues on reading vectors
            # So before that lets orthogonalise them
            UT = self._symmetric_orthogonalisation(UT)
            hessian = np.dot(np.dot(UT.T, f2), UT)
            # Make sure the dynamical matrix is real
            hessian = np.real(hessian)
            # We are going to store the non mass-weighted hessian, since the frequencies came from the MM/QM program
            # they are calculated using the program masses
            current_mass_dictionary = self.mass_dictionary()
            if self.debug:
                print("current mass dictionary", current_mass_dictionary)
            # There is a chance that the program_mass_dictionary hasn't been set - if it hasn't use the current masses
            if not self.program_mass_dictionary:
                self.program_mass_dictionary = current_mass_dictionary
            if self.debug:
                print("program mass dictionary", self.program_mass_dictionary)
            self.change_masses(self.program_mass_dictionary, {})
            masses = np.array(self.masses)*amu
            # remove the mass weighting from the hessian and store
            self.nomass_hessian = self._remove_mass_weighting(hessian,masses)
            # finally replace the masses with those set before we did this
            self.change_masses(current_mass_dictionary, {})
            if self.debug:
                print("non mass weighted hessian", self.nomass_hessian[0:4][0])
        # If the masses have been changed then alter the mass weighted hessian here
        masses = np.array(self.masses)*amu
        if self.debug:
            print("masses", masses)
            print("non mass weighted hessian", self.nomass_hessian[0:4][0])
        self.hessian = self._modify_mass_weighting(self.nomass_hessian, masses)
        if self.debug:
            print("mass weighted hessian", self.hessian[0:4][0])
        # Project out the translational modes if requested
        if self.eckart:
            self.hessian = self.project(self.hessian)
            if self.debug:
                print("projected hessian", self.hessian[0:4][0])
        # Find its eigenvalues and eigen vectors
        self.mass_weighted_normal_modes, self.frequencies = calculate_normal_modes_and_frequencies(self.hessian)
        if self.debug:
            print("calculated frequencies", self.frequencies)
            print("mass-weighted normal modes", self.mass_weighted_normal_modes)
        # end for i
        return self.mass_weighted_normal_modes

    def project(self, hessian):
        """Apply projection operators to remove translation.

        Parameters
        ----------
        hessian : np.array
            The given matrix.

        Returns
        -------
        np.array
            The matrix with translational modes projected out.

        """
        new_hessian = np.zeros_like(hessian)
        nmodes = self.nions*3
        unit = np.eye(nmodes)
        p1 = np.zeros(nmodes)
        p2 = np.zeros(nmodes)
        p3 = np.zeros(nmodes)
        for i in range(self.nions):
            p1[i*3+0] = math.sqrt(self.masses[i])
            p2[i*3+1] = math.sqrt(self.masses[i])
            p3[i*3+2] = math.sqrt(self.masses[i])
        # end for i
        # Normalise
        p1 = p1 / math.sqrt(np.dot(p1, p1))
        p2 = p2 / math.sqrt(np.dot(p2, p2))
        p3 = p3 / math.sqrt(np.dot(p3, p3))
        # Form the projection operators
        P1 = unit - np.outer(p1, p1)
        P2 = unit - np.outer(p2, p2)
        P3 = unit - np.outer(p3, p3)
        # Now project out
        new_hessian = np.dot(np.dot(P1.T, hessian), P1)
        new_hessian = np.dot(np.dot(P2.T, new_hessian), P2)
        return np.dot(np.dot(P3.T, new_hessian), P3)

    def _read_till_phrase(self,phrase):
        """Read lines from the current file until a match with phrase is found.

        Once a match is found, return the matching line.

        Parameters
        ----------
        phrase : str
            The phrase to match in each line of the file.

        Returns
        -------
        str
            The line from the file that matches the phrase.

        """
        line = ""
        while not phrase.match(line):
            line = self.file_descriptor.readline()
        return line

    def _dynamical_matrix(self, hessian):
        """Process the dynamical matrix.

        Processes the dynamical matrix by performing several steps: symmetrizing the Hessian, which is a nxn matrix of mass-weighted force constants, projecting out translational modes, diagonalizing the Hessian, and finally storing the frequencies and normal modes.

        Parameters
        ----------
        hessian : nxn array of floats
            The second derivative matrix

        Returns
        -------
        None

        Notes
        -----
        - Hessian is a nxn matrix of the mass weighted force constants.
        - The hessian is symmetrised.
        - Translational modes are projected out.
        - The hessian is diagonalised.
        - Finally, the frequencies and normal modes are stored.

        """
        if self.debug:
            print("_dynamical_matrix")
            print("hessian", hessian[0:4][0])
        masses = np.array(self.masses)*amu
        if self.debug:
            print("masses", self.masses, masses)
        if not self.nomass_hessian_has_been_set:
            # symmetrise the hessian and store it for later use
            if self.hessian_symmetrisation == "symm":
                hessian = 0.5 * (hessian + hessian.T)
            else:
                # Crystal only uses the upper triangle to calculate the eigenvalues
                hessian = hessian - np.tril(hessian) + np.triu(hessian).T
            self.nomass_hessian_has_been_set = True
            self.nomass_hessian = self._remove_mass_weighting(hessian,masses)
        if self.debug:
            print("non mass weighted hessian", self.nomass_hessian[0:4][0])
        hessian = self._modify_mass_weighting(self.nomass_hessian,masses)
        if self.debug:
            print("non mass weighted hessian", self.nomass_hessian[0:4][0])
        if self.debug:
            print("mass weighted hessian", hessian[0:4][0])
        # Project out the translational modes if requested
        if self.eckart:
            hessian = self.project(hessian)
        if self.debug:
            print("projected hessian", hessian[0:4][0])
        self.mass_weighted_normal_modes, self.frequencies = calculate_normal_modes_and_frequencies(hessian)
        if self.debug:
            print("non mass weighted hessian", self.nomass_hessian[0:4][0])
        return

    def reset_born_charges(self):
        """Reset the born charges to their original values if they are not currently being used.

        This method restores the `born_charges` attribute of the instance to the value stored in `original_born_charges`, but only if `original_born_charges_are_being_used` is `False`.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method changes the state of the reader by modifying its `born_charges` attribute to match `original_born_charges`, under the condition that `original_born_charges_are_being_used` is `False`.

        """        
        if not self.original_born_charges_are_being_used:
            self.born_charges = self.original_born_charges

    def neutralise_born_charges(self):
        """Neutralise Born charges within the object.

        Changes the state of `original_born_charges_are_being_used` 
        to False and saves the current `born_charges` as `original_born_charges` if
        `original_born_charges_are_being_used` is True, indicating that the 
        original Born charges are no longer being used directly. It then applies 
        the Born charge sum rule by calling the `_born_charge_sum_rule` method.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method is intended to be used within a context where Born charges 
        (representative of the polarization of ions in a solid under an electric field) 
        need to be neutralized or altered from their original state.

        """        
        if self.original_born_charges_are_being_used:
            self.original_born_charges = self.born_charges
            self.original_born_charges_are_being_used = False
        self._born_charge_sum_rule()
        return

    def _born_charge_sum_rule(self):
        """Apply a simple charge sum rule to all the elements of the born matrices.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This function applies a straightforward summation rule to the elements of the born matrices,

        """
        total = np.zeros((3, 3))
        born_charges = np.array(self.original_born_charges)
        new_born_charges = np.zeros_like(self.born_charges)
        total = np.sum(born_charges) / self.nions
        if self.debug:
            print("born charge sum", total)
        new_born_charges = born_charges - total
        self.born_charges = new_born_charges.tolist()
        return

    def _modify_mass_weighting(self,hessian,new):
        """Modify the Hessian matrix based on new mass weighting.

        This function iterates over the elements of the Hessian matrix (`hessian`) and adjusts each element based on the square root of the product of elements from a new weighting (`new`). 

        Parameters
        ----------
        hessian : numpy.ndarray
            The original Hessian matrix whose elements are to be modified. It is assumed to be a square matrix with dimensions [3*nions, 3*nions], where nions is the number of ions (or atoms) in the system.
        new : numpy.ndarray or list
            An array or list of new weighting factors, one per ion. Its length should be equal to `nions`. These factors are applied in the modification of the Hessian matrix elements.

        Returns
        -------
        numpy.ndarray
            The modified Hessian matrix with the same dimensions as the input matrix.

        """        
        ipos = -1
        new_hessian = np.empty_like(hessian)
        for i in range(self.nions):
            for _ix in range(3):
                ipos += 1
                jpos = -1
                for j in range(self.nions):
                    for _jx in range(3):
                        jpos += 1
                        new_hessian[ipos,jpos] = hessian[ipos,jpos] / math.sqrt(  new[i]*new[j] )
                    # end for jx
                # end for j
            # end for ix
        # end for i
        return new_hessian

    def _remove_mass_weighting(self,hessian,old):
        """Remove mass-weighting from a Hessian matrix.

        Parameters
        ----------
        hessian : ndarray
            The Hessian matrix to be un-weighted, typically a square, 2D array where
            the dimension correlates with 3 times the number of ions (nions*3) since
            it includes the x, y, and z coordinates for each ion.
        old : ndarray or list
            An array or list of mass weights for each ion.

        Returns
        -------
        new_hessian : ndarray
            The mass-weighting removed (un-weighted) Hessian matrix.

        Notes
        -----
        This function assumes that the Hessian matrix is provided in a form that
        has been mass-weighted. Each element of the `hessian` input is divided by
        the square root of the product of the mass weights of the corresponding ions,
        adjusting it to an un-weighted form.

        The size of `old` must match the number of ions (`nions`). Each entry in `old`
        is used to compute the square root of the product of mass weights for the 
        appropriate matrix element, thereby removing the mass weighting from the 
        original Hessian matrix.

        """        
        new_hessian = np.empty_like(hessian)
        ipos = -1
        for i in range(self.nions):
            for _ix in range(3):
                ipos += 1
                jpos = -1
                for j in range(self.nions):
                    for _jx in range(3):
                        jpos += 1
                        new_hessian[ipos,jpos] = hessian[ipos,jpos] * math.sqrt( old[i]*old[j] )
                    # end for jx
                # end for j
            # end for ix
        # end for i
        return new_hessian
