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
"""Read the contents of a directory containing Phonopy input and output files."""

import numpy as np

from PDielec.Constants import thz2cm1
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell


class PhonopyOutputReader(GenericOutputReader):
    """Read the contents of a directory containing Phonopy input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    names : list
        A list of file names to be used.

    """

    def __init__(self, names):
        """Initialise a Phonopy output reader.

        Parameters
        ----------
        names : list
            A list of file names to be used.

        This routine reads files with fixed names:
                  phonopy.yaml    contains details of the cells used
                  qpoints.yaml    contains the dynamical matrix
                  BORN_PDIELEC    generated from the BORN file using phonopy-pdielec-born

        """        
        GenericOutputReader.__init__(self, names)
        self.type                    = "Phonopy output"
        self.units                   = { "length"          : "angstrom",
                                         "mass"            : "AMU",
                                         "force"           : None,
                                         "force_constants" : None, }
        return

    def _read_output_files(self):
        """Read the Phonopy files in the directory.

        reads the dynamical matrix from phonopy
        reads the Born charges
        """
        # Calculate dynamical matrix
        qpoints_filename = self._outputfiles[0]
        phonopy_filename = self._outputfiles[1]
        born_filename    = self._outputfiles[2]
        self.read_dynamical_matrix(phonopy_filename,qpoints_filename)
        self.read_born_file(self.nions,born_filename)
        return

    def read_dynamical_matrix(self,phonopy_filename, qpoints_filename):
        """Read and process the dynamical matrix from output files.

        This method reads the dynamical matrix and other relevant data from the specified
        output files, calculates the Hessian matrix, converts its units, and computes the
        eigenvalues and eigenvectors to determine frequencies and mass-weighted normal modes.

        Parameters
        ----------
        phonopy_filename : str
            The phonopy.yaml file name
        qpoints_filename : str
            The qpoints.yaml file name

        Returns
        -------
        None
            This method directly modifies the instance attributes `_old_masses`, `frequencies`, 
            and `mass_weighted_normal_modes` based on the data read and calculations performed.

        Raises
        ------
        FileNotFoundError
            If any of the specified output files cannot be found or opened.
        YAMLError
            If there is an error parsing the YAML files.
        MemoryError
            If there is insufficient memory to compute the eigendecomposition.

        Notes
        -----
        - It first attempts to use the CLoader for parsing YAML files for performance reasons and falls
          back to the standard Loader if CLoader is not available, issuing a warning in such cases.
        - The dynamical matrix is processed to compute real-valued Hessian by considering only the real
          parts. It is then converted to the desired units before computing the eigenvalues and eigenvectors.
        - The method computes frequencies in THz by first converting the eigenvalues from the Hessian to
          the correct units, then calculates the square root of their absolute values, preserving the sign
          to distinguish between real and imaginary modes.
        - Mass-weighted normal modes are then computed for each ion based on the eigenvectors, and the 
          frequencies and normal modes are stored in corresponding instance attributes.

        """        
        #
        # Yaml imports of large files are really slow....
        # Attempt to use the PyYaml C parser, using yaml.CLoader
        #
        import yaml
        try:
            from yaml import CLoader as Loader
        except Exception:
            print("WARNING: Yaml CLoader is not avaiable, using fallback")
            from yaml import Loader as Loader
        # the first name has to be the qpoints file
        with open(qpoints_filename) as fd:
            data_q = yaml.load(fd, Loader=Loader)
        # the second name has to be the phonopy file
        with open(phonopy_filename) as fd:
            data_p = yaml.load(fd, Loader=Loader)
        #
        # process phonopy.yaml
        #
        try:
            conversion_factor_to_THz = data_p["phonopy"]["frequency_unit_conversion_factor"]
        except Exception:
            conversion_factor_to_THz = 15.633302
        conversion_factor_to_cm1 = conversion_factor_to_THz * thz2cm1
        #
        # Get units
        #
        if "physical_unit" in data_p:
            if "atomic_mass" in data_p["physical_unit"]:
                self.units["atomic_mass"]     = data_p["physical_unit"]["atomic_mass"]
            if "length" in data_p["physical_unit"]:
                self.units["length"]          = data_p["physical_unit"]["length"]
            if "force" in data_p["physical_unit"]:
                self.units["force"]           = data_p["physical_unit"]["force"]
            if "force_constants" in data_p["physical_unit"]:
                self.units["force_constants"] = data_p["physical_unit"]["force_constants"]
        #
        # Determine axes and cells
        #
        if "primitive_axes" in data_p["phonopy"]["configuration"]:
            primitive_axes = data_p["phonopy"]["configuration"]["primitive_axes"]
        else:
            primitive_axes = None
        if "primitive_matrix" in data_p:
            self.primitive_transformation = data_p["primitive_matrix"]
        else:
            self.primitive_transformation = None
        primitive_cell = self.read_cell(data_p["primitive_cell"]) if "primitive_cell" in data_p else None
        unit_cell = self.read_cell(data_p["unit_cell"]) if "unit_cell" in data_p else None
        #
        # Use the cell that is consistent with the primitive_axes
        #
        cell = primitive_cell
        self._old_masses = cell.get_atomic_masses()
        self.nions = len(self._old_masses)
        self.species = cell.get_species()
        self.nspecies = len(self.species)
        elements = cell.get_element_names()
        # Calculate the number of ions / type
        # Determine the atom type list
        self.ions_per_type = [ 0 for _species in self.species ]
        for element in elements:
            species_index = self.species.index(element)
            self.ions_per_type[species_index] += 1
            self.atom_type_list.append(species_index)
        # Determine the mass of each type
        self.mass_per_type = []
        for species in self.species:
             element_index = elements.index(species)
             self.masses_per_type.append(self._old_masses[element_index])
        # Proceed
        self.unit_cells = [ cell ]
        self.ncells = 1
        self.volume = cell.getVolume("Angstrom")
        #
        # Process qpoints.yaml
        #
        dynmat = []
        dynmat_data = data_q["phonon"][0]["dynamical_matrix"]
        for row in dynmat_data:
            vals = np.reshape(row, (-1, 2))
            dynmat.append(vals[:, 0] + vals[:, 1] * 1j)
        dynmat = np.array(dynmat)
        # Make sure the hessian is real
        hessian = np.real(dynmat)
        # We need to convert to cm-1
        hessian = hessian * conversion_factor_to_cm1 * conversion_factor_to_cm1
        # Find its eigenvalues and eigen vectors
        eig_val, eig_vec = np.linalg.eigh(hessian)
        self.mass_weighted_normal_modes = []
        nmodes = 3*self.nions
        # Store the new frequencies, using the negative convention for imaginary modes
        frequencies_a = np.sqrt(np.abs(eig_val.real)) * np.sign(eig_val.real)
        self.frequencies = frequencies_a.tolist()
        # Store the mass weighted normal modes
        for i in range(nmodes):
            mode = []
            n = 0
            for _j in range(self.nions):
                modea = [eig_vec[n][i], eig_vec[n+1][i], eig_vec[n+2][i]]
                n = n + 3
                mode.append(modea)
            self.mass_weighted_normal_modes.append(mode)
        # end for i
        return

    def read_cell(self, yaml):
        """Read and process yaml string and extract a unit cell.

        Parameters
        ----------
        yaml : str
            A yaml string containing a phonopy cell

        Returns
        -------
        cell : unitCell

        """        
        lattice = yaml["lattice"]
        points = yaml["points"]
        symbols = [ atom["symbol"] for atom in points ]
        coordinates = [ atom["coordinates"] for atom in points ]
        masses = [ atom["mass"] for atom in points ]
        cell = UnitCell(lattice[0], lattice[1], lattice[2], units=self.units["length"])
        cell.set_fractional_coordinates(coordinates)
        cell.set_element_names(symbols)
        cell.set_atomic_masses(masses)
        return cell

    def read_born_file(self,natoms,filename):
        """Read the BORN_PDIELEC file in the current directory.

        The BORN_PDIELEC file was created by a phonopy helper routine such as phonopy-pdielec-born
        It contains the optical permittivity and the Born charges

        Parameters
        ----------
        natoms : int
            The number of atoms
        filename : str
            The filename (probably "BORN_PDIELEC")

        Modifies
        --------
        zerof_optical_dielectric : the zero frequency optical permittivity
        born_charges             : the born charges

        """        
        with open(filename) as fd:
            #
            # Skip a line
            #
            line = fd.readline()
            #
            # Read dielectric constant
            #
            line = fd.readline().split()
            if len(line) != 9:
                print("BORN file format of line 2 is incorrect")
                return
            self.zerof_optical_dielectric = np.reshape([float(x) for x in line], (3, 3))
            #
            # Read Born effective charge
            #
            self.born_charges = np.zeros((natoms, 3, 3), dtype="double")
            for i in range(natoms):
                line = fd.readline().split()
                if len(line) == 0:
                    print("Number of lines for Born effect charge is not enough.")
                    return
                if len(line) != 9:
                    print("BORN file format of line %d is incorrect" % (i + 3))
                    return
                self.born_charges[i] = np.reshape([float(x) for x in line], (3, 3))
        return
