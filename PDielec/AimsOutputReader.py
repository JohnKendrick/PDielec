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
"""Read the contents of a directory containing FHI-AIMS input and output files."""

import os

import numpy as np

from PDielec import Constants
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell


class AimsOutputReader(GenericOutputReader):
    """Read the contents of a directory containing Aims input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    names : list
        A list of file names to be used.

    """

    def __init__(self, names):
        """Initialise a Aims output reader.

        Parameters
        ----------
        names : list
            A list of file names to be used.

        This routine reads files with fixed names:
                  geometry.in         contains details of the cells used
                  hessian.???.yaml    contains the dynamical matrix
                  BORN_PDIELEC        generated using aims-pdielec-born

        """        
        GenericOutputReader.__init__(self, names)
        self.type                    = "Aims output"
        self.units                   = { "length"          : "angstrom",
                                         "mass"            : "AMU",
                                         "force"           : None,
                                         "force_constants" : None, }
        return

    def _read_output_files(self):
        """Read the Phonopy files in the directory.

        reads the dynamical matrix from FHI-Aims
        reads the Born charges
        """
        # Calculate dynamical matrix
        energy_filename = self._outputfiles[0]
        masses_filename = self._outputfiles[1]
        geometry_filename = self._outputfiles[2]
        hessian_filename  = self._outputfiles[3]
        born_filename     = self._outputfiles[4]
        if os.path.exists(energy_filename):
            self.read_energy(energy_filename)
        if os.path.exists(masses_filename) and os.path.exists(geometry_filename):
            self.read_geometry(masses_filename,geometry_filename)
        if os.path.exists(hessian_filename):
            self.read_dynamical_matrix(hessian_filename)
        if os.path.exists(born_filename):
            self.read_born_file(self.nions,born_filename)
        return

    def read_energy(self, file):
        """Read and process the aims log file.

        Parameters
        ----------
        file : str
            The aims log file, usually aims.out

        Returns
        -------
        None

        Modifies
        --------
        self.final_energy_without_entropy float
        self.final_energies_without_entropy float

        """        
        # Read the aims.out file
        with open(file) as fd:
            for line in fd:
                if "  | Total energy of the DFT " in line:
                    self.final_energy_without_entropy = float(line.split()[-2])
                    self.final_energies_without_entropy.append(self.final_energy_without_entropy)
        return

    def read_geometry(self, masses_file, geometry_file):
        """Read and process the FHI-AIMS masses.*.dat and geometry.in file.

        Parameters
        ----------
        masses_file : str
            The filename containing the FHI-Aims masses, usually masses.*.dat
        geometry_file : str
            The filename containing the FHI-Aims geometry, usually geometry.in

        Returns
        -------
        None

        Modifies
        --------
        self.unit_cell as a list with a single unit cell in it

        """        
        # Read the masses
        masses = []
        with open(masses_file) as fd:
            for line in fd:
                masses.append(float(line.split()[0]))
        # Read the geometry
        with open(geometry_file) as fd:
            lines = fd.readlines()
        lattice = []
        coordinates = []
        symbols = []
        for line in lines:
            split_line = line.split()[1:]
            if "lattice_vector" in line:
                lattice.append([ float(f) for f in split_line[:3] ])
            elif "atom_frac" in line:
                cartesian = False
                coordinates.append([ float(f) for f in split_line[:3] ])
                symbols.append(split_line[3])
            elif "atom" in line:
                cartesian = True
                coordinates.append([ float(f) for f in split_line[:3] ])
                symbols.append(split_line[3])
        cell = UnitCell(lattice[0], lattice[1], lattice[2], units=self.units["length"])
        if cartesian:
            cell.set_xyz_coordinates(coordinates)
        else:
            cell.set_fractional_coordinates(coordinates)
        cell.set_element_names(symbols)
        cell.set_atomic_masses(masses)
        self.unit_cell = [ cell ]
        return

    def read_dynamical_matrix(self,filename):
        """Read and process the dynamical matrix from output files.

        This method reads the dynamical matrix from the specified
        output files calculates the Hessian matrix, converts its units, and computes the
        eigenvalues and eigenvectors to determine frequencies and mass-weighted normal modes.

        Parameters
        ----------
        filename : str
            The filename containing the hessian, usually of the form hessian.*.dat

        Returns
        -------
        None
            This method directly modifies the instance attributes `_old_masses`, `frequencies`, 
            and `mass_weighted_normal_modes` based on the data read and calculations performed.

        Raises
        ------
        FileNotFoundError
            If any of the specified output files cannot be found or opened.
        MemoryError
            If there is insufficient memory to compute the eigendecomposition.

        Notes
        -----
        - The dynamical matrix is processed to compute real-valued Hessian by considering only the real
          parts. It is then converted to the desired units before computing the eigenvalues and eigenvectors.
        - The method computes frequencies in THz by first converting the eigenvalues from the Hessian to
          the correct units, then calculates the square root of their absolute values, preserving the sign
          to distinguish between real and imaginary modes.
        - Mass-weighted normal modes are then computed for each ion based on the eigenvectors, and the 
          frequencies and normal modes are stored in corresponding instance attributes.

        """        
        hessian = np.loadtxt(filename)
        cell = self.unit_cell[-1]
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
        # Values of units in SI
        ev = 1.0/Constants.hartree2ev
        amu = Constants.amu
        angs = Constants.angs2bohr
        wavenumber = Constants.wavenumber
        #
        # Process the hessian
        # Mass weighting
        #
        mass_vector = np.zeros([0])
        for mass in self._old_masses:
            mass_vector = np.append(mass_vector, np.ones(3) * 1.0/np.sqrt(mass))
        mass_mat=mass_vector[:,np.newaxis]*mass_vector[np.newaxis,:]  #1/sqrt(mimj)
        hessian = hessian * mass_mat
        #
        # The FHI-Aims units for the hessian are ev/(amu*Ang*Ang)
        # PDielec works in atomic units 
        #
        conversion_factor = ev / (amu * angs * angs) 
        hessian *= conversion_factor
        # Symmetrise hessian
        hessian = (hessian.T + hessian)/2.0
        # Find its eigenvalues and eigen vectors
        eig_val, eig_vec = np.linalg.eigh(hessian)
        self.mass_weighted_normal_modes = []
        nmodes = 3*self.nions
        # Store the new frequencies, using the negative convention for imaginary modes
        frequencies_a = np.sqrt(np.abs(eig_val.real)) * np.sign(eig_val.real)
        self.frequencies = (frequencies_a/wavenumber).tolist()
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
