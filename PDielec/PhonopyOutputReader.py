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
"""Read the contents of a directory containing Phonopy input and output files.
"""

import numpy as np

from PDielec.GenericOutputReader import GenericOutputReader


class PhonopyOutputReader(GenericOutputReader):
    """Read the contents of a directory containing Phonopy input and output files.

    Inherits from :class:`~PDielec.GenericOutputReader.GenericOutputReader`

    Parameters
    ----------
    names : list
        A list of file names to be used.
    qmreader : object
        An instance of another class used for reading.

    Notes
    -----
    This method sets the type of the output to 'Phonopy output' and 
    associates a quantum mechanics (QM) reader object with the instance.

    """

    def __init__(self, names, qmreader):
        """Constructor for creating a new instance of the class.

        Parameters
        ----------
        names : list
            A list of file names to be used.
        qmreader : object
            An instance of another class used for reading.

        Notes
        -----
        This method sets the type of the output to 'Phonopy output' and 
        associates a quantum mechanics (QM) reader object with the instance.

        """        
        GenericOutputReader.__init__(self, names)
        # We have to use the qm reader to do the reading of the QM files
        self.type                    = 'Phonopy output'
        self.qmreader                = qmreader
        return

    def _read_output_files(self):
        """Read the Phonopy files in the directory.

        Uses the qmreader to read in the electronic structure information and copies it to the current object.
        Then it reads the dynamical matrix from phonopy
        """
        # Set the qm reader to have all the settings that phonopy reader has
        self.qmreader.eckart = self.eckart
        self.qmreader.debug = self.debug
        # trigger the reading of the qm files
        self.qmreader.read_output()
        # We don't call self._read_outputfile as this starts looking for keywords
        # for f in self._outputfiles:
        #     self._read_output_file(f)
        #
        # Instead copy anything useful from the QM calcs to self
        self.ncells                  = self.qmreader.ncells
        self.unit_cells              = self.qmreader.unit_cells
        self.volume                  = self.qmreader.volume
        self.spin                    = self.qmreader.spin
        self.energy_cutoff           = self.qmreader.energy_cutoff
        self.kpoints                 = self.qmreader.kpoints
        self.kpoint_grid             = self.qmreader.kpoint_grid
        self.nbands                  = self.qmreader.nbands
        self.nions                   = self.qmreader.nions
        self.ions_per_type           = self.qmreader.ions_per_type
        self.atom_type_list          = self.qmreader.atom_type_list
        self.electrons               = self.qmreader.electrons
        self.magnetization           = self.qmreader.magnetization
        self.final_energy_without_entropy = self.qmreader.final_energy_without_entropy
        self.final_free_energy       = self.qmreader.final_free_energy
        self.pressure                = self.qmreader.pressure
        self.masses_per_type         = self.qmreader.masses_per_type
        self.ions_per_type           = self.qmreader.ions_per_type
        self.masses                  = self.qmreader.masses
        self.nspecies                = self.qmreader.nspecies
        self.species                 = self.qmreader.getSpecies()
        self.born_charges            = self.qmreader.born_charges
        self.zerof_optical_dielectric= self.qmreader.zerof_optical_dielectric
        self.zerof_static_dielectric = self.qmreader.zerof_static_dielectric
        # Calculate dynamical matrix
        self.read_dynamical_matrix()
        return

    def read_dynamical_matrix(self):
        """Read and process the dynamical matrix from output files.

        This method reads the dynamical matrix and other relevant data from the specified
        output files, calculates the Hessian matrix, converts its units, and computes the
        eigenvalues and eigenvectors to determine frequencies and mass-weighted normal modes.

        Parameters
        ----------
        None explicitly, but relies on the instance attributes `_outputfiles` and `nions` 
        to determine the files to read and number of ions, respectively.

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
        except:
            print("WARNING: Yaml CLoader is not avaiable, using fallback")
            from yaml import Loader as Loader
        # the first name has to be the qpoints file
        fd = open(self._outputfiles[0])
        data_q = yaml.load(fd, Loader=Loader)
        fd.close
        # the second name has to be the phonopy file
        fd = open(self._outputfiles[1])
        data_p = yaml.load(fd, Loader=Loader)
        fd.close
        self._old_masses = []
        for i in range(self.nions):
            self._old_masses.append(data_p['primitive_cell']['points'][i]['mass'])
        #qpoints = data_q['phonon'][0]['q-position']
        # print('q-points',qpoints)
        #natom = data_q['natom']
        # print('natom:',natom)
        dynmat = []
        dynmat_data = data_q['phonon'][0]['dynamical_matrix']
        for row in dynmat_data:
            vals = np.reshape(row, (-1, 2))
            dynmat.append(vals[:, 0] + vals[:, 1] * 1j)
        dynmat = np.array(dynmat)
        # Make sure the hessian is real
        hessian = np.real(dynmat)
        # We need to convert to cm-1
        conversion_factor_to_THz = 15.633302
        conversion_factor_to_cm1 = conversion_factor_to_THz * 33.35641
        conv  = conversion_factor_to_cm1
        hessian = hessian * conv * conv
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
            for j in range(self.nions):
                modea = [eig_vec[n][i], eig_vec[n+1][i], eig_vec[n+2][i]]
                n = n + 3
                mode.append(modea)
            self.mass_weighted_normal_modes.append(mode)
        # end for i
        return

