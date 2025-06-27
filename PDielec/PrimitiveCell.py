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
"""The PrimitiveCell module."""

import numpy as np
import spglib

from PDielec.UnitCell import UnitCell


class PrimitiveCell(UnitCell):
    """Hold primitive-cell information; the unit cell and the primitive-cell size.
    
    Inherits from UnitCell, but adds a new variable, the original cell
    The initialisation is always performed with a unit cell.
    A transformation of the original cell can be applied using the transformaton matrix

    Parameters
    ----------
    aUnitCell : :class:`~PDielec.UnitCell.UnitCell`
        The unit cell associated with this instance.

    Attributes
    ----------
    unitCell : :class:`~PDielec.UnitCell.UnitCell`
        The primitive cell associated with this instance. (see :class:`~PDielec.UnitCell.UnitCell`)
    originalCell : :class:`~PDielec.UnitCell.UnitCell`
        The original cell associated with this instance. (see :class:`~PDielec.UnitCell.UnitCell`)
    method   : str
        A string with either "auto" or "given"
    transformation : 3x3 list
        A transformation matrix to create the primitive cell

    See Also
    --------
    :class:`~PDielec.UnitCell.UnitCell`

    """

    def __init__(self, aUnitCell, transformation=None):
        """Initialize the primitive cell with a unit cell.

        Parameters
        ----------
        aUnitCell : :class:`~PDielec.UnitCell.UnitCell`
            The unit cell associated with this instance.
        transformation : 3x3 list (optional)
             This is an optional parameter giving the 3x3 transformation matrix
             If it is None then the primitive cell is guessed using spglib

        Attributes
        ----------
        unitCell : :class:`~PDielec.UnitCell.UnitCell`
            The primitive cell associated with this instance.
        transformation : 3x3 list (optional)
            A transformation matrix to create the primitive cell

        """
        #
        # Initialise everything using the parent class
        #
        super().__init__()
        #
        # Keep a copy of the original cell
        #
        self.originalCell   = aUnitCell         # the unit cell on which everything is based
        # 
        # If there are any atoms outside the cell boundaries move them inside
        # NB. this undoes whole molecules
        #
        if isinstance(transformation,np.ndarray):
            #
            # Store the transformation
            #
            self.transformation = transformation
        elif isinstance(transformation,list):
            self.transformation = np.array(transformation)
        elif isinstance(transformation,str):
            #
            # Create an spglib cell and analyse the symmetry
            #
            spglib_cell = (aUnitCell.lattice, 
                           aUnitCell.fractional_coordinates, 
                           aUnitCell.get_atomic_numbers())
            (cell, coords, atnos) = spglib.find_primitive(spglib_cell)
            #
            # The transformation is defined by a mapping between the primitive and standard cell
            # In this routine we assume that the standard cell is given by self.originalCell
            #
            self.transformation = np.dot(cell, np.linalg.inv(self.originalCell.lattice))
            #
            # Create a UnitCell with the primitive cell dimension
            # but with the contents from the full cell
        else:
            self.transformation = np.eye(3)
        #
        # Now we are using self to store the primitive cell
        #
        new_lattice = np.dot(self.transformation.T, aUnitCell.lattice)
        self.set_lattice(new_lattice)
        #
        # Add the xyz coordinates and therefore calculate the fractional ones
        #
        self.set_xyz_coordinates(aUnitCell.xyz_coordinates)
        self.set_atomic_numbers(aUnitCell.get_atomic_numbers())
        self.set_atom_labels(aUnitCell.get_atom_labels())
        self.set_atomic_masses(aUnitCell.get_atomic_masses())
        self.set_element_names(aUnitCell.get_element_names())
        #jk print('JK1 ---------------------------------------')
        #jk print(new_lattice)
        #jk self.print()
        #jk print('JK1 ---------------------------------------')
        #
        # Trim this new cell so that any atoms outside the cell boundaries are removed
        # Atoms, labels, masses and element names are removed if the atom is not unique
        # The map lists help remember the relationship between everything.
        #
        #self.fold_cell()
        self.map_old_to_new, self.map_new_to_old = self.fold_cell()
        #jk print('JK2 ----After trim-------------------------')
        #jk self.print()
        #jk print('JK2 ---------------------------------------')
