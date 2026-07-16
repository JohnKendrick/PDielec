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
"""Finite-field Raman JSON output reader."""

import json
import logging
import math

import numpy as np

from PDielec.Calculator import calculate_normal_modes_and_frequencies, cleanup_symbol
from PDielec.Constants import amu, angs2bohr, average_masses, hartree2ev
from PDielec.GenericOutputReader import GenericOutputReader
from PDielec.UnitCell import UnitCell

logger = logging.getLogger(__name__)


class FiniteFieldOutputReader(GenericOutputReader):
    """Read PDielec finite-field Raman JSON datasets.

    The reader is code-agnostic: VASP, QE, CASTEP, or other finite-field
    workflows can write the same consolidated JSON schema.  Program/DFT masses
    are retained as ``program_mass_dictionary``, but the active masses default
    to PDielec's built-in average masses.
    """

    _SUPPORTED_SCHEMAS = {
        "pdielect-finite-field-raman-v1",
        "pdielect-r_epsilon-v1",
        "pdielect-vASP-consolidated-raman-results-v1",
    }

    def __init__(self, names):
        """Initialise the finite-field JSON reader."""
        GenericOutputReader.__init__(self, names)
        self.type = "Finite field Raman JSON"
        self._deps_dr = None
        return

    def _read_output_files(self):
        """Read the finite-field JSON file."""
        with open(self._outputfiles[0]) as fd:
            data = json.load(fd)

        schema = data.get("schema", "")
        if schema not in self._SUPPORTED_SCHEMAS:
            logger.warning(f"FiniteFieldOutputReader: unknown schema {schema}; attempting to read required fields")

        self._read_structure(data)
        epsilon_infinity = self._reference_value(data, "epsilon_infinity")
        born_charges = self._reference_value(data, "born_effective_charges")
        self.zerof_optical_dielectric = self._array_or_default(epsilon_infinity, (3, 3), as_list=True)
        self.zerof_static_dielectric = self.zerof_optical_dielectric
        self.born_charges = self._array_or_default(born_charges, (self.nions, 3, 3), as_list=True)

        chi2 = data.get("chi2")
        if chi2 is not None:
            self.nonlinear_optical_susceptibility = np.asarray(chi2, dtype=float)

        self._deps_dr = self._read_deps_dr(data)
        self._read_modes(data)
        self._read_force_constants(data)

        if self.nomass_hessian_has_been_set:
            self.calculate_mass_weighted_normal_modes()
        elif self._deps_dr is not None:
            self._recalculate_raman_tensors_from_deps_dr()
        return

    @staticmethod
    def _array_or_default(values, shape, as_list=False):
        """Return an array-like value or zeros if it is absent."""
        array = np.zeros(shape, dtype=float) if values is None else np.asarray(values, dtype=float)
        return array.tolist() if as_list else array

    @staticmethod
    def _reference_value(data, key):
        """Return a top-level value or its lower-level reference-block equivalent."""
        if data.get(key) is not None:
            return data[key]
        reference = data.get("reference", {})
        return reference.get(key) if isinstance(reference, dict) else None

    def _read_structure(self, data):
        """Read the structure, species, and active/program masses."""
        structure = data.get("structure")
        if structure is None:
            raise ValueError("Finite-field JSON reader requires a structure block")

        lattice = np.asarray(structure["lattice_angstrom"], dtype=float)
        symbols = [cleanup_symbol(symbol) for symbol in structure["symbols"]]
        self.nions = len(symbols)
        self.ncells = 1
        self.volume = float(structure.get("volume_angstrom3", abs(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))))

        cell = UnitCell(lattice[0].tolist(), lattice[1].tolist(), lattice[2].tolist(), units="Angstrom")
        cell.set_element_names(symbols)
        if structure.get("positions_fractional") is not None:
            cell.set_fractional_coordinates(structure["positions_fractional"])
        elif structure.get("positions_cartesian_angstrom") is not None:
            cell.set_xyz_coordinates(structure["positions_cartesian_angstrom"], units="Angstrom")
        else:
            raise ValueError("Finite-field JSON reader requires fractional or Cartesian positions")

        self.species = []
        self.atom_type_list = []
        self.ions_per_type = []
        for symbol in symbols:
            if symbol not in self.species:
                self.species.append(symbol)
                self.ions_per_type.append(0)
            species_index = self.species.index(symbol)
            self.atom_type_list.append(species_index)
            self.ions_per_type[species_index] += 1
        self.nspecies = len(self.species)
        self.species_list = symbols

        active_mass_by_species = {}
        for symbol in self.species:
            mass = average_masses.get(symbol)
            if mass is None:
                raise ValueError(f"No built-in PDielec average mass is available for {symbol}")
            active_mass_by_species[symbol] = mass

        program_masses = self._program_masses(data, symbols)
        if program_masses is None:
            program_masses = [active_mass_by_species[symbol] for symbol in symbols]
        self.program_mass_dictionary = {}
        for symbol, mass in zip(symbols, program_masses):
            self.program_mass_dictionary.setdefault(symbol, float(mass))

        self.masses_per_type = [active_mass_by_species[symbol] for symbol in self.species]
        self.masses = [self.masses_per_type[index] for index in self.atom_type_list]
        cell.set_atomic_masses(self.masses)
        self.unit_cells = [cell]
        self.volumes = [self.volume]
        return

    def _program_masses(self, data, symbols):
        """Return per-atom program masses from the JSON, if present."""
        candidates = [
            data.get("program_masses_amu"),
            data.get("masses_amu"),
        ]
        atoms = data.get("atoms")
        if atoms is not None:
            candidates.append([atom.get("mass_amu") for atom in atoms])
        for candidate in candidates:
            if candidate is None:
                continue
            masses = [candidate.get(symbol) for symbol in symbols] if isinstance(candidate, dict) else candidate
            if len(masses) == len(symbols) and all(mass is not None for mass in masses):
                return [float(mass) for mass in masses]
        return None

    def _read_modes(self, data):
        """Read frequencies, JSON normal modes, and projected Raman tensors."""
        modes = data.get("modes", [])
        if len(modes) == 0:
            self.frequencies = []
            self.mass_weighted_normal_modes = []
            self.raman_tensors = None
            return

        def mode_frequency(mode):
            frequency = float(mode.get("frequency_cm-1", 0.0))
            return -frequency if mode.get("imaginary", False) and frequency > 0.0 else frequency

        modes = sorted(modes, key=mode_frequency)
        self.frequencies = [mode_frequency(mode) for mode in modes]
        if all("mass_weighted_eigenvector" in mode for mode in modes):
            self.mass_weighted_normal_modes = [
                np.asarray(mode["mass_weighted_eigenvector"], dtype=float).tolist()
                for mode in modes
            ]
        else:
            self.mass_weighted_normal_modes = []
        if all("raman_tensor" in mode for mode in modes):
            self.raman_tensors = [
                np.asarray(mode["raman_tensor"], dtype=float)
                for mode in modes
            ]
        else:
            self.raman_tensors = None
        return

    def _read_force_constants(self, data):
        """Read unweighted force constants and store the atomic-unit Hessian."""
        matrix = data.get("force_constants_matrix_sym_eV_per_angstrom2")
        if matrix is None:
            matrix = data.get("force_constants_matrix_eV_per_angstrom2")
        if matrix is None and data.get("force_constants_eV_per_angstrom2") is not None:
            matrix = self._flatten_force_constants(data["force_constants_eV_per_angstrom2"])
        if matrix is None:
            matrix = self._force_constants_from_displacements(data)
        if matrix is None:
            return

        hessian = np.asarray(matrix, dtype=float)
        expected = 3 * self.nions
        if hessian.shape != (expected, expected):
            raise ValueError(f"Force-constant matrix has shape {hessian.shape}, expected {(expected, expected)}")
        hessian = 0.5 * (hessian + hessian.T)
        ev_per_angstrom2_to_hartree_per_bohr2 = (1.0 / hartree2ev) / (angs2bohr * angs2bohr)
        self.nomass_hessian = hessian * ev_per_angstrom2_to_hartree_per_bohr2
        self.nomass_hessian_has_been_set = True
        return

    def _flatten_force_constants(self, force_constants):
        """Flatten Phi[a,k,b,l] to a (3N, 3N) matrix."""
        phi = np.asarray(force_constants, dtype=float)
        if phi.shape != (self.nions, 3, self.nions, 3):
            raise ValueError(f"Force constants have shape {phi.shape}, expected {(self.nions, 3, self.nions, 3)}")
        matrix = np.zeros((3 * self.nions, 3 * self.nions), dtype=float)
        for atom_a in range(self.nions):
            for cart_a in range(3):
                row = 3 * atom_a + cart_a
                for atom_b in range(self.nions):
                    for cart_b in range(3):
                        matrix[row, 3 * atom_b + cart_b] = phi[atom_a, cart_a, atom_b, cart_b]
        return matrix

    def _force_constants_from_displacements(self, data):
        """Build a force-constant matrix from central finite-difference forces."""
        displacements = data.get("cartesian_displacements")
        if displacements is None:
            return None
        matrix = np.zeros((3 * self.nions, 3 * self.nions), dtype=float)
        seen = set()
        for item in displacements:
            if "forces_plus_eV_per_angstrom" not in item or "forces_minus_eV_per_angstrom" not in item:
                continue
            atom = int(item["atom"])
            cart = int(item["cart"])
            delta = float(item["delta_angstrom"])
            forces_plus = np.asarray(item["forces_plus_eV_per_angstrom"], dtype=float)
            forces_minus = np.asarray(item["forces_minus_eV_per_angstrom"], dtype=float)
            if forces_plus.shape != (self.nions, 3) or forces_minus.shape != (self.nions, 3):
                raise ValueError(f"Force arrays for atom={atom} cart={cart} have inconsistent shapes")
            row = 3 * atom + cart
            matrix[row, :] = -((forces_plus - forces_minus) / (2.0 * delta)).reshape(3 * self.nions)
            seen.add((atom, cart))
        if len(seen) == 0:
            return None
        missing = [(atom, cart) for atom in range(self.nions) for cart in range(3) if (atom, cart) not in seen]
        if missing:
            raise ValueError(f"Missing force finite differences for Cartesian displacements: {missing}")
        return matrix

    def _read_deps_dr(self, data):
        """Read Cartesian dielectric derivatives, if present."""
        deps_dr = data.get("deps_dr")
        if deps_dr is None:
            diagnostics = data.get("diagnostics", {})
            deps_dr = diagnostics.get("deps_dr") if isinstance(diagnostics, dict) else None
        if deps_dr is None:
            deps_dr = self._deps_dr_from_displacements(data)
        if deps_dr is None:
            return None
        deps_dr = np.asarray(deps_dr, dtype=float)
        if deps_dr.shape != (self.nions, 3, 3, 3):
            raise ValueError(f"deps_dr has shape {deps_dr.shape}, expected {(self.nions, 3, 3, 3)}")
        return deps_dr

    def _deps_dr_from_displacements(self, data):
        """Build d epsilon / d r from central finite-difference displacement records."""
        displacements = data.get("cartesian_displacements")
        if displacements is None:
            return None
        deps_dr = np.zeros((self.nions, 3, 3, 3), dtype=float)
        seen = set()
        for item in displacements:
            if "epsilon_plus" not in item or "epsilon_minus" not in item:
                continue
            atom = int(item["atom"])
            cart = int(item["cart"])
            delta = float(item["delta_angstrom"])
            epsilon_plus = np.asarray(item["epsilon_plus"], dtype=float)
            epsilon_minus = np.asarray(item["epsilon_minus"], dtype=float)
            if epsilon_plus.shape != (3, 3) or epsilon_minus.shape != (3, 3):
                raise ValueError(f"Dielectric arrays for atom={atom} cart={cart} have inconsistent shapes")
            deps_dr[atom, cart] = (epsilon_plus - epsilon_minus) / (2.0 * delta)
            seen.add((atom, cart))
        if len(seen) == 0:
            return None
        missing = [(atom, cart) for atom in range(self.nions) for cart in range(3) if (atom, cart) not in seen]
        if missing:
            raise ValueError(f"Missing dielectric finite differences for Cartesian displacements: {missing}")
        return deps_dr

    def _update_raman_tensors_after_mode_recalculation(self, old_modes, old_raman_tensors):
        """Rebuild Raman tensors from raw derivatives when available."""
        if self._deps_dr is not None:
            self._recalculate_raman_tensors_from_deps_dr()
            return
        super()._update_raman_tensors_after_mode_recalculation(old_modes, old_raman_tensors)
        return

    def calculate_mass_weighted_normal_modes(self):
        """Regenerate modes from finite-field force constants using active masses."""
        if self.nomass_hessian_has_been_set:
            old_modes = np.array(self.mass_weighted_normal_modes, dtype=float, copy=True)
            old_raman_tensors = None
            if self.raman_tensors is not None:
                old_raman_tensors = [np.array(tensor, dtype=complex, copy=True) for tensor in self.raman_tensors]
            masses = np.asarray(self.masses, dtype=float) * amu
            self.hessian = self._modify_mass_weighting(self.nomass_hessian, masses)
            if self.eckart:
                self.hessian = self.project(self.hessian)
            self.mass_weighted_normal_modes, self.frequencies = calculate_normal_modes_and_frequencies(self.hessian)
            if self._deps_dr is not None:
                self._recalculate_raman_tensors_from_deps_dr()
            elif old_raman_tensors is not None and old_modes.size != 0:
                self.raman_tensors = self._transform_raman_tensors_between_mode_bases(
                    old_modes,
                    self.mass_weighted_normal_modes,
                    old_raman_tensors,
                )
            return self.mass_weighted_normal_modes
        return super().calculate_mass_weighted_normal_modes()

    def _recalculate_raman_tensors_from_deps_dr(self):
        """Project d epsilon / d r onto the active mass-weighted modes."""
        if self._deps_dr is None:
            return
        if not isinstance(self.mass_weighted_normal_modes, np.ndarray) and not self.mass_weighted_normal_modes:
            return
        sqrt_volume = math.sqrt(self.volume)
        masses = np.asarray(self.masses, dtype=float)
        modes = np.asarray(self.mass_weighted_normal_modes, dtype=float)
        tensors = []
        for mode in modes:
            tensor = np.zeros((3, 3), dtype=float)
            for atom, mass in enumerate(masses):
                inv_sqrt_mass = 1.0 / math.sqrt(mass)
                for cart in range(3):
                    tensor += sqrt_volume * mode[atom, cart] * inv_sqrt_mass * self._deps_dr[atom, cart]
            tensors.append(tensor)
        self.raman_tensors = tensors
        return
