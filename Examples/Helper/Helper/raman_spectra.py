#!/usr/bin/python
# ruff: noqa: E402
"""Calculate powder and crystal Raman spectra using helper routines."""
#
# Copyright 2026 John Kendrick
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
# You should have received a copy of the MIT License
# along with this program, if not see https://opensource.org/licenses/MIT
#
import sys
from pathlib import Path

import numpy as np

home_directory = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(home_directory))

from PDielec.GUI.SingleCrystalLayer import SingleCrystalLayer
from PDielec.HelperRoutines import (
    calculate_crystal_raman_spectrum,
    calculate_dft_permittivity_object,
    calculate_powder_raman_spectrum,
    get_material,
)
from PDielec.Materials import External
from PDielec.Utilities import get_reader


def powder_raman_test():
    """Calculate a powder Raman spectrum from a CASTEP Raman example."""
    filename = home_directory / "Examples" / "Powder_Raman" / "Castep" / "ramand3.castep"
    reader = get_reader(str(filename), "Castep")
    reader.read_output()

    frequencies_cm1 = np.arange(1.0, 200.0, 0.2)
    matrix = get_material("ptfe")
    spectrum = calculate_powder_raman_spectrum(
        frequencies_cm1,
        reader,
        matrix=matrix,
        volume_fraction=0.1,
        sigma=5.0,
        laser_wavelength_nm=785.0,
        polarisation="HV",
    )

    print("")
    print("Powder Raman results")
    print("      freq", " intensity")
    for frequency, intensity in zip(frequencies_cm1[::20], spectrum[::20]):
        print(f"{frequency:10.2f}{intensity:12.5e}")


def crystal_raman_test():
    """Calculate a crystal Raman spectrum from a CRYSTAL Raman example."""
    filename = home_directory / "Examples" / "Crystal_Raman" / "Rutile" / "raman.log"
    reader = get_reader(str(filename), "Crystal")
    reader.read_output()
    permittivity = calculate_dft_permittivity_object(reader, sigma=5.0)
    rutile = External("Rutile", permittivity_object=permittivity, cell=reader.get_unit_cell())

    air = get_material("air")
    layers = [
        SingleCrystalLayer(air, thickness=1.0, thickness_unit="um"),
        SingleCrystalLayer(rutile, hkl=[0, 0, 1], thickness=1.0, thickness_unit="um", dielectricFlag=True),
        SingleCrystalLayer(air, thickness=1.0, thickness_unit="um"),
    ]

    frequencies_cm1 = np.arange(1.0, 800.0, 1.0)
    spectrum, mode_frequencies, mode_intensities, mode_sigmas = calculate_crystal_raman_spectrum(
        frequencies_cm1,
        reader,
        layers,
        incident_angle=45.0,
        global_azimuthal_angle=0.0,
        sigma=5.0,
        laser_wavelength_nm=532.0,
        incident_polarisation="p",
        detected_polarisation="unpolarised",
        approximate_es=True,
    )

    print("")
    print("Crystal Raman active modes")
    print("      freq", " intensity", " sigma")
    for frequency, intensity, sigma in zip(mode_frequencies, mode_intensities, mode_sigmas):
        print(f"{frequency:10.2f}{intensity:12.5e}{sigma:10.2f}")
    print("")
    print(f"Crystal Raman spectrum points: {len(spectrum)}")


def main():
    """Run the powder and crystal Raman helper examples."""
    powder_raman_test()
    crystal_raman_test()


if __name__ == "__main__":
    main()
