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
"""Helper Routines - Useful for scripting and in Jupyter Notebooks.

These routines are a useful starting point for investigating the way the code operates
"""

import numpy as np

from PDielec import Calculator, DielectricFunction, Utilities
from PDielec.Constants import amu, average_masses, isotope_masses, wavenumber
from PDielec.GUI.SingleCrystalScenarioTab import solve_single_crystal_equations
from PDielec.Materials import External, MaterialsDataBase


def calculateDFTPermittivityObject(reader,sigma=5.0,eckart=True,mass_definition="Average"):
    """Define a permittivity object for the DFT calculation from the given reader.

    This routine reads in the details of the DFT calculation
    Calculates the hessian and normal modes
    Calculates the oscillator strength
    Finally it defines a permittivity object from which the permittivity can be calculated

    Parameters
    ----------
    reader : outputReader object
        An output reader object, see :class:`~PDielec.GenericOutputReader` and its sub-classes.
    sigma  : float
        Default Lorentzian widths in cm-1
    eckart  : boolean
        If true apply eckart conditions in the reader
    mass_definition : string or list of floats
        Default is 'Average'
        If a list of floats it is a list of the masses to be used.
        If a string it must be one of 'Average', 'Isotope' or 'Program' (case is irrelevant)
        - Average causes the average weight for the element to be used
        - Isotope causes the most common isotope weight to be used
        - Program means use the mass used by the program

    Returns
    -------
    permittivityObject (see :class:`~PDielec.DielectricFunction.DielectricFunction` and its sub-classes)

    """
    reader.eckart = eckart
    epsilon_inf = reader.zerof_optical_dielectric
    cell = reader.get_unit_cell()
    volume_au = cell.getVolume("Bohr")
    mass_dictionary = []
    reader.reset_masses()
    if isinstance(mass_definition, (list, np.ndarray)):
        reader.set_edited_masses(mass_definition)
    else:
        mass_definition = mass_definition.lower()
        if mass_definition == "average":
            reader.change_masses(average_masses, mass_dictionary)
        elif mass_definition == "program":
            pass
        elif mass_definition == "isotope":
            reader.change_masses(isotope_masses, mass_dictionary)
        else:
            print("Helper: Error unkown mass definition", mass_definition )
    masses = np.array(reader.masses)*amu
    # The reader uses the internal masses to calculate the massweighted normal modes
    mass_weighted_normal_modes = reader.calculate_mass_weighted_normal_modes()
    masses = np.array(reader.masses)*amu
    # The frequencies will change if the masses change, so get them after changing the masses
    frequencies_cm1 = np.array(reader.frequencies)
    frequencies_au = wavenumber*frequencies_cm1
    sigmas_cm1 = [ sigma for i in frequencies_cm1 ]
    sigmas_au = wavenumber*np.array(sigmas_cm1)
    born_charges = np.array(reader.born_charges)
    if reader.type == "Experimental output":
        # Obtain oscillator strength from reader
        oscillator_strengths = np.array(reader.oscillator_strengths)
    else:
        # Calculate oscillator strengths
        # calculate normal modes in xyz coordinate space
        normal_modes = Calculator.normal_modes(masses, mass_weighted_normal_modes)
        # from the normal modes and the born charges calculate the oscillator strengths of each mode
        oscillator_strengths = Calculator.oscillator_strengths(normal_modes, born_charges)
    # calculate the intensities from the trace of the oscillator strengths
    intensities = Calculator.infrared_intensities(oscillator_strengths)
    # Decide which modes to select
    modes_selected = []
    mode_list = []
    for f,intensity in zip(frequencies_cm1,intensities):
        if f > 10.0 and intensity > 1.0E-6:
            modes_selected.append(True)
        else:
            modes_selected.append(False)
    mode_list = [i for i,mode in enumerate(modes_selected) if mode]
    drude = False
    drude_plasma_au = 0
    drude_sigma_au = 0
    permittivityObject = DielectricFunction.DFT(
                                  mode_list, frequencies_au, sigmas_au, oscillator_strengths,
                                  volume_au, drude, drude_plasma_au, drude_sigma_au )
    # Add the optical permittivity to the dielctric function
    permittivityObject.setEpsilonInfinity(epsilon_inf)
    return permittivityObject

def getMaterial(name,dataBaseName="MaterialsDataBase.xlsx",eckart=True,mass_definition="Average",debug=False):
    """Get a material with the given name.

    If the name is a file name, it is treated as a DFT (Density Functional Theory) or experimental file.
    If the name is a material name in the material database, this is used instead.
    Information about the material class can be found here (:class:`~PDielec.Materials.Material`)

    Parameters
    ----------
    name : str
        The name of the material, which can be a file name or a material name existing in the database.
    dataBaseName : str
        The name of the database which will be searched from the material
    eckart  : boolean, optional
        If true apply eckart conditions in the reader. Defaults to true.
    mass_definition : string or list of floats, optional
        Default is 'Average'
        If a list of floats it is a list of the masses to be used.
        If a string it must be one of 'Average', 'Isotope' or 'Program' (case is irrelevant)
        - Average causes the average weight for the element to be used
        - Isotope causes the most common isotope weight to be used
        - Program means use the mass used by the program
    debug : boolean
        If true debugging output is provided

    Returns
    -------
    Material (see :class:`~PDielec.Materials.Material`)
        The material obtained based on the given name.

    Examples
    --------
    ```
    material1 = getMaterial("example.dft")
    material2 = getMaterial("gold")
    ```

    """
    # Let's see if the name is a file name that can be read
    program = Utilities.find_program_from_name(name)
    if debug:
        print(f"getMaterial: program = {program}")
    if len(program) > 1:
        reader = Utilities.get_reader(name,program,debug)
        reader.read_output()
        if debug:
            print("getMaterial: reader.print()")
            reader.print()
        permittivityObject=calculateDFTPermittivityObject(reader,sigma=5.0,eckart=eckart,mass_definition=mass_definition)
        cell = reader.get_unit_cell()
        if debug:
            print("getMaterial: cell.print()")
            cell.print()
        material = External("Dielectric layer",permittivityObject=permittivityObject,cell=cell)
    else:
        dataBase = MaterialsDataBase(dataBaseName)
        sheets = dataBase.getSheetNames()
        if debug:
            print("getMaterial: sheets",sheets)
        if name in sheets:
            material = dataBase.getMaterial(name)
        else:
            print("Material name not valid: ",name)
            print("Available materials:     ",sheets)
            material = None
    return material

def calculateSingleCrystalSpectrum(frequencies_cm1, layers, incident_angle, global_azimuthal_angle, method="Scattering matrix"):
    """Calculate a single crystal spectrum.

    Calculate a single crystal spectrum from the frequencies, a list of layers (:class:`~PDielec.GUI.SingleCrystalLayer.SingleCrystalLayer`), the incident angle
    the global azimuthal angle and optional specification of the method of calculation.
    Information about the layer class can be found here (:class:`~PDielec.GUI.SingleCrystalLayer.SingleCrystalLayer`)

    Parameters
    ----------
    frequencies_cm1 : list
        A list of frequencies in cm-1.
    layers : list
        A list of layers ( :class:`~PDielec.GUI.SingleCrystalLayer.SingleCrystalLayer`)
    incident_angle : float
        The incident angle in degrees.
    global_azimuthal_angle : float
        The global azimuthal angle in degrees.
    method : str
        The method for solving Maxwell's equation either 'Scattering matrix' (the default) or 'Transfer matrix. See (:class:`~PDielec.GTMcore.ScatteringMatrixSystem` or :class:`~PDielec.GTMcore.TransferMatrixSystem`) for more information.

    Returns
    -------
    ((reflp,refls) (tranp, trans) (absop, absos)) : list of nd.arrrays
        - reflp the p- polarised reflection at each frequency as an np.array
        - refls the s- polarised reflection at each frequency as an np.array
        - tranp the p- polarised transmittance at each frequency as an np.array
        - trans the s- polarised transmittance at each frequency as an np.array
        - absop the p- polarised absorptance at each frequency as an np.array
        - absos the s- polarised absorptance at each frequency as an np.array

    Example
    -------

    ::

        air = getMaterial('air')
        ptfe = getMaterial('ptfe')
        Sapphire = getMaterial('Sapphire')
        layers = []
        # Add a substrate
        layers.append(SingleCrystalLayer(ptfe,thickness=1.0,thicknessUnit='um'))
         # Add the dielectric layer
        layers.append(SingleCrystalLayer(Sapphire,hkl=[0,0,1],azimuthal=0.0,thickness=1.0,thicknessUnit='um',
                                         incoherentOption='Coherent'))
         # Add the superstrate
        layers.append(SingleCrystalLayer(air,thickness=1.0,thicknessUnit='um'))
        frequencies_cm1 = np.arange( 0, 200, 0.2 )
        incident_angle = 80.0
        global_azimuthal_angle = 0.0
        (reflectance, transmittance, absorptance) = calculateSingleCrystalSpectrum(frequencies_cm1,layers,incident_angl e, global_azimuthal_angle, method='Scattering matrix')

    """
    theta = np.radians(0.0)
    phi = np.radians(0.0)
    psi = np.radians(global_azimuthal_angle)
    incident_angle = np.radians(incident_angle)
    exponent_threshold = 700
    reflp = []
    refls = []
    tranp = []
    trans = []
    absop = []
    absos = []
    for v_cm1 in frequencies_cm1:
        (v,r,R,t,T,epsilon,errors,largest_exponent) = solve_single_crystal_equations(layers,method,
                                                              theta,phi,psi,incident_angle,
                                                              exponent_threshold,v_cm1)
        reflp.append( R[0] + R[2] )
        refls.append( R[1] + R[3] )
        tranp.append( T[0] )
        trans.append( T[1] )
        absop.append( 1.0 - R[0] - R[2] - T[0] )
        absos.append( 1.0 - R[1] - R[3] - T[1] )
    return np.array([reflp,refls]), np.array([tranp,trans]), np.array([absop,absos])

def calculatePowderSpectrum(frequencies_cm1, dielectric, matrix, volume_fraction, method="Maxwell-Garnett"):
    """Calculate the powder IR spectrum of a mixture of spherical dielectric particles in a matrix with a given volume fraction.

    Parameters
    ----------
    frequencies_cm1 : list
        A list of frequencies in cm-1.
    dielectric : material
        A material with a tensor permittivity.
    matrix : material
        A material with a scalar permittivity.
    volume_fraction : float
        The volume fraction of dielectric.
    method : str
        One of 'Maxwell-Garnett', 'Bruggeman'.

    Returns
    -------
    absorption,permittivity : tuple of np.arrays
        absorption: Absorption coefficients in cm-1.
        permittivity: Complex permittivities

    Example
    -------
    ```
    frequencies_cm1 = np.arange( 0, 200, 0.2 )
    matrix = getMaterial('ptfe')
    dielectric = getMaterial('Sapphire')
    method = 'Maxwell-Garnett' 
    volume_fraction = 0.1
    absorption,permittivity = calculatePowderSpectrum(frequencies_cm1,dielectric, matrix, volume_fraction)
    ```

    """
    method = method.lower()
    particle_size_mu = 0
    particle_sigma_mu = 0
    matrixPermittivityFunction = matrix.getPermittivityFunction()
    shape = "Sphere"
    depolarisation = Calculator.initialise_sphere_depolarisation_matrix()
    concentration = 1.0
    atr_refractive_index = 5.0
    atr_theta = 45
    atr_spolfraction  = 1.0
    bubble_vf = 0.0
    bubble_radius = 0.0
    previous_solution_shared = np.eye( 3 )
    crystalPermittivityFunction = dielectric.getPermittivityFunction()
    results = []
    unitMatrix = np.eye(3)
    for v_cm1 in frequencies_cm1:
        crystalPermittivity = crystalPermittivityFunction(v_cm1)
        #
        # Cope with the case that the dielectric is isotropic
        #
        if dielectric.isScalar():
            crystalPermittivity = crystalPermittivity*unitMatrix
        result = Calculator.solve_effective_medium_equations(method,volume_fraction,
                       particle_size_mu,particle_sigma_mu,matrixPermittivityFunction,
                       shape,depolarisation,concentration,
                       atr_refractive_index,atr_theta,atr_spolfraction,
                       bubble_vf,bubble_radius,previous_solution_shared,(v_cm1,crystalPermittivity))
        results.append(result)
    # Prepare lists for results
    permittivity = []
    absorptionCoefficient = []
    molarAbsorptionCoefficient = []
    sp_atr = []
    for _v,_method,_size_mu,_size_sigma,_shape,_data,trace, absorption_coefficient,molar_absorption_coefficient,spatr in results:
         permittivity.append(trace)
         absorptionCoefficient.append(absorption_coefficient)
         molarAbsorptionCoefficient.append(molar_absorption_coefficient)
         sp_atr.append(spatr)
    return np.array(absorptionCoefficient), np.array(permittivity)

def maxwell_garnett(em, ei, f):
    """Calculate the dielectric constant of a mixture using Maxwell-Garnett.

    Parameters
    ----------
    em : float
        The dielectric constant of the host
    ei : float
        The dielectric constant of the inclusion
    f : float
        The volume fraction of the inclusion

    Returns
    -------
    float
        The dielectric constant of the mixture

    """
    return em*( 2*f*(ei-em)+ ei + 2*em) / ( 2*em + ei - f*(ei-em))
