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

import logging

import numpy as np

import PDielec.GTMcore as GTM
from PDielec import Calculator, DielectricFunction, Utilities
from PDielec.Constants import amu, average_masses, boltzmann_si, isotope_masses, planck_si, speed_light_si, wavenumber
from PDielec.GUI.CrystalScenarioTab import solve_single_crystal_equations
from PDielec.LayeredRamanCalculator import LayeredRamanCalculator, RamanLayer, lorentzian_broaden
from PDielec.Materials import External, MaterialsDataBase

logger = logging.getLogger(__name__)
gtm_methods = {"Coherent": GTM.CoherentLayer,
               "Incoherent (intensity)": GTM.IncoherentIntensityLayer,
               "Incoherent (phase cancelling)": GTM.IncoherentPhaseLayer,
               "Incoherent (phase averaging)": GTM.IncoherentAveragePhaseLayer,
               "Incoherent (non-reflective)": GTM.IncoherentThickLayer}



def calculate_dft_permittivity_object(reader,sigma=5.0,eckart=True,mass_definition="Average"):
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
    permittivity_object (see :class:`~PDielec.DielectricFunction.DielectricFunction` and its sub-classes)

    """
    reader.eckart = eckart
    epsilon_inf = reader.zerof_optical_dielectric
    cell = reader.get_unit_cell()
    volume_au = cell.get_volume("Bohr")
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
            logger.error(f"Helper: Error unkown mass definition {mass_definition}")
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
    if reader.type == "Experimental output" and reader.oscillator_strengths is not None and np.any(reader.oscillator_strengths):
        # Obtain oscillator strength from reader, or from ExperimentOutputReader's Born-charge calculation.
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
    permittivity_object = DielectricFunction.DFT(
                                  mode_list, frequencies_au, sigmas_au, oscillator_strengths,
                                  volume_au, drude, drude_plasma_au, drude_sigma_au )
    # Add the optical permittivity to the dielctric function
    permittivity_object.set_epsilon_infinity(epsilon_inf)
    return permittivity_object

def get_material(name,dataBaseName="MaterialsDataBase.xlsx",eckart=True,mass_definition="Average",debug=False):
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
    material1 = get_material("example.dft")
    material2 = get_material("gold")
    ```

    """
    # Let's see if the name is a file name that can be read
    program = Utilities.find_program_from_name(name)
    if debug:
        logger.debug(f"get_material: program = {program}")
    if len(program) > 1:
        reader = Utilities.get_reader(name,program)
        reader.read_output()
        if debug:
            logger.debug("get_material: reader.print()")
            reader.print()
        permittivity_object=calculate_dft_permittivity_object(reader,sigma=5.0,eckart=eckart,mass_definition=mass_definition)
        cell = reader.get_unit_cell()
        if debug:
            logger.debug("get_material: cell.print()")
            cell.print()
        material = External("Dielectric layer",permittivity_object=permittivity_object,cell=cell)
    else:
        dataBase = MaterialsDataBase(dataBaseName)
        sheets = dataBase.get_sheet_names()
        if debug:
            logger.debug(f"get_material: sheets {sheets}")
        if name in sheets:
            material = dataBase.get_material(name)
        else:
            logger.debug(f"Material name not valid: {name}")
            logger.warning(f"Available materials: {sheets}")
            material = None
    return material

def calculate_crystal_infrared_spectrum(frequencies_cm1, layers, incident_angle, global_azimuthal_angle,
                                        method="Scattering matrix"):
    """Calculate a crystal infrared spectrum.

    Calculate a single crystal spectrum from the frequencies, a list of layers
    (:class:`~PDielec.GUI.SingleCrystalLayer.SingleCrystalLayer`), the incident angle the global azimuthal angle and
    optional specification of the method of calculation. Information about the layer class can be found here
    (:class:`~PDielec.GUI.SingleCrystalLayer.SingleCrystalLayer`)

    Parameters
    ----------
    frequencies_cm1 : list
        A list of frequencies in cm-1.
    layers : list
        A list of layers (:class:`~PDielec.GUI.SingleCrystalLayer.SingleCrystalLayer`).
    incident_angle : float
        The incident angle in degrees.
    global_azimuthal_angle : float
        The global azimuthal angle in degrees.
    method : str
        The method for solving Maxwell's equation, either 'Scattering matrix' (the default) or 'Transfer matrix'.
        See (:class:`~PDielec.GTMcore.ScatteringMatrixSystem` or :class:`~PDielec.GTMcore.TransferMatrixSystem`)
        for more information.

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

        air = get_material('air')
        ptfe = get_material('ptfe')
        Sapphire = get_material('Sapphire')
        layers = []
        # Add a substrate
        layers.append(SingleCrystalLayer(ptfe,thickness=1.0,thickness_unit='um'))
         # Add the dielectric layer
        layers.append(SingleCrystalLayer(Sapphire,hkl=[0,0,1],azimuthal=0.0,thickness=1.0,thickness_unit='um',
                                         incoherentOption='Coherent'))
         # Add the superstrate
        layers.append(SingleCrystalLayer(air,thickness=1.0,thickness_unit='um'))
        frequencies_cm1 = np.arange( 0, 200, 0.2 )
        incident_angle = 80.0
        global_azimuthal_angle = 0.0
        (reflectance, transmittance, absorptance) = calculate_crystal_infrared_spectrum(
            frequencies_cm1,layers,incident_angle, global_azimuthal_angle, method='Scattering matrix')

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

def calculate_powder_infrared_spectrum(frequencies_cm1, dielectric, matrix, volume_fraction, method="Maxwell-Garnett"):
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
    matrix = get_material('ptfe')
    dielectric = get_material('Sapphire')
    method = 'Maxwell-Garnett' 
    volume_fraction = 0.1
    absorption,permittivity = calculate_powder_infrared_spectrum(frequencies_cm1,dielectric, matrix, volume_fraction)
    ```

    """
    method = method.lower()
    particle_size_mu = 0
    particle_sigma_mu = 0
    matrixPermittivityFunction = matrix.get_permittivity_function()
    shape = "Sphere"
    depolarisation = Calculator.initialise_sphere_depolarisation_matrix()
    concentration = 1.0
    bubble_vf = 0.0
    bubble_radius = 0.0
    previous_solution_shared = np.eye( 3 )
    crystalPermittivityFunction = dielectric.get_permittivity_function()
    results = []
    unit_matrix = np.eye(3)
    for v_cm1 in frequencies_cm1:
        crystalPermittivity = crystalPermittivityFunction(v_cm1)
        #
        # Cope with the case that the dielectric is isotropic
        #
        if dielectric.is_scalar():
            crystalPermittivity = crystalPermittivity*unit_matrix
        result = Calculator.solve_effective_medium_equations(method,volume_fraction,
                       particle_size_mu,particle_sigma_mu,matrixPermittivityFunction,
                       shape,depolarisation,concentration,
                       bubble_vf,bubble_radius,previous_solution_shared,(v_cm1,crystalPermittivity))
        results.append(result)
    # Prepare lists for results
    permittivity = []
    absorptionCoefficient = []
    molarAbsorptionCoefficient = []
    for _v,_method,_size_mu,_size_sigma,_shape,_data,trace,absorption_coefficient,molar_absorption_coefficient in results:
         permittivity.append(trace)
         absorptionCoefficient.append(absorption_coefficient)
         molarAbsorptionCoefficient.append(molar_absorption_coefficient)
    return np.array(absorptionCoefficient), np.array(permittivity)

def _raman_linewidths(reader, sigma):
    """Return per-mode Raman linewidths in cm-1."""
    frequencies_cm1 = np.asarray(reader.frequencies, dtype=float)
    if np.isscalar(sigma):
        return np.full(len(frequencies_cm1), float(sigma))
    sigmas_cm1 = np.asarray(sigma, dtype=float)
    if len(sigmas_cm1) != len(frequencies_cm1):
        raise ValueError("sigma must be a scalar or have one value per reader frequency")
    return sigmas_cm1

def _raman_modes_selected(reader, raman_tensors, modes_selected):
    """Return a boolean mask selecting Raman-active optical modes."""
    frequencies_cm1 = np.asarray(reader.frequencies, dtype=float)
    if modes_selected is not None:
        selected = np.asarray(modes_selected, dtype=bool)
        if len(selected) != len(frequencies_cm1):
            raise ValueError("modes_selected must have one value per reader frequency")
        return selected
    activities = Calculator.raman_intensities(raman_tensors, reader.volume)[:, 0]
    return np.asarray([frequency > 10.0 and activity > 1.0e-6
                       for frequency, activity in zip(frequencies_cm1, activities)], dtype=bool)

def _matrix_optical_permittivity(matrix):
    """Return the scalar optical permittivity of a matrix material."""
    if matrix is None or matrix == "none":
        return None
    eps = matrix.get_optical_permittivity()
    if np.isscalar(eps) or (hasattr(eps, "ndim") and eps.ndim == 0):
        return float(np.real(eps))
    return float(np.real(np.trace(eps))) / 3.0

def calculate_powder_raman_spectrum(frequencies_cm1, reader, matrix=None, volume_fraction=1.0, sigma=5.0,
                                    laser_wavelength_nm=785.0, polarisation="HV", temperature=298.0,
                                    modes_selected=None):
    """Calculate a powder Raman spectrum for spherical particles.

    Parameters
    ----------
    frequencies_cm1 : array_like
        Raman-shift axis in cm-1.
    reader : outputReader object
        Reader containing phonon frequencies, optical permittivity and Raman tensors.
    matrix : material or None
        Matrix material.  If ``None`` or ``"none"``, no particle-field correction is applied.
    volume_fraction : float
        Crystal volume fraction in the sample.
    sigma : float or array_like
        Lorentzian half-width(s) in cm-1.
    laser_wavelength_nm : float
        Laser wavelength in nm.
    polarisation : {"VV", "VH", "HV", "Unpolarised"}
        Raman polarisation geometry.
    temperature : float
        Temperature in K.
    modes_selected : array_like of bool, optional
        Optional per-mode selection mask.  If omitted, Raman-active optical modes are selected automatically.

    Returns
    -------
    spectrum : np.ndarray
        Raman intensity on ``frequencies_cm1`` in arbitrary units.

    """
    raman_tensors = reader.get_raman_tensors()
    if raman_tensors is None or len(raman_tensors) == 0:
        return np.zeros(len(frequencies_cm1))

    frequencies_cm1 = np.asarray(frequencies_cm1, dtype=float)
    mode_frequencies = np.asarray(reader.frequencies, dtype=float)
    sigmas_cm1 = _raman_linewidths(reader, sigma)
    selected = _raman_modes_selected(reader, raman_tensors, modes_selected)

    epsilon_e = _matrix_optical_permittivity(matrix)
    no_matrix = epsilon_e is None
    epsilon_inf_i = np.asarray(reader.zerof_optical_dielectric, dtype=complex)
    I3 = np.eye(3, dtype=complex)
    L = I3 / 3.0
    N = I3 if no_matrix else Calculator.compute_internal_field_tensor(L, epsilon_inf_i, epsilon_e)

    nu_L = 1.0e7 / laser_wavelength_nm
    hc_over_k = planck_si * speed_light_si * 100.0 / boltzmann_si
    spectrum = np.zeros(len(frequencies_cm1))

    for frequency, linewidth, is_selected, raman_tensor in zip(
            mode_frequencies, sigmas_cm1, selected, raman_tensors):
        if not is_selected or abs(frequency) < 1.0:
            continue
        R_eps = np.asarray(raman_tensor, dtype=complex)
        R_eff = R_eps if no_matrix else Calculator.compute_effective_raman_tensor(R_eps, N, N)

        if polarisation in ("VV", "VH", "HV"):
            vv, vh = Calculator.compute_powder_raman_intensities(R_eff)
            intensity_factor = vv if polarisation == "VV" else vh
        elif polarisation == "Unpolarised":
            alpha = np.trace(R_eff) / 3.0
            gamma_t = 0.5 * (R_eff + R_eff.T) - alpha * I3
            kappa_t = 0.5 * (R_eff - R_eff.T)
            alpha2 = float(np.real(alpha * np.conj(alpha)))
            gamma2 = 3.0 / 2.0 * float(np.real(np.sum(gamma_t * np.conj(gamma_t))))
            kappa2 = 3.0 / 2.0 * float(np.real(np.sum(kappa_t * np.conj(kappa_t))))
            intensity_factor = 45.0 * alpha2 + 7.0 * gamma2 + 5.0 * kappa2
        else:
            raise ValueError("polarisation must be one of 'VV', 'VH', 'HV' or 'Unpolarised'")

        x = hc_over_k * frequency / temperature if temperature > 0 else 1.0e18
        n_bose = 1.0 / np.expm1(x) if x > 1.0e-6 else 1.0 / x
        nu_s = nu_L - frequency
        if nu_s <= 0.0:
            continue
        strength = (nu_s ** 4) * (n_bose + 1.0) * intensity_factor / frequency
        spectrum += strength * linewidth / ((frequencies_cm1 - frequency) ** 2 + linewidth ** 2)

    return spectrum * volume_fraction

def calculate_crystal_raman_spectrum(frequencies_cm1, reader, layers, incident_angle, global_azimuthal_angle,
                                     method="Scattering matrix", sigma=5.0, laser_wavelength_nm=532.0,
                                     incident_polarisation="p", detected_polarisation="unpolarised",
                                     temperature=298.0, n_gauss=20, collection_side="superstrate",
                                     collection_angle=-1.0, coherent_layers=False, approximate_es=False):
    """Calculate a layered crystal Raman spectrum.

    Parameters are intentionally close to :func:`calculate_crystal_infrared_spectrum`, with Raman-specific laser,
    polarisation and quadrature options added.

    Returns
    -------
    spectrum, mode_frequencies, mode_intensities, mode_sigmas : tuple of np.ndarray
        Broadened Raman spectrum and the active-mode data used to generate it.

    """
    raman_tensors = reader.get_raman_tensors()
    if raman_tensors is None or len(raman_tensors) == 0:
        zeros = np.zeros(len(frequencies_cm1))
        return zeros, np.array([]), np.array([]), np.array([])

    theta = 0.0
    phi = 0.0
    psi = np.radians(global_azimuthal_angle)
    incident_angle_rad = np.radians(incident_angle)
    exponent_threshold = 700

    for layer in layers:
        if layer.is_tensor():
            hkl = layer.get_hkl()
            if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
                raise ValueError("Cannot calculate crystal Raman spectrum for a tensor layer with hkl=[0, 0, 0]")
            layer.calculate_euler_matrix()

    superstrate = GTM.SemiInfiniteLayer(layers[0], exponent_threshold=exponent_threshold)
    substrate = GTM.SemiInfiniteLayer(layers[-1], exponent_threshold=exponent_threshold)
    selected_layers = layers[1:-1]
    gtm_layers = []
    for layer in selected_layers:
        incoherent_option = layer.get_incoherent_option()
        gtm_layers.append(gtm_methods[incoherent_option](layer, exponent_threshold=exponent_threshold))
    if method == "Scattering matrix":
        system = GTM.ScatteringMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtm_layers)
    else:
        system = GTM.TransferMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtm_layers)

    system.superstrate.set_euler(theta, phi, psi)
    system.substrate.set_euler(theta, phi, psi)
    for gtm_layer in system.layers:
        gtm_layer.set_euler(theta, phi, psi)

    G_psi = np.array([
        [np.cos(psi), -np.sin(psi), 0.0],
        [np.sin(psi),  np.cos(psi), 0.0],
        [0.0,          0.0,         1.0],
    ])
    raman_tensors_physical = [np.asarray(raman_tensor, dtype=float) for raman_tensor in raman_tensors]
    phonon_frequencies = np.asarray(reader.frequencies, dtype=float)
    sigmas_cm1 = _raman_linewidths(reader, sigma)

    raman_layers = []
    for sys_idx, layer in enumerate(selected_layers):
        if not layer.is_dielectric():
            continue
        raman_layers.append(RamanLayer(
            layer_index=sys_idx,
            phonon_frequencies_cm1=phonon_frequencies,
            raman_tensors=raman_tensors_physical,
            rotation_matrix=G_psi @ layer.euler,
        ))
    if not raman_layers:
        raise ValueError("No Raman-active dielectric layers were found in the crystal layer stack")

    collection_angle_rad = incident_angle_rad if collection_angle < 0.0 else np.radians(collection_angle)
    calculator = LayeredRamanCalculator(
        system=system,
        raman_layers=raman_layers,
        laser_frequency_cm1=1.0e7 / laser_wavelength_nm,
        incident_angle_rad=incident_angle_rad,
        incident_pol=incident_polarisation,
        detected_pol=detected_polarisation,
        temperature_K=temperature,
        linewidths_cm1=sigmas_cm1,
        n_gauss=n_gauss,
        collection_side=collection_side,
        collection_angle_rad=collection_angle_rad,
        coherent_layers=coherent_layers,
        approximate_es=approximate_es,
    )
    mode_frequencies, mode_intensities, mode_sigmas = calculator.calculate_mode_intensities()
    if len(mode_frequencies) > 0:
        spectrum = lorentzian_broaden(mode_frequencies, mode_intensities, mode_sigmas, np.asarray(frequencies_cm1))
    else:
        spectrum = np.zeros(len(frequencies_cm1))
    return spectrum, mode_frequencies, mode_intensities, mode_sigmas

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
