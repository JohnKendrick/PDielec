#!/usr/bin/python
'''
Helper Routines - Useful for scripting and Jupyter Notebooks

Copyright 2024 John Kendrick

This file is part of PDielec

This program is free software; you can redistribute it and/or modify
it under the terms of the MIT License

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

You should have received a copy of the MIT License
along with this program, if not see https://opensource.org/licenses/MIT
'''

import numpy                      as np
import PDielec.Utilities          as Utilities
import PDielec.Calculator         as Calculator
import PDielec.DielectricFunction as DielectricFunction
import PDielec.GTMcore            as GTMcore

from PDielec.GUI.SingleCrystalLayer        import SingleCrystalLayer
from PDielec.GUI.SingleCrystalScenarioTab  import solve_single_crystal_equations
from PDielec.Materials                     import Material,MaterialsDataBase,External
from PDielec.Constants                     import wavenumber, angstrom, amu

def calculateDFTPermittivityObject(reader,sigma=5.0):
    '''
    Define a permittivity object for the DFT calculation from the given reader

    This routine reads in the details of the DFT calculation
    Calculates the hessian and normal modes
    Calculates the oscillator strength
    Finally it defines a permittivity object from which the permittivity can be calculated

    Parameters
    ----------
    reader : outputReader object
        An output reader object
    sigma  : float
        Default Lorentzian widths in cm-1

    Returns
    -------
    permittivityObject

    '''
    frequencies_cm1 = np.array(reader.frequencies)
    frequencies_au = wavenumber*frequencies_cm1
    sigmas_cm1 = [ sigma for i in frequencies_cm1 ]
    sigmas_au = wavenumber*np.array(sigmas_cm1)
    epsilon_inf = reader.zerof_optical_dielectric
    cell = reader.get_unit_cell()
    volume_au = cell.getVolume('Bohr')
    masses = np.array(cell.get_atomic_masses())*amu
    mass_weighted_normal_modes = reader.calculate_mass_weighted_normal_modes()
    born_charges = np.array(reader.born_charges)
    if reader.type == 'Experimental output':
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
    for index,(f,intensity) in enumerate(zip(frequencies_cm1,intensities)):
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

def getMaterial(name,dataBaseName='MaterialsDataBase.xlsx',qmprogram='vasp'):
    '''
    Get a material with the given name.

    If the name is a file name, it is treated as a DFT (Density Functional Theory) or experimental file.
    If the name is a material name in the material database, this is used instead.

    Parameters
    ----------
    name : str
        The name of the material, which can be a file name or a material name existing in the database.
    dataBaseName : str
        The name of the database which will be searched from the material
    qmprogram : str
        The name of the QM program used in the case of Phonopy DFT calculations

    Returns
    -------
    Material
        The material obtained based on the given name.

    Examples
    --------
    ```
    material1 = get_material("example.dft")
    material2 = get_material("gold")
    ```
    '''
    # Let's see if the name is a file name that can be read
    program = Utilities.find_program_from_name(name)
    if len(program) > 1:
        reader = Utilities.get_reader(name,program,qmprogram)
        reader.read_output()
        permittivityObject=calculateDFTPermittivityObject(reader,sigma=5.0)
        cell = reader.get_unit_cell()
        material = External('Dielectric layer',permittivityObject=permittivityObject,cell=cell)
    else:
        dataBase = MaterialsDataBase(dataBaseName)
        sheets = dataBase.getSheetNames()
        if name not in sheets:
            print('Material name not valid ',name)
        material = dataBase.getMaterial(name)
    return material

def calculateSingleCrystalSpectrum(frequencies_cm1, layers, incident_angle, global_azimuthal_angle, method='Scattering matrix'):
    '''
    Calculate a single crystal spectrum.

    Calculate a single crystal spectrum from the frequencies, a list of layers, the incident angle
    the global azimuthal angle and optional specification of the method of calculation.

    Parameters
    ----------
    frequencies_cm1 : list
        A list of frequencies in cm-1.
    layers : list
        A list of layers.
    incident_angle : float
        The incident angle in degrees.
    global_azimuthal_angle : float
        The global azimuthal angle in degrees.
    method : str
        The method for solving Maxwell's equation ('Scattering matrix' or 'Transfer matrix').

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

    '''
    theta = np.radians(0.0)
    phi = np.radians(0.0)
    psi = np.radians(global_azimuthal_angle)
    incident_angle = np.radians(incident_angle)
    slice_thickness = 0
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
                                                              slice_thickness,exponent_threshold,v_cm1)
        reflp.append( R[0] + R[2] )
        refls.append( R[1] + R[3] )
        tranp.append( T[0] )
        trans.append( T[1] )
        absop.append( 1.0 - R[0] - R[2] - T[0] )
        absos.append( 1.0 - R[1] - R[3] - T[1] )
    return np.array([reflp,refls]), np.array([tranp,trans]), np.array([absop,absos])

def calculatePowderSpectrum(frequencies_cm1, dielectric, matrix, volume_fraction, method='Maxwell-Garnett'):
    '''
    Calculate the powder IR spectrum of a mixture of spherical dielectric particles in a matrix with a given volume fraction.

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
    float
        Absorption coefficient in cm-1.

    Example
    -------
    ```
    frequencies_cm1 = np.arange( 0, 200, 0.2 )
    matrix = getMaterial('ptfe')
    dielectric = getMaterial('Sapphire')
    method = 'Maxwell-Garnett' 
    volume_fraction = 0.1
    absorption = calculatePowderSpectrum(frequencies_cm1,dielectric, matrix, volume_fraction)
    ```
    '''
    method = method.lower()
    particle_size_mu = 0
    particle_sigma_mu = 0
    matrixPermittivityFunction = matrix.getPermittivityFunction()
    shape = 'Sphere'
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
    for v_cm1 in frequencies_cm1:
        crystalPermittivity = crystalPermittivityFunction(v_cm1)
        result = Calculator.solve_effective_medium_equations(method,volume_fraction,
                       particle_size_mu,particle_sigma_mu,matrixPermittivityFunction,
                       shape,depolarisation,concentration,
                       atr_refractive_index,atr_theta,atr_spolfraction,
                       bubble_vf,bubble_radius,previous_solution_shared,(v_cm1,crystalPermittivity))
        results.append(result)
    # Prepare lists for results
    realPermittivity = []
    imagPermittivity = []
    absorptionCoefficient = []
    molarAbsorptionCoefficient = []
    sp_atr = []
    for v,method,size_mu,size_sigma,shape,data,trace, absorption_coefficient,molar_absorption_coefficient,spatr in results:
         realPermittivity.append(np.real(trace))
         imagPermittivity.append(np.imag(trace))
         absorptionCoefficient.append(absorption_coefficient)
         molarAbsorptionCoefficient.append(molar_absorption_coefficient)
         sp_atr.append(spatr)
    return np.array(absorptionCoefficient)

def powderTest():
    '''
    Test powder calculation
    '''
    frequencies_cm1 = np.arange( 0, 200, 0.2 )
    matrix = getMaterial('ptfe')
    dielectric = getMaterial('Sapphire')
    method = 'Maxwell-Garnett' 
    shape = 'Sphere'
    volume_fraction = 0.1
    absorption = calculatePowderSpectrum(frequencies_cm1,dielectric, matrix, volume_fraction)
    print('Powder absorption results')
    print('      freq', ' absorption', '      r_s', '      t_p', '      t_s', '      a_p', '      a_s')
    for f,a in zip(frequencies_cm1,absorption):
        print('{:10.2f}{:12.5f}'.format(f,a))
    
def singleCrystalTest():
    '''
    Test single crystal calculation
    '''
    # Define the materials
    air = getMaterial('air')
    ptfe = getMaterial('ptfe')
    Sapphire = getMaterial('Sapphire')
    # Prepare the layers
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
    (reflectance, transmittance, absorptance) = calculateSingleCrystalSpectrum(frequencies_cm1,layers,incident_angle, global_azimuthal_angle, method='Scattering matrix')
    # Print single crystal results
    print('Single Crystal Results')
    print('      freq', '      r_p', '      r_s', '      t_p', '      t_s', '      a_p', '      a_s')
    for r_p,r_s,t_p,t_s,a_p,a_s,f in zip(reflectance[0],reflectance[1],transmittance[0],transmittance[1],absorptance[0],absorptance[1],frequencies_cm1):
        print('{:10.2f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}'.format(f,r_p , r_s , t_p , t_s , a_p , a_s ))

def main():
    """
    Execute the main sequence of operations.

    This function runs single crystal and powder tests in sequence.

    Returns
    -------
    None
    """    
    singleCrystalTest()
    powderTest()
    return

if __name__ == '__main__':
    main()