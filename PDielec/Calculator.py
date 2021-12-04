#!/usr/bin/python
#
# Copyright 2015 John Kendrick
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
"""Calculate useful properties, used to be DielectricConstant.py """
from __future__ import print_function
import math
import sys
import os
import cmath
import random
import numpy as np
import scipy.optimize as sc
import PDielec.GTMcore as GTM
import string
from   PDielec.Constants import PI, d2byamuang2, speed_light_si, wavenumber
from   scipy.integrate import trapz
from   scipy.stats import lognorm

#
# Modify the crossover used in the PyMieScatt Mie routines
# The Mie routine is taken from PyMieScatt by B. Sumlin and can be found on github
#
from   PDielec import Mie
Mie.crossover = 0.01

points_on_sphere = None


def initialise_unit_tensor():
    '''Initialise a 3x3 tensor, the argument is a list of 3 real numbers for the diagonals, the returned tensor is an array'''
    x = np.zeros((3, 3), dtype=np.float)
    x[0, 0] = 1.0
    x[1, 1] = 1.0
    x[2, 2] = 1.0
    return x

def initialise_complex_diagonal_tensor(reals):
    '''Initialise a 3x3 tensor, the argument is a list of 3 real numbers for the diagonals, the returned tensor is an array'''
    x = np.zeros((3, 3), dtype=np.complex)
    x[0, 0] = reals[0]
    x[1, 1] = reals[1]
    x[2, 2] = reals[2]
    return x

def initialise_diagonal_tensor(reals):
    '''Initialise a 3x3 tensor, the argument is a list of 3 real numbers for the diagonals, the returned tensor is an array'''
    x = np.zeros((3, 3), dtype=np.float)
    x[0, 0] = reals[0]
    x[1, 1] = reals[1]
    x[2, 2] = reals[2]
    return x

def calculate_distance(a,b):
    '''Calculate the distance between a and b'''
    d = 0.0
    for c1, c2 in zip(a,b):
        d += (c2 - c1)*(c2 - c1)
    return math.sqrt(d)

def initialise_random_tensor(scale):
    '''Initialise a 3x3 complex tensor, the argument gives the maximum absolute value'''
    print("Error initialise_random_tensor not working", file=sys.stderr)
    exit(1)
    return

def initialise_sphere_depolarisation_matrix():
    '''Initialise a 3x3 tensor with the sphere depolarisation matrix, returns a tensor'''
    athird = 1.0 / 3.0
    tensor = initialise_diagonal_tensor([athird, athird, athird])
    tensor = tensor / np.trace(tensor)
    return tensor

def initialise_plate_depolarisation_matrix(normal):
    '''Initialise a 3x3 tensor with the plate depolarisation matrix, returns a tensor'''
    normal = normal / np.linalg.norm(normal)
    tensor = np.outer(normal, normal)
    tensor = tensor / np.trace(tensor)
    return tensor

def initialise_needle_depolarisation_matrix(unique):
    '''Initialise a 3x3 tensor with the needle depolarisation matrix'''
    # unique is the unique direction of the needle
    # depolarisation matrix is therefore half the sum of the
    # two outer products of the directions perpendicular to it
    unique = unique / np.linalg.norm(unique)
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])
    xdotn = np.dot(x, unique)
    ydotn = np.dot(y, unique)
    zdotn = np.dot(z, unique)
    absxdotn = np.abs(xdotn)
    absydotn = np.abs(ydotn)
    abszdotn = np.abs(zdotn)
    # choose the direction with the smallest projection along the needle
    # then project out any needle direction
    if absxdotn <= absydotn and absxdotn <= abszdotn:
        dir1 = x - (xdotn * unique)
    elif absydotn <= absxdotn and absydotn <= abszdotn:
        dir1 = y - (ydotn * unique)
    else:
        dir1 = z - (zdotn * unique)
    dir1 = dir1/np.linalg.norm(dir1)
    # now find the orthogonal direction
    dir2 = np.cross(unique, dir1)
    dir2 = dir2/np.linalg.norm(dir2)
    # compute the complex tensors from the outer product of each direction
    tensor = np.outer(dir1, dir1) + np.outer(dir2, dir2)
    tensor = tensor / np.trace(tensor)
    return tensor

def initialise_ellipsoid_depolarisation_matrix(unique, aoverb):
    '''Initialise a 3x3 tensor with the ellipsoid depolarisation matrix'''
    # unique is the unique direction of the ellipsoid
    # a and b are the sizes of the ellipsoid along the unique axis (a) and perpendicular to it b
    # if a > b then we have a prolate ellipsoid
    unique = unique / np.linalg.norm(unique)
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])
    xdotn = np.dot(x, unique)
    ydotn = np.dot(y, unique)
    zdotn = np.dot(z, unique)
    absxdotn = np.abs(xdotn)
    absydotn = np.abs(ydotn)
    abszdotn = np.abs(zdotn)
    # choose the direction with the smallest projection along the ellipsoid
    # then project out any ellipsoid direction
    if absxdotn <= absydotn and absxdotn <= abszdotn:
        dir1 = x - (xdotn * unique)
    elif absydotn <= absxdotn and absydotn <= abszdotn:
        dir1 = y - (ydotn * unique)
    else:
        dir1 = z - (zdotn * unique)
    dir1 = dir1/np.linalg.norm(dir1)
    # now find the orthogonal direction
    dir2 = np.cross(unique, dir1)
    dir2 = dir2/np.linalg.norm(dir2)
    bovera = 1.0 / aoverb
    small = 1.0E-8
    if bovera < 1.0-small:
        e = math.sqrt(1.0  - bovera*bovera)
        nz = (1 - e*e) * (np.log(((1+e) / (1-e))) - 2*e) / (2*e*e*e)
    elif bovera > 1.0+small:
        e = math.sqrt(bovera*bovera - 1.0)
        nz = (1 + e*e) * (e - np.arctan(e)) / (e*e*e)
    else:
        nz = 1.0/3.0
    nxy = (1 - nz) * 0.5
    #
    # compute the tensors from the outer product of each direction
    tensor = nz*np.outer(unique, unique) + nxy*np.outer(dir1, dir1) + nxy*np.outer(dir2, dir2)
    tensor = tensor / np.trace(tensor)
    return tensor

def fibonacci_sphere(samples=1,randomize=True):
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points.append([x,y,z])

    return points

def ionic_permittivity(mode_list, oscillator_strengths, frequencies, volume):
    """Calculate the low frequency permittivity or zero frequency permittivity
       oscillator_strengths are in atomic units
       frequencies are in atomic units
       volume is in atomic units"""
    permittivity = np.zeros((3, 3))
    for imode in mode_list:
        permittivity = permittivity + oscillator_strengths[imode] / (frequencies[imode] * frequencies[imode])
    # end for
    return permittivity * (4 * PI / volume)

def infrared_intensities(oscillator_strengths):
    """Calculate the IR intensities from the trace of the oscillator strengths,
       The intensities are returned in units of (D/A)^2/amu"""
    # Each mode has a 3x3 oscillator strength
    nmodes = np.size(oscillator_strengths, 0)
    intensities = np.zeros(nmodes)
    for imode, strength in enumerate(oscillator_strengths):
        # We calculate the intensity from the trace of the strengths
        intensities[imode] = intensities[imode] + strength[0, 0] + strength[1, 1] + strength[2, 2]
    # end for
    maxintensity = np.max(intensities)
    convert = maxintensity
    # convert the intensities to Castep units (D/A)**2/amu
    convert = d2byamuang2
    intensities = intensities / convert
    return intensities

def longitudinal_modes(frequencies, normal_modes, born_charges, masses, epsilon_inf, volume, qlist, reader):
    """Apply the nonanalytic correction to the dynamical matrix and calculate the LO frequencies
       frequencies are the frequencies (f) in atomic units
       normal_modes are the mass weighted normal modes (U)
       born_charges are the born charges (Z) stored as
          [Z1x Z1y Z1z] [Z2x Z2y Z2z] [Z3x Z3y Z3z]]
          where 1, 2, 3 are the directions of the field and x, y, z are the coordinates of the atom
       qlist is a list of direction vectors
       The subroutine returns a list of (real) frequencies in atomic units
       Any imaginary frequencies are set to 0
       If projection was requested in the reader, the correction is modified ensure translational invariance"""
    # Use a sqrt that returns a complex number
    # from numpy.lib.scimath import sqrt
    # First step is to reconstruct the dynamical matrix (D) from the frequencies and the eigenvectors
    # f^2 = UT . D . U
    # and U is a hermitian matrix so U-1 = UT
    # D = (UT)-1 f^2 U-1 = U f UT
    # Construct UT from the normal modes
    n = np.size(normal_modes, 0)
    m = np.size(normal_modes, 1)*3
    UT = np.zeros((n, m))
    for imode, mode in enumerate(normal_modes):
        n = 0
        for atom in mode:
            # in python the first index is the row of the matrix, the second is the column
            UT[imode, n+0] = atom[0]
            UT[imode, n+1] = atom[1]
            UT[imode, n+2] = atom[2]
            n = n + 3
        # end for atom
    # end for imode
    # zero the nonanalytical correction
    Wm = np.zeros((n, n))
    # convert the frequencies^2 to a real diagonal array
    # Warning we have to make sure the sign is correct here
    f2 = np.diag(np.sign(frequencies)*np.real(frequencies*frequencies))
    Dm = np.dot(np.dot(UT.T, f2), UT)
    # Make sure the dynamical matrix is real
    Dm = np.real(Dm)
    # Find its eigenvalues
    eig_val, eig_vec = np.linalg.eigh(Dm)
    # Store the results for returning to the main program
    results = []
    # Loop over q values
    for q in qlist:
        # Now calculate the nonanalytic part
        constant = 4.0 * PI / (np.dot(np.dot(q, epsilon_inf), q) * volume)
        # Loop over atom a
        for a, za in enumerate(born_charges):
            # atom is the atom index
            # born contains the polarisability tensor [z1x z1y z1z] [z2x z2y z2z] [z3x z3y z3z]]
            # where 1, 2, 3 are the directions of the field and x, y, z are the coordinates of the atom
            za = np.dot(q, za)
            # Loop over atom b
            for b, zb in enumerate(born_charges):
                zb = np.dot(q, zb)
                terms = np.outer(za, zb) * constant / math.sqrt(masses[a]*masses[b])
                i = a*3
                for termi in terms:
                    j = b*3
                    for term in termi:
                        Wm[i, j] = term
                        j = j + 1
                    # end for term
                    i = i + 1
                # end for i
            # end loop over b
        # end loop over a
        # Construct the full dynamical matrix with the correction
        Dmq = Dm + Wm
        # If projection was requested when the matrix was read, project out translation
        if reader.eckart:
            reader.project(Dmq)
        eig_val, eig_vec = np.linalg.eigh(Dmq)
        # If eig_val less than zero we set it to zero
        values = []
        for eig in eig_val:
            if eig >= 0:
                val = math.sqrt(eig)
            else:
                val = -math.sqrt(-eig)
            values.append(val)
        # end of for eig
        # Sort the eigen values in ascending order and append to the results
        results.append(np.sort(values))
    # end loop over q
    return results

def oscillator_strengths(normal_modes, born_charges):
    """Calculate oscillator strengths from the normal modes and the born charges
       normal_modes are in the mass weighted coordinate system and normalised
       born charges are in electrons, so atomic units"""
    # Each mode has a 3x3 oscillator strength
    nmodes = np.size(normal_modes, 0)
    oscillator_strengths = np.zeros((nmodes, 3, 3))
    for imode, mode in enumerate(normal_modes):
        # We calculate the dipole induced by displacement of each atom along the normal mode
        z_imode = np.zeros(3)
        for atom, born in enumerate(born_charges):
            # atom is the atom index
            # born contains the polarisability tensor [a1x a1y a1z] [a2x a2y a2z] [a3x a3y a3z]]
            # where 1, 2, 3 are the directions of the field and x, y, z are the coordinates of the atom
            z_imode = z_imode + np.dot(born, mode[atom])  # the displacement is an array [x, y, z]
        # end for
        # The oscillator strength matrix is the outer product of z
        oscillator_strengths[imode] = np.outer(z_imode, z_imode)
    # end for
    return oscillator_strengths

def normal_modes(masses, mass_weighted_normal_modes):
    """ Transform from mass weighted coordinates to xyz. Note this returns an array object.
        The returned normal modes have NOT been renormalised.
        The input masses are in atomic units
        the output normal modes are in atomic units """
    list_m = []
    normal_modes = np.zeros_like(mass_weighted_normal_modes)
    nions = np.size(masses)
    for a in range(nions):
        x = 1.0 / math.sqrt(masses[a])
        atom = [x, x, x]
        list_m.append(atom)
    # end of loop of ions
    array_m = np.array(list_m)
    for index, mode in enumerate(mass_weighted_normal_modes):
        normal_modes[index] = mode * array_m
    return normal_modes

def project_field(shape, shape_data, projection, efield):
    """Take the field directions in efield and apply shape projection."""
    if projection == "random":
        return np.array(efield)
    elif projection == "parallel":
        if shape == "sphere":
            return np.array(efield)
        else:
            return np.array(shape_data)
    elif projection == "perpendicular":
        if shape == "sphere":
            return np.array(efield)
        else:
            proj_field = []
            data = np.array(shape_data)
            for field in efield:
                proj = field - (np.dot(data, field) * data)
                size = np.linalg.norm(proj_field)
                if size > 0.0001:
                    proj_field.append(proj)
                # end if
            # end for
            return np.array(proj_field)
    else:
        print("Error in project_field, projection unkown: ", projection, file=sys.stderr)
        exit(1)
    return

def rogridgues_rotations(efield):
    """Take the field directions in efield and use each direction to calculate a random rotation about that axis.
       Use the field (which a random unit vector in xyz space) to generate an orthogonal rotation matrix
       Rodrigues rotation formula A = I3. cos(theta) + (1-cos(theta)) e . eT + ex sin(theta0
       I3 is a unit matrix
       e is the direction
       ex is the cross product matrix
       ( 0  -e3   e2)
       ( e3  0   -e1)
       (-e2  e1   0 )
       Input field is real
       Output is a list of rotations
    """
    rotations = []
    for field in efield:
        # Calculate a random angle between 0 and 180
        theta = PI*np.random.rand()
        cos   = np.cos(theta)
        sin   = np.sin(theta)
        rotation = np.zeros((3, 3))
        e1 = field[0]
        e2 = field[1]
        e3 = field[2]
        rotation[0, 0] = cos + (1-cos)*e1*e1
        rotation[0, 1] = (1-cos)*e2*e1 + sin*e3
        rotation[0, 2] = (1-cos)*e3*e1 - sin*e2
        rotation[1, 0] = (1-cos)*e1*e2 - sin*e3
        rotation[1, 1] = cos + (1-cos)*e2*e2
        rotation[1, 2] = (1-cos)*e3*e2 + sin*e1
        rotation[2, 0] = (1-cos)*e1*e3 + sin*e2
        rotation[2, 1] = (1-cos)*e2*e3 - sin*e1
        rotation[2, 2] = cos + (1-cos)*e3*e3
        rotations.append(rotation)
    return rotations

def absorption_from_mode_intensities(f, modes, frequencies, sigmas, intensities):
    """Calculate the absorption from the frequencies and intensities using a Lorentzian
       f is the frequency of the absorption in cm-1
       modes are a list of the modes
       frequencies(cm-1), sigmas(cm-1) and intensities (D2/A2/amu)
       The number 4225.6 converts the units of D2/A2/amu to L mole-1 cm-1 cm-1
       The output from this is the molar absorption coefficient at f, in L/mol/cm"""
    absorption = 0.0
    for mode in modes:
        v = np.real(frequencies[mode])
        sigma = sigmas[mode]
        icastep = intensities[mode]
        absorption = absorption + 2.0 * 4225.6 * icastep / PI * (sigma / (4.0 * (f - v)*(f - v) + sigma*sigma))
    return absorption

def calculate_size_factor (x):
    """
    Calculate a size effect using Equations 10.38 and 10.39 in Sihvola
    """
    if x < 1.0E-12:
        result = 1.0
    else:
        ix = complex(0,x)
        g1 = 2.0 / 3.0 * ( ( 1.0 + ix ) * np.exp(-ix) - 1.0 )
        g2 = ( 1 + ix - (7.0/15.0) * x*x - complex(0,2.0*x*x*x/15.0) ) * np.exp(-ix) -1.0
        result = 1 - g1 - g2
    return result

def averaged_permittivity(dielectric_medium, crystal_permittivity, shape, L, vf, size):
    """Calculate the effective constant permittivity using the averaged permittivity method
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape
       L is the shapes depolarisation matrix
       size is the dimensionless size parameter for the frequency under consideration (not used)
       The routine returns the effective dielectric constant"""
    effd = vf * crystal_permittivity + (1.0-vf) * dielectric_medium
    trace = np.trace(effd) / 3.0
    effdielec = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    return effdielec

def balan(dielectric_medium, crystal_permittivity, shape, L, vf, size):
    """Calculate the effective constant permittivity using the method of balan
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape
       L is the shapes depolarisation matrix
       size is the dimensionless size parameter for the frequency under consideration (not used)
       The routine returns the effective dielectric constant"""
    unit = initialise_unit_tensor()
    dielecvm1 = (crystal_permittivity - unit)
    deformation  = np.dot(dielectric_medium, np.linalg.inv(dielectric_medium + np.dot(L, (crystal_permittivity - dielectric_medium))))
    effd          = unit + np.dot(deformation, dielecvm1)
    trace = vf * np.trace(effd) / 3.0
    effdielec = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    return effdielec

def spherical_averaged_mie_scattering(dielectric_medium, crystal_permittivity, shape, L, vf, size, size_mu, size_distribution_sigma):
    """Calculate the effective constant permittivity using a Mie scattering approach
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape (NOT USED)
       L is the shapes depolarisation matrix (NOT USED)
       size is the dimensionless size parameter for the frequency under consideration
       size_distribution_sigma is the log normal value of sigma
       vf is the volume fraction of filler
       Mie only works for spherical particles, so shape, and L parameters are ignored
       The anisotropy of the permittivity is accounted for by sampling many directions
       and calculating the scattering in each direction
       The routine returns the effective dielectric constant"""
    global points_on_sphere
    # define i as a complex number
    i = complex(0,1)
    # We need to taken account of the change in wavelength and the change in size parameter due to the
    # None unit value of the dielectric of the embedding medium
    # The size parameter is 2pi r / lambda
    # The effective lambda in the supporting medium is lambda / sqrt(emedium)
    # Where the refractive index is taken to be sqrt(emedium) (non magnetic materials)
    emedium = np.trace(dielectric_medium) / 3.0
    refractive_index_medium = np.real(np.sqrt(emedium))
    lambda_vacuum_mu = 2 * PI * size_mu / size
    wavelength_nm = lambda_vacuum_mu * 1000 / refractive_index_medium
    radius_nm = size_mu * 1000
    # The wavevector in nm-1
    k_nm = 2 * PI / wavelength_nm
    # volume of a particle in nm^3
    V_nm = 4.0/3.0 * PI * radius_nm * radius_nm * radius_nm
    # Number density of particles (number / nm^3)
    N_nm = vf / V_nm
    # If there is a size distribution set up to use it
    if size_distribution_sigma:
        lower,upper = lognorm.interval(0.9999,size_distribution_sigma,scale=size_mu)
        numberOfBins = 40
        dp = np.logspace(np.log(lower),np.log(upper),numberOfBins,base=np.e)
        # The definitions used are confusing;
        # dp is the log of the variable
        # s is the standard deviation (shape function) of the log of the variate
        # scale is the mean of the underlying normal distribution
        ndp = lognorm.pdf(dp,s=size_distribution_sigma,scale=size_mu)
        #print("DP",dp)
        #print("NDP",ndp)
        #print("Upper lower",lower,upper)
        #print("Size_MU",size_mu)
    # Calculate the sampling points on a fibonacci sphere
    if points_on_sphere == None:
        points_on_sphere = fibonacci_sphere(samples=500,randomize=True)
    trace = 0.0
    # Now take the average of each direction on the sphere
    for point in points_on_sphere:
        point = np.array(point)
        rotated_dielec = np.dot(point, np.dot(point, crystal_permittivity))
        refractive_index = calculate_refractive_index_scalar(rotated_dielec) / refractive_index_medium
        # print('refractive_index', refractive_index)
        # print('refractive_index_medium', refractive_index_medium)
        # print('rotated_dielec', rotated_dielec)
        if size_distribution_sigma:
            # Calculate the integral of the forward scattering factors over the distribution
            s1_factors = []
            for r in dp:
                # The size parameter is 2pi r / lambda
                x = 2 * PI * r / lambda_vacuum_mu
                # Calculate the S1 and S2 scattering factors, and store in a list
                s1,s2 = Mie.MieS1S2(refractive_index, x*refractive_index_medium, 1)
                s1_factors.append(s1)
            # Now integrate
            s1 = trapz(s1_factors*ndp,dp)
            normal = trapz(ndp,dp)
            #mean = trapz(ndp*dp,dp)
            #true_mean = np.exp( np.log(size_mu) + size_distribution_sigma*size_distribution_sigma/2.0)
            #v_cm1 = 1.0E4/lambda_vacuum_mu
            #print("Frequency,normal,mean",v_cm1,normal,true_mean,mean)
            if np.abs(normal - 1.0) > 1.0E-2:
                print("Warning integration of log-normal distribution in error", normal)
                print("Stopping calculation - likely problem is too large a sigma for log-normal distribution")
                effdielec = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                return effdielec
        else:
            # Calculate the scattering factors at 0 degrees
            #jk print("refractive_index, size, refractive_index_medium", refractive_index, size, refractive_index_medium)
            s1,s2 = Mie.MieS1S2(refractive_index, size*refractive_index_medium, 1)
        # See van de Hulst page 129, 130
        # Refractive index of material is
        # the sign of the imaginary component has changed for compatibility with MG/Bruggeman
        refractive_index = refractive_index_medium * ( 1.0 + i * s1 * 2 * PI * N_nm / ( k_nm * k_nm * k_nm ) )
        trace += refractive_index
    # return an isotropic tensor
    trace = trace / len(points_on_sphere)
    eff = trace * trace
    effdielec = np.array([[eff, 0, 0], [0, eff, 0], [0, 0, eff]])
    #jk print(radius_nm, lambda_vacuum_mu*1000.0, qext,qsca,qabs,g,qpr,qback,qratio,np.real(s1),np.imag(s1),np.real(trace),np.imag(trace),np.real(eff),np.imag(eff))
    #print ("radius_nm, eff", radius_nm, eff)
    return effdielec


def mie_scattering(dielectric_medium, crystal_permittivity, shape, L, vf, size, size_mu, size_distribution_sigma):
    """Calculate the effective constant permittivity using a Mie scattering approach
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape (NOT USED)
       L is the shapes depolarisation matrix (NOT USED)
       size is the dimensionless size parameter for the frequency under consideration
       size_distribution_sigma is the log normal value of sigma
       vf is the volume fraction of filler
       In this method the MG method is used to calculate the averaged effective permittivity
       Then the permittivity of the isodtropic sphere that would give the same average permittivity is calculated
       Then the Mie scattering of that sphere is calculated
       The routine returns the effective dielectric constant"""
    #
    # Calculate the MG permittivity with no size effects
    #
    mg_permittivity = maxwell(dielectric_medium, crystal_permittivity, shape, L, vf, 0.0000001)
    # Use scalar quantities to calculate the dielectric constant of the equivalent isotropic sphere
    ef = np.trace(mg_permittivity) / 3.0
    em = np.trace(dielectric_medium) / 3.0
    # print('crystal_permittivity',crystal_permittivity)
    # print('Maxwell',ef)
    # print('EMedium',em)
    # Calculate the permittivity of an isotropic sphere that has the same effective permittivity
    einclusion = ( -3*vf*em*em - (ef - em)*em*(2+vf) ) / ((ef-em)*(1-vf) - 3*vf*em)
    # print('E inclusion',einclusion)
    dielecv = einclusion*np.eye(3)
    # print('New dielecv',dielecv)
    # define i as a complex number
    i = complex(0,1)
    # We need to taken account of the change in wavelength and the change in size parameter due to the
    # None unit value of the dielectric of the embedding medium
    # The size parameter is 2pi r / lambda
    # The effective lambda in the supporting medium is lambda / sqrt(emedium)
    # Where the refractive index is taken to be sqrt(emedium) (non magnetic materials)
    emedium = np.trace(dielectric_medium) / 3.0
    refractive_index_medium = np.real(np.sqrt(emedium))
    lambda_vacuum_mu = 2 * PI * size_mu / size
    wavelength_nm = lambda_vacuum_mu * 1000 / refractive_index_medium
    radius_nm = size_mu * 1000
    # The wavevector in nm-1
    k_nm = 2 * PI / wavelength_nm
    # volume of a particle in nm^3
    V_nm = 4.0/3.0 * PI * radius_nm * radius_nm * radius_nm
    # Number density of particles (number / nm^3)
    N_nm = vf / V_nm
    # If there is a size distribution set up to use it
    if size_distribution_sigma:
        lower,upper = lognorm.interval(0.9999,size_distribution_sigma,scale=size_mu)
        numberOfBins = 40
        dp = np.logspace(np.log(lower),np.log(upper),numberOfBins,base=np.e)
        # The definitions used are confusing;
        # dp is the log of the variable
        # s is the standard deviation (shape function) of the log of the variate
        # scale is the mean of the underlying normal distribution
        ndp = lognorm.pdf(dp,s=size_distribution_sigma,scale=size_mu)
        #print("DP",dp)
        #print("NDP",ndp)
        #print("Upper lower",lower,upper)
        #print("Size_MU",size_mu)
    refractive_index = calculate_refractive_index_scalar(einclusion) / refractive_index_medium
    #jk print('refractive_index', refractive_index)
    #jk print('refractive_index_medium', refractive_index_medium)
    #jk print('einclusion', einclusion)
    if size_distribution_sigma:
        # Calculate the integral of the forward scattering factors over the distribution
        s1_factors = []
        for r in dp:
            # The size parameter is 2pi r / lambda
            x = 2 * PI * r / lambda_vacuum_mu
            # Calculate the S1 and S2 scattering factors, and store in a list
            s1,s2 = Mie.MieS1S2(refractive_index, x*refractive_index_medium, 1)
            s1_factors.append(s1)
        # Now integrate
        s1 = trapz(s1_factors*ndp,dp)
        normal = trapz(ndp,dp)
        #mean = trapz(ndp*dp,dp)
        #true_mean = np.exp( np.log(size_mu) + size_distribution_sigma*size_distribution_sigma/2.0)
        #v_cm1 = 1.0E4/lambda_vacuum_mu
        #print("Frequency,normal,mean",v_cm1,normal,true_mean,mean)
        if np.abs(normal - 1.0) > 1.0E-2:
            print("Warning integration of log-normal distribution in error", normal)
            print("Stopping calculation - likely problem is too large a sigma for log-normal distribution")
            effdielec = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            return effdielec
    else:
        # Calculate the scattering factors at 0 degrees
        #jk print("refractive_index, size, refractive_index_medium", refractive_index, size, refractive_index_medium)
        s1,s2 = Mie.MieS1S2(refractive_index, size*refractive_index_medium, 1)
    # See van de Hulst page 129, 130
    # Refractive index of material is
    # the sign of the imaginary component has changed for compatibility with MG/Bruggeman
    refractive_index = refractive_index_medium * ( 1.0 + i * s1 * 2 * PI * N_nm / ( k_nm * k_nm * k_nm ) )
    eff = refractive_index * refractive_index
    effdielec = np.array([[eff, 0, 0], [0, eff, 0], [0, 0, eff]])
    #jk print(radius_nm, lambda_vacuum_mu*1000.0, qext,qsca,qabs,g,qpr,qback,qratio,np.real(s1),np.imag(s1),np.real(refractive_index),np.imag(refractive_index),np.real(eff),np.imag(eff))
    #print ("radius_nm, eff", radius_nm, eff)
    return effdielec

def anisotropic_mie_scattering(dielectric_medium, crystal_permittivity, shape, L, vf, size, size_mu, size_distribution_sigma):
    """Calculate the effective constant permittivity using a Mie scattering approach
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape (NOT USED)
       L is the shapes depolarisation matrix (NOT USED)
       size is the dimensionless size parameter for the frequency under consideration
       size_distribution_sigma is the log normal value of sigma
       vf is the volume fraction of filler
       Mie only works for spherical particles, so shape, and L parameters are ignored
       The routine returns the effective dielectric constant"""
    # define i as a complex number
    i = complex(0,1)
    # We need to taken account of the change in wavelength and the change in size parameter due to the
    # None unit value of the dielectric of the embedding medium
    # The size parameter is 2pi r / lambda
    # The effective lambda in the supporting medium is lambda / sqrt(emedium)
    # Where the refractive index is taken to be sqrt(emedium) (non magnetic materials)
    emedium = np.trace(dielectric_medium) / 3.0
    refractive_index_medium = np.real(np.sqrt(emedium))
    lambda_vacuum_mu = 2 * PI * size_mu / size
    wavelength_nm = lambda_vacuum_mu * 1000 / refractive_index_medium
    radius_nm = size_mu * 1000
    # To account for anisotropy we diagonalise the real part of the dielectric matrix and transform
    # the full matrix with the eigenvectors, U
    # Find U and E, such that UT. D. U = E (where D is the real part of crystal_permittivity)
    E,U = np.linalg.eig(np.real(crystal_permittivity))
    # Transform the full dielectric matrix
    rotated_dielec = np.matmul(U.T,np.matmul(crystal_permittivity,U))
    # The wavevector in nm-1
    k_nm = 2 * PI / wavelength_nm
    # volume of a particle in nm^3
    V_nm = 4.0/3.0 * PI * radius_nm * radius_nm * radius_nm
    # Number density of particles (number / nm^3)
    N_nm = vf / V_nm
    # If there is a size distribution set up to use it
    if size_distribution_sigma:
        lower,upper = lognorm.interval(0.9999,size_distribution_sigma,scale=size_mu)
        numberOfBins = 40
        dp = np.logspace(np.log(lower),np.log(upper),numberOfBins,base=np.e)
        # The definitions used are confusing;
        # dp is the log of the variable
        # s is the standard deviation (shape function) of the log of the variate
        # scale is the mean of the underlying normal distribution
        ndp = lognorm.pdf(dp,s=size_distribution_sigma,scale=size_mu)
        #print("DP",dp)
        #print("NDP",ndp)
        #print("Upper lower",lower,upper)
        #print("Size_MU",size_mu)
    # We are now going to ignore any off-diagonal elements
    trace = 0.0
    # Now take the average of each direction
    for index in [0,1,2]:
        refractive_index = calculate_refractive_index_scalar(rotated_dielec[index,index]) / refractive_index_medium
        #jk print('refractive_index', refractive_index)
        #jk print('refractive_index_medium', refractive_index_medium)
        #jk print('rotated_dielec', rotated_dielec[index,index])
        if size_distribution_sigma:
            # Calculate the integral of the forward scattering factors over the distribution
            s1_factors = []
            for r in dp:
                # The size parameter is 2pi r / lambda
                x = 2 * PI * r / lambda_vacuum_mu
                # Calculate the S1 and S2 scattering factors, and store in a list
                s1,s2 = Mie.MieS1S2(refractive_index, x*refractive_index_medium, 1)
                s1_factors.append(s1)
            # Now integrate
            s1 = trapz(s1_factors*ndp,dp)
            normal = trapz(ndp,dp)
            #mean = trapz(ndp*dp,dp)
            #true_mean = np.exp( np.log(size_mu) + size_distribution_sigma*size_distribution_sigma/2.0)
            #v_cm1 = 1.0E4/lambda_vacuum_mu
            #print("Frequency,normal,mean",v_cm1,normal,true_mean,mean)
            if np.abs(normal - 1.0) > 1.0E-2:
                print("Warning integration of log-normal distribution in error", normal)
                print("Stopping calculation - likely problem is too large a sigma for log-normal distribution")
                effdielec = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                return effdielec
        else:
            # Calculate the scattering factors at 0 degrees
            #jk print("refractive_index, size, refractive_index_medium", refractive_index, size, refractive_index_medium)
            s1,s2 = Mie.MieS1S2(refractive_index, size*refractive_index_medium, 1)
        # See van de Hulst page 129, 130
        # Refractive index of material is
        # the sign of the imaginary component has changed for compatibility with MG/Bruggeman
        refractive_index = refractive_index_medium * ( 1.0 + i * s1 * 2 * PI * N_nm / ( k_nm * k_nm * k_nm ) )
        trace += refractive_index
    # return an isotropic tensor
    trace = trace / 3.0
    eff = trace * trace
    effdielec = np.array([[eff, 0, 0], [0, eff, 0], [0, 0, eff]])
    #jk print(radius_nm, lambda_vacuum_mu*1000.0, qext,qsca,qabs,g,qpr,qback,qratio,np.real(s1),np.imag(s1),np.real(trace),np.imag(trace),np.real(eff),np.imag(eff))
    #print ("radius_nm, eff", radius_nm, eff)
    return effdielec

def maxwell(dielectric_medium, crystal_permittivity, shape, L, vf, size):
    """Calculate the effective constant permittivity using the maxwell garnett method
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape
       L is the shapes depolarisation matrix
       size is the dimensionless size parameter for the frequency under consideration
       vf is the volume fraction of filler
       The routine returns the effective dielectric constant"""
    unit = initialise_unit_tensor()
    emedium = np.trace(dielectric_medium) / 3.0
    # If appropriate calculate a size effect using Equations 10.38 and 10.39 in Sihvola
    size_factor = calculate_size_factor(size)
    # Equation 5.78 in Sihvola
    nalpha = emedium*vf*np.dot((crystal_permittivity - dielectric_medium), np.linalg.inv(dielectric_medium + size_factor * np.dot(L, (crystal_permittivity-dielectric_medium))))
    nalphal = np.dot((nalpha/emedium), L)
    # average the polarisability over orientation
    nalpha = np.trace(nalpha) / 3.0 * unit
    # average the polarisability*L over orientation
    nalphal = np.trace(nalphal) / 3.0 * unit
    polarisation = np.dot(np.linalg.inv(unit - nalphal), nalpha)
    effd         = dielectric_medium + polarisation
    trace = np.trace(effd) / 3.0
    effdielec = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    return effdielec

def maxwell_sihvola(dielectric_medium, crystal_permittivity, shape, L, vf, size):
    """Calculate the effective constant permittivity using the maxwell garnett method
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape
       L is the shapes depolarisation matrix
       size is the dimensionless size parameter for the frequency under consideration
       vf is the volume fraction of filler
       The routine returns the effective dielectric constant"""
    unit = initialise_unit_tensor()
    # Equation 6.29 on page 123 of Sihvola
    # Equation 6.40 gives the averaging over the orientation function
    # See also equation 5.80 on page 102 and equation 4.31 on page 70
    Me = dielectric_medium
    # assume that the medium is isotropic calculate the inverse of the dielectric
    Mem1 = 3.0 / np.trace(Me)
    Mi = crystal_permittivity
    # If appropriate calculate a size effect using Equations 10.38 and 10.39 in Sihvola
    size_factor = calculate_size_factor(size)
    # calculate the polarisability matrix x the number density of inclusions
    nA = vf*np.dot((Mi-Me), np.linalg.inv(unit + (size_factor * Mem1 * np.dot(L, (Mi - Me)))))
    nAL = np.dot((nA), L)
    # average the polarisability over orientation
    nA = np.trace(nA) / 3.0 * unit
    # average the polarisability*L over orientation
    nAL = np.trace(nAL) / 3.0 * Mem1 * unit
    # Calculate the average polarisation factor which scales the average field
    # based on equation 5.80
    # <P> = pol . <E>
    pol = np.dot(np.linalg.inv(unit - nAL), nA)
    # Meff . <E> = Me . <E> + <P>
    # Meff . <E> = Me. <E> + pol . <E>
    # Meff = Me + pol
    effd         = dielectric_medium + pol
    # Average over orientation
    trace = np.trace(effd) / 3.0
    effdielec = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    return effdielec

def coherent(dielectric_medium, crystal_permittivity, shape, L, vf, size, dielectric_apparent):
    """Driver for coherent2 method"""
    for i in range(10):
        dielectric_apparent = 0.1 * dielectric_apparent + 0.9 * coherent2(dielectric_medium, dielectric_apparent, crystal_permittivity, shape, L, vf, size)
    return dielectric_apparent

def coherent2(dielectric_medium, dielectric_apparent, crystal_permittivity, shape, L, vf, size):
    """Calculate the effective constant permittivity using the maxwell garnett method
       dielectric_medium is the dielectric constant tensor of the medium
       crystal_permittivity is the total frequency dielectric constant tensor at the current frequency
       shape is the name of the current shape
       L is the shapes depolarisation matrix
       size is the dimensionless size parameter for the frequency under consideration
       vf is the volume fraction of filler
       The routine returns the effective dielectric constant"""
    unit = initialise_unit_tensor()
    emedium = np.trace(dielectric_medium) / 3.0
    eapparent = np.trace(dielectric_apparent) / 3.0
    # If appropriate calculate a size effect using Equations 10.38 and 10.39 in Sihvola
    size_factor = calculate_size_factor(size)
    # Equation 5.78 in Sihvola
    nalpha = emedium*vf*np.dot((crystal_permittivity - dielectric_medium), np.linalg.inv(dielectric_medium + size_factor * np.dot(L, (crystal_permittivity-dielectric_medium))))
    nalphal = np.dot((nalpha/eapparent), L)
    # average the polarisability over orientation
    nalpha = np.trace(nalpha) / 3.0 * unit
    # average the polarisability*L over orientation
    nalphal = np.trace(nalphal) / 3.0 * unit
    polarisation = np.dot(np.linalg.inv(unit - nalphal), nalpha)
    effd         = dielectric_medium + polarisation
    trace = np.trace(effd) / 3.0
    effdielec = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    return effdielec

def bruggeman_minimise( eps1, eps2, shape, L, f2, size, epsbr):
    """Calculate the effective constant permittivity using the method of bruggeman
       eps1 is the dielectric constant tensor of 1 (The medium)
       eps2 is the dielectric constant tensor of 2 (The inclusion)
       shape is the name of the current shape
       L is the shapes depolarisation matrix
       f2 is the volume fraction of component 2
       size is the dimensionless size parameter for the frequency under consideration
       epsbr is an initial guess at the solution
       The routine returns the effective dielectric constant
       On the application of homogenization formalisms to active dielectric composite materials
       Tom G. Mackay, Akhlesh Lakhtakia """
    f1 = 1.0 - f2
    # we need to fool the optimiser into thinking that it has two real variables
    # in fact the second is imaginary and reconstructed in the _brug_minimise routine
    trace = np.trace(epsbr) / 3.0
    variables = np.array([np.real(trace), np.log(1.0 + np.abs(np.imag(trace)))])
    options = {'xtol': 1.0e-4,
               'ftol': 1.0E-4}
    sol = sc.minimize(_brug_minimise_tensor, variables, method='Powell', args=(eps1, eps2, shape, L, f1, size), options=options)
    if not sol.success:
        print("A Bruggeman solution was not found at this frequency")
    variables = sol.x
    # transform the imaginary variable back
    trace = complex(variables[0], np.exp(variables[1])-1.0)
    epsbr = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    return epsbr

def bruggeman_iter( eps1, eps2, shape, L, f2, size, epsbr):
    """Calculate the effective constant permittivity using the method of bruggeman
       eps1 is the dielectric constant tensor of 1 (The medium)
       eps2 is the dielectric constant tensor of 2 (The inclusion)
       shape is the name of the current shape
       L is the shapes depolarisation matrix
       f2 is the volume fraction of component 2
       size is the dimensionless size parameter for the frequency under consideration
       epsbr is an initial guess at the solution
       The routine returns the effective dielectric constant
       On the application of homogenization formalisms to active dielectric composite materials
       Tom G. Mackay, Akhlesh Lakhtakia """
    f1 = 1.0 - f2
    # perform an iteration
    converged = False
    niters = 0
    while not converged:
        niters += 1
        epsbr, error = _brug_iter_error(epsbr, eps1, eps2, shape, L, f1, size)
        if abs(error) < 1.0E-8:
            converged = True
        if niters > 3000:
            print("Bruggeman iterations failed, error=", error)
            converged = True
    epsbr = average_tensor(epsbr)
    return epsbr

def average_tensor(t):
    """Return the averaged tensor"""
    a = np.trace(t) / 3.0
    return np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])

def _brug_minimise_scalar(variables, eps1, eps2, shape, L, f1, size):
    """Bruggeman method using scalar quantities"""
    # unpack the complex number from the variables
    # two things going on here.
    # 1. the two variables refer to the real and imaginary components
    # 2. we require the imaginary component to be positive
    trace = complex(variables[0], np.exp(variables[1])-1.0)
    epsbr = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    f2 = 1.0 - f1
    # If appropriate calculate a size effect using Equations 10.38 and 10.39 in Sihvola
    size_factor = calculate_size_factor(size)
    b1 = np.dot(L, (eps1 - epsbr))
    b2 = size_factor * np.dot(L, (eps2 - epsbr))
    tb1 = np.trace(b1)/3.0
    tb2 = np.trace(b2)/3.0
    ta1 = 1.0/(1.0 + tb1)
    ta2 = 1.0/(1.0 + tb2)
    c1 = eps1-epsbr
    c2 = eps2-epsbr
    tc1 = np.trace(c1)/3.0
    tc2 = np.trace(c2)/3.0
    # alpha1 and 2 are the polarisabilities of 1 and 2 in the effective medium
    talpha1 = tc1 * ta1
    talpha2 = tc2 * ta2
    error = f1*talpha1 + f2*talpha2
    error = np.abs(error.conjugate() * error)
    # Nasty issue in the powell method, the convergence on tol is given
    # relative to the solution (0.0).  Only a small number is added.
    # So we shift the solution by 1.0, the tol is now relative to 1.0
    return 1.0+error

def _brug_minimise_tensor(variables, eps1, eps2, shape, L, f1, size):
    """Bruggeman method using tensor quantities"""
    # unpack the complex number from the variables
    # two things going on here.
    # 1. the two variables refer to the real and imaginary components
    # 2. we require the imaginary component to be positive
    trace = complex(variables[0], np.exp(variables[1])-1.0)
    epsbr = np.array([[trace,  0, 0], [0, trace, 0], [0, 0, trace]])
    f2 = 1.0 - f1
    # If appropriate calculate a size effect using Equations 10.38 and 10.39 in Sihvola
    size_factor = calculate_size_factor(size)
    b1 = np.dot(L, (eps1 - epsbr))
    b2 = size_factor * np.dot(L, (eps2 - epsbr))
    b1 = average_tensor(b1)
    b2 = average_tensor(b2)
    a1 = np.linalg.inv(epsbr + b1)
    a2 = np.linalg.inv(epsbr + b2)
    c1 = eps1-epsbr
    c2 = eps2-epsbr
    # c1 = average_tensor(eps1-epsbr)
    # c2 = average_tensor(eps2-epsbr)
    alpha1 = np.dot(c1, a1)
    alpha2 = np.dot(c2, a2)
    alpha1 = average_tensor(alpha1)
    alpha2 = average_tensor(alpha2)
    error = f1*alpha1 + f2*alpha2
    fr  = 0.0
    fi  = 0.0
    for c in error:
        for f in c:
            fr += np.real(f)**2
            fi += np.imag(f)**2
    error = np.linalg.norm(error)
    # Nasty issue in the powell method, the convergence on tol is given
    # relative to the solution (0.0).  Only a small number is added.
    # So we shift the solution by 1.0, the tol is now relative to 1.0
    return 1.0+error

def _brug_iter_error(epsbr, eps1, eps2, shape, L, f1, size):
    """Routine to calculate the error in the Bruggeman method"""
    f2 = 1.0 - f1
    # If appropriate calculate a size effect using Equations 10.38 and 10.39 in Sihvola
    size_factor = calculate_size_factor(size)
    leps1 = np.dot(L, (eps1 - epsbr))
    leps2 = size_factor * np.dot(L, (eps2 - epsbr))
    leps1 = average_tensor(leps1)
    leps2 = average_tensor(leps2)
    a1 = np.linalg.inv(epsbr + leps1)
    a2 = np.linalg.inv(epsbr + leps2)
    # alpha1 and 2 are the polarisabilities of 1 and 2 in the effective medium
    eps1av = average_tensor(eps1)
    eps2av = average_tensor(eps2)
    alpha1 = np.dot((eps1av-epsbr), a1)
    alpha2 = np.dot((eps2av-epsbr), a2)
    # the error or residual matrix should be zero for a bruggeman solution
    error = f1*alpha1 + f2*alpha2
    error = np.linalg.norm(error)
    m1 = f1*np.dot(eps1, a1)+f2*np.dot(eps2, a2)
    m2 = np.linalg.inv(f1*a1 + f2*a2)
    damp = 0.0
    epsbr = (1.0 - damp)*np.dot(m1, m2) + damp*epsbr
    trace = np.trace(epsbr) / 3.0
    epsbr = np.array([[trace, 0, 0], [0, trace, 0], [0, 0, trace]])
    return epsbr, error

def calculate_refractive_index(dielectric, debug=False):
    ''' Calculate the refractive index from the dielectric constant.
        Calculate the trace of the dielectric and calculate both square roots.
        The choose the root with the largest imaginary component This obeys the Konig Kramer requirements'''
    trace = np.trace(dielectric)/3.0
    solution = calculate_refractive_index_scalar(trace, debug)
    return solution

def calculate_refractive_index_scalar(dielectric_scalar, debug=False):
    ''' Calculate the refractive index from the dielectric constant.
        Calculate the trace of the dielectric and calculate both square roots.
        The choose the root with the largest imaginary component This obeys the Konig Kramer requirements'''
    solution1 = np.sqrt(dielectric_scalar)
    r, phase = cmath.polar(solution1)
    solution2 = cmath.rect(-r, phase)
    imag1 = np.imag(solution1)
    imag2 = np.imag(solution2)
    if imag1 > imag2:
        solution = solution1
    else:
        solution = solution2
    if np.abs(solution*solution-dielectric_scalar)/(1+np.abs(dielectric_scalar)) > 1.0E-8 or debug:
        print("There is an error in refractive index")
        print("Dielectric = ", dielectric_scalar)
        print("solution*solution = ", solution*solution, np.abs(solution*solution-dielectric_scalar))
        print("solution    = ", solution, solution*solution)
        print("solution1   = ", solution1, solution1*solution1)
        print("solution2   = ", solution2, solution2*solution2)
    return solution

def direction_from_shape(data, reader):
    """ Determine the unique direction of the shape from data
    data may contain a miller indices which defines a surface eg. (1,1,-1)
    or a direction as a miller direction vector eg. [1,0,-1] """
    surface = False
    # original = data
    i = data.find(",")
    commas = False
    if i >= 0:
        commas = True
    if data[0] == "{":
        surface = True
        data = data.replace("{", "")
        data = data.replace("}", "")
    elif data[0] == "(":
        surface = True
        data = data.replace("(", "")
        data = data.replace(")", "")
    elif data[0] == "[":
        surface = False
        data = data.replace("[", "")
        data = data.replace("]", "")
    else:
        print("Error encountered in interpretting the miller surface / vector", data)
        exit(1)
    if commas:
        data = data.replace(",", " ")
        hkl = ([int(f) for f in data.split()])
    else:
        i = 0
        hkl = [0, 0, 0]
        for k in [0, 1, 2]:
            sign = 1
            if data[i] == "-":
                i += 1
                sign = -1
            hkl[k] = sign * int(data[i])
            i += 1
    # end of handling no commas
    if not len(hkl) == 3:
        print("Error encountered in interpretting the miller surface / vector", data)
        exit(1)
    cell = reader.unit_cells[-1]
    if surface:
        direction = cell.convert_hkl_to_xyz(hkl)
    else:
        direction = cell.convert_abc_to_xyz(hkl)
    direction = direction / np.linalg.norm(direction)
    data = data.replace('"', '')
    data = data.replace("", '')
    # if surface:
    #     print("The miller indices for the surface ", original, "has a normal", direction, "in xyz")
    # else:
    #    print("The miller direction ", original, "is ", direction, "in xyz")
    return direction

def solve_effective_medium_equations( 
        method                   ,
        vf                       ,
        size_mu                  ,
        size_distribution_sigma  ,
        dielectric_medium        ,
        shape                    ,
        L                        ,
        concentration            ,
        atrPermittivity          ,
        atrTheta                 ,
        atrSPol                  ,
        bubble_vf                ,
        bubble_radius            ,
        previous_solution_shared ,
        params                  ,
        ):
    # call_parameters is an index into the frequency and dielectric arrays
    # In the case of Bruggeman and coherent we can use the previous result to start the iteration/minimisation
    # However to do this we need some shared memory, this allocated in previous_solution_shared
    v,crystal_permittivity = params
    vau = v * wavenumber
    # convert the size to a dimensionless number which is 2*pi*size/lambda
    lambda_mu = 1.0E4 / (v + 1.0e-12)
    if size_mu < 1.0e-12:
        size_mu = 1.0e-12
    size = 2.0*np.pi*size_mu / lambda_mu
    data = ''
    # Calculate the effect of bubbles in the matrix by assuming they are embedded in an effective medium defined above
    refractive_index = math.sqrt(np.trace(dielectric_medium)/3.0)
    if refractive_index.imag < 0.0:
        refractive_index = refractive_index.conjugate()
    if bubble_vf > 0:
        effdielec,refractive_index = calculate_bubble_refractive_index(v, refractive_index, bubble_vf, bubble_radius)
        dielectric_medium = effdielec
    if method == "balan":
        effdielec = balan(dielectric_medium, crystal_permittivity, shape, L, vf, size)
    elif method == "ap" or method == "averagedpermittivity" or method == "averaged permittivity" or method == "average permittivity":
            effdielec = averaged_permittivity(dielectric_medium, crystal_permittivity, shape, L, vf, size)
    elif method == "maxwell" or method == "maxwell-garnett":
        effdielec = maxwell(dielectric_medium, crystal_permittivity, shape, L, vf, size)
    elif method == "maxwell_sihvola":
        effdielec = maxwell_sihvola(dielectric_medium, crystal_permittivity, shape, L, vf, size)
    elif method == "coherent":
        eff  = previous_solution_shared
        if np.abs(np.trace(eff)) < 1.0e-8:
            eff = maxwell(dielectric_medium, crystal_permittivity, shape, L, vf, size)
        effdielec = coherent(dielectric_medium, crystal_permittivity, shape, L, vf, size, eff)
        previous_solution_shared = effdielec
    elif method == "bruggeman_minimise":
        eff  = previous_solution_shared
        if np.abs(np.trace(eff)) < 1.0e-8:
            eff = maxwell(dielectric_medium, crystal_permittivity, shape, L, vf, size)
        effdielec = bruggeman_minimise(dielectric_medium, crystal_permittivity, shape, L, vf, size, eff)
        previous_solution_shared = effdielec
    elif method == "bruggeman" or method == "bruggeman_iter":
        eff  = previous_solution_shared
        if np.abs(np.trace(eff)) < 1.0e-8:
            eff = maxwell(dielectric_medium, crystal_permittivity, shape, L, vf, size)
        effdielec = bruggeman_iter(dielectric_medium, crystal_permittivity, shape, L, vf, size, eff)
        previous_solution_shared = effdielec
    elif method == "anisotropic-mie":
        effdielec = anisotropic_mie_scattering(dielectric_medium, crystal_permittivity, shape, L, vf, size, size_mu, size_distribution_sigma)
    elif method == "mie":
        effdielec = mie_scattering(dielectric_medium, crystal_permittivity, shape, L, vf, size, size_mu, size_distribution_sigma)
    else:
        print('Unkown dielectric method: {}'.format(method))
        exit(1)
    # Average over all directions by taking the trace
    trace = (effdielec[0, 0] + effdielec[1, 1] + effdielec[2, 2]) / 3.0
    refractive_index = calculate_refractive_index(effdielec)
    #
    # absorption coefficient is calculated from the imaginary refractive index
    # see H.C. van de Hulst Light Scattering by Small Particles , page 267
    # This is different but related to Genzel and Martin Equation 16, Phys. Stat. Sol. 51(1972) 91-
    # I've add a factor of log10(e) because we need to assume a decadic Beer's law
    # units are cm-1
    absorption_coefficient = v * 4*PI * np.imag(refractive_index) * math.log10(math.e)
    # units are cm-1 L moles-1
    molar_absorption_coefficient = absorption_coefficient / concentration / vf
    # calculate the ATR reflectance
    spatr = reflectance_atr(refractive_index,atrPermittivity,atrTheta,atrSPol)
    return v,method,size_mu,size_distribution_sigma,shape,data,trace,absorption_coefficient,molar_absorption_coefficient,spatr

def calculate_bubble_refractive_index(v_cm1, ri_medium, vf, radius_mu):
    """Calculate the scattering from bubbles embedded in a possibly, complex dielectric at v_cm1
       v_cm1 is the frequency in cm-1
       ri_medium is the refractive index of the medium
       vf is the volume fraction of bubbles
       radius_mu is the radius of the bubbles in microns
       The routine returns the effective dielectric constant and the associated refractive index"""
    #
    # We need to taken account of the change in wavelength and the change in size parameter due to the
    # None unit value of the dielectric of the embedding medium
    # The size parameter is 2pi r / lambda
    # The effective lambda in the supporting medium is lambda / sqrt(emedium)
    # Where the refractive index is taken to be sqrt(emedium) (non magnetic materials)
    #
    if v_cm1 > 0:
        lambda_vacuum_nm = 1.0E3*1.0E4 / v_cm1
    else:
        lambda_vacuum_nm = 1.0E99
    # Treat the bubble as though it is air in matrix
    # The effective wave number k = sqrt(emedium)*2pi*v/c (complex!)
    radius_nm = radius_mu * 1000
    # volume of a bubble in nm^3
    V_nm = 4.0/3.0 * PI * radius_nm * radius_nm * radius_nm
    # Number density of bubbles (number / nm^3)
    N_nm = vf / V_nm
    k_nm = waterman_truell_scattering(lambda_vacuum_nm, N_nm, radius_nm, ri_medium)
    # k_nm = foldy_scattering(lambda_vacuum_nm, N_nm, radius_nm, ri_medium)
    ri_medium = k_nm *lambda_vacuum_nm / (2*PI)
    eff_medium = ri_medium * ri_medium
    effdielec = np.array([[eff_medium, 0, 0], [0, eff_medium, 0], [0, 0, eff_medium]])
    return effdielec,ri_medium

def foldy_scattering(lambda_vacuum_nm, N_nm,radius_nm,ri_medium):
    #
    # Solve the foldy equation for scattering of an air bubble embedded in a lossy medium
    # N_nm is the number density of bubbles in the volume nm^3
    # k_nm is the wavenumber of the incoming wave (complex)
    # radius_nm is the radius of the bubble
    # ri_medium is the refractive index of the medium the bubble is in
    #
    #
    k_nm = 2*PI*ri_medium/lambda_vacuum_nm
    # The size parameter is now also complex and dimensionless
    size = k_nm*radius_nm
    refractive_index = 1.0 / ri_medium
    # Calculate the forward and backward scattering amplitude
    s10,s20 = Mie.MieS1S2(refractive_index, size*ri_medium, 1)
    i = complex(0,1)
    f0 = i * s10 / k_nm
    new_k = np.sqrt( k_nm*k_nm + 4*PI*N_nm*f0 )
    if new_k.imag < 0.0:
        new_k = new_k.conjugate()
    return new_k

def waterman_truell_scattering(lambda_vacuum_nm, N_nm,radius_nm,ri_medium):
    #
    # Solve the waterman truell equation for scattering of an air bubble embedded in a lossy medium
    # N_nm is the number density of bubbles in the volume nm^3
    # k_nm is the wavenumber of the incoming wave (complex)
    # radius_nm is the radius of the bubble
    # ri_medium is the refractive index of the medium the bubble is in
    #
    #
    k_nm = 2*PI*ri_medium/lambda_vacuum_nm
    # The size parameter is now also complex and dimensionless
    size = k_nm*radius_nm
    refractive_index = 1.0 / ri_medium
    # Calculate the forward and backward scattering amplitude
    s10,s20 = Mie.MieS1S2(refractive_index, size*ri_medium, 1)
    s11,s21 = Mie.MieS1S2(refractive_index, size*ri_medium,-1)
    # the normalisation by 1/k_nm is performed when f is calculated
    i = complex(0,1)
    f0 = i*s10
    f1 = i*s11
    # print('Waterman_truell',abs(f0+f1))
    k2 = k_nm*k_nm
    f = 2*PI*N_nm/(k_nm*k_nm*k_nm)
    new_k = np.sqrt( k2 * ( (1+f*f0)*(1+f*f0) - f*f1*f*f1 ) )
    if new_k.imag < 0.0:
        new_k = new_k.conjugate()
    return new_k

def calculate_centre_of_mass(xyzs, masses):
   '''Calculate centre of mass'''
   cm = np.zeros(3)
   mass = 0.0
   for m,xyz in zip(masses,xyzs):
       mass += m
       cm   += m*xyz
   cm /= mass
   return mass,cm

def orthogonalise_projection_operator(ps):
   '''Orthogonalise the projection operators ps'''
   # The projection operator has dimension [6,natoms*3]
   maxcyc = 10
   cycle = 0
   max_overlap = 1.0
   while cycle < maxcyc and max_overlap > 1.0E-8:
       cycle += 1
       max_overlap = 0.0
       for i,p in enumerate(ps):
           # Normalise the projection operator
           anorm = np.dot(p,p)
           shift = 0.0
           if anorm < 1.0e-12:
               shift = 1.0E-12
           p = p / np.sqrt(np.dot(p,p) + shift )
           ps[i] = p
           for j,q in enumerate(ps):
               if j > i:
                   # Gramm Schmidt orthogonoalisation
                   dotprod = np.dot(p,q)
                   ps[j] = q - dotprod*p
                   if max_overlap < dotprod:
                       max_overlap = dotprod
   if cycle >= maxcyc:
       print('WARNING Schmidt Orthogonalisation Failed', max_overlap)
       exit()
   return ps

def construct_projection_operator(atoms, xyzs, masses, nats):
   '''Construct the projection operator for the molecule defined by atoms, xyz, masses'''
   mass,cm = calculate_centre_of_mass(xyzs,masses)
   # The projection operator has dimension number_of_constraints*natoms*3
   ps = np.zeros( (6,nats*3) )
   x = 0
   y = 1
   z = 2
   for i,mass,xyz in zip(atoms,masses,xyzs):
       ps[0,i*3+x] = math.sqrt(mass)
       ps[1,i*3+y] = math.sqrt(mass)
       ps[2,i*3+z] = math.sqrt(mass)
       # coordinates relative to the centre of mass
       relxyz = math.sqrt(mass) * (xyz - cm)
       # First rotations about x in the y/z plane
       # zz is really r*sin(theta), and sin(theta) = zz/r
       # yy is really r*cos(theta), and cos(theta) = yy/r
       ps[3,i*3+x] = 0.0
       ps[3,i*3+y] = -relxyz[z]
       ps[3,i*3+z] = +relxyz[y]
       # Next rotations about y in the x/z plane
       ps[4,i*3+x] = +relxyz[z]
       ps[4,i*3+y] = 0.0
       ps[4,i*3+z] = -relxyz[x]
       # Next rotations about z in the x/y plane
       ps[5,i*3+x] = -relxyz[y]
       ps[5,i*3+y] = +relxyz[x]
       ps[5,i*3+z] = 0.0
   return ps

def calculate_energy_distribution(cell, frequencies, normal_modes, debug=False):
   '''Calculate energy distribution in the phonon modes
      frequencies are the frequencies in cm-1
      normal_modes are the mass weighted normal modes'''
   molecules = cell.molecules
   atomic_masses = cell.atomic_masses
   xyz = cell.xyz_coordinates
   nats = len(xyz)
   # Calculate the projections operators for each molecule
   molecular_projection_operators = []
   molecule_masks = []
   for atoms in molecules:
       mol_xyzs  = [ xyz[atom] for atom in atoms]
       mol_mask = np.zeros(nats*3)
       for atom in atoms:
           mol_mask[3*atom+0] = 1
           mol_mask[3*atom+1] = 1
           mol_mask[3*atom+2] = 1
       mol_masses = [ atomic_masses[atom] for atom in atoms]
       projection_operators = construct_projection_operator(atoms,mol_xyzs,mol_masses,nats)
       projection_operators = orthogonalise_projection_operator(projection_operators)
       molecular_projection_operators.append(projection_operators)
       #
       molecule_masks.append(mol_mask)
   # Calculate the contributions to the kinetic energy in each mode
   energies_in_modes = []
   for imode,mode in enumerate(normal_modes):
       mode_cm = mode
       centre_of_mass_energy = 0.0
       rotational_energy = 0.0
       molecular_energies = []
       for imol,(ps,mask) in enumerate(zip(molecular_projection_operators,molecule_masks)):
           # Calculate total kinetic energy
           total_energy = np.dot(mode,mode)
           # Project out centre of mass motion of each molecule
           mode_cm = mode - np.dot(mode,ps[0])*ps[0] - np.dot(mode,ps[1])*ps[1] - np.dot(mode,ps[2])*ps[2]
           centre_of_mass_energy += total_energy - np.dot(mode_cm,mode_cm)
           # Project out molecular rotation of each molecule
           mode_cm = mode - np.dot(mode,ps[3])*ps[3] - np.dot(mode,ps[4])*ps[4] - np.dot(mode,ps[5])*ps[5]
           rotational_energy += total_energy - np.dot(mode_cm,mode_cm)
           # Now work out the energy of the molecule
           mode_cm = mask * mode
           mol_energy = np.dot(mode_cm,mode_cm)
           molecular_energies.append(mol_energy)
       # end for imol,ps
       vibrational_energy = total_energy - centre_of_mass_energy - rotational_energy
       energies_in_modes.append( (total_energy,centre_of_mass_energy,rotational_energy,vibrational_energy,molecular_energies) )
   # end for mode in normal modes
   return energies_in_modes
# end def


def hodrick_prescott_filter(y,damping,lambda_value,niters):
    #
    # Apply a Hodrick Prescott filter to the spectrum in x, y
    # y is the experimental absorption
    # damping is used to damp the iterations
    # lambda_value is the chosen smoothing factor
    # Based on ideas in the thesis of Mayank Kaushik (University Adelaide)
    #
    from scipy import sparse
    #
    # Create a sparse 3rd order difference operator
    #
    n = len(y)
    diag = np.ones(n-3)
    D = sparse.spdiags( [-1*diag, 3*diag, -3*diag, 1*diag],
                        [0, -1, -2, -3], n, n-3 ).tocsc()
    w = np.ones(n)
    for it in range(10*niters):
        W = sparse.spdiags(w, 0, n, n).tocsc()
        # Problems with overflow if lambda is large
        try:
            Z = W + pow(10,lambda_value) * (D.dot(D.transpose()))
        except:
            pass
        z = sparse.linalg.spsolve(Z, w*y)
        # residuals = y - z
        # error = sum( r*r for r in residuals if r < 0.0 )
        # error = math.sqrt(error/n)
        w = damping*(y>z) + (1-damping)*(y<z)
        # print(it, error)
    return y-z

def reflectance_atr(ns,n0,theta,atrSPolFraction):
    #
    # Calculate the atr s and p reflectance
    # ns is the complex permittivity of the effective medium
    # n0 is the permittivity of atr material
    # theta is the angle of incidence in degrees
    # atrSPolFraction is the fraction of S wave to be considered
    #                 The amount of P wave is 1-atrSPolFraction
    #
    # rs is the s-wave Fresnel amplitude
    # rp is the p-wave

    # Convert theta to an angle in radians
    theta = math.radians(theta)
    costheta = math.cos(theta)
    sintheta = math.sin(theta)
    # Calculate the Fresnel amplitudes for reflection
    # Equations taken from p368 of APPLIED SPECTROSCOPY REVIEWS
    # Vol. 39, No. 3, pp. 365–384, 2004
    # by Milan Milosevic
    # DOI: 10.1081/ASR-200030195
    rs = -1.0* (n0*costheta - cmath.sqrt(ns*ns - n0*n0*sintheta*sintheta)) / (n0*costheta + cmath.sqrt(ns*ns - n0*n0*sintheta*sintheta))
    rp = (ns*ns/n0*costheta - cmath.sqrt(ns*ns - n0*n0*sintheta*sintheta)) / (ns*ns/n0*costheta + cmath.sqrt(ns*ns - n0*n0*sintheta*sintheta))
    # Calculate the reflectance from the amplitudes - store as a real
    RS = np.real(rs * rs.conjugate())
    RP = np.real(rp * rp.conjugate())
    RSP = atrSPolFraction*RS + (1.0-atrSPolFraction)*RP
    # Now return the extinction
    RSP = -math.log10(RSP)
    return RSP

def solve_single_crystal_equations( 
        superstrateDielectricFunction ,
        substrateDielectricFunction   ,
        crystalPermittivityFunction   ,
        superstrateDepth              ,
        substrateDepth                ,
        crystalDepth                  ,
        mode                          ,
        theta                         ,
        phi                           ,
        psi                           ,
        angleOfIncidence              ,
        v ):
    """ This is a parallel call to the single crystal equation solver,
    system is a GTM system"""
    # Create 3 layers, thickness is converted from microns to metres
    superstrate      = GTM.Layer(thickness=superstrateDepth*1e-6,epsilon1=superstrateDielectricFunction)
    substrate        = GTM.Layer(thickness=substrateDepth*1e-6,  epsilon1=substrateDielectricFunction)
    crystal          = GTM.Layer(thickness=crystalDepth*1e-9,    epsilon=crystalPermittivityFunction)
    # Creat the system with the layers 
    if mode == 'Thick slab':
        system = GTM.System(substrate=crystal, superstrate=superstrate, layers=[])
    elif mode == 'Coherent thin film':
        system = GTM.System(substrate=substrate, superstrate=superstrate, layers=[crystal])
    else:
        system = GTM.System(substrate=substrate, superstrate=superstrate, layers=[crystal])
    # Rotate the dielectric constants to the laboratory frame
    system.substrate.set_euler(theta, phi, psi)
    system.superstrate.set_euler(theta, phi, psi)
    for layer in system.layers:
        layer.set_euler(theta, phi, psi)
    # 
    # convert cm-1 to frequency
    #
    freq = v * speed_light_si * 1e2
    system.initialize_sys(freq)
    zeta_sys = np.sin(angleOfIncidence)*np.sqrt(system.superstrate.epsilon[0,0])
    Sys_Gamma = system.calculate_GammaStar(freq, zeta_sys)
    r, R, t, T = system.calculate_r_t(zeta_sys)
    if len(system.layers) > 0:
        epsilon = system.layers[0].epsilon
    else:
        epsilon = system.substrate.epsilon
    return v,r,R,t,T,epsilon


def cleanup_symbol(s):
    """Return a true element from the symbol"""
    s = s.capitalize()
    s = s.replace('_','')
    for i in string.digits:
        s = s.replace(i,'')
    return s

def euler_rotation(vector, theta, phi, psi):
     """Apply a passive Euler rotation to the vector"""
     euler = np.zeros( (3,3) )
     euler[0, 0] =  np.cos(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.sin(psi)
     euler[0, 1] = -np.sin(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.cos(psi)
     euler[0, 2] =  np.sin(theta) * np.sin(phi)
     euler[1, 0] =  np.cos(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.sin(psi)
     euler[1, 1] = -np.sin(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.cos(psi)
     euler[1, 2] = -np.sin(theta) * np.cos(phi)
     euler[2, 0] =  np.sin(theta) * np.sin(psi)
     euler[2, 1] =  np.sin(theta) * np.cos(psi)
     euler[2, 2] =  np.cos(theta)
     result = np.matmul(euler, vector)
     return result

def get_pool(ncpus, threading, initializer=None, initargs=None, debugger=None ):
     """Return a pool of processors given the number of cpus and whether threading is requested"""
     if debugger is not None:
         debugger.print('get_pool ncpus = ',ncpus)
         debugger.print('get_pool threading = ',threading)
         debugger.print('get_pool initializer = ',initializer)
     # Switch off mkl threading
     try:
         import mkl
         mkl.set_num_threads(1)
     except:
         pass
     # see if threading has been requested
     if threading:
         from multiprocessing.dummy import Pool
         pool = Pool(ncpus, initializer=initializer, initargs=initargs)
     else:
         from multiprocessing import Pool
         pool = Pool(ncpus, initializer=initializer, initargs=initargs )
     return pool

def set_affinity_on_worker():
    '''When a new worker process is created, the affinity is set to all CPUs'''
    #JK print('I'm the process %d, setting affinity to all CPUs.' % os.getpid())
    #JK Commented out for the time being
    #JK os.system('taskset -p 0xff %d > /dev/null' % os.getpid())
    return


