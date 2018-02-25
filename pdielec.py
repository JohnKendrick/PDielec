#!/usr/bin/env python
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
"""PDielec driver program to calculate dielectric response at infrared and THz frequencies"""
from __future__ import print_function
import math
import os
import sys
import numpy as np
import ctypes
from multiprocessing import Pool, cpu_count, Array
from Python.Constants import amu, PI, avogadro_si, wavenumber, angstrom, isotope_masses, average_masses
from Python.Constants import support_matrix_db
from Python.VaspOutputReader import VaspOutputReader
from Python.PhonopyOutputReader import PhonopyOutputReader
from Python.CastepOutputReader import CastepOutputReader
from Python.GulpOutputReader import GulpOutputReader
from Python.CrystalOutputReader import CrystalOutputReader
from Python.AbinitOutputReader import AbinitOutputReader
from Python.QEOutputReader import QEOutputReader
from Python.ExperimentOutputReader import ExperimentOutputReader
from Python.Plotter import Plotter, print3x3, print_reals
import Python.Calculator as Calculator
from Python.Utilities import get_reader

def set_affinity_on_worker():
    """When a new worker process is created, the affinity is set to all CPUs"""
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK Commented out for the time being
    #JK os.system("taskset -p 0xff %d > /dev/null" % os.getpid())


def main():
    """Main Driver routine for PDielec - reads the input and directs the calculation"""
    def show_usage():
        """Show pdielec usage"""
        print("pdielec [-method ap/balan/maxwell/bruggeman/mie] [-plot molar_absorption]")
        print("        [-sphere] [-needle h k l] [-plate h k l] [-ellipsoid h k l aoverb]")
        print("        [-print] [-sigma {:.1f}] [-mode_sigma mode {:.1f}] ".format(sigma, sigma))
        print("        [-optical xx yy zz] [-optical_tensor xx yy zz xy xz yz] [-vmin {:.1f}] [-vmax {:.1f}]".format(vmin, vmax))
        print("        [-vf {:.3f}] ".format(volume_fraction))
        print("        [-mf {:.3f}] ".format(mass_fraction))
        print("        [-masses program/average/isotope] ")
        print("        [-mass element mass] ")
        print("        [-dielectric {:.4f}] [-density {:.4f}] ".format(matrix_dielectric, matrix_density))
        print("        [-matrix nujol/ptfe/kbr/ldpe/mdpe/hdpe] ")
        print("        [-LO h k l] [-LO_cart x y z]")
        print("        [-mode index] [-csv file] [-csv_ext file")
        print("        [-eckart] [-neutral]")
        print("        [-hessian crystal/symm]")
        print("        [-h -help]")
        print("        [-program name] file1 [file2]")
        print("pdielec: Calculates an IR spectrum from VASP, CASTEP, PWscf, Abinit, Phonopy or GULP")
        print("             At least one shape must be specified (-sphere -needle -plate)")
        print("         -program castep/crystal/vasp/abinit/qe/phonopy/gulp/experiment")
        print("             Specifies the program that created the files to be processed")
        print("             some programs require two or more files to be specified")
        print("             In the case of VASP file1 must be OUTCAR")
        print("             In the case of CASTEP file1 must be the seed name of the calculation")
        print("             In the case of GULP file1 must be the output file name of the calculation")
        print("             In the case of Crystal file1 must be the output file name of the calculation")
        print("             In the case of Quantum Espresso file1 must be the output file name of the calculation")
        print("             In the case of Abinit file1 must be the output file name of the abinit calculation")
        print("             In the case of phonopy file1 must the qpoints.yaml file")
        print("             In the case of experiment the file should be a file containing experimental information")
        print("         -method ap/balan/maxwell/bruggeman/bruggeman_iter/bruggeman_minimise/maxwell_sihvola/mie")
        print("             Choose the method to be used to, mie, balan, maxwell or bruggeman")
        print("             maxwell_sihvola should give the same results as maxwell but could be used to deal with chirality")
        print("             bruggeman and brugemann_iter solve the bruggeman equations iteratively..")
        print("             brugemann_minimise solves the bruggeman equations by minimisation..")
        print("             Note that mie only applies to spherical particles..")
        print("             The option can be used more than once")
        print("         -plot real/imaginary/absorption/molar_absorption/extinction/molar_extinction")
        print("             The real, imaginary or the absorption coefficent of ")
        print("             the dielectric are plotted the default is absorption/extinction coefficient")
        print("             The units of absorption and extinction coefficients are m-1")
        print("             The option can be used more than once")
        print("         -vf 0.1")
        print("             Include the following number as a volume fraction")
        print("             The option can be used more than once")
        print("         -mf 0.1")
        print("             Include the following number as a mass fraction")
        print("             The option can be used more than once")
        print("         -masses program")
        print("             By default the atomic masses are defined by the average weight")
        print("             The options are program/average/isotope           ")
        print("             average means the weight is averaged using the natural abundance")
        print("             isotope means the most abundant isotopic mass is used")
        print("             program means that the QM/MM programs masses are used")
        print("         -mass element symbol")
        print("             Allows individual elements to have their masses defined")
        print("             the -masses mass definition is applied first,     ")
        print("             then any masses defined by -mass are applied,     ")
        print("         -size a [width]")
        print("             specifies the radius of spherical particle to be considered in microns")
        print("             only the mie, maxwell and bruggeman methods incorporation size effects")
        print("             in the case of the mie method a log normal size distribution can be specified")
        print("             by giving the mean of the log distribution (in microns) and its width (log microns)")
        print("         -print")
        print("           print additional information about the calculation")
        print("         -csv file.csv")
        print("           print a csv file")
        print("         -csv_ext file.csv")
        print("           extended print option creating 3 csv files")
        print("           file_command.csv, file_frequencies.csv file_spectrum.csv")
        print("         -excel  file.xlsx")
        print("           output the results to an excel file")
        print("         -sigma gives the Lorentzian sigma width in wavenumbers default is %f" % sigma)
        print("         -mode_sigma mode sigma")
        print("             A given mode is assigned its own width")
        print("             all other modes not specified use the width specified by -sigma")
        print("             The directive can be used many times")
        print("             Warning a mode_sigma command is needed for all degenerate components")
        print("         -matrix ptfe")
        print("             Defines material, density and permittivity of the supporting medium")
        print("             Possible materials are;")
        for f in support_matrix_db:
            (rho, eps) = support_matrix_db[f]
            print("                 " + "{:6s} density ".format(f) + "{:.4f}".format(rho) + " permittivity " + " {:.4f}".format(eps))
        print("         -dielectric is the dielectric constant of the medium")
        print("             For kbr the value is 2.25")
        print("             For ptfe the value is 2.0")
        print("             For nujol the value is 2.155 (refractive index is 1.468) ")
        print("             For vacuum or air the value is 1.0")
        print("         -vmax vmax  report all frequencies from vmin to vmax")
        print("         -vmin vmin  report all frequencies from vmin to vmax")
        print("         -sphere")
        print("             Adds sphere shape to the shapes to be processed")
        print("             The option can be used more than once")
        print("         -needle h k l")
        print("             Adds the needle shape to the shapes to be processed")
        print("             hkl (integers) define the unique needle direction as [hkl]")
        print("             The option can be used more than once")
        print("         -plate h k l")
        print("             Adds the plate shape to the shapes to be processed")
        print("             hkl (integers) are miller indices, the plate lies in the (hkl) plane")
        print("             for non-orthogonal cells [hkl] and the normal to (hkl) are not necessarily the same")
        print("             The option can be used more than once")
        print("         -ellipsoid h k l aoverb")
        print("             Adds the ellipsoid shape to the list of shapes to be processed")
        print("             the unique ellipsoid direction is [hkl] ")
        print("             aoverb (a/b) is the eccentricity of the ellipsoid")
        print("             oblate ellipsoid eccentricity (a/b) < 1")
        print("             prolate ellipsoid eccentricity (a/b) > 1")
        print("             The option can be used more than once")
        print("         -optical xx yy zz are the optical permitivitties (diagonal only) ")
        print("             If the tensor is not diagonal use the -optical_tensor flag")
        print("         -optical_tensor xx yy zz xy xz yz is the full optical permitivitty tensor")
        print("             If either optical commands are given on the command line they overide")
        print("             the optical permittivity")
        print("         -i  step   gives the increment of frequency (0.2 cm-1)")
        print("         -ignore mode")
        print("             Ignore a mode in the construction of the dielectric")
        print("             by default all modes with a frequency less than 5cm-1 are ignored")
        print("             -ignore (can be used more than once) is used only modes specified are ignored")
        print("         -mode  index  only include this mode in the sum")
        print("             This option can be included several times")
        print("             By default all frequencies are included")
        print("         -eckart")
        print("             The translational modes will be projected from the dynamical matrix")
        print("             This option only applies when the dynmical matrix is read and used")
        print("             to calculate the frequencies and normal modes.")
        print("         -neutral")
        print("             Charge neutrality of the Born charges is enforced")
        print("         -hessian [crystal|symm]")
        print("             By default the hessian is symmetrised by forming 0.5 * ( Ht + H )")
        print("             however the crystal package symmetrises by making the upper triangle the same as the lower")
        print("             For compatibility -hessian crystal will symmetrise in the same way as the Crystal package")
        print("         -LO h k l ")
        print("             Define a wavector direction (hkl) for which the LO frequencies ")
        print("             will be calculated using the nonanalytical correction")
        print("             This option can be included several times")
        print("         -LO_cart x y z ")
        print("             As above but the directions are defined in cartesian coordinates")
        print("         -thresholds 1.0E-10 5")
        print("             Set the thresholds for treating a phonon mode")
        print("             Only modes with an intensity > 1.0E-10 and a frequency above 5cm-1")
        print("             Will be conidered in the calcualtion of the absorption spectrum")
        print("         -molesof [atoms|cells|molecules natoms]")
        print("             Definition of moles for the calculation of concentration.  The default is 'cells'")
        print("             if 'molecules' is given then the number of atoms in the molecules is needed")
        print("         -processors 0")
        print("             Use the specificed number of processors, the default is to use all available processors")
        print("         -drude drude_input_plasma drude_input_sigma")
        print("             Tells the program to include a Drude for metals in the calculation")
        print("             of the dielectric.")
        print("             drude_input_plasma is the plasma frequency and ")
        print("             drude_input_sigma is its sigma parameter")
        print("         -help to print out this help information")
        print("         -h    to print out this help information")
        return

    # Create a plotter to plot/print the final data
    plotter = Plotter()
    # define some constants which may change due to the parameters on the command line
    debug          = False
    increment      = 0.2
    program        = ""
    qmprogram      = ""
    optical        = []
    optical_tensor = []
    my_modes       = []
    qlist          = []
    qlist_input    = []
    qdata          = []
    ignore_modes   = []
    mode_sigmas    = {}
    mass_dictionary= {}
    vmax           = 300.0
    vmin           = 0.0
    shapes         = []
    shape_data     = []
    names          = []
    methods        = []
    volume_fractions = []
    mass_fractions = []
    mass_fractions_string = []
    mass_definition = 'average'
    fractional_types = []
    plot_types      = []
    sizes = []
    size_distribution_sigmas = []
    volume_fraction = 0.1
    mass_fraction = 0.0
    sigma          = 5.0
    eckart = False
    neutral = False
    matrix_dielectric  = 1.00  # air
    matrix_dielectric  = 2.25  # kbr
    matrix_dielectric  = 2.0   # ptfe
    matrix_density = 2.2
    print_info      = False
    csvfile        = ""
    excelfile      = ""
    csv_extended   = False
    drude = False
    drude_input_plasma = 0.0
    drude_input_sigma = 0.0
    threshold = {'intensity': 1.0E-10, 'frequency': 5}
    hessian_symmetrisation = "symm"
    number_of_processors = 0
    moles_of = 'cells'

    # check usage
    if len(sys.argv) <= 1:
        show_usage()
        exit()

    # Begin processing of command line
    command_line = ' '.join(sys.argv)
    tokens = sys.argv[1:]
    ntokens = len(tokens)
    itoken = 0
    while itoken < ntokens:
        token = tokens[itoken]
        if token == "-program":
            itoken += 1
            program = tokens[itoken]
            if program == "phonopy":
                itoken += 1
                qmprogram = tokens[itoken]
        elif token == "-sigma":
            itoken += 1
            sigma = float(tokens[itoken])
        elif token == "-print":
            print_info = True
        elif token == "-mode_sigma":
            itoken += 1
            mode = int(tokens[itoken])
            itoken += 1
            mode_sigmas[mode] = int(tokens[itoken])
        elif token == "-method":
            itoken += 1
            methods.append(tokens[itoken])
        elif token == "-plot":
            itoken += 1
            plot_types.append(tokens[itoken])
        elif token == "-vf":
            itoken += 1
            volume_fractions.append(float(tokens[itoken]))
            fractional_types.append("vf=%s" % (tokens[itoken]))
        elif token == "-mass":
            itoken += 1
            element = tokens[itoken]
            itoken += 1
            mass = float(tokens[itoken])
            mass_dictionary[element] = mass
        elif token == "-masses":
            itoken += 1
            if tokens[itoken][0] == "p":  mass_definition = "program"
            if tokens[itoken][0] == "a":  mass_definition = "average"
            if tokens[itoken][0] == "i":  mass_definition = "isotope"
        elif token == "-mf":
            itoken += 1
            mass_fractions.append(float(tokens[itoken]))
            mass_fractions_string.append(tokens[itoken])
        elif token == "-size":
            itoken += 1
            sizes.append(float(tokens[itoken]))
            if tokens[itoken+1][0] != '-':
                itoken += 1
                size_distribution_sigmas.append(float(tokens[itoken]))
            else:
                size_distribution_sigmas.append(None)
        elif token == "-dielectric":
            itoken += 1
            matrix_dielectric = float(tokens[itoken])
        elif token == "-density":
            itoken += 1
            matrix_density = float(tokens[itoken])
        elif token == "-matrix":
            itoken += 1
            token = tokens[itoken]
            if token in support_matrix_db:
                (matrix_density, matrix_dielectric) = support_matrix_db[token]
            else:
                print("Error, unkown material {:s}", token)
                exit(1)
        elif token == "-optical":
            optical = [complex(f) for f in tokens[itoken+1:itoken+4]]
            itoken = itoken + 3
            # Check that the optical diagonal is real or complex
            odc = np.array(optical)
            odi = np.absolute(np.imag(odc))
            sumi = np.sum(odi)
            if sumi < 1.0e-12:
                odr = np.real(odc)
                optical = odr.tolist()
        elif token == "-optical_tensor":
            temp = [complex(f) for f in tokens[itoken+1:itoken+10]]
            itoken = itoken + 9
            optical_tensor.append(temp[0:3])
            optical_tensor.append(temp[3:6])
            optical_tensor.append(temp[6:9])
            # Check that the optical is real or complex
            odc = np.array(optical_tensor)
            odi = np.absolute(np.imag(odc))
            sumi = np.sum(odi)
            if sumi < 1.0e-12:
                odr = np.real(odc)
                optical_tensor = odr.tolist()
        elif token == "-i":
            itoken += 1
            increment = float(tokens[itoken])
        elif token == "-mode":
            itoken += 1
            my_modes.append(int(tokens[itoken]))
        elif token == "-ignore":
            itoken += 1
            ignore_modes.append(int(tokens[itoken]))
        elif token == "-vmin":
            itoken += 1
            vmin = float(tokens[itoken])
        elif token == "-vmax":
            itoken += 1
            vmax = float(tokens[itoken])
        elif token == "-csv_ext":
            itoken += 1
            csvfile = tokens[itoken]
            csv_extended = True
        elif token == "-csv":
            itoken += 1
            csvfile = tokens[itoken]
            csv_extended = False
        elif token == "-excel":
            itoken += 1
            excelfile = tokens[itoken]
        elif token == "-sphere":
            shapes.append("sphere")
            shape_data.append(" ")
        elif token == "-plate":
            shapes.append("plate")
            itoken += 1
            token = tokens[itoken]
            if token[0] == "[" or token[0] == "(" or token[0] == "{":
                data = token
            else:
                data = "(" + token + "," + tokens[itoken+1] + "," + tokens[itoken+2] + ")"
                itoken += 2
            shape_data.append(data)
        elif token == "-ellipsoid":
            shapes.append("ellipsoid")
            itoken += 1
            token = tokens[itoken]
            if token[0] == "[" or token[0] == "(" or token[0] == "{":
                data = token
            else:
                data = "[" + token + "," + tokens[itoken+1] + "," + tokens[itoken+2] + "]"
                itoken += 2
            itoken += 1
            aoverb = float(tokens[itoken])
            shape_data.append((data, aoverb))
        elif token == "-needle":
            shapes.append("needle")
            itoken += 1
            token = tokens[itoken]
            if token[0] == "[" or token[0] == "(" or token[0] == "{":
                data = token
            else:
                data = "[" + token + "," + tokens[itoken+1] + "," + tokens[itoken+2] + "]"
                itoken += 2
            shape_data.append(data)
        elif token == "-eckart":
            eckart = True
        elif token == "-hessian":
            itoken += 1
            token = tokens[itoken]
            if token == "crystal":
                hessian_symmetrisation = token
            elif token == "symm":
                hessian_symmetrisation = token
            else:
                print("The -hessian directive must be qualified with \"symm\" or \"crystal\" {:s}".format(token))
                exit(1)
        elif token == "-neutral":
            neutral = True
        elif token == "-processors":
            itoken += 1
            number_of_processors = int(tokens[itoken])
        elif token == "-molesof":
            itoken += 1
            moles_of = tokens[itoken]
            if moles_of == 'molecules':
                itoken += 1
                number_of_atoms_per_molecule = int(tokens[itoken])
            elif moles_of == 'atoms':
                pass
            elif moles_of == 'cells':
                pass
            else:
                print("The -molesof directive has a wrong argument \"molecules\",\"atoms\" or \"cells\" {:s}".format(moles_of))
                exit(1)
        elif token == "-LO_cart":
            itoken += 1
            qlist_input.append([float(q) for q in tokens[itoken:itoken+3]])
            itoken += 2
        elif token == "-LO":
            itoken += 1
            token = tokens[itoken]
            if token[0] == "[" or token[0] == "(" or token[0] == "{":
                data = token
            else:
                data = "(" + token + "," + tokens[itoken+1] + "," + tokens[itoken+2] + ")"
                itoken += 2
            qdata.append(data)
        elif token == "-h":
            show_usage()
            exit()
        elif token == "-help":
            show_usage()
            exit()
        elif token == "-drude":
            drude = True
            itoken += 1
            drude_input_plasma = float(tokens[itoken])
            itoken += 1
            drude_input_sigma  = float(tokens[itoken])
        elif token == "-debug":
            debug = True
        elif token == "-thresholds":
            itoken += 1
            threshold['intensity'] = float(tokens[itoken])
            itoken += 1
            threshold['frequency'] = float(tokens[itoken])
        elif token[0] == "-":
            print("Error on input unkown option: {:s}".format(token))
            exit(1)
        else:
            names.append(token)
        itoken += 1
        # end loop over tokens

    # Look for obvious errors
    programs = ["castep", "abinit", "qe", "castep", "vasp", "crystal", "gulp", "phonopy", "experiment"]
    if program is not "":
        if program not in programs:
            print("program specified is: {:s}".format(program))
            print("Needs to be one of: " + " ".join(p for p in programs))
            exit(1)
    if qmprogram is not "":
        if qmprogram not in programs:
            print("QM program specified is: {:s}".format(program))
            print("Needs to be one of: " + " ".join(p for p in programs))
            exit(1)
    if len(names) <= 0:
        print("No files were specified")
        show_usage()
        exit(1)
    if len(shapes) <= 0:
        shapes.append("sphere")
        shape_data.append(" ")
        print("No shapes were specified, using a sphere")

    # Set default method
    if len(methods) == 0:
        methods.append('maxwell')

    # print out information
    print("The following methods will be used: " + " ".join(m for m in methods))

    # Process size information, if no sizes have been given the use 0.0
    if len(sizes) > 0:
        print("Particle sizes have been specified (in microns)")
        for a,sig in zip(sizes,size_distribution_sigmas):
            if sig:
                print("Mean radius of spherical particles and log normal sigma: ", a, sig)
            else:
                print("Radius of spherical particles: ", a)
    else:
        sizes.append(0.0)
        size_distribution_sigmas.append(None)
    

    fd_excelfile = None
    fd_csvfile = None
    fd_csvfiles = [None, None, None]
    if csvfile != "":
        if csv_extended:
            head = csvfile.replace(".csv", "")
            print("Creating three csv files: \n   {}, \n   {} \n   {}".format(head+"_command.csv", head+"_frequencies.csv", head+"_spectrum.csv"))
            fd_csvfiles[0] = open(head+"_command.csv", 'w')
            fd_csvfiles[1] = open(head+"_frequencies.csv", 'w')
            fd_csvfiles[2] = open(head+"_spectrum.csv", 'w')
        else:
            print("Creating a csv file: "+csvfile)
            fd_csvfiles[0] = open(csvfile, 'w')
            fd_csvfiles[1] = fd_csvfiles[0]
            fd_csvfiles[2] = fd_csvfiles[0]
        fd_csvfile = fd_csvfiles[0]
    if excelfile != "":
        print("Creating an excel file: "+excelfile)
        fd_excelfile = 1
    if len(plot_types) == 0:
        print("No plotting has been requested")
    else:
        print("Plotting types requested are: " + " ".join(p for p in plot_types))
    print("The frequency increment is {:.2f} cm-1".format(increment))
    print("The default width factor (-sigma) has been set to {:3f} cm-1".format(sigma))
    if len(mode_sigmas) > 0:
        for mode in mode_sigmas:
            print("Mode {:3d} has a width factor of: {:5.2f}".format(mode, mode_sigmas[mode]))
    print("The permittivity of the medium (-dielectric) has been set to {:.4f}".format(matrix_dielectric))
    print("The density of the medium (-density) has been set to {:.4f} g/cc".format(matrix_density))
    if len(optical_tensor) > 0:
        print3x3("The optical permitivity tensor has been specified: ", optical_tensor)
    if len(optical) > 0:
        print(" ")
        print("Optical permitivity diagonal elements specified: " + " ".join("{:.4f}".format(p) for p in optical))
        optical_tensor = Calculator.initialise_complex_diagonal_tensor(optical)
        # Check that the optical is real or complex
        odc = np.array(optical_tensor)
        odi = np.absolute(np.imag(odc))
        sumi = np.sum(odi)
        if sumi < 1.0e-12:
            odr = np.real(odc)
            optical_tensor = odr.tolist()
    # end if len(optical)
    print("Vmin is: {:7.1f} cm-1".format(vmin))
    print("Vmax is: {:7.1f} cm-1".format(vmax))
    if len(my_modes) > 0:
        print("Only consider contributions from phonon modes:" + " ".join(str(m) for m in my_modes))
    else:
        if len(ignore_modes) > 0:
            print("Ignoring the following phonon modes:")
            print("     " + " ".join("{:d}".format(m) for m in ignore_modes))
        else:
            print("All phonon modes will be considered")
    drude_sigma = 0.0
    drude_plasma = 0.0
    if drude:
        print("A Drude model will be used as well as the phonon modes")
        print("The Drude sigma parameter  is; {:6.1f} cm-1".format(drude_input_sigma))
        print("The Drude plasma frequency is; {:6.1f} cm-1".format(drude_input_plasma))
        drude_sigma = drude_input_sigma * wavenumber
        drude_plasma = drude_input_plasma * wavenumber
    # end if
    # Parallel options handled here
    if number_of_processors == 0:
        number_of_processors = cpu_count()
        print("The number of processors ({}) has been determined automatically".format(number_of_processors))
    else:
        print("The number of processors has been set to {} ".format(number_of_processors))
    #
    # Arrays will be used where possible from now on, also some arrays will contain complex numbers
    #

    # 
    # Get the reader
    #
    reader = get_reader(program,names,qmprogram)

    if eckart:
        print("The translational modes will be projected from the Dynamical matrix where possible")
    else:
        print("No projection of the dynamical matrix will be performed")
    if neutral:
        print("The charge neutrality of the Born atomic charges will be enforced")
    else:
        print("The charge neutrality of the Born atomic charges will not be enforced")
    # Modify the default settings of the reader
    reader.debug = debug
    # The order these settings is applied is important - the hessian symmetry is determined by the program so apply
    # The non mass weighted hesian stored in reader will be symmetrised
    reader.hessian_symmetrisation = hessian_symmetrisation
    # Initiate reading of the files - at this point any masses are set by the program
    reader.read_output()
    # Now apply the eckart and neutral conditions
    reader.eckart = eckart
    if neutral:
        reader.neutralise_born_charges()
    if debug or print_info:
        reader.print_info()
    # Calculate the depolarisation matrices.  Put here because we need the lattice vectors
    depolarisations = []
    for shape, data in zip(shapes, shape_data):
        label = ""
        printd = data
        if shape == "ellipsoid":
            (printd, label) = data
        title = "Shape depolarisation matrix: {} {} {}".format(shape, printd, label)
        if shape == "sphere":
            direction = []
            depolarisations.append(Calculator.initialise_sphere_depolarisation_matrix())
        if shape == "plate":
            direction = Calculator.direction_from_shape(data,reader)
            depolarisations.append(Calculator.initialise_plate_depolarisation_matrix(direction))
        if shape == "needle":
            direction = Calculator.direction_from_shape(data,reader)
            depolarisations.append(Calculator.initialise_needle_depolarisation_matrix(direction))
        if shape == "ellipsoid":
            ellipsez, aoverb = data
            direction = Calculator.direction_from_shape(ellipsez,reader)
            depolarisations.append(Calculator.initialise_ellipsoid_depolarisation_matrix(direction, aoverb))
        # Print out the depolarisation matrix
        if len(direction) > 0:
            title = title + "xyz=[" + ",".join("{:f}".format(d) for d in direction) + "]"
        print3x3(title, depolarisations[-1])
    if len(optical) > 0:
        print()
    # What mass definition are we using?
    if not mass_definition == "program" or mass_dictionary:
        print(" ")
        print("Atomic masses will be defined using: ", mass_definition)
        if mass_dictionary:
            print("Additional elemental masses have been specifically specified:")
            for element in mass_dictionary:
                print("    ",element,": ",mass_dictionary[element])
        print_reals("Old masses by atom types", reader.masses_per_type,format="{:10.6f}")
        print_reals("Old atomic mass list", reader.masses,format="{:10.6f}")
        if mass_definition == "average":
            reader.change_masses(average_masses, mass_dictionary)
        elif mass_definition == "isotope":
            reader.change_masses(isotope_masses, mass_dictionary)
        else:
            print("Failed to change mass definition", mass_definition)
        print_reals("New masses by atom types", reader.masses_per_type,format="{:10.6f}")
        print_reals("New atomic mass list", reader.masses,format="{:10.6f}")
    else:
        print(" ")
        print("Atomic masses will be defined using: ", mass_definition)
        print_reals("Masses by atom types", reader.masses_per_type,format="{:10.6f}")
        print_reals("Atomic mass list", reader.masses,format="{:10.6f}")
    # Define some of the constants
    dielectric_medium = Calculator.initialise_diagonal_tensor([matrix_dielectric, matrix_dielectric, matrix_dielectric])
    # access the information as numpy arrays.  Use atomic units for these arrays
    masses = np.array(reader.masses) * amu
    print("")
    print("Volume of unit cell is {:14.8f} Angstrom^3".format(reader.volume))
    mtotal = 0.0
    for mass in reader.masses:
        mtotal = mtotal + mass
    print("Total unit cell mass is: {:14.8f} g/mol".format(mtotal))
    crystal_density = mtotal/(avogadro_si * reader.volume * 1.0e-24)
    print("Crystal density is: {:9.5f} g/cc".format(crystal_density))
    rho1 = crystal_density
    rho2 = matrix_density
    print("")
    if fd_excelfile is not None:
        import xlsxwriter as xlsx
        workbook = xlsx.Workbook(excelfile)
        worksheet = workbook.add_worksheet('info')
        excel_row = 0;  worksheet.write(excel_row,0,command_line)
        excel_row += 2; worksheet.write(excel_row,0,'Volume of unit cel (Angstrom^3)')
        excel_row += 1; worksheet.write(excel_row,0,reader.volume)
        excel_row += 2; worksheet.write(excel_row,0,'Unit cell mass (g/mol)')
        excel_row += 1; worksheet.write(excel_row,0,mtotal)
        excel_row += 2; worksheet.write(excel_row,0,'Crystal density (g/cc)')
        excel_row += 1; worksheet.write(excel_row,0,crystal_density)
        excel_row += 2; worksheet.write(excel_row,0,'Matrix density (g/cc)')
        excel_row += 1; worksheet.write(excel_row,0,matrix_density)
    # Convert any mass fractions to volume fractions
    if len(mass_fractions) > 0:
        print("The density of the supporting matrix is {:9.4f} g/cc".format(matrix_density))
        print("Converting the following mass fractions to volume fractions:")
        for mf1, string in zip(mass_fractions, mass_fractions_string):
            mf2 = 1.0 - mf1
            vf1 = 1.0 / (1 + (mf2/mf1) * (rho1/rho2))
            volume_fractions.append(vf1)
            fractional_types.append("mf={}".format(string))
            print("   The equivalent volume fraction for mass fraction {:7.4f} is {:9.6f}".format(mf1, vf1))

    # Set default volume fraction if we need to
    if len(volume_fractions) == 0:
        volume_fractions.append(volume_fraction)
        fractional_types.append("vf={:f}".format(volume_fraction).rstrip("0"))

    # the concentration of unit cells per 1000 cc can be calculated in moles / L
    concentration = 1000.0 / (avogadro_si * reader.volume * 1.0e-24)
    print("")
    if moles_of == 'atoms':
        concentration = 1000.0 / (avogadro_si * reader.volume * 1.0e-24 / reader.nions )
        print("Concentration: {:9.4f} moles of atoms/Litre".format(concentration))
    elif moles_of == 'molecules':
        concentration = 1000.0 / (avogadro_si * reader.volume * 1.0e-24 * number_of_atoms_per_molecule / reader.nions )
        print("Concentration: {:9.4f} moles of molecules/Litre (1 molecule = {} atoms)".format(concentration,number_of_atoms_per_molecule))
    else:
        print("Concentration for unit volume fraction {:9.4f} moles/Litre".format(concentration))

    # Get the zero frequency optical tensor from the output
    epsilon_inf = np.array(reader.zerof_optical_dielectric)
    # Initialise the zero frequency optical tensor if it was supplied on the command
    if len(optical_tensor) > 0:
        epsilon_inf = np.array(optical_tensor)
    if np.max(epsilon_inf) < 0.001:
        print("WARNING! no epsilon infinity has been provided")
        print("WARNING! Please supply with -optical or -optical_tensor")
        exit()
    # Get the born charges
    born_charges = np.array(reader.born_charges)
    #
    # Ask the reader to calculate the massweighted normal modes
    # This call allows the projection to be performed if necessary
    #
    mass_weighted_normal_modes = reader.calculate_mass_weighted_normal_modes()
    frequencies_cm1 = np.array(reader.frequencies)
    frequencies = frequencies_cm1 * wavenumber
    volume = reader.volume*angstrom*angstrom*angstrom
    mode_list = []
    sigmas = []
    # set the widths and define the list of modes to me considered
    for imode, frequency in enumerate(reader.frequencies):
        mode_list.append(imode)
        sigmas.append(sigma)
    for mode in mode_sigmas:
        sigmas[mode] = mode_sigmas[mode]

    if debug:
        print("Complete mode list is: " + ",".join("{:4d}".format(m) for m in mode_list))
        print("Sigmas are:            " + ",".join("{:4f}".format(s) for s in sigmas))

    # convert sigmas to wavenumbers
    sigmas_cm1 = np.array(sigmas)
    sigmas = np.array(sigmas)*wavenumber

    if program == 'experiment':
        oscillator_strengths = reader.oscillator_strengths
    else:
        #
        # calculate normal modes in xyz coordinate space
        # should they be re-normalised or not?  According to Balan the mass weighted coordinates should be normalised
        normal_modes = Calculator.normal_modes(masses, mass_weighted_normal_modes)
        # from the normal modes and the born charges calculate the oscillator strengths of each mode
        oscillator_strengths = Calculator.oscillator_strengths(normal_modes, born_charges)
    # calculate the intensities from the trace of the oscillator strengths
    intensities = Calculator.infrared_intensities(oscillator_strengths)
    if fd_excelfile is not None:
        excel_row += 2; worksheet.write(excel_row,0,'Transverse optical frequencies (cm-1)')
        excel_row += 1; [ worksheet.write(excel_row,col,f) for col,f in enumerate(np.real(frequencies_cm1)) ]
    #
    # If LO frequencies have been requested lets calculate them now
    # First calculate the cartesian directions from any wavevectors
    if qdata:
        for q in qdata:
            qlist.append(Calculator.direction_from_shape(q,reader))
        # end for
        print(" ")
    # No add any cartesian directions which were specified
    if qlist_input:
        qlist.extend(qlist_input)
    # Loop over all the cartesian directions to find the splitting
    if qlist:
        lo_freqs = Calculator.longitudinal_modes(frequencies, mass_weighted_normal_modes, born_charges, masses, epsilon_inf, volume, qlist, reader)
        print_reals("Transverse frequencies (cm-1)", np.real(frequencies_cm1))
        print(" ")
        index = 0
        for lo_freq, qxyz in zip(lo_freqs, qlist):
            lo_freq = lo_freq / wavenumber
            qstring = " "
            if index < len(qdata):
                qstring = " " + qdata[index]
            index += 1
            qxyz_string = "[" + ",".join("{:7.4f}".format(xyz) for xyz in qxyz) + "]"
            print_reals("Longitudinal frequencies (cm-1) q->0 q(cart)=" + qxyz_string + qstring, np.real(lo_freq))
            print(" ")
            if fd_excelfile is not None:
                excel_row += 2; worksheet.write(excel_row,0,'Logitudinal optical frequencies (cm-1) q->0 q(cart)='+qxyz_string+qstring)
                excel_row += 1; [ worksheet.write(excel_row,col,f) for col,f in enumerate(np.real(lo_freq)) ]
        # end for qxyz, lo_freq
    # end if qlist
    #
    # Calculate degenerate lists
    #
    degeneracy_threshold = 1.0E-8
    degenerate_lists = {}
    mmax = len(frequencies)
    for m1 in range(len(frequencies)):
        degenerate_lists[m1] = []
    for m1 in range(len(frequencies)):
        f1 = frequencies[m1]
        for m2 in range(m1+1,min(mmax,m1+3)):
            f2 = frequencies[m2]
            if abs(f2-f1) < degeneracy_threshold:
                degenerate_lists[m1].append(m2)
                degenerate_lists[m2].append(m1)
    #
    # Only modes with non-zero oscillator strengths contribute to the dielectric
    # so calculate those modes which we can safely ignore and store them in ignore_modes
    #
    if len(ignore_modes) == 0:
        for mode, intensity in enumerate(intensities):
            # ignore modes with a low oscillator strength
            if intensity < threshold['intensity']:
                # If any of its degenerate modes have intensity then we shouldn't ignore it
                ignore = True
                for m in degenerate_lists[mode]:
                    if intensities[m] > threshold['intensity']:
                        ignore = False
                if ignore:
                    ignore_modes.append(mode)
            # ignore modes with low real frequency
            elif np.real(frequencies[mode])/wavenumber < threshold['frequency']:
                ignore_modes.append(mode)
            # ignore modes with imaginary frequency
            elif abs(np.imag(frequencies[mode]))/wavenumber > 1.0e-6:
                ignore_modes.append(mode)
            # end if intensity
        # end for
    # end if len()
    # Remove any unwanted modes
    ignore_modes = list(set(ignore_modes))
    if len(ignore_modes) > 0:
        for mode in ignore_modes:
            if mode in mode_list:
                mode_list.remove(mode)
        # end loop over modes to be ignored
    # end of if ignore_modes
    # If a selected list of frequencies has been selected then use these
    if len(my_modes) > 0:
        mode_list = list(set(my_modes))
    if debug:
        print("Selected mode list is: " + ",".join("{:4d}".format(m) for m in mode_list))
        print("")
    #
    print('Calculated Integrated IR band intensities (espilon_max assumes a fwhm line width)')
    print('Only IR active modes are shown')
    print()
    #      mode
    print("{:^4s} {:^8s} {:^17s} {:^27s} {:^11s} {:^8s}".format("Mode", "Frequency", "Intensity       ", "Integrated Molar Absorption", "Epsilon_max", "FWHM"))
    print("{:^4s} {:^8s} {:^17s} {:^27s} {:^11s} {:^8s}".format("    ", "cm-1     ", "Debye2/Angs2/amu", "Coefficient L/mole/cm/cm", "L/mole/cm", "cm-1"))
    print("{:^4s} {:^8s} {:^17s} {:^27s} {:^11s} {:^8s}".format("----", "---------", "----------------", "---------------------------", "-----------", "----"))
    if fd_csvfile is not None:
        print(command_line, file=fd_csvfile)
    fd_csvfile = fd_csvfiles[1]
    if fd_csvfile is not None:
        print("{:^4s},{:^8s},{:^17s},{:^27s},{:^11s},{:^8s}".format("Mode", "Frequency", "Intensity       ", "Integrated Molar Absorption", "Epsilon_max", "FWHM"), file=fd_csvfile)
        print("{:^4s},{:^8s},{:^17s},{:^27s},{:^11s},{:^8s}".format("    ", "cm-1     ", "Debye2/Angs2/amu", "Coefficient L/mole/cm/cm", "L/mole/cm", "cm-1"), file=fd_csvfile)
        print("{:^4s},{:^8s},{:^17s},{:^27s},{:^11s},{:^8s}".format("----", "---------", "----------------", "---------------------------", "---------", "----"), file=fd_csvfile)
    for mode in mode_list:
        freq = reader.frequencies[mode]
        intn = intensities[mode]
        widt = sigmas_cm1[mode]
        print("{:>4d} {:>8.2f} {:>17.6f} {:>27.6f} {:>11.4f} {:>8.1f}".format(mode, freq, intn, 4225.6*intn, 2*4225.6*intn/widt/PI, widt))
        if fd_csvfile is not None:
            print("{:>4d},{:>8.2f},{:>17.6f},{:>27.6f},{:>11.4f},{:>8.1f}".format(mode, freq, intn, 4225.6*intn, 2*4225.6*intn/widt/PI, widt), file=fd_csvfile)
    # from the oscillator strengths calculate the low frequency permittivity
    epsilon_ionic = Calculator.ionic_permittivity(mode_list, oscillator_strengths, frequencies, volume)
    print(" ")
    print3x3('Permittivity contribution from ionic displacements', epsilon_ionic)
    print3x3('Optical permittivity at zero frequency', epsilon_inf)
    epsilon_total = epsilon_inf + epsilon_ionic
    print3x3('Total permittivity', epsilon_total)
    # If we have been requested to write an xslx file do so now
    if fd_excelfile is not None:
        excel_row += 2
        for col,string in enumerate(["Mode", "Frequency", "Intensity", "Integrated Molar Absorption", "Epsilon Max", "FWHM"]):
            worksheet.write(excel_row,col,string)
        excel_row += 1
        for col,string in enumerate(["    ", "cm-1", "Debye2/Angs2/amu", "Coefficient L/mol/cm/cm", "L/mol/cm", "cm-1"]):
            worksheet.write(excel_row,col,string)
        for mode in mode_list:
            freq = reader.frequencies[mode]
            intn = intensities[mode]
            widt = sigmas_cm1[mode]
            excel_row += 1
            for col,result in enumerate([mode, freq, intn, 4225.6*intn, 2*4225.6*intn/widt/PI, widt]):
                worksheet.write(excel_row,col,result)
        excel_row += 2; 
        worksheet.write(excel_row,0,'Permittivity contribution from ionic displacements')
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_ionic[0]) ]
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_ionic[1]) ]
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_ionic[2]) ]
        #
        excel_row += 2;  worksheet.write(excel_row,0,'Optical permittivity at zero frequency')
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_inf[0]) ]
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_inf[1]) ]
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_inf[2]) ]
        #
        excel_row += 2; worksheet.write(excel_row,0,'Total permittivity')
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_total[0]) ]
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_total[1]) ]
        excel_row += 1; [ worksheet.write(excel_row,col,perm) for col,perm in enumerate(epsilon_total[2]) ]
        # Finished processing the first sheet of the xlsx file
    #
    # Prepare parallel call parameters for a Loop over frequencies
    call_parameters = []
    for v in np.arange(float(vmin), float(vmax)+0.5*float(increment), float(increment)):
        vau = v * wavenumber
        call_parameters.append( (v, vau, mode_list, frequencies, sigmas, oscillator_strengths, volume, epsilon_inf, drude, drude_plasma, drude_sigma) )
    # Create a pool of processors and get them to calculate the initial dielectric contribution from the phonon modes
    p = Pool(number_of_processors,initializer=set_affinity_on_worker,maxtasksperchild=10)
    results = p.map(Calculator.parallel_dielectric, call_parameters)
    # Prepare the share memory so we can pass previous solutions around
    # the number of plots will be less than or equal to nplots
    nplots = len(results)*len(methods)*len(volume_fractions)*len(shapes)
    # Allocate space for the shared memory, we need twice as much as we have a complex data type
    shared_array_base = Array(ctypes.c_double, 2*nplots*9)
    previous_solution_shared = np.ctypeslib.as_array(shared_array_base.get_obj())
    # Convert the space allocated to complex
    previous_solution_shared.dtype = np.complex128
    # Reshape the array and fill everything with zero's
    previous_solution_shared = previous_solution_shared.reshape(nplots, 3,3)
    previous_solution_shared.fill(0.0+0.0j)
    # Prepare parallel call parameters for the loop over frequencies, methods, volume fractions
    call_parameters = []
    for v,vau,dielecv in results:
        nplot = -1
        # Loop over methods
        for method in methods:
            # Loop over volume_fractions
            for vf, vf_type in zip(volume_fractions, fractional_types):
              # loop over sizes
              for size_mu,size_sigma in zip(sizes,size_distribution_sigmas):
                # convert the size to a dimensionless number which is 2*pi*size/lambda
                lambda_mu = 1.0E4 / (v + 1.0e-12)
                if size_mu < 1.0e-12:
                    size_mu = 1.0e-12
                size = 2.0*PI*size_mu / lambda_mu
                # Loop over shapes
                ap_count = 0  # Use a counter to only calculate the AveragedPermittivity once for all shapes
                for shape, data, L in zip(shapes, shape_data, depolarisations):
                    nplot += 1
                    if method == "ap" or method == "averagedpermittivity":
                        if ap_count == 0:
                            shape = "noshape"
                            ap_count += 1
                        else:
                            nplot -= 1
                            continue
                    call_parameters.append( (v,vau,dielecv,method,vf,vf_type,size_mu,size_sigma,size,nplot,dielectric_medium,shape,data,L,concentration,previous_solution_shared) )
                # end of loop over shapes
              # end of loop over sizes
            # end of loop over volume fractions
        # end of loops over methods
    # end of look over frequencies
    # Allocate some shared memory
    nplots = nplot
    # Now perform the parallel calculation of the effective medium equations
    results = p.map(Calculator.solve_effective_medium_equations, call_parameters)
    #  Close all the parallel workers
    p.close()
    p.join()
    # Add the results to the plotter function
    for v,nplot,method,vf_type,size_mu,size_sigma,shape,data,trace,absorption_coefficient,molar_absorption_coefficient in results:
        plotter.add_dielectric(nplot, method, vf_type, size_mu, size_sigma, shape, data, v, trace, absorption_coefficient, molar_absorption_coefficient)
    # Tidy up before plotting
    if fd_excelfile is not None:
        plotter.excel(workbook)
        workbook.close()
    if fd_csvfile is not None:
        fd_csvfile = fd_csvfiles[2]
        plotter.printout(fd_csvfile, csv_extended)
        for file_descriptor in fd_csvfiles:
            file_descriptor.close()
    # Call the plotter
    plotter.plot(plot_types)

if __name__ == "__main__":
    main()