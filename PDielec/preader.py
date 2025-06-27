#!/usr/bin/env python
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
"""Read the contents of a directory containing DFT output and create a csv style file of information.

::

    preader -program program [-eckart] [-neutral] [-nocalculation] [-masses average] [-pickle name] [-version] filenames .....
      "program" must be one of "abinit", "castep", "crystal", "gulp"
               "phonopy", "qe", "vasp", "experiment", "auto"
               The default is auto, so the program tries to guess the package from
               the contents of the directory.  However this is not fool-proof!
      -masses [average|isotope|program]  chooses the atomic mass definition average
              The default is "average"
      -mass  element mass
              Change the mass of the element specified
      -eckart projects out the translational components. Default is no projection
      -neutral imposes neutrality on the Born charges of the molecule.
               Default is no neutrality is imposed
      -nocalculation requests no calculations are performed
               A single line is output with results obtained by reading the output
               any of -mass -masses -eckart -neutral or -crystal are ignored
      -debug   to switch on more debug information
      -pickle  write each file reader to a pickled dump file for later processing
               the name of the file to hold all the pickled readers is given
               If the file exists it is not overwritten
      -version print the version of PDielec library being used
"""
import os
import sys
from multiprocessing.dummy import Pool

import dill as pickle
import numpy as np
import psutil

import PDielec.__init__
from PDielec import Calculator, DielectricFunction, Utilities
from PDielec.Constants import amu, angstrom, average_masses, isotope_masses, wavenumber
from PDielec.Utilities import find_program_from_name

version = PDielec.__init__.__version__

def set_affinity_on_worker():
    """When a new worker process is created, the affinity is set to all CPUs."""
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK for the time being this is simply commented out, but might be useful at some point
    #os.system("taskset -p 0xff %d > /dev/null" % os.getpid())

def read_a_file( calling_parameters):
    """Read data from a file and process it based on the given parameters.

    Parameters
    ----------
    calling_parameters : tuple
        A tuple that contains various parameters used for reading and processing the file. This includes:
        - name : str (Name of the file.)
        - eckart : bool (Specifies whether Eckart conditions should be applied.)
        - neutral : bool (Specifies whether to neutralise Born charges.)
        - mass_definition : str (Specifies the method of mass calculation ('average', 'isotope', or 'program').)
        - mass_dictionary : dict (A dictionary containing mass information when `mass_definition` is not 'program'.)
        - global_no_calculation : bool (A flag to indicate whether any calculations should be performed.)
        - program : str (Name of the computational chemistry program used to generate the file.)
        - debug : bool (If set to True, additional debug information will be printed.)

    Returns
    -------
    tuple
        A tuple containing:
        - reader : object (An object capable of reading and processing the data from the specified file.)
        - name : str (Name of the file that was processed.)
        - results_string : list (A list of strings that represent the processed data ready for output. Depending on conditions, it includes the initial and modified frequencies, intensities, and optionally calculated molar absorption rates.)

    """    
    name, eckart, neutral, mass_definition, mass_dictionary, global_no_calculation, program, debug = calling_parameters
    reader = Utilities.get_reader(name,program)
    # The order that the settings are applied is important
    # Eckart and neutral are applied after the file has been read, this way the original frequencies are those before any calculations
    reader.debug = debug
    reader.read_output()
    frequencies_cm1 = reader.frequencies
    frequencies = np.array(frequencies_cm1)
    no_calculation = global_no_calculation
    if len(frequencies) < 3:
        no_calculation = True
    mode_list = []
    ignore_modes = []
    sigmas = []
    epsinf = np.array(reader.zerof_optical_dielectric)
    if not no_calculation:
        # apply the eckart conditions before we change the masses
        reader.eckart = eckart
        # Get the born charges
        if neutral:
            reader.neutralise_born_charges()
        # What mass definition are we using?
        if mass_definition != "program" or mass_dictionary:
            if mass_definition == "average":
                reader.change_masses(average_masses, mass_dictionary)
            elif mass_definition == "isotope":
                reader.change_masses(isotope_masses, mass_dictionary)
        born_charges = np.array(reader.born_charges)
        # Calculate the mass weighted normal modes.  This just forces projection
        mass_weighted_normal_modes = reader.calculate_mass_weighted_normal_modes()
        # The masses might have changed in the calculation of the mass weighted normal modes
        masses = np.array(reader.masses) * amu
        normal_modes = Calculator.normal_modes(masses, mass_weighted_normal_modes)
        # from the normal modes and the born charges calculate the oscillator strengths of each mode
        oscillator_strengths = Calculator.oscillator_strengths(normal_modes, born_charges)
        modified_frequencies_cm1 = reader.frequencies
        modified_frequencies_cm1.sort()
        modified_frequencies = np.array(modified_frequencies_cm1)
        # if the frequency is less than 5 cm-1 assume that the oscillator strength is zero
        for imode,f in enumerate(modified_frequencies_cm1):
            mode_list.append(imode)
            sigmas.append(5.0)
            if f < 5.0:
                oscillator_strengths[imode]=np.zeros((3,3))
        # calculate the intensities from the trace of the oscillator strengths
        intensities = Calculator.infrared_intensities(oscillator_strengths)
        # Calculate eps0
        volume = reader.volume*angstrom*angstrom*angstrom
        #
        # Calculate degenerate lists
        #
        degeneracy_threshold = 1.0E-8
        threshold_intensity = 1.0E-6
        threshold_frequency = 5.0
        degenerate_lists = {}
        mmax = len(modified_frequencies_cm1)
        for m1 in range(len(modified_frequencies_cm1)):
            degenerate_lists[m1] = []
        for m1 in range(len(modified_frequencies_cm1)):
            f1 = modified_frequencies[m1]
            for m2 in range(m1+1,min(mmax,m1+3)):
                f2 = modified_frequencies_cm1[m2]
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
                if intensity < threshold_intensity:
                    # If any of its degenerate modes have intensity then we shouldn't ignore it
                    ignore = True
                    for m in degenerate_lists[mode]:
                        if intensities[m] > threshold_intensity:
                            ignore = False
                    if ignore:
                        ignore_modes.append(mode)
                # ignore modes with low real frequency
                elif np.real(modified_frequencies_cm1[mode])/wavenumber < threshold_frequency or abs(np.imag(modified_frequencies_cm1[mode]))/wavenumber > 1.0e-6:
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
        crystalPermittivity = DielectricFunction.DFT(mode_list, modified_frequencies*wavenumber, sigmas, oscillator_strengths, volume, False, 0.0, 0.0) 
        crystalPermittivity.setEpsilonInfinity(epsinf)
        ionicv = crystalPermittivity.calculate(0.0) - epsinf
    # absorption units here are L/mole/cm-1
    # Continue reading any data from the output file
    frequencies_cm1.sort()
    unitCell = reader.get_unit_cell()
    if unitCell is None:
        a,b,c,alpha,beta,gamma = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    else:
        a,b,c,alpha,beta,gamma = unitCell.convert_unitcell_to_abc()
    eps0   = np.array(reader.zerof_static_dielectric) if no_calculation else np.real(ionicv)
    eps0_xx = str(eps0[0,0])
    eps0_yy = str(eps0[1,1])
    eps0_zz = str(eps0[2,2])
    eps0_xy = str(( eps0[0,1] + eps0[1,0] ) /2.0)
    eps0_xz = str(( eps0[0,2] + eps0[2,0] ) /2.0)
    eps0_yz = str(( eps0[1,2] + eps0[2,1] ) /2.0)
    epsinf_xx = str(epsinf[0,0])
    epsinf_yy = str(epsinf[1,1])
    epsinf_zz = str(epsinf[2,2])
    epsinf_xy = str(( epsinf[0,1] + epsinf[1,0] ) /2.0)
    epsinf_xz = str(( epsinf[0,2] + epsinf[2,0] ) /2.0)
    epsinf_yz = str(( epsinf[1,2] + epsinf[2,1] ) /2.0)
    volume = str(reader.volume)
    carray = np.array(reader.elastic_constants)
    c11 = str(carray[0,0])
    c22 = str(carray[1,1])
    c33 = str(carray[2,2])
    c44 = str(carray[3,3])
    c55 = str(carray[4,4])
    c66 = str(carray[5,5])
    c12 = str(carray[0,1])
    c13 = str(carray[0,2])
    c23 = str(carray[1,2])
    # Assemble the output line for the unprojected frequencies and all the other data
    header = name + ",Read from program"
    string = ""
    string = string + "," + str(reader.electrons)
    string = string + "," + str(reader.magnetization)
    string = string + "," + str(reader.kpoints)
    string = string + "," + str(reader.kpoint_grid[0])
    string = string + "," + str(reader.kpoint_grid[1])
    string = string + "," + str(reader.kpoint_grid[2])
    string = string + "," + str(reader.energy_cutoff)
    string = string + "," + str(reader.final_free_energy)
    string = string + "," + str(reader.final_energy_without_entropy)
    string = string + "," + str(reader.pressure)
    string = string + "," + str(a) + "," + str(b) + "," + str(c) + "," + str(alpha) + "," + str(beta) + "," + str(gamma) + "," + volume
    string = string + "," + eps0_xx + "," + eps0_yy + "," + eps0_zz + "," + eps0_xy + "," + eps0_xz + "," + eps0_yz
    string = string + "," + epsinf_xx + "," + epsinf_yy + "," + epsinf_zz + "," + epsinf_xy + "," + epsinf_xz + "," + epsinf_yz
    string = string + "," + c11 + "," + c22 + "," + c33 + "," + c44 + "," + c55 + "," + c66
    string = string + "," + c12 + "," + c13 + "," + c23
    common_output=string
    string = header+common_output
    for f in frequencies_cm1:
        string = string + "," + str(f)
    results_string = [ string ]
    # Assemble the next line if any of eckart/neutral have been used
    option_string = ""
    if not no_calculation:
        option_string = "Calculated frequencies (cm-1)"
        if eckart :
            option_string = option_string+" eckart"
        if neutral :
            option_string = option_string+" neutral"
        option_string = option_string+" mass_definition="+mass_definition
        header = name+","+option_string
        string = header + common_output
        for f in modified_frequencies_cm1:
            string = string + "," + str(f)
        results_string.append(string)
        option_string = "Calculated Intensities (Debye2/Angs2/amu)"
        header = name+","+option_string
        string = header + common_output
        for f in intensities:
            string = string + "," + str(f)
        results_string.append(string)
        option_string = "Calculated Integrated Molar Absorption (L/mole/cm/cm)"
        header = name+","+option_string
        string = header + common_output
        for f in intensities:
            string = string + "," + str(f*4225.6)
        results_string.append(string)
    # End if not no_calculation
    return reader,name,results_string

def print_help():
    """Display the help message for the preader program.

    This function prints out the usage of the preader program, including the supported flags and their descriptions.

    Notes
    -----
    The function internally directs its output to `sys.stderr` instead of the default `sys.stdout`, to separate the help message from other outputs. It also exits the program after displaying the help message by calling `exit()`. 

    Examples
    --------
    To display the help message, simply call the function without any arguments:
    
        >>> print_help()

    No output is returned by this function as it directs its message to `sys.stderr` and exits the program.

    """    
    print("preader -program program [-eckart] [-neutral] [-nocalculation] [-masses average] [-pickle name] [-version] filenames .....", file=sys.stderr)
    print('  "program" must be one of "abinit", "castep", "crystal", "gulp"       ', file=sys.stderr)
    print('           "phonopy", "qe", "vasp", "experiment", "auto"               ', file=sys.stderr)
    print("           The default is auto, so the program tries to guess the package from   ", file=sys.stderr)
    print("           the contents of the directory.  However this is not fool-proof!       ", file=sys.stderr)
    print("  -masses [average|isotope|program]  chooses the atomic mass definition average  ", file=sys.stderr)
    print('          The default is "average"                                             ', file=sys.stderr)
    print("  -mass  element mass                                                            ", file=sys.stderr)
    print("          Change the mass of the element specified                               ", file=sys.stderr)
    print("  -eckart projects out the translational components. Default is no projection    ", file=sys.stderr)
    print("  -neutral imposes neutrality on the Born charges of the molecule.               ", file=sys.stderr)
    print("           Default is no neutrality is imposed                                   ", file=sys.stderr)
    print("  -nocalculation requests no calculations are performed                          ", file=sys.stderr)
    print("           A single line is output with results obtained by reading the output   ", file=sys.stderr)
    print("           any of -mass -masses -eckart -neutral or -crystal are ignored         ", file=sys.stderr)
    print("  -debug   to switch on more debug information                                   ", file=sys.stderr)
    print("  -pickle  write each file reader to a pickled dump file for later processing    ", file=sys.stderr)
    print("           the name of the file to hold all the pickled readers is given         ", file=sys.stderr)
    print("           If the file exists it is not overwritten                              ", file=sys.stderr)
    print("  -version print the version of PDielec library being used                       ", file=sys.stderr)
    print("  Version ",version,file=sys.stderr)
    sys.exit()


def main():
    """Process command line arguments.

    ::

        preader -program program [-eckart] [-neutral] [-nocalculation] [-masses average] [-pickle name] [-version] filenames .....
          "program" must be one of "abinit", "castep", "crystal", "gulp"
                   "phonopy", "qe", "vasp", "experiment", "auto"
                   The default is auto, so the program tries to guess the package from
                   the contents of the directory.  However this is not fool-proof!
          -masses [average|isotope|program]  chooses the atomic mass definition average
                  The default is "average"
          -mass  element mass
                  Change the mass of the element specified
          -eckart projects out the translational components. Default is no projection
          -neutral imposes neutrality on the Born charges of the molecule.
                   Default is no neutrality is imposed
          -nocalculation requests no calculations are performed
                   A single line is output with results obtained by reading the output
                   any of -mass -masses -eckart -neutral or -crystal are ignored
          -debug   to switch on more debug information
          -pickle  write each file reader to a pickled dump file for later processing
                   the name of the file to hold all the pickled readers is given
                   If the file exists it is not overwritten
          -version print the version of PDielec library being used

    """    
    # Start processing the directories
    if len(sys.argv) <= 1 :
        print_help()
    files = []
    eckart = False
    neutral = False
    mass_definition = "average"
    mass_dictionary = {}
    global_no_calculation = False
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    program = "auto"
    debug = False
    picklefile = None
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        token = token.replace("--","-")
        if token == "-debug":
            debug = True
        elif token == "-eckart":
            eckart = True
        elif token == "-neutral":
            neutral = True
        elif token == "-masses":
            itoken += 1
            token = tokens[itoken]
            if token[0] == "a":
                mass_definition = "average"
            elif token[0] == "i":
                mass_definition = "isotope"
            elif token[0] == "p":
                mass_definition = "program"
            else:
                print('-masses must be followed by "average" "isotope" or "program"',file=sys.stderr)
                sys.exit()
        elif token == "-mass":
            itoken += 1
            element = tokens[itoken]
            itoken += 1
            mass = float(tokens[itoken])
            mass_dictionary[element] = mass
        elif token == "-help":
            print_help()
        elif token == "-version":
            print("  Version ",version,file=sys.stderr)
            sys.exit()
        elif token == "-nocalculation":
            global_no_calculation = True
        elif token == "-pickle":
            itoken += 1
            picklefile = tokens[itoken]
            if os.path.isfile(picklefile):
                print(f"Error: pickle file {picklefile} already exists",file=sys.stderr)
                sys.exit()
        elif token == "-program":
            itoken += 1
            program = tokens[itoken]
        else:
            files.append(token)

    if len(program) < 1:
        print("Please use -program to define the package used to generate the output files",file=sys.stderr)
        sys.exit()

    if program == "qe":
        program = "quantum espresso"

    if program not in ["auto","abinit","castep","crystal","gulp","quantum espresso","vasp","phonopy","experiment"]:
        print("Program is not recognised: ",program,file=sys.stderr)
        sys.exit()

    print("  Program is ",program,file=sys.stderr)

    for f in files:
        if not os.path.isfile(f):
            print("Error file requested for analysis does not exist",f,file=sys.stderr)
            sys.exit()

    if global_no_calculation:
        print("  No calculations has been requested",file=sys.stderr)
        print("  The frequencies (if present) will come from the QM/MM program    ",file=sys.stderr)
        if eckart:
            print("  The -eckart flag will be ignored",file=sys.stderr)
            eckart = False
        if neutral:
            print("  The -neutral flag will be ignored",file=sys.stderr)
            neutral = False
        if mass_definition != "program":
            print("  The average and isotope mass definitions will be ignored",file=sys.stderr)
            mass_definition = "program"
        if mass_dictionary:
            print("  The masses specified by the -mass directive will be ignored",file=sys.stderr)
            mass_dictionary = {}
    else:
        print("  Eckart is ",eckart,file=sys.stderr)
        print("  Neutral is ",neutral,file=sys.stderr)
        print("  Mass definition is ",mass_definition,file=sys.stderr)
    #
    # Create a pool of processors to handle reading the files
    #
    number_of_processors = psutil.cpu_count(logical=False)
    p = Pool(number_of_processors,initializer=set_affinity_on_worker)
    # Create a tuple list of calling parameters
    calling_parameters = []
    files.sort()
    for name in files:
        prog = program
        if program == "auto":
            prog = find_program_from_name(name)
            print(f"  Analysing {name} generated by {prog}",file=sys.stderr)
        calling_parameters.append( (name, eckart, neutral, mass_definition, mass_dictionary, global_no_calculation, prog, debug) )
    # Calculate the results in parallel
    results_map_object = p.map_async(read_a_file,calling_parameters)
    results_map_object.wait()
    results = results_map_object.get()
    # Convert the results into a dictionary
    results_dictionary = {}
    for reader,name,strings in results:
        results_dictionary[name] = strings
        # If pickling then open the pickle file
        if picklefile:
            with open(picklefile,"wb") as picklefd:
                pickle.dump(reader,picklefd)
    # Print out the results
    print("directory,information,electrons,magnetization,kpnts,kpnt_1,kpnt_2,kpnt_3,energy_cutoff_eV,final_free_energy_eV,final_energy_without_entropy_eV,pressure_GPa,a_A,b_A,c_A,alpha,beta,gamma,volume,eps0_xx,eps0_yy,eps0_zz,eps0_xy,eps0_xz,eps0_yz,epsinf_xx,epsinf_yy,epsinf_zz,epsinf_xy,epsinf_xz,epsinf_yz, c11_gpa, c22_gpa, c33_gpa, c44_gpa, c55_gpa, c66_gpa, c12_gpa,  c13_gpa, c23_gpa,f1,f2,f3,f4,f5,f6....")
    for name in files:
        for string in results_dictionary[name]:
            print(string)
    #
    # Close the pool down
    #
    p.close()
    p.join()
# end of def main

if __name__ == "__main__":
    main()
