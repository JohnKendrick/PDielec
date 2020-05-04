#!/usr/bin/env python
"""Read the contents of a directory containing DFT output and create a csv style file of information"""
from __future__ import print_function
import string
import re
import numpy as np
import os, sys
from PDielec.Constants import amu, PI, avogadro_si, wavenumber, angstrom, isotope_masses, average_masses
from PDielec.VaspOutputReader import VaspOutputReader
from PDielec.CastepOutputReader import CastepOutputReader
from PDielec.GulpOutputReader import GulpOutputReader
from PDielec.CrystalOutputReader import CrystalOutputReader
from PDielec.AbinitOutputReader import AbinitOutputReader
from PDielec.QEOutputReader import QEOutputReader
from PDielec.PhonopyOutputReader import PhonopyOutputReader
from multiprocessing import Pool, cpu_count
import PDielec.Calculator as Calculator

def set_affinity_on_worker():
    """When a new worker process is created, the affinity is set to all CPUs"""
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK for the time being this is simply commented out, but might be useful at some point
    #os.system("taskset -p 0xff %d > /dev/null" % os.getpid())

def read_a_file( calling_parameters):
    name, eckart, neutral, mass_definition, mass_dictionary, global_no_calculation, program, hessian_symmetrisation, qmprogram, debug = calling_parameters
    fulldirname = name
    head,tail = os.path.split(fulldirname)
    root,ext = os.path.splitext(tail)
    if program == "castep":
        names = [ name ]
        reader = CastepOutputReader( names )
    elif program == "vasp":
        name1 = os.path.join(head,'OUTCAR')
        name2 = os.path.join(head,'KPOINTS')
        names = [ name1, name2 ]
        reader = VaspOutputReader( names )
    elif program == "gulp":
        names = [ name ]
        reader = GulpOutputReader( names )
    elif program == "crystal":
        names = [ name ]
        reader = CrystalOutputReader( names )
    elif program == "abinit":
        names = [ name ]
        reader = AbinitOutputReader( names )
    elif program == "qe":
        tail1 = root+'.dynG'
        tail2 = root+'.log'
        tail3 = root+'.out'
        name1 = os.path.join(head,tail2)
        name2 = os.path.join(head,tail3)
        name3 = os.path.join(head,tail)
        # We want the dynG entry last rounding causes problems otherwise
        name4 = os.path.join(head,tail1)
        names = []
        for n in [ name1, name2, name3, name4 ]:
            if os.path.isfile(n):
                if not n in names:
                    names.append(n)
        reader = QEOutputReader( names )
    elif program == "phonopy":
        # The order is important
        pname1 = os.path.join(head,'qpoints.yaml')
        pname2 = os.path.join(head,'phonopy.yaml')
        # Only works for VASP at the moment
        vname1 = os.path.join(head,'OUTCAR')
        vname2 = os.path.join(head,'KPOINTS')
        pnames = [ pname1, pname2 ]
        vnames = [ vname1, vname2 ]
        pnames.extend(vnames)
        names = pnames
        # Which QM program was used by PHONOPY?
        if qmprogram == "castep":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = CastepOutputReader(names)
        elif qmprogram == "vasp":
            qmreader = VaspOutputReader(vnames)
        elif qmprogram == "gulp":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = GulpOutputReader(names)
        elif qmprogram == "crystal":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = CrystalOutputReader(names)
        elif qmprogram == "abinit":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = AbinitOutputReader(names)
        elif qmprogram == "qe":
            print("Error in qmreader",qmprogram)
            exit()
            qmreader = QEOutputReader(names)
        # The QM reader is used to get info about the QM calculation
        reader = PhonopyOutputReader(pnames,qmreader)
    else:
        print('Program name not recognized',program,file=sys.stderr)
        exit()

    print('  Analysing ',names, file=sys.stderr)
    # The order that the settings are applied is important
    # Symmetry is applied before the file is read, it is an integral part of the hessian
    # Eckart and neutral are applied after the file has been read, this way the original frequencies are those before any calculations
    reader.debug = debug
    reader.hessian_symmetrisation = hessian_symmetrisation
    reader.read_output()
    frequencies_cm1 = reader.frequencies
    frequencies = np.array(frequencies_cm1)
    no_calculation = global_no_calculation
    if len(frequencies) < 3:
        no_calculation = True
    mode_list = []
    ignore_modes = []
    sigmas = []
    if not no_calculation:
        # apply the eckart conditions before we change the masses
        reader.eckart = eckart
        # Get the born charges
        if neutral:
            reader.neutralise_born_charges()
        # What mass definition are we using?
        if not mass_definition == "program" or mass_dictionary:
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
                elif np.real(modified_frequencies_cm1[mode])/wavenumber < threshold_frequency:
                    ignore_modes.append(mode)
                # ignore modes with imaginary frequency
                elif abs(np.imag(modified_frequencies_cm1[mode]))/wavenumber > 1.0e-6:
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
        ionicv = Calculator.dielectric_contribution(0.0, mode_list, modified_frequencies*wavenumber, sigmas, oscillator_strengths, volume)
    # absorption units here are L/mole/cm-1
    # Continue reading any data from the output file
    frequencies_cm1.sort()
    unitCell = reader.unit_cells[-1]
    a,b,c,alpha,beta,gamma = unitCell.convert_unitcell_to_abc()
    if not no_calculation:
        eps0 = np.real(ionicv)
    else:
        eps0   = np.array(reader.zerof_static_dielectric)
    epsinf = np.array(reader.zerof_optical_dielectric)
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
    header = fulldirname + ',Read from program'
    string = ''
    string = string + ',' + str(reader.electrons)
    string = string + ',' + str(reader.magnetization)
    string = string + ',' + str(reader.kpoints)
    string = string + ',' + str(reader.kpoint_grid[0])
    string = string + ',' + str(reader.kpoint_grid[1])
    string = string + ',' + str(reader.kpoint_grid[2])
    string = string + ',' + str(reader.energy_cutoff)
    string = string + ',' + str(reader.final_free_energy)
    string = string + ',' + str(reader.final_energy_without_entropy)
    string = string + ',' + str(reader.pressure)
    string = string + ',' + str(a) + ',' + str(b) + ',' + str(c) + ',' + str(alpha) + ',' + str(beta) + ',' + str(gamma) + ',' + volume
    string = string + ',' + eps0_xx + ',' + eps0_yy + ',' + eps0_zz + ',' + eps0_xy + ',' + eps0_xz + ',' + eps0_yz
    string = string + ',' + epsinf_xx + ',' + epsinf_yy + ',' + epsinf_zz + ',' + epsinf_xy + ',' + epsinf_xz + ',' + epsinf_yz
    string = string + ',' + c11 + ',' + c22 + ',' + c33 + ',' + c44 + ',' + c55 + ',' + c66
    string = string + ',' + c12 + ',' + c13 + ',' + c23
    common_output=string
    string = header+common_output
    for f in frequencies_cm1:
        string = string + ',' + str(f)
    results_string = [ string ]
    # Assemble the next line if any of eckart/neutral/symm have been used
    option_string = ''
    if not no_calculation:
        option_string = 'Calculated frequencies (cm-1)'
        if eckart :
            option_string = option_string+' eckart'
        if neutral :
            option_string = option_string+' neutral'
        option_string = option_string+' hessian_symmetrisation='+hessian_symmetrisation
        option_string = option_string+' mass_definition='+mass_definition
        header = fulldirname+','+option_string
        string = header + common_output
        for f in modified_frequencies_cm1:
            string = string + ',' + str(f)
        results_string.append(string)
        option_string = 'Calculated Intensities (Debye2/Angs2/amu)'
        header = fulldirname+','+option_string
        string = header + common_output
        for f in intensities:
            string = string + ',' + str(f)
        results_string.append(string)
        option_string = 'Calculated Integrated Molar Absorption (L/mole/cm/cm)'
        header = fulldirname+','+option_string
        string = header + common_output
        for f in intensities:
            string = string + ',' + str(f*4225.6)
        results_string.append(string)
    # End if not no_calculation
    return name,results_string

def main():
    # Start processing the directories
    if len(sys.argv) <= 1 :
        print('preader -program program [-eckart] [-neutral] [-nocalculation] [-hessian symm] [-masses average] filenames .....', file=sys.stderr)
        print('  \"program\" must be one of \"abinit\", \"castep\", \"crystal\", \"gulp\"       ', file=sys.stderr)
        print('           \"phonopy\", \"qe\", \"vasp\"                                         ', file=sys.stderr)
        print('           If phonopy is used it must be followed by the QM package              ', file=sys.stderr)
        print('  -masses [average|isotope|program]  chooses the atomic mass definition average  ', file=sys.stderr)
        print('          The default is \"average\"                                             ', file=sys.stderr)
        print('  -mass  element mass                                                            ', file=sys.stderr)
        print('          Change the mass of the element specified                               ', file=sys.stderr)
        print('  -eckart projects out the translational components. Default is no projection    ', file=sys.stderr)
        print('  -neutral imposes neutrality on the Born charges of the molecule.               ', file=sys.stderr)
        print('           Default is no neutrality is imposed                                   ', file=sys.stderr)
        print('  -hessian [symm|crystal] defines the symmetrisation of the hessian              ', file=sys.stderr)
        print('           \"crystal\" imposes Crystal14 symmetrisation                          ', file=sys.stderr)
        print('           \"symm\" symmetrises by averaging the hessian with its transpose      ', file=sys.stderr)
        print('           \"symm\" is the default                                               ', file=sys.stderr)
        print('  -nocalculation requests no calculations are performed                          ', file=sys.stderr)
        print('           A single line is output with results obtained by reading the output   ', file=sys.stderr)
        print('           any of -mass -masses -eckart -neutral or -crystal are ignored         ', file=sys.stderr)
        print('  -debug   to switch on more debug information                                   ', file=sys.stderr)
        exit()

    observables=False
    files = []
    eckart = False
    neutral = False
    hessian_symmetrisation = "symm"
    mass_definition = "average"
    mass_dictionary = {}
    global_no_calculation = False
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    program = ''
    qmprogram = ''
    debug = False
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token == "-observables":
            observables = True
        elif token == "-debug":
            debug = True
        elif token == "-eckart":
            eckart = True
        elif token == "-neutral":
            neutral = True
        elif token == "-hessian":
            itoken += 1
            hessian_symmetrisation = tokens[itoken]
            if (not hessian_symmetrisation == "symm") and ( not hessian_symmetrisation == "crystal" ):
                print('-hessian must be followed by \"symm\" or \"crystal\"',file=sys.stderr)
                exit()
        elif token == "-masses":
            itoken += 1
            token = tokens[itoken]
            if token[0] == 'a':
                mass_definition = 'average'
            elif token[0] == 'i':
                mass_definition = 'isotope'
            elif token[0] == 'p':
                mass_definition = 'program'
            else:
                print('-masses must be followed by \"average\" \"isotope\" or \"program\"',file=sys.stderr)
                exit()
        elif token == "-mass":
            itoken += 1
            element = tokens[itoken]
            itoken += 1
            mass = float(tokens[itoken])
            mass_dictionary[element] = mass
        elif token == "-nocalculation":
            global_no_calculation = True
        elif token == "-program":
            itoken += 1
            program = tokens[itoken]
            if program == 'phonopy':
                itoken += 1
                qmprogram = tokens[itoken]
        else:
            files.append(token)

    if len(program) < 1:
        print('Please use -program to define the package used to generate the output files',file=sys.stderr)
        exit()

    if not program in ['abinit','castep','crystal','gulp','qe','vasp','phonopy']:
        print('Program is not recognised: ',program,file=sys.stderr)
        exit()

    if program == 'phonopy':
        if not qmprogram in ['abinit','castep','crystal','gulp','qe','vasp']:
            print('Phonopy QM program is not recognised: ',qmprogram,file=sys.stderr)
            exit()
        print('  QM program used by Phonopy is: ',qmprogram,file=sys.stderr)

    print('  Program is ',program,file=sys.stderr)

    for f in files:
        if not os.path.isfile(f):
            print('Error file requested for analysis does not exist',f,file=sys.stderr)
            exit()

    if global_no_calculation:
        print('  No calculations has been requested',file=sys.stderr)
        print('  The frequencies (if present) will come from the QM/MM program    ',file=sys.stderr)
        print('  No hessian symmetrisation be performed',file=sys.stderr)
        if eckart:
            print('  The -eckart flag will be ignored',file=sys.stderr)
            eckart = False
        if neutral:
            print('  The -neutral flag will be ignored',file=sys.stderr)
            neutral = False
        if not mass_definition == "program":
            print('  The average and isotope mass definitions will be ignored',file=sys.stderr)
            mass_definition = "program"
        if mass_dictionary:
            print('  The masses specified by the -mass directive will be ignored',file=sys.stderr)
            mass_dictionary = {}
    else:
        print('  Eckart is ',eckart,file=sys.stderr)
        print('  Neutral is ',neutral,file=sys.stderr)
        print('  Mass definition is ',mass_definition,file=sys.stderr)
        print('  Hessian symmetrisation is ',hessian_symmetrisation,file=sys.stderr)
    #
    # Create a pool of processors to handle reading the files
    #
    p = Pool(initializer=set_affinity_on_worker)
    # Create a tuple list of calling parameters
    calling_parameters = []
    files.sort()
    for name in files:
        calling_parameters.append( (name, eckart, neutral, mass_definition, mass_dictionary, global_no_calculation, program, hessian_symmetrisation, qmprogram, debug) )
    # Calculate the results in parallel
    results_map_object = p.map_async(read_a_file,calling_parameters)
    results_map_object.wait()
    results = results_map_object.get()
    # Convert the results into a dictionary
    results_dictionary = {}
    for name,strings in results:
        results_dictionary[name] = strings
    # Print out the results
    print('directory,information,electrons,magnetization,kpnts,kpnt_1,kpnt_2,kpnt_3,energy_cutoff_eV,final_free_energy_eV,final_energy_without_entropy_eV,pressure_GPa,a_A,b_A,c_A,alpha,beta,gamma,volume,eps0_xx,eps0_yy,eps0_zz,eps0_xy,eps0_xz,eps0_yz,epsinf_xx,epsinf_yy,epsinf_zz,epsinf_xy,epsinf_xz,epsinf_yz, c11_gpa, c22_gpa, c33_gpa, c44_gpa, c55_gpa, c66_gpa, c12_gpa,  c13_gpa, c23_gpa,f1,f2,f3,f4,f5,f6....')
    for name in files:
        for string in results_dictionary[name]:
            print(string)
    #
    # Process the observables switch
    # ( for the time being this is commented out )
    exit()
    for name in files :
        if observables:
            filename = directory+'.observables'
            if directory == '.' :
                filename = 'dot.observables'
            filename = filename.replace("/","_")
            print('filename', filename, file=sys.stderr)
            fd = open(filename,'w')
            print('observables', file=fd)
            if eps0[0,0] != 0.0 :
                print('sdlc', file=fd)
                print('1 1 ', eps0[0,0], file=fd)
                print('sdlc', file=fd)
                print('2 2 ', eps0[1,1], file=fd)
                print('sdlc', file=fd)
                print('3 3 ', eps0[2,2], file=fd)
                print('hfdlc', file=fd)
                print('1 1 ', epsinf[0,0], file=fd)
                print('hfdlc', file=fd)
                print('2 2 ', epsinf[1,1], file=fd)
                print('hfdlc', file=fd)
                print('3 3 ', epsinf[2,2], file=fd)
            if carray[0,0] != 0.0:
                print('elastic', file=fd)
                print('1 1 ', carray[0,0], file=fd)
                print('elastic', file=fd)
                print('2 2 ', carray[1,1], file=fd)
                print('elastic', file=fd)
                print('3 3 ', carray[2,2], file=fd)
                print('elastic', file=fd)
                print('4 4 ', carray[3,3], file=fd)
                print('elastic', file=fd)
                print('5 5 ', carray[4,4], file=fd)
                print('elastic', file=fd)
                print('6 6 ', carray[5,5], file=fd)
                print('elastic', file=fd)
                print('1 2 ', carray[0,1], file=fd)
                print('elastic', file=fd)
                print('1 3 ', carray[0,2], file=fd)
                print('elastic', file=fd)
                print('2 3 ', carray[1,2], file=fd)
            nfreq = len(frequencies_cm1)
            if nfreq > 4 :
                for n in range(nfreq):
                    nlast = n
                    if frequencies_cm1[n].real > 1.0E-6:
                        break
                    # end of if
                # end of for
                if frequencies_cm1[nlast].real > 1.0E-8 :
                    n = nlast
                    print('frequency', nfreq - nlast, file=fd)
                    for f in frequencies_cm1[nlast:] :
                        n = n + 1
                        print(n, f, file=fd)
            print('end', file=fd)
            fd.close()
    # End of for loop over files
# end of def main

if __name__ == "__main__":
    main(sys)
