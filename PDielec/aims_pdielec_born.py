#!/usr/bin/env python3
"""Read the FHI Aims output after a calculation of intensities using get_vibrations.py.

Write the BORN charge file in PDielec format
"""
import optparse
import sys

import numpy as np
from scipy.constants import value

USAGE = """%prog [options] <name> 
<name> will be used as prefix for all output files.
Calculate Born charges and store in BORN_PDIELEC
The file structure is set up by get_vibrations.py.
The electronic permittivity is calculated in the Dielectric/ directory
"""

def find_number_of_atoms():
    """Find the number of atoms in the geometry.in file.

    Returns
    -------
    n_atoms : int
        The number of atoms in a geometry.in file

    """
    # Determine the number of atoms
    n_atoms = 0
    with open("geometry.in") as geo:
        for line in geo:
            if line.startswith("atom") :
                n_atoms += 1
    return n_atoms

def read_aims_output(folder, aimsoutname, n_atoms):
    """Read an Aims output file and find the volume and polarisation.

    Parameters
    ----------
    folder : str
        The directory holding the Aims file
    aimsoutname : str
         The name of the Aims file, usually aims.out
    n_atoms : int
         The number of atoms

    Returns
    -------
    volume : float
        The volume in angstrom^3
    pol_jump : np array of 3
        The polarizability
    polr_jump : np array of 3
        The polarization
    dielectric : np array of 3x3
        The permitivity epsilon infinity

    """
    volume = -9999
    pol_jump = -9999
    polr_jump = -9999
    dielectric = np.zeros( (3,3) )
    dielectric_row = -9999
    
    with open(folder+"/"+aimsoutname) as data:
        
        for line in data:
            if "| Unit cell volume " in line:
                volume = np.float64(line.split()[-2])
                continue

            if "Polarizability" in line:
                pol_jump = np.float64(line.split()[-6:])  # Periodic/cluster
                continue

            if "Cartesian Polarization" in line:
                polr_jump = np.float64(line.split()[-3:])  # Periodic
                continue

            if "DFPT for dielectric_constant:" in line:
                dielectric_row = 0
                continue

            if dielectric_row == 0:
                dielectric[0,0] = np.float64(line.split()[0])  # Periodic
                dielectric[0,1] = np.float64(line.split()[1])  # Periodic
                dielectric[0,2] = np.float64(line.split()[2])  # Periodic
                dielectric_row = 1
                continue

            if dielectric_row == 1:
                dielectric[1,0] = np.float64(line.split()[0])  # Periodic
                dielectric[1,1] = np.float64(line.split()[1])  # Periodic
                dielectric[1,2] = np.float64(line.split()[2])  # Periodic
                dielectric_row = 2
                continue

            if dielectric_row == 2:
                dielectric[2,0] = np.float64(line.split()[0])  # Periodic
                dielectric[2,1] = np.float64(line.split()[1])  # Periodic
                dielectric[2,2] = np.float64(line.split()[2])  # Periodic
                dielectric_row = -9999
                continue

    # Symmetrise the dielectric constant
    dielectric = 0.5 * ( dielectric + dielectric.T ) 
    return volume, pol_jump, polr_jump, dielectric

def main():
  """Generate Born Charge file from Aims to PDielec."""
  # Parse command line
  parser = optparse.OptionParser(usage=USAGE)
  parser.add_option("-d", "--delta", action="store", type="float",
                    help="Displacement (default 0.0025)", default=0.0025)
  parser.add_option("-p", "--permittivity", action="store", type="str",
                    help="The directory containing the electronic permittivity calculation", default="Dielectric")
  parser.add_option("-i", "--info", action="store_true",
                      help="Calculate Born charges and store in BORN_PDIELEC")
  options, args = parser.parse_args()
  if options.info:
      print(__doc__)
      sys.exit(0)
  if len(args) != 1:
      parser.error("Need exactly one arguments")

  name=args[0]
  delta=options.delta

  # Constants
  Ang = 1.0e-10
  elementary_charge = value("elementary charge")

  # File names
  bornchargefile = "BORN_PDIELEC"
  aimsoutname = "aims.out"

  # Define deltas and coefficients
  deltas=np.array([-delta,delta])
  coeff=np.array([-1,1])
  # define denominator in metres
  c_zero = 1. / (2. * Ang * delta)
  n_atoms = find_number_of_atoms()
  # Read in the electronic dielectric constant
  _volume, _pol_jump, _polr_jump, electronic_dielectric = read_aims_output(options.permittivity, aimsoutname, n_atoms)
  print("Electronic permittivity")
  print(electronic_dielectric)
  
  # Initialise arrays
  polr = np.zeros([n_atoms * 3, 3])
  pol = np.zeros([n_atoms*3,6])
  pol_jump=np.zeros([])
  polr_jump=np.zeros([])
  index=0
  for atom in np.arange(n_atoms):
    for coord in np.arange(3):
      for delta in deltas:
        folder=name+".i_atom_"+str(atom)+".i_coord_"+str(coord)+".displ_"+str(delta)
        volume, pol_jump, polr_jump, dielectric = read_aims_output(folder, aimsoutname, n_atoms)
        pol[index,:]=pol[index,:]+pol_jump*coeff[deltas==delta]*c_zero                  # grad polar using finite difference
        polr[index, :] = (polr[index, :] + polr_jump * coeff[deltas == delta] * c_zero) # C/m2  -> C/m3
      index=index+1

  # Calculation of Born charges
  volume = volume * Ang * Ang * Ang            # convert to SI units
  grad_polr= polr * volume / elementary_charge # convert charge units to electrons
  print("Grad POLR")
  print(grad_polr)

  with open(bornchargefile,"w") as out:
        d = electronic_dielectric
        print( "#epsilon and Z* PDielec/PDGui file format (FHI-Aims results)" ,file=out)
        print( f"{d[0,0]:15.9f} {d[0,1]:15.9f} {d[0,2]:15.9f} {d[1,0]:15.9f} {d[1,1]:15.9f} {d[1,2]:15.9f} {d[2,0]:15.9f} {d[2,1]:15.9f} {d[2,2]:15.9f}",file=out)
        j = 0
        d = grad_polr
        for _i in range(n_atoms):
            # print( f"{d[j+0,0]:15.9f} {d[j+0,1]:15.9f} {d[j+0,2]:15.9f}", end=' ',file=out)
            # print( f"{d[j+1,0]:15.9f} {d[j+1,1]:15.9f} {d[j+1,2]:15.9f}", end=' ',file=out)
            # print( f"{d[j+2,0]:15.9f} {d[j+2,1]:15.9f} {d[j+2,2]:15.9f}",file=out)
            print( f"{d[j+0,0]:15.9f} {d[j+1,0]:15.9f} {d[j+2,0]:15.9f}", end=" ",file=out)
            print( f"{d[j+0,1]:15.9f} {d[j+1,1]:15.9f} {d[j+2,1]:15.9f}", end=" ",file=out)
            print( f"{d[j+0,2]:15.9f} {d[j+1,2]:15.9f} {d[j+2,2]:15.9f}",file=out)
            j += 3
          
if __name__ == "__main__":
    main()
    print("\n Done. ")
