# Phonopy with Quantum-Espresso

The source of the QE input and the source of the cell and coo-ordinates comes from urea.scf.in.
This file contains the optimised geometry, the pseudo-potentials and the various cutoffs used, but it is a single point calculations.

There were some changes needed though to make the phonopy program work as the phonopy qe reader cannot cope with all the options available to pwscf input.
The most significant change I found was that the units on my CELL_PARAMETERS line were in curly brackets and Phonopy did not accept those.  Simply removing the curly brackets should be sufficient.

1x1x1 supercells should be sufficient for infrared intensity calculations.
These were generated based on the geometry in uread.scf.in with;

```
phonopy --qe -d --dim="1 1 1" -c urea.scf.in
```

This generated a set of 16 files containing the displaced geometries, supercell-001.in to supercell-016.in.
These files do not have pseudo-potential, cutoffs etc and so it is necessary to create the final input files by concatenating a header, the displaced geometry, and a kpoints description.

The header.in is based on urea.scf.in and is shown below:
The celldm(1) statement urea.scf.in has been removed, as the cell dimensions given by Phonopy are in Bohr in the supercell file.  
Otherwise it can be created from urea.scf.in by removing everything after and including the CELL_PARAMETERS line.

```
&CONTROL
  calculation = 'scf' 
  restart_mode = 'from_scratch'
  outdir = '.' 
  pseudo_dir = '/home/john/Data/PseudoPotentials/SSSP_1.3.0_Efficiency'
  nstep = 1000
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = 1.0d-7
  forc_conv_thr = 1.0d-6
/
&SYSTEM
  ibrav       = 0
  nat         = 16
  ntyp        = 4
  occupations = 'fixed'
  ecutwfc     = 90.0
  ecutrho     = 760.0
  input_dft   = 'pbe'
  vdw_corr    = 'grimme-d3'
  dftd3_version = 6
  dftd3_threebody = .false.
/
&ELECTRONS
  conv_thr = 1.0D-10
  mixing_beta = 0.3
  diagonalization = 'rmm-paro'
/
&IONS
  ion_dynamics = 'bfgs'
  pot_extrapolation = 'second_order'
/
&CELL
  cell_dynamics = 'bfgs'
  cell_dofree   = 'all'
/
```

The kpoints.in file starts with an empty line and is show below.
The blank line at the beginning is necessary because the supercell files do not end with a carriage return.
The need to use a kpoints.in file is not mentioned in the Phonopy documentation.

```

K_POINTS automatic
5 5 6 0 0 0
```

The final inputs for each supercell calculation can now be created in urea-001.in to urea-016.in.

```
for f in {001..016}
do
 cat header.in supercell-$f.in kpoints.in > urea-$f.in
done
```

The supercell files can now be deleted if necessary and the jobs submitted to a batch queuing system.
Assuming that only one job is running at a time....

```
for f in {001..016}
do
  runpwscf urea-$f
done
```

After these have completed the FORCE_SETS can be calculate by Phonopy, followed by the dynamical matrix.

```
phonopy  --qe -f urea-{001..016}.log
phonopy --dim="1 1 1" --qpoints="0 0 0" --writedm
```

To calculate the Born charges it is necessary to run the QE phonon code (ph.x) after a single point calculation of the wavefunction.
The directory Born/ holds the single ploint calculation using urea.scf.in and the wavefunction is stored in pwscf.save/.   
In this directory there is an xml file (data-file-schema.xml), the charge density and the wavefunction.
The xml file is sometimes reproduced as pwscf.xml in the directory where the job was run.
The file holds copies of the input and output information for the job.
The Born charges can be calculated running ph.x together with the following input (see file urea.born.in):

```
Born charges
 &inputph
! Calculate dynamic polarisabilities
  fpol=.false.
! Calculate the dielectric constant
  epsil=.true.
! Calculate the Born charges using the forces as a derivative of the field
  zeu = .true.
! Do not calculate Raman spectra
  lraman=.false.
! Do not calculate electron phonon coupling
  elop=.false.
! Do not calculate dyanmical matrix
  trans=.false.
! A very stringent criteria for self consistency
  tr2_ph=1.0d-14
! Increase the number of iterations used for the potential mixing
  nmix_ph=8
! Set the output directory
  outdir='.'
! Calculate the Grimme D3 hessian
  dftd3_hess='automatic.hess'
! Acoustic sum rule
  asr=.false.
! Give the name of the dyanmic matrix file
  fildyn='urea.dynG'
/
0.0 0.0 0.0
```

This calculates the Born charges and the zero frequency optical permittivity at the gamma point using the combinattion of commands, zeu=.true., trans=.false., and epsil=.true..
The first (zeu=.true.) ensures the Born charges are calculated.
The second (trans=.false.) stops the calculation of the full dynamical matrix.
The last (epsil=.true.) causes the calculation of the optical permittivity.
The only difference between this and a full infrared intensity calculation is the use of (trans=.false.).
However, unfortunately this does not cause urea.dynG to be written as there is no dynamical matrix, so the information regarding the charges and the permittivity must be obtained elsewhere.
Whilst the information is available in the standard output file, the precision to which it is written makes it unsuitable to be used.  Even the alat parameter is not written to a sufficient number of decimal places.
Instead it is necessary to process two xml files, pwscf.xml and tensor.xml.
The first is a copy of that found in the pwscf.save directory, the second is a copy of the _ph0/_pwscf.phsave/tensor.xml file.  These two files should be in the same directory as the phonopy files.
Once copies of these files have been made the _ph and pwscf.save directories and there contents can be deleted.
In this example the files have been copied to the Born/ directory and links have been made in the directory where the Phonopy calculations were performed.

It should be noted that pdgui only recognises the names pwscf.xml and tensors.xml.
