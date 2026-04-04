# Examples/FHI-Aims/Na2SO42

## Na2(SO4)2 FHI-Aims calculation of powder IR

The FHI-Aims calculation of vibrational intensities were performed using the get_vibrations.py script supplied with the FHI-Aims.  FHI-Aims calculations were performed using the "tight" basis set.  Atomic positions and cell dimension were optimised. The get_vibrations.py script (from now on called get_vibrations) was used to generate displacements for the minimum and the displaced forces and displaced polarisations were calculated.

A typical protocol for generating the required data is given below.  In the directory where the geometry and cell optimisation has been performed create a directory /Vibrations and copy the optimised geometry and the control file into it.

```bash
mkdir Vibrations
cp geometry.in.next_step Vibrations/geometry.in
cp control.in Vibrations/control.in
cd Vibrations
```

The geometry needs to be in cartesian coordinates, otherwise the get_vibrations.py script fails.  There is a utilities script called `clims-convert-geometry`, which is part of the clims distribution, which can be used to perform the conversion.

```bash
pip install clims
# clims requires some additional configuration....
clims-convert-geometry --filein geometry.in --cartesian
mv geometry.in geometry.in.fractional
mv geometry-converted.in geometry.in
```

The control.in file which was used for the optimisation needs to be modified to allow single point calculations of the forces and the polarisation at the displace geometries, keeping the basis set specification.  The header of  typical control.in file would look like this:

```text
#===============================================================================
xc                                 pbe
d3
relativistic                       atomic_zora scalar
k_grid                             7 6 5
sc_accuracy_rho                    1E-6
charge_mix_param                   0.4
KS_method                          serial
final_forces_cleaned               .true.
compute_forces                     .true.
output                             polarization   1   20  5  5 
output                             polarization   2    5 20  5 
output                             polarization   3    5  5 20
#===============================================================================
```

The instructions in this file are followed by the basis set specification used by the optimisation step.  The instructions must use the same density functional and dispersion correction and the same k-point sampling.  The self consistent convergence threshold is strigent and the charge mixing parameter set to a value which improves convergence.  The KS_method is set to serial as the Berry Phase method for calculating the polarisation is not fully parallelized.  The calculation of the polarisation is performed on a k-point grid.  The grid needs to be converged, the numbers given here are illustrative.  It is common to use more k-points in the direction of the polarisation.

The perturbed geometries are generated using the Python script get_vibrations.py distributed as part of the FHI-Aims distribution in the utilities directory.  This script needs to be available for execution.  The script generates a directory for each perturbed geometry.

```
python get_vibrations.py na2so42 0
```

The first parameter in the command line is used as a prefix to the directory names containing the perturbed geometries.  The second (a zero) indicates that the geometries are to be generated and the control.in file copied in each directory.

The single point calculations can be submitted to the batch queuing system.  Assuming that there is a command, ``runaims`` which submits an FHI-Aims job to the batch queue and generates an output file ``aims.out``, the following will submit all the single point jobs to the queue.

```bash
for f in na2so42.i* do
  (cd $f; runaims)
done
```

While these calculations are running, the job to calculate the electronic contribution to the permittivity, $\epsilon_{\infty}$, can be setup and run.  Create a directory called Dielectric and copy in to it the control.in and geometry.in files.    The header on the control.in files needs to be modified to calculate the dielectric constant using DFPT.  An example is shown below.  Submit this job to the job queue.

```textile
#===============================================================================
xc                                 pbe
d3
relativistic                       atomic_zora scalar
k_grid                             7 6 5
sc_accuracy_rho                    1E-6
charge_mix_param                   0.4
DFPT                               dielectric
#===============================================================================
```

Once the all the calculations have finished succesfully, get_vibrations.py can be used to calculate the intensities and frequencies of the vibrational modes.  

```bash
python get_vibrations.py --IR na2so42 2
```

This generates several files including; hessian.na2so42.dat, masses.na2so42.dat and na2so42.ir.  The first two are numpy save files for the hessian and masses.  The last is a list of frequencies and intensities for the phonon modes.

The PDielec script aims-pdielec-born can now be used to generate the Born charges which will be stored in BORN_PDIELEC.  The aims-pdielec-born command needs a Dielectric/ directory as well as the displacement directories.  In this example only the Dielectric/ directory is included.  

For the PDGui calculations, a 10% volume fraction of Na2(SO4)2 is supported in a ptfe matrix.  The Maxwell-Garnett effective medium method is used with powdered spheres.  A comparison is given of the results from VASP and from FHI-Aims in Vasp_FHI-Aims_Comparison.pdf

| **File**                     | **Description**                                           |
| ---------------------------- | --------------------------------------------------------- |
| results.ref.xlsx             | Reference results for the test suite                      |
| script.py                    | script file using for testing                             |
| command.pdmake               | pdmake command file for testing                           |
| control.in                   | AIMS control file for vibrational calculation             |
| geometry.in                  | AIMS geometry for vibrational calculation                 |
| BORN_PDIELEC                 | Born charges calculated using aims-pdielec-born script    |
| hessian.na2so42.dat          | Hessian file from the get_vibrations script               |
| masses.na2so42.dat           | Masses file from the get_vibrations script                |
| na2so42.ir                   | IR frequencies and intensities from get_vibrations        |
| na2so42.xyz                  | Normal modes in cartesian coordinates from get_vibrations |
| Vasp_FHI-Aims_Comparison.qti | qtiplot file for comparing Vasp and FHI-AIMS              |
| Vasp_FHI-Aims_Comparison.pdf | pdf file for comparing Vasp and FHI-AIMS                  |
| Dielectric/control.in        | The AIMS control file for the calculation epsilon inf     |
| Dielectric/geometry.in       | The AIMS geometry file for the calculation epsilon inf    |
| Dielectric/aims.out          | The AIMS output containing epsilon infinity               |
| Vasp_FHI-Aims_Comparison.qti |                                                           |

[Back](../README.md)
