# L-Analine Single Crystal Calculations (Crystal DFT)

In an experimental and computational study of the the vibrational modes of single crystals of l-alanine the experimental transmittance was reported for polarized radiation along each of the principal axes, a-, b- and c-. (See reference below for the published paper)

This directory contains the results of Crystal17 calculations similar to those reported in the paper below. Details of the calculation can be found in phonon.dl2.  The output can be found in phonon.out.  Other files produced by the Crystal17 calculation are BORN.DAT and HESSFREQ.DAT.  These are used by PDGui to calculate the theoretical spectrum.

The experimental results for the measurement of terahertz transmittance were kindly provided by the authors of the experimental paper; J. L. Allen, T. J. Sanders, J. Horvat, R. A. Lewis and K. C. Rule.
These measurements can be found in the Excel spreadsheets a-axis-experimental.xlsx, b-axis-experimental and c-axis-experimental.xlsx.   Each of these spreadsheets has two columns, the first is a column of frequencies in cm-1 and the second is the transmittance.

The scripts used to generate the crystal transmittance from the calculated phonon modes can be found in a-axis.py, b-axis.py and c-axis.py.  These scripts use the Fitter Tab to compare the calculated spectrum with the experimental one.  Changes have been made to the film thickness and to the widths of the peaks to improve the agreement between the two.
The a-axis example is particularly sensitive to numerical problems and may not run succesfully on Windows machines.  Where possible, the Linux operating system should be used to run these examples.  If this is not feasible then the thickness of the l-alanine film should be reduced.

To run the scripts use a command such as;

```
pdgui -script a-axis.py
```

There are also some scripts for comparing the scattering and transfer matrix methods.  These are; transfer_matrix.py, scattering_vs_transfer.py and scattering_matrix.py.

We would like to thank the authors of the paper below for providing the experimental data and for helpful discussions.


JL Allen, TJ Sanders, J Horvat, RA Lewis, KC Rule;  
Physical Review Letters Vol130, 226901 (2023)


| File                  | Description                                       |
| --------------------- | ------------------------------------------------- |
| a-axis-calculated.xlsx | Calculated results for a-axis field |
| a-axis-experimental.xlsx | Experimental results for a-axis field |
| a-axis.py | PDGui script for a-axis field |
| b-axis-calculated.xlsx | Calculated results for b-axis field |
| b-axis-experimental.xlsx | Experimental results for b-axis field |
| b-axis.py | PDGui script for b-axis field |
| c-axis-calculated.xlsx | Calculated results for c-axis field |
| c-axis-experimental.xlsx | Experimental results for c-axis field |
| c-axis.py | PDGui script for c-axis field |
| BORN.DAT | Crystal Born file |
| HESSFREQ.DAT | Crystal Hessian file |
| command.pdmake | pdmake command for test suite |
| optimised.cif | A cif file of the optimised structure |
| phonon.d12 | Crystal input for phonon calculation |
| phonon.out | Crystal output for phonon calculation |
| results.ref.xlsx | Reference results for test suite |
| scattering_matrix.py | Calculations using Scattering Matrix theory |
| scattering_vs_transfer.py | Script for comparison of Transfer and Scattering matrix theories |
| script.py | Script to create reference results |
| transfer_matrix.py | Calculations using Transfer Matrix theory |
| unpol-experimental.xlsx | Unpolarised experimental results |
