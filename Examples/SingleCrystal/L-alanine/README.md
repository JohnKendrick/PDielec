In an experimental and computational study of the the vibrational modes of single crystals of l-alanine the experimental transmittance was reported for polarized radiation along each of the principal axes, a-, b- and c-. (See reference below for the published paper)

This directory contains the results of Crystal17 calculations similar to those reported in the paper below. Details of the calculation can be found in phonon.dl2.  The output can be found in phonon.out.  Other files produced by the Crystal17 calculation are BORN.DAT and HESSFREQ.DAT.  These are used by PDGui to calculate the theoretical spectrum.

The experimental results for the measurement of terahertz transmittance were kindly provided by the authors of the experimental paper; J. L. Allen, T. J. Sanders, J. Horvat, R. A. Lewis and K. C. Rule.
These measurements can be found in the Excel spreadsheets a-axis-experimental.xlsx, b-axis-experimental and c-axis-experimental.xlsx.   Each of these spreadsheets has two columns, the first is a column of frequencies in cm-1 and the second is the transmittance.

The scripts used to generate the crystal transmittance from the calculated phonon modes can be found in a-axis.py, b-axis.py and c-axis.py.  These scripts use the Fitter Tab to compare the calculated spectrum with the experimental one.  Changes have been made to the film thickness and to the widths of the peaks to improve the agreement between the two.

To run the scripts use a command such as;

```
pdgui -script a-axis.py
```

We would like to thank the authors of the paper below for providing the experimental data and for helpful discussions.


JL Allen, TJ Sanders, J Horvat, RA Lewis, KC Rule;  
Physical Review Letters Vol130, 226901 (2023)
