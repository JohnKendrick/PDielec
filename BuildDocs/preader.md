author: John Kendrick
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom.
eMail: j.kendrick@leeds.ac.uk
author: Andrew Burnett
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom
title: PREADER
email: a.d.burnett@leeds.ac.uk
Bibliography: ./pdielec.bib

[INCLUDE="style"]

<!-- Comment out some of the options -->
<!-- Csl Style: ieee -->
<!-- Math Mode: static -->
<!-- [INCLUDE="style"] -->
<!-- Colorizer: javascript -->
<!-- Doc class: [10pt]article -->

~ MathDefs
\newcommand{\water}{H_{2}O}
\newcommand{\tensor}[1]{\bar{\bar{#1}}}
\newcommand{\tensorbs}[1]{\bar{\bar{\boldsymbol{#1}}}}
\newcommand{\tensorbf}[1]{\bar{\bar{\mathbf{#1}}}}
\newcommand{\fieldbf}[1]{\bar{\mathbf{#1}}}
~

[TITLE]

[TOC]

# INTRODUCTION

preader is a 'helper' program which uses the underlying modules of PDielec to read output files and summarise the results of various MM/QM packages.  The program can be used to perform some straightforward calculations.  For instance projection of any remaining centre of mass motion of the crystal can be performed to make sure that there are three zero frequencies.  Also the masses used in the calculation of the dynamical matrix can be altered.

# MM/QM PACKAGES

The MM / QM packages supported are summarised below.  Unlike PDielec it is not necessary to have performed a full calculation of the dynamical matrix.  In the majority of cases preader will read geometry optimisation runs.

**VASP** -program vasp OUTCAR
The name provided on the command line is an OUTCAR file. The OUTCAR is read by preader to determine the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. 

**CASTEP** -program castep seedname
The name provided on the command line is the seedname for the calculation. The corresponding seedname.castep file in the current directory is read and processed to determine the unit cell, atomic masses, optical permittivity and born charge tensors. The normal modes and their frequencies are determined from the seedname.phonon file.

**CRYSTAL** -program crystal outputfilename
The name on the command line is a file ending in .out, containing the output of a CRYSTAL14 run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory where preader is run from , it uses these files to calculate the Born charge tensors, frequencies and normal modes. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the preader package. Small differences in the calculated frequencies between the CRYSTAL program and preader have been observed. These have been found to be due to a slightly different method for symmetrising the 2^nd^ derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that preader should use the same symmetrisation as CRYSTAL14.

**ABINIT** -program abinit outputfilename
The output file should come from a run containing three datasets. One to calculate the wavefunction at the optimized geometry, one to calculate the field perturbations and one to calculate the second derivatives. Examples of input files and output files are available with the distribution.

**QE** -program qe outputfilename
The output file is the dynamical matrix file, specified by "filedyn" in a run of the quantum espresso phonon package. Examples of input and output files are given in the preader distribution

**GULP** -program gulp outputfilename
The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit cell. The key words; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it on the command line (see -optical and -optical\_tensor options below).

# PROGRAM OPTIONS

Examples of data sets for these packages are included with the distribution and can be found in the Examples/'Package'/preader directory. The program is run from the command line. There are several command options and these are summarized below. Some of the options may be repeated.

~TableFigure {#tab-options; caption: "Command line options for preader <br>^1^A tick indicates an option can be used more than once "; breakable:true}
| Option                     |  Purpose                                  | R^1^ |
| :------------------------- | :---------------------------------------| :-------: |
| -program string            |  string can be “abinit”,  “castep”, “crystal”, “gulp”, “qe”, “experiment” or “vasp” and specifies the program which generated the results to be analysed |           |
| -neutral                    | Impose neutrality on the Born charge matrices |           |
| -nocalculation             | Requests that no calculations are performed.  This results in a single line of output with just information from the program.  If -eckart, -mass, -masses, -neutral or -crystal have -hessian crystal have been specified they will be ignored |           |
| -eckart                    | The translational modes are  projected out of the hessian before diagonalisation |           |
| -hessian string            | string may be either “crystal” or  “symm”.  In the case of “crystal” the  hessian is symmetrised using the same algorithm as Crystal14.  “symm” is the default |           |
| -masses string             | string can be either “program”, “average” or“isotopic”, meaning that the masses used in the calculation of the frequenciesare either taken from the QM program or are the average of the isotopeabundances or are the most abundant isotope mass. |           |
| -mass string real          | The atomic mass of string is set to real.  This can be used to explore the effect ofisotope substitution on the calculated frequencies |     ✔     |
~


# EXAMPLES

         preader -program vasp `find . -name OUTCAR` > results.csv

This reads all the VASP OUTCAR files in the current and any of its subdirectories and summarises the results to results.csv.

         preader -program castep -eckart `find . -name \*.castep` > results.csv

This reads all the castep output files in the current and any of its subdirectories and summarises the results to results.csv.  For each file the centre of mass motion of the crystal is projected.  The results file contains both the unprojected and the projected results.

There are examples of preader being used in the Examples/'Package'/preader subdirectories of the distribution of PDielec.



