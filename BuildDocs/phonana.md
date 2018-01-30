author: John Kendrick
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom.
eMail: j.kendrick@leeds.ac.uk
author: Andrew Burnett
address: School of Chemistry, University of Leeds, Leeds, LS2 9JT United Kingdom
title: PHONANA
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
\newcommand{\tensorbs}[1]{\bar{\bar{\bm{#1}}}}
\newcommand{\tensorbf}[1]{\bar{\bar{\bm{#1}}}}
\newcommand{\fieldbf}[1]{\bar{\bm{#1}}}
~

[TITLE]

[TOC]


# INTRODUCTION

Phonana is a 'helper' program which uses the underlying modules of PDielec to read an output file from an MM or QM calculation of the phonon modes in a crystal.  The phonon modes are analysed and the percentage of molecular centre-of-mass and rigid body rotational motion in each mode is calculated.  Using the results of the analysis it is possible to differentiate between inter- and intra-molecular phonon modes (see for instance Jepsen et al [@Jepsen2007]).

# THEORY

## Molecular Systems

It is common in molecular calculations of
vibrational properties to construct a force constant matrix which enforces the
requirements of zero energy change for centre-of-mass motion and rigid-body rotation. This can be achieved by using projection operators to transform
the second derivative matrix to a set of coordinates which no longer include centre-of-mass motion or rigid-body rotation.

Defining a projection operator as;

~ Equation {#eq-projectionop}
\tensorbf{P} = \tensorbf{1} - \fieldbf{V} \fieldbf{V}^{T}
~

where $\fieldbf{V}$ is an orthonormal column vector with length $3N$ ($N$ is the number of atoms), the centre-of-mass motion
can be described in Cartesian space as all atoms moving along the x, y or
z axis with the same displacement. So considering the projection of the centre-of-mass motion along the x-axis, for each atom $a$ we can write;


~ Equation {#eq-cmx1}
\fieldbf{V}^{x}_{a} = \left(\begin{matrix} 1 \\ 0 \\ 0 \end{matrix}\right) 
~

For the projection operator which will project out all components of translation along the x-axis for every atom in the molecule we have;

~ Equation {#eq-cmx2}
\fieldbf{V}^{x} = \begin{pmatrix} \fieldbf{V}^{x}_{1} \\ \fieldbf{V}^{x}_{2}\\ \vdots \end{pmatrix} 
~

There are two other projection operators describing translation along the y- and z-axis.  

In a similar fashion it is possible to describe an infinitesimal molecular rotation using a vector V and therefore constructing a projection operator to remove rigid-body rotation.  If the coordinates of atom $a$ relative to the centre-of-mass of the molecule are $x, y \text{ and } z$, the component of the projection vector, $\tensorbf{V}$ representing rotation about the x-axis in the yz-plane is;


~ Equation {#eq-rotyz}
\fieldbf{V}^{yz}_{a} = \left(\begin{matrix} 0 \\ -z \\ +y \end{matrix}\right) 
~

In a similar fashion rotations about the y- and z-axis are respectively;


~ Equation {#eq-rotxz}
\fieldbf{V}^{xz}_{a} = \left(\begin{matrix} +z \\ 0 \\ -x \end{matrix}\right)
~

and

~ Equation {#eq-rotxy}
\fieldbf{V}^{xy}_{a} = \left(\begin{matrix} -y \\ +x \\ 0 \end{matrix}\right)
~

In a similar fashion to Equation [#eq-cmx2] the complete projection operators for the 3 rotational modes can be assembled (using rotation about the x-axis as an example)

~ Equation {#eq-rotyz2}
\fieldbf{V}^{yz} = \begin{pmatrix} \fieldbf{V}^{yz}_{1} \\ \fieldbf{V}^{yz}_{2}\\ \vdots \end{pmatrix} 
~

In practice the projection operators will be defined using mass-weighted cartesian coordinates as this simplifies the expressions used later.  As an example the transformation to mass-weighted coordinates for atom $a$ translating along the x-axis is;


~ Equation {#eq-mass-weight}
\sqrt{m_a} \fieldbf{V}^{x}_{a}
~

  
## Periodic Systems

The invariants of the energy in a periodic system are the three
translational modes corresponding to motion of all the atoms in the same
direction.  The rotation of an infinite lattice is not invariant.  However,  it is proposed to construct projection operators for each
molecule in the unit cell, which will project out the motion of
each *molecular* translation and rotation. Such projection operators will
be used to separate the external from the internal modes and there will therefore
be a projection operator for each external mode in the system.

Hug and Haesler have shown [@Hug2005c] that the vibrational kinetic energy can be decomposed into single centre atomic contributions by considering the kinetic energy, $T_p$, of a normal mode $Q_p$;


~ Equation {#eq-kep}
T_p = \frac{1}{2} \dot{Q_p} = 
\sum\limits_{a}{T_{a,p}} = 
\frac{1}{2}\sum\limits_{a}{m_a \fieldbf{\Delta\dot{x}}^T_{a,p}\fieldbf{\Delta \dot{x}}_{a,p} } =  
\frac{1}{2} \dot{Q^2_p} \sum\limits_{a} { \fieldbf{L}^T_{a,p} \fieldbf{L}_{a,p} }
~

Where the sum is over atoms, $a$, with mass $m_a$, $\fieldbf{\Delta \dot{x}}_{a,p}$ is the time
derivative of the Cartesian displacement vector of atom $a$ along the $p^{th}$ normal mode and $\fieldbf{L}_p$ are the
components of the orthogonal transformation relating the $p^{th}$ normal
mode and mass-weighted Cartesian displacements. Thus $\fieldbf{L}$ are the eigenvectors of the 
mass-weighted second derivative matrix.  The atom subscript is used to indicate
only those atomic coordinates involving atom $a$ are being considered.
Derivation of this expression has made use of the following
relationship between the mass weighted normal mode;

~  Equation {#eq-kinetic-energy}
\sqrt{ m_a } \fieldbf{\Delta x}_{a,p} = \fieldbf{L}_{a,p} Q_p
~

As proposed by Hug and Haesler [@Hug2005c], since the kinetic energy can be expressed as
atomic contributions, consideration of the virial theorem indicates that
this is true also for the potential energy. Using mass weighted
Cartesian coordinates for convenience, an analysis of the phonon modes
may therefore be constructed in the following way. A given phonon
mode, $p$,  will have its total kinetic energy partitioned between the
molecules in the unit cell according to;

~  Equation {#eq-kinetic-energy_total}
E^{total}_p = \sum\limits_{mol} {E^{mol}_p } 
~
~  Equation {#eq-kinetic-energy_molecular}
E^{mol}_p = \dot{Q^2_p} \sum\limits_{a \in mol} { \fieldbf{L}^T_{a,p} \fieldbf{L}_{a,p} }
~

To calculate the contribution a particular mode has to the centre-of-mass kinetic energy we use the projection operators given in Equation [#eq-cmx2] (now in mass-weighted Cartesian coordinates) where only the atoms in that molecule are used to construct the operator;

~ Equation {#eq-finaloperator}
\tensorbf{P}^x_{mol} = \tensorbf{1} - \fieldbf{V}^x_{mol} (\fieldbf{V}^x_{mol})^T
~

This approach is used to define  6 projection operators for each molecule. The centre-of-mass energy contained in the $p_{th}$ mode can be written;

~ Equation {#eq-energy-cm}
E^{cm}_p = \dot{Q^2_p} \sum\limits_{i=x,y,z}{\sum\limits_{mol} { (\tensorbf{P}_{mol}^i\fieldbf{L}_{p})^T (\tensorbf{P}_{mol}^i \fieldbf{L}_{p} ) }}
~

and the rigid-body rotational energy associated with the mode is;

~ Equation {#eq-energy-rot}
E^{rot}_p = \dot{Q^2_p} \sum\limits_{i=xy,yz,xz}{\sum\limits_{mol} { (\tensorbf{P}_{mol}^i\fieldbf{L}_{p})^T (\tensorbf{P}_{mol}^i \fieldbf{L}_{p} ) }}
~

The molecular vibrational contribution, which can be used to classify the internal modes of the system, can be obtained by subtracting the centre-of-mass and rigid-body rotational; energies from the total.

~ Equation {#eq-energy-vib}
E^{vib}_p = E^{total}_p - E^{cm}_p - E^{rot}_p
~

Since all of the energy terms depend in the same way on $\dot{Q}_p$ , it is not required in the calculation of the *relative* contributions to the energy coming from the external (molecular centre-of-mass and rigid-body rotation) and the internal modes (vibrational contributions)

# IMPLEMENTATION

In order to calculate the relative internal and external contributions to each phonon mode it is first necessary to identify the molecules in the crystal.  In many cases the unit cell is packed with atoms in such a way that the cell is filled, rather than in a way reflecting the bondedness of the molecules.  Phonana therefore first replicates the atoms in all cells neighbouring the central unit cell.  Within this supercell the bonds are determined by calculating the distances between all atoms in the supercell.  In practice an order N method is used whereby instead of searching all the supercell only the space around each atom is searched for potentially bonded partners.  The criterion of the presence of a bond between atoms $i$ and $j$ is given by the requirement that the distance between the atoms $r{ij}$ is less than the bonding requirement;


~ Equation {#eq-bond}
r_{ij} < scale (radius_i+radius_j)+toler
~

Here $radius_i$ is the covalent radius of atom $i$, $scale$ and $toler$ are factors which can be altered when running the program.

Once the bonding in the supercell has been determined the program starts with the first atom in the central cell and determines all atoms which are connected to the molecule that it is in.  If there are any remaining atoms in the central cell which are not bonded yet then further molecules are added until all atoms in the central cell have been allocated to a molecule.

Finally a new cell with the same dimensions as the original is constructed.  Where necessary each molecule is shifted into the cell so that the molecules centre-of-mass lies within the it.

The projection operators, Equation [#eq-finaloperator], can now be constructed and the relative energies calculated using Equations [#eq-energy-cm], [#eq-energy-rot] and [#eq-energy-vib]. 

## QM/MM PROGRAMS

The MM / QM packages supported are summarised below.

* VASP (-program vasp OUTCAR)
  : The name provided on the command line is an OUTCAR file. The OUTCAR is read by phonana to determine the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. The VASP run can be a DFPT or numerical calculation of the response.

* CASTEP (-program castep seedname)
  : The name provided on the command line is the seedname for the calculation. The corresponding seedname.castep file in the current directory is read and processed to determine the unit cell, atomic masses, optical permittivity and born charge tensors. The normal modes and their frequencies are determined from the seedname.phonon file.  The CASTEP run needs to be a DFPT (phonon+efield) task.

* CRYSTAL (-program crystal outputfilename)
  : The name on the command line is a file ending in .out, containing the output of a CRYSTAL14 run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes and Born charge tensors. However, the number of significant figures for the normal modes is not sufficient for an accurate calculation and it is therefore recommended that the HESSFREQ.DAT and BORN.DAT files are also made available. If they are present in the directory where phonana is run from , it uses these files to calculate the Born charge tensors, frequencies and normal modes. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the phonana package. The CRYSTAL calculation needs to be a frequency calculation (FREQCALC) with the infrared intensity (INTENS) selected. The default algorithm does not calculate the optical permittivity, so this needs to be provided on the command line. However, if the CPHF or CPKS algorithm is used for the frequency calculation, the optical permittivity is calculated and PDielec will automatically read it from the output file. By default CRYSTAL projects out the pure translational modes of the system before calculating the frequencies, this can also done by the PDielec package.Small differences in the calculated frequencies between the CRYSTAL program and phonana have been observed. These have been found to be due to a slightly different method for symmetrising the 2^nd^ derivative matrix, because of this an optional directive "-hessian crystal" can be used to indicate that phonana should use the same symmetrisation as CRYSTAL14.

* ABINIT (-program abinit outputfilename)
  : The output file should come from a run containing three datasets. One to calculate the wavefunction at the optimized geometry, one to calculate the field perturbations and one to calculate the second derivatives. Examples of input files and output files are available with the distribution.

* QE (-program qe outputfilename)
  : The output file is the dynamical matrix file, specified by "filedyn" in a run of the quantum espresso phonon package. Examples of input and output files are given in the phonana distribution

* GULP (-program gulp outputfilename)
  : The name on the command line is a file ending in .gout, containing the output of a GULP run. The contents of this file alone are sufficient to provide the unit cell, atomic masses, frequencies, normal modes, Born charge tensors and optical permittivity. Because GULP only writes out the Born charge matrices for the asymmetric unit, it is necessary to run a frequency calculation using P1 symmetry and a complete unit cell. The key words; nosymm, phonon, intensity, eigen and cart are recommended for the GULP calculation. In the case that no shells are used in the calculation the optical permittivity is not available in the output and it is necessary to provide it on the command line (see -optical and -optical\_tensor options below).

## COMMAND OPTIONS

There are several command options and these are summarized in below. Some of the options may be repeated and these are indicated by a ✔.  Where there is a default its value is show.

* -program program
  : program can be “abinit”,  “castep”, “crystal”, “gulp”, “qe”, “experiment” or “vasp” and specifies the program which generated the results to be analysed.
* -vmin vmin (0.0)
  : vmin is the starting wavenumber (cm^-1^)  for the frequency range.
* -vmax vmax (9000.0)
  : vmax is the final wavenumber (cm^-1^)  for the frequency range
* -radius element value (✔)
  : The covalent radius of 'element' is set to value
* -toler toler (0.1)
  : Set the tolerance value used for calculating bonds. A bond is created between two atoms whose distance apart is less than scale*(radi+radj)+toler, where radii is the covalent radius of atom i
* -scale scale (1.1)
  : Set the scale factor for calculating bonds
* -excel filename
  : Write the results to an excel  spread sheet with the name  specified.
* -csv filename
  : Output is sent to a comma separated  file specified.
* -print
  : Additional output is provided from  the QM or MM calculation
* -ignore k (✔)
  : Ignore the kth mode
* -mode k (✔)
  : Only use the kth mode in the analysis
* -eckart
  : The translational modes are  projected out of the hessian before diagonalisation
* -hessian string (symm)
  : string may be either “crystal” or  “symm”.  In the case of “crystal” the  hessian is symmetrised using the same algorithm as Crystal14
* -masses mass_definition (average)
  : mass_definition can be either “program”, “average” or “isotopic”, meaning that the masses used in the calculation of the frequenciesare either taken from the QM program or are the average of the isotopeabundances or are the most abundant isotope mass.
* -mass element mass (✔)
  : The atomic element set to the value mass.  This can be used to explore the effect of isotope substitution on the calculated frequencies
* -viewer
  : Opens a graphics window showing the molecular structure of the unit cell and an the displacement of the phonon modes as arrows


## INSTALLATION
Phonana doesn't need any extra installation above what is needed to run PDielec, except if the the -viewer option is used.  The viewer is a Mayavi application and needs the following installed.

    sudo pip install vtk
    sudo pip install pyvtk
    sudo pip install  vtkinterface
    sudo pip install  mayavi

The installation of these modules can be performed as root, if all users on the machine require access to them, or they can be installed in the user's file system using 'pip --user' instead of 'pip'

On some Linux systems they installation of VTK can cause problems.  The most common of which can be an error message indicating that the vtkOpenGLKitPython module cannot be found.

This can be fixed by adding the vtk/ directory from your Python site-packages directory to your LD_LIBRARY_PATH.  If your LD_LIBRARY_PATH is otherwise empty and (for this example) the VTK pacakge was installed for all users;

     export LD_LIBRARY_PATH=/usr/lib/python3.6/site-packages/vtk

# EXAMPLES


## Usage
         phonana -program vasp OUTCAR 

This reads the results of the VASP calculation and prints a summary of the molecular centre-of-mass and molecular rotation contributions to each phonon mode.

         phonana -program vasp OUTCAR -excel results.xlsx

As above but with the results written to a spreadsheet

         phonana -program vasp OUTCAR -excel results.xlsx -eckart -radius S 1.0

As above but with the eckardt projection of the crystal centre-of-mass motion.  The covalent radius of the sulphur atom is set to 1 Angstrom.

## Output

### Isoleucine

The example of a calculation on isoleucine using Castep is available in the Examples/Castep/Isoleucine directory of the distribution.  The analysis was performed using;

        phonana -program castep phonon.castep -excel results.xlsx

The program first finds the four molecules of isoleucine in the unit cell of the crystal structure.  After summarising the the unit cell information, the elements of the atoms are given.  The ordering of the atoms at this point is different to the ordering in the original cell and in the QM calculation.  The program lists the cartesian and fractional coordinates.  Again these will not be the same as in the original calculation.  For each molecule the program lists the index numbers of the atoms in the molecule, the molecular mass and the centres of mass of the molecules. The results are summarised in Table [#tab-isoleucine-mols].

~TableFigure {#tab-isoleucine-mols; caption: "Molecular constituents of the isoleucine unit cell (centre-of-masses are in fractional coordinates)"; page-align:forcehere}

| Molecule | Mass    | CM~a~     | CM~b~     | CM~c~     |
|----------|---------|----------|----------|----------|
|0         | 131.172 | 0.293600 | 0.051683 | 0.304906 |
|1         | 131.172 | 0.706400 | 0.551683 | 0.695094 |
|2         | 131.172 | 0.793231 | 0.588919 | 0.287125 |
|3         | 131.172 | 0.206769 | 0.088919 | 0.712875 |
~

Table [#tab-isoleucine-results] shows that results of the analyis of the phonon modes with frequencies below 100 cm^-1^.  The first 3 modes should have zero frequency as they are the translationally invariant modes of the lattice.  The negative number shown here actually indicates that the mode has an imaginary frequency.  As can be see under the column %col-cme they are almost completely associated with molecular centre-of-mass motion.  This is true also of the lowest non-zero frequency at 32,81 cm^-1^ and to some extent of the next mode at 39.73 cm^-1^.  Although for this latter mode there is considerable contribution from rigid-body rotational motion (see the %mol-rot column).  As the frequency of the phonon mode increases the contribution from rigid-body motion generally decreases and contribution from vibrational mode (see the %vib column) increases,  The total contribution from molecular motion is summarised in the last four columns. 

~ TableFigure {#tab-isoleucine-results; caption: "Isoleucine: percentage contribtions of the centre-of-mass and the rotational molecular modes to each phonon mode"; page-align:forcehere}

| Freq(cm-1) | %mol-cme | %mol-rot | %vib | %mol-0 | %mol-1 | %mol-2 | %mol-3 |
|------------|----------|----------|------|-------|-------|-------|-------|
| -0.07      | 99.5     | 0.5      | 0.0  | 26.9  | 26.9  | 23.1  | 23.1  |
| -0.06      | 98.6     | 1.1      | 0.3  | 23.6  | 23.6  | 26.4  | 26.4  |
| -0.04      | 99.7     | 0.2      | 0.1  | 23.8  | 23.8  | 26.2  | 26.2  |
| 32.81      | 92.0     | 5.9      | 2.1  | 26.8  | 26.8  | 23.2  | 23.2  |
| 39.73      | 77.9     | 19.6     | 2.5  | 22.1  | 22.1  | 27.9  | 27.9  |
| 49.42      | 10.4     | 81.9     | 7.7  | 19.0  | 19.0  | 31.0  | 31.0  |
| 52.12      | 80.7     | 17.3     | 2.0  | 17.6  | 17.6  | 32.4  | 32.4  |
| 54.15      | 73.4     | 19.0     | 7.5  | 21.2  | 21.2  | 28.8  | 28.8  |
| 59.59      | 37.7     | 56.6     | 5.7  | 23.9  | 23.9  | 26.1  | 26.1  |
| 62.05      | 77.3     | 17.2     | 5.5  | 23.4  | 23.4  | 26.6  | 26.6  |
| 69.91      | 9.6      | 69.4     | 21.0 | 17.6  | 17.6  | 32.4  | 32.4  |
| 73.19      | 7.3      | 61.1     | 31.6 | 6.6   | 6.6   | 43.4  | 43.4  |
| 76.00      | 13.6     | 61.7     | 24.7 | 31.6  | 31.6  | 18.4  | 18.4  |
| 76.25      | 3.2      | 48.1     | 48.7 | 10.6  | 10.6  | 39.4  | 39.4  |
| 80.49      | 31.5     | 43.8     | 24.7 | 19.3  | 19.3  | 30.7  | 30.7  |
| 87.54      | 37.9     | 35.4     | 26.7 | 26.8  | 26.8  | 23.2  | 23.2  |
| 89.31      | 12.2     | 58.1     | 29.6 | 27.5  | 27.5  | 22.5  | 22.5  |
| 90.49      | 16.9     | 53.9     | 29.2 | 22.4  | 22.4  | 27.6  | 27.6  |
| 94.51      | 49.1     | 36.9     | 14.0 | 39.4  | 39.4  | 10.6  | 10.6  |
| 94.79      | 47.0     | 16.9     | 36.1 | 23.0  | 23.0  | 27.0  | 27.0  |

~
 ￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅￅ
### BaTiO~3~

The example for BaTiO~3~ can be found in Example/AbInit/BaTiO3-phonana from the main directory of the PDielec distribution. The results of the Abinit calculation can be analysed using the following command;

        phonana -program abinit BaTiO3.out -radius Ba 0.1 -excel results.xlsx
        
The -radius option has been used to set the covalent radius of Ba to a small value so this atom will be treated as though it is not bonded to anything else in the cell.  This results in 2 'molecules' being found in the cell; a TiO~3~ moiety and the Ba^2+^ ion. 

The results of the analysis are shown in Table [#tab-batio3-results] below.  In this example the lowest 3 modes (which have not been projected) show that there is a problem with translational invariance of the calculation.  The frequencies should be 0.0.  The fact that they are mainly centre-of-mass modes is shown by the high percentage in the %mol-cme column.  The energy in these three modes seems to be mainly in molecule 1, which is the Ba^2+^ ion.  The modes at 198.01 cm^-1^ also have a large centre-of-mass component, this time mainly coming from molecule 0, which is the TiO~3~ moeity. Above 200 cm^-1^ there is little centre-of-mass contribution to the energy and all the modes are dominated by the TiO~3~ group.

~ TableFigure {#tab-batio3-results; caption: "BaTiO~3~: percentage contribtions of the centre-of-mass and the rotational molecular modes to each phonon mode"; page-align:forcehere}

| Freq(cm-1) | %mol-cme | %mol-rot | %vib | %mol-0 | %mol-1 |
|------------|----------|----------|------|--------|--------|
| 65.74      | 95.9     | 1.0      | 3.1  | 26.1   | 73.9   |
| 65.74      | 95.9     | 1.2      | 2.9  | 26.1   | 73.9   |
| 65.75      | 95.9     | 1.1      | 3.0  | 26.1   | 73.9   |
| 198.01     | 70.2     | 6.6      | 23.2 | 75.3   | 24.7   |
| 198.01     | 70.2     | 8.3      | 21.5 | 75.3   | 24.7   |
| 198.02     | 70.2     | 7.4      | 22.3 | 75.3   | 24.7   |
| 270.59     | 0.0      | 55.6     | 44.4 | 100.0  | 0.0    |
| 270.59     | 0.0      | 34.4     | 65.6 | 100.0  | 0.0    |
| 270.59     | 0.0      | 76.7     | 23.3 | 100.0  | 0.0    |
| 279.39     | 33.8     | 18.5     | 47.7 | 98.6   | 1.4    |
| 279.39     | 33.8     | 19.3     | 47.0 | 98.6   | 1.4    |
| 279.39     | 33.8     | 17.8     | 48.4 | 98.6   | 1.4    |
| 493.63     | 0.1      | 17.4     | 82.5 | 100.0  | 0.0    |
| 493.63     | 0.1      | 17.4     | 82.5 | 100.0  | 0.0    |
| 493.63     | 0.1      | 17.4     | 82.5 | 100.0  | 0.0    |
~

### Isoleucine Visualisation
Using the same example as used in [#sec-isoleucine] the command;

    phonana -program castep phonon castep -viewer
    
was used to raise a basic 3D molecular viewer showing the molecular structure in the unit-cell and the displacements corresponding to a phonon mode as shown in Figure [#fig-viewer].

~ Figure { #fig-viewer; caption: "The Phonan 3D Viewer"; page-align:here }
![img-viewer]
~

[img-viewer]: Figures/phonana_viewer.png { width:90%; }

The graphical window responds to mouse presses in the usual way.  Rotation can be carried out with the left mouse button .  The right or wheel mouse buttons control the zoom.  The middle button controls the panning.  Below the graphical window is a mode selector, which indicates the mode selected and the range of modes which can be selected using the slider. The frequency of the selected mode is indicated below the slider.  On changing the selected mode the displacements (shown as arrows) also change.

There is an animate button which will start an animation of the selected mode.  The animations take some time to start (and finish).  Once started there is new small window which pops up, this can be ignored.  It is best to start and stop the animation using the 'Animate' button on the main window.  Another press of the 'Animate' button will stop the animation and will show the original displacement vector view.  The time taken for the change of view is substantial and the user needs to be patient.

[BIB]

