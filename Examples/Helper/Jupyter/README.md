# Jupyter Notebook Examples
## Bonds_Angles.ipynb
Transforms the unit cell aspartic acid to a cell which has whole molecules in it.
Calculates the bond lengths, bond angles and torsion angles.
## calculatePowderSpectrum.ipynb
Calculate the powder spectrum of a VASP calculation on ZnO.
Compares absorption for Maxwell-Garnett and Bruggeman methods at two volume fractions.
## crossCorrelation.ipynb
Calculate the cross correlation of a couple of very similar spectra.
The spectra are taken from Castep calculations on Na2(SO4)2.  
The 'experimental' calculation is one with PTFE as the matrix support and the calculated spectrum is using NUJOL.
## MaterialPoperties.ipynb
Illustrate the use of the Materials class using a helper function.
Sapphire is used as the example.
Information about Sapphire is printed.
The permittivity function is used to calculate the real and imaginary comonents of the permittivities.
## ModifyAtomicMasses.ipynb
Modify the calculated spectrum by changing the masses.
A Vasp calculation is used to illustrate how to alter the masses and calculate the spectrum.
## modifyDrudeLorentzPermittivity.ipynb
Modify the drude-Lorentz permittivity function by changing the sigmas and by scaling the frequencies.
## ReadDFTFiles.ipynb
Illustrates the use of readers to read DFT output files.
The permittivity of a Vasp calculation on ZnO is calculated using the permittivity function of the material
## singleCrystalSpectrum.ipynb
Calculate a single crystal spectrum of a multilayered system.  
This uses the calculateSingleCrystalSpectrum() helper function.
## SuperCell.ipynb
Illustrate the properties of the SuperCell class.
## UnitCell.ipynb
Illustrate the properties of the UnitCell class.
## compareDFTAbsorptions
Compares calculations on Na2(SO4)2 using different DFT methods.  
The results of the calculations are collected together as separate directories in Na2SO42/
