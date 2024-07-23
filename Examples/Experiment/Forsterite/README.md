# Forsterite FPSQ model

The FPSQ model is  taken from Pierre, M. D. La, De La Pierre, M., Carteret, C., Orlando, R., & Dovesi, R. (2013). Use of ab initio methods for the interpretation of the experimental IR reflectance spectra of crystalline compounds. Journal of Computational Chemistry, 34(2), 1476–1485. https://doi.org/10.1002/jcc.23283

Single crystal calculations are performed for Forsterite films, which can be compared to the results in the above paper.  The Transfer Matrix method is used and a 'Thick Film' model us used for the Forsterite film.

| **File**              | **Description**                                            |
| --------------------- | ---------------------------------------------------------- |
| 9000534.cif           | Forsterite cif file from the Crystallography Open Database |
| command.pdmake        | pdmake script for the testing suite |
| forsterite.exp        | Experimental file format for FPSQ model of Forsterite |
| POSCAR                | VASP file use to generate the cell information in forsterite.exp |
| results.ref.xlsx      | Reference results for the test suite |
| script.py             | Script to generate the reference results |
