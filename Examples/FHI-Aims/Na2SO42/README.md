# Examples/FHI-Aims/Na2SO42

## Na2(SO4)2 FHI-Aims calculation of powder IR

The FHI-Aims calculation of vibrational intensities were performed using the get_vibrations.py script supplied with the FHI-Aims.  FHI-Aims calculations were performed using the "tight" basis set.  Atomic positions and cell dimension were optimised. The get_vibrations.py script (from now on called get_vibrations) was used to generate displacements for the minimum and the displaced forces and displaced polarisations were calculated.


Once the calculations had finished, get_vibrations was used to calculate the intensities and frequencies of the vibrational modes.  The PDielec script aims-pdielec-born was used to generate the Born charges which are stored in BORN_PDIELEC.  The aims-pdielec-born command needs a Dielectric/ directory as well as the displacement directories.  However, in this example only the Dielectric/ directory is included.  The calculation of $\epsilon_{\infin}$ is performed in the Dielectric/ directory.

For the PDielec calculations, a 10% volume fraction of Na2(SO4)2 is supported in a ptfe matrix.  The Maxwell-Garnett effective medium method is used with powdered spheres.  A comparison is given of the results from VASP and from FHI-Aims in Vasp_FHI-Aims_Comparison.pdf

| **File**                     | **Description**                                        |
| ---------------------------- | ------------------------------------------------------ |
| results.ref.xlsx             | Reference results for the test suite                   |
| script.py                    | script file using for testing                          |
| command.pdmake               | pdmake command file for testing                        |
| control.in                   | AIMS control file for vibrational calculation          |
| geometry.in                  | AIMS geometry for vibrational calculation              |
| BORN_PDIELEC                 | Born charges calculated using aims-pdielec-born script |
| hessian.na2so42.dat          | Hessian file from the get_vibrations script            |
| masses.na2so42.dat           | Masses file from the get_vibrations script             |
| na2so42.ir                   | IR frequencies and intensities from get_vibrations     |
| Vasp_FHI-Aims_Comparison.qti | qtiplot file for comparing Vasp and FHI-AIMS           |
| Vasp_FHI-Aims_Comparison.pdf | pdf file for comparing Vasp and FHI-AIMS               |
| Dielectric/control.in        | The AIMS control file for the calculation epsilon inf  |
| Dielectric/geometry.in       | The AIMS geometry file for the calculation epsilon inf |
| Dielectric/aims.out          | The AIMS output containing epsilon infinity            |
| Vasp_FHI-Aims_Comparison.qti |                                                        |

[Back](..)
