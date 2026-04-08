# Examples/Raman/ZnO/VASP
## VASP 5.4.4 calculation on ZnO

The PDGui calculation in script.py calculates the powder Raman spectrum using a DFT calculation on ZnO 
The calculation uses VASP as the DFT package to calculate the Raman tensor numerically.
The Raman calculation itself is driven by the Skelton group's scripts for calculating the Raman efficiency

The following are valid ways of running a raman powder analysis:

```
pdgui --spectroscopy powder_raman OUTCAR        # Reads OUTCAR and Raman-Tensors.yaml
pdgui --spectroscopy powder_raman vasprun.xml   # Reads vasprun.xml and Raman-Tensors.yaml
```

| **File**              | **Description**                                   |
| --------------------- | ------------------------------------------------- |
| INCAR                 | An INCAR file for an IBRION=6 VASP calculation    |
| KPOINTS               | KPOINTS                                           |
| POSCAR                | The optimised structure of ZnO                    |
| POTCAR.index          | The PAW potentials used                           |
| OUTCAR                | The output from VASP |
| Raman-Tensors.yaml    | The Raman tensors calculated by Skeltons' scripts |
| Raman.yaml            | A yaml used by Skeltons' scripts |
| vasprun.xml           | A Vasp XML file                  |
| command.pdmake        | pdmake script for the test suite |
| results.ref.xlsx      | reference results for the test suite |
| script.py             | script used by pdgui for the test suite |

[Back](../README.md)
