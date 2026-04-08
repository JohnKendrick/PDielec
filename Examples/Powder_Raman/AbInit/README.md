# Examples/Raman/ZnO/AbInit
## AbInit 10.6 calculation on ZnO

The PDGui calculation in script.py calculates the powder Raman spectrum using a DFT calculation on ZnO 
The calculation uses an LDA functional with PAW pseudopotentials.
To calculate the Raman spectrum run `pdgui --spectroscopy powder_raman ZnO_raman.abo`

| **File**              | **Description**                                   |
| --------------------- | ------------------------------------------------- |
| ZnO_opt.abi           | Input for the optimisation (results used for Raman calculation |
| ZnO_raman.abi         | Input for the Raman calculation                                |
| ZnO_raman.abo         | Output from the Raman calculation                              |
| Zn.LDA_PW-JTH.xml     | Pseudopotential for Zn                                         |
| O.LDA_PW-JTH.xml      | Pseudopotential for O                                          |
| command.pdmake        | pdmake script for the test suite |
| results.ref.xlsx      | reference results for the test suite |
| script.py             | script used by pdgui for the test suite |

[Back](../README.md)

