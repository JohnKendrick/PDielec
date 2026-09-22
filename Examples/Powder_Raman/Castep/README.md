# Examples/Raman/ZnO/Castep
## CASTEP 19.0 calculation on ZnO

The PDGui calculation in script.py calculates the powder Raman spectrum using a DFT calculation on ZnO 
The calculation uses the G06 dispersion correction and the PBE functional
To calculate the Raman spectrum run `pdgui --spectroscopy powder_raman ramand3.castep`

| **File**              | **Description**                                   |
| --------------------- | ------------------------------------------------- |
| ramand3.param         | Castep input parameters |
| ramand3.cell          | Castep geometry         |
| ramand3.castep        | Castep output           |
| ramand3.efield        | Castep efield output    |
| ramand3.phonon        | Castep phonon output    |
| command.pdmake        | pdmake script for the test suite |
| results.ref.xlsx      | reference results for the test suite |
| script.py             | script used by pdgui for the test suite |

[Back](../README.md)
