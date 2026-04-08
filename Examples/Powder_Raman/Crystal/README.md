# Examples/Raman/ZnO/Crystal
## Crystal 17 calculation on ZnO

The PDGui calculation in script.py calculates the powder Raman spectrum using a DFT calculation on ZnO 
The DFT calculation uses an TZVP basis and PBE functional with the D3 dispersion correction.
To calculate the Raman spectrum run `pdgui --spectroscopy powder_raman opt_raman.out`

| **File**              | **Description**                                   |
| --------------------- | ------------------------------------------------- |
| opt_raman.d12         | Input for a Raman calculation |
| opt_raman.out         | Output for a Raman calculation |
| optc020               | The optimised geometry |
| TENS_IR.DAT           | The IR tensor          |
| BORN.DAT              | The Born charge        |
| TENS_RAMAN.DAT        | The Raman tensor       |
| command.pdmake        | pdmake script for the test suite |
| results.ref.xlsx      | reference results for the test suite |
| script.py             | script used by pdgui for the test suite |

[Back](../README.md)

