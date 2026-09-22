# Examples/Raman/ZnO/Crystal
## Crystal 17 calculation on ZnO

The PDGui calculation in script.py calculates the powder Raman spectrum using a DFT calculation on ZnO 
The DFT calculation uses an TZVP basis and PBE functional with the D3 dispersion correction.
To calculate the Raman spectrum run `pdgui --spectroscopy powder_raman opt_raman.out`
`TENS_RAMAN.DAT` supplies the polarizability derivatives needed for Raman tensors.
For the electro-optic correction, PDielec reads the nonlinear susceptibility table
embedded in `opt_raman.out`. An optional `CHI2.DAT` in the same directory
overrides that table; the reader converts its `d(MKS)` values using chi(2) = 2d.

Here `CHI2.DAT` duplicates the embedded susceptibility table. It is retained as
an example of the optional companion format, but can be omitted without
changing the parsed susceptibility or the EO correction.

| **File**              | **Description**                                   |
| --------------------- | ------------------------------------------------- |
| opt_raman.d12         | Input for a Raman calculation |
| opt_raman.out         | Output for a Raman calculation |
| optc020               | The optimised geometry |
| TENS_IR.DAT           | The IR tensor          |
| BORN.DAT              | The Born charge        |
| TENS_RAMAN.DAT        | The Raman tensor       |
| CHI2.DAT              | The non-linear susceptibility tensor       |
| command.pdmake        | pdmake script for the test suite |
| results.ref.xlsx      | reference results for the test suite |
| script.py             | script used by pdgui for the test suite |

[Back](../README.md)

