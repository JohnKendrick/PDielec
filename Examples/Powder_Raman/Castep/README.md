# Examples/Powder_Raman/Castep
## CASTEP 24.1 calculation on ZnO

The PDGui calculation in script.py calculates the powder Raman spectrum using a DFT calculation on ZnO 
The calculation uses the G06 dispersion correction and the PBE functional.
The saved scenario is a sphere in PTFE, averaged analytically.
To calculate the Raman spectrum run `pdgui --spectroscopy powder_raman raman.castep`

| **File**              | **Description**                                   |
| --------------------- | ------------------------------------------------- |
| raman.param         | Castep input parameters |
| raman.cell          | Castep geometry         |
| raman.castep        | Castep output           |
| raman.efield        | Castep efield output    |
| raman.phonon        | Castep phonon output    |
| command.pdmake        | pdmake script for the test suite |
| results.ref.xlsx      | reference results for the test suite |
| script.py             | script used by pdgui for the test suite |

Run the example from this directory:

```sh
../../../pdmake command.pdmake
../../../pdmake --view command.pdmake
```

The first command compares `results.xlsx` with `results.ref.xlsx`; the second
opens the saved scenarios interactively. Loading only the raw input with
`pdgui` does not restore all settings in `script.py`.

[Back to Powder Raman](../README.md)
