# Examples/Crystal_Raman/Castep
## CASTEP 24.1 Crystal Raman calculation on ZnO

The PDGui calculation in `script.py` calculates the single-crystal Raman spectrum using a DFT
calculation on ZnO (wurtzite structure, 4 atoms per cell, 12 phonon modes).
The DFT calculation uses the G06 dispersion correction and the PBE functional.

To run the calculation interactively:
```
pdgui --spectroscopy crystal_raman raman.castep
```

To run the test suite:
```
pdmake command.pdmake
```

To regenerate the reference results after changing `script.py`:
```
pdmake --regenerate command.pdmake
```

## Test scenarios

The `script.py` defines six scenarios that exercise the key coherent/incoherent choices
available in the Crystal Raman scenario:

| **Scenario** | **Legend** | **What is tested** |
| ------------ | ---------- | ------------------ |
| 1 | Coherent depth - thin film 1um p-in | Per-layer `Coherent` mode; depth `Coherent amplitude`; p-in / unpolarised out; backscattering (superstrate collection); normal incidence |
| 2 | Incoherent depth - thick slab 100um | Per-layer `Incoherent (non-reflective)` mode; depth `Incoherent intensity`; p-in / unpolarised out; backscattering |
| 3 | s-in s-out polarisation 1um | Per-layer `Coherent`; depth `Coherent amplitude`; s-in / s-out polarisation channel |
| 4 | Forward scattering - substrate collection 1um | Per-layer `Coherent`; depth `Coherent amplitude`; collection side = `substrate` (transmission geometry) |
| 5 | Oblique incidence 45deg p-in 1um | Per-layer `Coherent`; depth `Coherent amplitude`; angle of incidence = 45° |
| 6 | NAC geometry correction 1um | Per-layer `Coherent`; depth `Coherent amplitude`; `Layer NAC mode = geometry` (Snell's law phonon wavevector) |

All scenarios use:
- ZnO with HKL = [0, 0, 1] (c-axis surface normal)
- Laser wavelength 532 nm, temperature 298 K
- Scattering matrix GTM method
- GL point density 20 points/µm

## Files

| **File**         | **Description**                                      |
| ---------------- | ---------------------------------------------------- |
| raman.param      | CASTEP input parameters                              |
| raman.cell       | CASTEP geometry                                      |
| raman.castep     | CASTEP output (main input file for PDielec)          |
| raman.efield     | CASTEP electric field / Born charge output           |
| raman.phonon     | CASTEP phonon output                                 |
| command.pdmake   | pdmake script for the test suite                     |
| results.ref.xlsx | Reference results for the test suite                 |
| script.py        | PDGui script defining the six test scenarios         |

[Back](../README.md)
