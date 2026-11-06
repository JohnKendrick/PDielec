# Powder Raman examples

These ZnO datasets exercise powder Raman calculations with several readers.

| Directory | Input |
| --- | --- |
| [AbInit](./AbInit/README.md) | ABINIT Raman and nonlinear response |
| [Castep](./Castep/README.md) | CASTEP Raman response |
| [Crystal23](./Crystal23/README.md) | CRYSTAL17 output and companion tensor files (legacy directory name) |
| [QE](./QE/README.md) | Quantum ESPRESSO log, dynG and tensors.xml |
| [Vasp](./Vasp/README.md) | VASP and Raman-Tensors.yaml |
| [Particle_Shape](./Particle_Shape/README.md) | CASTEP ZnO: sphere, needle, plates and sphere-limit ellipsoid |
| [Finite_difference](./Finite_difference/README.md) | VASP finite-difference JSON, with EO off/on scenarios |

All supported particle shapes use analytical rotational averaging in an isotropic
host, with shape axes fixed relative to the crystal. Particle frequency shifts,
mode mixing and available electro-optic (EO) corrections are retained. There is
no orientation-sample convergence setting; legacy `Raman orientation samples`
entries in saved scripts are ignored.

Run the complete seven-example suite from the repository root:

```sh
./pdmake test-powder_raman
```

Run `pdmake command.pdmake` within the selected directory. Use
`pdmake --view command.pdmake` to open its saved scenarios interactively.
Review differences before replacing references with
`pdmake --regenerate command.pdmake`.

[Back to Examples](../README.md)
