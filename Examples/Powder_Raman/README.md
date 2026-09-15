# Powder Raman examples

These ZnO datasets exercise powder Raman calculations with several readers.

| Directory | Input |
| --- | --- |
| [AbInit](./AbInit/README.md) | ABINIT Raman and nonlinear response |
| [Castep](./Castep/README.md) | CASTEP Raman response |
| [Crystal23](./Crystal23/README.md) | CRYSTAL23 and companion tensor files |
| [QE](./QE/README.md) | Quantum ESPRESSO log, dynG and tensors.xml |
| [Vasp](./Vasp/README.md) | VASP and Raman-Tensors.yaml |
| [Finite_difference](./Finite_difference/README.md) | VASP finite-difference JSON, with EO off/on scenarios |

Run `pdmake command.pdmake` within the selected directory. Use
`pdmake --view command.pdmake` to open its saved scenarios interactively.
Review differences before replacing references with
`pdmake --regenerate command.pdmake`.

[Back to Examples](../README.md)
