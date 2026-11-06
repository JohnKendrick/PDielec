# Examples/Powder_Raman/Finite_difference

[Back to Powder_Raman](../README.md)

VASP finite-difference ZnO powder Raman example comparing electro-optic contributions off and on. The local `finite_difference.json` is also used by the infrared and Raman reader examples. See the [example notes](example.md) and [calculation provenance](../../Finite_difference/finite_difference_examples.md).

Run the supplied example from this directory:

```sh
../../../pdmake command.pdmake
```

To open the PDGui example interactively:

```sh
../../../pdmake --view command.pdmake
```

Both scenarios use ZnO spheres in PTFE (volume fraction 0.1), a 785 nm laser,
VV polarisation and 298 K, with analytical orientation averaging. The test
writes `results.xlsx` and compares it with `results.ref.xlsx`; it is included
in `pdmake test-powder_raman`.
