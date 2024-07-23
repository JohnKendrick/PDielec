# F-Apatite phonon vibanalysis using VASP

This example used auto-selection of the internal coordinates.

```
vibanalysis -autosel vasp OUTCAR
```
&nbsp;
  
 | **File**          | **Description**               |
 | ----------------- | ----------------------------- |
 | INCAR             | VASP input settings           |
 | KPOINTS           | VASP k-point settings         |
 | POSCAR            | VASP cell and atom geometry   |
 | POTCAR.index      | VASP pseudo potentials        |
 | OUTCAR            | VASP output                   |
 | command.pdmake    | pdmake commands for testing   |
 | phonon.nma.ref    | reference vibanlysis output   |
