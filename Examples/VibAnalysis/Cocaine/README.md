# Cocaine phonon vibanalysis using Quantum Espresso

This example used auto-selection of the internal coordinates, along with vibrational mode decomposition analysis (--vmd).  ARD decomposition is switched off (--novmard)

```
vibanalysis -autosel --novmard -vmd qe phonon.dynG
```
&nbsp;
  
 | **File**          | **Description**               |
 | ----------------- | ----------------------------- |
 | phonon.dynG       | QE output file                |
 | phonon.log        | QE output file                |
 | phonon.nma        | QE output file                |
 | phonon.ph.in      | QE phonnon input file         |
 | phonon.scf.in     | QE phonon input file          |
 | command.pdmake    | pdmake commands for testing   |
 | phonon.nma.ref    | reference vibanlysis output   |
