# Na2(SO4)2 powder IR using CRYSTAL14

The CRYSTAL14 run by default projects translation from  the hessian.
The presence of the BORN and HESSFREQ.DAT files allows a full precision calculation by PDielec.
PDielec symmetrises the hessian using the same method as CRYSTAL14.  The frequencies and the dielectric constants calculated by
PDielec in in very good agreement with those provided by CRYSTAL14.
This calculation of the frequencies and intensities is a CPHF calculation, so CRYSTAL as calculated ε<sub>∞</sub> as well as the frequencies.

 | Files               | Description                   |
 | ------------------- | ----------------------------- |
 | Na2SO42.d12         | Input to CRYSTAL14 |
 | Na2SO42.out         | Output from the CRYSTAL14 run |
 | Na2SO42_crystal.cif | Optimised crystal structure   |
 | BORN.DAT            | Contains the Born charge tensors from the calculation |
 | HESSFREQ.DAT        | Contains the hessian from the calculation |
 | command.pdmake      | pdmake command for test suite |
 | results.ref.xlsx    | Reference results for test suite |
 | script.py           | pdgui script to create reference results |
