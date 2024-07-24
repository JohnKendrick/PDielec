# ZnO powder IR using CRYSTAL14 (Default settings)

The CRYSTAL14 run by default projects translation from  the hessian.
The presence of the BORN and HESSFREQ.DAT files allows a full precision calculation by PDielec.
The hessian is symmetrised using the same method as CRYSTAL14.
The frequencies and the dielectric constants calculated by PDielec in in very good agreement with those provided by CRYSTAL14.
This calculation of the frequencies and intensities is not a CPHF calculation, so ε<sub>∞</sub> is provided by the script.

 | Files               | Description                   |
 | ------------------- | ----------------------------- |
 | ZnO_default.d12     | Input to Crystal14 |
 | ZnO_default.out     | Output from the Crystal14 run |
 | BORN.DAT            | Contains the Born charge tensors from the calculation |
 | HESSFREQ.DAT        | Contains the hessian from the calculation |
 | command.pdmake      | pdmake command for test suite |
 | results.ref.xlsx    | Reference results for test suite |
 | script.py           | pdgui script to create reference results |



