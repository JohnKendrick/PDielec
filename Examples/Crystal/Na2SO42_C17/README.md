# Na2(SO4)2 powder IR using CRYSTAL17

The CRYSTAL17 run by default projects translation from  the hessian.
The presence of the BORN and HESSFREQ.DAT files allows a full precision calculation by PDielec.
In command.sh "-hessian crystal" is specified.  This makes PDielec symmetrise the hessian using
the same method as CRYSTAL17.  The frequencies and the dielectric constants calculated by
PDielec in in very good agreement with those provided by CRYSTAL17.
This calculation of the frequencies and intensities is a CPHF calculation, so CRYSTAL has calculated ε<sub>∞</sub> as well as the frequencies.

 | Files               | Description                   |
 | ------------------- | ----------------------------- |
 | Na2SO42.d12         | Input to Crystal17 |
 | Na2SO42.out         | Output from the Crystal17 run |
 | Na2SO42_crystal.cif | Optimised crystal structure   |
 | BORN.DAT            | Contains the Born charge tensors from the calculation |
 | HESSFREQ.DAT        | Contains the hessian from the calculation |
 | command.pdmake      | pdmake command for test suite |
 | results.ref.xlsx    | Reference results for test suite |
 | script.py           | pdgui script to create reference results |

