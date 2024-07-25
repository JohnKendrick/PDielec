# Examples/Crystal/Quartz
## Quartz powder IR using CRYSTAL23

The presence of the BORN and HESSFREQ.DAT files allows a full precision calculation by PDielec.
The frequencies and the dielectric constants calculated by PDielec in in very good agreement with those provided by CRYSTAL23.
This calculation of the frequencies and intensities does not use CRYSTAL's CPHF method, so no ε<sub>∞</sub> is calculated.  
Instead, an estimate is provided in script.py.

 | Files               | Description                   |
 | ------------------- | ----------------------------- |
 | quartz.d12          | Input to Crystal23 |
 | quartz.out          | Output from the Crystal23 run |
 | BORN.DAT            | Contains the Born charge tensors from the calculation |
 | HESSFREQ.DAT        | Contains the hessian from the calculation |
 | command.pdmake      | pdmake command for test suite |
 | results.ref.xlsx    | Reference results for test suite |
 | script.py           | pdgui script to create reference results |

[Back](..)
