# Examples/Raman/AbInit
## ZnO Raman Calculation using AbInit
This calculation uses the PEAD method.  The raman.abi input files produces several databases.
The databases for steps 4 and 5 were merged before processing with anaddb.

```
mrgddb analysis.ddb ramano_DS4_DDB ramano_DS5_DDB
anaddb < analysis.files
```

The output in analysis.abo reports the Raman tensor in units of (amu_electron/bohr)^(1/2)
To convert to amu/Angstrom^(1/2) requires multiplying the analysis.abo figures by 31.0585. 
This gives:

|  **Group**  |   **Raw/PDielec** | **Sum-rule corrected** | **Printed rsus converted** |
| ------------|------------------ | -----------------------|--------------------------- |
|  E2 low     | 0.0730969         | 0.0730968              | 0.0730969                  |
|  A1 TO      | 0.0560663         | 0.0356425              | 0.0356426                  |
|  E1 TO      | 0.0002049         | 0.00004965             | 0.00004965                 |
|  E2 high    | 0.112456          | 0.112456               | 0.112456                   |


| **Directory**                   | **Description**                                       |
| ------------------------------- | ----------------------------------------------------- |
| opt.abi                         | AbInit optimisation input                             |
| opt.abo                         | AbInit optimisation output                            |
| raman.abi                       | AbInit Raman output                                   |
| raman.abo                       | AbInit Raman output                                   |
| analysis.abi                    | anaddb input for the analysis step                    |
| analysis.abo                    | anaddb output                                         |
| analysis.ddb                    | the merged database                                   |
| analysis.files                  | anaddb input                                          |
| command.pdmake                  | part of the pdmake tests                              |
| script.py                       | part of the pdmake test system                        |
| O.psp8,Zn.psp8                  | the pseudo-potentials used                            |

[Back](../README.md)

