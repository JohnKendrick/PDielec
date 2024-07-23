# PDielec examples using GULP


| **File**  | **Description**                                  |
| ----------| ------------------------------------------------ |
| CaCO3     | Powder example, calculated using GULP 6.2 |
| Na2SO42   | Powder example, calculated using GULP 4.4 |
| preader -  |Check that preader is working on all Gulp examples |

## Comments
In order to perform gulp calculations it is necessary to either run without symmetry or to modify the borncharge.F90 code to print out the full born charge matrix, rather than that just for the asymmetric unit.
In any case it is important to increase the format specification for the printing in this routine from f10.4 to f14.8
