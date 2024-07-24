# Leucine powder IR using CRYSTAL14

Not all CRYSTAL methods for calculating the frequencies and intensities also calculate ε<sub>∞</sub>.
In this example ε<sub>∞</sub> is calculated using a CPHF calculation (no frequencies) and then provided to PDGui in the script, script.py.
There is no BORN.DAT or HESSFREQ.DAT files, so the precision of the calculation by PDGui is limited by the number of significant figures used to express the normal modes in the output file.


 | Files                                              | Description                   |
 | -------------------------------------------------- | ----------------------------- |
 | LEUCINE_DIELECTRC_PBED_6311Gdp.out                 | CPHF calculation of ε<sub>∞</sub> |
 | LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out | Calculation of frequencies    |
 | command.pdmake                                     | pdmake command for test suite |
 | results.ref.xlsx                                   | Reference results for test suite |
 | script.py                                          | pdgui script to create reference results |

Thanks to Tim Kortner's group for providing this example for CRYSTAL14
