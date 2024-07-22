# Mie Scattering Example for MgO using Castep

A log normal distribution of sizes is used in this example.  The particle size is 1 micron.  Three calculations are performed.  One for a single particle size and the other two for lognormal width distributions of 0.1 and 0.5 microns.

| File                 | Description                                                 |
| -------------------- | ----------------------------------------------------------- |
| application_note.py  | Script used in generating data for the application note     |
| application_note.xlsx| Spread-sheet used in the application note     |
| command.pdmake       | pdmake command for testing suite |
| phonon.castep        | Castep output from phonon calculation |
| phonon.cell          | Castep unit cell for phonon calculation |
| phonon.param         | Castep input for phonon calculation |
| phonon.phonon        | Castep output from phonon calculation |
| results.ref.xlsx     | Reference results for testing suite |
| script.py            | Script to generate results for testing suite |
