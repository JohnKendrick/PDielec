# Fluoro-apatite VASP calculation of powder IR

10% volume fraction of fluoro-apatite is supported in a ptfe matrix.  Maxwell-Garnett effective medium theory is used to examine the IR absorption of spheres, plates and needles.


| Files                 | Description                                               |
| --------------------- | --------------------------------------------------------- |
| apatite_opt.cif       | Crystal structure from the VASP optimised calculation     |
| INCAR                 | VASP input settings                                       |
| KPOINTS               | VASP k-point settings                                     |
| POSCAR                | VASP cell and atom geometry                               |
| POTCAR.index          | A pointer to the pseudo potentials used in the calculation|
| OUTCAR                | VASP output, analysed by PDielec                          |
| application_note.py   | script used in the application notes                      |
| application_note.xlsx | results used in the application notes                     |
| command.pdmake        | pdmake command for the test suite                         |
| results.ref.xlsx      | Reference results for the test suite                      |
| script.py             | pdgui script to create the reference results              |
