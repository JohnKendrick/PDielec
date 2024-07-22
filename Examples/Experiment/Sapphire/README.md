# Sapphire Single Crystal

This directory holds the files that were used to write the Sapphire application note for single crystal reflectivity.
The use of an experimental file for air is simply a trick to allow the rest of the materials to be specified through the database functionality of the code.

The permittivity of Sapphire is described by an FPSQ model taken from Lutterotti and Scardi, Journal of Applied Crystallography 1990, 23, 246-252.  Further details are available in the database.

| File                        | Description                                                                    |
| --------------------------- | ------------------------------------------------------------------------------ |
| air.exp                     | Experimental file for air                                                      |
| a-planes.py                 | Script for the reflectance of all a-planes                                     |
| r-planes.py                 | Script for the reflectance of all r-planes                                     |
| c-planes.py                 | Script for the reflectance of all c-planes                                     |
| reflectance-a-planes.py     | Script for calculating a-plane reflectance as a function of azimuthal angle    |
| reflectance-c-plane.py      | Script for calculating c-plane reflectance as a function of azimuthal angle    |
| reflectance-r-planes.py     | Script for calculating r-plane reflectance as a function of azimuthal angle    |
| DataBase.xlsx               | Database with the FPSQ model for sapphire                                      |
| command.pdmake              | pdmake file for the testing suite                                              |
| script.py                   | pdgui script for the testing suite                                             |
