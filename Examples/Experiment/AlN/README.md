# Multi-layered single crystal

This is an example of a multi layer system for AlN/SiC/Si
The system is described in the application notes, where comparison is made with experiment.

The permittivity of AlN is described by Drude-Lorentz model taken from Akasaki & Hashimoto, Solid State Comm. 5, 851 (1967) and  MacMillan, Devaty, Choyke, Applied Phys. Lett. 62, 750-752 (1993).
The permittivity of Si is given in DataBase.xlsx and is a constant refractive index of 3.42.
The permittvity of SiC is also given in DataBase.xlsx.  It is a Drude-Lorentz model details are given in the database.

  | File      		| Comment                                               |
  | ------------------- | ----------------------------------------------------- |
  | 6H-SiC.cif		| crystal structure of SiC				|
  | AlN.cif		| crystal structure of AlN				|
  | DataBase.xlsx	| Materials database with SiC and Si			|
  | AlN.exp		| Experimental file with the dielectric for AlN 	|
  | AlN_on_Si.py	| Script for AlN on Si as in the application notes 	|
  | SiC-on-AlN-on-Si.py	| Script for SiC/AlN/Si as in the application notes 	|
  | command.pdmake	| Command for testing 					|
  | script.py		| Script for testing 					|
  | results.ref.xlsx	| Reference xlsx file for testing 			|
