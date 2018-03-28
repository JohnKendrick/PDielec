Title: PDielec: The Calculation of Infrared and Terahertz Absorption for Powdered Crystals
Author: John Kendrick
email:john@kendrick.me.uk
author: Andrew Burnett
email:a.d.burnett@leeds.ac.uk

[INCLUDE="style"]
[TITLE]

The Python package PDielec is described, which calculates the infrared absorption characteristics of a crystalline material supported in a non-absorbing medium. PDielec post processes solid state quantum mechanical and molecular mechanical calculations of the phonons and dielectric response of the crystalline material. Packages supported are Abinit, Castep, Crystal14, Gulp, QuantumEspresso, Phonopy and VASP. 

~ Center
![#img-1]
[#img-1]: Figures/Polarisation_Schematic.png {width:60%}
~

Using an effective medium method, the package calculates the internal electric field arising from different particle morphologies and calculates the resulting shift in absorption frequency and intensity arising from the coupling between a phonon and the internal field. The theory of the approach is described, followed by a description of the implementation within PDielec. For the specific case of a low concentration of spherical particles, calculations based on Mie scattering allow the exploration of particle size effects. A section providing several examples of its application is given.

# Theory

~ TexOnly
~~ Snippet
\noindent The documentation includes a section on the underlying \href{Theory.pdf}{theory}. 
~~
~

~ HtmlOnly
<div>
The documentation includes a section on the underlying <a href="Theory.html">theory</a>.
</div>
~

# Applications Notes

~ TexOnly
~~ Snippet
\noindent There is a set of \href{Application_Notes_1.pdf}{"applications notes"}
and examples, which illustrate the use of the program
~~
~

~ HtmlOnly
<div>
There is set of 
<a href="Application_Notes_1.html">application notes</a>
and examples, which illustrate the use of the program
</div>
~

# Supporting programs
A couple of supporting commands have been developed with use the same underlying libraries to read the output files from the various QM/MM packages.
~ TexOnly
~~ Snippet
\noindent \href{preader.pdf}{"preader"} reads a list of input files and summarises the results in output which can be written to a csv file.
\noindent \href{pdgui.pdf}{"pdgui"} provides a gui interface to the pdielec package
\noindent \href{phonana.pdf}{"phonana"} Has been replaced by pdgui.  It reads an output file and calculates the contribution of molecular centre of mass and rotation to each phonon mode
~~
~

~ HtmlOnly
<div>
<a href="preader.html">preader</a>
reads a list of input files and summarises the results in output which can be written to a csv file.
<a href="pdgui.html">pdgui</a>
provides a gui interface to the pdielec package
<a href="phonana.html">phonana</a>
Has been replaced by pdgui.  It reads an output file and calculates the contribution of molecular centre of mass and rotation to each phonon mode
</div>
~
