Title: PDielec The Calculation of Infrared and Terahertz Absorption for Powdered Crystals
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

~ TexOnly
~~ Snippet
\noindent The documentation includes a section on the underlying \href{Theory.pdf}{theory}.  There is also a set of \href{Application_Notes_1.pdf}{"applications notes"}
~~
~

~ HtmlOnly
<div>
The documentation includes a section on the underlying <a href="Theory.html">theory</a>.

There is also a set of 
<a href="Application_Notes_1.html">application notes</\>
</div>
~
