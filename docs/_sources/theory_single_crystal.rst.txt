.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE

.. _Single-Crystal-Theory:

================================
Theory for Single Crystal Optics
================================

The *SingleCrystal tab* enables the calculation of the optical behaviour of a single crystal (thick slab or film).  The method uses the pyGTM code available on GitHub :cite:`pygtm`.  This code implements a generalised transfer matrix method described by Passler et al. :cite:`Passler2020` and which builds on previous publications :cite:`Passler2017,Passler2017a`.
The code has been modified to model incoherent systems as well as coherent ones.
A brief summary of the theory underlying the transfer matrix method is given to aid in the understanding of the methods implemented and the range of their applications.

In PDGui the number of media to consider is only three, a superstrate (usually air) through which the incident light passes, a dielectric medium which has a frequency dependent complex permittivity and optionally a substrate (which is also usually air).  
The media are indexed by 0,1 and 2 for the superstrate, dielectric and substrate respectively.
A schematic illustrating this is shown in :numref:`fig-definition-of-RTA`.
To aid understanding of the input and output and the limitations of the methods used, a brief summary of the transfer matrix method is given here, specific to its application in PDGui.
There are 4 different modes of operation of the matrix transfer method in PDGui and these are described below under the headings 'coherent thin film', 'incoherent thin film', 'partially incoherent thin film' and 'thick slab'


.. _fig-definition-of-RTA:

.. figure:: ./_static/Figures/Definition_of_Transmittance.png
   :scale: 90%

   Layers used to calculate Transmittance, Reflectance and Absorptance

Coherent thin film
------------------

At an interface the tangential (in-plane) electric and magnetic fields (given by a vector :math:`\fieldbf{F_i}`) have to match in both media.  
This requires that the 4 amplitudes (vector :math:`\fieldbf{A_i}`) of the electric and magnetic fields, forward and backward, s and p polarised, in medium :math:`i` are related through a 4x4 dynamical matrix :math:`\tensorbf{D_i}`.

.. math::
   :label: eq-dynamical-matrix

   \fieldbf{F}_i  = \tensorbf{D}_i \fieldbf{A}_i

Requiring that the in-plane fields match in both media at an interface gives;

.. math::
   :label: eq-field-matching
    
   \fieldbf{F}_{i-1}                    &= \fieldbf{F}_i                \\
   \tensorbf{D}_{i-1} \fieldbf{A}_{i-1} &= \tensorbf{D}_i \fieldbf{A}_i \\
   \fieldbf{A}_{i-1}                    &= \tensorbf{D}_{i-1}^{-1} \tensorbf{D}_i \fieldbf{A}_i

Within a given medium the mode amplitudes will change according to the propagation through the medium described by an diagonal exponential matrix :math:`\tensorbf{P}_i`.  Thus the amplitudes on the left side of a layer are related to those on the right side by;
      
.. math::
   :label: eq-propagation

   \fieldbf{A}_{i,left}  = \tensorbf{P}_i \fieldbf{A}_{i,right}

The total transfer matrix for the 3 media considered by PDGui (superstrate, dielectric and substrate; 0, 1 and 2 respectively) is therefore;

.. math::
   :label: eq-transfer-matrix

   \tensorbf{M}  = 
                 \tensorbf{D}_0^{-1}
                 (\tensorbf{D}_1
                 \tensorbf{P}_1
                 \tensorbf{D}_1^{-1})
                 \tensorbf{D}_2

The elements of the total 4x4 transfer matrix :math:`\tensorbf{M}` can be used to determine the total reflectance and transittance for each s and p mode.

Incoherent thin film
--------------------
Because the above approach uses amplitudes there is always the possibility that a mismatch in phase for a reflected or transmitted wave will result in interference and therefore oscillations in intensity as the frequency is changed.
Experimentally, for thick crystals, it is known that the presence of defects in the crystal and at the interfaces, and differences in thickness cause this coherence to be lost and instead incoherent light transmission is observed and the oscillations in intensity due to inteference are lost.
The "Incoherent thin film" mode creates a transfer matrix like that in :eq:`eq-transfer-matrix`, but instead of amplitudes, uses intensities for the transfer matrix of the dielectric. 

Such an approach has been implemented in FSRStools :cite:`FSRStools` a Python implementation of Yeh's transfer matrix method :cite:`Yeh1980`.
The theory is described by Beck :cite:`Beck2012` and by Katsidis and Siapkas :cite:`Katsidis2002`.
This approach has been incorporated in PDGui by modifying the total transfer matrix as follows;

.. math::
   :label: eq-transfer-matrix-int

   \tensorbf{M}_{int} =  \mid \tensorbf{D}_0^{-1} \tensorbf{D}_1 \mid^2 \mid \tensorbf{P}_1 \mid ^2 \mid \tensorbf{D}_1^{-1} \tensorbf{D}_2 \mid ^2

The use of intensities instead of amplitudes ensures that the phase information is lost.


Partially incoherent thin film
------------------------------

PDgui caters for partial incoherence, where some incoherence is introduced by lack of planarity, or uncertainties in the orientation angle of the crystal, by randomly averaging over spectra produced by varying the thickness, the orientation angles of the crystal and the incident light.
The sampling is done over a uniform distribution.  A single percentage is given for all parameters.  For the thickness the percentage is a percentage of the initial thicknes.  For angles the percentage is a percentage of 90 :math:`^\circ`.

Because of the random sampling in the partially incoherent case it has sometimes been found necessary to smooth the calculated spectra using a Savitzky-Golay filter. 
The parameters of the filter are the kernel size (the number of points to averaged) and a polynomial degree for the fitting of those points.  
If a kernel size of less than 2 is given, then no filtering is performed.  The kernel size should be and odd number.
The settings in PDGui typically start with a sample of 20 calculations and no filter.


Thick Slab
----------
For those cases where the partially incoherent thin film and incoherent thin film modes are not appropiate and only reflectance is of interest, the calculation mode referred to as a "thick slab" assumes that there are only two (sem-infinite) media.
The media through which the incident light travels and the crystalline media by which it is reflected.
Transmission is not of interest in this case, as in the case of a thick crystal it assumed that total absorption will take place and that there will be no internal reflection within the crystal.
This mode is similar to a standard Fresnel calculation of reflectance (R) but is appropriate for general permittivity tensors.
As a result of the assumptions associated with a "thick slab" it is assumed that transmittance (T) through a semi-infinite dielectric with some absorption will be zero and therefore the absorptance (A) will be;

.. math::
   :label: thick_absorptance

   A = 1 - R

For all other modes, the absorptance is defined as below;

.. math::
   :label: full_absorptance

   A = 1 - R - T

Comparison of computational approaches
--------------------------------------

:numref:`fig-mode-comparison` shows a comparison of the approaches discussed above to MgO.
The DFT calculations were performed by Castep.  The incident angle is 45\ :superscript:`o` and the film thickness is 1\ |micron|.
The 'Thick slab' mode can be regarded as a limiting case as the simulated slab becomes thicker, the calculations tend to follow the 'thick slab' curve more closely.
The 'Coherent thin film' curve shows large amplitude oscillations before and after the restrahlung region of absorption.  These are damped by including partial incoherence but are not completely removed.  
The partial incoherence settings in this case were; 100 samples and up to 10% deviation in incident angle, thickness and crystal orientation, with no smoothing.
The incoherent curve shows no oscillations and follows the thick mode curve near the restrahlung region.
It is common to see numerical problems in regions where the extinction coefficient is high and the dielectric is thick.  The numnerical problesm manifest themselves with sudden spikes in reflectance and indicate that only thinner samples should be considered.



.. _fig-mode-comparison:
.. figure:: ./_static/Figures/single_crystal_mode_comparison.png
   :scale: 90%

   Comparison of methods

.. _crystal-and-laboratory-coordinates:

Crystal & Laboratory Coordinate Frames
--------------------------------------
There are two coordinate systems to consider.  The first is the laboratory coordinate system (X,Y,Z).  
It is assumed that the normal to the crystal surface is aligned with the Z-axis and that the incident radiation is in the X-Z plane.  P-polarised light is polarised parallel to this plane and S-polarised perpendicular (senkrecht) to this plane.  
The laboratory frame is illustrated in :numref:`fig-lab-coords2`.

The other coordinate system of interest is that of the crystal.  The calculation of the permittivity tensor is performed in the crystal coordinate system (x,y,z).  The crystal plane is defined by a set of Miller indices (hkl) and the normal to the crystal surface is rotated to align with the laboratory Z-axis.  Finally in the laboratory frame the crystal may be rotated around the Z-axis by the azimuthal angle :math:`\phi`.
The precise definition of the azimuthal angle is somewhat arbitrary as it depends on the crystal unit-cell definition, but PDGui provides information of the details of the relationship between the crystal axes and the laboratory frame.
The arrangement described is illustrated in the Figure below where :math:`E_s` and :math:`E_p` are the field directions of the S- and P polarised incident light respectively. 
Information about the crystal coordinates relative to the laboratory frame is given in the GUI.


.. _fig-lab-coords2:
.. figure:: ./_static/Figures/SingleCrystalGeometry.png
   :scale: 90%

   Definition of single crystal laboratory coordinates in PDGui

