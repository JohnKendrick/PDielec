.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


===========================
Single Crystal Applications
===========================

The transmission, reflection and absorption of infrared light through single crystals can be simulated using PDielec and PDGui.  Partial incoherence can be modelled by averaging over a sample of thicknesses and orientations.

Single Crystal Study of Thickness Effects in MgO
================================================

This example looks at an isotropic system.  The results of a CASTEP calculation of MgO have provided the frequency dependent permittivity and this is used to look at the effect of crystal thickness on the reflectance, transmittance and absorptance.i
The examples shown here were calculated using the scripts mgo_coherent_size_effects.py and mgo_partially_incoherent_size_effects in the Examples/SingleCrystal/MgO directory of the GitHub distribution.
A 45\ :superscript:`o` incident beam on the (001) surface of MgO is used to calculate the transmittance, reflectance and absorbtance of thin films with thickness of 0.1, 1, 10 and 100 microns.
The results of the coherent calculations are shown in :numref:`fig-mgo-coherent`

.. _fig-mgo-coherent:
.. figure:: ./_static/Figures/MgO_Coherent_Size_Effects.svg
   :scale: 90%

   Effects of Film Thickness of MgO in Coherent Light, a 45\ :superscript:`o` incident beam with p- polarisation

The results for the 'thick slab' case are shown only for reflectance and absorbtance.  It can be clearly seen that the reflection in the restrahlen region between the TO and LO frequencies (388 and 693 |cm1|) is almost total with also very little absorbtance.  
The 0.1\ |micron| thickness shows peaks in reflectance and corresponding drops in transmittance at both the TO and LO frequencies, whilst showing only absorbtance at the TO frequency.
The 1\ |micron| thickness shows considerable broadening of these signals, with the TO reflectance splitting into a doublet. 
At 10\ |micron| below the TO frequency there are oscillations in the reflectance, transmittance and absorptance which are a result of interference between the internally reflected waves.  They occur because the light is coherent.  The oscillations continue until they meet the 'thick slab' limit and they then follow this limit until frequencies above the LO frequency.
At 100\ |micron| the oscillations below the TO frequency are very fast and the reach the 'thick slab' limiting case at lower frequencmes than for the 10\ |micron| case.  Similarly at higher frequencies the follow the 'thick slab' limiting case to higher frequency.

The introduction of partial incoherence is acheived in the program by averaging over a range of orientation parameters and is therefore quite crude, but it gives an indication of the type of effects that partial incoherence may have.
The term partial incoherence in this case is used to refer to incoherence introduced by sample defects such as planarity and varying thickness.

:numref:`fig-mgo-partially-incoherent` shows the smoothing effect that averaging over the orientation and thickness parameters has.
For these calculations the following parameters were used.

.. table:: Partially Incoherent Parameters
   :name: tab-partially-incoherent-parameters
   :column-alignment:  left right  right
   :header-alignment:  left center center
   :column-dividers:   none single single single none

   +-----------------+------------------------+-------------------+
   | Thickness       | Percentage Incoherence | Number of samples |
   +-----------------+------------------------+-------------------+
   | 0.1\ |micron|   |         20             |       10          |
   +-----------------+------------------------+-------------------+
   | 1.0\ |micron|   |         20             |       10          |
   +-----------------+------------------------+-------------------+
   | 10.0\ |micron|  |         60             |       1000        |
   +-----------------+------------------------+-------------------+
   | 100.0\ |micron| |         30             |       1000        |
   +-----------------+------------------------+-------------------+


It is clear, especially from the absorbtance figure, that the thicker samples are a better fit to the 'thick slab' limiting case.  
However there are issues when using this facility.  In particular the high frequency side of the reflectance spectra is showing shoulders which are artefacts of the quite large perturbations in the orientation and thickness parameters which have been imposed.
Also, for the 10\ |micron| sample it is clear that these perturbations are not sufficient by themselves to remove the oscillations that occur at low frequencies.  
This is because at these frequencies the changes in thickness are not sufficent to move a peak by the period of the oscillations.

.. _fig-mgo-partially-incoherent:
.. figure:: ./_static/Figures/MgO_Partially_Incoherent_Size_Effects.svg
   :scale: 90%

   Effects of Partial Incoherence and Film Thickness of MgO, a 45\ :superscript:`0` incident beam with p-polarisation

The effect of thickness on the "Incoherent thin film Mode" is show in :numref:`fig-mgo-incoherent`.
This method uses intensities rather than amplitudes for the transfer matrix and as a result, all the phase information is lost, resulting in no interferences after internal reflections.
The calculations shown are for p- polarised incident radiation.


The incoherent method deviates from the thick slab mode at very low frequencies, but as the film thickness increases the deviations become smaller.
Similarly at frequencies above the LO frequency the deviations from the thick slab mode are less for the thicker films.
The user needs to be aware that for large thicknesses the calculation of the propagation matrix can lead to numerical overflows.  These manifest themselves as sudden changes in transmittance.  Reflectance calculations are badly affected by this and it is recommended to use the 'thick slab' mode for such cases.


.. _fig-mgo-incoherent:

.. figure:: ./_static/Figures/MgO_Incoherent_Size_Effects.svg
   :scale: 90%

   The Incoherent Thin Film Mode - a 45\ :superscript:`o` incident beam with p- polarisation




Komandin et al :cite:`Komandin2009` have reported terahertz spectra of 0.5\ |micron| thick MgO crystals which are dominated by interference at low frequencies.  
In :numref:`fig-mgo-experiment` a digitised spectrum obtained from this publication is compared with calculations.
It can be seen that despite the rather poor digitisation the agreement is very good, given that we are using only calculated information (refractive index, LO and TO frequencies etc.).
Agreement with the number of interference bands and the amplitude of the bands is very good.

.. _fig-mgo-experiment:
.. figure:: ./_static/Figures/MgO_Single_Crystal_0p5_Micron_Experiment.svg
   :scale: 90%

   Comparison of calculated and experimental transmittance terahertz spectra


Finally we compare a digitised spectrum of the reflectance measured over the range 0 - 800 |cm1| with the spectrum calculated using PDGui and Castep. 
In this the 'thick slab' approximation has been used as no oscillations are detected in the low frequency experimental results.
The calculated and experimental results are in excellent agreement.
The shoulder on the experimental result is due to two phonon processes which are not taken account of in the present theoretical work.


.. _fig-mgo-experiment-reflectance:
.. figure:: ./_static/Figures/MgO_Single_Crystal_Experiment_Reflectance.svg
   :scale: 90%

   Comparison of calculated and experimental reflectance infrared spectra



Single Crystal Study of Forsterite
==================================

In a publication on the interpretation of the experimental spectrum of forsterite :cite:`Pierre2013a` a Four Parameters Semi-Quantum (FPSQ) model for the single crystal infrared spectrum of this material was presented.   In this application note that we use the experimental file format of PDielec to read in the FPSQ model parameters from the paper and use it to calculate the reflectance spectrum for comparison with that published in the paper.
The FPSQ model is defined by the following equation.


.. math::
   :label: eq-fpsqa

    \epsilon (\omega )=\epsilon _{\infty}\prod_{j} \frac{\Omega^2_{LO_j}-\omega ^2-i\gamma _{LO_j}\omega }{\Omega^2_{TO_j}-\omega ^2-i\gamma _{TO_j}\omega}

Forsterite is an orthosilicate belonging to an orthorhombic space group.  The permittivity matrix is diagonal when the x, y and z axes of the crystal are aligned with a, b and c axes of the unit cell.  The parameters for the FPSQ model are summarised in the tables below.

.. table:: Optical Permittivity
   :name: tab-optical
   :column-alignment:  center right  right  right
   :header-alignment:  center center center center 
   :column-dividers:   none single single single none

   +------------+------------+--------+--------+
   |            |     X      | Y      |  Z     |
   +------------+------------+--------+--------+
   |    X       | 2.83       | 0.0    | 0.0    |
   +------------+------------+--------+--------+
   |    Y       | 0.0        | 2.69   | 0.0    |
   +------------+------------+--------+--------+
   |    Z       | 0.0        | 0.0    | 2.76   |
   +------------+------------+--------+--------+



.. table:: FPSQ Model Parameters for Forsterite
   :name: tab-fpsqa
   :widths:                   1      1      1      1      1      1      1      1      1      1      1      1
   :column-dividers:   single single single single single single single single single single single single single single
   :column-alignment:         center center center center center center center center center center center center
   :header-alignment:         right  right  right  right  right  right  right  right  right  right  right  right 

   +----------------------------------------+----------------------------------------+----------------------------------------+
   | EPS(xx)                                | EPS(yy)                                | EPS(zz)                                |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | TO Freq | TO sigma| LO Freq | LO sigma | TO Freq | TO sigma| LO Freq | LO sigma | TO Freq | TO sigma| LO Freq | LO sigma |
   +=========+=========+=========+==========+=========+=========+=========+==========+=========+=========+=========+==========+
   | 202.6   |0.18     |202.9    | 0.28     | 144.9   |0.09     |145.5    | 0.21     | 278.3   |0.50     | 278.9   | 0.54     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 276.6   |0.59     |277.1    | 0.84     | 278.6   |0.68     |278.9    | 0.79     | 292.7   |0.90     | 305.6   | 0.36     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 296.1   |0.79     |300.6    | 1.18     | 290.8   |1.28     |311.3    | 0.68     | 306.3   |0.39     | 316.7   | 0.53     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 321.1   |0.60     |322.5    | 0.88     | 351.7   |2.77     |376.7    | 1.02     | 411.6   |2.03     | 414.6   | 1.97     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 383.2   |1.52     |390.5    | 1.52     | 397.3   |2.15     |412.4    | 2.30     | 418.5   |3.19     | 456.4   | 1.88     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 405.4   |1.58     |468.9    | 2.90     | 419.8   |2.17     |444.1    | 2.75     | 425.0   |4.22     | 424.5   | 3.83     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 476.6   |4.46     |478.0    | 3.69     | 458.1   |4.39     |489.2    | 2.59     | 478.9   |3.22     | 487.2   | 3.99     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 504.7   |3.81     |532.6    | 5.65     | 507.8   |2.77     |513.7    | 4.16     | 506.2   |2.71     | 580.8   | 8.79     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 532.3   |5.67     |552.3    | 9.41     | 529.5   |4.62     |575.0    | 8.33     | 874.4   |4.17     | 996.8   | 4.32     |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 606.6   |4.77     |651.6    | 12.10    | 838.6   |6.56     |844.5    | 7.63     |         |         |         |          |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 840.8   |6.75     |841.5    | 7.03     | 872.9   |5.02     |965.2    | 3.33     |         |         |         |          |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 959.9   |2.41     |965.5    | 3.30     | 986.3   |4.80     |994.8    | 2.96     |         |         |         |          |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+
   | 978.5   |3.83     |1081.1   | 5.76     |         |         |         |          |         |         |         |          |
   +---------+---------+---------+----------+---------+---------+---------+----------+---------+---------+---------+----------+


The publication presents experimental polarised reflectance infrared spectra of single crystals of forsterite.  The incident beam angle to the sample was :math:`10^{\circ}`.  The s-polarised geometry was used and measurements made with the electric field parallel to the a-axis, b-axis and c-axis.

In PDGui the laboratory frame is defined by XYZ, with the incident and reflected light in the surface in the XZ plane, with the surface normal aligned with the Z-axis.
The surface is therefore in the XY plane.
The azimuthal angle can be used to rotate the crystal about the normal until the crystal axis of interest is parallel to the Y axis of the laboratory frame.
With this geometry the s-polarised wave field will be oscillating parallel to the chosen axis.
The arrangement described is shown in the Figure below where :math:`E_s` and :math:`E_p` are the field directions of the incident light. 

.. _fig-lab-coords:

.. figure:: ./_static/Figures/SingleCrystalGeometry.png
   :scale: 90%

   Definition of single crystal laboratory coordinates in PDGui



As can be seen :math:`E_s` is parallel to the laboratory Y-axis, so it is necessary to line up the crystal unit cell so that the cell direction being investigated is parallel to the Y-axis.
In the single crystal scenario tab PDGui shows the crystal axes in terms of the laboratory coordinates so it is relatively straightforward to make sure that the cell direction of interest is aligned along the Y-axis.  It is possible to set each axis in one of two ways, depending on the surface being used.  For example to make sure the a-axis is aligned with laboratory Y-axis, we can use either the [001] or the [010] surfaces and simply choose the appropriate azimuthal angle.

The following geometries can be used by PDGui to ensure the relevant axis is aligned with the laboratory Y-axis;

- For the a-axis; 
    | the [001] surface with azimuthal angle =  90
    | the [010] with azimuthal angle = 90
- For the b-axis
    | the [001] surface and azimuthal angle = 0
    | the [100] surface and azimuthal = 90
- for the c-axis
    | the [100] surface azimuthal angle = 0
    | the [010] surface azimuthal angle = 0


The two figures below show the reflectance on s-polarised light along the a-axis and these figures show very good agreement with Figure 1 in the published experimental work :cite:`Pierre2013a`.
The data for this calculation can be found by looking in the :ref:`Examples` section of the installation guide under the Experimental/forsterite example.


.. _forsterite:

.. figure:: ./_static/Figures/forsterite-100-1200.png
   :scale: 90%
.. figure:: ./_static/Figures/forsterite-350-600.png
   :scale: 90%

   Forsterite a-axis reflectance (upper figure 100 - 1200 |cm-1|, lower figure 350 - 600 |cm-1|)


Single Crystal Study of L-Alanine
=================================
In an experimental and computational study of the vibrational modes of single crystals of L-alanine, the experimental transmittance was reported for polarized radiation along each of the principal axes, a-, b- and c-. :cite:`Allen2023`

The directory Examples/SingleCrystal/L-alanine contains files for the calculation of the transmittance and reflectance of single crystals of L-alanine based on similar calculations to those reported in the paper, using Crystal17 to perform the DFT calculations, and using PDGui to calculate the transmittance and reflectance properties of the crystal.  
Details of the DFT calculation using the B97-3c functional can be found in phonon.dl2.  The output can be found in phonon.out.  Other files produced by the Crystal17 calculation are BORN.DAT and HESSFREQ.DAT.  These are used by PDGui to calculate the theoretical spectrum.

The experimental results for the measurement of terahertz transmittance were kindly provided by the authors of the experimental paper: J. L. Allen, T. J. Sanders, J. Horvat, R. A. Lewis and K. C. Rule :cite:`Allen2023`.
We would like to thank the authors of the paper for providing the experimental data and for helpful discussions.
Their measurements can be found in the spreadsheets a-axis-experimental.xlsx, b-axis-experimental.xlsx and c-axis-experimental.xlsx.   Each of these spreadsheets has two columns, the first is a column of frequencies in |cm-1| and the second is the transmittance.

Scripts used to generate the crystal transmittance for the electric field along the respective crystal axes from the calculated phonon modes can be found in a-axis.py, b-axis.py and c-axis.py.  These scripts use the *Fitter Tab* to compare the calculated spectrum with the experimental one.  Changes have been made to the film thickness and to the widths of the peaks to improve the agreement between the two.

To run the scripts use a command such as,::

    pdgui -script a-axis.py

In the following, the steps used to get to the scripts are explained.
We will concentrate on the a-axis spectrum but the approach is similar for each axis considered.

L-alanine crystallises in an orthorhombic form, with 3 non-equivalent axes. 
To initialise PDGui with the results of the DFT calculation and starting the GUI interface with the *Scenario tab* in its single crystal mode use the command below,::

    pdgui phonon.out -scenario crystal

The *Main Tab* shows the calculated unit cell and the frequencies as calculated by Crystal17.

.. _fig-lalanine-maintab:
.. figure:: ./_static/Figures/lalanine-maintab.png
   :scale: 90%

The *Settings Tab* (see below) shows those settings which affect the calculation of the permittivity.  By default the frequencies are modified by projecting out any translational degrees of freedom.  
This is shown by the three zero frequencies in the table of frequencies and is a result of selecting *Apply Eckart conditions?*
The atomic masses are selected by PDGui and maybe different to those used by the DFT calculation.
A Lorentzian line width is applied to all phonon vibrations by default.
In this case we have chosen a large value by default and applied smaller and more realistic values for those phonon modes important to transmission for the a-polarised electric field.

The contribution to the permittivity arising from electronic polarisation is shown in this Tab.
In the case of Crystal17 and the B9703c functional the CPHF calculations is not possible because the exchange/correlation functional is not supported and the :math:`\epsilon_{\infty}` calculation is not performed and a zero matrix will be presented here.
In this example we have typed in the contribution from a Crystal17 calculation using PBE-D3.  
The resulting spectrum is not very sensitive to these values.

.. _fig-lalanine-settingstab:
.. figure:: ./_static/Figures/lalanine-settingstab.png
   :scale: 90%

The next tab is the *Scenario tab* which is used to reflect the experimental parameters that are possible.  
In this case we are looking at the single crystal scenario.
The mode of operation in the example shown below is the *Coherent thin film* mode.
In this case is the surface that the THz beam is incident on is the (001) surface, so the c-axis will be pointing along the Z-axis of the laboratory coordinate frame.
The permittivity of the superstrate and substrate is set to 1.0, consistent with the experimental conditions.
The geometry of the incident beam relative to the crystal is controlled by the azimuthal angle.  A zero value for the incident angle means the beam is at 90\ :superscript:`o` to the surface.  A non-zero angle defines the laboratory coordinates as the beam lies in the XZ plane.
For this example the incident beam angle will be set to 0\ :superscript:`o`.
Once this is done information about the position of the crystal axes relative to the crystal frame can be seen in the text box of the *Settings tab*.
Changes to the azimuthal angle will result in alterations to the definitions of the crystal axes in terms of the laboratory coordinates.
We are interested in radiation polarised with the electric field active along the a-axis of the crystal, which, according to the *Lab frame information* lies along  the Y-axis of the laboratory coordinate.
This corresponds to s-polarised light, so we will be looking at the s-component of the transmittance.

An important parameter to consider is the thickness of the film and the method by which we will calculate the transmission and reflectivity.   
To start with we will use the *coherent thin film* mode and a thickness of 0.2mm.  Experimentally :cite:`Allen2023`, the crystal thickness was measured to be 0.5±0.1mm, but with this thickness the calculation showed little transmission in the chosen frequency range. 


.. _fig-lalanine-scenariotab:
.. figure:: ./_static/Figures/lalanine-scenariotab.png
   :scale: 90%


The resulting spectrum is shown below.  
It can be seen that there are regular oscillations in the transmittance, especially at low frequencies.
These are due to interference of the directly transmitted light and light that has undergone internal reflection and generating a series of decaying etalons.

.. _fig-lalanine-plottertab-coherent:
.. figure:: ./_static/Figures/lalanine-plottingtab-coherent.png
   :scale: 90%


To see if we can remove this we will change the mode in the *Scenario tab* to *Incoherent thin film*.
Re-plotting the spectrum and looking at the results in the *Fitter tab* is shown in the plot below.
The calculated spectrum is compared directly with the experimental one.  
Both the calculated and experimental spectra show little transmittance above about 140\ |cm-1|.
The strong, very sharp peaks at 225\ |cm-1| and above are probably numerical problems that occur with thicker films whentransmission tends to zero.  These numerical problems are exacerbated by the *Incoherent film* mode as it uses intensities in the calculation of the transfer matrix, which requires the square the electric field and the propagation matrices.
The fact that these are numerical problems can be verified by looking at the *Coherent thin film* mode and seeing if they are present there.
Even if the experimental measurements are designed to minimise etalons, it is generally sensible to compare both incoherent and coherent calculations before comparing to experiment.


A major advantage of using the *Fitter tab* to understand the calculated spectrum is that the user can modify the Lorentzian widths and see the result immediately, while also and comparing them with the experimental spectrum.
This also shows the cross-correlation between experiment and calculation, in this we get a value of 0.9379 (where 1 is a perfect match).


.. _fig-lalanine-plottertab-coherent-aaxis:
.. figure:: ./_static/Figures/lalanine-fittertab-incoherent-aaxis.png
   :scale: 90%


Plots for the b-axis and c-axis electric fields are shown below.  The b-axis example is given by using the same (001) surface but using an azimuthal angle of 0\ :superscript:`o`.  The thickness used for the b-axis example is 0.1mm, smaller than the 0.2mm used for the a-axis example.
The c-axis field transmittance example uses the (010) surface with an azimuthal angle of 0\ :superscript:`o` and a thickness of 0.1mm.

.. _fig-lalanine-plottertab-coherent-baxis:
.. figure:: ./_static/Figures/lalanine-fittertab-incoherent-baxis.png
   :scale: 90%


.. _fig-lalanine-plottertab-coherent-caxis:
.. figure:: ./_static/Figures/lalanine-fittertab-incoherent-caxis.png
   :scale: 90%

AlN on Silicon and Silicon Carbide
==================================

This is an example of a multilayer system where two layers are not isotropic.
In an experimental and theoretical study of the infrared reflactance of aluminium nitride on silicon and silicon carbide, MacMillan, Devaty and Choyke :cite:`MacMillan1993` used four paramater semi-quantum  (FPSQ) model to describe their experimental results on aluminium nitride on various substrates.
In the Examples/Experimental/AlN directory the file AlN.exp provides an experimental file which describes their FPSQ model in a format the PDGui can use.
In addition parameters for the permittivities of silicon and 6H-SiC are given in the TestDataBase.xlsx file in the same directory.  
The 6H-SiC permittivities are provided by a Drude-Lorentz model and the silicon permittivities is calculated from an experimental refractive index.
References for the origins of the models and experimental data are given in the spreadsheet.

The first experimental/calculated system considered :cite:`MacMillan1993` was the reflectance of a 0.92μ AlN on a Si substrate.  
This was modelled using PDGui after reading in the experimental file, AlN.exp, the *Single Crystal Scenario Tab* was used to define the system as described in the paper.
The *Thick slab* mode was specified, the superstrate material was chosen as air, an AlN dielectric film of 0.92μ was specified on top of a 1.0μm film of silicon.
An angle of incidence of 7.2 \textsuperscript{o} was used.
The (001) surface of AlN was defined, so the perpendicular to the surface aligns with the laboratory Z-axis.
Because the program is operating in *Thick slab* mode, the bottom layer is treated as a semi-infinite layer, so the size of the silicon layer specified in the GUI is irrelevant.
The scenario is shown below:

.. _fig-aln-on-si-etalons-scenario:
.. figure:: ./_static/Figures/AlN-on-Si-Etalons-Scenario.png
   :scale: 90%

Two comparisons are made with the calculated results of the published paper after digitising figures 1 and 2 of the paper.  
In the first the range of the plot extends for 0 to 6000 |cm-1| and clearly shows the etalons associated with the interference.
In the second the plot examines the region between 400 and 1200 |cm-1|.  
The agreement is excellent.

.. _fig-aln-on-si-etalons:
.. figure:: ./_static/Figures/AlN-on-Si-Etalons.svg
   :scale: 90%

.. _fig-aln-on-si-close-up:
.. figure:: ./_static/Figures/AlN-on-Si-Close-Up.svg
   :scale: 90%

Figure 3 of the paper shows the reflectance of of a 0.98μm film of AlN on the (0001) surface of 6H_SiC. 
To model this with PDGui two layers have been created using the *Layer Editor* as shown below.
As with the previous example the *Thick slab* mode has been specified, so the SiC layer is semi-infinite.
The layer editor shows two (001) layers with the c-axis of the crystals pointing along the Z-laboratory axis.
The top layer, next to the superstrate (air) is the *Dielectric film* (AlN in this case) and the semi-infinite layer beneath is has been specified by adding a new layer from the materials database, TestDataBase.xlsx.

.. _fig-aln-on-sic-layereditor:
.. figure:: ./_static/Figures/AlN_on_SiC_LayerEditor.png
   :scale: 90%

A comparison is given below between the calculated spectrum reported in the paper that produced by PDGui.

.. _fig-aln-on-sic:
.. figure:: ./_static/Figures/AlN-on-6H-SiC.svg
   :scale: 90%

A final comparison is made between Figure 4 of the paper and the results of PDGui.   
In this case a 0.24μm film of SiC is supported on a 0.56μm film of AlN, which is deposited on a substrate of Silicon.
For this example the *Coherent thin film* mode is used.
The layer editor for this system is shown below.

.. _fig-sic-on-aln-layereditor:
.. figure:: ./_static/Figures/SiC_on_AlN-LayerEditor.png
   :scale: 90%

A comparison is given below between the calculated spectrum reported in the paper that produced by PDGui.

.. _fig-sic-on-aln:
.. figure:: ./_static/Figures/6H-SiC-on-AlN.svg
   :scale: 90%

