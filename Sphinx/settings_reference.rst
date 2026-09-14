.. include:: preamble.txt

.. _pdgui-settings-reference:

PDGui Settings Reference
========================

This reference describes the keys in each tab's ``settings`` dictionary. Use
these exact names, including spaces and capitalisation, in a saved PDGui
Python script. The :doc:`pdgui` guide describes the corresponding windows.

*Save settings* on the Main Tab provides a working starting point. For example::

    tab = self.notebook.settingsTab
    tab.settings['Sigma value'] = 5.0

    tab = self.notebook.scenarios[0]
    tab.settings['Legend'] = 'ZnO, parallel polarisations'

    tab = self.notebook.plottingTab
    tab.settings['Minimum frequency'] = 100.0
    tab.settings['Maximum frequency'] = 800.0
    tab.settings['Frequency increment'] = 1.0

Run the saved script with ``pdgui -script script.py``. Scenario indices start
at zero. Set up the required spectroscopy and scenario type before editing
its settings; changing a type string alone does not reconstruct a tab.

Values below are Python values: use ``True`` and ``False`` for switches,
quoted strings for choices, and lists for arrays. Defaults refer to a newly
created tab and can be replaced by the input file or a restored session.
Entries marked *derived*, *legacy* or *reserved* are included because they
can occur in the dictionary, but are not independent calculation controls.

Main Tab
--------

Use ``self.notebook.mainTab.settings``.

``Program``
    Input reader: ``'Abinit'``, ``'Aims'``, ``'Castep'``, ``'Crystal'``,
    ``'Experiment'``, ``'Finite_field'``, ``'Gulp'``, ``'Phonopy'``,
    ``'Quantum espresso'``, ``'Vasp'`` or ``'PDGui'``. The last choice reads
    a saved session. See :doc:`software` for supported inputs, including
    :ref:`experimental-file-format` and :ref:`finite-field-json-format`.

    Menu labels, saved reader names and CLI options are distinct:

    .. list-table:: Reader names
       :header-rows: 1

       * - Menu label
         - Saved ``settings['Program']``
         - CLI ``-program`` value (pdgui and preader)
       * - Quantum Espresso
         - ``'Quantum espresso'``
         - ``qe``
       * - Finite field
         - ``'Finite_field'``
         - ``finite_field``

``Output file name``
    Name of the calculation output to read. Saved sessions use the notebook's
    working directory for associated files; retain companion input files
    required by the reader.

``Excel file name``
    Destination spreadsheet name, ending in ``.xlsx``. An empty string means
    no spreadsheet has been selected.

``Script file name``
    Filename used by *Save settings* for the Python session script.

``Compatibility mode``
    Platform metadata, normally ``'Linux'`` or ``'Window'``. The current tab
    records this value but does not use it to select a calculation method.

Settings Tab
------------

Use ``self.notebook.settingsTab.settings``. These settings are common to all
scenarios. See :doc:`theory_powder` for the dielectric response from normal
modes and :doc:`PowderRaman` for Raman tensors and activities.

``Spectroscopy type``
    ``'Powder Infrared'``, ``'Powder ATR'``, ``'Crystal Infrared'``,
    ``'Powder Raman'`` or ``'Crystal Raman'``. Selects the family of calculations
    and available scenario controls. Use a saved session of the required type
    or the Settings Tab spectroscopy selector to establish matching scenarios.

``Eckart flag``
    Boolean, default ``True``. Projects out translational motion when
    rebuilding the dynamical matrix, helping to enforce the acoustic modes.

``Neutral Born charges``
    Boolean, default ``False``. Enforces charge neutrality on the Born
    effective charges before calculating the dielectric response.

``Sigma value``
    Common Lorentzian half-width in cm⁻¹, default ``5``. Supplies the mode
    broadening; individual linewidths are held in ``tab.sigmas_cm1``.

``Mass definition``
    ``'average'`` uses average atomic masses; ``'program'`` uses masses from
    the calculation; ``'isotope'`` uses the isotope mass table; ``'gui'`` uses
    user-entered masses from ``tab.masses_dictionary``. Default ``'average'``.
    Mass changes affect the normal modes and their response.

``Optical permittivity``
    The dimensionless high-frequency dielectric tensor as a 3 × 3 nested
    list. Initially ``None``; normally obtained from the reader. Used for the
    dielectric background, local fields and NAC calculations. See
    :doc:`theory_powder` and :doc:`CrystalRaman`.

``Optical permittivity edited``
    Boolean recording a manual tensor override. Set ``True`` when supplying
    ``Optical permittivity`` explicitly so that refresh and session saving
    preserve the edited tensor.

``Symmetrise optical permittivity``
    Boolean, default ``True``. Applies the available crystal symmetry to the
    optical tensor; ``False`` retains the unsymmetrised reader tensor.

``Raman activity units``
    ``'polarizability'`` displays activities in Å⁴/amu; ``'epsilon'`` uses
    the internal ``R_epsilon`` activity convention in Å/amu. This changes
    the Settings table and Raman polar display, not the physical scenario
    spectrum. See :doc:`PowderRaman` for the tensor conventions.

The selected-mode mask, individual linewidths, masses and NAC diagnostic
controls are not additional ``settings[]`` keys. Saved scripts may also
contain assignments to tab attributes such as ``modes_selected``,
``sigmas_cm1`` and ``masses_dictionary``; retain those assignments when
reproducing a session.

Scenario Tabs: Common Settings
------------------------------

Use ``self.notebook.scenarios[index].settings``. The following keys are
inherited by powder and crystal scenarios.

``Legend``
    Text used to identify the scenario in plot legends and output.

``Scenario type``
    Records the scenario family. This is session metadata; use the notebook's
    scenario creation and spectroscopy controls to change the type.

``Materials database``
    Path to the Excel materials database, normally ``MaterialsDataBase.xlsx``.
    Provides support matrices and optical layer materials. See the materials
    database description in :doc:`software`.

Powder Infrared and Powder ATR
------------------------------

Use a powder scenario's ``settings`` dictionary. Most particle and support
settings are also present in Powder Raman, with the qualifications below.
See :doc:`theory_powder` for effective-medium and particle-shape theory.

``Matrix``
    Support-material name in the materials database, initially ``'ptfe'``.
    ``'Material defined manually'`` selects a user-defined constant support.
    Powder Raman also accepts ``'none'``: this uses the bulk TO Raman response
    without support local-field or particle NAC/EO corrections.

``Matrix density``
    Support density in g/cm³. Normally supplied by the database; a positive
    manual value is used to convert between mass and volume fractions.

``Matrix permittivity``
    Dimensionless constant support permittivity for a manually defined
    material; may be complex. Database materials can instead provide a
    frequency-dependent response.

``Mass fraction``
    Crystallite mass fraction between zero and one, not a percentage.
    The corresponding volume fraction is calculated using the densities.

``Volume fraction``
    Crystallite volume fraction between zero and one, not a percentage.
    Controls the effective-medium mixture and scales the powder Raman
    response. Default ``0.1``.

``Mass or volume fraction``
    ``'mass'`` or ``'volume'`` identifies which fraction is the independent
    input when converting the other; default ``'volume'``.

``Bubble volume fraction``
    Fraction of additional air voids, between zero and one; default ``0``.
    Used in the powder infrared mixture calculation.

``Bubble radius``
    Air-void radius in µm, default ``30``. Used by the infrared bubble
    correction when the bubble fraction is non-zero.

``Effective medium method``
    ``'Maxwell-Garnett'`` or ``'Bruggeman'`` selects an effective-medium
    model; ``'Averaged Permittivity'`` averages the crystallite response;
    ``'Mie'`` includes finite-size spherical scattering. These choices govern
    powder infrared/ATR, not the separate Raman tensor calculation.
    See :doc:`theory_powder`.

``Particle shape``
    ``'Sphere'``, ``'Needle'``, ``'Plate'`` or ``'Ellipsoid'``. Determines
    the depolarisation tensor and hence local fields and particle response.

``Particle size(mu)``
    Particle radius in µm, default ``0.0001``. Relevant to finite-size
    infrared scattering, particularly Mie; not a Raman mode linewidth.

``Particle size distribution sigma(mu)``
    Width parameter for the infrared particle-size distribution, default
    ``0`` for a single size. This is the lognormal distribution parameter,
    not a Lorentzian frequency width.

``Ellipsoid a/b``
    Axial aspect ratio for an ellipsoid; default ``1`` gives a sphere.

``Unique direction - h``
    First component of the particle's unique direction.

``Unique direction - k``
    Second component of the particle's unique direction.

``Unique direction - l``
    Third component of the particle's unique direction. Together the three
    keys specify a direct-lattice [abc] axis for a needle/ellipsoid or an
    (hkl) surface normal for a plate. The default is ``[0, 0, 1]``;
    spheres do not use this direction. See :doc:`theory_powder`.

``ATR material refractive index``
    Real refractive index of the ATR prism, default ``4``. Powder ATR only;
    see :ref:`ATR_theory`.

``ATR theta``
    ATR incidence angle in degrees, default ``45``.

``ATR S polarisation fraction``
    Fraction of s-polarised illumination: ``0`` is pure p, ``1`` is pure s,
    and ``0.5`` is an equal mixture. Used to combine ATR responses.

Powder Raman
------------

These keys supplement the common scenario and powder particle settings.
See :doc:`PowderRaman` for orientation averaging, internal fields,
particle-corrected phonons and the electro-optic contribution.

``Raman laser wavelength``
    Laser wavelength in nm. A positive
    number, default ``785``; determines the excitation and Stokes frequency
    factors. Do not supply a frequency in Hz or cm⁻¹.

``Raman laser polarisation``
    ``'VV'`` selects parallel polarisations, ``'VH'`` and ``'HV'`` select
    crossed polarisations, and ``'Unpolarised'`` selects the unpolarised
    response. Default ``'VV'``.

``Raman temperature``
    Temperature in kelvin, default ``298``. Controls the Stokes thermal
    population factor.

``Raman orientation samples``
    Number of numerical orientations for non-spherical particles. The GUI
    offers ``4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096``; default
    ``512``. Increase to check convergence. Spheres use an analytic average.

``Raman electro-optic term``
    Boolean, default ``False`` for a new powder scenario. Includes the
    electro-optic contribution when particle correction data and χ⁽²⁾ are
    available. Has no effect for ``Matrix = 'none'``. See :doc:`PowderRaman`.

Crystal Infrared and Crystal Raman: Layers
------------------------------------------

Use a crystal scenario's ``settings`` dictionary. See
:doc:`theory_single_crystal` for the optical stack, orientation conventions
and coherent/incoherent propagation, and :doc:`CrystalRaman` for Raman fields.

``Global azimuthal angle``
    Rotation of the stack about its surface normal, in degrees; default ``0``.

``Angle of incidence``
    Incident angle from the surface normal in the superstrate, in degrees;
    default ``0``. See :ref:`crystal-and-laboratory-coordinates`.

``Mode``
    ``'Scattering matrix'`` (default) or ``'Transfer matrix'`` for the
    infrared optical solver. Scattering matrices avoid numerical problems
    that transfer matrices can encounter in strongly absorbing thick layers.
    See :doc:`theory_single_crystal`.

``Layer material names``
    List of database material names, in order from superstrate to substrate.
    ``'Dielectric layer'`` denotes the material read from the calculation.
    The initial stack is ``['air', 'Dielectric layer', 'air']``.

``Layer hkls``
    One three-component (hkl) surface normal per layer. Used for anisotropic
    layers; isotropic layers do not need a crystal orientation.

``Layer azimuthals``
    One in-plane rotation in degrees per layer, applied with its surface
    orientation and the global azimuthal angle.

``Layer thicknesses``
    One thickness per layer, in the corresponding ``Layer thickness units``.
    The first and last layers are semi-infinite ports, so their stored
    thicknesses do not represent finite films.

``Layer thickness units``
    List of length-unit strings: ``'ang'`` (Å), ``'nm'``, ``'um'``, ``'mm'``
    or ``'cm'``. The GUI offers the latter four; scripts can also use
    ``'ang'``. Keep the list aligned with the thicknesses.

``Layer incoherent options``
    One optical propagation choice per layer: ``'Coherent'`` retains phase;
    ``'Incoherent (intensity)'`` propagates intensities;
    ``'Incoherent (phase averaging)'`` averages sampled phases;
    ``'Incoherent (phase cancelling)'`` removes propagation phase; and
    ``'Incoherent (non-reflective)'`` treats a non-reflecting thick layer.
    See :doc:`theory_single_crystal`. These are optical layer choices,
    distinct from the Raman ``Depth coherence`` and ``Layer combination``.

``Layer dielectric flags``
    Boolean list marking layers that use the loaded dielectric material.
    Keep these flags consistent with ``Layer material names``; the initial
    three-layer stack uses ``[False, True, False]``.

All seven layer lists must have the same length. Saved sessions restore the
layer objects from these lists; live layer edits should use the layer editor
so that both the objects and their serialised settings remain consistent.

``Percentage average incoherence``
    Phase range as a percentage of a full 2π cycle, default ``100``.
    Used for layers with phase averaging.

``Number of average incoherence samples``
    Number of phase samples per phase-averaged layer, default ``10``.
    Multiple such layers require combinations of their sampled phases.

``Percentage partial incoherence``
    Strength of the infrared partial-incoherence sampling, default ``0``
    (disabled). Randomly perturbs the active crystal thickness, orientation
    and incidence angle; this is separate from phase averaging.

``Partially incoherent samples``
    Number of random samples for that infrared average, default ``20``.

``Filter kernel size``
    Savitzky–Golay smoothing window for the partially incoherent infrared
    spectrum. Default ``1`` disables smoothing; an enabled window must be
    odd, exceed the polynomial order, and fit within the spectrum length.

``Filter polynomial size``
    Polynomial order for that smoothing filter, default ``3``. Relevant only
    when the kernel size is greater than two.

``Frequency units``
    Legacy stored value, initially ``'wavenumber'``. The current crystal tab
    does not consume it; select display units with Plotting's singular
    ``Frequency unit`` key.

Crystal Raman: Scattering and Coherence
---------------------------------------

These settings supplement the crystal layer settings. See
:doc:`CrystalRaman` for the scattering geometry, reciprocal fields,
NAC/EO corrections and the different coherence choices.

``Laser wavelength nm``
    Positive laser wavelength in nm, default ``532``.

``Incident polarisation``
    ``'p'`` or ``'s'`` selects the incident channel; default ``'p'``.
    Mixed coherent polarisations can be represented by a two-component
    complex Jones vector in p/s order, as produced by a saved Porto setup.

``Detected polarisation``
    ``'p'`` or ``'s'`` selects an analyser, while ``'unpolarised'`` sums the
    two detected intensities. Default ``'p'``. A two-component complex
    p/s Jones vector specifies a coherent analyser combination. Summing
    unpolarised intensities is different from summing field amplitudes.

``Temperature K``
    Temperature in kelvin, default ``298``, used in the Stokes population factor.

``GL point density``
    Gauss–Legendre integration-point density per µm of active material,
    default ``20``. Increase to check depth-integration convergence.

``Collection side``
    ``'superstrate'`` collects backscattered light on the incident side;
    ``'substrate'`` collects forward-scattered light. Default ``'superstrate'``.

``Collection angle``
    Collection angle in degrees. ``-1`` (default) requests automatic geometry;
    other values specify a signed angle explicitly. Automatic collection uses
    minus the incidence angle on the superstrate side (retro-backscattering)
    and plus the incidence angle on the substrate side (forward collection).
    An explicit positive angle equal to incidence selects the historical
    specular-reflection geometry on the superstrate side.

``Layer combination``
    ``'Incoherent intensities'`` (default) sums the layer intensities;
    ``'Coherent amplitudes'`` sums complex layer amplitudes before squaring.
    Coherent layer combination is unavailable with incoherent depth integration.

``Depth coherence``
    ``'Coherent amplitude'`` (default) integrates amplitudes over depth before
    squaring. ``'Incoherent intensity'`` integrates local intensities. This is
    a Raman emission choice, not the optical layer's propagation coherence.

``Approximate ES``
    Script-only boolean, default ``False``. If ``True``, substitutes the laser
    field for the separately calculated scattered field. This changes the
    optical approximation; see :doc:`CrystalRaman`.

``Coalesce equivalent layers``
    Script-only boolean, default ``True``. Merges adjacent equivalent coherent
    layers before Raman integration. ``False`` retains the original layer
    boundaries and is useful for subdivision checks.

``Raman electro-optic term``
    Boolean, default ``True`` for crystal scenarios. Includes χ⁽²⁾ in the
    NAC-corrected Raman tensors when those data are available. Disabling it
    retains the NAC frequency correction without the EO tensor contribution.

``Layer NAC mode``
    ``'none'`` uses bulk TO modes; ``'geometry'`` (default, displayed as
    *Snell's law*) obtains phonon momentum from the macroscopic geometry;
    ``'dominant_mode'`` uses the dominant internal Berreman mode;
    ``'modal_pairs'`` (*All modes*) resolves the internal optical mode pairs.
    GUI labels may append *(EO)*; that suffix is not part of the stored value.

``Modal pair combination``
    Used with ``Layer NAC mode = 'modal_pairs'``. ``'Group q channels'``
    (default) combines amplitudes sharing a phonon momentum and detector
    channel. ``'Incoherent pairs'`` squares pairs separately for comparison
    with older calculations. ``'Coherent all pairs'`` is a diagnostic option
    that mixes distinct phonon final states. See :doc:`CrystalRaman`.

``Modal pair final-state model``
    For incoherent-depth modal-pair calculations, ``'Bulk phase matched'``
    (default) selects internal momenta relative to the external transfer.
    ``'Local incoherent'`` assigns local field components to the externally
    selected phonon state without the q-angle filter.

``Modal pair q-angle tolerance``
    Maximum angle in degrees between internal and external momentum transfer,
    from ``0`` to ``90``. Default ``90`` accepts the forward hemisphere in
    the bulk phase-matched model. Inactive for local incoherent depth.

``Modal pair include zero q``
    Reserved/legacy entry, initially ``None``. It is retained in the dictionary
    but is not read by the current scenario calculation; it does not control
    zero-momentum channel selection.

``Azimuthal sweep points``
    Positive integer number of angles over 0–360°, excluding the repeated
    endpoint; default ``36``. Used when opening an azimuthal sweep.

``Porto notation``
    Stored normal-incidence geometry text, initially empty. For example,
    ``'x(yy)x'`` denotes the backscattering shorthand and ``'x(yz)+x'``
    explicit forward scattering. Applying the Porto entry in the GUI sets
    the layer normal, angles, collection side and polarisations. Preserve
    these accompanying values in a script; assigning the text alone is not
    a replacement for applying the geometry. See :ref:`crystal-raman-gui`.

``Phonon boundary correction``
    Legacy key from older sessions. With ``Layer NAC mode = 'none'``, the
    old values ``'NAC'``, ``'slab-environment'`` and ``'slab'`` map to
    ``'geometry'``. Use ``Layer NAC mode`` in new scripts.

Plotting Tab
------------

Use ``self.notebook.plottingTab.settings``.

``Minimum frequency``
    Lower calculation limit in cm⁻¹, default ``1``. Stored in wavenumbers
    regardless of the chosen display unit.

``Maximum frequency``
    Upper calculation limit in cm⁻¹, default ``200``.

``Frequency increment``
    Calculation-grid spacing in cm⁻¹, default ``0.2``. Use sufficiently small
    spacing to resolve the selected mode linewidths.

``Frequency unit``
    Display unit: ``'wavenumber'`` (cm⁻¹), ``'THz'``, ``'GHz'``, or the
    wavelength units ``'ang'``, ``'nm'``, ``'um'``, ``'mm'``, ``'cm'``, ``'m'``.
    Wavelength is reciprocal to frequency; the calculation-grid settings
    above remain in cm⁻¹.

``Molar definition``
    ``'Unit cells'`` (default), ``'Atoms'`` or ``'Molecules'`` sets the amount
    used to normalise powder molar absorption. See :doc:`theory_powder`.

``Number of atoms``
    Number of atoms per molecule for ``Molar definition = 'Molecules'``;
    default ``1``. Used with the unit-cell atom count and volume to calculate
    the molecular concentration.

``Plot type``
    Exact curve name. Powder infrared offers ``'Powder Molar Absorption'``,
    ``'Powder Absorption'``, ``'Powder Real Permittivity'`` and
    ``'Powder Imaginary Permittivity'``. Powder ATR additionally offers
    ``'Powder ATR'``. Raman uses ``'Powder Raman'`` or ``'Crystal Raman'``.
    Crystal infrared offers ``'Crystal Reflectance (P polarisation)'``,
    ``'Crystal Reflectance (S polarisation)'``,
    ``'Crystal Transmittance (P polarisation)'``,
    ``'Crystal Transmittance (S polarisation)'``,
    ``'Crystal Absorbtance (P polarisation)'`` and
    ``'Crystal Absorbtance (S polarisation)'``. Preserve the historical
    spelling *Absorbtance* in scripts. See :doc:`theory_powder`,
    :ref:`ATR_theory`, :doc:`theory_single_crystal`, :doc:`PowderRaman`
    and :doc:`CrystalRaman` for the corresponding observables.

``Spectrum renormalisation``
    ``'none'`` retains the spectrum's units; ``'max=1'`` divides each curve
    by its peak; ``'area=1'`` divides by its integrated area.
    ``'layer depth'`` is available for Crystal Raman and divides by the
    total active-layer depth. These are plotting normalisations.

``Scenarios to plot``
    List of booleans, one per scenario in notebook order, for example
    ``[True, False, True]``. Controls visibility. An initially empty list is
    expanded by the tab to match the scenarios.

``concentration``
    Derived molar concentration in mol/L for the selected molar definition.
    Recalculated from the cell volume, atom count and molar definition;
    set those inputs rather than this value directly.

``cell concentration``
    Derived unit-cell concentration in mol/L, used in molar normalisation.

Analysis Tab
------------

Use ``self.notebook.analysisTab.settings``. See :doc:`analysis` for the
internal/external and molecular decomposition of modes.

``Radii``
    List of covalent radii in Å, in the tab's element order. Initially
    ``None``; populated from the element table and editable thereafter.

``Minimum frequency``
    Lower frequency shown in the analysis plot, in cm⁻¹; default ``-1``.

``Maximum frequency``
    Upper frequency shown in the analysis plot, in cm⁻¹; default ``400``.

``title``
    Analysis plot title, default ``'Analysis'``. This key is lower case.

``Covalent radius scaling``
    Multiplier on the sum of the atomic covalent radii when deciding bonds;
    default ``1.1``. Changes molecular grouping and therefore mode analysis.

``Bonding tolerance``
    Additional distance in Å added to the scaled sum of radii when deciding
    bonds; default ``0.1``. Bonding information is also used by the Viewer.

``Bar width``
    Width of the mode bars in the analysis plot, default ``0.5``.

3D Viewer Tab
-------------

Use ``self.notebook.viewerTab.settings``. These are display settings; they do
not change the calculated spectrum. Colours use ``[red, green, blue, alpha]``
with integer components from zero to 255.

``Selected mode``
    One-based phonon number to display/animate; initially ``4``. This differs
    from zero-based Python array and scenario indices.

``Atom scaling``
    Multiplier on covalent radii for atom spheres, default ``0.5``.

``Maximum displacement``
    Maximum displayed displacement in Å, default ``1``. Scales the vibration
    for viewing rather than changing the normal-mode eigenvector.

``Bond colour``
    Bond-cylinder RGBA colour, default ``[80, 80, 80, 255]``.

``Bond radius``
    Bond-cylinder radius in Å, default ``0.1``.

``Cell colour``
    Unit-cell edge RGBA colour, default ``[255, 0, 0, 255]``.

``Cell radius``
    Radius of the cell-edge cylinders in Å, default ``0.1``.

``Text colour``
    Cell-label RGBA colour, default ``[255, 255, 0, 255]``.

``Text size``
    Cell-label font size, stored as a string; default ``'24'``.

``Background colour``
    Viewport background RGBA colour, default ``[120, 120, 120, 255]``.

``Arrow colour``
    Displacement-arrow RGBA colour, default ``[0, 255, 0, 255]``.

``Arrow radius``
    Displacement-arrow radius in Å, default ``0.07``.

``Number of phase steps``
    Number of displacement steps for the animation, default ``41``. Must be
    odd; the Viewer increases an even value by one.

``Super Cell``
    Three positive integer repetitions along the cell vectors, for example
    ``[2, 2, 1]``. Default ``[1, 1, 1]``.

``Transform``
    3 × 3 lattice transformation, normally stored as nested lists of strings.
    Default identity. Changes the displayed cell; see
    :ref:`viewer-transform-window`.

``Primitive transform``
    Derived transformation obtained from the reader when the primitive-cell
    action is used. Normally preserve the value saved by that action.

``hkl``
    Three-component surface normal used to orient the view, default
    ``(0, 0, 1)``. See :ref:`crystal-and-laboratory-coordinates`.

``uvw``
    Three-component direct-lattice direction used to orient the view,
    default ``(1, 0, 0)``.

``Element palette``
    ``'Jmol'`` (default) or ``'Vesta'`` selects the element colour table;
    other names fall back to the built-in elemental palette.

``Element colours``
    Optional list of RGBA colours in the Viewer species order, overriding
    the palette for those species. Initially ``None``.

``Toggle states``
    Five booleans in this order: cell labels, cell edges, orientation axes,
    bonds, atoms. Default ``[True, True, True, True, True]``.

Fitter Tab
----------

Use ``self.notebook.fitterTab.settings``. The Fitter compares an experimental
spectrum with the selected calculated spectrum; see the Fitter section in
:doc:`pdgui` for the fitting workflow and experimental spectrum columns.

``Experimental file name``
    Filename of the measured spectrum to compare, initially empty. This is
    a frequency/intensity spectrum, distinct from an ``.exp`` material input.

``Plot title``
    Comparison plot title, initially
    ``'Experimental and Calculated Spectral Comparison'``.

``Fitting type``
    ``'Minimise x-correlation'`` selects the correlation-based objective;
    ``'Minimise spectral difference'`` selects the spectral RMSE objective.
    The optimisation adjusts selected linewidths and, optionally, frequency scale.

``Number of iterations``
    Optimisation iteration control, default ``20``.

``Frequency scaling factor``
    Multiplier on calculated frequencies in the comparison, default ``1``.
    Does not rescale the underlying normal modes in the input reader.

``Optimise frequency scaling``
    Boolean, default ``False``. Includes the frequency scaling factor among
    the variables adjusted during fitting.

``Independent y-axes``
    Boolean, default ``True``. Plots experimental and calculated spectra on
    independent y axes; ``False`` uses a shared intensity axis.

``Spectral difference threshold``
    Threshold applied after intensity normalisation when calculating the
    spectral difference, default ``0.05``. Values below it are set to zero
    in the comparison objective.

``Baseline removal``
    Boolean, default ``False``. Applies a Hodrick–Prescott background removal
    to the resampled experimental spectrum.

``HPFilter lambda``
    Base-10 logarithm of the Hodrick–Prescott smoothing parameter, default
    ``7`` (a coefficient of 10⁷). Used when baseline removal is enabled.

``Scenario index``
    Zero-based index of the scenario being fitted; initially the last
    scenario in the notebook.

``Spectrum scaling``
    Legacy boolean, initially ``False``. The current Fitter stores this
    value but does not apply it in the plotting or fitting calculation.

``Spectrum scaling factor``
    Legacy multiplier, initially ``1``. Retained with ``Spectrum scaling``;
    it is not an active intensity-scaling control in the current Fitter.
