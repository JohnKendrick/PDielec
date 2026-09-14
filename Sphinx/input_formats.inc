Tabulated permittivity
----------------------

The *interpolate* models read a frequency followed by real and imaginary parts
of the permittivity. Frequencies are in |cm-1| and should increase strictly;
provide enough points for interpolation and cover the range to be plotted.
Finish each block with ``end`` or ``&end``. Commas or spaces separate values.

.. list-table:: Interpolation columns
   :header-rows: 1
   :widths: 1 3

   * - Directive
     - Columns after the frequency
   * - ``interpolate``
     - Real and imaginary parts of an isotropic permittivity.
   * - ``interpolate_3``
     - Real/imaginary pairs for xx, yy and zz, in that order.
   * - ``interpolate_6``
     - Real/imaginary pairs for xx, yy, zz, xy, xz and yz, in that order.
       The off-diagonal entries are copied to their transposed positions.

For example, an isotropic table is written as::

    interpolate
    100.0  2.10  0.01
    200.0  2.15  0.02
    300.0  2.20  0.03
    400.0  2.25  0.02
    end

A tabulation supplies the permittivity itself, not an oscillator strength.
The *constant* model instead takes three rows of three complex numbers, using
Python's ``j`` notation for the imaginary unit. The *fpsq* and *drude-lorentz*
models above require three consecutive sections, xx, yy and zz, including a
count for each. Their frequencies, damping parameters and Drude-Lorentz
strength parameters are in |cm-1|. The Drude-Lorentz strength is squared in
:ref:`the model expression <drude>`; it is not an integrated intensity.

Vibrational and Raman data in experimental files
------------------------------------------------

An experimental file can also supply frequencies, normal modes, Born charges
and Raman tensors. This makes it useful for importing results from a program
which PDielec does not read directly. The contents determine the available
features: a dielectric model alone supports optical calculations, whereas
normal-mode analysis and animation require eigenvectors. Raman spectra require
Raman tensors; polar-mode corrections additionally need the structure, masses,
complete normal modes, Born charges and optical permittivity.

Put the lattice, species and unit cell first, followed by ``epsinf`` (or its
synonym ``static``), then the following blocks. Keep identical atom and mode
ordering in all blocks. Use a complete set of :math:`3N` modes for a cell of
:math:`N` atoms, including acoustic and silent modes. List modes in ascending
frequency order and use negative frequencies for imaginary modes.

.. list-table:: Vibrational blocks
   :header-rows: 1
   :widths: 1 3

   * - Header
     - Following data
   * - ``frequencies M``
     - M lines, each starting with a frequency in |cm-1|. An optional scalar
       strength supplies all three diagonal oscillator strengths; three
       strengths supply xx, yy and zz. Omitting strengths sets them to zero.
       These are oscillator-strength tensor entries in atomic units (squared
       electron charge divided by electron mass), not Raman
       activities. When complete modes and Born charges are supplied the
       reader calculates oscillator strengths from them.
   * - ``born_charges N``
     - Three rows of three real numbers per atom, giving the Born effective
       charge tensor in electron-charge units. Use the atom order of
       ``unitcell``. An optional atom label before each tensor is ignored.
   * - ``normal_modes M``
     - N rows of x, y, z components per mode. These are mass-weighted
       eigenvectors, normalised so that the sum of squared components is one.
       They must be mutually orthogonal and consistent with the supplied masses.
   * - ``raman_tensors M epsilon``
     - Three rows of three real numbers per mode, giving
       :math:`R_\epsilon=\sqrt{V}\,\partial\epsilon_\infty/\partial Q` in
       :math:`(\mathrm{Angstrom}/\mathrm{amu})^{1/2}`.
       There is no frequency line or mode-number line between tensors.

The optional unit word after ``raman_tensors M`` is significant.
``epsilon`` stores the tensor directly. ``castep``, the default when the word
is omitted, multiplies the input by :math:`4\pi`: its input is the
polarizability-volume derivative divided by :math:`\sqrt{V}`. ``skelton``
multiplies the polarizability-volume derivative in
:math:`\mathrm{Angstrom}^2/\sqrt{\mathrm{amu}}` by
:math:`4\pi/\sqrt{V}`. Use explicit ``epsilon`` when exporting PDielec tensors.
Do not insert Raman activities in place of tensor components.

Blank lines and comment lines beginning with ``#`` in the first column are
permitted between data rows. Supply exactly the declared number of rows;
interpolation blocks need their terminator. Lattice vectors are multiplied by
the scale on the line after ``lattice`` and the resulting lengths are in
Angstrom. ``unitcell`` coordinates are fractional. ``species`` masses are in
amu. The reader also accepts CP2K-style ``&CELL`` and ``&COORD`` blocks, as
illustrated by the existing ``Examples/Experiment/cp2k`` examples.

The experimental reader preserves supplied normal modes; changing masses does
not provide the same eigenvector recalculation as an ab initio dynamical
matrix. Use masses and eigenvectors from the same calculation. For general
mass/Eckart reprojection, use the JSON format below with force constants or a
complete reference mode set. The experimental format currently has no
``chi2`` input block; use JSON to supply the nonlinear response for EO.

Two downloadable examples are provided:
:download:`a constant dielectric <_static/inputs/constant.exp>` and
:download:`ZnO with frequencies, modes, Born charges and Raman tensors <_static/inputs/zno-raman.exp>`.
The latter was converted from the CASTEP ZnO example; it is not measured data.
To open it, select *Experiment* in the *Main Tab*, or run::

    pdgui experiment zno-raman.exp -spectroscopy powder_raman

The ``.exp`` format describes a material. It is different from the two-column
experimental spectrum read by the *Fitter Tab*.

.. _finite-field-json-format:

Finite-Difference and Finite-Field JSON Format
==============================================

The *Finite_field* input program reads a consolidated JSON description of a
structure and its vibrational and Raman response. A finite-displacement or
finite-field workflow using VASP, Quantum ESPRESSO, CASTEP or another program
can write this common format. PDielec reads the results; it does not launch
the displaced or finite-field electronic-structure calculations.

Choose *Finite_field* and *Finite field (*.json)* in the *Main Tab* file
selector, or use::

    pdgui finite_field zno-raman.json -spectroscopy crystal_raman

Recognised schema names are exactly ``pdielect-finite-field-raman-v1``,
``pdielect-r_epsilon-v1`` and ``pdielect-vASP-consolidated-raman-results-v1``.
The spelling and case are part of these identifiers. An unfamiliar schema
produces a warning and the reader attempts to read the fields below. Automatic
program detection recognises schema strings beginning with ``pdielect-`` and
containing ``raman`` or ``r_epsilon``; specify ``finite_field`` explicitly for
other schema names. JSON uses double-quoted keys, ordinary arrays and numbers,
and ``true``/``false`` booleans. Comments and complex-number literals are not
accepted.

Structure and reference response
--------------------------------

The following is a minimal structural/optical example. It contains no phonons
or Raman response and is intended to illustrate the required structure keys.

.. literalinclude:: _static/inputs/minimal.json
   :language: json

.. list-table:: Structure and reference fields
   :header-rows: 1
   :widths: 2 4

   * - Key
     - Meaning
   * - ``structure.lattice_angstrom``
     - Required 3 by 3 array: rows are the a, b and c lattice vectors in Angstrom.
   * - ``structure.symbols``
     - Required list of N element symbols, one per atom.
   * - ``structure.positions_fractional``
     - N by 3 fractional coordinates. Takes precedence over
       ``structure.positions_cartesian_angstrom``, the alternative N by 3
       Cartesian coordinates in Angstrom. One must be present.
   * - ``structure.volume_angstrom3``
     - Optional cell volume. Defaults to the absolute determinant of the lattice.
       If supplied it must agree with the lattice, since Raman normalisation
       and electrostatic corrections use this volume.
   * - ``program_masses_amu``
     - Optional N masses or a dictionary keyed by element symbol. Takes
       precedence over ``masses_amu``, then the ``mass_amu`` entries of an
       ``atoms`` list. Masses default to PDielec average masses when absent.
   * - ``epsilon_infinity``
     - 3 by 3 relative optical permittivity, dimensionless. If absent at the
       top level, read from ``reference.epsilon_infinity``; otherwise defaults
       to zeros. Supply physical values for optical/NAC calculations.
   * - ``born_effective_charges``
     - N by 3 by 3 Born effective charge tensors in electron-charge units,
       in structure order. The ``reference`` block is the fallback. Missing
       tensors default to zero, eliminating the corresponding polar coupling.
   * - ``chi2``
     - Optional real 3 by 3 by 3 nonlinear susceptibility in PDielec's internal,
       cell-dependent EO convention. This array is stored directly; it is
       **not** a raw pm/V tensor or a Voigt 3 by 6 d-tensor.

For a physical susceptibility in pm/V the internal ``chi2`` is

.. math::

   \widetilde\chi^{(2)}_{ijl} = \chi^{(2)}_{ijl}[\mathrm{pm/V}]
   \frac{\sqrt{\mathrm{amu}/m_e}\,e\,10^8}
        {4\pi\epsilon_0\sqrt{V[\mathrm{Angstrom}^3]}}.

Convert a reported :math:`d=\chi^{(2)}/2` tensor to :math:`\chi^{(2)}` first.
The meaning and limitations of this electronic response are explained in
:ref:`Crystal-Raman-Theory`. A non-zero array alone is insufficient for EO:
the polar correction also needs the Born charges, optical permittivity,
normal modes and the relevant scenario option enabled.

Modes, force constants and dielectric derivatives
-------------------------------------------------

``modes`` is a list of objects, one per mode. ``frequency_cm-1`` is a real
frequency; ``imaginary: true`` makes a positive supplied frequency negative.
Modes are sorted by this signed frequency. Each mode may contain an N by 3
``mass_weighted_eigenvector`` and a 3 by 3 ``raman_tensor``. Eigenvectors must
be orthonormal and refer to the declared program masses. Raman tensors use
:math:`R_\epsilon=\sqrt{V}\,\partial\epsilon_\infty/\partial Q`, including
:math:`\sqrt V`, in :math:`(\mathrm{Angstrom}/\mathrm{amu})^{1/2}`.
All modes must contain a given array for that array set to be loaded.
Provide the complete :math:`3N` mode set for recalculation and mode mixing.

Unweighted force constants may replace or supplement the supplied modes.
The first available representation in the following list is used:

1. ``force_constants_matrix_sym_eV_per_angstrom2``: a 3N by 3N matrix.
2. ``force_constants_matrix_eV_per_angstrom2``: the same shape.
3. ``force_constants_eV_per_angstrom2``: an N by 3 by N by 3 array,
   indexed by atom, Cartesian component, atom, Cartesian component.
4. ``cartesian_displacements``: central finite-difference force records.

Matrix row and column order is atom 0 x, y, z, atom 1 x, y, z, and so on.
The reader symmetrises the matrix, converts eV/Angstrom² to atomic units,
and diagonalises it using the active masses. This recalculates frequencies
and eigenvectors, rather than retaining conflicting supplied frequencies.
The active masses initially use PDielec average masses; the *Settings Tab*
can select the original program masses instead.

Cartesian dielectric derivatives are supplied as ``deps_dr``, an N by 3 by 3
by 3 array with indices atom, displacement component, dielectric row,
dielectric column, in inverse Angstrom. The fallbacks, in order, are
``diagnostics.deps_dr`` and dielectric central differences in
``cartesian_displacements``. With active modes :math:`u_{mak}` the reader
forms

.. math::

   R^{(m)}_{ij}=\sqrt V\sum_{ak}
       \frac{u_{mak}}{\sqrt{M_a}}\frac{\partial\epsilon_{ij}}{\partial r_{ak}}.

These derivatives take precedence over supplied mode Raman tensors during
reprojection. When only projected tensors are supplied, also include reference
eigenvectors so their phases and mode mixing can be tracked. Without reference
eigenvectors or derivatives, the reader retains tensors only when masses are
unchanged and Eckart projection is off; unsupported reprojection raises an
error. Such a tensor-only dataset is not a substitute for a complete mode set.

For central differences, each ``cartesian_displacements`` entry contains
zero-based ``atom`` and ``cart`` (0=x, 1=y, 2=z), a non-zero positive
``delta_angstrom``, and either or both pairs of arrays:

* ``forces_plus_eV_per_angstrom`` and ``forces_minus_eV_per_angstrom``:
  N by 3 forces at the positive and negative displacement. The force constant
  is the negative force difference divided by twice the displacement.
* ``epsilon_plus`` and ``epsilon_minus``: 3 by 3 optical permittivities at
  those displacements. Their difference divided by twice the displacement
  gives ``deps_dr``.

Supply all 3N Cartesian displacements for each response being reconstructed.
Missing displacements and inconsistent array shapes raise errors. Keep all
structures, responses, masses and array orderings consistent; schema recognition
is not a full validation of the physical dataset.

Download the :download:`minimal JSON <_static/inputs/minimal.json>` or the
:download:`complete EO-capable ZnO JSON <_static/inputs/zno-raman.json>`.
The ZnO file contains structure, masses, reference response, ``chi2``, modes,
Raman tensors and Cartesian dielectric derivatives reconstructed from the
CASTEP mode tensors. It is a format-conversion example, not a claim of a new
finite-difference calculation. Its provenance is stored in the file.
``Sphinx/tools/make_format_examples.py`` regenerates these examples and checks
that they can be read. Existing finite-difference workflow data can also be
found in ``Examples/Powder_Raman/Finite_difference`` and
``Examples/Crystal_Raman/Finite_field``.
