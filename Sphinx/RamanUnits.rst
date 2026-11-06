.. include:: preamble.txt

.. _Raman-Units:

Raman units and implementation conventions
==========================================

This is the shared numerical-conversion reference for the powder and crystal
Raman chapters. The physical derivations use SI quantities; numerical Raman
frequencies denoted by :math:`\nu` are wavenumbers in cm\ :sup:`−1`.
Use :math:`f` for ordinary frequency in Hz and :math:`\omega=2\pi f` for
angular frequency. The GTM optical solver takes Hz, whereas Raman interfaces
and linewidths use wavenumbers. Material properties are relative permittivities.

.. list-table:: Common optical and Raman conversion boundaries
   :header-rows: 1
   :widths: 28 36 36

   * - Quantity
     - Stored convention
     - Conversion or consequence
   * - Optical frequency
     - Wavenumber in cm\ :sup:`−1`
     - :math:`f_{Hz}=100c_{SI}\nu`, :math:`\omega=2\pi f_{Hz}`
   * - Laser wavelength
     - Nanometres
     - :math:`\nu_L=10^7/\lambda_{nm}`
   * - Optical depth and phase
     - Depth in metres; :math:`\zeta,q_z` dimensionless
     - :math:`k_0=2\pi100\nu`; propagation phase uses :math:`k_0q_z z`
   * - Raman depth quadrature
     - Weights in Å
     - Metre-valued quadrature weights multiplied by :math:`10^{10}`
   * - Raman tensor
     - :math:`R_\epsilon` in :math:`\sqrt{\mathrm{Å}/\mathrm{amu}}`
     - Multiply by :math:`\sqrt{10^{-10}/m_u}` for :math:`\sqrt{\mathrm m/\mathrm{kg}}`
   * - Raw coherent overlap squared
     - Å cubed per amu
     - Multiply by :math:`10^{-30}/m_u` for cubic metres per kg
   * - Local or normalized-profile weight
     - Å squared per amu
     - Multiply by :math:`10^{-20}/m_u` for square metres per kg
   * - Stokes prefactor
     - :math:`\nu_S^4(n+1)/\nu_m`, all wavenumbers in cm\ :sup:`−1`
     - Applied once; excludes absolute radiation and collection constants
   * - Broadened line
     - :math:`I_m\sigma_m/[(\nu-\nu_m)^2+\sigma_m^2]`
     - Full area :math:`\pi I_m`; width is HWHM

Here :math:`m_u` is the atomic mass unit in kg. These conversions do not equate
raw coherent and local observables or supply missing physical normalization.


Cell, mode and response conversions
-----------------------------------

.. list-table:: Lattice and response quantities
   :header-rows: 1
   :widths: 24 36 40

   * - Quantity
     - Stored convention
     - Boundary or consequence
   * - Cell volume
     - Reader volume in Å cubed
     - :math:`V_{SI}=10^{-30}V`; NAC uses :math:`V_{au}=V\,\mathrm{angs2bohr}^3`
   * - Atomic masses
     - Reader masses in amu
     - NAC uses :math:`M_{au}=M\,\mathrm{Constants.amu}`; SI uses :math:`M m_u`
   * - Born charge
     - Elementary-charge units, crystal Cartesian frame
     - Multiply by :math:`e` for coulombs; mass-weighted matrix divides by :math:`\sqrt{M_{au}}` in the NAC kernel
   * - Mode coordinate
     - :math:`Q` in Å times square root of amu; dimensionless unit-normalised mass-weighted eigenvector :math:`u`
     - :math:`Q_{SI}=10^{-10}\sqrt{m_u}Q`; displacement :math:`r_a=u_a Q/\sqrt{M_a}`
   * - Dynamical matrix
     - Hartree atomic units in NAC calculations
     - Signed wavenumber :math:`\nu=\operatorname{sgn}(\lambda)\sqrt{|\lambda|}/\mathrm{Constants.wavenumber}`; negative values flag instability
   * - Polarizability-volume Raman derivative
     - :math:`R_\alpha` in Å squared per square root of amu
     - :math:`R_\epsilon=4\pi R_\alpha/\sqrt{V}`; this source convention differs from SI dipole polarizability
   * - Activity display
     - Å to the fourth power per amu
     - Multiply the corresponding internal quadratic tensor activity by :math:`V/(16\pi^2)`; this is not a spectrum or cross section
   * - Nonlinear susceptibility
     - Physical input in pm/V, then cell-dependent internal value
     - :math:`\widetilde\chi^{(2)}=\chi^{(2)}[\mathrm{pm/V}]\sqrt{m_u/m_e}\,e10^8/(4\pi\epsilon_0\sqrt{V})`
   * - Local-field factors
     - Dimensionless tensors in the same Cartesian frame as the Raman tensor
     - :math:`R_{eff}=N_S^T R_\epsilon N_L` retains the bulk tensor units; transpose is not conjugate transpose

Here :math:`V` in reader conversions is the numerical volume in Å cubed,
:math:`m_e` and :math:`m_u` are masses in kg, and :math:`e` is in coulombs.
``Constants.amu`` is the mass ratio :math:`m_u/m_e`, not a mass in kg.
``Constants.wavenumber`` is :math:`100hc_{SI}/E_h`. For an SI dynamical
matrix, use :math:`\operatorname{sgn}(\lambda)\sqrt{|\lambda|}/(2\pi100c_{SI})`
instead. Do not insert SI charges or masses into the atomic-unit NAC kernel.
The internal EO contraction uses :math:`-8\pi\widetilde\chi^{(2)}` with
atomic-unit mass-weighted charges; it cannot take raw pm/V input.

Implementation and frames
-------------------------

* ``GenericOutputReader.get_raman_tensors`` returns bulk
  :math:`R_\epsilon=\sqrt{V_{cell}}\partial\epsilon/\partial Q` in crystal
  Cartesian axes. Reader-specific native conversions precede this boundary.
  ``chi2_pm_per_v_to_repsilon`` implements the separate nonlinear conversion.
* ``NACDiagnostics.compute_nac_direction_diagnostics`` accepts atomic-unit
  masses and volume. Its eigenvectors are dimensionless; its reported mode
  charges explicitly convert back to :math:`e/\sqrt{\mathrm{amu}}`.
* ``Calculator.compute_effective_raman_tensor`` applies the two optical
  factors. A crystal-to-laboratory rotation is :math:`G R G^T`; rotating a
  tensor changes its components, not its units or cell normalization.
* ``LayeredRamanCalculator`` uses optical positions in metres and Raman
  integration weights in Å. Its fields are reduced unit-electric-amplitude
  launch fields, not physical V/m inputs. GTM magnetic output is scaled by
  vacuum impedance relative to the same electric reference amplitude.
* ``RamanSpectrum.stokes_prefactor`` supplies the thermal/radiation weight
  once. ``lorentzian_broaden`` broadens the weighted line coefficients.
  The SI unit of a theoretical efficiency must not be attached to these
  uncalibrated reduced spectra.

Worked conversion examples
--------------------------

These examples check conventions with prescribed tensors and fields; they are
not material predictions or new DFT results.

**Reader tensor.** For :math:`V=100` Å cubed and one component
:math:`R_\alpha=2` Å squared per square root of amu,
:math:`R_\epsilon=0.8\pi` in :math:`\sqrt{\mathrm{Å}/\mathrm{amu}}`.
Its square is :math:`0.64\pi^2` Å/amu. Multiplication by
:math:`V/(16\pi^2)` recovers the squared source component, 4 Å to the fourth
power per amu. This scalar check is not the full invariant powder activity.
The SI Raman component is :math:`0.8\pi\sqrt{10^{-10}/m_u}`.

**Powder particle.** Prescribe :math:`R_\epsilon=I` in stored tensor units,
relative inclusion permittivity 4, host permittivity 1, and a sphere.
At both optical frequencies :math:`N=3/(4+2)I=I/2`, so
:math:`R_{eff}=I/4`. The implemented activities are VV = :math:`45/16`
and VH = 0 in Å/amu. The VV activity is 45 times the orientational mean
squared amplitude. For a stable mode at 400 cm\ :sup:`−1`, laser at
20000 cm\ :sup:`−1`, zero temperature, unit volume fraction and HWHM 5
cm\ :sup:`−1`, the reduced peak height is
:math:`(45/16)(19600)^4/(400\times5)`. It is not an absolute cross section.

**Crystal layer.** Let the reduced bilinear overlap be constant,
:math:`e_S^T R_\epsilon e_L=1` in :math:`\sqrt{\mathrm{Å}/\mathrm{amu}}`,
over a 100 nm layer. Its optical thickness is :math:`10^{-7}` m; its Raman
quadrature weights sum to 1000 Å. The raw coherent squared overlap is
:math:`10^6` Å cubed/amu, whereas the local-depth weight is 1000 Å squared/amu.
A single normalized uniform profile gives the latter weight for this constant
overlap. Neither expression includes the Stokes prefactor. Doubling thickness
multiplies the raw coherent weight by four and the local weight by two.

**Spectral axis.** With :math:`a=2\pi100c_{SI}` and
:math:`\omega=a\nu`, a spectral density obeys
:math:`I_\omega(\omega)=I_\nu(\omega/a)/a`. Merely relabelling the axis
or multiplying the line centre and width by :math:`a` does not preserve area.
The numerical Lorentzian area is :math:`\pi` times its supplied coefficient.

Validation scope
----------------

The conversion reference records implementation conventions and the bounded
checks described in the Raman theory chapters. It does not certify all
materials, optical regimes or physical final-state models. Broader NAC material
validation and general absorbing/anisotropic multilayer validation will be
carried out with the package paper's examples. Each example must establish the
assumptions it relies on before its spectrum is interpreted. Existing datasets
and controlled limits should be used first; additional DFT calculations require
an identified missing input or validation case. General extensions outside the
chosen examples remain separate work.
