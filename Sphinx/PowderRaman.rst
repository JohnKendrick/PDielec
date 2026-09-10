.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package powder Raman theory
   :keywords: Raman, Powder, Placzek, Local Field, Depolarisation, Effective Raman Tensor

.. _Powder-Raman-Theory:

========================
Theory for Powder Raman
========================

The *Powder Raman* calculation predicts non-resonant Raman spectra from Raman tensors
calculated by a solid-state quantum mechanical program.  The calculation follows the
Placzek approximation, in which the intensity of a first-order Stokes line is determined by
the derivative of the polarisability with respect to the normal-mode coordinate.  For a mode
:math:`m`, incident laser frequency :math:`\nu_L`, scattered frequency
:math:`\nu_S = \nu_L - \nu_m`, incident polarisation :math:`\mathbf{e}_L`, and detected
polarisation :math:`\mathbf{e}_S`, the scattering strength is proportional to

.. math::
   :label: eq-powder-raman-strength

   S_m =
   \frac{\nu_S^4}{c^4}
   \frac{n(\nu_m)+1}{8\pi^2\nu_m}
   \left|\mathbf{e}_S^{T}\tensorbf{R}^{(m)}\mathbf{e}_L\right|^2

where :math:`n(\nu_m)` is the Bose-Einstein occupation factor.  PDielec reports Raman
intensities in arbitrary units, so the most important outputs are the relative mode
strengths and the changes introduced by particle shape, local fields, polarisation geometry
and linewidth.

Raman Tensor
------------

Every Raman output reader returns the bulk dielectric-derivative tensor

.. math::
   :label: eq-powder-raman-tensor

   \tensorbf{R}^{(m)}_\epsilon = \sqrt{V_{cell}}
   \frac{\partial \tensorbs{\varepsilon}_\infty}{\partial Q_m}

where :math:`V_{cell}` is in Angstrom cubed and :math:`Q_m` is in
Angstrom times the square root of amu.  The tensor units are
:math:`(\mathrm{Angstrom}/\mathrm{amu})^{1/2}`.  The reader already includes
:math:`\sqrt{V_{cell}}`; the powder and crystal calculations apply no additional
cell-volume factor.  A source polarizability-volume derivative
:math:`R_\alpha=\partial[V_{cell}(\varepsilon-I)/(4\pi)]/\partial Q_m`
is converted on read using :math:`R_\epsilon=4\pi R_\alpha/\sqrt{V_{cell}}`.
Activities in :math:`\mathrm{Angstrom}^4/\mathrm{amu}` are display quantities,
obtained by multiplying the internal activities by :math:`V_{cell}/(16\pi^2)`.

For particles whose Raman tensor does not depend on shape, the powder average can be
written in terms of rotational invariants of :math:`\tensorbf{R}`.  For one mode,

.. math::
   :label: eq-powder-raman-invariants

   \alpha &= \frac{1}{3}\left(R_{xx}+R_{yy}+R_{zz}\right) \\
   \tensorbf{\gamma} &= \frac{1}{2}\left(\tensorbf{R}+\tensorbf{R}^{T}\right) -
                      \alpha\tensorbf{1} \\
   \tensorbf{\kappa} &= \frac{1}{2}\left(\tensorbf{R}-\tensorbf{R}^{T}\right)

and

.. math::
   :label: eq-powder-raman-invariant-squares

   \alpha^2 &= \alpha\alpha^* \\
   \gamma^2 &= \frac{3}{2}\sum_{ij}\gamma_{ij}\gamma_{ij}^* \\
   \kappa^2 &= \frac{3}{2}\sum_{ij}\kappa_{ij}\kappa_{ij}^*

The usual parallel (VV), crossed (VH), and unpolarised powder strengths are then

.. math::
   :label: eq-powder-raman-vv-vh

   I_{VV}     &\propto 45\alpha^2 + 4\gamma^2 \\
   I_{VH}     &\propto 3\gamma^2 + 5\kappa^2 \\
   I_{total}  &\propto 45\alpha^2 + 7\gamma^2 + 5\kappa^2

For non-chiral, non-magnetic Raman scattering the antisymmetric term is normally zero.

Local Field Correction
----------------------

For powder infrared calculations PDielec treats each crystallite as a small inclusion in a
non-absorbing matrix.  The powder Raman calculation uses the same macroscopic picture
when local-field effects are requested.  The particle is assumed to be much smaller than the
laser wavelength, so the field inside the particle is uniform.  For an inclusion with optical
permittivity :math:`\tensorbs{\varepsilon}_i`, a matrix permittivity
:math:`\varepsilon_e`, and a depolarisation tensor :math:`\tensorbf{L}`, the internal field is

.. math::
   :label: eq-powder-raman-internal-field

   \mathbf{E}_i = \tensorbf{N}\mathbf{E}_e

with

.. math::
   :label: eq-powder-raman-N

   \tensorbf{N} =
   \left[
   \tensorbf{1} +
   \frac{1}{\varepsilon_e}
   \tensorbf{L}
   \left(\tensorbs{\varepsilon}_i-\varepsilon_e\tensorbf{1}\right)
   \right]^{-1}

For a sphere :math:`\tensorbf{L}` has diagonal elements of :math:`1/3`.  Other ellipsoidal
shapes use the same depolarisation tensor as the powder infrared effective-medium theory.

For a reciprocal dielectric, the effective bulk Raman tensor includes one
incident and one scattered optical local-field factor:

.. math::
   :label: eq-powder-raman-particle-tensor

   \tensorbf{R}^{(m)}_{eff} =
   \tensorbf{N}(\nu_S)^T
   \left[\tensorbf{R}^{(m)}_{\epsilon,0}+\Delta\tensorbf{R}^{(m)}_{\epsilon,EO}\right]
   \tensorbf{N}(\nu_L)

The transpose is ordinary, including for complex reciprocal permittivities.
The current powder implementation uses the same optical permittivity at both
frequencies.  :math:`R_{eff}` retains the bulk :math:`R_\epsilon` normalization;
it is not an extensive particle dipole-polarizability derivative.  The former
three-factor expression double-counted the scattered-field response.

Particle Frequencies
--------------------

In ionic crystals the particle boundary conditions can also shift the vibrational
frequencies.  The phonon-induced polarisation is

.. math::
   :label: eq-powder-raman-phonon-polarisation

   \mathbf{P}_{ph} = \frac{1}{V}\tensorbf{Z}^{mw}\mathbf{x}

where :math:`\tensorbf{Z}^{mw}` is the mass-weighted Born charge tensor and
:math:`\mathbf{x}` is the vector of mass-weighted atomic displacements.  The
depolarisation field generated by this polarisation applies a restoring force to the ions.
This gives an effective particle dynamical matrix

.. math::
   :label: eq-powder-raman-particle-dynamical

   \tensorbf{D}^{particle} =
   \tensorbf{D}^{TO} +
   \frac{1}{\epsilon_0\varepsilon_e V}
   \left(\tensorbf{Z}^{mw}\right)^T
   \tensorbf{N}_{bg}\tensorbf{L}\tensorbf{Z}^{mw}

where :math:`\tensorbf{D}^{TO}` is the transverse-optic dynamical matrix and
:math:`\tensorbf{N}_{bg}` is the internal-field tensor constructed from the background
permittivities.  Diagonalising :math:`\tensorbf{D}^{particle}` gives particle-mode
frequencies and eigenvectors appropriate to the selected crystallite shape.

In the SI expression above Born charges carry coulombs and masses carry kilograms.
The implementation uses charges in electrons, masses in electron-mass units,
volume in Bohr cubed, and the prefactor :math:`4\pi/(\varepsilon_e V)`.
The same electrostatic kernel
:math:`K=N_{bg}L/\varepsilon_e` determines the particle EO correction:

.. math::

   \Delta R^{(m)}_{\epsilon,EO,ij}
   =-8\pi\sum_l\widetilde\chi^{(2)}_{ijl}(K Z^{mw}u_m)_l .

Here :math:`\widetilde\chi^{(2)}` is the cell-dependent internal reader value,
not the raw susceptibility in pm/V.  Its conversion is specified in
:ref:`Crystal-Raman-Theory`.  Both the mechanical Raman tensor and the mode
charge must use the same particle eigenvector.  No bulk propagation direction
is needed for this ellipsoid boundary-value model.  It describes a quasistatic
particle in a lossless scalar host; it is not a bulk directional LO/TO powder
average.  ``Matrix=none`` bypasses particle and optical-field corrections;
use a host with permittivity one to describe an isolated particle in vacuum.

Orientation Averaging
---------------------

For a crystallite in orientation :math:`g`, the particle Raman tensor in the laboratory frame
is

.. math::
   :label: eq-powder-raman-rotation

   \tensorbf{R}^{(m)}_{lab}(g) =
   g\tensorbf{R}^{(m)}_{eff}g^T

The powder strength is the orientational average over all rotations,

.. math::
   :label: eq-powder-raman-average

   \left<S_m\right>_{powder} =
   \int_{SO(3)}
   S_m(g)\,dg

For spherical particles this average reduces to the invariant expressions above.  For more
general particle tensors, PDielec samples the rotation group and averages the fixed
laboratory polarisation geometry over many crystallite orientations.

The broadened Raman spectrum is formed as a sum over active modes,

.. math::
   :label: eq-powder-raman-spectrum

   I(\Delta\nu) =
   \sum_m
   \left<S_m\right>_{powder}
   \frac{\sigma_m}{(\Delta\nu-\nu_m)^2+\sigma_m^2}

where :math:`\sigma_m` is the Lorentzian half-width at half maximum.
