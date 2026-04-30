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

For a particle of volume :math:`V`, the dipole polarisability tensor
:math:`\tensorbf{\alpha}` relates the electric field to the induced dipole moment.  The Raman
tensor of mode :math:`m` is

.. math::
   :label: eq-powder-raman-tensor

   \tensorbf{R}^{(m)} =
   \frac{\partial \tensorbf{\alpha}}{\partial Q_m}
   =
   \epsilon_0 V
   \frac{\partial \tensorbs{\chi}}{\partial Q_m}
   =
   \epsilon_0 V
   \frac{\partial \tensorbs{\varepsilon}}{\partial Q_m}

where :math:`Q_m` is the mass-weighted normal-mode coordinate.  The final equality follows
because the electric susceptibility and relative permittivity differ only by the identity
tensor.  Output readers may store a volume-normalised tensor; where required PDielec
rescales it before calculating activities.

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

   I_{VV}     &\propto 45\alpha^2 + 4\gamma^2 + 5\kappa^2 \\
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

The particle Raman tensor connects the external incident field to the external scattered
field.  Assuming that the optical permittivities are constant across the laser and scattered
frequencies, the effective particle tensor is

.. math::
   :label: eq-powder-raman-particle-tensor

   \tensorbf{R}^{(m)}_{particle} =
   \tensorbf{N}
   \left[
   \tensorbf{R}^{(m)}_{\varepsilon}
   -
   \frac{1}{\varepsilon_e}
   \left(\tensorbs{\varepsilon}_i-\varepsilon_e\tensorbf{1}\right)
   \tensorbf{N}\tensorbf{L}\tensorbf{R}^{(m)}_{\varepsilon}
   \right]
   \tensorbf{N}

where :math:`\tensorbf{R}^{(m)}_{\varepsilon}` is the Raman tensor derived from
:math:`\partial\tensorbs{\varepsilon}/\partial Q_m`.  The first and last
:math:`\tensorbf{N}` factors account for the field entering and leaving the particle.  The
middle correction accounts for the dependence of the particle polarisability on dielectric
contrast.

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

Orientation Averaging
---------------------

For a crystallite in orientation :math:`g`, the particle Raman tensor in the laboratory frame
is

.. math::
   :label: eq-powder-raman-rotation

   \tensorbf{R}^{(m)}_{lab}(g) =
   g\tensorbf{R}^{(m)}_{particle}g^T

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
