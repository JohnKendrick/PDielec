.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package crystal Raman theory
   :keywords: Raman, Crystal, Layer, Transfer Matrix, Scattering Matrix, Reciprocity

.. _Crystal-Raman-Theory:

================================
Theory for Crystal Layer Raman
================================

The *Crystal Raman* calculation models Raman scattering from a crystal slab or a multilayer
stack.  It combines the Raman tensors and phonon frequencies from a lattice-dynamical
calculation with the optical electric fields calculated by the transfer matrix or scattering
matrix methods described in :ref:`Single-Crystal-Theory`.

The method is intended for non-resonant Raman scattering from homogeneous Raman-active
layers.  The incident laser field and the reciprocal scattered field may be modified by
refraction, reflection, absorption and interference within the stack.  These effects are
included by evaluating the electric fields throughout each Raman-active layer and
integrating the Raman source over the layer thickness.

Geometry
--------

The crystal layer uses the same laboratory frame as the single-crystal infrared calculation.
The layer surface lies in the :math:`XY` plane and the surface normal is the laboratory
:math:`Z` direction.  The incident radiation lies in the :math:`XZ` plane; p-polarised light
is polarised in this plane and s-polarised light is polarised along :math:`Y`.

.. _fig-crystal-raman-layers:

.. figure:: ./_static/Figures/Raman_Multi_Layers.png
   :scale: 80%

   Raman scattering from a multilayer system.

The laser is incident with angle :math:`\theta_L`.  The detector may be on the superstrate
side for backscattering or on the substrate side for forward scattering.  The detected
scattered radiation has frequency

.. math::
   :label: eq-crystal-raman-scattered-frequency

   \nu_S = \nu_L - \nu_m

for Stokes scattering from mode :math:`m`.

Layer Raman Amplitude
---------------------

At a point :math:`z_s` inside a Raman-active layer the scattering amplitude of mode
:math:`m` is

.. math::
   :label: eq-crystal-raman-local-amplitude

   A_m(z_s) \propto
   \mathbf{E}_S(z_s,\nu_S)^T
   \tensorbf{R}^{(m)}_{lab}
   \mathbf{E}_L(z_s,\nu_L)

where :math:`\mathbf{E}_L` is the local electric field produced by the incident laser,
:math:`\mathbf{E}_S` is the reciprocal field associated with the detected scattered channel,
and :math:`\tensorbf{R}^{(m)}_{lab}` is the Raman tensor rotated into the laboratory frame.
The transpose on :math:`\mathbf{E}_S` is not a Hermitian conjugate; this is the usual
reciprocity overlap used for Raman scattering from layered optical systems.

For a homogeneous layer :math:`\ell`, the layer amplitude is obtained by integrating through
the layer thickness,

.. math::
   :label: eq-crystal-raman-layer-amplitude

   A_{\ell m} =
   \int_{z_\ell}^{z_{\ell+1}}
   \mathbf{E}_S(z,\nu_S)^T
   \tensorbf{R}^{(m)}_{\ell,lab}
   \mathbf{E}_L(z,\nu_L)\,dz

PDielec evaluates this integral numerically using Gauss-Legendre quadrature points inside
each Raman-active layer.

Optical Fields
--------------

The fields :math:`\mathbf{E}_L` and :math:`\mathbf{E}_S` are calculated from the same
transfer matrix or scattering matrix machinery used for single-crystal infrared spectra.
At optical laser frequencies the phonon contribution to the permittivity is not included;
the high-frequency optical permittivity of each material is used for propagation through
the stack.

The incident field :math:`\mathbf{E}_L(z,\nu_L)` is calculated in the original multilayer
system at the laser frequency.  The scattered field is calculated at the shifted frequency
:math:`\nu_S`.  For backscattering the reciprocal field is launched from the superstrate
side.  For forward scattering it is launched from the substrate side, which is implemented
by using the reversed stack and mapping each integration point to the corresponding
position in that reversed system.

The detected channel can be p-polarised, s-polarised, or unpolarised.  In the unpolarised
case the p- and s-detected intensities are summed incoherently.

Raman Tensor Orientation
------------------------

The Raman tensor from the quantum mechanical calculation is defined in the crystal
coordinate frame.  If :math:`\tensorbf{G}` maps vectors from the crystal frame to the
laboratory frame,

.. math::
   :label: eq-crystal-raman-vector-rotation

   \mathbf{v}_{lab} = \tensorbf{G}\mathbf{v}_{crystal}

then the laboratory Raman tensor is

.. math::
   :label: eq-crystal-raman-tensor-rotation

   \tensorbf{R}^{(m)}_{lab} =
   \tensorbf{G}\tensorbf{R}^{(m)}_{crystal}\tensorbf{G}^{T}

The rotation includes the surface orientation, specified by the layer Miller indices, and the
azimuthal rotation about the laboratory :math:`Z` axis.  This allows the Raman intensity to
be calculated as a function of azimuthal angle.

Combining Layers and Modes
--------------------------

For one mode, the layer amplitudes may be combined coherently or incoherently.  Coherent
combination sums amplitudes before squaring,

.. math::
   :label: eq-crystal-raman-coherent

   I_m \propto \left|\sum_\ell A_{\ell m}\right|^2

whereas incoherent combination sums the layer intensities,

.. math::
   :label: eq-crystal-raman-incoherent

   I_m \propto \sum_\ell \left|A_{\ell m}\right|^2

This layer summation option is distinct from the *Coherent* or *Incoherent*
choice in the layer table.  The layer-table option controls optical propagation
through that layer: whether the transfer/scattering matrix calculation retains
phase interference from internal reflections, averages over phase, or treats a
thick layer approximately.  It therefore changes the optical fields
:math:`\mathbf{E}_L` and :math:`\mathbf{E}_S` used in the Raman overlap.

The *Coherent layer summation* option controls only how Raman amplitudes from
different Raman-active layers are combined after those optical fields have been
calculated.  Leave it unchecked when separate active layers should contribute
independent intensities, and enable it only when Raman emission from those layers
is expected to retain a fixed phase relationship and interfere.

The *Raman depth integration* option is a third, more local coherence choice.  It
controls how Raman sources are combined through the thickness of each individual
Raman-active layer.  *Coherent amplitude* integrates the complex source amplitude
before squaring, which is appropriate for thin coherent films.  *Incoherent
intensity* integrates the local intensity, which is more appropriate for thick or
bulk samples where long-range Raman phase coherence is not physical.

The mode intensity is then multiplied by the Stokes thermal factor,

.. math::
   :label: eq-crystal-raman-bose-factor

   \frac{n(\nu_m)+1}{\nu_m}

where :math:`n(\nu_m)` is the Bose-Einstein occupation factor.  The final spectrum is
formed by applying Lorentzian broadening to each active mode,

.. math::
   :label: eq-crystal-raman-spectrum

   I(\Delta\nu) =
   \sum_m
   I_m
   \frac{\sigma_m}{(\Delta\nu-\nu_m)^2+\sigma_m^2}

where :math:`\sigma_m` is the Lorentzian half-width at half maximum.

Layer Phonon Frequencies
------------------------

The optical field calculation and the active-layer vibrational calculation are separate
parts of the model.  The transfer or scattering matrix calculation determines the laser
field and the reciprocal scattered field at optical frequencies.  It does not, by itself,
change the phonon frequencies.  Any shift of an infrared-active Raman mode away from the
bulk transverse-optic value must therefore be included in the phonon dynamical matrix used
for the Raman-active layer.

For each Raman-active layer PDielec starts from the transverse-optic dynamical matrix,
the mass-weighted Born charge tensor and the Raman tensor of the layer material.  The
phonon-induced polarisation is written as

.. math::
   :label: eq-crystal-raman-phonon-polarisation

   \mathbf{P}_{ph} = \frac{1}{V}\tensorbf{Z}^{mw}\mathbf{x}

where :math:`\mathbf{x}` is the mass-weighted normal-coordinate displacement vector and
:math:`V` is the unit-cell volume used for the Born charges.  A macroscopic electrostatic
field generated by this polarisation gives an additional restoring-force term,

.. math::
   :label: eq-crystal-raman-layer-dynamical

   \tensorbf{D}^{eff} =
   \tensorbf{D}^{TO} +
   \frac{1}{\epsilon_0 V}
   \left(\tensorbf{Z}^{mw}\right)^T
   \tensorbf{S}_{ph}
   \tensorbf{Z}^{mw}

where :math:`\tensorbf{S}_{ph}` is the phonon-field screening tensor.  Different choices
of :math:`\tensorbf{S}_{ph}` define the phonon-frequency model for the active layer.
After the correction has been applied, :math:`\tensorbf{D}^{eff}` is diagonalised and the
resulting frequencies and eigenvectors are used to form the Raman spectrum.  The Raman
tensors are transformed consistently with the corrected eigenvectors when mode mixing is
introduced by the macroscopic field correction.

PDielec provides three levels of treatment for the active-layer phonon frequencies.
These levels should be regarded as separate approximations for the phonon problem, not as
alternatives for solving the optical propagation problem.

**Level 1: bulk TO frequencies**
   No macroscopic phonon-field correction is applied:

   .. math::
      :label: eq-crystal-raman-level1-screening

      \tensorbf{S}_{ph}=\tensorbf{0},
      \qquad
      \tensorbf{D}^{eff}=\tensorbf{D}^{TO}

   This is appropriate for modes that are not infrared active, for calculations where the
   long-range electric field is deliberately neglected, or as a reference calculation.  It
   is also the most direct comparison with a conventional Placzek Raman calculation using
   the bulk transverse-optic phonon frequencies.

**Level 2: travelling-wave NAC correction**
   The phonon is assigned the momentum carried by the Raman process.  In the laboratory
   frame,

   .. math::
      :label: eq-crystal-raman-qph

      \mathbf{q}_{ph}=\mathbf{k}_{L}-\mathbf{k}_{S}

   where :math:`\mathbf{k}_{L}` is the optical wave vector of the laser field in the
   active layer and :math:`\mathbf{k}_{S}` is the wave vector of the emitted Stokes photon
   in the same layer.  The direction
   :math:`\hat{\mathbf{q}}_{ph}=\mathbf{q}_{ph}/|\mathbf{q}_{ph}|` is used in the usual
   non-analytic correction.  After rotating the layer background permittivity into the
   laboratory frame, the screening tensor is

   .. math::
      :label: eq-crystal-raman-nac-screening

      \tensorbf{S}^{NAC}_{ph} =
      \frac{\hat{\mathbf{q}}_{ph}\hat{\mathbf{q}}_{ph}^{T}}
           {\hat{\mathbf{q}}_{ph}^{T}
            \tensorbs{\varepsilon}^{b}_{i,lab}
            \hat{\mathbf{q}}_{ph}}

   This is the bulk non-analytic correction evaluated for the Raman scattering geometry.
   It is useful when the phonon is treated as a propagating long-wavelength bulk mode in
   the active material.  In a multilayer calculation the optical fields may contain several
   forward and backward components.  The implementation therefore uses the wave-vector
   branch associated with the selected incident and detected optical channels; if more than
   one branch is retained, the corresponding Raman contributions are evaluated and combined
   using the same coherent or incoherent choices as the optical Raman amplitudes.

**Level 3: slab electrostatic correction**
   The active layer is treated as a laterally infinite slab for the phonon electrostatic
   problem.  The relevant macroscopic phonon field is then the depolarisation field created
   by the surface-normal component of the phonon-induced polarisation.  The slab
   depolarisation tensor is

   .. math::
      :label: eq-crystal-raman-slab-depolarisation

      \tensorbf{L}_{slab} = \hat{\mathbf{n}}\hat{\mathbf{n}}^T

   where :math:`\hat{\mathbf{n}}` is the surface normal in the laboratory frame.  The
   simplest isolated-slab form uses the active layer's own background permittivity,

   .. math::
      :label: eq-crystal-raman-slab-internal-screening

      \tensorbf{S}^{slab}_{ph} =
      \frac{\hat{\mathbf{n}}\hat{\mathbf{n}}^{T}}
           {\hat{\mathbf{n}}^{T}
            \tensorbs{\varepsilon}^{b}_{i,lab}
            \hat{\mathbf{n}}}

   which is equivalent to the NAC expression with
   :math:`\hat{\mathbf{q}}_{ph}` replaced by :math:`\hat{\mathbf{n}}`.
   When the active layer is embedded between other media, PDielec may instead include the
   dielectric screening of the adjacent layers.  A scalar external background permittivity
   is formed from the normal components of the media above and below the active layer,

   .. math::
      :label: eq-crystal-raman-slab-env-epsilon

      \epsilon^b_e = \tfrac{1}{2}\hat{\mathbf{n}}^T
      \left(\tensorbs{\varepsilon}^{b}_{above,lab}+
            \tensorbs{\varepsilon}^{b}_{below,lab}\right)
      \hat{\mathbf{n}}

   and the internal phonon-field factor is

   .. math::
      :label: eq-crystal-raman-slab-env-nbg

      \tensorbf{N}^{slab}_{ph} =
      \left[\tensorbf{I}+
      \frac{1}{\epsilon^b_e}\tensorbf{L}_{slab}
      \left(\tensorbs{\varepsilon}^{b}_{i,lab}-
      \epsilon^b_e\tensorbf{I}\right)\right]^{-1}

   giving

   .. math::
      :label: eq-crystal-raman-slab-env-sbg

      \tensorbf{S}^{slab-env}_{ph} =
      \frac{1}{\epsilon^b_e}\tensorbf{N}^{slab}_{ph}\tensorbf{L}_{slab}

   This level is intended for phonon modes whose macroscopic field is controlled mainly by
   the planar boundaries of the active layer, rather than by the photon momentum transfer.
   It is therefore the natural correction for polar modes in a thin slab or in a strongly
   dielectric multilayer environment.

For an infrared-inactive mode :math:`\tensorbf{Z}^{mw}` gives no macroscopic restoring-force
correction, so all three levels reduce to the same transverse-optic frequency.  For polar
modes the corrected frequencies can depend on the surface orientation, the optical
scattering geometry and, in the slab-environment model, the background permittivities of
the neighbouring layers.  The optical field enhancement and interference factors are still
calculated separately at :math:`\nu_L` and :math:`\nu_S`; they should not be interpreted as
replacing the phonon-frequency correction described here.
