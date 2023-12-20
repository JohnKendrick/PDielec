.. include:: preamble.txt

..
    .. contents::
       :local:
..


.. meta::
   :description: PDielec package for the calculation of infrared and terahertz absorption from QM calculations
   :keywords: Quantum Mechanics, Effective Field Theory, Maxwell, Garnett, Mie, Infrared, Terahertz, Castep, Abinit, VASP, GULP. QE


===============
Phonon Analysis
===============

PDGui performs an analysis of a phonon mode in terms of the percentage 
of molecular centre-of-mass and rigid body rotational motion in each mode.  Using the results of the analysis it is possible to differentiate between inter- and intra-molecular phonon modes (see for instance Jepsen et al :cite:`Jepsen2007`).  In addition, for those systems with more than one formula unit in the unit-cell it is possible to break down the contribution of each molecule to a particular phonon mode.

Molecular Systems
-----------------

It is common in molecular calculations of
vibrational properties to construct a force constant matrix which enforces the
requirements of zero energy change for centre-of-mass motion and rigid-body rotation. This can be achieved by using projection operators to transform
the second derivative matrix to a set of coordinates that no longer include centre-of-mass motion or rigid-body rotation.

Defining a projection operator as;

.. math::
   :label: eq-projectionop

   \tensorbf{P} = \tensorbf{1} - \fieldbf{V} \fieldbf{V}^{T}


where :math:`\fieldbf{V}` is an orthonormal column vector with length :math:`3N` (:math:`N` is the number of atoms), the centre-of-mass motion
can be described in Cartesian space as all atoms moving along the x, y or
z axis with the same displacement. So considering the projection of the centre-of-mass motion along the x-axis, for each atom :math:`a` we can write;


.. math::
   :label: eq-cmx1

   \fieldbf{V}^{x}_{a} = \left(\begin{matrix} 1 \\ 0 \\ 0 \end{matrix}\right) 


For the projection operator which will project out all components of translation along the x-axis for every atom in the molecule we have;

.. math::
   :label: eq-cmx2

   \fieldbf{V}^{x} = \begin{pmatrix} \fieldbf{V}^{x}_{1} \\ \fieldbf{V}^{x}_{2}\\ \vdots \end{pmatrix} 


There are two other projection operators describing translation along the y- and z-axis.  

In a similar fashion, it is possible to describe an infinitesimal molecular rotation using a vector V and therefore constructing a projection operator to remove rigid-body rotation.  If the coordinates of atom :math:`a` relative to the centre-of-mass of the molecule are :math:`x, y \text{ and } z`, the component of the projection vector, :math:`\tensorbf{V}` representing rotation about the x-axis in the yz-plane is;


.. math::
   :label: eq-rotyz

    \fieldbf{V}^{yz}_{a} = \left(\begin{matrix} 0 \\ -z \\ +y \end{matrix}\right) 

In a similar fashion rotations about the y- and z-axis are respectively;


.. math::
   :label: eq-rotxz

   \fieldbf{V}^{xz}_{a} = \left(\begin{matrix} +z \\ 0 \\ -x \end{matrix}\right)


and

.. math::
   :label: eq-rotxy

   \fieldbf{V}^{xy}_{a} = \left(\begin{matrix} -y \\ +x \\ 0 \end{matrix}\right)


In a similar fashion to Equation :eq:`eq-cmx2` the complete projection operators for the 3 rotational modes can be assembled (using rotation about the x-axis as an example)

.. math::
   :label: eq-rotyz2

    \fieldbf{V}^{yz} = \begin{pmatrix} \fieldbf{V}^{yz}_{1} \\ \fieldbf{V}^{yz}_{2}\\ \vdots \end{pmatrix} 

In practice, the projection operators will be defined using mass-weighted cartesian coordinates as this simplifies the expressions used later.  As an example the transformation to mass-weighted coordinates for atom :math:`a` translating along the x-axis is;


.. math::
   :label: eq-mass-weight

   \sqrt{m_a} \fieldbf{V}^{x}_{a}


  
Periodic Systems
----------------

The invariants of the energy in a periodic system are the three
translational modes corresponding to motion of all the atoms in the same
direction.  The rotation of an infinite lattice is not invariant.  However,  it is proposed to construct projection operators for each
molecule in the unit-cell, which will project out the motion of
each molecular translational and rotational degree of freedom.
Such projection operators will be used to separate the external from the internal modes and there will therefore
be a projection operator for each external mode in the system.

Hug and Haesler have shown :cite:`Hug2005c` that the vibrational kinetic energy can be decomposed into single centre atomic contributions by considering the kinetic energy, :math:`T_p`, of a normal mode :math:`Q_p`;


.. math::
   :label: eq-kep

   T_p = \frac{1}{2} \dot{Q_p} = 
   \sum\limits_{a}{T_{a,p}} = 
   \frac{1}{2}\sum\limits_{a}{m_a \fieldbf{\Delta\dot{x}}^T_{a,p}\fieldbf{\Delta \dot{x}}_{a,p} } =  
   \frac{1}{2} \dot{Q^2_p} \sum\limits_{a} { \fieldbf{L}^T_{a,p} \fieldbf{L}_{a,p} }


Where the sum is over atoms, :math:`a`, with mass :math:`m_a`, :math:`\fieldbf{\Delta \dot{x}}_{a,p}` is the time
derivative of the Cartesian displacement vector of atom :math:`a` along the :math:`p^{th}` normal mode and :math:`\fieldbf{L}_p` are the
components of the orthogonal transformation relating the :math:`p^{th}` normal
mode and mass-weighted Cartesian displacements. Thus :math:`\fieldbf{L}` are the eigenvectors of the 
mass-weighted second derivative matrix.  The atom subscript is used to indicate that
only those atomic coordinates involving atom :math:`a` are being considered.
Derivation of this expression has made use of the following
relationship between the mass-weighted normal mode;

.. math::
   :label: eq-kinetic-energy

   \sqrt{ m_a } \fieldbf{\Delta x}_{a,p} = \fieldbf{L}_{a,p} Q_p

As proposed by Hug and Haesler :cite:`Hug2005c`, since the kinetic energy can be expressed as
atomic contributions, consideration of the Virial theorem indicates that
this is true also for the potential energy. Using mass-weighted
Cartesian coordinates for convenience, an analysis of the phonon modes
may therefore be constructed in the following way. A given phonon
mode, :math:`p`,  will have its total kinetic energy partitioned between the
molecules in the unit-cell according to;

.. math::
  :label: eq-kinetic-energy_total

   E^{total}_p = \sum\limits_{mol} {E^{mol}_p } 

.. math::
   :label: eq-kinetic-energy_molecular

   E^{mol}_p = \dot{Q^2_p} \sum\limits_{a \in mol} { \fieldbf{L}^T_{a,p} \fieldbf{L}_{a,p} }


To calculate the contribution a particular mode has to the centre-of-mass kinetic energy we use the projection operators given in Equation :eq:`eq-cmx2` (now in mass-weighted Cartesian coordinates) where only the atoms in that molecule are used to construct the operator;

.. math::
   :label: eq-finaloperator

   \tensorbf{P}^x_{mol} = \tensorbf{1} - \fieldbf{V}^x_{mol} (\fieldbf{V}^x_{mol})^T


This approach is used to define  6 projection operators for each molecule. The centre-of-mass energy contained in the :math:`p_{th}` mode can be written as;

.. math::
   :label: eq-energy-cm

   E^{cm}_p = \dot{Q^2_p} \sum\limits_{i=x,y,z}{\sum\limits_{mol} { (\tensorbf{P}_{mol}^i\fieldbf{L}_{p})^T (\tensorbf{P}_{mol}^i \fieldbf{L}_{p} ) }}


and the rigid-body rotational energy associated with the mode is;

.. math::
   :label: eq-energy-rot

   E^{rot}_p = \dot{Q^2_p} \sum\limits_{i=xy,yz,xz}{\sum\limits_{mol} { (\tensorbf{P}_{mol}^i\fieldbf{L}_{p})^T (\tensorbf{P}_{mol}^i \fieldbf{L}_{p} ) }}


The molecular vibrational contribution, which can be used to classify the internal modes of the system, can be obtained by subtracting the centre-of-mass and rigid-body rotational; energies from the total;

.. math::
   :label: eq-energy-vib

   E^{vib}_p = E^{total}_p - E^{cm}_p - E^{rot}_p


Since all of the energy terms depend in the same way on :math:`\dot{Q}_p`, it is not required in the calculation of the relative contributions to the energy coming from the external (molecular centre-of-mass and rigid-body rotation) and the internal modes (vibrational contributions)

Implementation
--------------

In order to calculate the relative internal and external contributions to each phonon mode it is first necessary to identify the molecules in the crystal.  In many cases, the unit-cell is packed with atoms in such a way that the cell is filled, rather than in a way reflecting the bondedness of the molecules.  The program first replicates the atoms in all cells neighbouring the central unit-cell.  Within this supercell the bonds are determined by calculating the distances between all atoms in the supercell.  In practice, an order N method is used whereby instead of searching all the supercell only the space around each atom is searched for potentially bonded partners.  The criterion of the presence of a bond between atoms :math:`i` and :math:`j` is given by the requirement that the distance between the atoms :math:`r{ij}` is less than the bonding requirement;


.. math::
   :label: eq-bond

   r_{ij} < scale (radius_i+radius_j)+toler


Here :math:`radius_i` is the covalent radius of atom :math:`i`, :math:`scale` and :math:`toler` are factors that can be altered when running the program.

Once the bonding in the supercell has been determined the program starts with the first atom in the central cell and determines all atoms which are connected to the molecule that it is in.  If there are any remaining atoms in the central cell which are not bonded yet then further molecules are added until all atoms in the central cell have been allocated to a molecule.

Finally a new cell with the same dimensions as the original is constructed.  Where necessary each molecule is shifted into the cell so that the molecule's centre-of-mass lies within it.

The projection operators, Equation :eq:`eq-finaloperator`, can now be constructed and the relative energies calculated using Equations :eq:`eq-energy-cm`, :eq:`eq-energy-rot` and :eq:`eq-energy-vib`. 

Examples
========

Isoleucine
----------

The example of a calculation on isoleucine using Castep is available in the Examples/Castep/Isoleucine directory of the distribution.  The analysis was performed using PDGui.
The program first finds the four molecules of isoleucine in the unit-cell of the crystal structure.  The ordering of the atoms at this point is different to the ordering in the original cell and in the QM calculation.  The visualiser shows the atoms in the new positions.  These will not be the same as in the original calculation.  The molecules and their centres of mass are summarised in :numref:`tab-isoleucine-mols`.

.. table:: Molecular constituents of the isoleucine unit-cell (centre-of-mass is in fractional coordinates)
   :name: tab-isoleucine-mols
   :column-dividers:  single single single single single single single 
   :widths:                  1      1      1      1      1      
   :column-alignment:        center center center center center
   :header-alignment:        center center center center center
   :align: center

   +----------+---------+--------------------+--------------------+---------------------+
   | Molecule | Mass    | CM\ :subscript:`a` | CM\ :subscript:`b` | CM\ :subscript:`c`  |
   +==========+=========+====================+====================+=====================+
   | 0        | 131.172 | 0.293600           | 0.051683           | 0.304906            |
   +----------+---------+--------------------+--------------------+---------------------+
   | 1        | 131.172 | 0.706400           | 0.551683           | 0.695094            |
   +----------+---------+--------------------+--------------------+---------------------+
   | 2        | 131.172 | 0.793231           | 0.588919           | 0.287125            |
   +----------+---------+--------------------+--------------------+---------------------+
   | 3        | 131.172 | 0.206769           | 0.088919           | 0.71287 5           |
   +----------+---------+--------------------+--------------------+---------------------+


:numref:`tab-isoleucine-results` shows the results of the analysis of the phonon modes with frequencies below 100 |cm-1|.  The first 3 modes should have zero frequency as they are the translationally invariant modes of the lattice.  The negative number shown here actually indicates that the mode has an imaginary frequency.  As can be seen under the column %mol-cme they are almost completely associated with molecular centre-of-mass motion.  This is true also of the lowest non-zero frequency at 32.81 |cm-1| and to some extent of the next mode at 39.73 |cm-1|.  
However, for this latter mode there is a considerable contribution from rigid-body rotational motion (see the %mol-rot column).  As the frequency of the phonon mode increases the contribution from rigid-body motion generally decreases and the contribution from vibrational modes (see the %vib column) increases.  The total contribution from molecular motion is summarised in the last four columns. 

.. table:: Isoleucine: percentage contributions of the centre-of-mass and the rotational molecular modes to each phonon mode
   :name: tab-isoleucine-results
   :column-dividers:  single single single single single single single single single single
   :column-alignment:        right  right  right  right  right  right  right  right
   :header-alignment:        center center center center center center center center
   :widths:                  1      1      1      1      1      1      1      1    
   :align: center

   +--------------+----------+----------+------+--------+--------+--------+--------+
   | Freq(|cm-1|) | %cm      | %rot     | %vib | %mol-0 | %mol-1 | %mol-2 | %mol-3 |
   +==============+==========+==========+======+========+========+========+========+
   | -0.07        | 99.5     | 0.5      | 0.0  | 26.9   | 26.9   | 23.1   | 23.1   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | -0.06        | 98.6     | 1.1      | 0.3  | 23.6   | 23.6   | 26.4   | 26.4   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | -0.04        | 99.7     | 0.2      | 0.1  | 23.8   | 23.8   | 26.2   | 26.2   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 32.81        | 92.0     | 5.9      | 2.1  | 26.8   | 26.8   | 23.2   | 23.2   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 39.73        | 77.9     | 19.6     | 2.5  | 22.1   | 22.1   | 27.9   | 27.9   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 49.42        | 10.4     | 81.9     | 7.7  | 19.0   | 19.0   | 31.0   | 31.0   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 52.12        | 80.7     | 17.3     | 2.0  | 17.6   | 17.6   | 32.4   | 32.4   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 54.15        | 73.4     | 19.0     | 7.5  | 21.2   | 21.2   | 28.8   | 28.8   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 59.59        | 37.7     | 56.6     | 5.7  | 23.9   | 23.9   | 26.1   | 26.1   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 62.05        | 77.3     | 17.2     | 5.5  | 23.4   | 23.4   | 26.6   | 26.6   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 69.91        | 9.6      | 69.4     | 21.0 | 17.6   | 17.6   | 32.4   | 32.4   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 73.19        | 7.3      | 61.1     | 31.6 | 6.6    | 6.6    | 43.4   | 43.4   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 76.00        | 13.6     | 61.7     | 24.7 | 31.6   | 31.6   | 18.4   | 18.4   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 76.25        | 3.2      | 48.1     | 48.7 | 10.6   | 10.6   | 39.4   | 39.4   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 80.49        | 31.5     | 43.8     | 24.7 | 19.3   | 19.3   | 30.7   | 30.7   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 87.54        | 37.9     | 35.4     | 26.7 | 26.8   | 26.8   | 23.2   | 23.2   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 89.31        | 12.2     | 58.1     | 29.6 | 27.5   | 27.5   | 22.5   | 22.5   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 90.49        | 16.9     | 53.9     | 29.2 | 22.4   | 22.4   | 27.6   | 27.6   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 94.51        | 49.1     | 36.9     | 14.0 | 39.4   | 39.4   | 10.6   | 10.6   |
   +--------------+----------+----------+------+--------+--------+--------+--------+
   | 94.79        | 47.0     | 16.9     | 36.1 | 23.0   | 23.0   | 27.0   | 27.0   |
   +--------------+----------+----------+------+--------+--------+--------+--------+


BaTiO\ :subscript:`3`
---------------------

The example for 
BaTiO\ :subscript:`3`
can be found in Example/AbInit/BaTiO3 in the main directory of the PDielec distribution. The results of the Abinit calculation can be analysed using PDGui.
For this example, the default covalent radius of the Barium atoms has been changed to 0.3 Å, so this atom will be treated as though it is not bonded to anything else in the cell.  This results in 2 'molecules' being found in the cell; a TiO\ :subscript:`3` moiety and the Ba\ :superscript:`2+` ion. 

The results of the analysis are shown in :numref:`tab-batio3-results` below.  In this example, the lowest 3 modes have been projected so there is no problem with translational invariance.  They are purely centre-of-mass modes.  The energy in these three modes seems to be mainly in molecule 1, which is the Ba\ :superscript:`2+` ion.  The modes at 197.73 |cm-1| also have a large centre-of-mass component, this time mainly coming from molecule 0, which is the TiO\ :subscript:`3` moiety. Above 200 |cm-1| there is a much smaller centre-of-mass contribution to the energy and all the modes are dominated by the TiO\ :subscript:`3` group.

.. table:: BaTiO\ :subscript:`3`: percentage contributions of the centre-of-mass and the rotational molecular modes to each phonon mode
   :name: tab-batio3-results
   :column-dividers:  single single single single single single single single 
   :column-alignment:        right  right  right  right  right  right  
   :header-alignment:        center center center center center center 
   :widths:                  1      1      1      1      1      1      
   :align: center

   +--------------+----------+----------+------+--------+--------+
   | Freq(|cm-1|) | %cm      | %rot     | %vib | %mol-0 | %mol-1 |
   +==============+==========+==========+======+========+========+
   | 0.00         | 100.0    | 0.0      | 0.0  | 41.1   | 58.9   |
   +--------------+----------+----------+------+--------+--------+
   | 0.00         | 100.0    | 0.0      | 0.0  | 41.1   | 58.9   |
   +--------------+----------+----------+------+--------+--------+
   | 0.00         | 100.0    | 0.0      | 0.0  | 41.1   | 58.9   |
   +--------------+----------+----------+------+--------+--------+
   | 197.73       | 72.9     | 6.7      | 20.4 | 70.0   | 30.0   |
   +--------------+----------+----------+------+--------+--------+
   | 197.73       | 72.9     | 6.7      | 20.4 | 70.0   |  30.0  |
   +--------------+----------+----------+------+--------+--------+
   | 197.74       | 72.9     | 6.8      | 20.4 | 70.0   | 30.0   |
   +--------------+----------+----------+------+--------+--------+
   | 269.34       | 27.1     | 20.3     | 52.7 | 88.9   | 11.1   |
   +--------------+----------+----------+------+--------+--------+
   | 269.34       | 27.0     | 20.3     | 52.7 | 88.9   | 11.1   |
   +--------------+----------+----------+------+--------+--------+
   | 269.34       | 27.0     | 20.3     | 52.7 | 88.9   | 11.1   |
   +--------------+----------+----------+------+--------+--------+
   | 270.59       | 0.0      | 55.6     | 44.4 | 100.0  | 0.0    |
   +--------------+----------+----------+------+--------+--------+
   | 270.59       | 0.0      | 55.6     | 44.4 | 100.0  | 0.0    |
   +--------------+----------+----------+------+--------+--------+
   | 270.59       | 0.0      | 55.6     | 44.4 | 100.0  | 0.0    |
   +--------------+----------+----------+------+--------+--------+
   | 493.61       | 0.1      | 17.4     | 82.5 | 100.0  | 0.0    |
   +--------------+----------+----------+------+--------+--------+
   | 493.61       | 0.1      | 17.4     | 82.5 | 100.0  | 0.0    |
   +--------------+----------+----------+------+--------+--------+
   | 493.61       | 0.1      | 17.4     | 82.5 | 100.0  | 0.0    |
   +--------------+----------+----------+------+--------+--------+

The visualisation of the internal and external contributions to the phonon modes is shown in the :numref:`fig-batio3-analysis`

.. _fig-batio3-analysis:

.. figure::  ./_static/Figures/batio3_internal_external.png
   :scale: 90%

   Internal and external contributions to the phonon modes of BaTiO\ :subscript:`3`
