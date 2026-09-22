============
Installation
============

.. meta::
   :description: Installing PDielec for infrared, terahertz and Raman calculations
   :keywords: PDielec, installation, Python, PyQt6, conda, pip, Raman

Requirements and installation choices
=====================================

PDielec 10 requires **Python 3.10 or later**. Install it in a dedicated conda
or Python virtual environment so that its dependencies are kept together.
The GUI uses **PyQt6** through QtPy; the 3D viewer also requires working
OpenGL support from the operating system and graphics driver.

Choose a published package for normal use, or a Git checkout to run the
examples, develop the code, or use changes not yet published. The version
available on `PyPI <https://pypi.org/project/pdielec/>`_ or
`conda-forge <https://anaconda.org/conda-forge/pdielec>`_ may differ from the
version described by the development documentation. Check the installed
version with::

    python -c "from PDielec import __version__; print(__version__)"

Package installations provide the Python modules, GUI and command-line
entry points. Obtain the example datasets and regression tests from the Git
repository; they are not included in the installed wheel. Use a checkout
matching the installed version when comparing with reference results.

Conda
=====

With conda installed, create and activate an environment. On Windows, use a
terminal configured for conda, such as the Miniconda prompt::

    conda create --name pdielec --channel conda-forge pdielec
    conda activate pdielec
    pdgui

Update to the version available on that channel with::

    conda activate pdielec
    conda update --channel conda-forge pdielec

For development from source, create an environment containing Python and pip,
then follow the Git checkout instructions below. For example::

    conda create --name pdielec-dev --channel conda-forge python=3.12 pip
    conda activate pdielec-dev

Python 3.12 is an example environment choice, not the minimum requirement.
The package installation step below installs its declared dependencies.

PyPI
====

Use Python 3.10 or later to create a virtual environment. On Linux or macOS::

    python3 -m venv pdielec-env
    source pdielec-env/bin/activate

On Windows, the equivalent commands in Command Prompt are::

    python -m venv pdielec-env
    pdielec-env\Scripts\activate.bat

If you already activated a conda environment, use that environment instead.
In the activated environment, install and start PDielec::

    python -m pip install --upgrade pip
    python -m pip install pdielec
    pdgui

To update a PyPI installation, use ``python -m pip install --upgrade pdielec``.
Always run pip through the Python interpreter of the intended environment.
See the `Python Packaging installation guide
<https://packaging.python.org/en/latest/tutorials/installing-packages/>`_
for environment setup on other shells.

If a dependency cannot be installed from a wheel, check that it provides a
wheel for your operating system, processor architecture and Python version.
Use a supported Python version with matching wheels, or install the build
tools required by that dependency. PDielec itself is a Python package, but
several scientific and GUI dependencies contain compiled components.

Installing from GitHub
======================

Install Git and activate a Python environment as above. From a directory
where you want to keep the source, run::

    git clone https://github.com/JohnKendrick/PDielec.git
    cd PDielec
    python -m pip install -e .
    pdgui

The editable installation reads dependencies from ``pyproject.toml`` and
creates commands such as ``pdgui``, ``preader`` and ``pdmake`` in the active
environment. It works from a checkout on Linux, macOS or Windows without
manually copying scripts or relying on Windows symbolic links. Activate the
same environment whenever using those commands.

To install a fixed copy instead of an editable checkout, use
``python -m pip install .``. With an editable installation, Python source
changes take effect directly. After updating a checkout with ``git pull``,
repeat ``python -m pip install -e .`` to pick up dependency or entry-point
changes as well.

Checking an installation
------------------------

Check dependency consistency and launch the GUI::

    python -m pip check
    pdgui

From the repository root, run the pytest suites with::

    pdmake test-pytests

Run pytest followed by the example regression suites with::

    pdmake tests

On Linux or macOS, ``./pdmake`` explicitly selects the wrapper in the current
checkout. These tests need the checkout's ``Examples`` and ``PDielec/Tests``
directories. Normal regression runs compare generated results with references;
``--regenerate`` replaces references and is not an installation check.

Package dependencies
====================

The authoritative list is ``[project].dependencies`` in ``pyproject.toml``.
Pip installs these automatically when installing PDielec. The accompanying
``requirements.txt`` lists the same dependencies for environments that need
an explicit requirements file::

    dill
    imageio
    imageio-ffmpeg
    matplotlib
    numpy
    openpyxl
    psutil
    PyOpenGL
    PyQt6
    pytest
    PyYAML
    QtPy
    scikit-learn
    scipy
    setuptools
    spglib
    termcolor
    XlsxWriter

QtPy is an abstraction layer; it does not replace the PyQt6 binding.
Pytest is currently a declared dependency because PDielec includes test-running
commands. A separate ``mkl`` installation is not required by PDielec.

Optional tools
--------------

**Phonopy:** install ``phonopy`` to use ``phonopy-pdielec-born``, which calls
the Phonopy API to expand Born charges. It is not in the core dependency list::

    python -m pip install phonopy

**Documentation:** from the repository root, install the documentation tools::

    python -m pip install -r Sphinx/requirements.txt

This includes Sphinx, its configured extensions and theme, and
``sphinx-autobuild`` for live previews. With Make available, ``make livehtml``
in ``Sphinx/`` starts a preview that updates when sources change. A direct
cross-platform alternative, run from the repository root, is::

    sphinx-autobuild -b html Sphinx Sphinx/_build/html

For a static preview, use ``sphinx-build -b html Sphinx Sphinx/_build/html``.
The custom ``make html`` target also copies built HTML into ``docs/``.
Screenshot regeneration needs a working Qt/OpenGL display; see
``Sphinx/tools/README.md`` for the capture and input-validation tools.

**PDF manual:** install the additional Python converter, together with a TeX
installation providing ``pdflatex`` and ``latexmk``, and the system
``rsvg-convert`` command (provided by librsvg). From the repository root, run::

    python -m pip install -r Sphinx/requirements-pdf.txt
    make -C Sphinx pdf

Alternatively, if you are already inside ``Sphinx/``, run::

    make pdf

This builds the manual and copies it to ``PDielec.pdf`` in the repository root
after a successful build. The intermediate manual remains at
``Sphinx/_build/latex/PDielec.pdf``; ``make latexpdf`` builds only that copy.
SVG figures are converted to vector PDF files automatically in the build
directory. HTML continues to use the original
SVG files and does not require the PDF converter or TeX. Unicode mappings for
scientific symbols are applied only to the LaTeX output.

**Release builds:** install the build frontend before creating archives::

    python -m pip install build
    python -m build

The build environment installs the backend declared in ``pyproject.toml``,
including ``setuptools>=77.0.3`` for SPDX license metadata. ``pdmake pypi``
also builds archives on supported platforms, but removes previous build
output first. Building archives does not upload them to PyPI. Maintainers
can install ``twine`` separately to check or upload a release.

GitHub - PDielec directory structure
====================================

Paths below are relative to the repository root:

* ``PDielec/`` contains the readers and numerical routines; ``PDielec/GUI/`` contains PDGui.
* ``PDielec/Tests/`` contains pytest regression tests.
* ``Examples/`` contains input/output datasets, regression recipes and helper scripts.
* ``Sphinx/`` contains the documentation sources; ``docs/`` contains published HTML.
* ``pyproject.toml`` defines dependencies, package data and installed command entry points.

.. _Examples:

Examples
========

The catalogue below describes the tracked example families. A directory with
``command.pdmake`` contains a regression recipe. PDGui recipes normally use
``script.py`` and compare a generated ``results.xlsx`` with
``results.ref.xlsx``. Reader and analysis recipes instead use CSV, CIF or
normal-mode-analysis references. Supporting files and notebooks are not all
standalone regression tests, and ``pdmake tests`` runs the configured suites,
not every file below ``Examples/``.

From an example directory containing a PDGui recipe, run::

    pdmake command.pdmake
    pdmake --view command.pdmake

The second command opens the saved scenario interactively. From the repository
root, ``pdmake tests`` runs pytest first and then the configured example suites.
``pdmake benchmarks`` runs the benchmark selection without reference comparison;
see :ref:`performance`. Use ``--regenerate`` only after reviewing differences.

The small ``Crystal_Raman/AbInit`` and ``Crystal_Raman/Finite_field`` directories
currently contain provenance notes only. Use ``Powder_Raman/AbInit`` and
``Powder_Raman/Finite_difference`` for the corresponding complete datasets.
The CASTEP 25.12 directory includes both a crystal Raman recipe and data used
by the atomic Raman reader tests.

.. list-table:: Examples available in the repository
   :header-rows: 1
   :widths: 28 16 56

   * - Directory
     - Input or workflow
     - Purpose
   * - ``ATR/AlAs``
     - AbInit
     - Maxwell Garnett calculation of the ATR spectrum of an ellipsoid along [001]. The incident angle varies from 0 to 80 degrees.
   * - ``ATR/Na2SO42``
     - Vasp
     - Maxwell-Garnett calculation of the ATR spectrum, changes the S polarisation component from 0 to 100%
   * - ``ATR/Na2SO42_fit``
     - Vasp
     - Maxwell-Garnett calculation of the ATR spectrum, an example of fitting the spectrum to experiment
   * - ``AbInit/AlAs``
     - AbInit
     - Average permittivity and Maxwell-Garnett calculation of sphere, plate and ellipsoid
   * - ``AbInit/BaTiO3``
     - AbInit
     - Average permittivity and Maxwell-Garnett calculations of sphere, plate and ellipsoid, using average isotope masses
   * - ``AbInit/Na2SO42``
     - AbInit
     - Average permittivity and Maxwell-Garnett calculations of Na2(SO4)2, sphere, plate and ellipsoid, using program-defined masses
   * - ``AbInit/preader``
     - AbInit
     - AbInit preader test
   * - ``Castep/AsparticAcid``
     - Castep
     - Average permittivity and Maxwell-Garnett calculations of sphere, plate and ellipsoid, using program-defined masses
   * - ``Castep/Bubbles``
     - Castep
     - Maxwell-Garnett calculation showing the effect of air bubbles at 24% volume fraction and 30 micron radius
   * - ``Castep/Castep17``
     - Castep
     - Castep 17, Maxwell-Garnett sphere and plates with 3 surfaces
   * - ``Castep/Isoleucine``
     - Castep
     - Maxwell-Garnett sphere
   * - ``Castep/MgO``
     - Castep
     - Comparison of MG, Bruggeman and AP methods changing shapes and volume fractions
   * - ``Castep/Na2SO42``
     - Castep
     - Comparison of MG and Bruggeman, for needle, ellipsoid and plate shapes
   * - ``Castep/preader``
     - Castep
     - Castep preader test
   * - ``Crystal/Leucine``
     - Crystal
     - Comparison of MG, plates and ellipsoids
   * - ``Crystal/Na2SO42``
     - Crystal
     - Comparison of MG for needle, ellipsoid and plate shapes
   * - ``Crystal/Na2SO42_C17``
     - Crystal
     - Comparison of MG for needle, ellipsoid and plate shapes, reading output from Crystal 17
   * - ``Crystal/Quartz``
     - Crystal
     - Comparison of MG for needle, ellipsoid and plate shapes
   * - ``Crystal/ZnO``
     - Crystal
     - PDGui Crystal of ZnO (CPHF); PDGui Crystal of ZnO (default settings); PDGui Crystal of ZnO (NOECKART settings)
   * - ``Crystal/preader``
     - Crystal
     - Crystal preader test
   * - ``Crystal_Raman/Castep``
     - Crystal Raman
     - PDGui Castep Crystal Raman of ZnO
   * - ``Crystal_Raman/Castep-25.12``
     - Crystal Raman
     - ZnO crystal Raman scenarios and CASTEP 25.12 atomic polar tensor reader checks.
   * - ``Crystal_Raman/Crystal23``
     - Crystal Raman
     - CRYSTAL23 ZnO crystal Raman scenarios, including optical geometry and coherence.
   * - ``Experiment/AlN``
     - Experiment
     - Aluminium Nitride multi-layer system including SiC and Si
   * - ``Experiment/Forsterite``
     - Experiment
     - Single crystal calculations of a thick slab, for a, b and c axis alignments with polarisation direction. Uses FPSQ model for permittivity.
   * - ``Experiment/Mayerhofer``
     - Experiment
     - Example of a Drude Lorentz model permittivity
   * - ``Experiment/Sapphire``
     - Experiment
     - Sapphire example and test of the materials database
   * - ``Experiment/constant``
     - Experiment
     - Example of a constant permittivity with loss
   * - ``Experiment/cp2k``
     - Experiment
     - PDGui experiment Single Crystal for cysteine form3
   * - ``Experiment/drude-lorentz``
     - Experiment
     - A Drude-Lorentz model for MgO, varying the angle of incidence
   * - ``Experiment/fpsq``
     - Experiment
     - An FPSQ model for Quartz, showing polarisation on along different axes and different incident angles.
   * - ``Experiment/interpolation``
     - Experiment
     - An example of an interpolation model
   * - ``FHI-Aims/Na2SO42``
     - FHI-Aims
     - Powder calculation for spheres, including a comparison with VASP results
   * - ``FHI-Aims/preader``
     - FHI-Aims
     - FHI-Aims preader test
   * - ``Finite_difference/preader``
     - Finite-field JSON
     - Finite-difference VASP ZnO preader infrared (Eckart)
   * - ``Finite_difference/preader_raman``
     - Finite-field JSON
     - Finite-difference VASP ZnO preader Raman
   * - ``Gulp/Na2SO42``
     - Gulp
     - Maxwell-Garnett and Bruggeman on needle, ellipsoid and plate
   * - ``Gulp/calcite``
     - Gulp
     - Maxwell-Garnett method on Sphere and Plate
   * - ``Gulp/preader``
     - Gulp
     - Gulp preader test
   * - ``Helper/Helper``
     - Python API
     - Python helper API demonstration; inspect test_helper.py for its workflow.
   * - ``Helper/Jupyter``
     - Python API
     - Jupyter notebooks and supporting datasets for API workflows.
   * - ``Mie/MgO``
     - Castep
     - Mie method with varying volume fractions and sphere sizes
   * - ``Mie/MgO_lognormal``
     - Castep
     - Mie method with varying volume fractions and sphere size distributions
   * - ``P2Cif``
     - CIF export
     - Testing p2cif
   * - ``Phonopy/Al2O3``
     - Phonopy
     - Primitive and standard cell example of using Phonopy
   * - ``Phonopy/Na2SO42``
     - Phonopy
     - Maxwell-Garnett and Bruggeman method for needle, ellipsoid and plate shapes, with varying volume fractions
   * - ``Phonopy/QE``
     - Phonopy/QE
     - Powder and single crystal Phonopy example using QE
   * - ``Phonopy/Vasp``
     - Phonopy/Vasp
     - Powder and single crystal Phonopy example using Vasp
   * - ``Phonopy/ZnO``
     - Phonopy
     - Maxwell-Garnett and Bruggeman method for needle, ellipsoid and plate shapes
   * - ``Phonopy/preader``
     - Phonopy
     - Phonopy preader test
   * - ``Powder_Raman/AbInit``
     - Powder Raman
     - PDGui Abinit powder Raman of ZnO
   * - ``Powder_Raman/Castep``
     - Powder Raman
     - PDGui Castep powder Raman of ZnO
   * - ``Powder_Raman/Crystal23``
     - Powder Raman
     - PDGui Crystal powder Raman of ZnO
   * - ``Powder_Raman/Finite_difference``
     - Powder Raman
     - VASP ZnO finite-difference JSON; powder Raman with EO disabled and enabled.
   * - ``Powder_Raman/QE``
     - Powder Raman
     - PDGui QE powder Raman of ZnO
   * - ``Powder_Raman/Vasp``
     - Powder Raman
     - PDGui Vasp powder Raman of ZnO
   * - ``QE/Cocaine``
     - Quantum Espresso
     - Maxwel-Garnett sphere, using QE 4.1
   * - ``QE/Na2SO42``
     - Quantum Espresso
     - Maxwell-Garnett and Bruggeman on needle, ellipsoid and plate, using QE 5.1
   * - ``QE/Na2SO42-v7``
     - Quantum Espresso
     - Maxwell-Garnett and Bruggeman on needle, ellipsoid and plate, using QE 7.3.1
   * - ``QE/Urea``
     - Quantum Espresso
     - Maxwell-Garnett and single crystal, using QE 7.3.1
   * - ``QE/ZnO``
     - Quantum Espresso
     - Maxwell-Garnett and Bruggeman on needle, ellipsoid and plate, using QE 5.4.0
   * - ``QE/preader``
     - QE
     - QE preader test
   * - ``SingleCrystal/Bi2Se3``
     - Vasp
     - Single crystal example of thick slab, angle of incidence varies from 0 to 90
   * - ``SingleCrystal/Bi2Se3_film``
     - Vasp
     - Single crystal example of thin film, angle of incidence varies from 0 to 90
   * - ``SingleCrystal/L-alanine``
     - SingleCrystal
     - PDGui Single crystal L-alanine test
   * - ``SingleCrystal/MgO``
     - SingleCrystal
     - PDGui Castep of MgO Single Crystal coherent
   * - ``SingleCrystal/ScatteringMatrix``
     - SingleCrystal
     - PDGui Single crystal scattering matrix test
   * - ``SizeEffects/BaTiO3``
     - Abinit
     - Exploration of size effects in Bruggeman effective medium theory
   * - ``SizeEffects/MgO``
     - Castep
     - Exploration of size effects in Bruggeman and Maxwell-Garnett effective medium theories
   * - ``SizeEffects/ZnO``
     - Vasp
     - Exploration of size effects in Maxwell-Garnett effective medium theory
   * - ``Vasp/F-Apatite``
     - Vasp
     - Maxwell-Garnett, sphere plates and needles, using Vasp 5.3.5
   * - ``Vasp/Na2SO42``
     - Vasp
     - Maxwell-Garnett and Bruggeman, needle, plate and needle, using Vasp 5.3.5
   * - ``Vasp/Na2SO42_v``
     - Vasp
     - PDGui Vasp for Na2(SO4)2 v dependent matrix permittivities
   * - ``Vasp/Urea``
     - Vasp
     - Powder and single crystal exampl, using Vasp 5.4.4
   * - ``Vasp/Vasp6``
     - Vasp
     - The DFT calculation were performed using VASP 6.4.2, with thanks to Dr David Santos-Carballal for performing the DFT calculation
   * - ``Vasp/ZnO``
     - Vasp
     - Maxwell-Garnett and Bruggeman, needle, plate and needle, mass fraction, using Vasp 5.3.5
   * - ``Vasp/preader``
     - Vasp
     - Vasp preader test
   * - ``VibAnalysis/AsparticAcid``
     - Mode analysis
     - Vibanalysis Castep for L-Aspartic acid
   * - ``VibAnalysis/BaTiO3``
     - Mode analysis
     - Vibanalysis Abinit for BaTiO3
   * - ``VibAnalysis/Cocaine``
     - Mode analysis
     - Vibanalysis QE for cocaine salt
   * - ``VibAnalysis/F-Apatite``
     - Mode analysis
     - Vibanalysis Vasp for F-Apatite
   * - ``VibAnalysis/Isoleucine``
     - Mode analysis
     - Vibanalysis Castep for Isoleucine
   * - ``VibAnalysis/Leucine``
     - Mode analysis
     - Vibanalysis Crystal for Leucine
