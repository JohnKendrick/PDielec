# PDGui documentation refresh plan

Status: implemented following the September 2026 GUI redesign. See the
implementation record below for delivered files and verification.

The deliverable is an updated `Sphinx/pdgui.rst` with current screenshots and a
description of every visible control in each PDGui tab and application-specific
window. Include the four infrared/Raman scenario types and Powder ATR and the additional states that expose
conditional controls. Retain the established writing and figure style.

## Selected screenshot example

Use the wurtzite ZnO CASTEP calculation in
`Examples/Crystal_Raman/Castep/raman.castep` as the main dataset for both Raman
scenarios and, where useful, the common tabs and infrared scenarios. Use the
existing input and companion files in place; prepare separate documentation
scripts instead of changing the regression example scripts or spreadsheets.

The current reader was exercised during planning. It returned four atoms,
twelve phonon modes, twelve 3 by 3 Raman tensors, four Born-charge tensors,
a 3 by 3 optical permittivity tensor and a non-zero 3 by 3 by 3 nonlinear
susceptibility tensor. The CASTEP output contains the block headed
`Nonlinear Optical Susceptibility (pm/V)`. This establishes that the example
has the required EO input data; an EO-on/off spectrum comparison remains part
of preparation before capture.

The ZnO CRYSTAL example in `Examples/Powder_Raman/Crystal23/`, with `CHI2.DAT`,
was also read successfully and is a fallback. Prefer the CASTEP example for
the guide because its nonlinear susceptibility is included in the main output.
Do not reuse the existing multi-scenario scripts without checking their settings;
several deliberately select TO modes or older polarisation choices for tests.

### Raman capture settings

Use a consistent laser wavelength of 532 nm and temperature of 298 K. Start
with a 0–800 cm⁻¹ spectrum and 5 cm⁻¹ Lorentzian half-widths, then check that
all displayed peaks and axis labels are legible. These are example settings,
which must be distinguished from program defaults in the text.

* **Crystal Raman:** air / ZnO / air, (001) surface, 1 µm ZnO thickness,
  scattering matrix, normal incidence, p-in/p-out, superstrate collection and
  automatic collection angle. Select Snell's law (`geometry`) and enable the
  electro-optic term. Use incoherent depth integration and incoherent layer
  combination for the introductory example. Verify numerical convergence of
  the integration density before fixing it in the capture script.
* **Powder Raman:** the same ZnO input, a spherical particle in air and VV
  polarisation, with the electro-optic term explicitly enabled. Record the
  matrix and concentration settings. Add a non-spherical view to expose the
  direction, aspect-ratio and orientation-sampling controls that apply.
* **EO comparison:** duplicate each Raman scenario, changing only the EO
  setting. Check finite results and a non-zero change in a polar-mode
  intensity. Keep identical scales and avoid independent normalisation that
  would hide the effect. If the first geometry suppresses the relevant mode,
  choose and record an oblique geometry before capture.
* **Visible confirmation:** the main crystal screenshot must show the checked
  EO control and `Snell's law (EO)`. Show `All modes (EO)` in an additional
  view to expose its conditional controls. The powder guide should show its
  checked EO control and the corresponding calculated spectrum.

For the azimuthal sweep, use a ZnO (100) orientation if the (001) response is
uninformative under rotation. State this change explicitly in its caption.
The Raman polar window must be described according to its own calculation:
do not imply that it inherits scenario NAC/EO corrections without checking.

## Window and screenshot inventory

The first implementation step is a live GUI walkthrough to complete a control
inventory. The following is the starting capture list, based on the GUI source.

| Tab or window | Required views and description |
| --- | --- |
| Main | Loaded ZnO calculation; program and input selection, script and spreadsheet paths, actions, cell and frequency outputs. |
| Settings | Infrared and Raman views; masses, sum rules, permittivity, linewidths, mode selection, table columns and units, Raman-specific controls and actions. |
| Powder Infrared scenario | Material database, matrix, concentration, method, particle geometry and ATR controls; additional method-dependent states. |
| Crystal Infrared scenario | Optical method, incidence and azimuth, layer table, layer actions and applicable options. |
| Powder Raman scenario | EO-enabled example, VV polarisation, wavelength, temperature, particle and matrix settings; non-spherical options. |
| Crystal Raman scenario | EO-enabled main view, conditional All modes view, polarisation/Porto controls, collection geometry, integration/coherence settings, layer table and sweep controls. |
| Layer editor | Crystalline layer, orientation and laboratory-frame information; isotropic-layer view where controls differ. |
| Plotting | Infrared and Raman output views, frequency range and units, scenario selection, normalisation, plot types and export controls. |
| Analysis | Mode decomposition, frequency selection, bonding/radii settings, outputs and plot controls. |
| 3D Viewer | Phonon displacement view plus each settings subtab, display toggles, supercell/orientation, animation and file output. |
| Transformation window | Matrix entries and all actions available from the viewer. |
| Fitter | Main comparison and each settings subtab, input format, fitting choices, adjustable variables, diagnostics and actions. |
| Raman polar window | Mode selection, axes/polarisations, numerical readouts and plot/export controls. |
| Azimuthal sweep window | Polar, colour-map, waterfall and line views as needed; mode selection, grouping tolerance and plot/export controls. |

Include any further application-specific dialogs discovered in the walkthrough.
Describe standard file/colour selectors through their invoking controls; do not
create separate reference sections for platform-native dialogs.

Use ZnO for the common tabs where it provides an informative illustration.
For the Fitter, retain an established example with actual experimental input
after verifying that input, rather than presenting invented experimental ZnO
data. A molecular example may supplement Analysis if necessary to explain
molecular decomposition clearly. Identify every change of example in captions.

## Description style

Follow the existing Main, Settings and Layer Editor sections in
`Sphinx/pdgui.rst`: introduce the purpose of the tab, give a labelled figure
with a short caption, then explain its controls in connected prose in their
screen order. Use British spelling, italicised GUI labels, existing unit
macros, mathematical notation and Sphinx cross-references. Retain existing
figure anchors where their meaning is unchanged.

For each variable explain its meaning, units, editable or read-only status,
available choices, verified default and effect on the calculation or display.
Explain when a control is enabled or hidden and which other setting governs
it. Describe table columns, buttons and output readouts as well as input fields.
Use compact option tables where there are several parallel choices, as in the
existing Fitter section; avoid turning every section into a widget inventory.

Keep the scientific derivations in `Sphinx/PowderRaman.rst` and
`Sphinx/CrystalRaman.rst`, with links from the practical GUI explanations.
Explain EO prerequisites, polarisation conventions, NAC choices and the
distinct coherence settings sufficiently for users to make a choice.

Reflect the completed redesign: wavelength wording, powder VV default,
crystal p/p and Snell's law defaults, Raman controls before layer information,
and azimuthal sweep controls below it. Explain `Coalesce equivalent layers`
and `Approximate ES` in the scripting reference, including their preserved
defaults and effects, rather than listing them as GUI controls.

## Work sequence

### Input-format reference updates

Extend the scope to `Sphinx/software.rst`, with links from the Main Tab
description and examples of selecting and opening each format.

* **Experimental (`.exp`) format:** refresh the existing
  `experimental-file-format` section against `ExperimentOutputReader`.
  Cover structure, species/masses, permittivity tensors, all supported
  dielectric models and interpolation variants, plus the newer frequencies,
  normal modes, Born charges and Raman tensor blocks. Specify syntax, counts,
  tensor ordering, units, normalisation, defaults and required combinations
  of blocks. Replace the blanket statement that normal modes cannot be known
  with an explanation of which GUI features are available for each supplied
  dataset. Include a minimal dielectric example and a runnable vibrational/
  Raman example. Distinguish this format from the experimental spectrum
  spreadsheet imported by the Fitter.
* **Finite-difference / finite-field Raman JSON format:** add a dedicated
  reference for the consolidated JSON consumed by `FiniteFieldOutputReader`.
  The current GUI calls this input program `Finite_field` and its file
  selector says `Finite field (*.json)`; explain the relationship to
  finite-difference workflows without inventing a `finite_difference`
  program selector. Document accepted schema identifiers, required and
  optional keys, structure and coordinate alternatives, masses, reference
  permittivity and Born charges, mode/force-constant representations,
  dielectric derivatives and Raman tensors, and optional `chi2` for EO.
  Verify array shapes, ordering, units and normalisation directly from the
  reader, especially the internal convention expected for `chi2`.
  Explain precedence where alternative representations coexist, mass changes,
  missing-data behaviour and limitations. Include valid minimal and complete
  examples and instructions for opening them in PDGui and from the command
  line.

Validate the documented examples with the current readers and check the
resulting GUI capabilities. Include an EO-capable JSON example if a suitable
existing dataset is available; otherwise prepare a clearly attributed example
from verified data. Keep the CASTEP ZnO dataset as the principal screenshot
example. Both format references and their runnable examples are required
deliverables, not optional follow-up work.

### Implementation order

1. **Audit the final interface.** Build a checklist mapping every visible
   control and conditional state to a paragraph and screenshot. Verify labels,
   defaults and behaviour against current source and the running GUI. Update
   the stale introductory scenario list and related launch instructions.
2. **Prepare reproducible sessions.** Create small documentation-only PDGui
   scripts and a capture manifest recording dataset, settings, program version,
   window size, Qt style, font and display scaling. Validate the EO comparisons
   and useful mode/geometry choices before capturing images.
3. **Capture actual windows.** Save PNGs under `Sphinx/_static/Figures/` using
   consistent dimensions and styling. Wait for calculations and rendering to
   finish. Use additional views for hidden controls or long panels rather
   than reducing text to an unreadable size. Check the 3D view with a working
   OpenGL display; an empty headless render is not an acceptable screenshot.
4. **Rewrite the GUI chapter.** Refresh existing sections and add the two
   Raman scenario descriptions and missing auxiliary windows. Add scripting
   descriptions for the removed options. Make only the related overview and
   theory cross-reference updates needed for a coherent guide. Update both
   input-format references and validate their examples as specified above.
5. **Review and build.** Check every screenshot against the text and control
   checklist. Build Sphinx HTML in a preview directory, review the rendered
   pages, and fix new warnings, missing images and broken references. Report
   pre-existing build issues separately. Generate `docs/` only through the
   established Sphinx publication process; never edit generated HTML directly.

The implementation is complete when every application window type has a current,
legible snapshot, every visible variable has an accurate explanation, conditional
states are covered, both Raman scenarios have working EO examples, and the
rendered guide matches the established style. Preserve the existing example
datasets, reference spreadsheets and unrelated working-tree changes. The nested
`Examples/Crystal_Raman/Castep_project/` repository is outside this work.

## Implementation record

Delivered:

* `Sphinx/pdgui.rst`: refreshed GUI guide, Raman scenarios, auxiliary windows,
  scripting settings, current launch instructions and preserved main figure anchors.
* `Sphinx/gui_options.rst`: individual viewer and fitter settings views.
* `Sphinx/software.rst` and `Sphinx/input_formats.inc`: experimental format
  extensions and the new JSON reference, including units and representation precedence.
* `Sphinx/_static/Figures/gui/`: 51 actual Qt/OpenGL snapshots and a capture
  manifest recording visible controls, choices, versions and calculation checks.
* `Sphinx/_static/inputs/`: minimal `.exp` and JSON files plus complete ZnO
  vibrational `.exp` and EO-capable JSON examples converted from the existing
  CASTEP input. These are documentation fixtures, not new DFT calculations.
* `Sphinx/tools/`: reproducible capture, example generation and GUI validation
  scripts, with instructions in `README.md`.

The walkthrough also covered Powder ATR and the scenario-selection and Raman
polar-description dialogs. Viewer/Fitter subtab snapshots show the settings
strip at readable size rather than repeating the entire window. The Fitter
uses the existing sodium-sulphate ATR experiment.

Checks completed:

* The current CASTEP reader supplies non-zero Raman, Born-charge, permittivity
  and nonlinear susceptibility data for the ZnO example.
* Both Raman scenarios show a non-zero EO-on/off spectrum difference with
  identical other settings. The maximum pointwise change divided by the
  EO-on peak is approximately 0.981 for the crystal and 0.824 for powder in air.
* Increasing crystal depth integration from 40 to 80 points/µm changes the
  spectrum by 8.2e-10 relative to the peak. The manifest records the values.
* The OpenGL viewer has a valid context and its actual image was inspected.
* All four downloadable inputs parse. Both full vibrational inputs run a
  PDGui Crystal Raman calculation; the JSON example displays Snell's law (EO)
  and changes its spectrum when EO is disabled.
* The finite-field reader and experimental-default tests pass (6 tests).
* Sphinx HTML builds successfully in `/tmp/pdielec-docs-preview`. The updated
  pages have no broken local links or missing images. The final incremental
  build is warning-free; the initial full build also exposed existing API
  docstring formatting warnings in unrelated package/test sources.
* Screenshots and a rendered HTML page were visually reviewed. `git diff
  --check` passes for the documentation changes.

The codebase graph inventory used Tier 2 evidence and direct reads of changed
source. Its recorded generation was 2026-08-25T18:20:12Z, so runtime capture
and reader checks supplied current verification. Generated published `docs/`
was not hand-edited or republished; the HTML preview is available for review.
Existing regression datasets/spreadsheets and unrelated working-tree changes
were not edited by these documentation scripts.
