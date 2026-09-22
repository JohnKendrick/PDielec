# Reproducing the PDGui documentation

Run from the repository root with PDielec and the Sphinx requirements installed:

```sh
python Sphinx/tools/make_format_examples.py
python Sphinx/tools/validate_doc_inputs.py
python Sphinx/tools/capture_pdgui.py
sphinx-build -b html Sphinx /tmp/pdielec-docs-preview
```

The capture script needs a working Qt/OpenGL display. It opens application
windows, reads the existing CASTEP ZnO and sodium-sulphate ATR examples, and
saves PNGs under `Sphinx/_static/Figures/gui/`. It uses Fusion, DejaVu Sans 10,
one worker thread and explicit example settings. Set Qt display scaling to 1
when reproducing these screenshots. No regression input or spreadsheet is
rewritten. The Fitter figures use the existing experimental data, without
running a new fit or changing its saved settings.

`manifest.json` records the PDielec/Qt versions, dimensions, visible controls,
selected values and EO/convergence checks. The script calculates EO-on/off
spectra with identical geometry and axes. Its relative-change metric is
`max(abs(on - off)) / max(abs(on))`; the integration-density check compares
40 and 80 points/µm. These checks verify that the illustrated settings affect
the calculation, not absolute agreement with an experimental cross section.

The format generator exports the existing CASTEP response into small `.exp`
and JSON files under `_static/inputs/`. The JSON `deps_dr` is reconstructed
from the complete reference mode/tensor basis. It is clearly labelled as a
format conversion, not a new finite-displacement electronic calculation. The
JSON includes the reader's internal `chi2` without an extra unit conversion.
The experimental format has no EO input block. `validate_doc_inputs.py`
exercises both vibrational examples through a PDGui Raman calculation and
checks that disabling EO changes the JSON example's spectrum.

The default Sphinx build includes generated API documentation. Existing
malformed source docstrings can produce warnings independently of these user
guide pages. Inspect the build log and require the changed guide pages to be
free of new warnings. `make html` inside `Sphinx/` also publishes generated
HTML into `docs/`; use the preview command above for review before publishing.

The OpenGL overview is captured from the displayed native window, not through
`QWidget.grab()`: off-screen widget rendering can produce incorrect colours,
depth and text overlays. Keep the viewer window visible and unobscured during
capture. To refresh only this image and its manifest entry, run:

```sh
python Sphinx/tools/capture_pdgui.py --viewer-only
```
