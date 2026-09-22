# GaN Raman application note

Examples for Crystal Raman, Powder Raman and incoherent final-state models.
The full walkthrough is in
[the application note](../../../Sphinx/application_notes_gan.rst).

- `GaN_Schulz_Thiemann_1977.cif`: experimental starting structure.
- `Optimisation/`: CRYSTAL23 geometry optimisation inputs and results.
- `Raman/`: CRYSTAL23 Raman calculation and companion data required by PDGui.
- `Powder/`: parallel/crossed powder spectra with and without particle corrections.
- `ModalPairs/`: Mina-inspired near-forward crystal comparison and checks.
- `ModalPairs/Irmer2013/`: aperture scans and complementary finite-q response.
- `ModalPairs/IncoherentModels/`: bulk and local incoherent final-state comparisons.
- `regenerate.py`: generates the workbooks, summaries and figures.

From the repository root, in a PDielec Python environment:

```bash
python Examples/ApplicationNotes/GaN/regenerate.py
```

Regenerate before building the documentation or running the workbook checks.
Generated outputs are ignored by Git. The driver does not rerun CRYSTAL23.
