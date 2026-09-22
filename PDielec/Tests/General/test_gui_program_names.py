"""Reader identifiers remain independent of their GUI labels."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from qtpy.QtWidgets import QApplication

from PDielec.GUI.MainTab import MainTab


@pytest.mark.parametrize(
    "program,label,stored",
    [
        ("qe", "Quantum Espresso", "Quantum espresso"),
        ("Quantum espresso", "Quantum Espresso", "Quantum espresso"),
        ("finite_field", "Finite field", "Finite_field"),
        ("Finite_field", "Finite field", "Finite_field"),
    ],
)
def test_program_labels_preserve_reader_identifiers(program, label, stored):
    """CLI and saved-session names select readable labels and survive refresh."""
    app = QApplication.instance() or QApplication([])
    tab = MainTab(None, program, "", "")
    try:
        assert tab.program_cb.currentText() == label
        assert tab.settings["Program"] == stored
        index = tab.program_cb.findText(label)
        tab.on_program_cb_activated(index)
        assert tab.settings["Program"] == stored
        tab.calculation_required = False
        tab.refresh(force=True)
        assert tab.program_cb.currentText() == label
    finally:
        tab.deleteLater()
        app.processEvents()
