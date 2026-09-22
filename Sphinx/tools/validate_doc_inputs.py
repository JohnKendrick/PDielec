"""Exercise the downloadable vibrational examples through PDGui calculations."""
import os
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/pdielec-docs-matplotlib')
import numpy as np
from qtpy.QtWidgets import QApplication
from PDielec.GUI.App import App

app = QApplication([])
for filename, program in [('zno-raman.exp', 'experiment'), ('zno-raman.json', 'finite_field')]:
    window = App(['pdgui', program, str(ROOT / 'Sphinx/_static/inputs' / filename), '-spectroscopy', 'crystal_raman', '-cpus', '1', '-threading'], progressbar=None)
    nb = window.notebook
    tab = nb.scenarios[0]
    tab.refresh(force=True)
    tab.calculate(np.arange(0., 801., 2.))
    result = np.asarray(tab.raman_spectrum)
    assert np.isfinite(result).all() and np.max(result) > 0
    if filename.endswith('.json'):
        assert tab.phonon_bc_cb.currentText() == "Snell's law (EO)"
        tab.settings['Raman electro-optic term'] = False
        tab.refresh(force=True)
        tab.calculate(np.arange(0., 801., 2.))
        assert np.max(np.abs(result - np.asarray(tab.raman_spectrum))) > 0
    print(filename, 'GUI Raman calculation passed')
    nb.pool.close()
    nb.pool.join()
    window.deleteLater()
