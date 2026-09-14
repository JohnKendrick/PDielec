"""Capture real PDGui windows for the user guide; run from the repository root.

Requires a working Qt/OpenGL display. Does not modify regression examples.
"""
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
os.environ.setdefault('MPLCONFIGDIR', '/tmp/pdielec-docs-matplotlib')
import numpy as np
from qtpy import QT_VERSION as QT_VERSION_STR
from qtpy.QtCore import QTimer
from qtpy.QtGui import QFont
from qtpy.QtWidgets import QApplication, QWidget, QLabel, QComboBox, QAbstractButton, QLineEdit, QAbstractSpinBox, QTableWidget
from PDielec import __version__
from PDielec.GUI.App import App
from PDielec.GUI.OpenGLWidget import OpenGLWidget
from PDielec.GUI.SingleCrystalLayer import ShowLayerWindow
from PDielec.GUI.ViewerTab import TransformWindow

OUT = ROOT / 'Sphinx/_static/Figures/gui'
OUT.mkdir(parents=True, exist_ok=True)
app = QApplication([])
app.setStyle('Fusion')
app.setFont(QFont('DejaVu Sans', 10))
manifest = {'version': __version__, 'qt': QT_VERSION_STR, 'style': 'Fusion', 'font': 'DejaVu Sans 10', 'images': {}, 'checks': {}}


def settle():
    for _ in range(5):
        app.processEvents()
        time.sleep(0.05)


def capture(widget, name, size=(1100, 900)):
    widget.resize(*size)
    widget.show()
    settle()
    # QWidget.grab() re-renders OpenGL off screen and can corrupt colours,
    # depth and QPainter text overlays. Read the displayed native window
    # instead when it contains a visible OpenGL view.
    has_opengl = any(w.isVisibleTo(widget) for w in widget.findChildren(OpenGLWidget))
    if has_opengl:
        if not widget.isWindow():
            raise ValueError('Capture the top-level window containing the OpenGL view')
        widget.raise_()
        settle()
        pixmap = widget.screen().grabWindow(int(widget.winId()))
    else:
        pixmap = widget.grab()
    if pixmap.isNull() or not pixmap.save(str(OUT / (name + '.png'))):
        raise RuntimeError('Unable to capture ' + name)
    controls = []
    for w in widget.findChildren(QWidget):
        if not w.isVisibleTo(widget):
            continue
        d = {'type': type(w).__name__, 'tooltip': w.toolTip()}
        if isinstance(w, (QLabel, QAbstractButton, QLineEdit, QAbstractSpinBox)):
            d['text'] = w.text()
        if isinstance(w, QComboBox):
            d['choices'] = [w.itemText(i) for i in range(w.count())]
            d['value'] = w.currentText()
        if isinstance(w, QTableWidget):
            d['columns'] = [w.horizontalHeaderItem(i).text() if w.horizontalHeaderItem(i) else '' for i in range(w.columnCount())]
        if any(k in d for k in ('text', 'choices', 'columns')):
            d['enabled'] = w.isEnabled()
            controls.append(d)
    manifest['images'][name] = {'size': [widget.size().width(), widget.size().height()], 'controls': controls, 'capture_method': 'displayed native window' if has_opengl else 'Qt widget'}
    print('Captured', name, flush=True)


def tabshot(nb, tab, name, size=(1100, 900)):
    nb.tabs.setCurrentWidget(tab)
    tab.refresh(force=True)
    settle()
    capture(nb.app, name, size)


def spectrum(nb, tab):
    tab.refresh(force=True)
    tab.calculate(np.arange(0., 801., 1.))
    return np.asarray(tab.raman_spectrum, dtype=float).copy()


def capture_modal(name):
    dialog = app.activeModalWidget()
    capture(dialog, name, (800, 550))
    dialog.reject()


def switch(nb, kind):
    nb.settingsTab.settings['Spectroscopy type'] = kind
    nb.set_spectroscopy_type(kind)
    nb.settingsTab.refresh(force=True)


def main():
    window = App(['pdgui', 'castep', 'Examples/Crystal_Raman/Castep/raman.castep', '-spectroscopy', 'crystal_raman', '-cpus', '1', '-threading'], progressbar=None)
    nb = window.notebook
    nb.plottingTab.settings.update({'Minimum frequency': 0., 'Maximum frequency': 800., 'Frequency increment': 1., 'Spectrum renormalisation': 'none'})
    crystal = nb.scenarios[0]
    crystal.settings.update({'Legend': 'ZnO crystal — EO on', 'Laser wavelength nm': 532., 'Temperature K': 298., 'Depth coherence': 'Incoherent intensity', 'Layer combination': 'Incoherent intensities', 'GL point density': 40., 'Raman electro-optic term': True})
    crystal.refresh(force=True)
    on = spectrum(nb, crystal)
    crystal.settings['Raman electro-optic term'] = False
    off = spectrum(nb, crystal)
    crystal.settings['Raman electro-optic term'] = True
    assert np.isfinite(on).all() and np.max(np.abs(on-off)) > 0
    manifest['checks']['crystal_eo_relative_change'] = float(np.max(np.abs(on-off))/np.max(np.abs(on)))
    crystal.settings['GL point density'] = 80.
    fine = spectrum(nb, crystal)
    manifest['checks']['crystal_density_40_to_80_relative_change'] = float(np.max(np.abs(on-fine))/np.max(np.abs(fine)))
    crystal.settings['GL point density'] = 40.
    tabshot(nb, nb.mainTab, 'main')
    tabshot(nb, nb.settingsTab, 'settings-raman')
    nb.settingsTab._nac_apply_cb.setChecked(True)
    nb.settingsTab._lo_eo_cb.setChecked(True)
    capture(window, 'settings-raman-eo')
    nb.settingsTab._nac_apply_cb.setChecked(False)
    tabshot(nb, crystal, 'crystal-raman', (1100, 1100))
    assert crystal.phonon_bc_cb.currentText() == "Snell's law (EO)"
    tabshot(nb, nb.plottingTab, 'plotting-crystal-raman')
    QTimer.singleShot(200, lambda: capture_modal('select-scenarios'))
    nb.plottingTab.open_scenario_selection_dialog()
    crystal.settings['Layer NAC mode'] = 'modal_pairs'
    tabshot(nb, crystal, 'crystal-raman-all-modes', (1100, 1200))
    crystal.settings['Layer NAC mode'] = 'geometry'
    crystal.refresh(force=True)
    for i, name in [(1, 'layer-crystal'), (0, 'layer-isotropic')]:
        dialog = ShowLayerWindow(crystal.layers[i], message='ZnO layer' if i else 'Superstrate layer')
        capture(dialog, name, (800, 550))
        dialog.hide()
    nb.settingsTab.on_raman_polar_button_clicked()
    polar = nb.settingsTab._raman_polar_windows[-1]
    capture(polar, 'raman-polar', (1250, 800))
    QTimer.singleShot(200, lambda: capture_modal('raman-polar-description'))
    polar._show_description()
    polar.hide()
    crystal.settings['Layer hkls'][1] = [1, 0, 0]
    crystal.settings['Azimuthal sweep points'] = 36
    spectrum(nb, crystal)
    crystal.on_azimuthal_sweep_btn_clicked()
    sweep = crystal._sweep_windows[-1]
    for i, name in enumerate(['polar', 'colour-map', 'waterfall', 'line']):
        sweep._plot_type_cb.setCurrentIndex(i)
        capture(sweep, 'sweep-' + name)
    sweep._plot_type_cb.setCurrentIndex(0)
    sweep._degen_cb.setCurrentIndex(1)
    capture(sweep, 'sweep-grouped')
    sweep.hide()
    crystal.settings['Layer hkls'][1] = [0, 0, 1]
    switch(nb, 'Powder Raman')
    powder = nb.scenarios[0]
    powder.settings.update({'Legend': 'ZnO powder — EO on', 'Matrix': 'air', 'Mass or volume fraction': 'volume', 'Volume fraction': 0.1, 'Raman laser wavelength': 532., 'Raman laser polarisation': 'VV', 'Raman electro-optic term': True})
    pon = spectrum(nb, powder)
    powder.settings['Raman electro-optic term'] = False
    poff = spectrum(nb, powder)
    powder.settings['Raman electro-optic term'] = True
    assert np.isfinite(pon).all() and np.max(np.abs(pon-poff)) > 0
    manifest['checks']['powder_eo_relative_change'] = float(np.max(np.abs(pon-poff))/np.max(np.abs(pon)))
    tabshot(nb, powder, 'powder-raman')
    nb.plottingTab.settings['Plot type'] = 'Powder Raman'
    tabshot(nb, nb.plottingTab, 'plotting-powder-raman')
    powder.settings['Particle shape'] = 'Ellipsoid'
    tabshot(nb, powder, 'powder-raman-ellipsoid')
    switch(nb, 'Powder Infrared')
    tabshot(nb, nb.scenarios[0], 'powder-infrared')
    tabshot(nb, nb.settingsTab, 'settings-infrared')
    tabshot(nb, nb.plottingTab, 'plotting-powder-infrared')
    switch(nb, 'Powder ATR')
    tabshot(nb, nb.scenarios[0], 'powder-atr')
    switch(nb, 'Crystal Infrared')
    tabshot(nb, nb.scenarios[0], 'crystal-infrared')
    tabshot(nb, nb.plottingTab, 'plotting-crystal-infrared')
    tabshot(nb, nb.analysisTab, 'analysis')
    nb.viewerTab.settings['Selected mode'] = 7
    tabshot(nb, nb.viewerTab, 'viewer')
    nb.viewerTab.opengl_widget.molecule_rotate(1., 25., 35., 10.)
    nb.viewerTab.opengl_widget.zoom(1.5)
    capture(window, 'viewer')
    assert nb.viewerTab.opengl_widget.isValid(), 'OpenGL context is not valid'
    for i in range(nb.viewerTab.settingsTab.count()):
        nb.viewerTab.settingsTab.setCurrentIndex(i)
        capture(nb.viewerTab.settingsTab, 'viewer-option-' + str(i), (920, 100))
    dialog = TransformWindow(nb.viewerTab.settings['Transform'])
    capture(dialog, 'transform', (500, 280))
    dialog.hide()
    manifest['settings'] = {'crystal': crystal.settings, 'powder': powder.settings}
    nb.pool.close()
    nb.pool.join()
    window.hide()
    os.chdir(ROOT / 'Examples/ATR/Na2SO42_fit')
    fitwindow = App(['pdgui', '-script', 'script.py', '-cpus', '1', '-threading'], progressbar=None)
    fitnb = fitwindow.notebook
    fitnb.fitterTab.refresh(force=True)
    tabshot(fitnb, fitnb.fitterTab, 'fitter')
    for i in range(fitnb.fitterTab.settingsTab.count()):
        fitnb.fitterTab.settingsTab.setCurrentIndex(i)
        capture(fitnb.fitterTab.settingsTab, 'fitter-option-' + str(i), (900, 100))
    fitnb.pool.close()
    fitnb.pool.join()
    fitwindow.hide()
    os.chdir(ROOT)
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2, default=str) + '\n')

def capture_viewer_only():
    """Refresh the OpenGL overview without recalculating other screenshots."""
    saved = json.loads((OUT / 'manifest.json').read_text())
    manifest.update(saved)
    window = App(['pdgui', 'castep', 'Examples/Crystal_Raman/Castep/raman.castep', '-cpus', '1', '-threading'], progressbar=None)
    nb = window.notebook
    nb.viewerTab.settings['Selected mode'] = 7
    nb.tabs.setCurrentWidget(nb.viewerTab)
    window.resize(1100, 950)
    window.show()
    nb.viewerTab.refresh(force=True)
    settle()
    gl = nb.viewerTab.opengl_widget
    assert gl.isValid(), 'OpenGL context is not valid'
    gl.molecule_rotate(1., 25., 35., 10.)
    gl.zoom(1.5)
    capture(window, 'viewer', (1100, 950))
    nb.pool.close()
    nb.pool.join()
    window.hide()
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2, default=str) + '\n')


if __name__ == '__main__':
    if '--viewer-only' in sys.argv:
        capture_viewer_only()
    else:
        main()
