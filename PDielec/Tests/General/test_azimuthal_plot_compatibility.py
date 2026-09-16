"""Render azimuthal plots without Matplotlib's removed cm.get_cmap API."""

from types import SimpleNamespace

import matplotlib.cm
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from PDielec.GUI.AzimuthalSweepWindow import AzimuthalSweepWindow


def test_waterfall_renders_without_legacy_colormap_api(monkeypatch):
    """Keep the three spectra and their offsets when using current Matplotlib."""
    monkeypatch.delattr(matplotlib.cm, "get_cmap", raising=False)
    figure = Figure()
    canvas = FigureCanvasAgg(figure)
    window = SimpleNamespace(figure=figure)
    AzimuthalSweepWindow._plot_waterfall(
        window, np.array([0, 90, 180]), np.array([100, 200]), np.array([[1, 2], [3, 4], [5, 6]]),
    )
    canvas.draw()
    lines = figure.axes[0].lines
    assert len(lines) == 3
    np.testing.assert_allclose(lines[0].get_ydata(), [1, 2])
    np.testing.assert_allclose(lines[2].get_ydata(), [14.6, 15.6])
    assert len(figure.axes) == 2  # Spectrum axes and colour bar.


def test_all_mode_lines_render_without_legacy_colormap_api(monkeypatch):
    """Render each mode with its data and legend on current Matplotlib."""
    monkeypatch.delattr(matplotlib.cm, "get_cmap", raising=False)
    figure = Figure()
    canvas = FigureCanvasAgg(figure)
    window = SimpleNamespace(figure=figure)
    intensities = np.array([[1, 2], [3, 4], [5, 6]])
    AzimuthalSweepWindow._plot_line(
        window, np.array([0, 90, 180]), np.array([100, 200]), intensities, ["Mode A", "Mode B"], -1,
    )
    canvas.draw()
    lines = figure.axes[0].lines
    assert len(lines) == 2
    for line, expected in zip(lines, intensities.T):
        np.testing.assert_allclose(line.get_ydata(), expected)
    assert [text.get_text() for text in figure.axes[0].get_legend().get_texts()] == ["Mode A", "Mode B"]
