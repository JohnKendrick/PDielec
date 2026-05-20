#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""AzimuthalSweepWindow module — pop-up window showing Raman intensity vs azimuthal angle."""
import matplotlib.cm
import matplotlib.colors
import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget,
)


def group_degenerate_modes(frequencies, intensities_2d, tol_cm1=2.0):
    """Group near-degenerate modes by summing their intensities.

    Parameters
    ----------
    frequencies : ndarray, shape (M,)
        Mode frequencies in cm⁻¹.
    intensities_2d : ndarray, shape (N, M)
        Per-mode intensities at each azimuthal angle.
    tol_cm1 : float, optional
        Frequency tolerance in cm⁻¹.  Modes whose frequencies span no more
        than this range are merged into a single group.  Default is 2.0.

    Returns
    -------
    grouped_freqs : ndarray, shape (G,)
        Mean frequency for each group.
    grouped_ints : ndarray, shape (N, G)
        Summed intensities for each group.
    group_labels : list of str
        Human-readable label for each group.

    """
    M = len(frequencies)
    if M == 0:
        return frequencies.copy(), intensities_2d.copy(), []

    sort_idx = np.argsort(frequencies)
    sorted_freqs = frequencies[sort_idx]
    sorted_ints = intensities_2d[:, sort_idx]

    # Greedy consecutive grouping: extend the current group as long as the
    # span from its first element does not exceed tol_cm1.
    groups = []
    current = [0]
    for i in range(1, M):
        if sorted_freqs[i] - sorted_freqs[current[0]] <= tol_cm1:
            current.append(i)
        else:
            groups.append(current)
            current = [i]
    groups.append(current)

    grouped_freqs = np.array([np.mean(sorted_freqs[g]) for g in groups])
    grouped_ints = np.column_stack([sorted_ints[:, g].sum(axis=1) for g in groups])

    labels = []
    for g in groups:
        gf = sorted_freqs[g]
        if len(gf) == 1:
            labels.append(f"{gf[0]:.1f} cm\u207b\u00b9")
        else:
            labels.append(
                f"{gf[0]:.1f}\u2013{gf[-1]:.1f} cm\u207b\u00b9 ({len(gf)} modes)"
            )

    return grouped_freqs, grouped_ints, labels


class AzimuthalSweepWindow(QWidget):
    """Non-modal window displaying Raman intensity vs azimuthal angle in four plot styles.

    Parameters
    ----------
    sweep_data : dict
        Pre-computed sweep data with keys:

        ``psi_array`` : ndarray, shape (N,)
            Azimuthal angles in degrees.
        ``mode_frequencies`` : ndarray, shape (M,)
            Active mode frequencies in cm⁻¹.
        ``mode_intensities_2d`` : ndarray, shape (N, M)
            Per-mode intensities at each azimuthal angle.
        ``vs_cm1`` : ndarray, shape (K,)
            Wavenumber axis in cm⁻¹.
        ``spectra_2d`` : ndarray, shape (N, K)
            Broadened spectra at each azimuthal angle.

    title : str, optional
        Window title string.
    parent : QWidget or None, optional
        Parent widget (pass ``None`` for an independent top-level window).

    """

    def __init__(self, sweep_data, title="Azimuthal sweep", parent=None):
        """Initialise the window and draw the initial plot."""
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setWindowFlags(self.windowFlags() | Qt.Window)

        self._psi = np.asarray(sweep_data["psi_array"])
        self._mode_freqs = np.asarray(sweep_data["mode_frequencies"])
        self._mode_ints = np.asarray(sweep_data["mode_intensities_2d"])
        self._vs = np.asarray(sweep_data["vs_cm1"])
        self._spectra = np.asarray(sweep_data["spectra_2d"])

        # --- Controls row ---
        controls_hbox = QHBoxLayout()

        plot_type_label = QLabel("Plot type:")
        self._plot_type_cb = QComboBox(self)
        self._plot_type_cb.addItems(["Polar", "Colour map", "Waterfall", "Line"])
        controls_hbox.addWidget(plot_type_label)
        controls_hbox.addWidget(self._plot_type_cb)

        controls_hbox.addItem(QSpacerItem(20, 0, QSizePolicy.Expanding, QSizePolicy.Minimum))

        mode_label = QLabel("Mode:")
        self._mode_cb = QComboBox(self)
        self._mode_cb.addItem("All modes (integrated)")
        for idx, freq in enumerate(self._mode_freqs):
            self._mode_cb.addItem(f"Mode {idx + 1}: {freq:.1f} cm⁻¹")
        controls_hbox.addWidget(mode_label)
        controls_hbox.addWidget(self._mode_cb)

        degen_label = QLabel("Mode display:")
        self._degen_cb = QComboBox(self)
        self._degen_cb.addItems(["Individual eigenmodes", "Group near-degenerate"])
        self._degen_cb.setToolTip(
            "Individual eigenvectors inside a degenerate subspace are basis-dependent.\n"
            "Grouping near-degenerate modes shows the physically stable summed intensity."
        )
        controls_hbox.addWidget(degen_label)
        controls_hbox.addWidget(self._degen_cb)

        self._tol_label = QLabel("Tol (cm\u207b\u00b9):")
        self._tol_sb = QDoubleSpinBox(self)
        self._tol_sb.setRange(0.1, 50.0)
        self._tol_sb.setValue(2.0)
        self._tol_sb.setSingleStep(0.5)
        self._tol_sb.setDecimals(1)
        self._tol_sb.setToolTip("Frequency tolerance for grouping near-degenerate modes.")
        self._tol_label.setVisible(False)
        self._tol_sb.setVisible(False)
        controls_hbox.addWidget(self._tol_label)
        controls_hbox.addWidget(self._tol_sb)

        # --- Matplotlib canvas ---
        self.figure = matplotlib.figure.Figure(figsize=(7.2, 6.0), constrained_layout=False)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumSize(480, 360)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        # --- Main layout ---
        vbox = QVBoxLayout(self)
        vbox.addLayout(controls_hbox)
        vbox.addWidget(self.toolbar)
        vbox.addWidget(self.canvas)
        self.setLayout(vbox)

        # Connect signals after layout is built
        self._plot_type_cb.currentIndexChanged.connect(self._on_controls_changed)
        self._mode_cb.currentIndexChanged.connect(self._on_controls_changed)
        self._degen_cb.currentIndexChanged.connect(self._on_degen_changed)
        self._tol_sb.valueChanged.connect(self._on_tol_changed)

        self.resize(720, 660)
        self._update_plot()

    # ------------------------------------------------------------------
    # Signal handlers
    # ------------------------------------------------------------------

    def _on_controls_changed(self):
        """React to plot-type or mode combo-box changing."""
        plot_type = self._plot_type_cb.currentText()
        # Mode selection is irrelevant for colour map and waterfall
        disable_mode = plot_type in ("Colour map", "Waterfall")
        self._mode_cb.setEnabled(not disable_mode)
        self._update_plot()

    def _on_degen_changed(self):
        """React to degenerate-mode display combo changing."""
        grouping = self._degen_cb.currentText() == "Group near-degenerate"
        self._tol_label.setVisible(grouping)
        self._tol_sb.setVisible(grouping)
        self._rebuild_mode_combo()
        self._update_plot()

    def _on_tol_changed(self):
        """React to tolerance spinbox changing."""
        self._rebuild_mode_combo()
        self._update_plot()

    # ------------------------------------------------------------------
    # Mode grouping helpers
    # ------------------------------------------------------------------

    def _get_effective_modes(self):
        """Return (frequencies, intensities_2d, labels) with optional grouping.

        Returns
        -------
        freqs : ndarray, shape (G,)
        ints : ndarray, shape (N, G)
        labels : list of str, length G

        """
        if self._degen_cb.currentText() == "Group near-degenerate":
            tol = self._tol_sb.value()
            return group_degenerate_modes(self._mode_freqs, self._mode_ints, tol)
        labels = [f"{f:.1f} cm\u207b\u00b9" for f in self._mode_freqs]
        return self._mode_freqs, self._mode_ints, labels

    def _rebuild_mode_combo(self):
        """Repopulate the mode combo to match the current grouping."""
        freqs, _ints, labels = self._get_effective_modes()
        prev_idx = self._mode_cb.currentIndex()
        self._mode_cb.blockSignals(True)
        self._mode_cb.clear()
        self._mode_cb.addItem("All modes (integrated)")
        for idx, label in enumerate(labels):
            self._mode_cb.addItem(f"Mode {idx + 1}: {label}")
        # Restore position if still valid; else default to "All modes"
        self._mode_cb.setCurrentIndex(min(prev_idx, self._mode_cb.count() - 1))
        self._mode_cb.blockSignals(False)

    # ------------------------------------------------------------------
    # Dispatch
    # ------------------------------------------------------------------

    def _update_plot(self):
        """Clear the figure and redraw according to the current control state."""
        self.figure.clear()
        plot_type = self._plot_type_cb.currentText()
        mode_idx = self._mode_cb.currentIndex() - 1  # -1 means "all modes"
        freqs, ints, labels = self._get_effective_modes()

        if plot_type == "Polar":
            self._plot_polar(self._psi, freqs, ints, labels, mode_idx)
        elif plot_type == "Colour map":
            self._plot_colourmap(self._psi, self._vs, self._spectra)
        elif plot_type == "Waterfall":
            self._plot_waterfall(self._psi, self._vs, self._spectra)
        else:  # "Line"
            self._plot_line(self._psi, freqs, ints, labels, mode_idx)

        self.canvas.draw()

    # ------------------------------------------------------------------
    # Individual plot methods
    # ------------------------------------------------------------------

    def _plot_polar(self, psi, mode_freqs, mode_ints, labels, mode_idx):
        """Draw a polar plot of Raman intensity vs azimuthal angle.

        Parameters
        ----------
        psi : ndarray
            Azimuthal angles in degrees, shape (N,).
        mode_freqs : ndarray
            Mode frequencies in cm⁻¹, shape (M,).
        mode_ints : ndarray
            Per-mode intensities, shape (N, M).
        labels : list of str
            Display label for each mode/group, length M.
        mode_idx : int
            Index of selected mode, or -1 for the integrated sum over all modes.

        """
        ax = self.figure.add_subplot(111, projection="polar")
        theta = np.radians(psi)
        # Close the loop
        theta_full = np.append(theta, theta[0])

        if mode_idx < 0:
            r = mode_ints.sum(axis=1)
            label = "Integrated"
        else:
            r = mode_ints[:, mode_idx]
            label = labels[mode_idx]
        r_full = np.append(r, r[0])

        ax.plot(theta_full, r_full)
        ax.set_title(f"Raman intensity vs azimuthal angle\n{label}")

    def _plot_colourmap(self, psi, vs, spectra):
        """Draw a colour-map (image) of the broadened Raman spectrum vs azimuthal angle.

        Parameters
        ----------
        psi : ndarray
            Azimuthal angles in degrees, shape (N,).
        vs : ndarray
            Wavenumber axis in cm⁻¹, shape (K,).
        spectra : ndarray
            Broadened spectra, shape (N, K).

        """
        ax = self.figure.add_subplot(111)
        im = ax.imshow(
            spectra,
            aspect="auto",
            origin="upper",
            extent=[vs[0], vs[-1], psi[-1], psi[0]],
            cmap="viridis",
        )
        cbar = self.figure.colorbar(im, ax=ax)
        cbar.set_label("Raman intensity")
        ax.set_xlabel("Raman shift (cm⁻¹)")
        ax.set_ylabel("Azimuthal angle ψ (°)")
        ax.set_title("Raman spectrum vs azimuthal angle")

    def _plot_waterfall(self, psi, vs, spectra):
        """Draw a waterfall plot of the broadened Raman spectrum vs azimuthal angle.

        Parameters
        ----------
        psi : ndarray
            Azimuthal angles in degrees, shape (N,).
        vs : ndarray
            Wavenumber axis in cm⁻¹, shape (K,).
        spectra : ndarray
            Broadened spectra, shape (N, K).

        """
        ax = self.figure.add_subplot(111)
        stride = max(1, len(psi) // 36)
        shown_indices = list(range(0, len(psi), stride))
        n_shown = len(shown_indices)
        cmap = matplotlib.cm.get_cmap("viridis", n_shown)  # noqa: RUF100
        vmax = spectra.max() or 1.0
        offset_step = vmax * 0.8

        for plot_idx, i in enumerate(shown_indices):
            ax.plot(vs, spectra[i] + offset_step * plot_idx, color=cmap(plot_idx), lw=0.8)

        sm = matplotlib.cm.ScalarMappable(
            cmap=matplotlib.cm.get_cmap("viridis"),  # noqa: RUF100
            norm=matplotlib.colors.Normalize(vmin=psi[0], vmax=psi[-1]),
        )
        sm.set_array([])
        cbar = self.figure.colorbar(sm, ax=ax)
        cbar.set_label("Azimuthal angle ψ (°)")
        ax.set_xlabel("Raman shift (cm⁻¹)")
        ax.set_ylabel("Intensity + offset")
        ax.set_title("Waterfall: Raman spectrum vs azimuthal angle")

    def _plot_line(self, psi, mode_freqs, mode_ints, labels, mode_idx):
        """Draw a line plot of Raman mode intensity vs azimuthal angle.

        Parameters
        ----------
        psi : ndarray
            Azimuthal angles in degrees, shape (N,).
        mode_freqs : ndarray
            Mode frequencies in cm⁻¹, shape (M,).
        mode_ints : ndarray
            Per-mode intensities, shape (N, M).
        labels : list of str
            Display label for each mode/group, length M.
        mode_idx : int
            Index of selected mode, or -1 to plot all modes.

        """
        ax = self.figure.add_subplot(111)
        ax.set_xlabel("Azimuthal angle ψ (°)")
        ax.set_ylabel("Raman intensity")

        if mode_idx < 0:
            n_modes = len(mode_freqs)
            cmap = matplotlib.cm.get_cmap("tab10", max(n_modes, 1))  # noqa: RUF100
            for m_idx in range(n_modes):
                ax.plot(psi, mode_ints[:, m_idx], color=cmap(m_idx), label=labels[m_idx])
            if n_modes > 0:
                ax.legend(fontsize="small")
            ax.set_title("Raman mode intensities vs azimuthal angle")
        else:
            ax.plot(psi, mode_ints[:, mode_idx])
            ax.set_title(f"Mode {mode_idx + 1}: {labels[mode_idx]} vs azimuthal angle")
