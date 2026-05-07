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
"""Qt window for tensor-level Raman polar mode plots."""

import logging

import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QAbstractItemView,
    QDoubleSpinBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from PDielec.RamanPolarCalculator import POLAR_GEOMETRIES, calculate_polar_curves

logger = logging.getLogger(__name__)


class RamanPolarWindow(QWidget):
    """Non-modal Raman polar plot simulator window."""

    def __init__(self, polar_data, title="Raman polar mode simulator", parent=None):
        """Initialise the window."""
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setWindowFlags(self.windowFlags() | Qt.Window)

        self._raman_tensors = polar_data["raman_tensors"]
        self._modes = polar_data["modes"]
        self._angles_rad = np.linspace(0.0, 2.0 * np.pi, 361)

        main_hbox = QHBoxLayout(self)
        controls_vbox = QVBoxLayout()
        plot_vbox = QVBoxLayout()

        self._mode_table = QTableWidget(self)
        self._mode_table.setColumnCount(4)
        self._mode_table.setHorizontalHeaderLabels(["Include", "Mode", "Frequency (cm-1)", "Raman total"])
        self._mode_table.verticalHeader().setVisible(False)
        self._mode_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._mode_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._mode_table.setRowCount(len(self._modes))
        self._mode_table.itemChanged.connect(self._on_mode_item_changed)
        self._populate_mode_table()
        self._mode_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self._mode_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self._mode_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self._mode_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)
        self._mode_table.setMinimumWidth(430)

        mode_label = QLabel("Raman-active modes", self)
        controls_vbox.addWidget(mode_label)
        controls_vbox.addWidget(self._mode_table)

        buttons_hbox = QHBoxLayout()
        self._select_all_btn = QPushButton("Select all active modes", self)
        self._select_all_btn.clicked.connect(self._select_all_modes)
        self._clear_btn = QPushButton("Clear selection", self)
        self._clear_btn.clicked.connect(self._clear_modes)
        buttons_hbox.addWidget(self._select_all_btn)
        buttons_hbox.addWidget(self._clear_btn)
        controls_vbox.addLayout(buttons_hbox)

        controls_vbox.addWidget(QLabel("Light direction Rz", self))
        self._rz_spins = self._add_vector_controls(controls_vbox, (0.0, 0.0, 1.0))
        controls_vbox.addWidget(QLabel("Zero-angle polarisation Pza", self))
        self._pza_spins = self._add_vector_controls(controls_vbox, (1.0, 0.0, 0.0))
        self._orthogonality_label = QLabel("Orthogonality: 0.000000", self)
        controls_vbox.addWidget(self._orthogonality_label)

        self._raw_max_label = QLabel("Raw max intensity: 0.0", self)
        controls_vbox.addWidget(self._raw_max_label)

        self._ack_btn = QPushButton("Acknowledgement", self)
        self._ack_btn.clicked.connect(self._show_acknowledgement)
        controls_vbox.addWidget(self._ack_btn)
        controls_vbox.addStretch(1)

        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        plot_vbox.addWidget(self.toolbar)
        plot_vbox.addWidget(self.canvas)

        main_hbox.addLayout(controls_vbox)
        main_hbox.addLayout(plot_vbox, stretch=1)
        self.setLayout(main_hbox)
        self.resize(1050, 700)
        self._update_plot()

    def _populate_mode_table(self):
        """Populate the mode table."""
        self._mode_table.blockSignals(True)
        for row, mode in enumerate(self._modes):
            include_item = QTableWidgetItem("")
            include_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            include_item.setCheckState(Qt.Checked if mode["selected"] else Qt.Unchecked)
            self._mode_table.setItem(row, 0, include_item)

            mode_item = QTableWidgetItem(str(mode["display_mode"]))
            freq_item = QTableWidgetItem(f"{mode['frequency_cm1']:.4f}")
            activity_item = QTableWidgetItem(f"{mode['raman_total']:.6g}")
            for item in (mode_item, freq_item, activity_item):
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self._mode_table.setItem(row, 1, mode_item)
            self._mode_table.setItem(row, 2, freq_item)
            self._mode_table.setItem(row, 3, activity_item)
        self._mode_table.blockSignals(False)

    def _add_vector_controls(self, layout, values):
        """Add x/y/z spin boxes to a layout."""
        hbox = QHBoxLayout()
        spins = []
        for label_text, value in zip(("x", "y", "z"), values):
            hbox.addWidget(QLabel(label_text, self))
            spin = QDoubleSpinBox(self)
            spin.setRange(-1.0, 1.0)
            spin.setDecimals(6)
            spin.setSingleStep(0.1)
            spin.setValue(value)
            spin.valueChanged.connect(self._update_plot)
            hbox.addWidget(spin)
            spins.append(spin)
        layout.addLayout(hbox)
        return spins

    def _selected_mode_indices(self):
        """Return selected zero-based mode indices."""
        selected = []
        for row, mode in enumerate(self._modes):
            item = self._mode_table.item(row, 0)
            if item is not None and item.checkState() == Qt.Checked:
                selected.append(mode["index"])
        return selected

    def _select_all_modes(self):
        """Select every Raman-active mode."""
        self._mode_table.blockSignals(True)
        for row in range(self._mode_table.rowCount()):
            self._mode_table.item(row, 0).setCheckState(Qt.Checked)
        self._mode_table.blockSignals(False)
        self._update_plot()

    def _clear_modes(self):
        """Clear every selected mode."""
        self._mode_table.blockSignals(True)
        for row in range(self._mode_table.rowCount()):
            self._mode_table.item(row, 0).setCheckState(Qt.Unchecked)
        self._mode_table.blockSignals(False)
        self._update_plot()

    def _on_mode_item_changed(self, _item):
        """Redraw after a mode selection change."""
        self._update_plot()

    def _vector_from_spins(self, spins):
        """Return a vector from three spin boxes."""
        return np.array([spin.value() for spin in spins], dtype=float)

    def _update_plot(self):
        """Redraw the polar plot."""
        self.figure.clear()
        ax = self.figure.add_subplot(111, projection="polar")
        mode_indices = self._selected_mode_indices()
        if not mode_indices:
            self._raw_max_label.setText("Raw max intensity: 0.0")
            self._orthogonality_label.setText("Orthogonality: n/a")
            ax.set_title("No Raman-active modes selected")
            self.canvas.draw()
            return

        try:
            result = calculate_polar_curves(
                self._raman_tensors,
                mode_indices,
                rotation_axis=self._vector_from_spins(self._rz_spins),
                zero_angle_polarisation=self._vector_from_spins(self._pza_spins),
                angles_rad=self._angles_rad,
                geometries=POLAR_GEOMETRIES,
                normalise=True,
            )
        except ValueError as exc:
            self._raw_max_label.setText("Raw max intensity: n/a")
            self._orthogonality_label.setText(str(exc))
            ax.set_title("Invalid polarisation basis")
            self.canvas.draw()
            return

        colors = {"VV": "tab:red", "HV": "tab:blue", "theta": "tab:green"}
        for geometry in POLAR_GEOMETRIES:
            ax.plot(result["angles_rad"], result["curves"][geometry], label=geometry, color=colors[geometry])
        raw_max = result["raw_max_intensity"]
        self._raw_max_label.setText(f"Raw max intensity: {raw_max:.6g}")
        self._orthogonality_label.setText(f"Orthogonality: {result['orthogonality']:.6f}")
        if abs(result["orthogonality"]) > 1.0e-6:
            self._orthogonality_label.setText(self._orthogonality_label.text() + " (Pza adjusted for plot)")
        ax.set_title(f"Polar Raman intensities (max: {raw_max:.6g})")
        ax.set_ylim(0.0, 1.05)
        ax.legend(loc="upper right", bbox_to_anchor=(1.18, 1.12))
        self.canvas.draw()

    def _show_acknowledgement(self):
        """Show a compact citation acknowledgement."""
        QMessageBox.information(
            self,
            "Acknowledgement",
            "This polar Raman visualisation is inspired by the CRD/Oulu Raman simulator.\n\n"
            "Please cite:\n"
            'M. Bagheri and H.-P. Komsa, "High-throughput computation of Raman spectra '
            'from first principles", Scientific Data 10, 80 (2023).\n'
            "https://doi.org/10.1038/s41597-023-01988-5",
        )
