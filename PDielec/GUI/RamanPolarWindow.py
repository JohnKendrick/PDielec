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
    QCheckBox,
    QDoubleSpinBox,
    QGridLayout,
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

from PDielec.RamanPolarCalculator import (
    POLAR_GEOMETRIES,
    apply_eo_correction,
    calculate_polar_curves,
)

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
        self._activity_display = polar_data.get("raman_activity_display", {
            "factor": 1.0,
            "unit_label": "Å/amu",
            "tooltip": "Displayed in the internal R_epsilon convention.",
        })

        # EO correction data (None when χ^(2) is not available).
        self._chi2 = polar_data.get("chi2_repsilon", polar_data.get("chi2_pm_per_v"))
        self._Z_mat = polar_data.get("Z_mat")
        self._eigvecs = polar_data.get("eigvecs")
        self._eps_inf = polar_data.get("eps_inf")
        self._has_eo = self._chi2 is not None

        # LO frequency data (None when NAC data is unavailable).
        self._hessian        = polar_data.get("hessian")
        self._U_TO           = polar_data.get("U_TO")
        self._volume_au      = polar_data.get("volume_au")
        self._volume_angstrom = polar_data.get("volume_angstrom")
        self._has_lo = (self._hessian is not None and self._Z_mat is not None)

        main_hbox = QHBoxLayout(self)
        controls_vbox = QVBoxLayout()
        plot_vbox = QVBoxLayout()

        self._mode_table = QTableWidget(self)
        self._mode_table.setColumnCount(5)
        self._mode_table.setHorizontalHeaderLabels(
            [
                "Include",
                "Mode",
                "TO freq (cm-1)",
                "LO freq (cm-1)",
                f"Raman total\n({self._activity_display['unit_label']})",
            ]
        )
        self._mode_table.horizontalHeaderItem(0).setToolTip("Check to include this mode in the polar plot.")
        self._mode_table.horizontalHeaderItem(1).setToolTip("Mode number.")
        self._mode_table.horizontalHeaderItem(2).setToolTip("Transverse optical (TO) frequency in cm⁻¹.")
        self._mode_table.horizontalHeaderItem(3).setToolTip(
            "Longitudinal optical (LO) frequency in cm⁻¹\nfor the current q̂ direction."
        )
        self._mode_table.horizontalHeaderItem(4).setToolTip(
            f"Total Raman scattering activity. {self._activity_display['tooltip']}"
        )
        self._mode_table.verticalHeader().setVisible(False)
        self._mode_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._mode_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._mode_table.setRowCount(len(self._modes))
        self._mode_table.itemChanged.connect(self._on_mode_item_changed)
        self._populate_mode_table()
        self._mode_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self._mode_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self._mode_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self._mode_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self._mode_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self._mode_table.setMinimumWidth(400)

        mode_label = QLabel("Raman-active modes", self)
        controls_vbox.addWidget(mode_label)
        controls_vbox.addWidget(self._mode_table)

        buttons_hbox = QHBoxLayout()
        self._select_all_btn = QPushButton("Select all active modes", self)
        self._select_all_btn.setToolTip("Select all Raman-active modes for inclusion in the polar plot.")
        self._select_all_btn.clicked.connect(self._select_all_modes)
        self._clear_btn = QPushButton("Clear selection", self)
        self._clear_btn.setToolTip("Deselect all modes.")
        self._clear_btn.clicked.connect(self._clear_modes)
        buttons_hbox.addWidget(self._select_all_btn)
        buttons_hbox.addWidget(self._clear_btn)
        controls_vbox.addLayout(buttons_hbox)

        # All three direction rows share one grid so x/y/z columns line up.
        vec_grid = QGridLayout()
        vec_grid.setColumnStretch(2, 1)
        vec_grid.setColumnStretch(4, 1)
        vec_grid.setColumnStretch(6, 1)
        controls_vbox.addLayout(vec_grid)

        rz_tip = ("The light propagation direction (in and out of the sample),\n"
                  "or equivalently the crystal rotation axis.")
        self._rz_spins, _ = self._add_vector_controls(
            vec_grid, 0, (0.0, 0.0, 1.0), prefix_label="Light direction Rz:", tooltip=rz_tip
        )

        pza_tip = ("The incoming laser polarisation direction at rotation angle 0°.\n"
                   "Should be perpendicular to Rz.\n"
                   "When Rz changes, Pza is automatically made perpendicular to Rz.")
        self._pza_spins, _ = self._add_vector_controls(
            vec_grid, 1, (1.0, 0.0, 0.0), prefix_label="Zero-angle pol. Pza:", tooltip=pza_tip
        )

        q_init = polar_data.get("q_hat_init", (0.0, 0.0, 1.0))
        q_tip = ("The phonon wavevector direction in the crystal frame.\n"
                 "Used for LO frequency and electro-optic correction calculations.")
        self._q_spins, self._q_row_widgets = self._add_vector_controls(
            vec_grid, 2, q_init, prefix_label="Phonon q̂ direction:", tooltip=q_tip
        )

        # Hide the LO freq column when no NAC data is available.
        if not self._has_lo:
            self._mode_table.setColumnHidden(3, True)

        self._orthogonality_label = QLabel("Orthogonality: 0.000000", self)
        self._orthogonality_label.setToolTip(
            "Dot product Rz · Pza.\n"
            "Should be zero for physically meaningful results.\n"
            "A non-zero value means Pza is not perpendicular to Rz."
        )
        controls_vbox.addWidget(self._orthogonality_label)

        # EO correction controls — only shown when χ^(2) is available.
        self._eo_checkbox = QCheckBox("Include EO correction (χ⁽²⁾ available)", self)
        self._eo_checkbox.setChecked(True)
        self._eo_checkbox.stateChanged.connect(self._update_plot)
        self._eo_checkbox.setVisible(self._has_eo)
        self._eo_checkbox.setToolTip(
            "Apply electro-optic correction to the Raman tensors\n"
            "using the χ⁽²⁾ tensor (requires CHI2.DAT data)."
        )
        controls_vbox.addWidget(self._eo_checkbox)

        show_q_controls = self._has_eo or self._has_lo
        for w in self._q_row_widgets:
            w.setVisible(show_q_controls)
        self._sync_q_btn = QPushButton("Sync q̂ to Rz", self)
        self._sync_q_btn.setToolTip(
            "Copy the light direction Rz into q̂.\n"
            "Appropriate for backscattering geometry where q̂ ≈ Rz."
        )
        self._sync_q_btn.clicked.connect(self._sync_q_to_rz)
        self._sync_q_btn.setVisible(show_q_controls)
        controls_vbox.addWidget(self._sync_q_btn)

        self._desc_btn = QPushButton("Description", self)
        self._desc_btn.setToolTip("Show a description of the coordinate system and plot geometries.")
        self._desc_btn.clicked.connect(self._show_description)
        controls_vbox.addWidget(self._desc_btn)
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

    def _display_activity(self, activity):
        """Convert an R_epsilon activity to the selected display units."""
        return float(activity) * float(self._activity_display.get("factor", 1.0))

    def set_activity_display(self, activity_display):
        """Update Raman activity display units and refresh table/title text."""
        self._activity_display = activity_display
        header_item = self._mode_table.horizontalHeaderItem(4)
        if header_item is not None:
            header_item.setText(f"Raman total\n({self._activity_display['unit_label']})")
            header_item.setToolTip(f"Total Raman scattering activity. {self._activity_display['tooltip']}")
        self._update_plot()

    def _populate_mode_table(self):
        """Populate the mode table."""
        self._mode_table.blockSignals(True)
        ro_flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        for row, mode in enumerate(self._modes):
            include_item = QTableWidgetItem("")
            include_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            include_item.setCheckState(Qt.Checked if mode["selected"] else Qt.Unchecked)
            self._mode_table.setItem(row, 0, include_item)

            mode_item     = QTableWidgetItem(str(mode["display_mode"]))
            freq_item     = QTableWidgetItem(f"{mode['frequency_cm1']:.4f}")
            lo_item       = QTableWidgetItem("—")
            activity_item = QTableWidgetItem(f"{self._display_activity(mode['raman_total']):.4f}")
            for item in (mode_item, freq_item, lo_item, activity_item):
                item.setFlags(ro_flags)
            for item in (freq_item, lo_item, activity_item):
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self._mode_table.setItem(row, 1, mode_item)
            self._mode_table.setItem(row, 2, freq_item)
            self._mode_table.setItem(row, 3, lo_item)
            self._mode_table.setItem(row, 4, activity_item)
        self._mode_table.blockSignals(False)

    def _add_vector_controls(self, grid, row, values, prefix_label=None, tooltip=None):
        """Add x/y/z spin boxes into one row of a shared QGridLayout.

        Placing all three direction rows in the same grid ensures the x/y/z
        spin-box columns are aligned across rows.

        Parameters
        ----------
        grid : QGridLayout
            Shared grid for all vector rows.
        row : int
            Grid row index.
        values : sequence of float
            Initial (x, y, z) values.
        prefix_label : str, optional
            Text for the label placed in column 0.
        tooltip : str, optional
            Tooltip applied to every widget in the row.

        Returns
        -------
        tuple[list[QDoubleSpinBox], list[QWidget]]
            The three spin boxes and all widgets in the row (prefix label first
            if present, then alternating axis-label / spin-box for x, y, z).
        """
        row_widgets = []
        col = 0
        if prefix_label is not None:
            lbl_widget = QLabel(prefix_label, self)
            if tooltip:
                lbl_widget.setToolTip(tooltip)
            grid.addWidget(lbl_widget, row, col)
            row_widgets.append(lbl_widget)
        col += 1  # column 0 is always the prefix column; skip it even when absent
        spins = []
        for label_text, value in zip(("x:", "y:", "z:"), values):
            axis_lbl = QLabel(label_text, self)
            if tooltip:
                axis_lbl.setToolTip(tooltip)
            grid.addWidget(axis_lbl, row, col)
            row_widgets.append(axis_lbl)
            col += 1
            spin = QDoubleSpinBox(self)
            spin.setRange(-1.0, 1.0)
            spin.setDecimals(6)
            spin.setSingleStep(1.0)
            spin.setValue(value)
            spin.valueChanged.connect(self._update_plot)
            if tooltip:
                spin.setToolTip(tooltip)
            grid.addWidget(spin, row, col)
            row_widgets.append(spin)
            col += 1
            spins.append(spin)
        return spins, row_widgets

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

    def _sync_q_to_rz(self):
        """Copy the current Rz vector into the q̂ spinboxes."""
        for q_spin, rz_spin in zip(self._q_spins, self._rz_spins):
            q_spin.blockSignals(True)
            q_spin.setValue(rz_spin.value())
            q_spin.blockSignals(False)
        self._update_plot()

    def _update_lo_table(self, q_hat=None):
        """Refresh the LO freq column in the mode table for the given q̂.

        Parameters
        ----------
        q_hat : ndarray, shape (3,) or None
            Normalised phonon wavevector.  When *None* the column is filled
            with ``"—"``.

        """
        if not self._has_lo:
            return
        self._mode_table.blockSignals(True)
        if q_hat is None:
            for row in range(self._mode_table.rowCount()):
                item = self._mode_table.item(row, 3)
                if item is not None:
                    item.setText("—")
        else:
            from PDielec.RamanPolarCalculator import compute_lo_frequencies
            lo_freqs = compute_lo_frequencies(
                q_hat, self._hessian, self._U_TO,
                self._Z_mat, self._eps_inf, self._volume_au,
            )
            for row, mode in enumerate(self._modes):
                item = self._mode_table.item(row, 3)
                if item is not None:
                    idx = mode["index"]
                    item.setText(f"{lo_freqs[idx]:.4f}" if idx < len(lo_freqs) else "—")
        self._mode_table.blockSignals(False)

    def _vector_from_spins(self, spins):
        """Return a vector from three spin boxes."""
        return np.array([spin.value() for spin in spins], dtype=float)

    def _update_plot(self):
        """Redraw the polar plot."""
        self.figure.clear()
        ax = self.figure.add_subplot(111, projection="polar")
        mode_indices = self._selected_mode_indices()
        if not mode_indices:
            self._orthogonality_label.setText("Orthogonality: n/a")
            ax.set_title("No Raman-active modes selected")
            self.canvas.draw()
            return

        # Determine q̂ (used for both EO correction and LO frequencies).
        q_raw = self._vector_from_spins(self._q_spins)
        q_norm = float(np.linalg.norm(q_raw))
        q_hat = (q_raw / q_norm) if q_norm > 1e-12 else None

        # Update LO freq column (independent of EO).
        self._update_lo_table(q_hat)

        # Apply EO correction when χ^(2) is available and the checkbox is on.
        tensors = self._raman_tensors
        eo_active = False
        if self._has_eo and self._eo_checkbox.isChecked() and q_hat is not None:
            try:
                tensors = apply_eo_correction(
                    self._raman_tensors, self._chi2, q_hat,
                    self._Z_mat, self._eigvecs, self._eps_inf,
                )
                eo_active = True
            except Exception as exc:
                logger.warning("EO correction failed: %s", exc)

        # Update Raman total column with EO-corrected values (or restore TO values).
        self._mode_table.blockSignals(True)
        if eo_active and self._volume_angstrom is not None:
            from PDielec import Calculator
            acts = Calculator.raman_intensities(tensors, self._volume_angstrom)
            for row, mode in enumerate(self._modes):
                item = self._mode_table.item(row, 4)
                if item is not None:
                    idx = mode["index"]
                    item.setText(f"{self._display_activity(acts[idx, 0]):.4f}" if idx < len(acts) else "—")
        else:
            for row, mode in enumerate(self._modes):
                item = self._mode_table.item(row, 4)
                if item is not None:
                    item.setText(f"{self._display_activity(mode['raman_total']):.4f}")
        self._mode_table.blockSignals(False)

        try:
            result = calculate_polar_curves(
                tensors,
                mode_indices,
                rotation_axis=self._vector_from_spins(self._rz_spins),
                zero_angle_polarisation=self._vector_from_spins(self._pza_spins),
                angles_rad=self._angles_rad,
                geometries=POLAR_GEOMETRIES,
                normalise=True,
            )
        except ValueError as exc:
            self._orthogonality_label.setText(str(exc))
            ax.set_title("Invalid polarisation basis")
            self.canvas.draw()
            return

        colors = {"VV": "tab:red", "HV": "tab:blue", "Theta": "tab:green"}
        for geometry in POLAR_GEOMETRIES:
            ax.plot(result["angles_rad"], result["curves"][geometry], label=geometry, color=colors[geometry])
        raw_max = result["raw_max_intensity"]
        self._orthogonality_label.setText(f"Orthogonality: {result['orthogonality']:.6f}")
        if abs(result["orthogonality"]) > 1.0e-6:
            self._orthogonality_label.setText(self._orthogonality_label.text() + " (Pza adjusted for plot)")
        eo_suffix = " +EO" if eo_active else ""
        display_max = self._display_activity(raw_max)
        unit_label = self._activity_display["unit_label"]
        ax.set_title(f"Polar Raman intensities{eo_suffix}\n(max: {display_max:.4f} {unit_label})")
        ax.set_ylim(0.0, 1.05)
        ax.legend(loc="upper right", bbox_to_anchor=(1.05, 1.12))
        self.canvas.draw()

    def _show_description(self):
        """Show a description of the coordinate system and plot geometries."""
        QMessageBox.information(
            self,
            "Polar Raman Description",
            "COORDINATE SYSTEM\n"
            "─────────────────\n"
            "Three directions in crystal coordinates define the geometry:\n\n"
            "1. The light propagation direction Rz.\n\n"
            "2. The direction Pz showing the laser polarisation direction at zero-angle.\n"
            "   This should be perpendicular to Rz. \n"
            "   When Rz is changed orthogonality of Pz is ensured automatically. \n"
            "   Pza can also be changed, but orthogonality is not enforces\n"
            "   'Orthogonality' displays the orthogonality between Pz and Rz\n"
            "3. The remaining orthogonal direct Pp = Rz × Pz\n\n"
            "PLOT LEGEND\n"
            "───────────────\n"
            "The sample is rotated by angle θ (via rotation matrix M):\n\n"
            "• VV (vertical-vertical / parallel)\n"
            "    Both incident and scattered polarisations rotate together.\n"
            "    E_L = M·Pz,  E_S = M·Pz\n\n"
            "• HV (horizontal-vertical / perpendicular)\n"
            "    Incident and scattered polarisations are always perpendicular;\n"
            "    both rotate with the sample — equivalent to simultaneous\n"
            "    rotation around Pz and Pp.\n"
            "    E_L = M·Pza,  E_S = M·Pp\n\n"
            "• Theta (fixed incident, rotating analyser)\n"
            "    Incident polarisation is fixed; the analyser rotates.\n"
            "    E_L = Pz  (fixed),  E_S = M·Pz\n\n"
            "REFERENCE\n"
            "─────────\n"
            "This visualisation is inspired by the CRD/Oulu Raman simulator.\n"
            "https://ramandb.oulu.fi/simulator \n"
            "Please cite:\n"
            "M. Bagheri and H.-P. Komsa, \n "
            "\"High-throughput computation of Raman spectra from first principles\"\n"
            " Scientific Data 10, 80 (2023).\n"
            "https://doi.org/10.1038/s41597-023-01988-5",
        )
