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
"""SettingsTab module."""
# -*- coding: utf8 -*-

import logging

import numpy as np
from qtpy.QtCore import QCoreApplication, QSize, Qt
from qtpy.QtGui import QKeySequence
from qtpy.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from PDielec import Calculator, DielectricFunction
from PDielec.Constants import amu, angs2bohr, angstrom, average_masses, isotope_masses, wavenumber
from PDielec.RamanPolarCalculator import raman_active_mode_indices

logger = logging.getLogger(__name__)
class FixedQTableWidget(QTableWidget):
    """A custom QTableWidget with fixed number of rows and columns.

    This widget is designed to have a fixed size based on a specified number of rows and columns, rather than adjusting
    dynamically to the content it displays. If the number of rows and columns isn't specifically provided, it defaults
    to the current count of rows and columns in the widget.

    Parameters
    ----------
    args
        Variable length argument list.
    parent : Optional[QWidget]
        The parent widget. Default is None.
    rows : Optional[int]
        The fixed number of rows for the table. If None, the table's row count is used. Default is None.
    columns : Optional[int]
        The fixed number of columns for the table. If None, the table's column count is used. Default is None.

    Methods
    -------
    sizeHint() Calculates and returns the recommended size for the table widget based on the fixed number of rows and
    columns, or based on its current row and column counts if no fixed numbers are provided.

    See Also
    --------
    QTableWidget : The base class from which this class is derived.

    """

    def __init__(self, *args, parent=None, rows=None, columns=None):
        """Initialize a QTableWidget instance.

        Parameters
        ----------
        args
            Variable length argument list passed to the super class initializer.
        parent : optional
            The parent widget of the QTableWidget. Default is None.
        rows : int, optional
            The initial number of rows in the table. Default is None.
        columns : int, optional
            The initial number of columns in the table. Default is None.

        """        
        self.columns = columns
        self.rows = rows
        super(QTableWidget, self).__init__(*args)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectItems)

    def keyPressEvent(self, event):
        """Copy and paste rectangular table selections with the system clipboard."""
        if event.matches(QKeySequence.Copy):
            self.copy_selection_to_clipboard()
            event.accept()
            return
        if event.matches(QKeySequence.Paste):
            self.paste_clipboard_from_selection()
            event.accept()
            return
        super().keyPressEvent(event)

    def copy_selection_to_clipboard(self):
        """Copy selected cells as tab-separated text."""
        selected = {(index.row(), index.column()) for index in self.selectedIndexes() if index.isValid()}
        if not selected:
            return
        min_row = min(row for row, _col in selected)
        max_row = max(row for row, _col in selected)
        min_col = min(col for _row, col in selected)
        max_col = max(col for _row, col in selected)

        rows = []
        for row in range(min_row, max_row + 1):
            values = []
            for col in range(min_col, max_col + 1):
                if (row, col) in selected:
                    item = self.item(row, col)
                    values.append("" if item is None else item.text())
                else:
                    values.append("")
            rows.append("\t".join(values))
        QApplication.clipboard().setText("\n".join(rows))

    def paste_clipboard_from_selection(self):
        """Paste tab-separated text into editable cells starting at the selection."""
        text = QApplication.clipboard().text()
        if not text:
            return

        lines = text.splitlines()
        values = [line.split("\t") for line in lines]
        if not values:
            return

        ranges = self.selectedRanges()
        if ranges:
            start_row = min(selection.topRow() for selection in ranges)
            start_col = min(selection.leftColumn() for selection in ranges)
        else:
            start_row = self.currentRow() if self.currentRow() >= 0 else 0
            start_col = self.currentColumn() if self.currentColumn() >= 0 else 0

        for row_offset, row_values in enumerate(values):
            row = start_row + row_offset
            if row >= self.rowCount():
                break
            for col_offset, value in enumerate(row_values):
                col = start_col + col_offset
                if col >= self.columnCount():
                    break
                item = self.item(row, col)
                if item is None:
                    item = QTableWidgetItem(value)
                    item.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    self.setItem(row, col, item)
                elif item.flags() & Qt.ItemIsEnabled and item.flags() & Qt.ItemIsEditable:
                    item.setText(value)

    def sizeHint(self):
        """Calculate the suggested dimensions for the widget based on its content.

        This method calculates the preferred size of the widget by summing up the widths of all columns and the heights
        of all rows, including any scrollbars, headers, and frame widths.

        Parameters
        ----------
        None

        Returns
        -------
        QSize The calculated size that suggests the optimal dimensions for the widget. The width is determined by the
        total width of all columns, the vertical header width, the vertical scrollbar width, and twice the frame width.
        The height is determined by the total height of all rows, the horizontal header height, the horizontal scrollbar
        height, and twice the frame width.

        """        
        width = 0
        columns = self.columnCount() if self.columns is None else self.columns
        for i in range(columns):
            width += self.columnWidth(i)
        width += self.verticalHeader().sizeHint().width()
        width += self.verticalScrollBar().sizeHint().width()
        width += self.frameWidth()*2
        height = 0
        rows = self.rowCount() if self.rows is None else self.rows
        for i in range(rows):
            height += self.rowHeight(i)
        height += self.verticalHeader().sizeHint().width()
        height += self.horizontalScrollBar().sizeHint().height()
        height += self.frameWidth()*2
        return QSize(width,height)

class SettingsTab(QWidget):
    """A class for managing and displaying settings related to optical permittivity, mass calculation, and vibrational calculations for crystal structures within a GUI application.

    This class allows for the configuration of various simulation parameters, such as Eckart conditions, Born charge
    neutrality, atomic mass definitions, and optical permittivity values. It also enables the calculation of infrared
    intensities, dielectric functions, and permittivity values based on these settings.

    Parameters
    ----------
    parent : QWidget
        The parent widget (container) of this widget.
    debug : bool, optional
        Flag indicating whether debugging messages should be printed. Defaults to False.

    Attributes
    ----------
    notebook : parent The notebook widget that contains this settings tab. refresh_required : bool A flag indicating
    whether the displayed information needs to be refreshed. calculation_required : bool A flag indicating whether new
    calculations are needed based on the current settings. settings : dict A dictionary of settings related to
    calculations. Includes flags for Eckart conditions, Born charge neutrality, mass definition, and sigma values.
    mass_definition_options : list A list of options for defining atomic masses. masses_dictionary : dict A dictionary
    where keys are element symbols and values are the corresponding atomic masses. modes_selected : list A list
    indicating which vibrational modes are selected for calculations. frequencies_cm1 : list A list of vibrational
    frequencies in cm^-1. frequencies_have_been_edited : bool A flag indicating if frequencies have been manually
    edited. intensities : list A list of calculated infrared intensities for the vibrational modes. sigmas_cm1 : list A
    list of Lorentzian width factors (σ) in cm^-1 for the vibrational modes. oscillator_strengths : list A list of
    oscillator strengths calculated for the vibrational modes. mass_weighted_normal_modes : np.ndarray An array of
    mass-weighted normal modes. CrystalPermittivityObject : DielectricFunction An object for calculating dielectric
    functions. vs_cm1 : list A list of frequencies in cm^-1 at which the crystal permittivity is calculated.
    crystal_permittivity : list A list of calculated crystal permittivity values. recalculate_selected_modes : bool A
    flag indicating whether selected modes should be recalculated. reader : object A reader object from the parent
    notebook, used for accessing and manipulating input data.


    Methods
    -------
    set_element_mass(element, mass)
        Sets the mass for a specified element.
    create_intensity_table()
        Creates a table displaying calculated intensities for selected vibrational modes.
    request_refresh()
        Requests a refresh of displayed data.
    write_spreadsheet()
        Writes the current settings and calculated values to a spreadsheet.
    redraw_output_tw()
        Redraws the output table widget with updated values.
    on_sigma_changed()
        Handles changes to the sigma value.
    on_mass_cb_activated(index)
        Activates the specified mass definition option based on the selected index.
    set_masses_tw()
        Sets the masses in the table widget based on the current mass definition setting.
    on_output_tw_itemChanged(item)
        Handles changes to items in the output table widget.
    on_element_masses_tw_itemClicked(item)
        Handles item clicks in the element masses table widget.
    on_element_masses_tw_itemChanged(item)
        Handles changes to items in the element masses table widget.
    on_optical_tw_itemChanged(item)
        Handles changes to items in the optical permittivity table widget.
    on_optical_tw_itemClicked(item)
        Handles item clicks in the optical permittivity table widget.
    refresh(force=False)
        Refreshes the displayed information and recalculates values if necessary.
    refresh_optical_permittivity_tw()
        Refreshes the displayed optical permittivity values in the table widget.
    set_optical_permittivity_tw()
        Sets the optical permittivity values based on reader data.
    on_born_changed()
        Handles changes to the Born charge neutrality setting.
    on_eckart_changed()
        Handles changes to the Eckart condition setting.
    on_spectroscopy_type_cb_activated(index)
        Handles changes to the spectroscopy type
    calculate(vs_cm1)
        Calculates the permittivity of the crystal over the specified range of frequencies.
    get_crystal_permittivity(vs_cm1)
        Returns the crystal permittivity for the specified frequencies.
    get_crystal_permittivity_object()
        Returns the crystal permittivity object used for calculations.

    """

    def __init__(self, parent, debug=False):
        """Initialize the QWidget object with specific settings for an application.

        Parameters
        ----------
        parent : QWidget
            The parent widget (container) of this widget.
        debug : bool, optional
            Flag indicating whether debugging messages should be printed. Defaults to False.

        Notes
        -----
        This method sets up various elements of the interface, including checkboxes, comboBoxes, spinBoxes, and
        tableWidgets with predefined settings and connects various signals and slots. Settings related to the Eckart
        flag, neutral Born charges, sigma values, mass definition, and optical permittivity are initialized.
        Additionally, UI elements for adjusting these settings are created. It involves initializing debugging
        capabilities, setting up a notebook reference, and various flags and settings necessary for the functionality of
        the application.

        """        
        super(QWidget, self).__init__(parent)
        logger.debug("Start:: initialisation")
        self.notebook = parent
        self.refresh_required = True
        self.calculation_required = False
        self.settings = {}
        self.settings["Eckart flag"] = True
        self.settings["Neutral Born charges"] = False
        self.settings["Sigma value"] = 5
        self.mass_definition_options = ["average","program","isotope","gui"]
        self.settings["Mass definition"] = "average"
        self.settings["Optical permittivity"] = None
        self.settings["Optical permittivity edited"] = False
        self.settings["Symmetrise optical permittivity"] = True
        self.settings["Spectroscopy type"] = "Powder Infrared"
        self.settings["Raman activity units"] = "polarizability"
        self.spectroscopy_types = ["Powder Infrared", "Powder ATR", "Powder Raman", "Crystal Infrared", "Crystal Raman"]
        self.raman_activity_unit_options = [
            ("polarizability", "Polarizability volume (Å⁴/amu)"),
            ("epsilon", "R_epsilon (Å/amu)"),
        ]
        self.masses_dictionary = {}
        self.modes_selected = []
        self.frequencies_cm1 = []
        self.frequencies_have_been_edited = False
        self.intensities = []
        self.raman_intensities = []
        self.raman_intensities_par = []
        self.raman_intensities_perp = []
        self.sigmas_cm1 = []
        self.oscillator_strengths = []
        self.modal_pair_use_nac = []
        self.mass_weighted_normal_modes = None
        self.CrystalPermittivityObject = None
        self.vs_cm1 = []
        self.crystal_permittivity = []
        self.recalculate_selected_modes = True
        #
        # get the reader from the main tab
        self.reader = self.notebook.reader
        #
        # Create second tab - SettingsTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # Spectroscopy type selection
        self.spectroscopy_type_cb = QComboBox(self)
        self.spectroscopy_type_cb.setToolTip("Choose the spectroscopy type for all scenarios")
        self.spectroscopy_type_cb.addItems(self.spectroscopy_types)
        self.spectroscopy_type_cb.activated.connect(self.on_spectroscopy_type_cb_activated)
        form.addRow(QLabel("Spectroscopy type:", self), self.spectroscopy_type_cb)
        #
        # The eckart checkbox
        self.eckart_cb = QCheckBox(self)
        self.eckart_cb.setToolTip("Applying Eckart conditions ensures three zero translation mode)")
        self.eckart_cb.setText("")
        self.eckart_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings["Eckart flag"]:
            self.eckart_cb.setCheckState(Qt.Checked)
        else:
            self.eckart_cb.setCheckState(Qt.Unchecked)
        self.eckart_cb.stateChanged.connect(self.on_eckart_changed)
        form.addRow(QLabel("Apply Eckart conditions?", self), self.eckart_cb)
        #
        # Add the Born neutral condition
        self.born_cb = QCheckBox(self)
        self.born_cb.setToolTip("Applying Born charge neutrality ensures unit cell has zero charge")
        self.born_cb.setText("")
        self.born_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings["Neutral Born charges"]:
            self.born_cb.setCheckState(Qt.Checked)
        else:
            self.born_cb.setCheckState(Qt.Unchecked)
        self.born_cb.stateChanged.connect(self.on_born_changed)
        form.addRow(QLabel("Apply Born charge neutrality?",self),self.born_cb)
        #
        # The mass definition combobox
        self.mass_cb = QComboBox(self)
        self.mass_cb.setToolTip("The atomic masses used to calculate frequencies and intensities can be give here")
        self.mass_cb.addItem("Average natural abundance")
        self.mass_cb.addItem("Mass taken from QM/MM program")
        self.mass_cb.addItem("Most common isotope mass")
        self.mass_cb.addItem("Masses set individually")
        self.mass_cb.model().item(3).setEnabled(False)
        #
        # set default to average natural abundance
        self.mass_cb.activated.connect(self.on_mass_cb_activated)
        self.current_mass_definition_index = self.mass_definition_options.index(self.settings["Mass definition"])
        self.mass_cb.setCurrentIndex(0)
        form.addRow(QLabel("Atomic mass defintion:", self), self.mass_cb)
        #
        # Create Table containing the masses - block signals until the table is loaded
        self.element_masses_tw = FixedQTableWidget(parent=self)
        self.element_masses_tw.setToolTip("Individual element masses can be modified here")
        self.element_masses_tw.itemClicked.connect(self.on_element_masses_tw_itemClicked)
        self.element_masses_tw.itemChanged.connect(self.on_element_masses_tw_itemChanged)
        self.element_masses_tw.setRowCount(1)
        self.element_masses_tw.blockSignals(True)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        self.element_masses_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel("Atomic masses", self), self.element_masses_tw)
        #
        # The lorentzian width - sigma
        self.sigma_sb = QDoubleSpinBox(self)
        self.sigma_sb.setRange(0.1,100.0)
        self.sigma_sb.setSingleStep(0.1)
        self.sigma_sb.setDecimals(1)
        self.sigma_sb.setValue(self.settings["Sigma value"])
        self.sigma_sb.setToolTip("Set the default Lorentzian width factor (σ) in cm-1")
        self.sigma_sb.valueChanged.connect(self.on_sigma_changed)
        label = QLabel("Lorentzian width (σ):", self)
        label.setToolTip("Set the default Lorentzian width factor (σ) in cm-1")
        form.addRow(label, self.sigma_sb)
        #
        # Raman activity display units
        self.raman_activity_units_cb = QComboBox(self)
        self.raman_activity_units_cb.setToolTip(
            "Choose the units used to display Raman activities in this table and the Raman polar plot."
        )
        for _key, label_text in self.raman_activity_unit_options:
            self.raman_activity_units_cb.addItem(label_text)
        self.raman_activity_units_cb.activated.connect(self.on_raman_activity_units_cb_activated)
        self.raman_activity_units_label = QLabel("Raman activity units:", self)
        self.raman_activity_units_label.setToolTip(self.raman_activity_units_cb.toolTip())
        form.addRow(self.raman_activity_units_label, self.raman_activity_units_cb)
        #
        # Create the Optical permittivity table widget and block signals until a click on the widget
        self.optical_tw = FixedQTableWidget(3,3,parent=self)
        self.optical_tw.setToolTip("The optical permittivity is taken from the calculation where this is possible.  If it is not availble suitbale values should be provided here")
        #
        # Set the header names
        self.optical_tw.setHorizontalHeaderLabels(["x","y","z"])
        self.optical_tw.setVerticalHeaderLabels  (["x","y","z"])
        self.optical_tw.itemClicked.connect(self.on_optical_tw_itemClicked)
        self.optical_tw.itemChanged.connect(self.on_optical_tw_itemChanged)
        #
        # Block the widget until the optical permittivity is loaded
        self.optical_tw.blockSignals(True)
        self.optical_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel("Optical permittivity:", self), self.optical_tw)
        self.symmetrise_optical_cb = QCheckBox(self)
        self.symmetrise_optical_cb.setToolTip(
            "When checked, the optical permittivity tensor loaded from the DFT output is\n"
            "averaged over all point-group symmetry operations of the crystal:\n"
            "  ε_sym[i,j] = (1/N) Σ_R  R[i,a] R[j,b] ε[a,b]\n"
            "This projects the raw tensor onto the invariant subspace of the crystal\n"
            "symmetry, removing numerical noise that violates the lattice symmetry\n"
            "(e.g. spurious in-plane birefringence in a uniaxial crystal).\n"
            "Disable only if you need to use the raw, unsymmetrised DFT tensor."
        )
        self.symmetrise_optical_cb.setText("")
        self.symmetrise_optical_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings["Symmetrise optical permittivity"]:
            self.symmetrise_optical_cb.setCheckState(Qt.Checked)
        else:
            self.symmetrise_optical_cb.setCheckState(Qt.Unchecked)
        self.symmetrise_optical_cb.stateChanged.connect(self.on_symmetrise_optical_changed)
        form.addRow(QLabel("Symmetrise optical permittivity?", self), self.symmetrise_optical_cb)
        vbox.addLayout(form)
        #
        # output window
        # Create Table containing the IR active modes
        self.output_tw = FixedQTableWidget(parent=self)
        self.output_tw.setToolTip("Output showing the frequencies and strengths of the IR active modes only")
        self.output_tw.itemChanged.connect(self.on_output_tw_itemChanged)
        self.output_tw.setRowCount(1)
        self.output_tw.blockSignals(True)
        vbox.addWidget(self.output_tw)
        #
        # NAC Options row (always visible — applies to both IR and Raman)
        self._lo_q_label = QLabel("NAC Options:", self)
        self._lo_q_label.setToolTip(
            "Non-analytical correction (NAC) options.\n"
            "When enabled, LO frequencies are computed from the dynamical matrix\n"
            "with the long-range Coulomb correction at the given phonon wavevector q̂.\n"
            "Applies to both IR and Raman mode tables."
        )
        lo_q_hbox = QHBoxLayout()
        lo_q_hbox.setContentsMargins(0, 0, 0, 0)
        lo_q_hbox.setSpacing(2)
        self._nac_apply_cb = QCheckBox("Apply NAC", self)
        self._nac_apply_cb.setChecked(False)
        self._nac_apply_cb.setToolTip(
            "Apply non-analytical correction (NAC) to compute LO frequencies.\n"
            "When checked, the LO frequency column is filled for all modes using\n"
            "the phonon wavevector q̂ below.  Works for both IR and Raman spectroscopy."
        )
        self._nac_apply_cb.stateChanged.connect(self._refresh_lo_columns)
        lo_q_hbox.addWidget(self._nac_apply_cb)
        lo_q_hbox.addSpacing(8)
        self._lo_q_spins = []
        q_tooltips = (
            "x-component of the phonon wavevector q̂ in crystal coordinates",
            "y-component of the phonon wavevector q̂ in crystal coordinates",
            "z-component of the phonon wavevector q̂ in crystal coordinates",
        )
        for axis_label, default_val, tip in zip(("x:", "y:", "z:"), (0, 0, 1), q_tooltips):
            lbl = QLabel(axis_label, self)
            lbl.setToolTip(tip)
            lo_q_hbox.addWidget(lbl)
            spin = QSpinBox(self)
            spin.setRange(-100, 100)
            spin.setSingleStep(1)
            spin.setValue(default_val)
            spin.setToolTip(tip)
            spin.valueChanged.connect(self._refresh_lo_columns)
            lo_q_hbox.addWidget(spin)
            lo_q_hbox.addSpacing(8)
            self._lo_q_spins.append(spin)
        self._lo_eo_cb = QCheckBox("Include EO correction (χ⁽²⁾)", self)
        self._lo_eo_cb.setChecked(False)
        self._lo_eo_cb.setToolTip(
            "Apply the electro-optic (EO) correction to Raman tensors using the\n"
            "second-order susceptibility χ⁽²⁾ and the phonon wavevector q̂.\n"
            "Only available for Raman spectroscopy when χ⁽²⁾ data is present (e.g. CRYSTAL)."
        )
        self._lo_eo_cb.stateChanged.connect(self._refresh_lo_columns)
        lo_q_hbox.addWidget(self._lo_eo_cb)
        lo_q_hbox.addStretch(1)
        self._lo_q_widget = QWidget(self)
        self._lo_q_widget.setLayout(lo_q_hbox)
        self._lo_q_widget.setContentsMargins(0, 0, 0, 0)
        form.addRow(self._lo_q_label, self._lo_q_widget)
        # Raman polar mode simulator (Raman only)
        self.raman_polar_button = QPushButton("Raman polar plot", self)
        self.raman_polar_button.setToolTip(
            "Open a tensor-level Raman polar plot simulator for the Raman-active modes.\n"
            "Displays polar intensity curves for VV, HV, and theta geometries."
        )
        self.raman_polar_button.clicked.connect(self.on_raman_polar_button_clicked)
        self.raman_polar_label = QLabel("Raman polar simulator:", self)
        self.raman_polar_label.setToolTip(
            "Open a tensor-level Raman polar plot simulator for the Raman-active modes."
        )
        form.addRow(self.raman_polar_label, self.raman_polar_button)
        self.update_raman_polar_button()
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        logger.debug("Finished:: initialisation")

    def set_element_mass(self,element,mass):
        """Set the mass value of a specific element in the GUI.

        This method updates the mass value for a specified element in the graphical user interface (GUI). It also
        triggers various updates in the GUI to reflect this change, including updating the mass selection and
        recalculating modes if required.

        Parameters
        ----------
        element : str
            The element whose mass is to be set.
        mass : float
            The mass value to set for the specified element.

        Returns
        -------
        None

        Notes
        -----
        - This method assumes that `self.settings`, `self.masses_dict`, `self.mass_cb`, and other related attributes are
          already defined in the class. - It triggers a refresh of the GUI and potentially recalculates selected modes
          based on the new mass settings.

        """        
        logger.debug(f"Start::  set_element_mass {element} {mass}")
        self.settings["Mass definition"] = "gui"
        self.masses_dict[element] = mass
        self.mass_cb.setCurrentIndex(3)
        self.set_masses_tw()
        self.refresh_required = True
        self.recalculate_selected_modes = True
        self.refresh()
        logger.debug(f"Finished::  set_element_mass {element} {mass}")

    def create_intensity_table(self):
        """Generate the intensity table for spectroscopy analysis.

        This method is responsible for creating a table that includes calculated intensities, frequencies, oscillator
        strengths,  and other related parameters required for spectroscopic analysis. It adapts based on reader
        settings, mass definitions,  and whether to recalculate selected modes.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method relies on several external variables and settings, including: - reader settings for neutralizing or
        resetting born charges - mass definition settings to determine how masses are handled in calculations - a check
        on if frequencies have been edited, which affects reading directly from the reader or not - uses `Calculator`
        class methods to compute normal modes, oscillator strengths, infrared intensities, 
        Raman intensities and ionic permittivity - the
        method also determines the dielectric function based on the configuration and recalculates selected modes if
        necessary - it finally populates a provided table widget (`self.output_tw`) with computed spectroscopic
        properties

        See Also
        --------
        - `Calculator.normal_modes`: For calculation of normal modes. - `Calculator.oscillator_strengths`: For
          calculating oscillator strengths based on normal modes and born charges. - `Calculator.infrared_intensities`:
          For the calculation of infrared intensities from oscillator strengths. - `Calculator.ionic_permittivity`: For
          calculating the ionic permittivity based on mode list, oscillator strengths, and frequencies.

        """        
        logger.debug("Start:: create_intensity_table")
        self.reader = self.notebook.reader

        # Only calculate if the reader is set
        if self.reader is None:
            logger.debug("create_intensity_table aborting as now reader available")
            return

        # Deal with Born charges
        if self.settings["Neutral Born charges"]:
            self.reader.neutralise_born_charges()
        else:
            self.reader.reset_born_charges()

        # Set the Eckart flag in the reader
        self.reader.eckart = self.settings["Eckart flag"]

        # Set up the Mass definitions
        mass_dictionary = []
        self.reader.reset_masses()
        if self.settings["Mass definition"] == "average":
            self.reader.change_masses(average_masses, mass_dictionary)
        elif self.settings["Mass definition"] == "program":
            pass
        elif self.settings["Mass definition"] == "isotope":
            self.reader.change_masses(isotope_masses, mass_dictionary)
        elif self.settings["Mass definition"] == "gui":
            self.reader.change_masses(self.masses_dictionary, mass_dictionary)
        else:
            logger.error(f"Error unknown mass definition {self.settings['Mass definition']}")
        QCoreApplication.processEvents()

        # Calculate normal modes
        self.mass_weighted_normal_modes = self.reader.calculate_mass_weighted_normal_modes()
        if not isinstance(self.mass_weighted_normal_modes, np.ndarray) and not self.mass_weighted_normal_modes:
            logger.warning("create_intensity_table: no normal modes available (no phonon data in input file)")
            return

        # Find the sigmas for each frequency, the sigmas might have been edited by hand
        if not self.frequencies_have_been_edited:
            self.frequencies_cm1 = self.reader.frequencies
        if len(self.sigmas_cm1) == 0:
            self.sigmas_cm1 = [ self.settings["Sigma value"] for i in self.frequencies_cm1 ]

        # Get the Born charges from the reader
        born_charges = np.array(self.reader.born_charges)

        # If reader is Experimental then also get oscilator strengths, otherwise calculate them
        if self.reader.type == "Experimental output":
            self.oscillator_strengths = self.reader.oscillator_strengths
        else:
            #
            # calculate normal modes in xyz coordinate space
            masses = np.array(self.reader.masses) * amu
            normal_modes = Calculator.normal_modes(masses, self.mass_weighted_normal_modes)
            # from the normal modes and the born charges calculate the oscillator strengths of each mode
            self.oscillator_strengths = Calculator.oscillator_strengths(normal_modes, born_charges)
            # Apply point-group symmetrisation when enabled (same setting as ε_∞)
            if self.settings.get("Symmetrise optical permittivity", True):
                self.oscillator_strengths = self._symmetrise_oscillator_strengths(self.oscillator_strengths)

        # calculate the intensities from the trace of the oscillator strengths
        self.intensities = Calculator.infrared_intensities(self.oscillator_strengths)
        self.modal_pair_use_nac = self._classify_modal_pair_nac_modes(self.oscillator_strengths)

        # calculate Raman activities if Raman tensors are available
        raman_tensors = self.reader.get_raman_tensors() if self.reader else None

        # Set the Raman intensities for each mode, set to zero if they are not available
        show_raman_col = raman_tensors is not None and len(raman_tensors) == len(self.frequencies_cm1)
        if show_raman_col:
            _acts = Calculator.raman_intensities(raman_tensors, self.reader.volume)
            self.raman_intensities      = _acts[:, 0].tolist()
            self.raman_intensities_par  = _acts[:, 1].tolist()
            self.raman_intensities_perp = _acts[:, 2].tolist()
        else:
            self.raman_intensities      = [0.0] * len(self.frequencies_cm1)
            self.raman_intensities_par  = [0.0] * len(self.frequencies_cm1)
            self.raman_intensities_perp = [0.0] * len(self.frequencies_cm1)

        # Decide which modes to select based on IR/Raman activity
        if self.recalculate_selected_modes and len(self.intensities) > 0 and len(self.frequencies_cm1) > 0:
            logger.debug("create_intensity_table: recalculating selected modes")
            is_raman = "Raman" in self.settings["Spectroscopy type"]
            self.modes_selected = []
            self.mode_list = []
            for f, ir_intensity, raman_activity in zip(self.frequencies_cm1, self.intensities, self.raman_intensities):
                if is_raman:
                    self.modes_selected.append(f > 10.0 and raman_activity > 1.0e-6)
                else:
                    self.modes_selected.append(f > 10.0 and ir_intensity > 1.0e-6)
            self.mode_list = [i for i,mode in enumerate(self.modes_selected) if mode]
            logger.debug(f"Selected modes are; {self.mode_list}")
            self.recalculate_selected_modes = False

        # Calculate the ionic contribution to the permittivity
        frequencies_au = wavenumber*np.array(self.frequencies_cm1)
        volume_au = self.reader.volume*angstrom*angstrom*angstrom
        self.epsilon_ionic = Calculator.ionic_permittivity(self.mode_list, self.oscillator_strengths, frequencies_au, volume_au )
        # Make an np array of epsilon infinity
        epsilon_inf = np.array(self.settings["Optical permittivity"])

        # If the reader already has a Dielectric Constant then use this
        # Otherwise create one from the data that has been read in
        if self.reader.CrystalPermittivity:
            self.CrystalPermittivityObject = self.reader.CrystalPermittivity
        else:
            drude = False
            drude_plasma_au = 0
            drude_sigma_au = 0
            sigmas_au = np.array(self.sigmas_cm1)*wavenumber
            logger.debug(f"CreateIntensityTable: Calculating dielectric {self.mode_list}")
            self.CrystalPermittivityObject = DielectricFunction.DFT(
                                         self.mode_list, frequencies_au, sigmas_au, self.oscillator_strengths,
                                         volume_au, drude, drude_plasma_au, drude_sigma_au )

        # Add the optical permittivity to the dielctric function
        self.CrystalPermittivityObject.set_epsilon_infinity(epsilon_inf)

        # Prepare to finish and output table
        ncols = 9 if show_raman_col else 6
        self.output_tw.setRowCount(len(self.sigmas_cm1))
        self.output_tw.setColumnCount(ncols)

        # Set the table headers
        if show_raman_col:
            raman_display = self._raman_activity_display_metadata()
            self.output_tw.setHorizontalHeaderLabels([
                "   Sigma   \n(cm-1)", " TO freq \n(cm-1)", " LO freq \n(cm-1)",
                "  Intensity  \n(Debye²/Å²/amu)",
                "Integrated Molar Absorption\n(L/mole/cm²)",
                "Absorption maximum\n(L/mole/cm)",
                f"Raman total\n({raman_display['unit_label']})",
                f"Raman ∥\n({raman_display['unit_label']})",
                f"Raman ⟂\n({raman_display['unit_label']})",
            ])
            for col in (6, 7, 8):
                header_item = self.output_tw.horizontalHeaderItem(col)
                if header_item is not None:
                    header_item.setToolTip(raman_display["tooltip"])
        else:
            self.output_tw.setHorizontalHeaderLabels([
                "   Sigma   \n(cm-1)", " TO freq \n(cm-1)", " LO freq \n(cm-1)",
                "  Intensity  \n(Debye²/Å²/amu)",
                "Integrated Molar Absorption\n(L/mole/cm²)",
                "Absorption maximum\n(L/mole/cm)",
            ])

        # Draw the table then fill in the LO frequency column
        self.redraw_output_tw()
        self._refresh_lo_columns()
        QCoreApplication.processEvents()
        # LO freq header renders narrower than TO freq (font width of 'L' < 'T');
        # force column 2 to be at least as wide as column 1.
        to_width = self.output_tw.columnWidth(1)
        if self.output_tw.columnWidth(2) < to_width:
            self.output_tw.setColumnWidth(2, to_width)
        logger.debug("Finished:: create_intensity_table")
        return

    def _raman_activity_display_metadata(self):
        """Return display metadata for Raman activities."""
        units = self.settings.get("Raman activity units", "polarizability")
        if units == "epsilon":
            return {
                "key": "epsilon",
                "factor": 1.0,
                "unit_label": "Å/amu",
                "tooltip": "Displayed in the internal R_epsilon convention: R_epsilon = sqrt(Vcell) dε/dQ.",
            }
        factor = 1.0
        if self.reader is not None and getattr(self.reader, "volume", None) is not None:
            factor = self.reader.volume / (16.0 * np.pi * np.pi)
        return {
            "key": "polarizability",
            "factor": factor,
            "unit_label": "Å⁴/amu",
            "tooltip": (
                "Displayed in polarizability-volume units: "
                "R_alpha = d[Vcell(epsilon-I)/(4*pi)]/dQ = sqrt(Vcell) R_epsilon/(4*pi)."
            ),
        }

    def _raman_activity_display_factor(self):
        """Return the factor converting R_epsilon activities to display units."""
        return self._raman_activity_display_metadata()["factor"]

    def request_refresh(self):
        """Request to refresh the environment or state.

        This method sets a flag to indicate that a refresh is required, logging the start and end of its execution.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        logger.debug("Start:: request_refresh")
        self.refresh_required = True
        logger.debug("Finished:: request_refresh")
        return

    def write_spreadsheet(self):
        """Write the configuration and calculation results to the associated spreadsheet.

        This function iterates through various settings and results stored within the object, such as optical
        permittivity, mass definitions, and vibrational mode calculations. For each setting or result, it writes the
        relevant information to a new row in the spreadsheet. The spreadsheet is expected to have a 'Settings'
        worksheet, which this function selects and updates with the provided values.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        logger.debug("Start:: write_spreadsheet")

        # Set the spreadsheet and test if it exists in the Notebook
        sp = self.notebook.spreadsheet
        if sp is None:
            logger.debug("Finished:: write_spreadsheet - Aborting write of spreadsheet")
            return
        logger.debug("Writing of spreadsheet")

        # Select and delete the existing Settings tab of the worksheet
        sp.select_work_sheet("Settings")
        sp.delete()

        # Start writing the rows with information
        sp.write_next_row(["Settings and calculations of frequencies and absorption"], row=0, col=1)
        for item in sorted(self.settings):
            if item == "Optical permittivity" and self.settings[item] is not None:
               sp.write_next_row([item], col=1)
               sp.write_next_row(self.settings[item][0], col=2, check=1)
               sp.write_next_row(self.settings[item][1], col=2, check=1)
               sp.write_next_row(self.settings[item][2], col=2, check=1)
               eps = self.epsilon_ionic.tolist()
               sp.write_next_row(["Ionic Contribution to permittivity at zero frequency"], col=1)
               sp.write_next_row(eps[0], col=2, check=1)
               sp.write_next_row(eps[1], col=2, check=1)
               sp.write_next_row(eps[2], col=2, check=1)
               eps = self.epsilon_ionic + np.array(self.settings[item])
               eps = eps.tolist()
               sp.write_next_row(["Total Contribution to permittivity at zero frequency"], col=1)
               sp.write_next_row(eps[0], col=2, check=1)
               sp.write_next_row(eps[1], col=2, check=1)
               sp.write_next_row(eps[2], col=2, check=1)
            elif item == "Mass definition":
               sp.write_next_row([item,self.settings[item]], col=1, check=1)
               for c in self.masses_dictionary:
                   sp.write_next_row([c,self.masses_dictionary[c]],col=3,check=1)
            else:
               sp.write_next_row([item,self.settings[item]], col=1, check=1)
        sp.write_next_row([""], col=1)

        # Now write out the a spreadsheet of information
        sp.write_next_row(["Mode","Include?","Sigma(cm-1)","Frequency(cm-1)","Intensity(Debye2/Angs2/amu)",
                           "Integrated Molar Absorption(L/mole/cm2)","Absorption maximum (L/mole/cm)",
                           "Raman R_total (Å4/amu)","Raman R_parallel (Å4/amu)","Raman R_perp (Å4/amu)"], col=1)
        for mode,(f,intensity,raman_act,raman_par,raman_perp,sigma,selected) in enumerate(zip(
                self.frequencies_cm1, self.intensities, self.raman_intensities,
                self.raman_intensities_par, self.raman_intensities_perp,
                self.sigmas_cm1, self.modes_selected)):
            yn = "No"
            if selected:
                yn = "Yes"
            sp.write_next_row([mode, yn, sigma, f, intensity, 4225.6*intensity, 2*4225.6*intensity/sigma/np.pi,
                               raman_act, raman_par, raman_perp], col=1)
        logger.debug("Finished:: write_spreadsheet")

    def redraw_output_tw(self):
        """Redraws the output table widget with updated values.

        This method loops through frequencies, sigmas, and intensities, updating the table widget rows
        with new items representing these values and custom flags for item behavior. It also calculates
        and displays additional derived values in the table. Selection flags for each item are set
        based on whether the mode is selected or not. The table's columns are resized to content at the end.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - This method assumes that `self.frequencies_cm1`, `self.sigmas_cm1`, `self.intensities`, and
          `self.modes_selected` contain iterable values that are indexed correspondingly. - The method temporarily
          blocks signals from `self.output_tw` to prevent unwanted signal handling during updates. - Table item flags
          are set to control their selectability, checkability, editability, and enablement based on whether the
          corresponding mode is selected. - After updating, it ensures the table columns fit the content and re-enables
          the signals. - Uses `Qt` enumeration for setting check state, item flags, and text alignment.

        """        
        logger.debug("Start:: redraw_output_tw")
        ncols = self.output_tw.columnCount()
        show_lo_col    = ncols >= 6   # col 2 is LO freq for both IR (6-col) and Raman (9-col)
        show_raman_col = ncols == 9   # Raman-only cols (total/∥/⟂) at cols 6–8
        raman_display_factor = self._raman_activity_display_factor()
        self.output_tw.blockSignals(True)

        # Loop over modes, frequencies, intensities and raman activities
        for i,(f,sigma,intensity,raman_act,raman_par,raman_perp) in enumerate(zip(
                self.frequencies_cm1, self.sigmas_cm1, self.intensities,
                self.raman_intensities, self.raman_intensities_par, self.raman_intensities_perp)):
            #
            # Sigma and check / unchecked column
            items = []
            itemFlags = []
            item = QTableWidgetItem(f"{sigma:.1f}")
            if self.modes_selected[i]:

                # Set selection flags for active modes
                item.setCheckState(Qt.Checked)
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable )
                freqFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
                otherFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled | Qt.ItemIsSelectable
            else:

                # Set selection flags for inactive modes
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable )
                item.setCheckState(Qt.Unchecked)
                freqFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled | Qt.ItemIsSelectable
                otherFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled | Qt.ItemIsSelectable
            items.append(item)

            # TO frequency column cm-1 (col 1)
            items.append(QTableWidgetItem(f"{f:.4f}" ) )
            itemFlags.append( freqFlags )

            # LO frequency placeholder (col 2; filled by _refresh_lo_columns for both IR and Raman)
            if show_lo_col:
                lo_item = QTableWidgetItem("—")
                lo_item.setFlags(otherFlags)
                items.append(lo_item)
                itemFlags.append(otherFlags)

            # Intensity column Debye2/Angs2/amu (col 3 when LO present, col 2 otherwise)
            items.append(QTableWidgetItem(f"{intensity:.4f}" ) )
            itemFlags.append( otherFlags )

            # Integrated molar absorption L/mole/cm/cm
            items.append(QTableWidgetItem(f"{intensity*4225.6:.2f}" ) )
            itemFlags.append( otherFlags )

            # Maximum extinction L/mole/cm
            items.append(QTableWidgetItem(f"{2*intensity*4225.6/self.sigmas_cm1[i]/np.pi:.2f}" ) )
            itemFlags.append( otherFlags )

            # Raman columns displayed in the selected activity units.
            if show_raman_col:
                items.append(QTableWidgetItem(f"{raman_act * raman_display_factor:.6f}"))
                itemFlags.append(otherFlags)
                items.append(QTableWidgetItem(f"{raman_par * raman_display_factor:.6f}"))
                itemFlags.append(otherFlags)
                items.append(QTableWidgetItem(f"{raman_perp * raman_display_factor:.6f}"))
                itemFlags.append(otherFlags)

            # Set the text alignment
            for j,(item,flag) in enumerate(zip(items,itemFlags)):
                item.setFlags(flag)
                item.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                self.output_tw.setItem(i, j, item )

        # Resize the column widths to content
        self.output_tw.resizeColumnsToContents()

        # Release the block on signals for the frequency output table
        self.output_tw.blockSignals(False)
        QCoreApplication.processEvents()
        logger.debug("Finished:: redraw_output_tw")

    def _refresh_lo_columns(self):
        """Recompute and display LO frequencies, LO IR intensities, and LO/EO Raman strengths.

        Called whenever the NAC checkbox, q̂ spinboxes, or EO checkbox change, and
        at the end of :meth:`create_intensity_table`.  Acts when the modes table has
        6 columns (IR) or 9 columns (Raman).

        Four cases are handled:

        * NAC off, EO off  → TO frequencies (col 2 = "—"), TO IR intensities, TO Raman strengths.
        * NAC off, EO on   → TO frequencies, TO IR intensities, EO-corrected Raman strengths.
        * NAC on,  EO off  → LO frequencies, LO IR intensities, LO Raman strengths.
        * NAC on,  EO on   → LO frequencies, LO IR intensities, EO-corrected LO Raman strengths.
        """
        ncols = self.output_tw.columnCount()
        if ncols < 6:
            return
        show_raman_col = (ncols == 9)
        n_rows = self.output_tw.rowCount()
        raman_display_factor = self._raman_activity_display_factor()

        nac_on = hasattr(self, "_nac_apply_cb") and self._nac_apply_cb.isChecked()
        eo_on  = (show_raman_col
                  and hasattr(self, "_lo_eo_cb") and self._lo_eo_cb.isChecked()
                  and self.reader is not None
                  and getattr(self.reader, "nonlinear_optical_susceptibility", None) is not None)

        # Outputs — None means fall back to TO values
        lo_freqs         = None   # fills col 2; None → "—"
        lo_ir_intensities = None  # fills cols 3–5; None → self.intensities
        lo_raman_acts    = None   # fills cols 6–8; None → self.raman_intensities*

        # Check whether NAC/EO ingredients are available
        has_nac_data = (
            self.reader is not None
            and self.mass_weighted_normal_modes is not None
            and getattr(self.reader, "hessian", None) is not None
            and len(getattr(self.reader, "born_charges", [])) > 0
        )

        # Normalise q̂ (needed by both NAC and EO)
        q_hat = None
        if (nac_on or eo_on) and has_nac_data:
            q_raw  = np.array([spin.value() for spin in self._lo_q_spins])
            q_norm = float(np.linalg.norm(q_raw))
            if q_norm >= 1e-12:
                q_hat = q_raw / q_norm

        if q_hat is not None and has_nac_data:
            from PDielec.RamanPolarCalculator import (
                apply_eo_correction,
                build_eigvecs_from_normal_modes,
                build_Z_mat,
                compute_lo_modes,
            )

            # Build shared NAC ingredients
            masses_au    = np.array(self.reader.masses) * amu
            born_charges = np.array(self.reader.born_charges)
            eps_inf      = np.array(self.settings["Optical permittivity"], dtype=float)
            if eps_inf.ndim == 1:
                eps_inf = np.diag(eps_inf)
            Z_mat   = build_Z_mat(born_charges, masses_au)
            eigvecs = build_eigvecs_from_normal_modes(self.mass_weighted_normal_modes)  # (3N, n_to)
            U_TO    = eigvecs.T                                                          # (n_to, 3N)

            if nac_on:
                # --- LO frequencies and eigenvectors ---
                lo_freqs, lo_eigvecs = compute_lo_modes(
                    q_hat, np.array(self.reader.hessian, dtype=float),
                    U_TO, Z_mat, eps_inf, self.reader.volume * angs2bohr ** 3,
                )

                # --- LO IR intensities ---
                n_atoms  = len(self.reader.masses)
                n_to     = lo_eigvecs.shape[1]
                lo_mwnm  = [lo_eigvecs[:, n].reshape(n_atoms, 3) for n in range(n_to)]
                lo_nm    = Calculator.normal_modes(masses_au, lo_mwnm)
                lo_osc   = Calculator.oscillator_strengths(lo_nm, born_charges)
                lo_ir_intensities = Calculator.infrared_intensities(lo_osc)

                # --- LO Raman strengths (optionally EO-corrected) ---
                if show_raman_col:
                    raman_tensors = self.reader.get_raman_tensors()
                    if raman_tensors and len(raman_tensors) >= n_rows:
                        # Project TO Raman tensors onto LO mode basis
                        n_r = min(n_to, len(raman_tensors))
                        t_stack = np.array(
                            [np.asarray(raman_tensors[k], dtype=complex) for k in range(n_r)]
                        )  # (n_r, 3, 3)
                        lo_raman_tensors = []
                        for n in range(n_to):
                            C_n  = U_TO[:n_r, :] @ lo_eigvecs[:, n]   # (n_r,) overlaps
                            R_lo = np.einsum("k,kij->ij", C_n, t_stack)
                            lo_raman_tensors.append(R_lo)
                        if eo_on:
                            chi2 = self.reader.nonlinear_optical_susceptibility
                            lo_raman_tensors = apply_eo_correction(
                                lo_raman_tensors, chi2, q_hat, Z_mat, lo_eigvecs, eps_inf)
                        lo_raman_acts = Calculator.raman_intensities(lo_raman_tensors, self.reader.volume)

            elif eo_on:
                # NAC off but EO on: apply EO correction to TO Raman tensors.
                raman_tensors = self.reader.get_raman_tensors()
                if raman_tensors and len(raman_tensors) == n_rows:
                    chi2      = self.reader.nonlinear_optical_susceptibility
                    corrected = apply_eo_correction(raman_tensors, chi2, q_hat, Z_mat, eigvecs, eps_inf)
                    lo_raman_acts = Calculator.raman_intensities(corrected, self.reader.volume)

        # --- Update the table ---
        self.output_tw.blockSignals(True)
        for i in range(n_rows):
            # col 2: LO frequency
            lo_item = self.output_tw.item(i, 2)
            if lo_item is not None:
                lo_item.setText(
                    f"{lo_freqs[i]:.4f}" if lo_freqs is not None and i < len(lo_freqs) else "—"
                )

            # cols 3, 4, 5: IR intensities (LO or TO)
            ir_val  = (lo_ir_intensities[i]
                       if lo_ir_intensities is not None and i < len(lo_ir_intensities)
                       else (self.intensities[i] if i < len(self.intensities) else 0.0))
            sigma_i = self.sigmas_cm1[i] if i < len(self.sigmas_cm1) else 1.0
            for col, text in (
                (3, f"{ir_val:.4f}"),
                (4, f"{ir_val * 4225.6:.2f}"),
                (5, f"{2.0 * ir_val * 4225.6 / sigma_i / np.pi:.2f}"),
            ):
                item = self.output_tw.item(i, col)
                if item is not None:
                    item.setText(text)

            # cols 6, 7, 8: Raman strengths (LO/EO-corrected or TO)
            if show_raman_col:
                if lo_raman_acts is not None and i < len(lo_raman_acts):
                    for col_offset in range(3):
                        item = self.output_tw.item(i, 6 + col_offset)
                        if item is not None:
                            item.setText(f"{lo_raman_acts[i, col_offset] * raman_display_factor:.6f}")
                else:
                    to_vals = (
                        self.raman_intensities[i]      if i < len(self.raman_intensities)      else 0.0,
                        self.raman_intensities_par[i]  if i < len(self.raman_intensities_par)  else 0.0,
                        self.raman_intensities_perp[i] if i < len(self.raman_intensities_perp) else 0.0,
                    )
                    for col_offset, val in enumerate(to_vals):
                        item = self.output_tw.item(i, 6 + col_offset)
                        if item is not None:
                            item.setText(f"{val * raman_display_factor:.6f}")
        self.output_tw.blockSignals(False)

    def on_sigma_changed(self):
        """Update the sigma value in settings and apply it across the frequency range, then refreshes the output.

        This method is triggered when the sigma value has been changed. It updates the stored sigma value in the
        settings, applies the new sigma value across the pre-defined frequency range, redraws the output table widget,
        and then initiates a refresh to reflect the newly applied changes.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        logger.debug("Start:: redraw_output_tw")
        self.settings["Sigma value"] = self.sigma_sb.value()
        self.sigmas_cm1 = [ self.settings["Sigma value"] for i in self.frequencies_cm1 ]
        self.redraw_output_tw()
        logger.debug(f"on sigma change {self.settings['Sigma value']}")
        self.refresh_required = True
        self.refresh()
        logger.debug("Finished:: on_sigma_changed")

    def on_mass_cb_activated(self,index):
        """Handle activation of a combobox option in a mass-related setting.

        This method updates the current mass definition based on the selected option
        from a combobox, updates related settings, and refreshes the GUI accordingly.
        It disables the 4th option in the combobox if one of the first three options is selected.
        Additionally, it processes all pending events to ensure the GUI is updated immediately.

        Parameters
        ----------
        index : int
            The index of the selected option in the combobox.

        Returns
        -------
        None

        Notes
        -----
        - `refresh()` is called to update the UI based on the new selection.
        - It also makes use of `QCoreApplication.processEvents()` to ensure the UI is 
          responsive by processing all pending events.

        See Also
        --------
        set_masses_tw : A method to update the masses table widget based on the selection.
        refresh : A method to refresh the UI components.

        """        
        logger.debug(f"Start:: on_mass_combobox_activated {self.mass_cb.currentText()}")
        self.settings["Mass definition"] = self.mass_definition_options[index]
        self.current_mass_definition_index = index
        if index < 3:
            self.mass_cb.model().item(3).setEnabled(False)

        # Modify the element masses
        self.set_masses_tw()
        self.refresh_required = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        logger.debug(f"Finished:: on_mass_combobox_activated {self.mass_cb.currentText()}")

    def set_masses_tw(self):
        """Set the element masses in the table widget based on the mass_definition setting.

        This method updates the element masses table widget with masses according to the selected mass definition mode
        ('average', 'program', 'isotope', 'gui'). It retrieves species information from the reader attribute, calculates
        or fetches the corresponding masses, and updates an internal dictionary and the table widget accordingly.
        Signals from the table widget are temporarily blocked to prevent unwanted event triggers during updates.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        logger.debug("Start:: set_masses_tw")
        if self.reader:
            self.element_masses_tw.blockSignals(True)

            # Get the species list
            species = self.reader.get_species()

            # set the initial dictionary according to the mass_definition
            masses = []
            if self.settings["Mass definition"] == "average":

                # Set up the average masses 
                self.mass_cb.setCurrentIndex(0)
                for element in species:
                    mass = average_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings["Mass definition"] == "program":

                # Set up the masses as used by the DFT program
                self.mass_cb.setCurrentIndex(1)
                self.reader.reset_masses()
                masses = self.reader.masses_per_type
                for mass,element in zip(masses,species):
                    self.masses_dictionary[element] = mass
            elif self.settings["Mass definition"] == "isotope":

                # Set up the isotopic masess
                self.mass_cb.setCurrentIndex(2)
                for element in species:
                    mass = isotope_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings["Mass definition"] == "gui":

                # Set up the masses defined in the GUI
                for element in species:
                    mass = self.masses_dictionary[element]
                    masses.append(mass)
            else:
                 logger.debug(f"Error mass_definition not recognised {self.settings['Mass definition']}")
            self.element_masses_tw.setColumnCount(len(masses))
            self.element_masses_tw.setHorizontalHeaderLabels(species)
            self.element_masses_tw.setVerticalHeaderLabels([""])
            logger.debug(f"masses_dictionary {self.masses_dictionary}")
            logger.debug(f"masses {masses}")
            #
            # set masses of the elements in the table widget according to the mass definition
            for i,(mass,element) in enumerate(zip(masses,species)):
                logger.debug(f"set_masses_tw {self.settings['Mass definition']} {i} {mass} {element}")
                qw = QTableWidgetItem()
                if self.settings["Mass definition"] == "program":

                    # Set up the program masses 
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{mass:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    self.element_masses_tw.setItem(0,i, qw )
                elif self.settings["Mass definition"] == "average":

                    # Set up the average masses 
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{average_masses[element]:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    logger.debug(f"average {average_masses[element]}")
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings["Mass definition"] == "isotope":

                    # Set up the isotopic masses 
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{isotope_masses[element]:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    logger.debug(f"isotope {isotope_masses[element]}")
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings["Mass definition"] == "gui":

                    # Set up the gui masses 
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{self.masses_dictionary[element]:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    logger.debug(f"gui {self.masses_dictionary[element]}")
                    self.element_masses_tw.setItem(0,i, qw )
                else:
                    logger.warning(f"Mass definition not processed {self.settings['Mass definition']}")

            # unblock the table signals
            self.element_masses_tw.blockSignals(False)

        # Wait for any un processed events
        QCoreApplication.processEvents()
        logger.debug("Finished:: set_masses_tw")

    def on_output_tw_itemChanged(self, item):
        """Handle item changes in the output table widget.

        This method is called when an item in the output table widget (`output_tw`) is changed. It manages the selection
        state of modes based on the state of checkbox items in the first column, updates the "sigmas_cm1" and
        "frequencies_cm1" lists based on the user input in the respective columns, refreshes the widget, and processes
        pending GUI events.

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was changed.

        Returns
        -------
        None

        Notes
        -----
        - This method blocks signals from the output table widget at the beginning to prevent recursive calls during its
          execution. - It distinguishes between changes in the first column (selection state and mode sigma values) and
          the second column (frequency values). - The method updates internal state to reflect the changes made by the
          user, including which modes are selected and the values for "sigmas_cm1" and "frequencies_cm1". - After making
          changes, this method refreshes the output table widget and enforces a GUI update with
          `QCoreApplication.processEvents()`.

        See Also
        --------
        redraw_output_tw : Redraws the output table widget.
        refresh : Refreshes the internal state based on the latest changes.

        """        
        self.output_tw.blockSignals(True)
        logger.debug(f"Start:: on_output_tw_itemChanged {item.row()} {item.column()}")
        col = item.column()
        row = item.row()
        if col == 0:

            # If this is the first column alter the check status but reset the sigma value
            if item.checkState() == Qt.Checked:

                # Item was checked set selected to true and create and new mode_list
                logger.debug(f"on_output_tw_itemChanged setting selected mode to True {row}")
                self.modes_selected[row] = True
                self.mode_list = [i for i,mode in enumerate(self.modes_selected) if mode]
            else:

                # Item was not checked set selected to false and create and new mode_list
                logger.debug(f"on_output_tw_itemChanged setting selected mode to False {row}")
                self.modes_selected[row] = False
                self.mode_list = [i for i,mode in enumerate(self.modes_selected) if mode]

            # Check to see if a new value has been given for sigma
            new_value = float(item.text())
            if new_value != self.sigmas_cm1[row]:
                self.sigmas_cm1[row] = new_value
                self.redraw_output_tw()
        elif col == 1:
            # This is the frequency column
            if item.text() == "":
                self.frequencies_have_been_edited = False
            else:
                new_value = float(item.text())
                self.frequencies_cm1[row] = new_value
                self.frequencies_have_been_edited = True
            self.redraw_output_tw()
        else:
            self.redraw_output_tw()
        logger.debug(f"on_output_tw_itemChanged selected_modes {self.modes_selected}")
        self.refresh_required = True
        self.refresh()
        QCoreApplication.processEvents()
        logger.debug("Finished:: on_output_tw_itemChanged")

    def on_element_masses_tw_itemClicked(self, item):
        """Handle the item clicked event on the element mass table widget.

        The only thing this function does is to make sure that signal from the `element_masses_tw` table widget are not
        blocked,

        Parameters
        ----------
        item : QTableWidgetItem
            The item in the table widget that was clicked.

        Returns
        -------
        None

        """        
        logger.debug(f"Start:: on_element_masses_tw_itemClicked) {item.row()} {item.column()}")
        self.element_masses_tw.blockSignals(False)
        logger.debug("Finished:: on_element_masses_tw_itemClicked)" )

    def on_element_masses_tw_itemChanged(self, item):
        """Handle item changed event in the mass table widget.

        This function is triggered whenever a cell in the masses table widget is edited. It updates the internal masses
        dictionary based on the new value provided in the table cell. Additionally, it sets some UI components to
        reflect the change and flags the need for recalculation and refresh.

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was changed.

        Returns
        -------
        None

        Notes
        -----
        - `self.reader.get_species()` returns a list of elements corresponding to table columns.
        - `self.masses_dictionary` is a dict where keys are element symbols and values are their masses.
        - `self.refresh_required` and `self.recalculate_selected_modes` are boolean flags used to indicate that a
          refresh and a recalculation are necessary. - `self.refresh()` is a method that refreshes or updates the UI
          components as needed based on the new changes.

        """        
        logger.debug(f"Start:: on_element_masses_tw_itemChanged) {item.row()} {item.column()}")
        elements = self.reader.get_species()
        col = item.column()
        self.mass_cb.model().item(3).setEnabled(True)
        self.settings["Mass definition"] = "gui"
        self.mass_cb.setCurrentIndex(3)
        self.masses_dictionary[elements[col]] = float(item.text())
        logger.debug(f"masses_dictionary {self.masses_dictionary}")
        self.refresh_required = True
        self.recalculate_selected_modes = True
        self.refresh()
        logger.debug("Finished:: on_element_masses_tw_itemChanged)" )

    def on_optical_tw_itemChanged(self, item):
        """Handle item change events for an optical permittivity table.

        This method updates the optical permittivity settings based on changes made to the table widget items. It
        enforces symmetry in the optical permittivity matrix, marks the settings as edited, and triggers a series of
        updates and recalculations.

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was changed. The item's row, column, and text are used to update the settings.

        Returns
        -------
        None

        Notes
        -----
        - Updates the 'Optical permittivity' matrix in the settings dictionary based on the row, column, and text of the
          changed item. - Enforces symmetry by setting the [column][row] element equal to the [row][column] element. -
          Marks the 'Optical permittivity edited' flag as True. - Calls `refresh_optical_permittivity_tw` to refresh the
          optical permittivity table widget. - Flags `refresh_required`, `recalculate_selected_modes`, and `refresh`
          methods/attributes for further actions. - Processes any pending events with `QCoreApplication.processEvents()`
          to ensure the UI remains responsive.

        """        
        logger.debug(f"Start::on_optical_itemChanged) {item.row()} {item.column()} {item.text()}")
        self.settings["Optical permittivity"][item.row()][item.column()] = float(item.text())
        self.settings["Optical permittivity"][item.column()][item.row()] = float(item.text())
        self.settings["Optical permittivity edited"] = True
        self.refresh_optical_permittivity_tw()
        self.refresh_required = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        logger.debug("Finished::on_optical_itemChanged)")
        return

    def on_optical_tw_itemClicked(self, item):
        """Handle item clicked event for an optical widget's table item.

        This method is designed to be connected to the itemClicked signal of a table widget
        in a graphical user interface built with PyQt or PySide. When a table widget item
        is clicked, it first prints the start message along with the row and column of the
        clicked item. It then unblocks signals for the table widget and processes any pending
        GUI events, ensuring the application remains responsive. 

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was clicked.

        Returns
        -------
        None

        """        
        logger.debug(f"Start:: on_optical_itemClicked) {item.row()} {item.column()}")
        self.optical_tw.blockSignals(False)
        QCoreApplication.processEvents()
        logger.debug("Finished:: on_optical_itemClicked)")
        return

    def on_spectroscopy_type_cb_activated(self, index):
        """Handle changes to the spectroscopy type combo box.

        Updates the spectroscopy type setting and tells the notebook to replace
        all scenarios with a new one of the appropriate type.

        Parameters
        ----------
        index : int
            The index of the selected spectroscopy type.

        Returns
        -------
        None

        """
        logger.debug(f"Start:: on_spectroscopy_type_cb_activated {index}")
        spectroscopy_type = self.spectroscopy_types[index]
        self.settings["Spectroscopy type"] = spectroscopy_type
        self.recalculate_selected_modes = True
        self.refresh_required = True
        self.update_raman_polar_button()
        self.notebook.set_spectroscopy_type(spectroscopy_type)
        logger.debug(f"Finished:: on_spectroscopy_type_cb_activated {index}")

    def on_raman_activity_units_cb_activated(self, index):
        """Handle changes to the Raman activity display units combo box."""
        key = self.raman_activity_unit_options[index][0]
        if self.settings.get("Raman activity units") == key:
            return
        self.settings["Raman activity units"] = key
        self.create_intensity_table()
        if hasattr(self, "_raman_polar_windows"):
            display = self._raman_activity_display_metadata()
            for window in list(self._raman_polar_windows):
                if window is not None:
                    window.set_activity_display(display)
        QCoreApplication.processEvents()

    def update_raman_polar_button(self):
        """Show the Raman polar plot button for Raman spectroscopy; NAC options are always shown."""
        if not hasattr(self, "raman_polar_button"):
            return
        is_raman = self.settings.get("Spectroscopy type") in ("Powder Raman", "Crystal Raman")
        self.raman_polar_label.setVisible(is_raman)
        self.raman_polar_button.setVisible(is_raman)
        self.raman_polar_button.setEnabled(is_raman)
        if hasattr(self, "raman_activity_units_cb"):
            self.raman_activity_units_label.setVisible(is_raman)
            self.raman_activity_units_cb.setVisible(is_raman)
        # EO checkbox — requires Raman AND χ^(2) data
        if hasattr(self, "_lo_eo_cb"):
            has_chi2 = (is_raman
                        and self.reader is not None
                        and getattr(self.reader, "nonlinear_optical_susceptibility", None) is not None)
            self._lo_eo_cb.setVisible(has_chi2)

    def _collect_raman_polar_data(self):
        """Collect Raman tensor and SettingsTab mode data for the polar plot window."""
        if self.reader is None:
            return None, "No file has been read."
        raman_tensors = self.reader.get_raman_tensors()
        if raman_tensors is None or len(raman_tensors) == 0:
            return None, "No Raman tensors are available for the current file."
        if len(self.frequencies_cm1) == 0 or len(self.raman_intensities) == 0:
            self.create_intensity_table()
        active_indices = raman_active_mode_indices(
            self.frequencies_cm1,
            self.raman_intensities,
            raman_tensors,
        )
        if not active_indices:
            return None, "No Raman-active modes were found."

        modes = []
        for index in active_indices:
            selected = index < len(self.modes_selected) and self.modes_selected[index]
            modes.append({
                "index": index,
                "display_mode": index + 1,
                "frequency_cm1": float(self.frequencies_cm1[index]),
                "raman_total": float(self.raman_intensities[index]),
                "selected": selected,
            })
        if not any(mode["selected"] for mode in modes):
            for mode in modes:
                mode["selected"] = True

        data = {
            "raman_tensors": raman_tensors,
            "modes": modes,
            "raman_activity_display": self._raman_activity_display_metadata(),
        }

        # Pass current q̂ so the polar window initialises aligned with SettingsTab.
        if hasattr(self, "_lo_q_spins"):
            data["q_hat_init"] = tuple(spin.value() for spin in self._lo_q_spins)

        # NAC infrastructure — needed for LO frequencies (independent of EO).
        has_nac = (
            getattr(self.reader, "hessian", None) is not None
            and len(getattr(self.reader, "born_charges", [])) > 0
            and self.mass_weighted_normal_modes is not None
        )
        if has_nac:
            from PDielec.RamanPolarCalculator import build_eigvecs_from_normal_modes, build_Z_mat
            masses_au    = np.array(self.reader.masses) * amu
            born_charges = np.array(self.reader.born_charges)
            eps_inf      = np.array(self.settings["Optical permittivity"], dtype=float)
            if eps_inf.ndim == 1:
                eps_inf = np.diag(eps_inf)
            eigvecs = build_eigvecs_from_normal_modes(self.mass_weighted_normal_modes)
            data["Z_mat"]           = build_Z_mat(born_charges, masses_au)
            data["eigvecs"]         = eigvecs
            data["hessian"]         = np.array(self.reader.hessian, dtype=float)
            data["U_TO"]            = eigvecs.T
            data["volume_au"]       = self.reader.volume * angs2bohr ** 3
            data["volume_angstrom"] = self.reader.volume
            data["eps_inf"]         = eps_inf

        # EO correction — additionally requires χ^(2).
        chi2 = getattr(self.reader, "nonlinear_optical_susceptibility", None)
        if chi2 is not None and has_nac:
            data["chi2_repsilon"] = chi2

        return data, None

    def on_raman_polar_button_clicked(self):
        """Open the tensor-level Raman polar plot simulator."""
        data, message = self._collect_raman_polar_data()
        if data is None:
            QMessageBox.warning(self, "Raman polar plot", message)
            return
        from PDielec.GUI.RamanPolarWindow import RamanPolarWindow

        window = RamanPolarWindow(
            data,
            title="Raman polar mode simulator",
            parent=None,
        )
        window.show()
        if not hasattr(self, "_raman_polar_windows"):
            self._raman_polar_windows = []
        self._raman_polar_windows.append(window)

    def refresh(self, force=False):
        """Refresh the current state based on notebook content changes or user request.

        This method updates the GUI components such as tables for masses, optical permittivity, sigma value, Eckart term
        checkbox, and Born charges checkbox based on the changes made in the notebook or when forced by the user. It
        also triggers refreshes for any associated scenarios, plotting, analysis, and fitting tabs. The method blocks
        all signals from widgets during the update process to prevent unwanted calls to event handlers.

        Parameters
        ----------
        force : bool, optional
            If True, forces the refresh even if it's deemed not necessary. The default is False.

        Returns
        -------
        None

        Notes
        -----
        - Blocks all QWidget signals before starting the update process to prevent unintended event handling.
        - Updates various GUI components based on the current notebook's reader and settings information.
        - Triggers refreshes on associated tabs and scenarios within the notebook.
        - Unblocks all QWidget signals after updates are applied.
        - Manages internal flags to keep track of the need for refreshes and calculations.
        - Processes all pending Qt events with `QCoreApplication.processEvents()`.

        """        
        logger.debug(f"Start:: refresh {force}")
        if not self.reader and self.notebook.reader:
            self.refresh_required = True
        if not self.refresh_required and not force:
            logger.debug(f"Finished:: refresh not required {force}")
            return
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        self.reader = self.notebook.reader
        if self.reader:
            # Masses
            self.set_masses_tw()
            # Optical dielectric
            if self.settings["Optical permittivity"] is None:
                self.settings["Symmetrise optical permittivity"] = getattr(
                    self.reader,
                    "default_symmetrise_optical_permittivity",
                    True,
                )
                self.set_optical_permittivity_tw()
            else:
                self.refresh_optical_permittivity_tw()
        self.sigma_sb.setValue(self.settings["Sigma value"])
        index = self.spectroscopy_type_cb.findText(self.settings["Spectroscopy type"], Qt.MatchFixedString)
        if index >= 0:
            self.spectroscopy_type_cb.setCurrentIndex(index)
        if hasattr(self, "raman_activity_units_cb"):
            unit_keys = [key for key, _label in self.raman_activity_unit_options]
            unit_key = self.settings.get("Raman activity units", "polarizability")
            unit_index = unit_keys.index(unit_key) if unit_key in unit_keys else 0
            self.raman_activity_units_cb.setCurrentIndex(unit_index)
        self.update_raman_polar_button()
        if self.settings["Eckart flag"]:
            self.eckart_cb.setCheckState(Qt.Checked)
        else:
            self.eckart_cb.setCheckState(Qt.Unchecked)
        if self.settings["Neutral Born charges"]:
            self.born_cb.setCheckState(Qt.Checked)
        else:
            self.born_cb.setCheckState(Qt.Unchecked)
        if self.settings["Symmetrise optical permittivity"]:
            self.symmetrise_optical_cb.setCheckState(Qt.Checked)
        else:
            self.symmetrise_optical_cb.setCheckState(Qt.Unchecked)
        #
        # Unlock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.create_intensity_table()
        self.refresh_required = False
        # Refresh the scenarios that depend on this dielectric
        if self.notebook.scenarios is not None:
            for scenario in self.notebook.scenarios:
                scenario.request_refresh()
        if self.notebook.plottingTab is not None:
            self.notebook.plottingTab.request_refresh()
        if self.notebook.analysisTab is not None:
            self.notebook.analysisTab.request_refresh()
        if self.notebook.fitterTab is not None:
            self.notebook.fitterTab.request_refresh()
        QCoreApplication.processEvents()
        self.refresh_required = False
        self.calculation_required = True
        logger.debug(f"Finished:: refresh {force}")
        return

    def refresh_optical_permittivity_tw(self):
        """Refresh the table widget with optical permittivity values.

        This function takes no parameters explicitly, it operates on the instance's attributes.
        It updates the `optical_tw` table widget with optical permittivity values from the `settings` attribute.
        Formatting of the values is to four decimal places and aligned both horizontally and vertically to the center.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - Assumes: 
            - `self.settings` is a dict containing the key 'Optical permittivity' which maps to a 2D list (i.e., list of
              lists) of numeric values. - `self.optical_tw` is a Qt TableWidget or similar object with a `.setItem()`
              method accepting row, column, and QTableWidgetItem arguments. - The table widget's signals are temporarily
              blocked to prevent unintended side-effects during the update. - `Qt.AlignHCenter | Qt.AlignVCenter` is
              used to center-align the text, assuming `Qt` and `QTableWidgetItem` are properly imported and used. - The
              `QCoreApplication.processEvents()` call forces the application to process all pending events, ensuring the
              UI updates in real-time.

        """        
        logger.debug("Start:: refresh_optical_permittivity")
        optical = self.settings["Optical permittivity"]
        self.optical_tw.blockSignals(True)
        for i,row in enumerate(optical):
            for j, value in enumerate(row):
                qw = QTableWidgetItem(f"{value:.4f}" )
                qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                self.optical_tw.setItem(i,j,qw)
        self.optical_tw.blockSignals(False)
        QCoreApplication.processEvents()
        logger.debug("Finished:: refresh_optical_permittivity")
        return

    def _get_cartesian_rotations(self):
        """Return the list of Cartesian rotation matrices for the current unit cell.

        Fetches symmetry operations from spglib for the final unit cell in the reader,
        converts the integer fractional-coordinate rotation matrices to Cartesian
        rotation matrices, and returns them.  Returns an empty list if the reader has
        no unit-cell data or if spglib is unavailable.

        Returns
        -------
        list of ndarray, each shape (3, 3)
            Cartesian rotation matrices, or ``[]`` on failure.

        """
        if self.reader is None or not getattr(self.reader, "unit_cells", []):
            logger.debug("_get_cartesian_rotations: no unit cell")
            return []
        cell = self.reader.unit_cells[-1]
        if len(cell.fractional_coordinates) == 0 or len(cell.element_names) == 0:
            logger.debug("_get_cartesian_rotations: empty unit cell")
            return []
        try:
            import spglib
            numbers = cell.get_atomic_numbers()
            dataset = spglib.get_symmetry(
                (cell.lattice, cell.fractional_coordinates, numbers),
                symprec=1.0e-5,
            )
        except Exception as exc:
            logger.debug(f"_get_cartesian_rotations: spglib failed ({exc})")
            return []
        rotations = getattr(dataset, "rotations", None)
        if rotations is None and hasattr(dataset, "get"):
            rotations = dataset.get("rotations")
        if rotations is None or len(rotations) == 0:
            logger.debug("_get_cartesian_rotations: no rotations found")
            return []
        lattice_t = np.asarray(cell.lattice, dtype=float).T
        lattice_t_inv = np.linalg.inv(lattice_t)
        return [lattice_t @ np.asarray(R, dtype=float) @ lattice_t_inv for R in rotations]

    def _project_rank2_tensor(self, tensor, cart_rotations):
        """Project a real rank-2 tensor onto the invariant subspace of a point group.

        Computes ``T_sym[i,j] = (1/N) Σ_R  R[i,a] R[j,b] T[a,b]`` over the
        supplied Cartesian rotation matrices.  If the list is empty the raw tensor
        is returned unchanged.

        Parameters
        ----------
        tensor : array-like, shape (3, 3)
            Input tensor to symmetrise.
        cart_rotations : list of ndarray, each shape (3, 3)
            Cartesian rotation matrices from :meth:`_get_cartesian_rotations`.

        Returns
        -------
        ndarray, shape (3, 3)
            Symmetrised tensor.

        """
        t = np.array(tensor, dtype=float)
        if not cart_rotations:
            return t
        t_sym = np.zeros((3, 3), dtype=float)
        for R in cart_rotations:
            t_sym += np.einsum("ia,jb,ab->ij", R, R, t)
        t_sym /= float(len(cart_rotations))
        return t_sym

    def _symmetrise_optical_permittivity(self, tensor):
        """Return the optical permittivity averaged over the crystal point group.

        Each symmetry operation R (a Cartesian rotation matrix derived from the
        space-group rotations of the final unit cell) is applied as a rank-2 tensor
        transformation and the results are averaged:

            ε_sym[i,j] = (1/N) Σ_R  R[i,a] R[j,b] ε[a,b]

        This projects the raw DFT tensor onto the invariant subspace of the crystal
        point group, removing numerical noise that violates the lattice symmetry
        (e.g. spurious in-plane birefringence in a uniaxial crystal).

        Falls back to the raw tensor if the reader has no unit-cell information or
        if spglib is unavailable.

        Parameters
        ----------
        tensor : array-like, shape (3, 3)
            Raw optical permittivity tensor from the DFT reader.

        Returns
        -------
        list of list of float
            Point-group-symmetrised 3×3 tensor as a nested list (same format as
            stored in ``settings["Optical permittivity"]``).

        """
        cart_rotations = self._get_cartesian_rotations()
        t_sym = self._project_rank2_tensor(tensor, cart_rotations)
        t = np.array(tensor, dtype=float)
        if cart_rotations:
            diff = t - t_sym
            max_abs = float(np.max(np.abs(diff)))
            max_t = float(np.max(np.abs(t)))
            logger.info(
                f"Optical permittivity point-group symmetrisation: "
                f"max |ε - sym(ε)| = {max_abs:.3g}, "
                f"relative = {max_abs / max(max_t, 1e-30):.3g}, "
                f"operations = {len(cart_rotations)}"
            )
        return t_sym.tolist()

    def _symmetrise_oscillator_strengths(self, strengths):
        """Apply point-group symmetrisation to all oscillator strength tensors.

        Each oscillator strength ``S_mode`` is a real 3×3 tensor.  The same
        rank-2 point-group projection used for ε_∞ is applied to every mode:

            S_sym[mode, i, j] = (1/N) Σ_R  R[i,a] R[j,b] S[mode, a, b]

        For degenerate mode pairs (e.g. E₁ partners), individual mode tensors
        are modified but their sum is correctly symmetrised:
        ``S₇_sym + S₈_sym = (S₇ + S₈)_sym``.  This ensures that the full
        ionic contribution ``Σ S_mode/(ω_mode² − ω²)`` transforms correctly
        under the crystal point group at every frequency, including the laser
        wavelength where the Berreman propagation matrix is evaluated.

        Falls back to the raw list if no unit-cell data or spglib is unavailable.

        Parameters
        ----------
        strengths : list of array-like, each shape (3, 3)
            Oscillator strength tensors, one per normal mode, in the same units
            as returned by ``Calculator.oscillator_strengths``.

        Returns
        -------
        list of ndarray, each shape (3, 3)
            Point-group-symmetrised oscillator strength tensors.

        """
        cart_rotations = self._get_cartesian_rotations()
        if not cart_rotations:
            return strengths
        sym_strengths = [self._project_rank2_tensor(s, cart_rotations) for s in strengths]
        # Log the largest change
        max_change = max(
            float(np.max(np.abs(np.array(s_raw) - s_sym)))
            for s_raw, s_sym in zip(strengths, sym_strengths)
        )
        logger.info(
            f"Oscillator strength point-group symmetrisation: "
            f"max |S - sym(S)| = {max_change:.3g}, operations = {len(cart_rotations)}"
        )
        return sym_strengths

    def _classify_modal_pair_nac_modes(self, strengths):
        """Classify modes that should retain q-resolved modal-pair NAC treatment.

        Modes with negligible oscillator-strength tensor are treated as
        q-independent/non-polar and use the standard coherent total-field path
        in modal-pair calculations.  Modes with finite oscillator strength are
        polar/NAC-active and keep q-resolved modal-pair treatment.  This
        intrinsic classification is stable with respect to layer geometry,
        incidence angle, and polarisation.

        Parameters
        ----------
        strengths : list of array-like, each shape (3, 3)
            Oscillator strength tensors for all phonon modes.

        Returns
        -------
        list of bool
            ``True`` for polar/NAC-active modes; ``False`` for non-polar modes.

        """
        if strengths is None:
            return []
        norms = [float(np.linalg.norm(np.asarray(s, dtype=float))) for s in strengths]
        if len(norms) == 0:
            return []
        max_norm = max(norms)
        threshold = max(1.0e-12, 1.0e-10 * max_norm)
        return [norm > threshold for norm in norms]

    def on_symmetrise_optical_changed(self):
        """Handle a change in the 'Symmetrise optical permittivity' checkbox.

        Updates the setting, re-applies (or removes) symmetrisation to the currently
        loaded tensor, and triggers a full refresh so that all scenario tabs pick up
        the change.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        logger.debug(f"Start:: on_symmetrise_optical_changed {self.symmetrise_optical_cb.isChecked()}")
        self.settings["Symmetrise optical permittivity"] = self.symmetrise_optical_cb.isChecked()
        # Re-load from the reader so the (de)symmetrisation takes effect immediately,
        # but only if the user has not manually edited the tensor.
        if not self.settings["Optical permittivity edited"] and self.reader is not None:
            self.set_optical_permittivity_tw()
        self.refresh_required = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        logger.debug(f"Finished:: on_symmetrise_optical_changed {self.symmetrise_optical_cb.isChecked()}")
        return

    def set_optical_permittivity_tw(self):
        """Set the optical permittivity in the current settings.

        This method sets the 'Optical permittivity' value in the settings to the value of `zerof_optical_dielectric`
        from the reader (optionally Hermitian-symmetrised), triggers a refresh for the optical permittivity widget,
        and marks the selected modes for recalculation. It also signals that a refresh is required in the application
        and ensures the GUI events are processed to keep the application responsive.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        logger.debug("Start:: set_optical_permittivity_tw")
        raw = self.reader.zerof_optical_dielectric
        if self.settings.get("Symmetrise optical permittivity", True):
            self.settings["Optical permittivity"] = self._symmetrise_optical_permittivity(raw)
        else:
            self.settings["Optical permittivity"] = raw
        self.refresh_optical_permittivity_tw()
        self.recalculate_selected_modes = True
        self.refresh_required = True
        QCoreApplication.processEvents()
        logger.debug("Finished:: set_optical_permittivity_tw")
        return

    def on_born_changed(self):
        """Handle changes in the Born checkbox state.

        This method is triggered when the state of the Born checkbox changes. It updates the settings to reflect the new
        state, marks the session for refresh and recalculation, and then refreshes the UI.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        logger.debug(f"Start:: on_born_change {self.born_cb.isChecked()}")
        self.settings["Neutral Born charges"] = self.born_cb.isChecked()
        logger.debug(f"on born change {self.settings['Neutral Born charges']}")
        self.refresh_required = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        logger.debug(f"Finished:: on_born_change {self.born_cb.isChecked()}")
        return

    def on_eckart_changed(self):
        """Handle changes in the Eckart flag state.

        This function is triggered when the state of the Eckart checkbox changes.
        It logs the state change, updates the settings to reflect the new state,
        marks refresh and recalculations as needed, and triggers a UI refresh.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        logger.debug(f"Start:: on_eckart_change {self.eckart_cb.isChecked()}")
        self.settings["Eckart flag"] = self.eckart_cb.isChecked()
        logger.debug(f"on eckart change {self.settings['Eckart flag']}")
        self.refresh_required = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        logger.debug(f"Finished:: on_eckart_change {self.eckart_cb.isChecked()}")
        return

    def calculate(self,vs_cm1):
        """Calculate the permittivity of the crystal over the range of frequencies in vs_cm1.

        Parameters
        ----------
        vs_cm1 : list
            Frequencies in cm-1

        Returns
        -------
        None

        """
        logger.debug("Start:: calculate ")
        if len(vs_cm1) == 0:
            logger.debug("Finished:: calculate aborted vs_cm1 has not been set")
            return
        self.vs_cm1 = vs_cm1.copy()
        dielectricFunction = self.CrystalPermittivityObject.function()
        self.crystal_permittivity = []
        logger.debug("About to calculate settings crystal dielectric using pool")
        # Loop over the frequencies and calculate the crystal dielectric for each frequency
        if self.notebook.pool is None:
            self.notebook.start_pool()
        for dielecv in self.notebook.pool.imap(dielectricFunction, vs_cm1, chunksize=20):
            self.crystal_permittivity.append(dielecv)
            self.notebook.progressbars_update()
        QCoreApplication.processEvents()
        self.calculation_required = False
        logger.debug("Finished:: calculate ")
        return

    def get_crystal_permittivity(self,vs_cm1):
        """Return the crystal permittivity.

        Parameters
        ----------
        vs_cm1 : list of floats
            The frequencies in cm-1

        Returns
        -------
        self.crystal_permittivity : array_like
            The crystal permittivity (a 3x3 tensor) at each frequency

        """
        logger.debug(f"Start:: get_crystal_permittivity {self.refresh_required}")
        if self.calculation_required or self.refresh_required or  ( len(self.vs_cm1) != len(vs_cm1) ) or ( self.vs_cm1[0] != vs_cm1[0] ) or ( self.vs_cm1[1] != vs_cm1[1] ) :
            logger.debug("get_crystal_permittivity refreshing and recalculating" )
            self.refresh()
            self.calculate(vs_cm1)
        logger.debug(f"Finished:: get_crystal_permittivity {self.refresh_required}")
        return self.crystal_permittivity


    def get_crystal_permittivity_object(self):
        """Return the crystal permittivity object.

        Parameters
        ----------
        None

        Returns
        -------
        self.CrystalPermittivityObject : a permittivity object
            An instance of dielectric function

        """
        logger.debug(f"Start:: get_crystal_permittivity_object {self.refresh_required}")
        if self.calculation_required or self.refresh_required:
            logger.debug("get_crystal_permittivity_object refreshing and recalculating" )
            self.refresh()
            self.calculate(self.vs_cm1)
        logger.debug(f"Finished:: get_crystal_permittivity_object {self.refresh_required}")
        return self.CrystalPermittivityObject
