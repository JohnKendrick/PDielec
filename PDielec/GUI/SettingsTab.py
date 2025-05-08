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

import numpy as np
from qtpy.QtCore import QCoreApplication, QSize, Qt
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from PDielec import Calculator, DielectricFunction
from PDielec.Constants import amu, angstrom, average_masses, isotope_masses, wavenumber
from PDielec.Utilities import Debug


class FixedQTableWidget(QTableWidget):
    """A custom QTableWidget with fixed number of rows and columns.

    This widget is designed to have a fixed size based on a specified number of rows and columns, rather than adjusting dynamically to the content it displays. If the number of rows and columns isn't specifically provided, it defaults to the current count of rows and columns in the widget.

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
    sizeHint()
        Calculates and returns the recommended size for the table widget based on the fixed number of rows and columns, or based on its current row and column counts if no fixed numbers are provided.

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

    def sizeHint(self):
        """Calculate the suggested dimensions for the widget based on its content.

        This method calculates the preferred size of the widget by summing up the widths of all columns and the heights of all rows, including any scrollbars, headers, and frame widths.

        Parameters
        ----------
        None

        Returns
        -------
        QSize
            The calculated size that suggests the optimal dimensions for the widget. The width is determined by the total width of all columns, the vertical header width, the vertical scrollbar width, and twice the frame width. The height is determined by the total height of all rows, the horizontal header height, the horizontal scrollbar height, and twice the frame width.

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

    This class allows for the configuration of various simulation parameters, such as Eckart conditions, Born charge neutrality, atomic mass definitions, and optical permittivity values. It also enables the calculation of infrared intensities, dielectric functions, and permittivity values based on these settings.

    Parameters
    ----------
    parent : QWidget
        The parent widget (container) of this widget.
    debug : bool, optional
        Flag indicating whether debugging messages should be printed. Defaults to False.

    Attributes
    ----------
    notebook : parent
        The notebook widget that contains this settings tab.
    refreshRequired : bool
        A flag indicating whether the displayed information needs to be refreshed.
    calculationRequired : bool
        A flag indicating whether new calculations are needed based on the current settings.
    settings : dict
        A dictionary of settings related to calculations. Includes flags for Eckart conditions, Born charge neutrality, mass definition, and sigma values.
    mass_definition_options : list
        A list of options for defining atomic masses.
    masses_dictionary : dict
        A dictionary where keys are element symbols and values are the corresponding atomic masses.
    modes_selected : list
        A list indicating which vibrational modes are selected for calculations.
    frequencies_cm1 : list
        A list of vibrational frequencies in cm^-1.
    frequencies_have_been_edited : bool
        A flag indicating if frequencies have been manually edited.
    intensities : list
        A list of calculated infrared intensities for the vibrational modes.
    sigmas_cm1 : list
        A list of Lorentzian width factors (σ) in cm^-1 for the vibrational modes.
    oscillator_strengths : list
        A list of oscillator strengths calculated for the vibrational modes.
    mass_weighted_normal_modes : np.ndarray
        An array of mass-weighted normal modes.
    CrystalPermittivityObject : DielectricFunction
        An object for calculating dielectric functions.
    vs_cm1 : list
        A list of frequencies in cm^-1 at which the crystal permittivity is calculated.
    crystal_permittivity : list
        A list of calculated crystal permittivity values.
    recalculate_selected_modes : bool
        A flag indicating whether selected modes should be recalculated.
    reader : object
        A reader object from the parent notebook, used for accessing and manipulating input data.

    Methods
    -------
    setElementMass(element, mass)
        Sets the mass for a specified element.
    createIntensityTable()
        Creates a table displaying calculated intensities for selected vibrational modes.
    requestRefresh()
        Requests a refresh of displayed data.
    writeSpreadsheet()
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
    calculate(vs_cm1)
        Calculates the permittivity of the crystal over the specified range of frequencies.
    getCrystalPermittivity(vs_cm1)
        Returns the crystal permittivity for the specified frequencies.
    getCrystalPermittivityObject()
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
        This method sets up various elements of the interface, including checkboxes, comboBoxes, spinBoxes, and tableWidgets with predefined settings and connects various signals and slots. Settings related to the Eckart flag, neutral Born charges, sigma values, mass definition, and optical permittivity are initialized. Additionally, UI elements for adjusting these settings are created. It involves initializing debugging capabilities, setting up a notebook reference, and various flags and settings necessary for the functionality of the application.

        """        
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug, "SettingsTab:")
        debugger.print("Start:: initialisation")
        self.notebook = parent
        self.refreshRequired = True
        self.calculationRequired = False
        self.settings = {}
        self.settings["Eckart flag"] = True
        self.settings["Neutral Born charges"] = False
        self.settings["Sigma value"] = 5
        self.mass_definition_options = ["average","program","isotope","gui"]
        self.settings["Mass definition"] = "average"
        self.settings["Optical permittivity"] = None
        self.settings["Optical permittivity edited"] = False
        self.masses_dictionary = {}
        self.modes_selected = []
        self.frequencies_cm1 = []
        self.frequencies_have_been_edited = False
        self.intensities = []
        self.sigmas_cm1 = []
        self.oscillator_strengths = []
        self.mass_weighted_normal_modes = None
        self.CrystalPermittivityObject = None
        self.vs_cm1 = []
        self.crystal_permittivity = []
        self.recalculate_selected_modes = True
        # get the reader from the main tab
        self.reader = self.notebook.reader
        # Create second tab - SettingsTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The eckart checkbox
        #
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
        #
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
        #
        self.mass_cb = QComboBox(self)
        self.mass_cb.setToolTip("The atomic masses used to calculate frequencies and intensities can be give here")
        self.mass_cb.addItem("Average natural abundance")
        self.mass_cb.addItem("Mass taken from QM/MM program")
        self.mass_cb.addItem("Most common isotope mass")
        self.mass_cb.addItem("Masses set individually")
        self.mass_cb.model().item(3).setEnabled(False)
        # set default to average natural abundance
        self.mass_cb.activated.connect(self.on_mass_cb_activated)
        self.current_mass_definition_index = self.mass_definition_options.index(self.settings["Mass definition"])
        self.mass_cb.setCurrentIndex(0)
        form.addRow(QLabel("Atomic mass defintion:", self), self.mass_cb)
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
        #
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
        # Create the Optical permittivity table widget and block signals until a click on the widget
        #
        self.optical_tw = FixedQTableWidget(3,3,parent=self)
        self.optical_tw.setToolTip("The optical permittivity is taken from the calculation where this is possible.  If it is not availble suitbale values should be provided here")
        # Set the header names
        self.optical_tw.setHorizontalHeaderLabels(["x","y","z"])
        self.optical_tw.setVerticalHeaderLabels  (["x","y","z"])
        self.optical_tw.itemClicked.connect(self.on_optical_tw_itemClicked)
        self.optical_tw.itemChanged.connect(self.on_optical_tw_itemChanged)
        # Block the widget until the optical permittivity is loaded
        self.optical_tw.blockSignals(True)
        self.optical_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel("Optical permittivity:", self), self.optical_tw)
        vbox.addLayout(form)
        # output window
        # Create Table containing the IR active modes
        self.output_tw = FixedQTableWidget(parent=self)
        self.output_tw.setToolTip("Output showing the frequencies and strengths of the IR active modes only")
        self.output_tw.itemChanged.connect(self.on_output_tw_itemChanged)
        self.output_tw.setRowCount(1)
        self.output_tw.blockSignals(True)
        vbox.addWidget(self.output_tw)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        debugger.print("Finished:: initialisation")

    def setElementMass(self,element,mass):
        """Set the mass value of a specific element in the GUI.

        This method updates the mass value for a specified element in the graphical user interface (GUI). It also triggers various updates in the GUI to reflect this change, including updating the mass selection and recalculating modes if required.

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
        - This method assumes that `self.settings`, `self.masses_dict`, `self.mass_cb`, and other related attributes are already defined in the class.
        - It triggers a refresh of the GUI and potentially recalculates selected modes based on the new mass settings.

        """        
        debugger.print("Start::  setElementMass",element,mass)
        self.settings["Mass definition"] = "gui"
        self.masses_dict[element] = mass
        self.mass_cb.setCurrentIndex(3)
        self.set_masses_tw()
        self.refreshRequired = True
        self.recalculate_selected_modes = True
        self.refresh()
        debugger.print("Finished::  setElementMass",element,mass)

    def createIntensityTable(self):
        """Generate the intensity table for spectroscopy analysis.

        This method is responsible for creating a table that includes calculated intensities, frequencies, oscillator strengths, 
        and other related parameters required for spectroscopic analysis. It adapts based on reader settings, mass definitions, 
        and whether to recalculate selected modes.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method relies on several external variables and settings, including:
        - reader settings for neutralizing or resetting born charges
        - mass definition settings to determine how masses are handled in calculations
        - a check on if frequencies have been edited, which affects reading directly from the reader or not
        - uses `Calculator` class methods to compute normal modes, oscillator strengths, infrared intensities, and ionic permittivity
        - the method also determines the dielectric function based on the configuration and recalculates selected modes if necessary
        - it finally populates a provided table widget (`self.output_tw`) with computed spectroscopic properties

        See Also
        --------
        - `Calculator.normal_modes`: For calculation of normal modes.
        - `Calculator.oscillator_strengths`: For calculating oscillator strengths based on normal modes and born charges.
        - `Calculator.infrared_intensities`: For the calculation of infrared intensities from oscillator strengths.
        - `Calculator.ionic_permittivity`: For calculating the ionic permittivity based on mode list, oscillator strengths, and frequencies.

        """        
        debugger.print("Start:: createIntensityTable")
        self.reader = self.notebook.reader
        # Only calculate if the reader is set
        if self.reader is None:
            debugger.print("createIntensityTable aborting as now reader available")
            return
        if self.settings["Neutral Born charges"]:
            self.reader.neutralise_born_charges()
        else:
            self.reader.reset_born_charges()
        self.reader.eckart = self.settings["Eckart flag"]
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
            print("Error unkown mass definition", self.settings["Mass definition"] )
        QCoreApplication.processEvents()
        self.mass_weighted_normal_modes = self.reader.calculate_mass_weighted_normal_modes()
        # convert cm-1 to au
        if not self.frequencies_have_been_edited:
            self.frequencies_cm1 = self.reader.frequencies
        if len(self.sigmas_cm1) == 0:
            self.sigmas_cm1 = [ self.settings["Sigma value"] for i in self.frequencies_cm1 ]
        born_charges = np.array(self.reader.born_charges)
        if self.reader.type == "Experimental output":
            self.oscillator_strengths = self.reader.oscillator_strengths
        else:
            #
            # calculate normal modes in xyz coordinate space
            masses = np.array(self.reader.masses) * amu
            normal_modes = Calculator.normal_modes(masses, self.mass_weighted_normal_modes)
            # from the normal modes and the born charges calculate the oscillator strengths of each mode
            self.oscillator_strengths = Calculator.oscillator_strengths(normal_modes, born_charges)
        # calculate the intensities from the trace of the oscillator strengths
        self.intensities = Calculator.infrared_intensities(self.oscillator_strengths)
        # Decide which modes to select
        if self.recalculate_selected_modes and len(self.intensities) > 0 and len(self.frequencies_cm1) > 0:
            debugger.print("createIntensityTable: recalculating selected modes")
            self.modes_selected = []
            self.mode_list = []
            for f,intensity in zip(self.frequencies_cm1,self.intensities):
                if f > 10.0 and intensity > 1.0E-6:
                    self.modes_selected.append(True)
                else:
                    self.modes_selected.append(False)
            self.mode_list = [i for i,mode in enumerate(self.modes_selected) if mode]
            debugger.print("Selected modes are;",self.mode_list)
            self.recalculate_selected_modes = False
        # end if
        #
        # Calculate the ionic contribution to the permittivity
        frequencies_au = wavenumber*np.array(self.frequencies_cm1)
        volume_au = self.reader.volume*angstrom*angstrom*angstrom
        self.epsilon_ionic = Calculator.ionic_permittivity(self.mode_list, self.oscillator_strengths, frequencies_au, volume_au )
        # Make an np array of epsilon infinity
        epsilon_inf = np.array(self.settings["Optical permittivity"])
        #
        # If the reader already has a Dielectric Constant then use this
        #
        if self.reader.CrystalPermittivity:
            self.CrystalPermittivityObject = self.reader.CrystalPermittivity
        else:
            # Otherwise create one from the data that has been read in
            drude = False
            drude_plasma_au = 0
            drude_sigma_au = 0
            sigmas_au = np.array(self.sigmas_cm1)*wavenumber
            debugger.print("CreateIntensityTable: Calculating dielectric",self.mode_list)
            self.CrystalPermittivityObject = DielectricFunction.DFT(
                                         self.mode_list, frequencies_au, sigmas_au, self.oscillator_strengths,
                                         volume_au, drude, drude_plasma_au, drude_sigma_au )
        # Add the optical permittivity to the dielctric function
        self.CrystalPermittivityObject.setEpsilonInfinity(epsilon_inf)
        #
        # Prepare to finish
        #
        self.output_tw.setRowCount(len(self.sigmas_cm1))
        self.output_tw.setColumnCount(5)
        self.output_tw.setHorizontalHeaderLabels(["   Sigma   \n(cm-1)", " Frequency \n(cm-1)", "  Intensity  \n(Debye2/Å2/amu)", "Integrated Molar Absorption\n(L/mole/cm2)", "Absorption maximum\n(L/mole/cm)"])
        QCoreApplication.processEvents()
        self.redraw_output_tw()
        QCoreApplication.processEvents()
        # if self.notebook.spreadsheet is not None:
        #     self.writeSpreadsheet()
        QCoreApplication.processEvents()
        debugger.print("Finished:: createIntensityTable")
        return

    def requestRefresh(self):
        """Request to refresh the environment or state.

        This method sets a flag to indicate that a refresh is required, logging the start and end of its execution.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: requestRefresh")
        self.refreshRequired = True
        # self.refresh()
        debugger.print("Finished:: requestRefresh")
        return

    def writeSpreadsheet(self):
        """Write the configuration and calculation results to the associated spreadsheet.

        This function iterates through various settings and results stored within the object, such as optical permittivity, mass definitions, and vibrational mode calculations. For each setting or result, it writes the relevant information to a new row in the spreadsheet. The spreadsheet is expected to have a 'Settings' worksheet, which this function selects and updates with the provided values.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: writeSpreadsheet")
        sp = self.notebook.spreadsheet
        if sp is None:
            debugger.print("Finished:: writeSpreadsheet - Aborting write of spreadsheet")
            return
        debugger.print("Writing of spreadsheet")
        sp.selectWorkSheet("Settings")
        sp.delete()
        sp.writeNextRow(["Settings and calculations of frequencies and absorption"], row=0, col=1)
        for item in sorted(self.settings):
            if item == "Optical permittivity" and self.settings[item] is not None:
               sp.writeNextRow([item], col=1)
               sp.writeNextRow(self.settings[item][0], col=2, check=1)
               sp.writeNextRow(self.settings[item][1], col=2, check=1)
               sp.writeNextRow(self.settings[item][2], col=2, check=1)
               eps = self.epsilon_ionic.tolist()
               sp.writeNextRow(["Ionic Contribution to permittivity at zero frequency"], col=1)
               sp.writeNextRow(eps[0], col=2, check=1)
               sp.writeNextRow(eps[1], col=2, check=1)
               sp.writeNextRow(eps[2], col=2, check=1)
               eps = self.epsilon_ionic + np.array(self.settings[item])
               eps = eps.tolist()
               sp.writeNextRow(["Total Contribution to permittivity at zero frequency"], col=1)
               sp.writeNextRow(eps[0], col=2, check=1)
               sp.writeNextRow(eps[1], col=2, check=1)
               sp.writeNextRow(eps[2], col=2, check=1)
            elif item == "Mass definition":
               sp.writeNextRow([item,self.settings[item]], col=1, check=1)
               for c in self.masses_dictionary:
                   sp.writeNextRow([c,self.masses_dictionary[c]],col=3,check=1)
            else:
               sp.writeNextRow([item,self.settings[item]], col=1, check=1)
        sp.writeNextRow([""], col=1)
        sp.writeNextRow(["Mode","Include?","Sigma(cm-1)","Frequency(cm-1)","Intensity(Debye2/Angs2/amu","Integrated Molar Absorption(/L/mole/cm2","Absorption maximum (L/mole/cm)"], col=1)
        for mode,(f,intensity,sigma,selected) in enumerate(zip(self.frequencies_cm1, self.intensities, self.sigmas_cm1,self.modes_selected)):
            yn = "No"
            if selected:
                yn = "Yes"
            sp.writeNextRow([mode, yn, sigma, f, intensity, 4225.6*intensity, 2*4225.6*intensity/sigma/np.pi], col=1)
        debugger.print("Finished:: writeSpreadsheet")

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
        - This method assumes that `self.frequencies_cm1`, `self.sigmas_cm1`, `self.intensities`, and `self.modes_selected` contain iterable values that are indexed correspondingly.
        - The method temporarily blocks signals from `self.output_tw` to prevent unwanted signal handling during updates.
        - Table item flags are set to control their selectability, checkability, editability, and enablement based on whether the corresponding mode is selected.
        - After updating, it ensures the table columns fit the content and re-enables the signals.
        - Uses `Qt` enumeration for setting check state, item flags, and text alignment.

        """        
        debugger.print("Start:: redraw_output_tw")
        # If the frequencies haven't been set yet just don't try to do anything
        self.output_tw.blockSignals(True)
        for i,(f,sigma,intensity) in enumerate(zip(self.frequencies_cm1, self.sigmas_cm1, self.intensities)):
            # Sigma and check / unchecked column
            items = []
            itemFlags = []
            item = QTableWidgetItem(f"{sigma:.1f}")
            if self.modes_selected[i]:
                item.setCheckState(Qt.Checked)
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable )
                freqFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
                otherFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled
            else:
                #itemFlags.append( item.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled & ~Qt.ItemIsSelectable & ~Qt.ItemIsEditable )
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled )
                freqFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
                item.setCheckState(Qt.Unchecked)
                freqFlags = item.flags() & Qt.NoItemFlags 
                otherFlags = item.flags() & Qt.NoItemFlags
            items.append(item)
            # Frequency column cm-1
            items.append(QTableWidgetItem(f"{f:.4f}" ) )
            itemFlags.append( freqFlags )
            # Intensity column Debye2/Angs2/amu
            items.append(QTableWidgetItem(f"{intensity:.4f}" ) )
            itemFlags.append( otherFlags )
            # Integrated molar absorption L/mole/cm/cm
            items.append(QTableWidgetItem(f"{intensity*4225.6:.2f}" ) )
            itemFlags.append( otherFlags )
            # Maximum extinction L/mole/cm
            items.append(QTableWidgetItem(f"{2*intensity*4225.6/self.sigmas_cm1[i]/np.pi:.2f}" ) )
            itemFlags.append( otherFlags )
            for j,(item,flag) in enumerate(zip(items,itemFlags)):
                item.setFlags(flag)
                item.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                self.output_tw.setItem(i, j, item )
        # Release the block on signals for the frequency output table
        self.output_tw.resizeColumnsToContents()
        self.output_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print("Finished:: redraw_output_tw")

    def on_sigma_changed(self):
        """Update the sigma value in settings and apply it across the frequency range, then refreshes the output.

        This method is triggered when the sigma value has been changed. It updates the stored sigma value in the settings,
        applies the new sigma value across the pre-defined frequency range, redraws the output table widget, and then
        initiates a refresh to reflect the newly applied changes.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: redraw_output_tw")
        self.settings["Sigma value"] = self.sigma_sb.value()
        self.sigmas_cm1 = [ self.settings["Sigma value"] for i in self.frequencies_cm1 ]
        self.redraw_output_tw()
        debugger.print("on sigma change ", self.settings["Sigma value"])
        self.refreshRequired = True
        self.refresh()
        debugger.print("Finished:: on_sigma_changed")

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
        debugger.print("Start:: on_mass_combobox_activated", self.mass_cb.currentText())
        self.settings["Mass definition"] = self.mass_definition_options[index]
        self.current_mass_definition_index = index
        if index < 3:
            self.mass_cb.model().item(3).setEnabled(False)
        # Modify the element masses
        self.set_masses_tw()
        self.refreshRequired = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print("Finished:: on_mass_combobox_activated", self.mass_cb.currentText())

    def set_masses_tw(self):
        """Set the element masses in the table widget based on the mass_definition setting.

        This method updates the element masses table widget with masses according to the selected mass definition mode ('average', 'program', 'isotope', 'gui'). It retrieves species information from the reader attribute, calculates or fetches the corresponding masses, and updates an internal dictionary and the table widget accordingly. Signals from the table widget are temporarily blocked to prevent unwanted event triggers during updates.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: set_masses_tw")
        if self.reader:
            self.element_masses_tw.blockSignals(True)
            species = self.reader.getSpecies()
            # set the initial dictionary according to the mass_definition
            masses = []
            if self.settings["Mass definition"] == "average":
                self.mass_cb.setCurrentIndex(0)
                for element in species:
                    mass = average_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings["Mass definition"] == "program":
                self.mass_cb.setCurrentIndex(1)
                self.reader.reset_masses()
                masses = self.reader.masses_per_type
                for mass,element in zip(masses,species):
                    self.masses_dictionary[element] = mass
            elif self.settings["Mass definition"] == "isotope":
                self.mass_cb.setCurrentIndex(2)
                for element in species:
                    mass = isotope_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings["Mass definition"] == "gui":
                for element in species:
                    mass = self.masses_dictionary[element]
                    masses.append(mass)
            else:
                 debugger.print("Error mass_definition not recognised", self.settings["Mass definition"])
            self.element_masses_tw.setColumnCount(len(masses))
            self.element_masses_tw.setHorizontalHeaderLabels(species)
            self.element_masses_tw.setVerticalHeaderLabels([""])
            debugger.print("masses_dictionary",self.masses_dictionary)
            debugger.print("masses",masses)
            # set masses of the elements in the table widget according to the mass definition
            for i,(mass,element) in enumerate(zip(masses,species)):
                debugger.print("set_masses_tw", self.settings["Mass definition"],i,mass,element)
                qw = QTableWidgetItem()
                if self.settings["Mass definition"] == "program":
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{mass:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    self.element_masses_tw.setItem(0,i, qw )
                elif self.settings["Mass definition"] == "average":
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{average_masses[element]:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    debugger.print("average",average_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings["Mass definition"] == "isotope":
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{isotope_masses[element]:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    debugger.print("isotope",isotope_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings["Mass definition"] == "gui":
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(f"{self.masses_dictionary[element]:.6f}")
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    debugger.print("gui",self.masses_dictionary[element])
                    self.element_masses_tw.setItem(0,i, qw )
                else:
                    print("Mass definition not processed", self.settings["Mass definition"])
            # unblock the table signals
            self.element_masses_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print("Finished:: set_masses_tw")

    def on_output_tw_itemChanged(self, item):
        """Handle item changes in the output table widget.

        This method is called when an item in the output table widget (`output_tw`) is changed. It manages the selection state of modes based on the state of checkbox items in the first column, updates the "sigmas_cm1" and "frequencies_cm1" lists based on the user input in the respective columns, refreshes the widget, and processes pending GUI events.

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was changed.

        Returns
        -------
        None

        Notes
        -----
        - This method blocks signals from the output table widget at the beginning to prevent recursive calls during its execution.
        - It distinguishes between changes in the first column (selection state and mode sigma values) and the second column (frequency values).
        - The method updates internal state to reflect the changes made by the user, including which modes are selected and the values for "sigmas_cm1" and "frequencies_cm1".
        - After making changes, this method refreshes the output table widget and enforces a GUI update with `QCoreApplication.processEvents()`.

        See Also
        --------
        redraw_output_tw : Redraws the output table widget.
        refresh : Refreshes the internal state based on the latest changes.

        """        
        self.output_tw.blockSignals(True)
        debugger.print("Start:: on_output_tw_itemChanged", item.row(), item.column() )
        col = item.column()
        row = item.row()
        if col == 0:
            # If this is the first column alter the check status but reset the sigma value
            if item.checkState() == Qt.Checked:
                debugger.print("on_output_tw_itemChanged setting selected mode to True",row )
                self.modes_selected[row] = True
                self.mode_list = [i for i,mode in enumerate(self.modes_selected) if mode]
            else:
                debugger.print("on_output_tw_itemChanged setting selected mode to False",row )
                self.modes_selected[row] = False
                self.mode_list = [i for i,mode in enumerate(self.modes_selected) if mode]
            new_value = float(item.text())
            if new_value != self.sigmas_cm1:
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
        debugger.print("on_output_tw_itemChanged selected_modes",self.modes_selected )
        self.refreshRequired = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print("Finished:: on_output_tw_itemChanged")

    def on_element_masses_tw_itemClicked(self, item):
        """Handle the item clicked event on the element mass table widget.

        The only thing this function does is to make sure that signal from the `element_masses_tw` table widget are not blocked,

        Parameters
        ----------
        item : QTableWidgetItem
            The item in the table widget that was clicked.

        Returns
        -------
        None

        """        
        debugger.print("Start:: on_element_masses_tw_itemClicked)", item.row(),item.column() )
        self.element_masses_tw.blockSignals(False)
        debugger.print("Finished:: on_element_masses_tw_itemClicked)" )

    def on_element_masses_tw_itemChanged(self, item):
        """Handle item changed event in the mass table widget.

        This function is triggered whenever a cell in the masses table widget is edited. It updates the internal masses dictionary based on the new value provided in the table cell. Additionally, it sets some UI components to reflect the change and flags the need for recalculation and refresh.

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was changed.

        Returns
        -------
        None

        Notes
        -----
        - `self.reader.getSpecies()` returns a list of elements corresponding to table columns.
        - `self.masses_dictionary` is a dict where keys are element symbols and values are their masses.
        - `self.refreshRequired` and `self.recalculate_selected_modes` are boolean flags used to indicate that a refresh and a recalculation are necessary.
        - `self.refresh()` is a method that refreshes or updates the UI components as needed based on the new changes.

        """        
        debugger.print("Start:: on_element_masses_tw_itemChanged)", item.row(), item.column() )
        elements = self.reader.getSpecies()
        col = item.column()
        self.mass_cb.model().item(3).setEnabled(True)
        self.settings["Mass definition"] = "gui"
        self.mass_cb.setCurrentIndex(3)
        self.masses_dictionary[elements[col]] = float(item.text())
        debugger.print("masses_dictionary", self.masses_dictionary)
        self.refreshRequired = True
        self.recalculate_selected_modes = True
        self.refresh()
        debugger.print("Finished:: on_element_masses_tw_itemChanged)" )

    def on_optical_tw_itemChanged(self, item):
        """Handle item change events for an optical permittivity table.

        This method updates the optical permittivity settings based on changes made to the table widget items. It enforces symmetry in the optical permittivity matrix, marks the settings as edited, and triggers a series of updates and recalculations.

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was changed. The item's row, column, and text are used to update the settings.

        Returns
        -------
        None

        Notes
        -----
        - Updates the 'Optical permittivity' matrix in the settings dictionary based on the row, column, and text of the changed item.
        - Enforces symmetry by setting the [column][row] element equal to the [row][column] element.
        - Marks the 'Optical permittivity edited' flag as True.
        - Calls `refresh_optical_permittivity_tw` to refresh the optical permittivity table widget.
        - Flags `refreshRequired`, `recalculate_selected_modes`, and `refresh` methods/attributes for further actions.
        - Processes any pending events with `QCoreApplication.processEvents()` to ensure the UI remains responsive.

        """        
        debugger.print("Start::on_optical_itemChanged)", item.row(), item.column(), item.text() )
        self.settings["Optical permittivity"][item.row()][item.column()] = float(item.text())
        self.settings["Optical permittivity"][item.column()][item.row()] = float(item.text())
        self.settings["Optical permittivity edited"] = True
        self.refresh_optical_permittivity_tw()
        self.refreshRequired = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print("Finished::on_optical_itemChanged)")
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
        debugger.print("Start:: on_optical_itemClicked)", item.row(), item.column() )
        self.optical_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print("Finished:: on_optical_itemClicked)")
        return

    def refresh(self, force=False):
        """Refresh the current state based on notebook content changes or user request.

        This method updates the GUI components such as tables for masses, optical permittivity, sigma value, Eckart term checkbox, and Born charges checkbox based on the changes made in the notebook or when forced by the user. It also triggers refreshes for any associated scenarios, plotting, analysis, and fitting tabs. The method blocks all signals from widgets during the update process to prevent unwanted calls to event handlers.

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
        debugger.print("Start:: refresh", force )
        if not self.reader and self.notebook.reader:
            self.refreshRequired = True
        if not self.refreshRequired and not force:
            debugger.print("Finished:: refresh not required",force)
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
                self.set_optical_permittivity_tw()
            else:
                self.refresh_optical_permittivity_tw()
        self.sigma_sb.setValue(self.settings["Sigma value"])
        if self.settings["Eckart flag"]:
            self.eckart_cb.setCheckState(Qt.Checked)
        else:
            self.eckart_cb.setCheckState(Qt.Unchecked)
        if self.settings["Neutral Born charges"]:
            self.born_cb.setCheckState(Qt.Checked)
        else:
            self.born_cb.setCheckState(Qt.Unchecked)
        #
        # Unlock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.createIntensityTable()
        self.refreshRequired = False
        # Refresh the scenarios that depend on this dielectric
        if self.notebook.scenarios is not None:
            for scenario in self.notebook.scenarios:
                scenario.requestRefresh()
        if self.notebook.plottingTab is not None:
            self.notebook.plottingTab.requestRefresh()
        if self.notebook.analysisTab is not None:
            self.notebook.analysisTab.requestRefresh()
        if self.notebook.fitterTab is not None:
            self.notebook.fitterTab.requestRefresh()
        QCoreApplication.processEvents()
        self.refreshRequired = False
        self.calculationRequired = True
        debugger.print("Finished:: refresh", force )
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
            - `self.settings` is a dict containing the key 'Optical permittivity' which maps to a 2D list (i.e., list of lists) of numeric values.
            - `self.optical_tw` is a Qt TableWidget or similar object with a `.setItem()` method accepting row, column, and QTableWidgetItem arguments.
            - The table widget's signals are temporarily blocked to prevent unintended side-effects during the update.
            - `Qt.AlignHCenter | Qt.AlignVCenter` is used to center-align the text, assuming `Qt` and `QTableWidgetItem` are properly imported and used.
            - The `QCoreApplication.processEvents()` call forces the application to process all pending events, ensuring the UI updates in real-time.

        """        
        debugger.print("Start:: refresh_optical_permittivity")
        optical = self.settings["Optical permittivity"]
        self.optical_tw.blockSignals(True)
        for i,row in enumerate(optical):
            for j, value in enumerate(row):
                qw = QTableWidgetItem(f"{value:.4f}" )
                qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                self.optical_tw.setItem(i,j,qw)
        self.optical_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print("Finished:: refresh_optical_permittivity")
        return

    def set_optical_permittivity_tw(self):
        """Set the optical permittivity in the current settings.

        This method sets the 'Optical permittivity' value in the settings to the value of `zerof_optical_dielectric` from the reader, triggers a refresh for the optical permittivity widget, and marks the selected modes for recalculation. It also signals that a refresh is required in the application and ensures the GUI events are processed to keep the application responsive.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: set_optical_permittivity_tw")
        self.settings["Optical permittivity"] = self.reader.zerof_optical_dielectric
        self.refresh_optical_permittivity_tw()
        self.recalculate_selected_modes = True
        self.refreshRequired = True
        QCoreApplication.processEvents()
        debugger.print("Finished:: set_optical_permittivity_tw")
        return

    def on_born_changed(self):
        """Handle changes in the Born checkbox state.

        This method is triggered when the state of the Born checkbox changes. It updates the settings to reflect the new state, marks the session for refresh and recalculation, and then refreshes the UI.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: on_born_change ", self.born_cb.isChecked())
        self.settings["Neutral Born charges"] = self.born_cb.isChecked()
        debugger.print("on born change ", self.settings["Neutral Born charges"])
        self.refreshRequired = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print("Finished:: on_born_change ", self.born_cb.isChecked())
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
        debugger.print("Start:: on_eckart_change ", self.eckart_cb.isChecked())
        self.settings["Eckart flag"] = self.eckart_cb.isChecked()
        debugger.print("on eckart change ", self.settings["Eckart flag"])
        self.refreshRequired = True
        self.recalculate_selected_modes = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print("Finished:: on_eckart_change ", self.eckart_cb.isChecked())
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
        debugger.print("Start:: calculate ")
        if len(vs_cm1) == 0:
            debugger.print("Finished:: calculate aborted vs_cm1 has not been set")
            return
        self.vs_cm1 = vs_cm1.copy()
        dielectricFunction = self.CrystalPermittivityObject.function()
        self.crystal_permittivity = []
        debugger.print("About to calculate settings crystal dielectric using pool")
        # Loop over the frequencies and calculate the crystal dielectric for each frequency
        if self.notebook.pool is None:
            self.notebook.startPool()
        for dielecv in self.notebook.pool.imap(dielectricFunction, vs_cm1, chunksize=20):
            self.crystal_permittivity.append(dielecv)
            self.notebook.progressbars_update()
        QCoreApplication.processEvents()
        self.calculationRequired = False
        debugger.print("Finished:: calculate ")
        return

    def getCrystalPermittivity(self,vs_cm1):
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
        debugger.print("Start:: get_crystal_permittivity", self.refreshRequired)
        if self.calculationRequired or self.refreshRequired or  ( len(self.vs_cm1) != len(vs_cm1) ) or ( self.vs_cm1[0] != vs_cm1[0] ) or ( self.vs_cm1[1] != vs_cm1[1] ) :
            debugger.print("get_crystal_permittivity refreshing and recalculating" )
            self.refresh()
            self.calculate(vs_cm1)
        debugger.print("Finished:: get_crystal_permittivity", self.refreshRequired)
        return self.crystal_permittivity


    def getCrystalPermittivityObject(self):
        """Return the crystal permittivity object.

        Parameters
        ----------
        None

        Returns
        -------
        self.CrystalPermittivityObject : a permittivity object
            An instance of dielectric function

        """
        debugger.print("Start:: getCrystalPermittivityObject", self.refreshRequired)
        if self.calculationRequired or self.refreshRequired:
            debugger.print("getCrystalPermittivityObject refreshing and recalculating" )
            self.refresh()
            self.calculate(self.vs_cm1)
        debugger.print("Finished:: getCrystalPermittivityObject", self.refreshRequired)
        return self.CrystalPermittivityObject

