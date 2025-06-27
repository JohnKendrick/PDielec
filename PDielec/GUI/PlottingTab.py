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
"""PlottingTab module."""

# Import plotting requirements
import matplotlib
import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from qtpy.QtCore import QCoreApplication, Qt
from qtpy.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from PDielec.Constants import avogadro_si
from PDielec.Utilities import Debug

possible_frequency_units = ["wavenumber","THz","GHz","ang","nm","um","mm","cm","m"]

def isThisAFrequencyUnit(unit):
    """Return true if this is a frequency unit, false if a wavelength.

    Units of frequency are 'wavenumber','THz','GHz'
    Units of wavelength are 'ang','nm','um','mm','cm' or 'm'

    Returns
    -------
    bool
        True if this is a frequency unit, False if a wavelength.

    """
    index = possible_frequency_units.index(unit)
    return index <= 2

def convert_frequency_units( value, unit_in, unit_out ):
    """Convert between frequency and wavelength units.

    The input can be either a scalar value or a numpy array of values. The function will return the converted value(s) in the output units specified.
    The unit strings are turned into lower-case so case is irrelevant

    Parameters
    ----------
    value : scalar or numpy array
        The value(s) for which the conversion is to be made.
    unit_in : str
        The units of the input value(s). Can be one of 'cm-1' (or 'wavenumber'), 'GHz', 'THz', 'nm', 'um', 'mm', 'cm', 'm'.
    unit_out : str
        The units of the output value(s). Must be one of 'cm-1' (or 'wavenumber'), 'GHz', 'THz', 'nm', 'um', 'mm', 'cm', 'm'.

    Returns
    -------
    scalar or numpy array
        The converted value(s) in the output units specified.

    """
    unit_in  = unit_in.lower()
    unit_out = unit_out.lower()
    # the conversion dictionary has a tuple for each unit
    # the first entry is the multiplicative factor, the second is the operator
    # To convert wavelength they have to be scaled and then inverted
    wavenumber = { "cm-1"       : (1,"*"),
                   "wavenumber" : (1,"*"),
                   "thz"        : (33.356,"*"),
                   "ghz"        : (33.356E+3,"*"),
                   "ang"        : (1E-8,"/"),
                   "nm"         : (1E-7,"/"),
                   "um"         : (1E-4,"/"),
                   "mm"         : (1E-1,"/"),
                   "cm"         : (1.0 ,"/"),
                   "m"          : (1E2 ,"/"),
                 }
    # Convert the input unit to cm-1
    if isinstance(value,np.ndarray):
        for i,x in enumerate(value):
            if x <= 0:
                # Problem with negative or zero conversions between wavelength and frequency
                print("Warning a zero or negative frequency/wavelength is not permitted",i,x,unit_in,unit_out)
                x[i] = 1.0E-8
    elif value <= 0 :
        value = 1.0E-8
    scale,operator = wavenumber[unit_in]
    value = scale * value     
    if operator == "/":
        value = 1.0 / value  
    # convert the internal value from cm-1 to the output unit
    scale,operator = wavenumber[unit_out]
    if operator == "/":
        value = 1.0 / value  
    return value / scale
    
class PlottingTab(QWidget):
    """A class used for creating and managing a plotting tab in a graphical user interface. It inherits from QWidget.

    This class is responsible for handling plotting functionalities like choosing plot types, setting molar definitions and frequency values, managing interactions with UI components such as spin boxes, combo boxes, and buttons, and plotting the data using matplotlib.

    Parameters
    ----------
    parent : QWidget
        The parent widget.
    debug : bool, optional
        Enables debug mode which provides additional logging details. Defaults to False.

    Attributes
    ----------
    settings : dict
        A dictionary containing settings related to plotting such as minimum frequency, maximum frequency, frequency increment, molar definitions, number of atoms, plot type, and frequency unit.
    refreshRequired : bool
        A flag indicating whether the plot needs to be refreshed.
    subplot : matplotlib subplot object
        The subplot used for plotting data.
    vmin, vmax, vinc : float
        Variables to store current minimum frequency, maximum frequency, and frequency increment values for the plot.
    molar_cb_current_index : int
        Tracks the current index of the molar definition combo box to handle changes.
    notebook : QWidget
        A reference to the parent widget which contains the widget structure.
    reader : object
        An object to read and handle data from different file formats or sources.
    legends, vs_cm1 : list
        Lists used for storing legends for the plots and frequency values respectively.
    frequency_length : int
        Stores the length of the frequency array to check for changes.
    progressbar : QProgressBar
        A progress bar widget to display operation progress.

    Methods
    -------
    get_total_number_of_frequency_calculations()
        Calculates the total number of frequency calculations required.
    requestRefresh()
        Requests a refresh for the plotting.
    requestScenarioRefresh()
        Requests a refresh for all scenarios in the application.
    on_vmin_changed(), on_vmax_changed(), on_vinc_changed(value), on_funits_cb_activated(index), on_molar_cb_activated(index), on_natoms_changed(value), on_plot_type_cb_activated(index)
        Event handlers for UI component changes.
    refresh(force=False)
        Refreshes the plot based on current settings and data.
    plot()
        Generates and displays the plot based on the selected plot type and data.
    get_number_of_califications_required()
        Returns the number of calculations needed for updating the plot.
    greyed_out()
        Greys out options in the UI that are not available.
    writeSpreadsheet()
        Writes results to a spreadsheet.
    set_concentrations()
        Sets the concentrations for the plotting based on the molar definition.

    """

    def __init__(self, parent, debug=False ):
        """Initialise the Plotting tab within a given parent.

        Parameters
        ----------
        parent : QWidget
            The parent widget in which this plotting tab will be initialized.
        debug : bool, optional
            Determines whether debugging messages should be shown, by default False.

        Attributes
        ----------
        settings : dict
            A dictionary to store settings such as frequency range, molar definitions, etc.
        refreshRequired : bool
            Indicates whether the plot needs to be refreshed.
        subplot : NoneType or matplotlib subplot
            Placeholder for a future subplot object, initialized as None.
        molar_definitions : list
            A list of possible definitions for a mole, including 'Unit cells', 'Atoms', 'Molecules'.
        legends : list
            A list to store legend entries.
        vs_cm1 : list
            A list to store frequency data.
        frequency_length : int
            Stores the length of `vs_cm1` list.
        vmin, vmax, vinc : float
            Variables to store minimum, maximum, and increment values for frequency/wavelength after conversion.
        molar_cb_current_index : int
            Tracks the currently selected index in the molar definitions combo box.
        notebook : QWidget
            A reference to the parent widget, here it's assumed to be a notebook containing multiple tabs.
        reader : Object
            An object that deals with reading data, accessed through the `notebook`'s mainTab.

        Notes
        -----
        Widget elements such as QDoubleSpinBox, QComboBox, QLabel, QPushButton, QProgressBar, and the matplotlib figure are also initialized and configured within this method. This involves setting up signal-slot connections for interactive behavior, tooltips for user guidance, and incorporating the matplotlib canvas and toolbar for plotting functionality. Additionally, the method is designed for setup within a Qt layout structure, ensuring proper arrangement and display of the interface components.

        """        
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,"PlottingTab")
        debugger.print("Start:: Plotting tab initialisation")
        self.settings = {}
        self.refreshRequired = True
        self.subplot = None
        self.setWindowTitle("Plotting")
        self.settings["Minimum frequency"] = 1.0
        self.settings["Maximum frequency"] = 200
        self.settings["Frequency increment"] = 0.2
        self.molar_definitions = ["Unit cells","Atoms","Molecules"]
        self.settings["Molar definition"] = "Unit cells"
        self.settings["Number of atoms"] = 1
        self.settings["Plot type"] = "Powder Molar Absorption"
        self.settings["Frequency unit"] = "wavenumber"
        # self.settings['Plot title'] = 'Plot Title'
        self.legends = []
        self.vs_cm1 = []
        self.frequency_length = len(self.vs_cm1)
        self.vmin = 0.0
        self.vmax = 0.0
        self.vinc = 0.0
        self.molar_cb_current_index = 0
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Create last tab - PlottingTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        # Preamble for frequencies / wavelengths display
        vmin = self.settings["Minimum frequency"]
        vmax = self.settings["Maximum frequency"]
        vinc = self.settings["Frequency increment"]
        # n is the number of samples
        n = len(np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc)))
        # convert units
        vmin = convert_frequency_units(vmin,"wavenumber",self.settings["Frequency unit"])
        vmax = convert_frequency_units(vmax,"wavenumber",self.settings["Frequency unit"])
        # vinc is calculated using the number of samples
        vinc = (vmax-vmin) / n
        vinc = convert_frequency_units(vinc,"wavenumber",self.settings["Frequency unit"])
        # If dealing with wavelength then swap the order of vmin and vmax
        if not isThisAFrequencyUnit(self.settings["Frequency unit"]):
            vmin, vmax = vmax, vmin
        #
        # The minimum frequency
        #
        self.vmin_sb = QDoubleSpinBox(self)
        self.vmin_sb.setRange(0.00000001,900000000000)
        self.vmin_sb.setValue(vmin)
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            self.vmin_sb.setToolTip("Set the minimum frequency to be considered)")
        else:
            self.vmin_sb.setToolTip("Set the minimum wavelength to be considered)")
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        #
        # The maximum frequency
        #
        self.vmax_sb = QDoubleSpinBox(self)
        self.vmax_sb.setRange(0.00000001,900000000000)
        self.vmax_sb.setValue(vmax)
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            self.vmax_sb.setToolTip("Set the maximum frequency to be considered)")
        else:
            self.vmax_sb.setToolTip("Set the maximum wavelength to be considered)")
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        #
        # Choose a suitable increment
        #
        self.vinc_sb = QDoubleSpinBox(self)
        self.vinc_sb.setRange(0.0000001,5000000000.0)
        self.vinc_sb.setSingleStep(0.1)
        self.vinc_sb.setDecimals(4)
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            self.vinc_sb.setToolTip("Choose an increment for the frequency when plotting")
        else:
            self.vinc_sb.setToolTip("Choose an increment for the wavelength when plotting")
        self.vinc_sb.setValue(vinc)
        self.vinc_sb.valueChanged.connect(self.on_vinc_changed)
        #
        # Set the frequency units
        #
        self.funits_cb = QComboBox(self)
        self.funits_cb.addItems( possible_frequency_units )
        self.funits_cb.activated.connect(self.on_funits_cb_activated)
        index = possible_frequency_units.index(self.settings["Frequency unit"])
        self.funits_cb.setCurrentIndex(index)
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            self.funits_cb.setToolTip("Set the frequency unit")
            self.frequency_form_label = QLabel("Frequency min, max and increment", self)
            self.frequency_form_label.setToolTip("Choose minimum, maximum and increment for frequency")
        else:
            self.funits_cb.setToolTip("Set the wavelength unit")
            self.frequency_form_label = QLabel("Wavelength min, max and increment", self)
            self.frequency_form_label.setToolTip("Choose minimum, maximum and increment for wavelength")
        hbox = QHBoxLayout()
        hbox.addWidget(self.vmin_sb)
        hbox.addWidget(self.vmax_sb)
        hbox.addWidget(self.vinc_sb)
        hbox.addWidget(self.funits_cb)
        form.addRow(self.frequency_form_label, hbox)
        #
        # Define molar quantity
        #
        self.molar_cb = QComboBox(self)
        self.molar_cb.setToolTip("Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given")
        self.molar_cb.addItems(self.molar_definitions)
        try:
            self.molar_cb_current_index = self.molar_definitions.index(self.settings["Molar definition"])
        except Exception:
            self.molar_cb_current_index = 0
            self.settings["Molar definition"] = self.molar_definitions[self.molar_cb_current_index]
        self.molar_cb.setCurrentIndex(self.molar_cb_current_index)
        self.molar_cb.activated.connect(self.on_molar_cb_activated)
        label = QLabel("Molar definition", self)
        label.setToolTip("Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given")
        form.addRow(label, self.molar_cb)
        #
        # Number of atoms in a molecule
        #
        self.natoms_sb = QSpinBox(self)
        self.natoms_sb.setToolTip("Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed")
        self.natoms_sb.setRange(1,500)
        self.natoms_sb.setValue(self.settings["Number of atoms"])
        self.natoms_sb.valueChanged.connect(self.on_natoms_changed)
        self.natoms_sb.setEnabled(False)
        label = QLabel("Number of atoms per molecule", self)
        label.setToolTip("Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed")
        form.addRow(label, self.natoms_sb)
        #
        # Final button
        #
        self.plot_type_cb = QComboBox(self)
        self.plot_type_cb.setToolTip("Choose the which data to plot")
        self.plot_types = [
                            "Powder Molar Absorption",
                            "Powder Absorption",
                            "Powder Real Permittivity",
                            "Powder Imaginary Permittivity",
                            "Powder ATR",
                            "Crystal Reflectance (P polarisation)",
                            "Crystal Reflectance (S polarisation)",
                            "Crystal Transmittance (P polarisation)",
                            "Crystal Transmittance (S polarisation)",
                            "Crystal Absorbtance (P polarisation)",
                            "Crystal Absorbtance (S polarisation)",
                          ]
        self.plot_ylabels = {
                     "Powder Molar Absorption": r"Molar Absorption Coefficient $\mathdefault{(L mole^{-1} cm^{-1})}$",
                           "Powder Absorption": r"Absorption Coefficient $\mathdefault{(cm^{-1})}$",
                    "Powder Real Permittivity": r"Real Component of Permittivity",
               "Powder Imaginary Permittivity": r"Imaginary Component of Permittivity",
                                  "Powder ATR": r"ATR absorption",
        "Crystal Reflectance (P polarisation)": r"Fraction of p-polarised reflectance",
        "Crystal Reflectance (S polarisation)": r"Fraction of s-polarised reflectance",
      "Crystal Transmittance (P polarisation)": r"Fraction of p-polarised transmitted",
      "Crystal Transmittance (S polarisation)": r"Fraction of s-polarised transmitted",
        "Crystal Absorbtance (P polarisation)": r"Fraction of p-polarised absorbtance",
        "Crystal Absorbtance (S polarisation)": r"Fraction of s-polarised absorbtance",
                            }

        self.plot_type_cb.activated.connect(self.on_plot_type_cb_activated)
        self.plot_type_cb.addItems( self.plot_types )
        label = QLabel("Choose plot type", self)
        label.setToolTip("Choose the plot type")
        index = self.plot_type_cb.findText(self.settings["Plot type"], Qt.MatchFixedString)
        self.plot_type_cb.setCurrentIndex(index)
        plot_button = QPushButton("Update plot")
        plot_button.clicked.connect(self.refresh)
        plot_button.setToolTip("Update the plot")
        hbox = QHBoxLayout()
        hbox.addWidget(self.plot_type_cb)
        hbox.addWidget(plot_button)
        form.addRow(label, hbox)
        # Add a progress bar
        self.progressbar = QProgressBar(self)
        self.progressbar.setToolTip("Show the progress of any calculations")
        # Append the progress bar to the list of progress bars managed by the notebook
        self.notebook.progressbars_add(self.progressbar)
        self.notebook.progressbars_set_maximum(self.get_total_number_of_frequency_calculations())
        label = QLabel("Calculation progress", self)
        label.setToolTip("Show the progress of any calculations")
        form.addRow(label,self.progressbar)
        # Add the matplotlib figure to the bottom
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding))
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        form.addRow(self.canvas)
        form.addRow(self.toolbar)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        # Create the plot
        debugger.print("Finished:: Plotting tab initialisation")
        return

    def get_total_number_of_frequency_calculations(self):
        """Calculate and return the total number of frequency calculations required.

        This method computes the total number of frequency calculations based on the current settings for minimum frequency, maximum frequency, and frequency increment. It updates the scenario if the frequency range or the increment has changed since the last update.

        Parameters
        ----------
        None

        Returns
        -------
        int
            The total number of frequency calculations required, which is the product of the number of calculations required (as obtained from `get_number_of_calculations_required`) and the new length of frequency range computed.

        """        
        vmin = self.settings["Minimum frequency"]
        vmax = self.settings["Maximum frequency"]
        vinc = self.settings["Frequency increment"]
        self.vs_cm1 = np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc))
        new_length = len(self.vs_cm1)
        if self.frequency_length != new_length or self.vmin != vmin or self.vmax != vmax or self.vinc != vinc :
            self.requestScenarioRefresh()
            self.frequency_length = new_length
            self.vmin = vmin
            self.vmax = vmax
            self.vinc = vinc
        n = self.get_number_of_calculations_required()
        return n*new_length

    def requestRefresh(self):
        """Initiate a refresh request.

        This function flags that a refresh is required
        It doesn't return any value.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: requestRefresh")
        self.refreshRequired = True
        debugger.print("Finished:: requestRefresh")
        return

    def requestScenarioRefresh(self):
        """Request a refresh on all scenarios within a notebook.

        This function triggers a refresh process for the settings tab and all scenarios within the notebook. 

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: requestScenarioRefresh")
        self.notebook.settingsTab.requestRefresh()
        for scenario in self.notebook.scenarios:
            scenario.requestRefresh()
        debugger.print("Finished:: requestScenarioRefresh")
        return

    def on_vinc_changed(self,value):
        """Handle the change in frequency increment and update GUI accordingly.

        This function is triggered when there's a change in the frequency increment value. It adjusts the number of GUI elements based on 
        new values of minimum and maximum frequency and the changed frequency increment. It also updates the settings dict with the 
        new frequency increment value.

        Parameters
        ----------
        value : float
            The new value of the frequency increment.

        Returns
        -------
        None

        Side Effects
        ------------
        - Temporarily blocks signals from the vinc_sb spin box to prevent recursive calls.
        - Updates 'Frequency increment' in the settings dictionary with the new frequency increment value.
        - Requests a refresh of the fitterTab if needed.

        """        
        debugger.print("Start:: on_vinc_changed", value)
        if value <=0 :
            return
        self.vinc_sb.blockSignals(True)
        # Use existing units to calculate the number of samples
        vmin = self.vmin_sb.value()
        vmax = self.vmax_sb.value()
        ngui = int((vmax - vmin) / value) + 1
        # Work out the increment needed for cm-1
        vmin = self.settings["Minimum frequency"]
        vmax = self.settings["Maximum frequency"]
        if ngui <= 2:
            return
        vinc = (vmax - vmin) / (ngui - 1)
        self.settings["Frequency increment"] = vinc
        self.notebook.fitterTab.requestRefresh()
        self.refreshRequired = True
        debugger.print("on_vinc_change ", self.settings["Frequency increment"])
        self.vinc_sb.blockSignals(False)
        debugger.print("Finished:: on_vinc_changed", value)

    def on_vmin_changed(self):
        """Handle the change in the minimum frequency setting.

        This method adjusts the stored minimum frequency in the settings according to the new
        minimum frequency input by the user. It first blocks the signals of the minimum frequency spinbox
        to prevent recursive calls during the adjustment process. It then updates the setting based
        on whether the frequency units are compatible with a predefined list, converting the frequency
        accordingly. After adjusting the settings, it triggers a refresh request for the application's
        fitter tab and marks the need for a general refresh. Finally, it re-enables the signals for
        the minimum frequency spinbox and logs the completion of the process.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - Uses the `blockSignals` method on `self.vmin_sb` to prevent signal-slot recursion.
        - Adjusts settings based on a unit conversion utility function, `convert_frequency_units`.
        - The effect of this method extends beyond just the internal state changes; it influences the UI and possibly other components' states through the requested refresh.

        """        
        debugger.print("Start:: on_vmin_changed")
        self.vmin_sb.blockSignals(True)
        vmin = self.vmin_sb.value()
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            self.settings["Minimum frequency"] = convert_frequency_units(vmin,self.settings["Frequency unit"],"wavenumber")
        else:
            self.settings["Maximum frequency"] = convert_frequency_units(vmin,self.settings["Frequency unit"],"wavenumber")
        debugger.print("on_vmin_changed setting vmin to", self.settings["Minimum frequency"])
        self.notebook.fitterTab.requestRefresh()
        self.refreshRequired = True
        self.vmin_sb.blockSignals(False)
        debugger.print("Finished:: on_vmin_changed")

    def on_vmax_changed(self):
        """Handle the event when the maximum frequency setting is changed.

        This method is triggered whenever there is a change in the maximum frequency setting (`vmax`). It blocks signal emission from the `vmax_sb` spinner box (presumably, a GUI element for setting `vmax`), reads the current minimum (`vmin`) and maximum (`vmax`) frequency settings, and updates the corresponding setting based on the specified frequency unit. If the frequency unit denotes a frequency, `vmax` is converted and stored as the 'Maximum frequency' in the settings dictionary in the 'wavenumber' unit. Otherwise, it updates the 'Minimum frequency' with the converted value. This method also triggers any necessary refresh operations in the user interface and finally unblocks signal emissions from the `vmax_sb`.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - This method assumes the existence of `vmin_sb` and `vmax_sb` attributes, which should be spinner box GUI elements (or similar) for setting minimum and maximum frequencies, respectively.
        - The `settings` dictionary must have a 'Frequency unit' key, and possibly 'Maximum frequency' and 'Minimum frequency' keys which are updated based on the condition.

        """        
        debugger.print("Start:: on_vmax_changed")
        self.vmax_sb.blockSignals(True)
        vmax = self.vmax_sb.value()
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            self.settings["Maximum frequency"] = convert_frequency_units(vmax,self.settings["Frequency unit"],"wavenumber")
        else:
            self.settings["Minimum frequency"] = convert_frequency_units(vmax,self.settings["Frequency unit"],"wavenumber")
        debugger.print("on_vmax_changed setting vmax to ", self.settings["Maximum frequency"])
        self.notebook.fitterTab.requestRefresh()
        self.refreshRequired = True
        self.vmax_sb.blockSignals(False)
        debugger.print("Finished:: on_vmax_changed")

    def refresh(self,force=False):
        """Refresh the current state based on the changes in settings, forces refresh if needed.

        This function updates the frequency settings (minimum frequency, maximum frequency, and frequency increment) based on the current settings. It ensures that these settings are within logical limits. The function updates GUI elements (spin boxes and combo boxes) with these new calculations. Additionally, it handles the conversion of frequency units, updates GUI tooltips based on the frequency unit, and builds or refreshes the visualization plot based on the refreshed data. If the 'force' flag is set or if a refresh is deemed necessary due to significant changes in the settings, it proceeds with the refresh operation, else it skips the refresh to optimize performance.

        Parameters
        ----------
        force : bool, optional
            A boolean flag indicating whether to forcefully execute a refresh regardless of whether it was deemed necessary based on internal conditions. The default is False.

        Returns
        -------
        None

        Notes
        -----
        - The actual refresh operation involves several steps:
            - Checking if a refresh is required based on the current settings and the 'force' parameter.
            - Temporarily blocking signals from all child widgets to prevent unintended side effects during settings updates.
            - Verifying and updating frequency settings to maintain logical constraints.
            - Updating GUI components with new settings values and tooltips based on the current frequency unit.
            - Processing and plotting data based on the updated settings.
            - Re-enabling signals for child widgets after modifications are complete.
        - The processEvents call is used to ensure the UI remains responsive during long operations.

        """        
        debugger.print("Start:: refresh", force)
        if not self.refreshRequired and not force:
            self.plot()
            debugger.print("Finished:: refreshing widget not required")
            return
        #
        # Block signals during refresh
        #
        self.greyed_out()
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Now refresh values
        if self.settings["Maximum frequency"] < self.settings["Minimum frequency"]:
            self.settings["Maximum frequency"] = self.settings["Minimum frequency"]+1
        if self.settings["Frequency increment"] > self.settings["Maximum frequency"] - self.settings["Minimum frequency"]:
            self.settings["Frequency increment"] = (self.settings["Maximum frequency"] - self.settings["Minimum frequency"])/2
        # The GUI uses frequency and wavelength, convert the values in settings from cm-1
        # Work out the increment needed for cm-1
        vmin = self.settings["Minimum frequency"]
        vmax = self.settings["Maximum frequency"]
        vinc = self.settings["Frequency increment"]
        # Protect the code from over-exuberant choice of parameters
        if (vmax - vmin)/vinc > 90000:
            vinc = (vmax - vmin) / 90000
            print("Warning - the number data points in a plot has been limited to 9000")
            print("          this happens if a 0 wavelength or frequency is entered in the GUI")
            self.settings["Frequency increment"] = vinc
        ncm1 = len(np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc)))
        vmin = convert_frequency_units(self.settings["Minimum frequency"],"wavenumber",self.settings["Frequency unit"])
        vmax = convert_frequency_units(self.settings["Maximum frequency"],"wavenumber",self.settings["Frequency unit"])
        # Reorder the GUI values
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        # Set the GUI values
        vinc = (vmax - vmin) / (ncm1-1)
        self.vmin_sb.setValue(vmin)
        self.vmax_sb.setValue(vmax)
        self.vinc_sb.setValue(vinc)
        # Update the tool tips as the frequency unit could be a wavelength
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            self.vmin_sb.setToolTip("Set the minimum frequency to be considered)")
            self.vmax_sb.setToolTip("Set the maximum frequency to be considered)")
            self.vinc_sb.setToolTip("Choose an increment for the frequency when plotting")
            self.funits_cb.setToolTip("Set the frequency unit")
            self.frequency_form_label.setToolTip("Choose minimum, maximum and increment for frequency")
            self.frequency_form_label.setText("Frequency min, max and increment")
            self.frequency_form_label.setToolTip("Choose minimum, maximum and increment for frequency")
        else:
            self.vmin_sb.setToolTip("Set the minimum wavelength to be considered)")
            self.vmax_sb.setToolTip("Set the maximum wavelength to be considered)")
            self.vinc_sb.setToolTip("Choose an increment for the wavelength when plotting")
            self.funits_cb.setToolTip("Set the wavelength unit")
            self.frequency_form_label.setText("Wavelength min, max and increment")
            self.frequency_form_label.setToolTip("Choose minimum, maximum and increment for wavelength")
        index = possible_frequency_units.index(self.settings["Frequency unit"])
        self.funits_cb.setCurrentIndex(index)
        index = self.plot_type_cb.findText(self.settings["Plot type"], Qt.MatchFixedString)
        self.plot_type_cb.setCurrentIndex(index)
        try:
            self.molar_cb_current_index = self.molar_definitions.index(self.settings["Molar definition"])
        except Exception:
            self.molar_cb_current_index = 0
            self.settings["Molar definition"] = self.molar_definitions[self.molar_cb_current_index]
        self.molar_cb.setCurrentIndex(self.molar_cb_current_index)
        self.natoms_sb.setValue(self.settings["Number of atoms"])
        # Refresh the widgets that depend on the reader
        self.reader = self.notebook.reader
        if self.reader is not None:
            self.set_concentrations()
        # Reset the progress bar
        self.notebook.progressbars_set_maximum(self.get_total_number_of_frequency_calculations())
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print("calling plot from refresh")
        self.plot()
        self.refreshRequired = False
        debugger.print("Finished:: refresh", force)
        return

    def on_natoms_changed(self, value):
        """Handle the change in the number of atoms.

        This method is called when the number of atoms changes. It updates the relevant setting in the instance, recalculates the concentration based on the new number of atoms, flags that a refresh is required, and then triggers the refresh process.

        Parameters
        ----------
        value : int or float
            The new number of atoms. This value is used to update settings and recalculate concentrations.

        Notes
        -----
        - `self.reader.volume` and `self.reader.nions` are expected to be available and contain the volume of the container and the number of ions, respectively.
        - The method refreshes `self.notebook.fitterTab` object with a `requestRefresh` method to initiate the refresh process.

        """        
        debugger.print("Start:: on_natoms_changed", value)
        self.settings["Number of atoms"] = value
        debugger.print("on natoms changed ", self.settings["Number of atoms"])
        self.settings["concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 * self.settings["Number of atoms"] / self.reader.nions)
        debugger.print("The concentration has been set", self.settings["Molar definition"], self.settings["concentration"])
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refresh()
        debugger.print("Finished:: on_natoms_changed", value)

    def on_plot_type_cb_activated(self, index):
        """Handle plot type change from a combo box.

        This method gets triggered when the plot type in a combo box is changed.
        It updates the plot type in the settings, marks the current state as needing a refresh,
        requests the fitter tab to refresh, and finally refreshes the current view.

        Parameters
        ----------
        index : int
            The index of the selected item in the combo box.

        Returns
        -------
        None

        """        
        debugger.print("Start:: on_plot_type_cb_activated", index)
        self.settings["Plot type"] = self.plot_type_cb.currentText()
        debugger.print("Changed plot type to ", self.settings["Plot type"])
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refresh()
        debugger.print("Finished:: on_plot_type_cb_activated", index)

    def on_funits_cb_activated(self, index):
        """Handle the activation of a frequency unit combo box item.

        This method updates the active frequency unit in the settings, triggers necessary refresh processes and logs the change.

        Parameters
        ----------
        index : int
            The index of the activated item in the frequency units combo box.

        Returns
        -------
        None

        Notes
        -----
        This method is part of a GUI application where `self` refers to an instance of the application or a relevant widget. It manages updates to the application's settings and GUI components based on the user's selection of a frequency unit from a combo box.

        - `self.settings` is a dictionary where application settings are stored.
        - `self.refreshRequired` is a boolean flag used to indicate whether a refresh of certain GUI components is necessary.
        - `self.notebook` appears to be a widget container (like a tab widget), with `fitterTab` being one of its child tabs.
        - `self.vmin_sb` is another component (likely a spin box or similar input widget) that may need its signals blocked/unblocked during the process to avoid unwanted signal emission.
        - The `debugger.print` calls are used for logging and are not a standard Python function; they imply the existence of a custom logging or debugging utility named `debugger`.

        The actual refreshing of the GUI and handling of signal blocking is done within other methods not shown here, such as `self.refresh()` and `self.notebook.fitterTab.requestRefresh()`.

        Raises
        ------
        This function does not explicitly raise any exceptions but depends on the proper functioning of the methods it calls and the state of `self` and its attributes.

        """        
        debugger.print("Start:: on_funits_cb_activated", index)
        self.settings["Frequency unit"] = possible_frequency_units[index]
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refreshRequired = True
        self.vmin_sb.blockSignals(False)
        self.refresh()
        debugger.print("Frequency unit changed to ", self.settings["Frequency unit"])
        debugger.print("Finished:: on_funits_cb_activated", index)

    def on_molar_cb_activated(self, index):
        """Handle the activation of the molar combobox option.

        This method is tied to a GUI event where a selection from a molar combobox triggers various updates in the application state, including setting concentrations, refreshing UI elements, and potentially triggering further calculations or refreshes as needed.

        Parameters
        ----------
        index : int
            The index of the selected item in the molar combobox. This index corresponds to a specific molar definition and is used to update the application settings and state accordingly.

        Returns
        -------
        None

        See Also
        --------
        set_concentrations : A method to update the concentrations based on the selected molar definition.
        refresh : A method to refresh the UI elements.

        """        
        debugger.print("Start:: on_molar_cb_activated", index)
        self.molar_cb_current_index = index
        self.settings["Molar definition"] = self.molar_definitions[index]
        self.set_concentrations()
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refresh()
        debugger.print("The concentration has been set", self.settings["Molar definition"], self.settings["concentration"])
        debugger.print("Finished:: on_molar_cb_activated", index)
        return

    def set_concentrations(self):
        """Set the concentration based on the molar definition in settings.

        This method updates the concentration value in the settings dictionary based on the 'Molar definition' key. It calculates concentration values differently based on whether the molar definition is set to 'Molecules', 'Unit cells', or 'Atoms'. It also enables or disables the `natoms_sb` spin box accordingly.
        """        
        debugger.print("Start:: set_concentration")
        if self.settings["Molar definition"] == "Molecules":
            self.settings["concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 * self.settings["Number of atoms"] / self.reader.nions)
            self.natoms_sb.setEnabled(True)
        elif self.settings["Molar definition"] == "Unit cells":
            self.settings["concentration"]      = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
            self.settings["cell concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
            self.natoms_sb.setEnabled(False)
        elif self.settings["Molar definition"] == "Atoms":
            self.settings["concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 / self.reader.nions)
            self.natoms_sb.setEnabled(False)
        debugger.print("Finished:: set_concentration")
        return

    def writeSpreadsheet(self):
        """Update and write the results of powder and crystal scenarios to a spreadsheet.

        This function navigates through each scenario defined in the `notebook` attribute, extracts relevant data such as absorption coefficients, permittivities, and reflectances, and writes these along with scenario settings to a spreadsheet. It handles different types of scenarios (Powder or Crystal) and makes use of the spreadsheet object's methods for selecting worksheets, writing rows, and dealing with data transformation. The function also handles the calculation of molar absorption coefficients with unit conversion when necessary.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start::writeSpreadsheet")
        if self.notebook.spreadsheet is None:
            debugger.print("Finished::writeSpreadsheet spreadsheet is None")
            return
        # make sure the plottingTab is up to date
        self.refresh()
        # Handle powder plots
        molarAbsorptionCoefficients = []
        absorptionCoefficients      = []
        realPermittivities          = []
        imagPermittivities          = []
        sp_atrs                     = []
        R_ps                        = []
        R_ss                        = []
        T_ps                        = []
        T_ss                        = []
        A_ps                        = []
        A_ss                        = []
        powder_legends              = []
        crystal_legends             = []
        # Deal with Scenarios 
        sp = self.notebook.spreadsheet
        sp.selectWorkSheet("Scenarios")
        sp.delete()
        sp.writeNextRow(["A list of the scenarios used:"],col=1)
        for index,scenario in enumerate(self.notebook.scenarios):
            if scenario.scenarioType == "Powder":
                direction = scenario.direction
                depolarisation = scenario.depolarisation
                sp.writeNextRow([""],col=1)
                sp.writeNextRow(["Scenario "+str(index)],col=1,check=1)
                settings = scenario.settings
                for key in sorted(settings,key=str.lower):
                    sp.writeNextRow([key, settings[key]],col=1,check=1)
                sp.writeNextRow(["Normalised unique direction"]+direction.tolist(), col=1,check=1)
                sp.writeNextRow(["Depolarisation matrix"], col=1,check=1)
                sp.writeNextRow(depolarisation[0].tolist(), col=2, check=1)
                sp.writeNextRow(depolarisation[1].tolist(), col=2, check=1)
                sp.writeNextRow(depolarisation[2].tolist(), col=2, check=1)
                molarAbsorptionCoefficients.append( scenario.get_result(self.vs_cm1,self.plot_types[0] ) )
                absorptionCoefficients.append( scenario.get_result(self.vs_cm1,self.plot_types[1] ) )
                realPermittivities.append( scenario.get_result(self.vs_cm1,self.plot_types[2] ) )
                imagPermittivities.append( scenario.get_result(self.vs_cm1,self.plot_types[3] ) )
                sp_atrs.append( scenario.get_result(self.vs_cm1,self.plot_types[4] ) )
                powder_legends.append(scenario.settings["Legend"])
            else:
                sp.writeNextRow([""],col=1)
                sp.writeNextRow(["Scenario "+str(index)],col=1,check=1)
                settings = scenario.settings
                for key in sorted(settings,key=str.lower):
                    sp.writeNextRow([key, settings[key]],col=1,check=1)
                dielectricLayerIndex = scenario.getDielectricLayerIndex()
                if dielectricLayerIndex is not None and scenario.layers[dielectricLayerIndex].isTensor():
                    sp.writeNextRow("Dielectric layer laboratory frame:")
                    sp.writeNextRow(scenario.layers[dielectricLayerIndex].labframe[0].tolist(), col=2, check=1)
                    sp.writeNextRow(scenario.layers[dielectricLayerIndex].labframe[1].tolist(), col=2, check=1)
                    sp.writeNextRow(scenario.layers[dielectricLayerIndex].labframe[2].tolist(), col=2, check=1)
                # Store the reflectance and transmittance
                R_ps.append( scenario.get_result(self.vs_cm1,self.plot_types[5] ) )
                R_ss.append( scenario.get_result(self.vs_cm1,self.plot_types[6] ) )
                T_ps.append( scenario.get_result(self.vs_cm1,self.plot_types[7] ) )
                T_ss.append( scenario.get_result(self.vs_cm1,self.plot_types[8] ) )
                A_ps.append( scenario.get_result(self.vs_cm1,self.plot_types[9] ) )
                A_ss.append( scenario.get_result(self.vs_cm1,self.plot_types[10] ) )
                crystal_legends.append(scenario.settings["Legend"])
        # Single crystal Permittivity
        dielecv = self.notebook.settingsTab.getCrystalPermittivity(self.vs_cm1)
        # Powder results
        # Work out what molar units we are using
        if len(molarAbsorptionCoefficients) > 0:
            if self.settings["Molar definition"] == "Molecules":
                sheet_name = "Powder Molar Absorption (mols)"
            elif self.settings["Molar definition"] == "Unit cells":
                sheet_name = "Powder Molar Absorption (cells)"
            elif self.settings["Molar definition"] == "Atoms":
                sheet_name = "Powder Molar Absorption (atoms)"
            # Always write out the moles of cell
            self.write_powder_results(sp, "Powder Molar Absorption (cells)", self.vs_cm1, powder_legends, molarAbsorptionCoefficients)
            if self.settings["Molar definition"] != "Unit cells":
                # If some other molar definition has been used then write that out too
                molarAbsorptionCoefficients_mols = []
                molar_scaling = self.settings["cell concentration"]/self.settings["concentration"]
                for absorption in molarAbsorptionCoefficients:
                    molarAbsorptionCoefficients_mols.append(molar_scaling * np.array(absorption))
                self.write_powder_results(sp, sheet_name,                      self.vs_cm1, powder_legends, molarAbsorptionCoefficients_mols)
            # end if
            self.write_powder_results(sp, "Powder Absorption",             self.vs_cm1, powder_legends, absorptionCoefficients)
            self.write_powder_results(sp, "Powder Real Permittivity",      self.vs_cm1, powder_legends, realPermittivities)
            self.write_powder_results(sp, "Powder Imaginary Permittivity", self.vs_cm1, powder_legends, imagPermittivities)
            self.write_powder_results(sp, "Powder ATR Reflectance",        self.vs_cm1, powder_legends, sp_atrs)
        # Single Crystal results
        if len(R_ps) > 0:
            self.write_crystal_results(sp, "Crystal R_p", self.vs_cm1, crystal_legends, R_ps)
            self.write_crystal_results(sp, "Crystal R_s", self.vs_cm1, crystal_legends, R_ss)
            self.write_crystal_results(sp, "Crystal T_p", self.vs_cm1, crystal_legends, T_ps)
            self.write_crystal_results(sp, "Crystal T_s", self.vs_cm1, crystal_legends, T_ss)
            self.write_crystal_results(sp, "Crystal A_p", self.vs_cm1, crystal_legends, A_ps)
            self.write_crystal_results(sp, "Crystal A_s", self.vs_cm1, crystal_legends, A_ss)

        if len(dielecv) > 0:
            self.write_eps_results(sp, self.vs_cm1, dielecv)
        debugger.print("Finished::writeSpreadsheet")
        return

    def write_eps_results(self, sp, vs, dielecv):
        """Write real and imaginary parts of crystal permittivity to a spreadsheet.

        Parameters
        ----------
        sp : object
            An instance of a spreadsheet processing class with methods to manipulate data in a spreadsheet.
        vs : list
            A list containing frequency values.
        dielecv : numpy.ndarray
            A complex numpy array where the real parts represent the real permittivities and the imaginary parts represent the imaginary permittivities of a crystal. The array should have a shape of (N,3,3) where N is the number of frequency values, and the 3x3 inner arrays represent the permittivity tensor for each frequency.

        Returns
        -------
        None

        Notes
        -----
        The function does two main tasks:
        1. Selects the 'Real Crystal Permittivity' sheet, deletes its current content if any, and writes the real parts of the permittivity tensor for each frequency along with the frequency values themselves.
        2. Selects the 'Imag Crystal Permittivity' sheet, deletes its current content if any, and writes the imaginary parts of the permittivity tensor for each frequency along with the frequency values themselves.

        Both sections write data in the format: frequencies (cm-1), xx, yy, zz, xy, xz, yz, where xx, yy, zz, xy, xz, and yz are components of the permittivity tensor.
        The output data starts from the second column (index 1), and for each row written, a 'check' flag is set to 1.

        """        
        debugger.print("Start:: write_eps_results length vs",len(vs))
        sp.selectWorkSheet("Real Crystal Permittivity")
        sp.delete()
        headers = ["frequencies (cm-1)", "xx", "yy", "zz", "xy", "xz", "yz" ]
        sp.writeNextRow(headers,row=0, col=1)
        for v,eps in zip(vs,dielecv):
            eps_xx_r = np.real(eps[0][0])
            eps_yy_r = np.real(eps[1][1])
            eps_zz_r = np.real(eps[2][2])
            eps_xy_r = np.real(eps[0][1])
            eps_xz_r = np.real(eps[0][2])
            eps_yz_r = np.real(eps[1][2])
            output = [v, eps_xx_r, eps_yy_r, eps_zz_r, eps_xy_r, eps_xz_r, eps_yz_r ]
            sp.writeNextRow(output, col=1,check=1)
        sp.selectWorkSheet("Imag Crystal Permittivity")
        sp.delete()
        sp.writeNextRow(headers,row=0, col=1)
        for v,eps in zip(vs,dielecv):
            eps_xx_i = np.imag(eps[0][0])
            eps_yy_i = np.imag(eps[1][1])
            eps_zz_i = np.imag(eps[2][2])
            eps_xy_i = np.imag(eps[0][1])
            eps_xz_i = np.imag(eps[0][2])
            eps_yz_i = np.imag(eps[1][2])
            output = [v, eps_xx_i, eps_yy_i, eps_zz_i, eps_xy_i, eps_xz_i, eps_yz_i ]
            sp.writeNextRow(output, col=1,check=1)
        debugger.print("Finished:: write_eps_results length vs",len(vs))
        return

    def write_crystal_results(self, sp, name, vs, legends, yss):
        """Write single crystal results to a spread sheet.

        Parameters
        ----------
        sp : object
            The spreadsheet object.
        name : str
            The worksheet name used for writing.
        vs : np.array
            An array of the frequencies.
        legends : list
            The heading names for the yss.
        yss : list of np.arrays
            A list of numpy arrays of the reflections and transmittance.

        Returns
        -------
        None

        """
        debugger.print("Start:: write_crystal_results")
        debugger.print("write_crystal_results name",name)
        debugger.print("write_crystal_results legends",legends)
        debugger.print("write_crystal_results length vs",len(vs))
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ["frequencies (cm-1)"]
        headers.extend(legends)
        sp.writeNextRow(headers,row=0, col=1)
        for iv,v in enumerate(vs):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)
        debugger.print("Finished:: write_crystal_results")
        return

    def write_powder_results(self, sp, name, vs, legends, yss):
        """Write the powder simulation results to a worksheet.

        Parameters
        ----------
        sp : object
            An object with methods for manipulating a worksheet.
        name : str
            The name of the worksheet to select or create.
        vs : list
            A list of frequency values (in cm-1).
        legends : list of str
            The legends to be used as headers, alongside 'frequencies (cm-1)'.
        yss : list of lists
            A list where each element is a list of intensities corresponding to each frequency in `vs`.

        Returns
        -------
        None

        """        
        debugger.print("Start:: write powder results")
        debugger.print("write_powder_results name",name)
        debugger.print("write_powder_results legends",legends)
        debugger.print("write_powder_results length vs",len(vs))
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ["frequencies (cm-1)"]
        #for isc,ys in enumerate(yss):
        #    headers.append('Scenario'+str(isc))
        headers.extend(legends)
        sp.writeNextRow(headers,row=0, col=1)
        for iv,v in enumerate(vs):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)
        debugger.print("Finished:: write powder results")
        return

    def plot(self):
        # import matplotlib.pyplot as pl
        # mp.use('qtagg')
        """Plot the results based on specified settings and scenarios.

        This function generates a plot for the given scenarios under the analysis settings
        defined in the notebook object's tabs. It handles various conditions to properly format
        the plot, including checking if necessary parameters and objects are set, converting
        frequency units, and adjusting data for display based on plot type. The actual plotting
        is done using matplotlib, and various attributes like plot legends, titles, and axes labels
        are dynamically set based on the configuration.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - self.notebook.mainTab.settings: Dictionary containing program settings.
        - self.notebook.mainTab.getFullFileName(): Method that returns the currently selected filename.
        - self.notebook.mainTab.reader: Object used for reading data files.
        - self.notebook.settingsTab.CrystalPermittivityObject: Object containing crystal permittivity settings.
        - self.settings: Dictionary containing plot-related settings such as frequency range and plot type.
        - self.notebook.scenarios: List of scenarios to be plotted.

        """        
        debugger.print("Start:: plot")
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings["Program"]
        filename = self.notebook.mainTab.getFullFileName()
        reader = self.notebook.mainTab.reader
        if reader is None:
            debugger.print("Finished:: plot aborting because reader is NONE")
            return
        if program == "":
            debugger.print("Finished:: plot aborting because program is not set")
            return
        if filename == "":
            debugger.print("Finished:: plot aborting because filename is not set")
            return
        if self.notebook.settingsTab.CrystalPermittivityObject is None:
            debugger.print("Finished:: plot aborting because settingTab.CrystalPermittivityObject is not set")
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        vmin = self.settings["Minimum frequency"]
        vmax = self.settings["Maximum frequency"]
        vinc = self.settings["Frequency increment"]
        self.vs_cm1 = np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc))
        self.subplot = None
        self.figure.clf()
        x = np.array(self.vs_cm1)
        removeFirstElement = False
        reverseElements = False
        if x[0] <= 0:
            # We can't have a 0 frequency converted to a wavelength
            removeFirstElement = True
            x = x[1:]
        if isThisAFrequencyUnit(self.settings["Frequency unit"]):
            x = x[::-1]
            reverseElements = True
        xlabel = self.settings["Frequency unit"]
        if self.settings["Frequency unit"] == "wavenumber":
            xlabel = r"Frequency $\mathdefault{(cm^{-1})}$"
        x = convert_frequency_units( x, "wavenumber", self.settings["Frequency unit"] )
        self.subplot = self.figure.add_subplot(111)
        self.notebook.progressbars_set_maximum(self.get_total_number_of_frequency_calculations())
        self.legends = []
        plots = 0
        for scenario in self.notebook.scenarios:
            legend = scenario.settings["Legend"]
            self.legends.append(legend)
            y = scenario.get_result(self.vs_cm1,self.settings["Plot type"])
            if y is not None and len(y) > 0:
                y = np.array(y)
                if removeFirstElement:
                    y = y[1:]
                if reverseElements:
                     y = y[::-1]
                if self.settings["Plot type"] == "Powder Molar Absorption":
                    y = y * self.settings["cell concentration"]/self.settings["concentration"]
                plots += 1
                line, = self.subplot.plot(x,y,lw=2, label=legend )
        if plots > 0:
            self.subplot.set_xlabel(xlabel)
            self.subplot.set_ylabel(self.plot_ylabels[self.settings["Plot type"]])
            self.subplot.legend(loc="best")
            self.subplot.set_title(self.settings["Plot type"])
            #self.subplot.set_autoscaley_on(False)
            #self.subplot.set_ylim([0,1])
            self.canvas.draw_idle()
        QApplication.restoreOverrideCursor()
        debugger.print("Finished:: plot")

    def get_number_of_calculations_required(self):
        """Return the total number of spectra that need to be calculated.

        Only spectra that need a refresh are included.

        Parameters
        ----------
        None

        Returns
        -------
        int
            The total number of spectra that require calculation.

        """
        debugger.print("Start:: get_number_of_calculations_required")
        n = 0
        for scenario in self.notebook.scenarios:
            n += scenario.getNoCalculationsRequired()
        debugger.print("get_number_of_calculations_required",n)
        return n

    def greyed_out(self):
        """Handle items that should be greyed out if they are not needed.

        Parameters
        ----------
        None

        Returns
        -------
        int

        """
        debugger.print("Start:: greyed_out")
        powder_scenarios_present = False
        crystal_scenarios_present = False
        for scenario in self.notebook.scenarios:
            if scenario.scenarioType == "Powder":
                powder_scenarios_present = True
            else:
                crystal_scenarios_present = True
        # end of for loop
        # 
        # Disable any plot types that are not needed
        #
        self.plot_type_cb.model().item(0).setEnabled(True)
        self.plot_type_cb.model().item(1).setEnabled(True)
        self.plot_type_cb.model().item(2).setEnabled(True)
        self.plot_type_cb.model().item(3).setEnabled(True)
        self.plot_type_cb.model().item(4).setEnabled(True)
        self.plot_type_cb.model().item(5).setEnabled(True)
        self.plot_type_cb.model().item(6).setEnabled(True)
        self.plot_type_cb.model().item(7).setEnabled(True)
        self.plot_type_cb.model().item(8).setEnabled(True)
        self.plot_type_cb.model().item(9).setEnabled(True)
        self.plot_type_cb.model().item(10).setEnabled(True)
        index = self.plot_type_cb.findText(self.settings["Plot type"], Qt.MatchFixedString)
        if not powder_scenarios_present:
            self.plot_type_cb.model().item(0).setEnabled(False)
            self.plot_type_cb.model().item(1).setEnabled(False)
            self.plot_type_cb.model().item(2).setEnabled(False)
            self.plot_type_cb.model().item(3).setEnabled(False)
            self.plot_type_cb.model().item(4).setEnabled(False)
            if index < 5:
                self.plot_type_cb.setCurrentIndex(5)
                self.settings["Plot type"] = self.plot_type_cb.currentText()
        if not crystal_scenarios_present:
            self.plot_type_cb.model().item(5).setEnabled(False)
            self.plot_type_cb.model().item(6).setEnabled(False)
            self.plot_type_cb.model().item(7).setEnabled(False)
            self.plot_type_cb.model().item(8).setEnabled(False)
            self.plot_type_cb.model().item(9).setEnabled(False)
            self.plot_type_cb.model().item(10).setEnabled(False)
            if index >= 5:
                self.plot_type_cb.setCurrentIndex(0)
                self.settings["Plot type"] = self.plot_type_cb.currentText()
        debugger.print("Finished:: greyed_out")
