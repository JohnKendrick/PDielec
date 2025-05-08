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
"""AnalysisTab Module."""

# Import plotting requirements
import copy

import matplotlib
import matplotlib.figure
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.ticker import MaxNLocator
from qtpy.QtCore import QCoreApplication, Qt
from qtpy.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSizePolicy,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from PDielec import Calculator
from PDielec.Constants import covalent_radii
from PDielec.GUI.SettingsTab import FixedQTableWidget
from PDielec.Utilities import Debug


class AnalysisTab(QWidget):
    """A widget class for analyzing vibrational modes, molecular composition, and bonding configurations within a molecular dataset.

    This class inherits from QWidget
    Key features include plotting vibrational mode decompositions

    Parameters
    ----------
    parent : QObject
        The parent widget or object, typically the main application or main window in which this widget will be embedded.
    debug : bool, optional
        Indicates if debugging is enabled for this widget. The default is False.

    Attributes
    ----------
    subplot : matplotlib subplot
        A subplot attribute for plotting the analysis results. 
    settings : dict
        A dictionary for storing various settings related to analysis, like radii, frequencies, and plotting configurations.
    settings['Radii'] : list 
        Holds specific radii settings if provided; otherwise, set to None.
    settings['Minimum frequency'] : int
        The minimum frequency boundary for analysis, initialized to -1.
    settings['Maximum frequency'] : int
        The maximum frequency boundary for analysis, initialized to 400.
    settings['title'] : str
        The title for plots generated within the widget. Defaults to 'Analysis'.
    settings['Covalent radius scaling'] : float
        A scaling factor for adjusting covalent radii in bonding calculations. Defaults to 1.1.
    settings['Bonding tolerance'] : float
        The tolerance level used in bonding calculations, determining what is considered a bond. Defaults to 0.1.
    settings['Bar width'] : float
        Determines the width of bars in bar plots. Defaults to 0.5.
    refreshRequired : bool
        A flag indicating whether the widget needs to be refreshed to update the visuals or calculations. Defaults to True.
    plot_types : list
        A list of available plot types the user can select from.
    plot_type_index : int
        The index of the currently selected plot type in `plot_types` list.
    number_of_molecules : int
        The number of molecules identified during the last analysis. Defaults to 0.
    frequency_units : str or None
        String representing the units used for frequency measurements. Initially set to None.
    cell_of_molecules : object or None
        An object containing information about the cell of molecules, if analysis has been run. Initially set to None.
    frequencies_cm1 : list
        A list holding frequency values in cm-1 units.
    mode_energies : list
        Stores the calculated energies for different vibrational modes.
    element_radii : dict
        A dictionary mapping elements to their corresponding covalent radii used in the analysis.
    species : list
        A list of species involved in the current analysis.
    notebook : object
        A reference to the parent notebook-like structure which hosts this widget.
    reader : object
        An object responsible for reading and interpreting molecular data files.

    Methods
    -------
    on_element_radii_tw_itemClicked(self, item)
        Handler for click events on table widget items related to element radii.
    on_element_radii_tw_itemChanged(self, item)
        Handler for changes in the table widget items related to element radii.
    set_radii_tw(self)
        Updates the radii table widget based on current settings or data.
    setCovalentRadius(self, element, radius)
        Sets the covalent radius for a given element and updates the analysis.
    writeSpreadsheet(self)
        Writes the analysis results to a spreadsheet hosted by the parent structure.
    on_width_changed(self, value)
        Handler for changes in the bar width setting.
    on_scale_changed(self, value)
        Handler for changes in the covalent radius scaling setting.
    on_tolerance_changed(self, value)
        Handler for changes in the bonding tolerance setting.
    on_title_changed(self, text)
        Handler for changes in the plot title.
    on_vmin_changed(self)
        Handler for changes in the minimum frequency setting.
    on_vmax_changed(self)
        Handler for changes in the maximum frequency setting.
    requestRefresh(self)
        Requests a refresh of the analysis and visualization.
    refresh(self, force=False)
        Refreshes the widget based on current settings, optionally forcing a refresh.
    on_plottype_cb_changed(self, index)
        Handler for changes in the selected plot type.
    calculate(self)
        Performs calculations based on current settings, preparing data for visualization.
    plot(self)
        Invokes the plotting method corresponding to the current plot type selection.
    plot_molecular(self)
        Generates a plot based on molecular composition of vibrational energy.
    plot_internal_external(self)
        Generates a plot based on the decomposition of vibrational energy into internal and external movements.

    """

    def __init__(self, parent, debug=False ):
        """Initialise the analysis tab with configurable settings and plotting capabilities.

        Parameters
        ----------
        parent : QObject
            The parent widget or object, typically the main application or main window in which this widget will be embedded.
        debug : bool, optional
            Indicates if debugging is enabled for this widget. The default is False.

        """        
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,"AnalysisTab")
        self.settings = {}
        self.subplot = None
        self.setWindowTitle("Analysis")
        self.settings["Radii"] = None
        self.settings["Minimum frequency"] = -1
        self.settings["Maximum frequency"] = 400
        self.settings["title"] = "Analysis"
        self.settings["Covalent radius scaling"] = 1.1
        self.settings["Bonding tolerance"] = 0.1
        self.settings["Bar width"] = 0.5
        self.refreshRequired = True
        self.plot_types = ["Internal vs External","Molecular Composition"]
        self.plot_type_index = 0
        self.number_of_molecules = 0
        self.frequency_units = None
        self.cell_of_molecules = None
        self.frequencies_cm1 = []
        self.mode_energies = []
        self.element_radii = covalent_radii
        self.species = []
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Create last tab - AnalysisTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The minimum frequency
        #
        self.vmin_sb = QDoubleSpinBox(self)
        self.vmin_sb.setRange(-100,9000)
        self.vmin_sb.setValue(self.settings["Minimum frequency"])
        self.vmin_sb.setToolTip("Set the minimum frequency to be considered)")
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        label = QLabel("Minimum frequency:", self)
        label.setToolTip("Set the minimum frequency to be considered)")
        form.addRow(label, self.vmin_sb)
        #
        # The maximum frequency
        #
        self.vmax_sb = QDoubleSpinBox(self)
        self.vmax_sb.setRange(0,9000)
        self.vmax_sb.setValue(self.settings["Maximum frequency"])
        self.vmax_sb.setToolTip("Set the maximum frequency to be considered)")
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        label = QLabel("Maximum frequency:", self)
        label.setToolTip("Set the maximum frequency to be considered)")
        form.addRow(label, self.vmax_sb)
        #
        # The bonding tolerance and scaling of radii frequency
        #
        hbox = QHBoxLayout()
        self.scale_sp = QDoubleSpinBox(self)
        self.scale_sp.setRange(0.01,10.0)
        self.scale_sp.setSingleStep(0.01)
        self.scale_sp.setDecimals(2)
        self.scale_sp.setValue(self.settings["Covalent radius scaling"])
        self.scale_sp.setToolTip("Scale the covalent radii to determine bonding")
        self.scale_sp.valueChanged.connect(self.on_scale_changed)
        hbox.addWidget(self.scale_sp)
        self.tolerance_sp = QDoubleSpinBox(self)
        self.tolerance_sp.setRange(0.01,2.0)
        self.tolerance_sp.setSingleStep(0.01)
        self.tolerance_sp.setDecimals(2)
        self.tolerance_sp.setValue(self.settings["Bonding tolerance"])
        self.tolerance_sp.setToolTip("Tolerance for bonding is determined from scale*(radi+radj)+toler")
        self.tolerance_sp.valueChanged.connect(self.on_tolerance_changed)
        hbox.addWidget(self.tolerance_sp)
        label = QLabel("Bonding scale and tolerance", self)
        label.setToolTip("Bonding is determined from scale*(radi+radj)+toler")
        form.addRow(label, hbox)
        # Add a table of covalent radii
        self.element_radii_tw = FixedQTableWidget(parent=self)
        self.element_radii_tw.setToolTip("Individual covalent radii used to determine bonding can be set here")
        self.element_radii_tw.itemClicked.connect(self.on_element_radii_tw_itemClicked)
        self.element_radii_tw.itemChanged.connect(self.on_element_radii_tw_itemChanged)
        self.element_radii_tw.setRowCount(1)
        self.element_radii_tw.blockSignals(False)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        self.element_radii_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel("Atomic radii", self), self.element_radii_tw)
        # Add number of molecules found
        self.molecules_le = QLineEdit(self)
        self.molecules_le.setEnabled(False)
        self.molecules_le.setText(f"{self.number_of_molecules}")
        self.molecules_le.setToolTip("The bonding tolerances can change the number of molecules found")
        label = QLabel("Number of molecules found", self)
        label.setToolTip("The bonding tolerances can change the number of molecules found")
        form.addRow(label, self.molecules_le)
        #
        # The plotting width of bar
        #
        self.width_sp = QDoubleSpinBox(self)
        self.width_sp.setRange(0.01,2.0)
        self.width_sp.setSingleStep(0.01)
        self.width_sp.setDecimals(2)
        self.width_sp.setValue(self.settings["Bar width"])
        self.width_sp.setToolTip("Change the width of the bars - should be between 0 and 1")
        self.width_sp.valueChanged.connect(self.on_width_changed)
        label = QLabel("Bar width", self)
        label.setToolTip("Change the width of the bars - should be between 0 and 1")
        form.addRow(label, self.width_sp)
        #
        # Set the plot title
        #
        self.title_le = QLineEdit(self)
        self.title_le.setToolTip("Set the plot title")
        self.title_le.setText(self.settings["title"])
        self.title_le.textChanged.connect(self.on_title_changed)
        label = QLabel("Plot title", self)
        label.setToolTip("Set the plot title")
        form.addRow(label, self.title_le)
        #
        # Add a comb box to select which type of plot
        #
        self.plottype_cb = QComboBox(self)
        self.plottype_cb.setToolTip("The energy can be decomposed either according to internal vs external motion or into a molecular based decompostion")
        self.plottype_cb.addItems(self.plot_types)
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        self.plottype_cb.currentIndexChanged.connect(self.on_plottype_cb_changed)
        label = QLabel("Choose the plot type", self)
        label.setToolTip("The energy can be decomposed either according to internal vs external motion or into a molecular based decompostion")
        form.addRow(label, self.plottype_cb)
        #
        # Add the matplotlib figure to the bottom
        #
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
        #if self.notebook.spreadsheet is not None:
        #    self.writeSpreadsheet()
        #QCoreApplication.processEvents()

    def on_element_radii_tw_itemClicked(self,item):
        """Handle the item clicked event for a TableWidget related to element radii.

        Parameters
        ----------
        item : QTableWidgetItem
            The item from the TableWidget that was clicked.

        Returns
        -------
        None

        Notes
        -----
        This function unblocks signals for the `element_radii_tw` TableWidget after an item click
        event has occurred, allowing subsequent events to be processed.

        """        
        self.element_radii_tw.blockSignals(False)

    def on_element_radii_tw_itemChanged(self,item):
        """Respond to changes in element radii within a widget.

        Parameters
        ----------
        item : QTableWidgetItem
            The table widget item that was changed. This item contains the new radius value as its text, and its column indicates which element's radius was modified.

        Returns
        -------
        None

        """        
        if self.reader is None:
            return
        col = item.column()
        try:
            debugger.print("Changing the element radius",col,item.text())
            self.settings["Radii"][col] = float(item.text())
            self.calculate()
            self.plot()
            if self.notebook.viewerTab is not None:
                self.notebook.viewerTab.requestRefresh()
        except Exception:
            debugger.print("Failed Changing the element radius",col,item.text())
            pass

    def set_radii_tw(self):
        """Set or update the atomic radii in the GUI's table widget based on the current settings and active file.

        This method updates the radii configuration within the element radii table widget of the GUI. It first ensures that
        a reader object, a program, and a filename are set and available. If any of these are missing, the process is aborted
        with a message. It then checks if custom radii are set in the settings; if not, it defaults to using pre-defined radii.
        The species and their corresponding radii are then used to populate the table widget. After the table is populated,
        it updates the settings with the possibly new radii values.

        Parameters
        ----------
        None 

        Returns
        -------
        None

        """        
        self.reader = self.notebook.mainTab.reader
        program = self.notebook.mainTab.settings["Program"]
        filename = self.notebook.mainTab.getFullFileName()
        if self.reader is None:
            debugger.print("set_radii_tw aborting - no reader")
            return
        if program == "":
            debugger.print("set_radii_tw aborting - no program")
            return
        if filename == "":
            debugger.print("set_radii_tw aborting - no filename")
            return
        debugger.print("set_radii_tw starting")
        self.element_radii_tw.blockSignals(True)
        self.species = self.reader.getSpecies()
        if self.settings["Radii"] is None:
            radii = [ self.element_radii[el] for el in self.species ]
        else:
            radii = self.settings["Radii"]
            for sp,rad in zip(self.species,self.settings["Radii"]):
               self.element_radii[sp] = rad
        self.element_radii_tw.setColumnCount(len(self.species))
        self.element_radii_tw.setHorizontalHeaderLabels(self.species)
        self.element_radii_tw.setVerticalHeaderLabels([""])
        for i,radius in enumerate(radii):
            qw = QTableWidgetItem()
            qw.setText(f"{radius:.6f}")
            qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.element_radii_tw.setItem(0,i, qw )
        self.element_radii_tw.blockSignals(False)
        self.settings["Radii"] = radii
        debugger.print("set_radii_tw finishing")
        return

    def setCovalentRadius(self,element,radius):
        """Set the covalent radius for a given element and update the plot.

        Parameters
        ----------
        element : str
            The chemical symbol of the element for which to set the covalent radius.
        radius : float
            The new covalent radius for the element.

        Returns
        -------
        None

        """        
        self.element_radii[element] = radius
        self.set_radii_tw()
        self.calculate()
        self.plot()

    def writeSpreadsheet(self):
        """Write analysis data into a selected worksheet in the notebook's spreadsheet.

        This method assumes the existence of a spreadsheet object within the notebook
        attribute of the class instance. It selects a worksheet named 'Analysis',
        writes a header row and data rows containing calculated percentages of
        vibrational mode contributions for molecules and internal/external modes.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        if self.notebook.spreadsheet is None:
            return
        sp = self.notebook.spreadsheet
        sp.selectWorkSheet("Analysis")
        sp.delete()
        sp.writeNextRow(["Analysis of the vibrational modes into percentage contributions for molecules and internal/external modes"], row=0,col=1)
        headers = ["Mode","Frequency (cm-1)", "Centre of mass %","Rotational %", "Vibrational %"]
        for mol in range(self.number_of_molecules):
            headers.append("Molecule "+str(mol)+" %")
        sp.writeNextRow(headers,col=1)
        for imode,(freq,energies) in enumerate(zip(self.frequencies_cm1,self.mode_energies)):
           tote,cme,rote,vibe, molecular_energies = energies
           tote = max(tote,1.0E-8)
           output = [ imode+1, freq, 100*cme/tote, 100*rote/tote, 100*vibe/tote ]
           for e in molecular_energies:
               output.append(100*e/tote)
           sp.writeNextRow(output,col=1,check=1)


    def on_width_changed(self,value):
        """Handle changes to the width property.

        This method is called when the width property of an object is changed. It updates the stored width value in the object's settings and then re-plots the object to reflect the new width.

        Parameters
        ----------
        value : int or float
            The new value for the width.

        Returns
        -------
        None

        """        
        debugger.print("on width changed ", value)
        self.settings["Bar width"] = value
        self.plot()

    def on_scale_changed(self,value):
        """Handle the event when the scale setting is changed.

        This function is typically connected to a signal that is emitted when the
        scaling factor for the covalent radius is adjusted in a graphical user interface.
        It updates the stored 'Covalent radius scaling' value, marks that a refresh is
        required, calculates necessary data based on new settings, and finally triggers
        a plot update.

        Parameters
        ----------
        value : float
            The new scale value for the covalent radius.

        Returns
        -------
        None

        """        
        debugger.print("on scale_le changed ", value)
        self.settings["Covalent radius scaling"] = value
        self.refreshRequired = True
        self.calculate()
        self.plot()

    def on_tolerance_changed(self,value):
        """Handle the event when the tolerance value changes.

        This function updates the 'Bonding tolerance' setting based on the provided value, marks the system as requiring a refresh, and then recalculates and replots the data.

        Parameters
        ----------
        value : float
            The new tolerance value.

        Returns
        -------
        None

        """        
        debugger.print("on_tolerance_le changed ", value)
        self.settings["Bonding tolerance"] = value
        self.refreshRequired = True
        self.calculate()
        self.plot()

    def on_title_changed(self,text):
        """Handle title change events.

        Parameters
        ----------
        text : str
            The new title text.

        Returns
        -------
        None

        """        
        self.settings["title"] = text
        if self.subplot is not None:
            self.subplot.set_title(self.settings["title"])
            self.canvas.draw_idle()
        debugger.print("on title change ", self.settings["title"])

    def on_vmin_changed(self):
        """Handle the change in minimum value of frequency.

        Parameters
        ----------
        text : str
            The new title text.

        Returns
        -------
        None

        """        
        self.vmin_sb.blockSignals(True)
        vmin = self.vmin_sb.value()
        vmax = self.vmax_sb.value()
        if vmin < vmax:
            self.settings["Minimum frequency"] = vmin
            debugger.print("on_vmin_changed new value", self.settings["Minimum frequency"])
        else:
            self.vmin_sb.setValue(self.settings["Maximum frequency"]-1)
            self.settings["Minimum frequency"] = self.settings["Maximum frequency"]-1
            debugger.print("on_vmin_changed restricting value to", self.settings["Minimum frequency"])
        self.plot()
        self.vmin_sb.blockSignals(False)
        return

    def on_vmax_changed(self):
        """Handle the change in maximum value for frequency.

        Parameters
        ----------
        text : str
            The new title text.

        Returns
        -------
        None

        """        
        self.vmin_sb.blockSignals(True)
        vmin = self.vmin_sb.value()
        vmax = self.vmax_sb.value()
        if vmax > vmin:
            self.settings["Maximum frequency"] = vmax
            debugger.print("on_vmax_changed new value", self.settings["Maximum frequency"])
        else:
            self.vmin_sb.setValue(self.settings["Minimum frequency"]+1)
            self.settings["Maximum frequency"] = self.settings["Minimum frequency"]+1
            debugger.print("on_vmax_changed restricting value to", self.settings["Maximum frequency"])
        self.plot()
        self.vmin_sb.blockSignals(False)
        return

    def requestRefresh(self):
        """Mark the instance as requiring a refresh.

        Sets the instance attribute `refreshRequired` to True, indicating that a refresh is necessary.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.refreshRequired = True

    def refresh(self, force=False):
        """Refresh the widget state, optionally enforcing refresh.

        This method updates the widgets' states and values according to the current settings. It first checks if a refresh is required or if the `force` parameter is set to `True`. If neither condition is met, it exits early. Otherwise, it proceeds to block signals from all child QWidget instances to avoid unwanted signal emission during state update. It then updates various UI components with new settings values, calculates and plots according to the updated settings, and finally re-enables signals for all child QWidget instances.

        Parameters
        ----------
        force : bool, optional
            If set to `True`, the widget refresh is forced even if it is deemed not required. Defaults to `False`.

        Returns
        -------
        None

        """        
        if not self.refreshRequired and not force:
            debugger.print("return with no refresh", self.refreshRequired, force)
            return
        debugger.print("Refreshing widget")
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        self.vmin_sb.setValue(self.settings["Minimum frequency"])
        self.vmax_sb.setValue(self.settings["Maximum frequency"])
        self.scale_sp.setValue(self.settings["Covalent radius scaling"])
        self.tolerance_sp.setValue(self.settings["Bonding tolerance"])
        self.molecules_le.setText(f"{self.number_of_molecules}")
        self.width_sp.setValue(self.settings["Bar width"])
        self.title_le.setText(self.settings["title"])
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        self.set_radii_tw()
        self.calculate()
        self.plot()
        #
        # Unlock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        return

    def on_plottype_cb_changed(self, index):
        """Handle a change in plot type selection.

        Parameters
        ----------
        index : int
            The new index that indicates the selected plot type.

        Returns
        -------
        None

        """        
        self.plot_type_index = index
        debugger.print("Plot type index changed to ", self.plot_type_index)
        self.plot()

    def calculate(self):
        """Perform calculations for the current object state including molecular contents, normal modes, and energy distributions.

        This method orchestrates the calculation process, setting up necessary parameters, processing molecular contents
        based on the provided settings, calculating normal modes and their mass weighting, and finally computing the energy
        distribution among the modes. It updates the UI and internal state as necessary based on the calculation results.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        The function interacts with several attributes of the object it belongs to, including:
        - `self.notebook.mainTab.settings`: Access to settings specific to the main tab.
        - `self.notebook.mainTab.reader`: The reader associated with the main tab for accessing file data.
        - `self.species` and `self.settings['Radii']`: Used for setting up element radii.
        - `self.notebook.settingsTab.settings`, `self.notebook.settingsTab.frequencies_cm1`, and 
        - `self.notebook.settingsTab.mass_weighted_normal_modes`: Access to settings and data specific to the settings tab.
        - `self.cell_of_molecules`, `self.number_of_molecules`: Modifies and utilizes these attributes
        - to hold molecular content, and the number of molecules, respectively.

        """        
        debugger.print("calculate")
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings["Program"]
        filename = self.notebook.mainTab.getFullFileName()
        self.reader = self.notebook.mainTab.reader
        if self.reader is None:
            return
        if program == "":
            return
        if filename == "":
            return
        for sp,rad in zip(self.species,self.settings["Radii"]):
           self.element_radii[sp] = rad
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        mass_weighted_normal_modes = self.notebook.settingsTab.mass_weighted_normal_modes
        scale = self.settings["Covalent radius scaling"]
        tolerance = self.settings["Bonding tolerance"]
        # Find the last unit cell read by the reader and its masses
        # Take a copy of this cell, as we do not want the reader copy of the cell changing
        self.cell_of_molecules = copy.deepcopy(self.reader.get_unit_cell())
        if self.cell_of_molecules is None:
            return
        nmols = self.cell_of_molecules.calculate_molecular_contents(scale=scale,
                                                  tolerance=tolerance,
                                                  radii=self.element_radii)
        # if the number of molecules has changed then tell the viewerTab that the cell has changed
        if self.number_of_molecules != self.cell_of_molecules.get_number_of_molecules():
            if self.notebook.viewerTab is not None:
                self.notebook.viewerTab.requestRefresh()
            self.number_of_molecules = nmols
        self.molecules_le.setText(f"{self.number_of_molecules}")
        # Calulate the distribution in energy for the normal modes
        mode_energies = Calculator.calculate_energy_distribution(self.cell_of_molecules, 
                                                                 self.frequencies_cm1,
                                                                 mass_weighted_normal_modes)
        # Deal with degeneracies
        degenerate_list = [ [] for f in self.frequencies_cm1]
        for i,fi in enumerate(self.frequencies_cm1):
            for j,fj in enumerate(self.frequencies_cm1):
                if abs(fi-fj) < 1.0E-5:
                    degenerate_list[i].append(j)
        self.mode_energies = []
        for i in range(len(self.frequencies_cm1)):
            tote,cme,rote,vibe,mole = mode_energies[i]
            tote = max(tote,1.0E-8)
            sums = [0.0]*5
            sume = [0.0]*len(mole)
            degeneracy = len(degenerate_list[i])
            for j in degenerate_list[i]:
                tote,cme,rote,vibe,mole = mode_energies[j]
                tote = max(tote,1.0E-8)
                sums[0] += tote / degeneracy
                sums[1] += cme / degeneracy
                sums[2] += rote / degeneracy
                sums[3] += vibe / degeneracy
                for k,e in enumerate(mole):
                    sume[k] += e / degeneracy
                sums[4] = sume
            self.mode_energies.append(sums)
        # Flag that a recalculation is not needed
        self.refreshRequired = False
        QApplication.restoreOverrideCursor()

    def plot(self):
        """Plot the spectroscopic data based on the selected plot type.

        This method selects between two plotting strategies: 'internal_external'
        or 'molecular', depending on the state of `plot_type_index`. If no data
        reader is associated with the instance or if there are no frequencies
        to plot, the function will terminate without plotting.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method relies on `self.reader` to be non-`None`, `self.frequencies_cm1` to
        contain frequency data, and `self.plot_type_index` to determine the plotting
        strategy. It calls `self.plot_internal_external()` or `self.plot_molecular()`
        based on the value of `self.plot_type_index`.

        """        
        if self.reader is None:
            return
        if len(self.frequencies_cm1) <= 0:
            return
        if self.plot_type_index == 0:
            self.plot_internal_external()
        else:
            self.plot_molecular()

    def plot_molecular(self):
        """Plot the molecular composition of vibrational energy.

        This method visualizes the distribution of vibrational energy across different molecules
        within specified frequency limits. It creates a bar chart representing the percentage of 
        energy attributed to each molecule per mode within the set frequency range.

        Parameters
        ----------
        None

        Returns
        -------
        None


        Notes
        -----
        - `self.settings` should contain 'Minimum frequency', 'Maximum frequency', and 'Bar width' keys.
        - `self.mode_energies` must be a list where each element is a tuple containing mode energies
          (total energy, center of mass energy, rotational energy, vibrational energy, molecular energy).
        - `self.frequencies_cm1` should be a list of mode frequencies in cm^-1.
        - `self.number_of_molecules` indicates the number of different molecules considered in the analysis.
        - The method updates `self.subplot` with the new plot and requires `self.figure` and `self.canvas`
          for plotting and refreshing the plot display, respectively.
        - Utilizes `matplotlib.axes.Axes.bar` for plotting and `matplotlib.ticker.MaxNLocator` for setting
          tick locations on the x-axis.
        - Assumes matplotlib, specifically the `matplotlib.pyplot` and `matplotlib.ticker.MaxNLocator` classes,
          are appropriately imported or accessible within the scope.

        """        
        self.subplot = None
        self.figure.clf()
        self.subplot = self.figure.add_subplot(111)
        self.subplot.xaxis.set_major_locator(MaxNLocator(integer=True))
        xlabel = "Mode Number"
        ylabel = "Percentage energy"
        # Decide which modes to analyse
        mode_list = []
        mode_list_text = []
        vmin = self.settings["Minimum frequency"]
        vmax = self.settings["Maximum frequency"]
        tote,cme,rote,vibe,mole = self.mode_energies[0]
        tote = max(tote,1.0E-8)
        mol_energies = [ [] for _ in range(self.number_of_molecules) ]
        mol_bottoms  = [ [] for _ in range(self.number_of_molecules) ]
        for imode, frequency in enumerate(self.frequencies_cm1):
            if frequency >= vmin and frequency <= vmax:
                mode_list.append(imode+1)
                mode_list_text.append(str(imode+1))
                tote,cme,rote,vibe,mole = self.mode_energies[imode]
                tote = max(tote,1.0E-8)
                for i,mol in enumerate(mole):
                    mol_energies[i].append(100.0*mol/tote)
                    if i == 0:
                        mol_bottoms[i].append(0.0)
                    else:
                        mol_bottoms[i].append(mol_bottoms[i-1][-1]+mol_energies[i-1][-1])
        if len(mode_list) < 3:
            return
        width = self.settings["Bar width"]
        plots = []
        colours = ["y","b","r","c","m","k"]
        for i,(energies,bottoms) in enumerate(zip(mol_energies, mol_bottoms )):
            plots.append(self.subplot.bar(mode_list,energies,width, bottom=bottoms,color=colours[i%6]))
        legends = []
        for i in range(self.number_of_molecules):
            legends.append("Molecule "+str(i))
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        self.subplot.legend( plots, legends)
        self.subplot.set_title("Molecular Composition of Vibrational Energy")
        self.canvas.draw_idle()

    def plot_internal_external(self):
        """Plot the internal and external composition of vibrational energy for modes within specified frequency ranges.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - This method modifies the subplot attribute by adding the bar plots for the internal and external energy contributions of different vibrational modes within the specified frequency ranges.
        - This method utilizes the 'Minimum frequency' and 'Maximum frequency' values from the settings attribute to filter modes for plotting.
        - The method sets the X-axis labels to mode numbers and the Y-axis label to percentage energy.
        - If the number of modes within the specified frequency range is less than 3, the method returns early without performing any plotting.
        - The method updates the canvas with the newly plotted data.

        """        
        self.subplot = None
        self.figure.clf()
        self.subplot = self.figure.add_subplot(111)
        self.subplot.xaxis.set_major_locator(MaxNLocator(integer=True))
        xlabel = "Mode Number"
        ylabel = "Percentage energy"
        # Decide which modes to analyse
        mode_list = []
        mode_list_text = []
        cme_energy = []
        rot_energy = []
        vib_energy = []
        vib_bottom = []
        mol_energy = []
        vmin = self.settings["Minimum frequency"]
        vmax = self.settings["Maximum frequency"]
        colours = ["y","b","r","c","m","k"]
        for imode, frequency in enumerate(self.frequencies_cm1):
            if frequency >= vmin and frequency <= vmax:
                mode_list.append(imode+1)
                mode_list_text.append(str(imode+1))
                tote,cme,rote,vibe,molecular_energies = self.mode_energies[imode]
                molecular_energies = np.array(molecular_energies)
                tote = max(tote,1.0E-8)
                cme_energy.append(cme/tote*100.0)
                rot_energy.append(rote/tote*100.0)
                vib_energy.append(vibe/tote*100.0)
                vib_bottom.append( (cme+rote)/tote*100.0 )
                mol_energy.append(molecular_energies/tote*100)
        if len(mode_list) < 3:
            return
        width = self.settings["Bar width"]
        p1 = self.subplot.bar(mode_list,cme_energy,width,color=colours[0])
        p2 = self.subplot.bar(mode_list,rot_energy,width,bottom=cme_energy,color=colours[1])
        p3 = self.subplot.bar(mode_list,vib_energy,width,bottom=vib_bottom,color=colours[2])
        plots = ( p1[0], p2[0], p3[0] )
        legends = ("translation","rotation","vibration")
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        self.subplot.legend( plots, legends)
        self.subplot.set_title("Internal-External Composition of Vibrational Energy")
        self.canvas.draw_idle()

