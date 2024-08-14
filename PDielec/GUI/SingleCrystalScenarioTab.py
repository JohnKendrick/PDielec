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
"""SingleCrystalScenarioTab module."""
import copy
from functools import partial
from itertools import product

import numpy as np
from qtpy.QtCore import QCoreApplication, QSize, Qt
from qtpy.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QSpinBox,
    QStyle,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)
from scipy import signal

import PDielec.GTMcore as GTM
from PDielec import Materials
from PDielec.Constants import speed_light_si
from PDielec.GUI.ScenarioTab import ScenarioTab
from PDielec.GUI.SingleCrystalLayer import ShowLayerWindow, SingleCrystalLayer
from PDielec.Materials import MaterialsDataBase
from PDielec.Utilities import Debug

thickness_conversion_factors = {"ang":1.0E-10, "nm":1.0E-9, "um":1.0E-6, "mm":1.0E-3, "cm":1.0E-2}
thickness_units = list(thickness_conversion_factors.keys())
# incoherentOptions = ['Coherent','Incoherent (intensity)','Incoherent (phase cancelling)','Incoherent (phase averaging)','Incoherent (non-reflective)'] 
incoherentOptions = ["Coherent","Incoherent (intensity)","Incoherent (phase averaging)","Incoherent (non-reflective)"] 
gtmMethods ={"Coherent":GTM.CoherentLayer,
             "Incoherent (intensity)":GTM.IncoherentIntensityLayer,
             "Incoherent (phase cancelling)":GTM.IncoherentPhaseLayer,
             "Incoherent (phase averaging)":GTM.IncoherentAveragePhaseLayer,
             "Incoherent (non-reflective)":GTM.IncoherentThickLayer} 

def solve_single_crystal_equations( 
        layers                        ,
        mode                          ,
        theta                         ,
        phi                           ,
        psi                           ,
        angleOfIncidence              ,
        exponent_threshold            ,
        v                             ,
        ):
    """Solve single crystal equations, system is a GTM system.

    Parameters
    ----------
    layers : list
        A list of material layers (their permittivity functions).
    mode : {'Transfer matrix', 'Scattering matrix'}
        The mode of calculation.
    theta : float
        The theta angle of the slab.
    phi : float
        The phi angle of the slab.
    psi : float
        The psi angle of the slab.
    angleOfIncidence : float
        The angle of incidence.
    exponent_threshold : float
        The largest exponent allowed in the calculation of the propagation matrix.
    v : float
        The frequency of the light in cm^-1.

    Returns
    -------
        A tuple of :
            - v                 : the frequency
            - r                 : the reflectance coefficients
            - R                 : the reflectance 
            - t                 : the transmittance coefficients
            - T                 : the transmittance 
            - epsilon           : the DFT permittivity
            - errors            : a list of errors
            - largest_exponent  : the largest exponent

    """
    # Create superstrate from the first layer
    superstrate      = GTM.SemiInfiniteLayer(layers[0],exponent_threshold=exponent_threshold)
    # Create substrate from the last layer
    substrate        = GTM.SemiInfiniteLayer(layers[-1],exponent_threshold=exponent_threshold)
    selectedLayers = layers[1:-1]
    gtmLayers = []
    # Create layers from all the layers between first and last
    for layer in selectedLayers:
        incoherentOption = layer.getIncoherentOption()
        gtmLayers.append(gtmMethods[incoherentOption](layer, exponent_threshold=exponent_threshold))
    # Creat the system with the layers 
    if mode == "Scattering matrix":
        system = GTM.ScatteringMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtmLayers)
    else:
        system = GTM.TransferMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtmLayers)
    # Rotate the dielectric constants to the laboratory frame
    # This is a global rotation of all the layers.
    system.substrate.set_euler(theta, phi, psi)
    system.superstrate.set_euler(theta, phi, psi)
    for layer in system.layers:
        layer.set_euler(theta, phi, psi)
    # 
    # convert cm-1 to frequency
    #
    freq = v * speed_light_si * 1e2
    system.initialize_sys(freq)
    zeta_sys = np.sin(angleOfIncidence)*np.sqrt(system.superstrate.epsilon[0,0])
    system.calculate_GammaStar(freq, zeta_sys)
    r, R, t, T = system.calculate_r_t(zeta_sys)
    epsilon = system.layers[0].epsilon if len(system.layers) > 0 else system.substrate.epsilon
    errors,largest_exponent = system.overflowErrors()
    return v,r,R,t,T,epsilon,errors,largest_exponent

class SingleCrystalScenarioTab(ScenarioTab):
    """A tab class for handling a single crystal scenario.

    This class extends a generic scenario tab (:class:`~PDielec.GUI.ScenarioTab.ScenarioTab`) with functionalities specific to handling single crystal scenarios, including setting up the scenario environment, defining the crystal orientation, layer information, and calculation modes (e.g., using Transfer Matrix or Scattering Matrix methods). It supports incoherence effects through partial incoherence simulations and provides tools for fine-tuning simulation parameters such as incidence angles, layer materials, thicknesses, and orientations.

    Parameters
    ----------
    parent : QWidget
        The parent widget (typically a main window or a larger container) to which this tab will belong.
    debug : bool, optional
        A flag to indicate whether debugging information should be printed to the console. Default is False.

    Attributes
    ----------
    refreshRequired : bool
        Flag indicating whether the widget's data needs to be refreshed.
    calculationRequired : bool
        Flag indicating whether a new calculation is required based on changes in parameters or settings.
    scenarioType : str
        A string representing the type of scenario. For this class, it is set to 'Single crystal'.
    refreshRequired : bool
        Indicates whether the scenario settings have been changed and thus require the scenario to be redrawn or recalculated.
    noCalculationsRequired : int
        The number of calculations required for the simulation. This value may change depending on the specifics of the scenario configuration (e.g., the inclusion of incoherent effects requiring multiple sample calculations).
    settings : dict
        A dictionary holding various scenario parameters and their values, which include settings for the type of analysis, material layers, angles, and method-specific parameters.
    materialNames : list
        A list of names for materials that can be used in the scenario layers.
    p_reflectance, s_reflectance, p_transmittance, s_transmittance, p_absorbtance, s_absorbtance : list
        Lists holding computed values for reflectance, transmittance, and absorbance for both P and S polarization modes.
    epsilon : list
        A list of dielectric function values computed during simulation.
    layers : list
        A list of layer configurations used in the simulation, each containing information about the material, thickness, orientation, and other layer-specific properties.

    Methods
    -------
    There are several methods within the class for handling events (e.g., button clicks, combo box activations), performing calculations, redrawing tables, and managing layer settings. These include methods for adding, deleting, or altering layers; opening material databases; changing simulation settings; calculating and displaying results; and initializing the user interface components relevant to the single crystal scenario.

    angleOfIncidenceWidget
        Create a widget to set the angle of incidence
    average_incoherent_calculator
        Calculates the incoherent light reflectance and transmittance
    calculate
        Perform the calculations as define in the GUI
    coherent_calculator
        Calculate the coherent light reflectance and transmittance
    createToolBar
         Create the tool bar used for the material layer
    createToolBarDeleteButton
         Create the delete button in the tool bar
    createToolBarMoveDownButton
         Create the move down button in the tool bar
    createToolBarMoveUpButton
         Create the move up button in the tool bar
    deleteLayer
         Delete a layer from the list of layers
    drawLayerTable
         Draw the layer table widget, each layer has a line in the table
    generateLayerSettings
         Generate the settings dictionary with an item for each layer
    getDielectricLayerIndex
         Get the index of the dielectric layer
    getMaterialFromDataBase
         Get a material from the database, or define the material from the given permittivity
    get_result
         Return the request information
    get_results
         Calculate all the reflectance, transmittance and absorptance information
    globalAzimuthalWidget
         Set the global azimuthal angle
    greyed_out
         Grey out menu items according to the GUI settings
    moveLayerDown
         Move a layer up the layer table
    moveLayerUp
         Move a layer down the layer table
    newLayerWidget
         Create a combox widget for adding a new layer
    on_angle_of_incidence_sb_changed
         Handle a change to the angle of incidence spin box
    on_azimuthal_angle_sb_changed
         Handle a change to the azimuthal spin box
    on_film_thickness_sb_changed
         Handle a change to the film thickness spin box
    on_global_azimuthal_angle_sb_changed
         Handle a change to the global azimuthal angle spin box
    on_hkl_sb_changed
         Handle a change to the h,k or l spin boxes
    on_layerTable_itemChanged
         Handle a change to an item in the layer table
    on_layer_button_clicked
         Handle a click on a layer name
    on_mode_cb_activated
         Handle a click on the mode combobox
    on_newLayer_cb_activated
         Handle a click on the new combobox
    on_option_cb_activated
         Handle a change to the incoherent option comobox
    on_partially_incoherent_kernel_sb_changed
         Handle a change to the incoherent kernel spinbox
    on_partially_incoherent_polynomial_sb_changed
         Handle a change to the length of incoherent polynomial spinbox
    on_partially_incoherent_samples_sb_changed
         Handle a change to the number of samples in the averaged incoherent method
    on_percentage_partial_incoherence_sb_changed
         Handle a change to the maximum percentage incoherence
    on_print_button_clicked
         Handle a print button click
    on_thickness_units_cb_activated
         Activate the thickness units combobox
    openDB_button_clicked
         Handle a click on the opend database button
    partialIncoherenceWidget
         Create a partialIncoherence widget
    partially_incoherent_calculator
         Calculate incoherent light scattering using a partially incoherent method
    printLayerSettings
         Print out the layer settings dictionary
    redrawLayerTable
         Redraw the layer table
    redrawLayerTableRow
         Redraw a row of the layer table
    refresh
         Perform a refresh of the GUI
    setMaterialNames
         Read the material names from the database
    set_noCalculationsRequired
         Calculate the number of calculations required
    settings2Layers
         Read the settings dictionary and create the necessary layers
    smoothingWidget
         A widget to handle the smoothing information

    """

    def __init__(self, parent, debug=False ):
        """Initialize the Single Crystal Scenario Tab.

        This initializer sets up the GUI components, populates settings with default values, 
        and prepares the tab for interaction. It inherits from a parent scenario, initializes
        a debug mode if required, and sets various default parameters for the single crystal 
        scenario, including layer information and calculation modes.

        Parameters
        ----------
        parent : QWidget
            The parent widget or window, typically the main application window or a central widget
            that this tab will be a part of.
        debug : bool, optional
            A boolean flag to indicate whether debug messages should be printed to the console
            or log. Defaults to False if not specified.

        Attributes
        ----------
        refreshRequired : bool
            Indicates if the view needs refreshing.
        calculationRequired : bool
            Indicates if a new calculation is required based on changed parameters.
        scenarioType : str
            The type of scenario, hardcoded as 'Single crystal'.
        noCalculationsRequired : int
            Number of calculations required, initialized to 1.
        settings : dict
            A dictionary of settings for the simulation, initialized with default values.
        materialNames : list
            List of material names used in the layers.
        p_reflectance, s_reflectance, p_transmittance, s_transmittance, p_absorbtance, s_absorbtance : list
            Lists for storing calculations results.
        epsilon : list
            List for storing dielectric constants.
        layers : list
            List containing layer information.
        notebook : QWidget
            Reference to the parent notebook or main application window.
        reader : Object
            Reference to an object capable of reading input data (e.g., unit cells, materials).
        cell : Object or None
            The last read unit cell, if any.
        exponent_threshold : int
            Threshold for exponent values in calculations, initialized to 11000.
        DataBase : MaterialsDataBase
            Reference to a materials database object.
        number_of_average_incoherent_layers : int
            Counter for the average incoherent layers, initialized to 0.

        Notes
        -----
        The initializer also connects UI elements with corresponding event handlers, sets up
        a QVBoxLayout for the tab, and adds UI components systematically to form a structured 
        layout for the user interface.

        """        
        ScenarioTab.__init__(self,parent)
        global debugger
        debugger = Debug(debug,"SingleCrystalScenarioTab:")
        debugger.print("Start:: initialiser")
        self.refreshRequired = True
        self.calculationRequired = True
        self.scenarioType = "Single crystal"
        self.noCalculationsRequired = 1
        self.settings["Scenario type"] = "Single crystal"
        self.settings["Global azimuthal angle"] = 0.0
        self.settings["Angle of incidence"] = 0.0
        self.settings["Mode"] = "Transfer matrix"
        self.settings["Mode"] = "Scattering matrix"
        self.settings["Frequency units"] = "wavenumber"
        self.settings["Partially incoherent samples"] = 20
        self.settings["Percentage partial incoherence"] = 0
        self.settings["Filter kernel size"] = 1
        self.settings["Filter polynomial size"] = 3
        # Define a default superstrate/dielectric/substrate system
        self.settings["Layer material names"]     = ["air",     "Dielectric layer","air"     ]
        self.settings["Layer hkls"]               = [ [0,0,0],  [0,0,1],           [0,0,0]   ]
        self.settings["Layer azimuthals"]         = [ 0,         0,                 0        ]
        self.settings["Layer thicknesses"]        = [ 1,         1,                 1        ]
        self.settings["Layer thickness units"]    = ["um",      "um",              "um"      ]
        self.settings["Layer incoherent options"] = ["Coherent","Coherent",        "Coherent"]
        self.settings["Layer dielectric flags"]   = [ False,     True,              False    ]
        # The maximum allowed thickness of a layer in metres
        # used to subdivide thicker films into many thinner films
        # if zero no subdivision is performed
        self.settings["Percentage average incoherence"] = 100
        self.settings["Number of average incoherence samples"] = 10
        self.number_of_average_incoherent_layers = 0
        self.materialNames = []
        self.p_reflectance = []
        self.s_reflectance = []
        self.p_transmittance = []
        self.s_transmittance = []
        self.p_absorbtance = []
        self.s_absorbtance = []
        self.epsilon = []
        self.layers = []
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Get the last unit cell in the reader
        self.cell = None
        if self.reader is not None:
            self.cell = self.reader.get_unit_cell()
        # Set the exponent threshold to be used by GTM
        self.exponent_threshold = 11000    
        # Open the database and get the material names
        self.DataBase = MaterialsDataBase(self.settings["Materials database"],debug=debugger.state())
        self.settings["Materials database"] = self.DataBase.getFileName()
        self.materialNames = self.setMaterialNames()
        # Create the layers - superstrate / dielectric / substrate from the defaults layer settings
        if self.reader is not None:
            self.settings2Layers()
        # Create last tab - SingleCrystalTab
        vbox = QVBoxLayout()
        self.form = QFormLayout()
        #
        # Option to open a database of permittivities for the support
        # label and button are defined the parent class
        self.form.addRow(self.openDB_label, self.openDB_button)
        #
        # Add a name for the database
        # label and button are defined the parent class
        self.form.addRow(self.database_le_label, self.database_le)
        #
        # Chose mode of operation
        #
        self.mode_cb = QComboBox(self)
        self.mode_cb.setToolTip("Set the method for calculating light transmission and reflectanceb;\n Transfer matrix.  This method is fast but can be numerically unstable.\n Scattering matrix. This method is slow but is numerically stable")
        self.mode_cb.addItems( ["Transfer matrix","Scattering matrix"] )
        index = self.mode_cb.findText(self.settings["Mode"], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        self.mode_cb.activated.connect(self.on_mode_cb_activated)
        label = QLabel("Single crystal methodology", self)
        label.setToolTip("Set the method for calculating light transmission and reflectanceb;\n Transfer matrix.  This method is fast but can be numerically unstable.\n Scattering matrix. This method is slow but is numerically stable")
        self.form.addRow(label, self.mode_cb)
        #
        # Define the global azimuthal angle widget
        #
        label,layout = self.globalAzimuthalWidget()
        self.form.addRow(label, layout)
        #
        # Define the angle of incidence widget
        #
        label, layout = self.angleOfIncidenceWidget()
        self.form.addRow(label,layout)
        #
        # Layer information widget
        #
        label = QLabel("Layer information")
        line  = QFrame()
        line.setFrameShape(QFrame.HLine)
        hbox = QHBoxLayout()
        hbox.addWidget(line)
        hbox.setAlignment(Qt.AlignVCenter)
        self.form.addRow(label,hbox)
        self.form.addRow(self.drawLayerTable())
        label = QLabel("    ")
        line  = QFrame()
        line.setFrameShape(QFrame.HLine)
        hbox = QHBoxLayout()
        hbox.addWidget(line)
        hbox.setAlignment(Qt.AlignVCenter)
        self.form.addRow(label,hbox)
        #
        # Partial incoherence widget
        #
        label,layout = self.partialIncoherenceWidget()
        #jk self.form.addRow(label, layout)
        #
        # Smoothing widget
        #
        label,layout = self.smoothingWidget()
        #jk self.form.addRow(label, layout)
        #
        # Add a legend option
        #
        self.legend_le = QLineEdit(self)
        self.legend_le.setToolTip("The legend will be used to describe the results in the plot")
        self.legend_le.setText(self.settings["Legend"])
        self.legend_le.textChanged.connect(self.on_legend_le_changed)
        label = QLabel("Scenario legend")
        label.setToolTip("The legend will be used to describe the results in the plot")
        self.form.addRow(label, self.legend_le)
        #
        #
        # Final buttons for changing/deleting and switching the scenarios
        #
        hbox = self.add_scenario_buttons()
        self.form.addRow(hbox)
        #
        # Spacer item to soak up space
        #
        verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.form.addItem(verticalSpacer)
        vbox.addLayout(self.form)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        debugger.print("Finished:: initialiser")

    def redrawLayerTable(self):
        """Redraw the layer table widget.

        Adds a row for each layer in the gadget.
        The row is drawn using redrawLayerTableRow()
        """
        self.layerTable_tw.setRowCount(1)
        rowCount = 0
        for sequenceNumber,layer in enumerate(self.layers):
            rowCount += 1
            firstLayer = False
            lastLayer = False
            if layer == self.layers[0]:
                firstLayer = True
            elif layer == self.layers[-1]:
                lastLayer = True
            self.layerTable_tw.setRowCount(rowCount)
            self.redrawLayerTableRow(sequenceNumber,layer,rowCount,firstLayer,lastLayer)
        # Add a 'create new layer' button
        rowCount += 1
        newLayer_cb = self.newLayerWidget()
        newLayer_cb.setStyleSheet("Text-align:left")
        self.layerTable_tw.setRowCount(rowCount)
        self.layerTable_tw.setCellWidget(rowCount-1,0,newLayer_cb)

    def redrawLayerTableRow(self,sequenceNumber,layer,rowCount,firstLayer,lastLayer):
        """Draw a row of the layer table.

        Each row has a material name, a thickness (and unit), an h, k, l option and a toolbar
        to move or delete the layer in the list.
        Needed a separate routine for this as there are problems with the lambda code only keeping the last in a list.

        Parameters
        ----------
        sequenceNumber : int
            The sequence number of the layer in the list
        layer : layer object
            The layer itself
        rowCount : int
            The row count
        firstLayer : bool
            True if this is the first layer in the list
        lastLayer : bool
            True if this is the last layer in the list

        Returns
        -------
        None

        """
        # Create a layer button
        material = layer.getMaterial()
        materialName = material.getName()
        layer_button = QPushButton(materialName)
        layer_button.setToolTip("Show the material properties in a new window")
        layer_button.setStyleSheet("Text-align:left")
        layer_button.clicked.connect(lambda x: self.on_layer_button_clicked(x,layer,sequenceNumber))
        self.layerTable_tw.setCellWidget(sequenceNumber,0,layer_button)
        # Handle thickness 
        materialThickness = layer.getThickness()
        thicknessUnit = layer.getThicknessUnit()
        film_thickness_sb = QDoubleSpinBox(self)
        film_thickness_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        film_thickness_sb.setToolTip("Define the thin film thickness in the defined thickness units")
        film_thickness_sb.setRange(0,100000)
        film_thickness_sb.setDecimals(3)
        film_thickness_sb.setSingleStep(0.001)
        film_thickness_sb.setValue(materialThickness)
        film_thickness_sb.valueChanged.connect(lambda x: self.on_film_thickness_sb_changed(x,layer))
        self.layerTable_tw.setCellWidget(sequenceNumber,1,film_thickness_sb)
        # thickness unit
        thickness_unit_cb = QComboBox(self)
        thickness_unit_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thickness_unit_cb.setToolTip("Set the units to be used for thickness; either angs nm, um, mm or cm")
        thickness_unit_cb.addItems( thickness_units )
        index = thickness_unit_cb.findText(thicknessUnit, Qt.MatchFixedString)
        thickness_unit_cb.setCurrentIndex(index)
        thickness_unit_cb.activated.connect(lambda x: self.on_thickness_units_cb_activated(x, layer))
        self.layerTable_tw.setCellWidget(sequenceNumber,2,thickness_unit_cb)
        # define hkl
        h_sb = QSpinBox(self)
        h_sb.setToolTip("Define the h dimension of the unique direction")
        h_sb.setRange(-20,20)
        h_sb.setSingleStep(1)
        h_sb.setValue(layer.getHKL()[0])
        h_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,0,layer))
        h_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        k_sb = QSpinBox(self)
        k_sb.setToolTip("Define the k dimension of the unique direction")
        k_sb.setRange(-20,20)
        k_sb.setSingleStep(1)
        k_sb.setValue(layer.getHKL()[1])
        k_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,1,layer))
        k_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        l_sb = QSpinBox(self)
        l_sb.setToolTip("Define the l dimension of the unique direction")
        l_sb.setRange(-20,20)
        l_sb.setSingleStep(1)
        l_sb.setValue(layer.getHKL()[2])
        l_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,2,layer))
        l_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.layerTable_tw.setCellWidget(sequenceNumber,3,h_sb)
        self.layerTable_tw.setCellWidget(sequenceNumber,4,k_sb)
        self.layerTable_tw.setCellWidget(sequenceNumber,5,l_sb)
        # define azimuthal angle
        azimuthal = layer.getAzimuthal()
        azimuthal_angle_sb = QDoubleSpinBox(self)
        azimuthal_angle_sb.setToolTip("Define the slab azimuthal angle (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal in the laboratory frame can be seen in the laboratory frame information below")
        azimuthal_angle_sb.setRange(-180,360)
        azimuthal_angle_sb.setSingleStep(10)
        azimuthal_angle_sb.setValue(azimuthal)
        azimuthal_angle_sb.valueChanged.connect(lambda x: self.on_azimuthal_angle_sb_changed(x,layer))
        azimuthal_angle_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.layerTable_tw.setCellWidget(sequenceNumber,6,azimuthal_angle_sb)
        # Create a checkbox for coherence/incoherence
        option_cb = QComboBox(self)
        option_cb.setToolTip("Change optional settings for the layer")
        option_cb.addItems( incoherentOptions )
        # We can't use incoherent intensity method with the scattering matrix method
        if self.settings["Mode"] == "Scattering matrix" and layer.getIncoherentOption() == "Incoherent (intensity)":
            layer.setIncoherentOption("Coherent")
        index = option_cb.findText(layer.getIncoherentOption(), Qt.MatchFixedString)
        option_cb.setCurrentIndex(index)
        option_cb.activated.connect(lambda x: self.on_option_cb_activated(x,layer))
        option_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        # Disable the intensity incoherence option if it is using scattering
        if self.settings["Mode"] == "Scattering matrix":
            index = option_cb.findText("Incoherent (intensity)", Qt.MatchFixedString)
            option_cb.model().item(index).setEnabled(False)
        self.layerTable_tw.setCellWidget(sequenceNumber,7,option_cb)
        # Create a toolbar for up down delete
        toolbar = self.createToolBar(layer,sequenceNumber,len(self.layers))
        self.layerTable_tw.setCellWidget(sequenceNumber,8,toolbar)
        # Add a Print option if debug is on
        if debugger.state():
            printButton = QPushButton("Print")
            printButton.setToolTip("Print the permittivity")
            printButton.clicked.connect(lambda x: self.on_print_button_clicked(x,layer))
            self.layerTable_tw.setCellWidget(sequenceNumber,9,printButton)
        if layer.isScalar():
            h_sb.setEnabled(False)
            k_sb.setEnabled(False)
            l_sb.setEnabled(False)
            azimuthal_angle_sb.setEnabled(False)
        if firstLayer or lastLayer:
            film_thickness_sb.setEnabled(False)
            thickness_unit_cb.setEnabled(False)
        return

    def drawLayerTable(self):
        """Draw a table with all the layers in it.

        Each layer has a row in the table showing, name, thickness, h,k,l, azimuthal angle, options and a toolbar

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.layerTable_tw = QTableWidget()
        self.layerTable_tw.setToolTip("Define the layers in the system\nThe layer at the top of the list is the superstrate\nThe layer at the bottom is the substrate\nThe calculated, DFT, permittivity is called the 'Dielectric layer'")
        self.layerTable_tw.itemChanged.connect(self.on_layerTable_itemChanged)
        self.layerTable_tw.setStyleSheet("QTableWidget::item {padding-left: 0px; border; 0px}")
        self.layerTable_tw.verticalHeader().setVisible(False)
        self.layerTable_tw.setShowGrid(False)
        headers = ["Material", "Thickness", "Units", "H", "K", "L", "Azimuthal", "Options", "Move"]
        if debugger.state():
            headers.append("Print")
        self.layerTable_tw.setRowCount(1)
        self.layerTable_tw.setColumnCount(len(headers))
        self.layerTable_tw.setHorizontalHeaderLabels(headers)
        header = self.layerTable_tw.horizontalHeader()
        # Material 
        header.setSectionResizeMode(0,QHeaderView.Stretch)
        # Thickness and units 
        header.setSectionResizeMode(1,QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2,QHeaderView.ResizeToContents)
        # HKL
        header.setSectionResizeMode(3,QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4,QHeaderView.ResizeToContents)
        header.setSectionResizeMode(5,QHeaderView.ResizeToContents)
        # Azimuthal
        header.setSectionResizeMode(6,QHeaderView.ResizeToContents)
        # Options
        header.setSectionResizeMode(7,QHeaderView.ResizeToContents)
        # Move
        header.setSectionResizeMode(8,QHeaderView.ResizeToContents)
        if "Print" in headers:
            # Print
            header.setSectionResizeMode(9,QHeaderView.ResizeToContents)
        if self.reader is not None:
            self.redrawLayerTable()
        return self.layerTable_tw

    def deleteLayer(self,x,layer,layerIndex):
        """Handle a delete layer button press.

        Parameters
        ----------
        x : any
            A dummy variable
        layer : layer object
            the layer to be deleted
        layerIndex : int
            index of the layer to be deleted

        Returns
        -------
        None

        """
        new = layerIndex + 1
        if layerIndex == 0 and layer[new].isTensor():
            #  Only allow scalar materials as the superstrate
            print("New superstrate material must be a scalar dielectric")
            return
        # Delete the layer
        del self.layers[layerIndex]
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def moveLayerUp(self,x,layer,layerIndex):
        """Move a layer up (sequence number gets smaller by 1).

        Parameters
        ----------
        x : any
            A dummy variable
        layer : layer object
            the layer to be deleted
        layerIndex : int
            index of the layer to be deleted

        Returns
        -------
        None

        """
        if layerIndex < 1:
            return
        if layerIndex == 1 and layer.isTensor():
            #  Only allow scalar materials as the superstrate
            print("New superstrate material must be a scalar dielectric")
            return
        new = layerIndex - 1
        item = self.layers[layerIndex]
        self.layers.pop(layerIndex)
        self.layers.insert(new, item)
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def moveLayerDown(self,x,layer,layerIndex):
        """Move a layer down (sequence number gets larger by 1).

        Parameters
        ----------
        x : any
            A dummy variable
        layer : layer object
            the layer to be deleted
        layerIndex : int
            index of the layer to be deleted

        Returns
        -------
        None

        """
        last = len(self.layers) - 1
        if layerIndex >= last:
            return
        new = layerIndex + 1
        if layerIndex == 0 and self.layers[new].isTensor():
            #  Only allow scalar materials as the superstrate
            print("New superstrate material must be a scalar dielectric")
            return
        item = self.layers[layerIndex]
        self.layers.pop(layerIndex)
        self.layers.insert(new, item)
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def createToolBarMoveUpButton(self,layer,layerIndex,nLayers):
        """Create the move up button as part of the layer toolbar.

        Parameters
        ----------
        layer : layer object
            the layer to be deleted
        layerIndex : int
            index of the layer to be deleted
        nLayers : int
            The number of layers

        Returns
        -------
        None

        """
        moveUpButton = QPushButton()
        moveUpButton.setIcon(QApplication.style().standardIcon(QStyle.SP_ArrowUp))
        moveUpButton.clicked.connect(lambda x: self.moveLayerUp(x,layer,layerIndex))
        moveUpButton.setFixedSize(20,20)
        moveUpButton.setIconSize(QSize(20,20))
        moveUpButton.setStyleSheet("border: none;")
        moveUpButton.setToolTip("Move this layer up the list of layers")
        return moveUpButton

    def createToolBarMoveDownButton(self,layer,layerIndex,nLayers):
        """Create the move down button as part of the layer toolbar.

        Parameters
        ----------
        layer : layer object
            the layer to be deleted
        layerIndex : int
            index of the layer to be deleted
        nLayers : int
            The number of layers

        Returns
        -------
        None

        """
        moveDownButton = QPushButton()
        moveDownButton.setIcon(QApplication.style().standardIcon(QStyle.SP_ArrowDown))
        moveDownButton.clicked.connect(lambda x: self.moveLayerDown(x,layer,layerIndex))
        moveDownButton.setFixedSize(20,20)
        moveDownButton.setIconSize(QSize(20,20))
        moveDownButton.setStyleSheet("border: none;")
        moveDownButton.setToolTip("Move this layer down the list of layers")
        return moveDownButton

    def createToolBarDeleteButton(self,layer,layerIndex,nLayers):
        """Create the delete button as part of the layer toolbar.

        Parameters
        ----------
        layer : layer object
            the layer to be deleted
        layerIndex : int
            index of the layer to be deleted
        nLayers : int
            The number of layers

        Returns
        -------
        None

        """
        deleteButton = QPushButton()
        deleteButton.setIcon(QApplication.style().standardIcon(QStyle.SP_DialogCloseButton))
        deleteButton.clicked.connect(lambda x: self.deleteLayer(x,layer,layerIndex))
        deleteButton.setFixedSize(20,20)
        deleteButton.setIconSize(QSize(20,20))
        deleteButton.setStyleSheet("border: none;")
        deleteButton.setToolTip("Delete this layer")
        return deleteButton

    def createToolBar(self,layer,layerIndex,nLayers):
        """Create the tool bar used for the material layer.

        Parameters
        ----------
        layer : a layer object
            Layer is the layer concerned.
        layerIndex : int
            Its index in the list.
        nLayers : int
            The number of layers in the list.

        Returns
        -------
        None

        """
        if nLayers <= 1:
            return None
        frame = QFrame()
        frame_layout = QHBoxLayout()
        frame.setLayout(frame_layout)
        # Create the buttons in different routines because of the lambda function usage
        moveUpButton   =  self.createToolBarMoveUpButton(layer,layerIndex,nLayers)
        moveDownButton =  self.createToolBarMoveDownButton(layer,layerIndex,nLayers)
        deleteButton   =  self.createToolBarDeleteButton(layer,layerIndex,nLayers)
        nextIndex = layerIndex+1
        # disable any buttons that are irrelevant to the layer
        if layerIndex == 0:
            moveUpButton.setEnabled(False)
        if layerIndex == nLayers-1:
            moveDownButton.setEnabled(False)
        if layerIndex == 0 and self.layers[nextIndex].isTensor():
            moveDownButton.setEnabled(False)
            deleteButton.setEnabled(False)
        if layerIndex == 1 and layer.isTensor():
            moveUpButton.setEnabled(False)
        # Add the buttons to the frame and return the frame
        frame_layout.addWidget(moveUpButton)
        frame_layout.addWidget(moveDownButton)
        frame_layout.addWidget(deleteButton)
        return frame
        
    def newLayerWidget(self):
        """Create and return a new layer widget as a QComboBox.

        Returns a QComboBox widget configured for creating a new layer in a graphical user interface.
        The combo box is initialized with a default item for creating a new layer, followed by any existing layer names.

        Parameters
        ----------
        None

        Returns
        -------
        QComboBox
            A QComboBox widget configured with an initial item for creating a new layer and any additional existing material layer names.

        Notes
        -----
        - The first item in the combo box is a placeholder for creating a new layer, indicated by "New layer...".
        - The combo box is connected to the `on_newLayer_cb_activated` method, which should handle the action performed upon selecting an item from the combo box.

        """        
        newLayer_cb = QComboBox()
        newLayer_cb.setToolTip("Create a new layer")
        materialNames = ["New layer..."]
        materialNames += self.materialNames
        newLayer_cb.addItems(materialNames)
        newLayer_cb.setCurrentIndex(0)
        newLayer_cb.activated.connect(self.on_newLayer_cb_activated)
        return newLayer_cb

    def on_layerTable_itemChanged(self,item):
        """Handle a change to the layer table.

        An empty routine as changes are handled by other routines

        Parameters
        ----------
        item : the item changed
            The item which has changed

        Returns
        -------
        None

        """
        print("on_layerTable_itemChanged: ",item)
        return

    def on_newLayer_cb_activated(self,index):
        """Handle a new layer button click.

        Based on the index chosen a new material is created and added to the list of layers.
        The settings dictionary is updated and a refresh is forced.

        Parameters
        ----------
        index : int
            The index of the entry in the combobox

        Returns
        -------
        None

        """
        if index == 0:
            return
        # Subtract 1 from the index because the widget thinks the list includes 'New layer...' at the start
        newMaterialName = self.materialNames[index-1]
        if "manual" in newMaterialName:
            return
        newMaterial = self.getMaterialFromDataBase(newMaterialName)
        hkl = [0,0,0]
        if newMaterial.isTensor():
            hkl = [0,0,1]
        new_layer = SingleCrystalLayer(newMaterial,hkl=hkl,azimuthal=0.0,thickness=1.0,thicknessUnit="um")
        self.layers.append(new_layer)
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired = True
        return

    def on_print_button_clicked(self,x,layer):
        """Print the permittivity for the layer.

        The print button is only visible if in debug mode.
        The permittivity is printed to a csv file with a name containing the material involved.

        Parameters
        ----------
        x : any
            A dummy parameter
        layer : a layer object
            The layer object for which a print of the permittivity has been requested

        Returns
        -------
        None

        """
        material = layer.getMaterial()
        permittivityObject = material.getPermittivityObject()
        name = material.getName()
        name = name.replace(" ","_")
        name += "_permittivity.csv"
        print("Printing permittivity information to",name)
        permittivityObject.print(0.0,2000.0,1.0,file=name)
        return

    def on_option_cb_activated(self,index,layer):
        """Activate the incoherence option combob.

        - The layer's incoherent options are changed
        - The number of calculations requested is updated
        - The settings dictionary is regenerated
        - A refresh is requested

        Parameters
        ----------
        index : int
            The index of the incoherent option combobox
        layer : a layer object
            The layer which will have its incoherent property changed

        Returns
        -------
        None

        """
        debugger.print("on_incoherence_cb_activated", index,layer.getName())
        option = incoherentOptions[index]
        layer.setIncoherentOption(option)
        self.set_noCalculationsRequired()
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_film_thickness_sb_changed(self,value,layer):
        """Handle film thickness spin box change.

        - The thickness of the film in this layer is changed
        - The settings dictionary is regenerated
        - A refresh is requested

        Parameters
        ----------
        value : float
            The thickness of the film in the current thickness units
        layer : a layer object
            The layer which will have its thickness changed

        Returns
        -------
        None

        """
        debugger.print("on_film_thickness_sb_changed", value, layer.getName())
        layer.setThickness(value)
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_thickness_units_cb_activated(self, index, layer):
        """Handle the activation of a thickness unit option in a combo box.

        - the thickness unit of this layer is changed
        - the settings dictionary is updated
        - a refresh is requested

        Parameters
        ----------
        index : int
            The index of the selected thickness unit in the combo box.
        layer : a layer object
            The layer object whose thickness unit is to be updated.

        Returns
        -------
        None

        """        
        debugger.print("Start:: on_thickness_units_cb_activated",index,layer.getName())
        unit = thickness_units[index]
        layer.setThicknessUnit(unit)
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_azimuthal_angle_sb_changed(self,value,layer):
        """Handle azimuthal spin box change.

        - The azimuthal angle of a layer is changed to the value given.
        - The laboratory frame information is updated
        - The settings dictionary is refreshed 
        - A refresh of the GUI is requested

        Parameters
        ----------
        value : float
            The azimuthal angle in degress
        layer : a layer object
             The layer whose azimuthal angle has changed

        Returns
        -------
        None

        """
        debugger.print("on_azimuthal_angl_sb_changed", value, layer.getName())
        layer.setAzimuthal(value)
        layer.changeLabFrameInfo()
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_hkl_sb_changed(self,value,hkorl,layer):
        """Handle a change in h, k, or l.

        The h, k and l parameters determine the surface of the film
        - determine the value of h, k or l
        - set the layer's surface
        - generate the settings dictionary
        - flag that a refresh is necessary

        Parameters
        ----------
        value : int
            The value of h, k or l
        hkorl : int
            A int 0, 1 or 2 representing 'h', 'k', or 'l'
        layer :  a layer object
            The layer whose surface is to be changed

        Returns
        -------
        None

        """
        debugger.print("on_hkl_sb_changed", value)
        hkl = layer.getHKL()
        hkl[hkorl] = value
        layer.setHKL(hkl)
        layer.changeLabFrameInfo()
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def generateLayerSettings(self):
        """Generate the settings dictionary for every layer.

        Parameters
        ----------
        None

        Returns
        -------
        list
            A list of dictionary settings for each layer.

        """
        debugger.print(self.settings["Legend"],"generateLayerSettings")
        self.settings["Layer material names"]  = []
        self.settings["Layer hkls"]            = []
        self.settings["Layer azimuthals"]      = []
        self.settings["Layer thicknesses"]     = []
        self.settings["Layer thickness units"] = []
        self.settings["Layer dielectric flags"] = []
        self.settings["Layer incoherent options"] = []
        for layer in self.layers:
            self.settings["Layer material names"].append(layer.getMaterial().getName())
            self.settings["Layer hkls"].append(layer.getHKL())
            self.settings["Layer azimuthals"].append(layer.getAzimuthal())
            self.settings["Layer thicknesses"].append(layer.getThickness())
            self.settings["Layer thickness units"].append(layer.getThicknessUnit())
            self.settings["Layer dielectric flags"].append(layer.isDielectric())
            self.settings["Layer incoherent options"].append(layer.getIncoherentOption())
        return

    def setMaterialNames(self):
        """Set and append a material name to the list of material names.

        This function retrieves the current list of material names from the database using
        `self.DataBase.getSheetNames()`, appends 'Dielectric layer' to the list, and returns the updated list.

        Parameters
        ----------
        None

        Returns
        -------
        list
            The list of material names including the appended default 'Dielectric layer'.

        """        
        # Get the list of material names from the database
        materialNames = self.DataBase.getSheetNames()
        materialNames.append("Dielectric layer")
        return materialNames

    def printLayerSettings(self,message):
        """Print the settings related to the layers.

        Parameters
        ----------
        message : str
            Message to be printed before the layer settings.

        Returns
        -------
        None

        Notes
        -----
        This function prints the settings associated with layers, including material names, HKL values (Miller indices for the crystallographic planes), azimuthal angles, thicknesses with their units, flags indicating whether the layer is dielectric, and options for incoherent scattering, if applicable.

        This function does not return any value. It is primarily used for debugging or informative purposes to quickly display the current configuration of layer settings within a class that has `settings` as an attribute storing these properties.

        """        
        print(message)
        print(self.settings["Layer material names"])
        print(self.settings["Layer hkls"])
        print(self.settings["Layer azimuthals"])
        print(self.settings["Layer thicknesses"])
        print(self.settings["Layer thickness units"])
        print(self.settings["Layer dielectric flags"])
        print(self.settings["Layer incoherent options"])
        print()

    def settings2Layers(self):
        """Read the layer settings and generate a list of layers.

        This method reads all the entries in the settings dictionary
        to create a list of layers

        Parameters
        ----------
        None

        Returns
        -------
        list
            A list of layers generated from the layer settings.

        """
        debugger.print(self.settings["Legend"],"settings2Layers")
        self.layers = []
        self.materialNames = self.setMaterialNames()
        # Process the settings information and append each layer to the list
        for  name, hkl, azimuthal, thickness, thicknessUnit, dielectricFlag, incoherentOption in zip(
                          self.settings["Layer material names"],
                          self.settings["Layer hkls"],
                          self.settings["Layer azimuthals"],
                          self.settings["Layer thicknesses"],
                          self.settings["Layer thickness units"],
                          self.settings["Layer dielectric flags"],
                          self.settings["Layer incoherent options"]):
            if name not in self.materialNames:
                print("Error material ", name, " not available ", self.materialNames)
                name = "air"
            material = self.getMaterialFromDataBase(name)
            self.layers.append(SingleCrystalLayer(material,hkl=hkl,azimuthal=azimuthal,
                                     thickness=thickness,thicknessUnit=thicknessUnit,
                                     incoherentOption=incoherentOption,dielectricFlag=dielectricFlag))
        return

    def getMaterialFromDataBase(self,name,permittivity=None):
        """Get the given material from the database.

        - If name is 'DielectricLayer' then a material is defined from the crystal permittivity of the DFT information
        - If name is 'Material defined manually' the the permittivity parameter is used
        - Otherwise the name is used to get the material from the database

        Parameters
        ----------
        name : str
            The name of the material. It can also be 'Dielectric layer' or 'Material defined manually'.
        permittivity : float
            The permittivity of a 'Defined manually material' material.

        Returns
        -------
        None

        """
        if name == "Dielectric layer":
            # Create the dielectric material
            crystalPermittivityObject = self.notebook.settingsTab.getCrystalPermittivityObject()
            material = Materials.External(name,permittivityObject=crystalPermittivityObject,cell=self.cell)
        elif name == "Material defined manually":
            material = Materials.Constant("Material defined manually",permittivity=permittivity)
        else:
            # Get the material from the data base
            # set the units for frequency to Hz for all materials
            material = self.DataBase.getMaterial(name)
        return material

    def getDielectricLayerIndex(self):
        """Return the index of the dielectric layer in the list of layers.

        - Sorts through the list of layers and reports the index of the first layer that came from the provided DFT calculation
        - If there is no such layer then the value None is returned

        Parameters
        ----------
        None

        Returns
        -------
        int or None
            The index of the dielectric layer in the list.

        """
        for index,layer in enumerate(self.layers):
            if layer.isDielectric():
                self.dielectricLayer = layer
                return index
        return None

    def on_layer_button_clicked(self,x,layer,layerIndex):
        """Handle a click on the show layer widget.

        A new window is shown with additional information about material in the layer.
        For materials with a 3x3 tensor permittivity the window shows the relationship between the laboratory and crystal coordinates.

        Parameters
        ----------
        x : any
            A dummy variable
        layer : a layer object
            A layer object
        layerIndex : int
            The index of the layer object in the list

        Returns
        -------
        None

        """
        # Create the dialog box with all the information on the layer, work on a copy of the layer
        if layerIndex == 0:
            message = "Superstrate layer"
        elif layerIndex == len(self.layers)-1:
            message = "Substrate layer"
        else:
            message = "Device layer " + str(layerIndex)
        showLayerWindow = ShowLayerWindow(copy.copy(layer),message=message,debug=debugger.state())
        if showLayerWindow.exec():
            # The 'Ok' button was pressed
            # get the new Layer and replace the old one
            self.layers[layerIndex] = showLayerWindow.getLayer()
            self.generateLayerSettings()
            self.refresh(force=True)
            self.refreshRequired = True
        return

    def globalAzimuthalWidget(self):
        """Create a global azimuthal angle widget.

        This widget is shown in the main single crystal scenario tab.
        Any change in this angle affects all layers.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.global_azimuthal_angle_sb = QDoubleSpinBox(self)
        self.global_azimuthal_angle_sb.setToolTip("Define the global azimuthal angle for all layers (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal layers in the laboratory frame is shown in the edit layers window")
        self.global_azimuthal_angle_sb.setRange(-180,360)
        self.global_azimuthal_angle_sb.setSingleStep(10)
        self.global_azimuthal_angle_sb.setValue(self.settings["Global azimuthal angle"])
        self.global_azimuthal_angle_sb.valueChanged.connect(self.on_global_azimuthal_angle_sb_changed)
        label = QLabel("Global azimuthal angle")
        label.setToolTip("Define the global azimuthal angle for all layers (rotation of the layers about the lab Z-axis).\nThe orientation of the crystal layers in the laboratory frame is shown in the edit layers window")
        return label,self.global_azimuthal_angle_sb

    def angleOfIncidenceWidget(self):
        """Create an angle of incidence widget.

        This widget is used in the layer table widget to control the azimuthal angle of a single layer.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.angle_of_incidence_sb = QDoubleSpinBox(self)
        self.angle_of_incidence_sb.setToolTip("Define the angle of incidence, (normal incidence is 0 degrees)")
        self.angle_of_incidence_sb.setRange(0,90)
        self.angle_of_incidence_sb.setSingleStep(5)
        self.angle_of_incidence_sb.setValue(self.settings["Angle of incidence"])
        self.angle_of_incidence_sb.valueChanged.connect(self.on_angle_of_incidence_sb_changed)
        label = QLabel("Angle of incidence")
        label.setToolTip("Define the angle of incidence, (normal incidence is 0 degrees).")
        return label,self.angle_of_incidence_sb

    def partialIncoherenceWidget(self):
        """Create a partial incoherence widget.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        hbox = QHBoxLayout()
        self.percentage_partial_incoherence_sb = QSpinBox()
        self.percentage_partial_incoherence_sb.setToolTip("Define the maximum percentage changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees\nIf the value of the percentage incoherence is zero, no incoherence is calculated")
        self.percentage_partial_incoherence_sb.setRange(0,100)
        self.percentage_partial_incoherence_sb.setSingleStep(1)
        self.percentage_partial_incoherence_sb.setValue(self.settings["Percentage partial incoherence"])
        self.percentage_partial_incoherence_sb.valueChanged.connect(self.on_percentage_partial_incoherence_sb_changed)
        self.partially_incoherent_samples_sb = QSpinBox()
        self.partially_incoherent_samples_sb.setToolTip("Define the number of samples to be used in the calculation of an incoherent spectrum.\nA large number of samples will take a long time but will give smoother results.")
        self.partially_incoherent_samples_sb.setRange(0,10000)
        self.partially_incoherent_samples_sb.setSingleStep(1)
        self.partially_incoherent_samples_sb.setValue(self.settings["Partially incoherent samples"])
        self.partially_incoherent_samples_sb.valueChanged.connect(self.on_partially_incoherent_samples_sb_changed)
        hbox.addWidget(self.percentage_partial_incoherence_sb)
        hbox.addWidget(self.partially_incoherent_samples_sb)
        label = QLabel("Partial incoherence (percentage & samples)")
        label.setToolTip("Define the percentage partial incoherence and the number of samples to be used in the calculation of an incoherent spectrum.\nThe percentage reflects changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence.\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees.\nA large number of samples will take a long time but will give smoother results.\nIf the value of the percentage incoherence is zero, no incoherence is calculated")
        return label,hbox

    def smoothingWidget(self):
        """Create a smoothing widget used by the partial incoherence method.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        hbox = QHBoxLayout()
        self.partially_incoherent_kernel_sb = QSpinBox()
        self.partially_incoherent_kernel_sb.setRange(1,1001)
        self.partially_incoherent_kernel_sb.setSingleStep(2)
        self.partially_incoherent_kernel_sb.setValue(self.settings["Filter kernel size"])
        self.partially_incoherent_kernel_sb.valueChanged.connect(self.on_partially_incoherent_kernel_sb_changed)
        self.partially_incoherent_kernel_sb.setToolTip("Define the kernel size for the smoothing of incoherent spectra (must be an odd number)\nIf the kernel size is less than 3, no smoothing is done.\nThe larger the number, the smoother the spectrum but beware of too much smoothing.")
        self.partially_incoherent_polynomial_sb = QSpinBox()
        self.partially_incoherent_polynomial_sb.setToolTip("Define the maximum degree of polynomial to be used in the smoothing filter")
        self.partially_incoherent_polynomial_sb.setRange(2,10)
        self.partially_incoherent_polynomial_sb.setSingleStep(1)
        self.partially_incoherent_polynomial_sb.setValue(self.settings["Filter polynomial size"])
        self.partially_incoherent_polynomial_sb.valueChanged.connect(self.on_partially_incoherent_polynomial_sb_changed)
        hbox.addWidget(self.partially_incoherent_kernel_sb)
        hbox.addWidget(self.partially_incoherent_polynomial_sb)
        label = QLabel("Smoothing filter (kernel & polynomial size)")
        label.setToolTip("Define the kernel size for the smoothing of incoherent spectra (must be an odd number)\nIf the kernel size is less than 3, no smoothing is done.\nThe larger the number, the smoother the spectrum but beware of too much smoothing.\nAlso defines the polynomial size for the fitting of the points in the kernel")
        return label,hbox

    def openDB_button_clicked(self):
        """Open a new materials' database.

        After the database is opened a refresh is requested.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        debugger.print("Start:: openDB_button_clicked")
        self.openDataBase()
        self.refresh(force=True)
        self.refreshRequired = True
        return

    def on_partially_incoherent_kernel_sb_changed(self,value):
        """Handle the change in the partially incoherent kernel spinbox value.

        Parameters
        ----------
        value : int or float
            The new value of the partially incoherent kernel spinbox.

        Returns
        -------
        None

        """        
        debugger.print(self.settings["Legend"],"on_partially_incoherent_kernel_sb_changed", value)
        self.refreshRequired = True
        self.settings["Filter kernel size"] = value
        return

    def on_partially_incoherent_polynomial_sb_changed(self,value):
        """Handle changes to the partially incoherent polynomial spinbox value.

        This method updates the relevant setting based on the new value and marks
        the need for refreshment. It also logs the change using the debugger.

        Parameters
        ----------
        value : int
            The new value for the partially incoherent polynomial spinbox.

        Returns
        -------
        None

        """        
        debugger.print(self.settings["Legend"],"on_partially_incoherent_polynomial_sb_changed", value)
        self.refreshRequired = True
        self.settings["Filter polynomial size"] = value
        return

    def on_partially_incoherent_samples_sb_changed(self,value):
        """Handle the change in partially incoherent samples setting.

        Parameters
        ----------
        value : int or float
            The new value for the partially incoherent samples setting.

        Returns
        -------
        None

        """        
        debugger.print(self.settings["Legend"],"on_partially_incoherent_samples_sb_changed", value)
        self.refreshRequired = True
        self.settings["Partially incoherent samples"] = value
        self.noCalculationsRequired = value
        return

    def on_percentage_partial_incoherence_sb_changed(self,value):
        """Handle changes in percentage partial incoherence spinner box value.

        This method updates the 'Percentage partial incoherence' setting based on the new value provided by the relevant spinner box. It marks that a refresh is required and logs the change using the configured debugger.

        Parameters
        ----------
        value : int or float
            The new value from the percentage partial incoherence spinner box.

        Returns
        -------
        None

        """        
        debugger.print(self.settings["Legend"],"on_percentage_partial_incoherence_sb_changed", value)
        self.refreshRequired = True
        self.settings["Percentage partial incoherence"] = value
        return

    def on_global_azimuthal_angle_sb_changed(self,value):
        """Handle the change in global azimuthal angle setting.

        This method updates the 'Global azimuthal angle' in the settings 
        based on the given value and marks a refresh as required.

        Parameters
        ----------
        value : float
            The new value for the global azimuthal angle to be set.

        Returns
        -------
        None

        """        
        debugger.print(self.settings["Legend"],"on_global_azimuthal_angl_sb_changed", value)
        self.refreshRequired = True
        self.settings["Global azimuthal angle"] = value
        return

    def on_angle_of_incidence_sb_changed(self,value):
        """Update the 'Angle of incidence' setting and mark a refresh as required upon a change in the angle of incidence spinbox value.

        A refresh is requested

        Parameters
        ----------
        value : int or float
            The new value of the angle of incidence from the spinbox.

        Returns
        -------
        None

        """        
        debugger.print(self.settings["Legend"],"on_angle_of_incidence_sb_changed", value)
        self.refreshRequired = True
        self.settings["Angle of incidence"] = value
        return

    def refresh(self,force=False):
        """Refresh the current state, optionally forcing a refresh.

        This method updates the current state based on various settings and conditions. It checks
        if a refresh is required or if a forceful refresh is requested. It then updates several
        settings, including material database, layer settings, calculation requirements, UI elements states,
        and more to reflect the current state. It also processes UI events to ensure the UI is responsive
        during the update.

        Parameters
        ----------
        force : bool, optional
            A flag to force the refresh even if it's not required. Default is False.

        Returns
        -------
        None

        """        
        debugger.print(self.settings["Legend"],"Start:: refresh, force =", force)
        if not self.refreshRequired and not force :
            debugger.print(self.settings["Legend"],"Finished:: refreshing widget aborted", self.refreshRequired,force)
            return
        # Check to see if there is a new reader, if there is set up the cell
        self.reader = self.notebook.reader
        if self.reader is None:
            return
        #
        # Determine the exponent_threshold to be used by GTM
        #
        if "Linux" in self.notebook.mainTab.settings["Compatibility mode"]:
            self.exponent_threshold = 11000
        elif "Windows" in self.notebook.mainTab.settings["Compatibility mode"]:
            self.exponent_threshold = 700
        else:
            self.exponent_threshold = 11000
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Initialise the cell
        self.cell = self.reader.get_unit_cell()
        # Open database and get materials
        if self.settings["Materials database"] != self.DataBase.getFileName():
            self.DataBase = MaterialsDataBase(self.settings["Materials database"],debug=debugger.state())
            self.settings["Materials database"] = self.DataBase.getFileName()
            self.database_le.setText(self.settings["Materials database"])
            # Update the possible  material names from the database
            self.materialNames = self.setMaterialNames()
        # Generate the layers from the settings
        self.settings2Layers()
        self.generateLayerSettings()
        # Force recalculation
        self.calculationRequired = True
        # Change any greyed out items
        self.greyed_out()
        #
        # Now refresh values that need updating
        #
        # Check the mode 
        index = self.mode_cb.findText(self.settings["Mode"], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        # Work out how many calculations are going to be performed
        self.set_noCalculationsRequired()
        # Update the Legend widget
        self.legend_le.setText(self.settings["Legend"])
        # Update angle widgets
        self.global_azimuthal_angle_sb.setValue(self.settings["Global azimuthal angle"])
        self.angle_of_incidence_sb.setValue(self.settings["Angle of incidence"])
        # For partial incoherent case, set percentage variation of angles and thickness and the number of samples
        self.percentage_partial_incoherence_sb.setValue(self.settings["Percentage partial incoherence"])
        self.partially_incoherent_samples_sb.setValue(self.settings["Partially incoherent samples"])
        # For partial incoherent case, set the smoothing parameters
        self.partially_incoherent_kernel_sb.setValue(self.settings["Filter kernel size"])
        self.partially_incoherent_polynomial_sb.setValue(self.settings["Filter polynomial size"])
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        # Redraw the layer information widget
        self.redrawLayerTable()
        for _i in range(20):
            QCoreApplication.processEvents()
        self.refreshRequired = False
        debugger.print(self.settings["Legend"],"Finished:: refresh, force =", force)
        return

    def set_noCalculationsRequired(self):
        """Determine the number of calculations required.

        The routine first works out the number of layers needing phase averaging,
        then it looks at the number of partial incoherent calculations being performed.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # First see how many layers are using phase averaging
        self.number_of_average_incoherent_layers = 0
        for layer in self.layers:
            if layer.getIncoherentOption() == "Incoherent (phase averaging)":
                self.number_of_average_incoherent_layers += 1
        # First see how many layers are using phase averaging
        if self.number_of_average_incoherent_layers > 0:
            number_of_samples = self.settings["Number of average incoherence samples"]
            self.noCalculationsRequired = pow(number_of_samples,self.number_of_average_incoherent_layers)
        else:
            self.noCalculationsRequired = 1
        # Now see if partial incoherence is being used
        if self.settings["Percentage partial incoherence"] > 0:
            self.noCalculationsRequired *= self.settings["Partially incoherent samples"]
        return

    def on_mode_cb_activated(self, index):
        """Handle the activation of the mode combo box in the GUI.

        Parameters
        ----------
        index : int
            The index of the selected mode in the combo box. This controls what mode the system will use. 
            Index 0 corresponds to the 'Transfer matrix' mode, while index 1 corresponds to the 'Scattering matrix' mode.

        Returns
        -------
        None

        Notes
        -----
        This function performs several operations based on the selected mode:
        - Sets the mode in the settings to either 'Transfer matrix' or 'Scattering matrix'.
        - If 'Scattering matrix' mode is selected (index 1), it goes through all layers and changes any layer with the incoherent option set to 'Incoherent (intensity)' to be 'Coherent'.
        - Invokes set_noCalculationsRequired to indicate that new calculations are necessary based on the changed settings.
        - Calls generateLayerSettings to update the layer settings according to the new mode.
        - Forces a refresh of the display to reflect any changes.
        - Sets a flag indicating that a refresh is required.

        """        
        debugger.print(self.settings["Legend"],"Start:: on_mode_cb_activated")
        if index == 0:
            self.settings["Mode"] = "Transfer matrix"
        elif index == 1:
            self.settings["Mode"] = "Scattering matrix"
            # If any layers are using intensity incoherence move them to coherent
            for layer in self.layers:
                incoherentOption = layer.getIncoherentOption()
                if incoherentOption == "Incoherent (intensity)":
                    layer.setIncoherentOption("Coherent")
        self.set_noCalculationsRequired()
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired = True
        debugger.print(self.settings["Legend"],"Mode changed to ", self.settings["Mode"])
        debugger.print(self.settings["Legend"],"Finished:: on_mode_cb_activated")
        return

    def average_incoherent_calculator( self,
                            layers,
                            mode,
                            theta,
                            phi,
                            psi,
                            angleOfIncidence,
                            exponent_threshold):
        """Calculate the incoherent component of light reflectance and transmission by averaging over the phase shift.

        Parameters
        ----------
        layers : a list of layers
            A list of layers
        mode : str
            Specifies either 'Transfer matrix' or 'Scattering matrix'
        theta : float
            The angle theta
        phi : float
            The angle phi
        psi : float
            The angle psi
        angleOfIncidence : float
            The angle incidence
        exponent_threshold: float
            The exponent threshold

        Returns
        -------
        A tuple of lists or nd.arrays
            - p_reflectance
            - s_reflectance
            - p_transmittance
            - s_transmittance
            - p_absorbtance
            - s_absorbtance
            - epsilon

        """
        debugger.print(self.settings["Legend"],"Start:: partially_incoherent_calculator")
        #
        # Zero the arrays we will need
        #
        size = len(self.vs_cm1)
        av_p_reflectance = np.zeros(size) 
        av_s_reflectance = np.zeros(size) 
        av_p_transmittance = np.zeros(size) 
        av_s_transmittance = np.zeros(size) 
        av_s_absorbtance = np.zeros(size) 
        av_p_absorbtance = np.zeros(size) 
        av_epsilon = np.zeros((size,3,3),dtype=np.cdouble)
        # Work out which of the layers is the crystal dielectric
        averageList = []
        for layer in layers:
            if layer.getIncoherentOption() == "Incoherent (phase averaging)":
                averageList.append(layer)
        #
        # Calculate the list of phase shift combinations from the number of samples
        #
        frac = self.settings["Percentage average incoherence"] / 100.0
        number_of_samples = self.settings["Number of average incoherence samples"]
        beta = [frac *2 * np.pi * s / number_of_samples for s in range(number_of_samples)]
        betas = product(beta, repeat=len(averageList))
        # loop over the combination of phases that are possible
        for beta in betas:
            # set the phase shift for each layer
            for index,layer in enumerate(averageList):
                layer.setPhaseShift(beta[index])
            ( p_reflectance, 
            s_reflectance, 
            p_transmittance, 
            s_transmittance, 
            s_absorbtance, 
            p_absorbtance, 
            epsilon) = self.coherent_calculator(
                                       layers,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence,
                                       exponent_threshold)
            av_p_reflectance   += np.array(p_reflectance)   / number_of_samples
            av_s_reflectance   += np.array(s_reflectance)   / number_of_samples
            av_p_transmittance += np.array(p_transmittance) / number_of_samples
            av_s_transmittance += np.array(s_transmittance) / number_of_samples
            av_s_absorbtance   += np.array(p_absorbtance)   / number_of_samples
            av_p_absorbtance   += np.array(s_absorbtance)   / number_of_samples
            av_epsilon         += np.array(epsilon)         / number_of_samples
        return (  av_p_reflectance.tolist(), av_s_reflectance.tolist(), av_p_transmittance.tolist(), av_s_transmittance.tolist(), av_p_absorbtance.tolist(), av_s_absorbtance.tolist(), av_epsilon.tolist() )

    def partially_incoherent_calculator( self,
                            layers,
                            mode,
                            theta,
                            phi,
                            psi,
                            angleOfIncidence,
                            exponent_threshold):
        """Calculate the incoherent component of light reflectance and transmission by sampling the path length in the incident medium.

        Parameters
        ----------
        layers : a list of layers
            A list of layers
        mode : str
            Specifies either 'Transfer matrix' or 'Scattering matrix'
        theta : float
            The angle theta
        phi : float
            The angle phi
        psi : float
            The angle psi
        angleOfIncidence : float
            The angle incidence
        exponent_threshold: float
            The exponent threshold

        Returns
        -------
        A tuple of lists or nd.arrays
            - p_reflectance
            - s_reflectance
            - p_transmittance
            - s_transmittance
            - p_absorbtance
            - s_absorbtance
            - epsilon

        """
        debugger.print(self.settings["Legend"],"Start:: partially_incoherent_calculator")
        #
        # Zero the arrays we will need
        #
        size = len(self.vs_cm1)
        av_p_reflectance = np.zeros(size) 
        av_s_reflectance = np.zeros(size) 
        av_p_transmittance = np.zeros(size) 
        av_s_transmittance = np.zeros(size) 
        av_s_absorbtance = np.zeros(size) 
        av_p_absorbtance = np.zeros(size) 
        av_epsilon = np.zeros((size,3,3),dtype=np.cdouble)
        # Work out which of the layers is the crystal dielectric
        crystalLayer = None
        for layer in layers:
            if layer.isDielectric():
                crystalLayer = layer
        #
        # Loop over the number of samples requred
        #
        crystalDepth = crystalLayer.getThickness()
        keepCrystalDepth = crystalDepth
        d = crystalDepth
        t = theta
        p  = phi
        q  = psi
        a  = angleOfIncidence
        fractionalIncoherence = self.settings["Percentage partial incoherence"]/100.0
        for _s in range(self.settings["Partially incoherent samples"]):
            crystalDepth = d + d*( -1 + 2*np.random.rand())*fractionalIncoherence
            crystalLayer.setThickness(crystalDepth)
            theta = t + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            phi   = p + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            psi = q + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            angleOfIncidence = a + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            ( p_reflectance, 
            s_reflectance, 
            p_transmittance, 
            s_transmittance, 
            s_absorbtance, 
            p_absorbtance, 
            epsilon) = self.coherent_calculator(
                                       layers,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence,
                                       exponent_threshold)
            av_p_reflectance   += np.array(p_reflectance) / self.settings["Partially incoherent samples"]
            av_s_reflectance   += np.array(s_reflectance) / self.settings["Partially incoherent samples"]
            av_p_transmittance += np.array(p_transmittance) / self.settings["Partially incoherent samples"]
            av_s_transmittance += np.array(s_transmittance) / self.settings["Partially incoherent samples"]
            av_s_absorbtance   += np.array(p_absorbtance) / self.settings["Partially incoherent samples"]
            av_p_absorbtance   += np.array(s_absorbtance) / self.settings["Partially incoherent samples"]
            av_epsilon         += np.array(epsilon) / self.settings["Partially incoherent samples"]
        crystalLayer.setThickness(keepCrystalDepth)
        # Only apply the smoothing filter if the kernel is larger than 2
        k = self.settings["Filter kernel size"]
        if k > 2:
            p = self.settings["Filter polynomial size"]
            av_p_reflectance   = signal.savgol_filter(av_p_reflectance, k, p, mode="nearest")
            av_p_transmittance = signal.savgol_filter(av_p_transmittance, k, p, mode="nearest")
            av_p_absorbtance   = signal.savgol_filter(av_p_absorbtance, k, p, mode="nearest")
            av_s_reflectance   = signal.savgol_filter(av_s_reflectance, k, p, mode="nearest")
            av_s_transmittance = signal.savgol_filter(av_s_transmittance, k, p, mode="nearest")
            av_s_absorbtance   = signal.savgol_filter(av_s_absorbtance, k, p, mode="nearest")
        return (  av_p_reflectance.tolist(), av_s_reflectance.tolist(), av_p_transmittance.tolist(), av_s_transmittance.tolist(), av_p_absorbtance.tolist(), av_s_absorbtance.tolist(), av_epsilon.tolist() )

    def coherent_calculator( self,
                            layers,
                            mode,
                            theta,
                            phi,
                            psi,
                            angleOfIncidence,
                            exponent_threshold):
        """Calculate the coherent component of light reflectance and transmission.

        Parameters
        ----------
        layers : a list of layers
            A list of layers
        mode : str
            Specifies either 'Transfer matrix' or 'Scattering matrix'
        theta : float
            The angle theta
        phi : float
            The angle phi
        psi : float
            The angle psi
        angleOfIncidence : float
            The angle incidence
        exponent_threshold: float
            The exponent threshold

        Returns
        -------
        A tuple of lists or nd.arrays
            - p_reflectance
            - s_reflectance
            - p_transmittance
            - s_transmittance
            - p_absorbtance
            - s_absorbtance
            - epsilon

        """
        debugger.print(self.settings["Legend"],"Entering the coherent_calculator function")
        #
        # Initialise the partial function to pass through to the pool
        #
        partial_function = partial(solve_single_crystal_equations,
                                       layers,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence,
                                       exponent_threshold)
        results = []
        # About to call
        debugger.print(self.settings["Legend"],"About to calculate single crystal scenario using pool")
        if self.notebook.pool is None:
            self.notebook.startPool()
        for result in self.notebook.pool.imap(partial_function, self.vs_cm1, chunksize=20):
            self.notebook.progressbars_update()
            results.append(result)
        QCoreApplication.processEvents()
        # Initialise result variables
        p_reflectance = []
        s_reflectance = []
        p_transmittance = []
        s_transmittance = []
        p_absorbtance = []
        s_absorbtance = []
        epsilon = []
        debugger.print(self.settings["Legend"],"About to extract results for single crystal scenario")
        for v,_r,R,_t,T,eps,errors,largest_exponent in results:
            if self.settings["Mode"] == "Transfer matrix" and errors > 0:
                print("Warning exponential overflow occured at frequency",v,errors,largest_exponent)
            p_reflectance.append(R[0]+R[2])
            s_reflectance.append(R[1]+R[3])
            p_transmittance.append(T[0])
            s_transmittance.append(T[1])
            p_absorbtance.append(1.0 - R[0]-R[2]-T[0])
            s_absorbtance.append(1.0 - R[1]-R[3]-T[1])
            epsilon.append(eps)
        debugger.print(self.settings["Legend"],"Finished the coherent_calculator function")
        return ( p_reflectance, s_reflectance, p_transmittance, s_transmittance, p_absorbtance, s_absorbtance, epsilon )

    def calculate(self,vs_cm1):
        """Perform simulation for calculating various properties such as reflectance, transmittance, and absorbance for a given set of material layers and configurations.

        Parameters
        ----------
        vs_cm1 : list or ndarray
            List or array of wavenumbers (in cm^-1) for which the simulation is to be run.

        Returns
        -------
        None

        Notes
        -----
        This method updates the object's state by calculating and setting various optical properties (reflectance, transmittance, absorbtance) based on the input wavenumbers, layer configurations, and other settings such as the angle of incidence. It requires that settings, frequencies, and other parameters be previously defined and valid. The method also makes use of different calculator objects depending on the coherence conditions and partial incoherence percentage. It handles conditions such as missing program, file reader, or filename settings by aborting the calculations and logging appropriate messages through a `debugger`.

        """        
        debugger.print(self.settings["Legend"],"Start:: calculate - number of frequencies",len(vs_cm1))
        if not self.calculationRequired:
            debugger.print(self.settings["Legend"],"Finished:: calculate aborted because calculationRequired false")
            return
        QCoreApplication.processEvents()
        self.vs_cm1 = vs_cm1
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings["Program"]
        filename = self.notebook.mainTab.getFullFileName()
        if self.reader is None:
            debugger.print(self.settings["Legend"],"Finished:: Calculate aborting - no reader")
            return
        if program == "":
            debugger.print(self.settings["Legend"],"Finished:: Calculate aborting - no program")
            return
        if filename == "":
            debugger.print(self.settings["Legend"],"Finished:: Calculate aborting - no file")
            return
        # Make sure the filter kernel size is odd, if not make it so and update the GUI
        if self.settings["Filter kernel size"] % 2 == 0:
            self.settings["Filter kernel size"] += 1
            self.partially_incoherent_kernel_sb.setValue(self.settings["Filter kernel size"])
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        # The euler angles are set to zero, apart from the global azimuthal angle
        # the rotation of each layer is now handled by the layer class.
        theta = 0.0
        phi = 0.0
        # global azimuthal angle in radians
        psi = np.radians(self.settings["Global azimuthal angle"])
        # Set the angle of incidence in radians
        angleOfIncidence = np.radians(self.settings["Angle of incidence"])
        # Tell each layer to calculate the euler matrix for rotation to the lab frame
        for layer in self.layers:
            if layer.isTensor():
                hkl = layer.getHKL()
                if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
                    QMessageBox.about(self,"",f"Unable to calculate surface for scenario {self.settings['Legend']}, hkl=[0,0,0]")
                    return
                layer.calculate_euler_matrix()
        # Define the mode of calculation transfer or scattering matrix
        mode = self.settings["Mode"]
        # See if the partially incoherent method will be used or the average incoherent
        # A different calculator is selected depending on the settings
        if self.settings["Percentage partial incoherence"] > 0:
            calculator = self.partially_incoherent_calculator
        elif self.number_of_average_incoherent_layers > 0:
            calculator = self.average_incoherent_calculator
        else:
            calculator = self.coherent_calculator
        # Pass the exponent threshold
        exponent_threshold = self.exponent_threshold
        # Call the relevant calculator
        ( self.p_reflectance, 
        self.s_reflectance, 
        self.p_transmittance, 
        self.s_transmittance, 
        self.p_absorbtance, 
        self.s_absorbtance, 
        self.epsilon) = calculator( self.layers,
                                    mode,
                                    theta,
                                    phi,
                                    psi,
                                    angleOfIncidence,
                                    exponent_threshold)
        debugger.print(self.settings["Legend"],"Finished:: calculate - number of frequencies",len(vs_cm1))
        return

    def get_result(self, vs_cm1, plot_type):
        """Return a particular result.

        Parameters
        ----------
        vs_cm1 : a list of floats
            The frequencies
        plot_type : str
            The data needed for plotting one of:
            - 'Crystal Reflectance (P polarisation)'
            - 'Crystal Reflectance (S polarisation)'
            - 'Crystal Transmittance (P polarisation)'
            - 'Crystal Transmittance (S polarisation)'
            - 'Crystal Absorbtance (P polarisation)'

        Returns
        -------
        list of floats
            The results to be plotted

        """
        debugger.print(self.settings["Legend"],"Start:: get_result",len(vs_cm1),plot_type)
        self.get_results(vs_cm1)
        debugger.print(self.settings["Legend"],"Finished:: get_result",len(vs_cm1),plot_type)
        return{ "Crystal Reflectance (P polarisation)"  : self.p_reflectance,
                "Crystal Reflectance (S polarisation)"  : self.s_reflectance,
                "Crystal Transmittance (P polarisation)": self.p_transmittance,
                "Crystal Transmittance (S polarisation)": self.s_transmittance,
                "Crystal Absorbtance (P polarisation)"  : self.p_absorbtance,
                "Crystal Absorbtance (S polarisation)"  : self.s_absorbtance,
        }.get(plot_type)

    def get_results(self, vs_cm1):
        """Return the results of the single crystal calculation.

        If a refresh has been requested it is performed before the calculate method is called.

        Parameters
        ----------
        vs_cm1 : a list of floats
            The frequencies in cm-1

        Returns
        -------
        None

        """
        debugger.print(self.settings["Legend"],"Start:: get_results",len(vs_cm1),self.refreshRequired)
        if len(vs_cm1)>0 and (self.refreshRequired or len(self.vs_cm1) != len(vs_cm1) or self.vs_cm1[0] != vs_cm1[0] or self.vs_cm1[1] != vs_cm1[1]) :
            debugger.print(self.settings["Legend"],"get_results recalculating")
            self.refresh()
            self.calculate(vs_cm1)
        else:
            debugger.print(self.settings["Legend"],"get_results no need for recalculation")
            #self.notebook.progressbars_update(increment=len(vs_cm1))
        debugger.print(self.settings["Legend"],"Finished:: get_results",len(vs_cm1),self.refreshRequired)
        return

    def greyed_out(self):
        """Have a look through the settings and see if we need to grey anything out.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # At the moment it appears there is nothing to do.
        debugger.print(self.settings["Legend"],"Finished:: greyed_out")
        return

