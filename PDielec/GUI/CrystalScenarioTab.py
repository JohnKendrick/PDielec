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
"""CrystalScenarioTab module."""
import copy
import logging
import math
from functools import partial
from itertools import product

import numpy as np
from qtpy.QtCore import QCoreApplication, QSize, Qt
from qtpy.QtWidgets import (
    QApplication,
    QCheckBox,
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
from PDielec.Constants import amu, angs2bohr, speed_light_si, wavenumber
from PDielec.GUI.ScenarioTab import ScenarioTab
from PDielec.GUI.SingleCrystalLayer import ShowLayerWindow, SingleCrystalLayer
from PDielec.LayeredRamanCalculator import LayeredRamanCalculator, RamanLayer, lorentzian_broaden
from PDielec.Materials import MaterialsDataBase

logger = logging.getLogger(__name__)
thickness_conversion_factors = {"ang":1.0E-10, "nm":1.0E-9, "um":1.0E-6, "mm":1.0E-3, "cm":1.0E-2}
thickness_units = list(thickness_conversion_factors.keys())
# incoherent_options = ['Coherent','Incoherent (intensity)','Incoherent (phase cancelling)','Incoherent (phase averaging)','Incoherent (non-reflective)'] 
incoherent_options = ["Coherent","Incoherent (intensity)","Incoherent (phase averaging)","Incoherent (non-reflective)"] 
gtm_methods ={"Coherent":GTM.CoherentLayer,
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
    gtm_layers = []
    # Create layers from all the layers between first and last
    for layer in selectedLayers:
        incoherentOption = layer.get_incoherent_option()
        gtm_layers.append(gtm_methods[incoherentOption](layer, exponent_threshold=exponent_threshold))
    # Creat the system with the layers 
    if mode == "Scattering matrix":
        system = GTM.ScatteringMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtm_layers)
    else:
        system = GTM.TransferMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtm_layers)
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
    errors,largest_exponent = system.overflow_errors()
    return v,r,R,t,T,epsilon,errors,largest_exponent

class CrystalScenarioTab(ScenarioTab):
    """A tab class for handling a crystal infrared  scenario.

    This class extends a generic scenario tab (:class:`~PDielec.GUI.ScenarioTab.ScenarioTab`) with functionalities
    specific to handling crystal infrared scenarios, including setting up the scenario environment, defining the crystal
    orientation, layer information, and calculation modes (e.g., using Transfer Matrix or Scattering Matrix methods). It
    supports incoherence effects through partial incoherence simulations and provides tools for fine-tuning simulation
    parameters such as incidence angles, layer materials, thicknesses, and orientations.

    Parameters
    ----------
    parent : QWidget
        The parent widget (typically a main window or a larger container) to which this tab will belong.
    debug : bool, optional
        A flag to indicate whether debugging information should be printed to the console. Default is False.

    Attributes
    ----------
    refresh_required : bool
        Flag indicating whether the widget's data needs to be refreshed.
    calculation_required : bool
        Flag indicating whether a new calculation is required based on changes in parameters or settings.
    refresh_required : bool
       Indicates whether the scenario settings have been changed and thus require the scenario to be redrawn or recalculated.
    no_calculations_required : int
       The number of calculations required for the simulation.
       This value may change depending on the specifics of the scenario
       configuration (e.g., the inclusion of incoherent effects requiring multiple sample calculations).
    settings : dict
       A dictionary holding various scenario parameters and their values,
       which include settings for the type of analysis,
       material layers, angles, and method-specific parameters.
    material_names : list
       A list of names for materials that can be used in the scenario layers.
    p_reflectance, s_reflectance, p_transmittance, s_transmittance, p_absorbtance, s_absorbtance : list
       Lists holding computed values for reflectance, transmittance, and absorbance for both P and S polarization modes.
    epsilon : list
       A list of dielectric function values computed during simulation.
    layers : list
       A list of layer configurations used in the simulation,
       each containing information about the material, thickness, orientation, and other layer-specific properties.


    Methods
    -------
    There are several methods within the class for handling events (e.g., button clicks, combo box activations),
    performing calculations, redrawing tables, and managing layer settings. These include methods for adding, deleting,
    or altering layers; opening material databases; changing simulation settings; calculating and displaying results;
    and initializing the user interface components relevant to the crystal infrared scenario.

    angle_of_incidence_widget
        Create a widget to set the angle of incidence
    average_incoherent_calculator
        Calculates the incoherent light reflectance and transmittance
    calculate
        Perform the calculations as define in the GUI
    coherent_calculator
        Calculate the coherent light reflectance and transmittance
    create_tool_bar
         Create the tool bar used for the material layer
    create_tool_bar_delete_button
         Create the delete button in the tool bar
    create_tool_bar_move_down_button
         Create the move down button in the tool bar
    create_tool_bar_move_up_button
         Create the move up button in the tool bar
    delete_layer
         Delete a layer from the list of layers
    draw_layer_table
         Draw the layer table widget, each layer has a line in the table
    generate_layer_settings
         Generate the settings dictionary with an item for each layer
    get_dielectric_layer_index
         Get the index of the dielectric layer
    get_material_from_data_base
         Get a material from the database, or define the material from the given permittivity
    get_result
         Return the request information
    get_results
         Calculate all the reflectance, transmittance and absorptance information
    global_azimuthal_widget
         Set the global azimuthal angle
    move_layer_down
         Move a layer up the layer table
    move_layer_up
         Move a layer down the layer table
    new_layer_widget
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
    on_new_layer_cb_activated
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
    open_db_button_clicked
         Handle a click on the opend database button
    partial_incoherence_widget
         Create a partialIncoherence widget
    partially_incoherent_calculator
         Calculate incoherent light scattering using a partially incoherent method
    print_layer_settings
         Print out the layer settings dictionary
    redraw_layer_table
         Redraw the layer table
    redraw_layer_table_row
         Redraw a row of the layer table
    refresh
         Perform a refresh of the GUI
    set_material_names
         Read the material names from the database
    set_no_calculations_required
         Calculate the number of calculations required
    settings2Layers
         Read the settings dictionary and create the necessary layers
    smoothing_widget
         A widget to handle the smoothing information

    """

    def __init__(self, parent, spectroscopy="Crystal Infrared", debug=False ):
        """Initialize the crystal infrared Scenario Tab.

        This initializer sets up the GUI components, populates settings with default values, 
        and prepares the tab for interaction. It inherits from a parent scenario, initializes
        a debug mode if required, and sets various default parameters for the crystal infrared 
        scenario, including layer information and calculation modes.

        Parameters
        ----------
        parent : QWidget
            The parent widget or window, typically the main application window or a central widget
            that this tab will be a part of.
        spectroscopy : str
            Spectroscopy can be one of "Crystal Infrared" or "Crystal Raman"
        debug : bool, optional
            A boolean flag to indicate whether debug messages should be printed to the console
            or log. Defaults to False if not specified.

        Attributes
        ----------
        refresh_required : bool
            Indicates if the view needs refreshing.
        calculation_required : bool
            Indicates if a new calculation is required based on changed parameters.
        no_calculations_required : int
            Number of calculations required, initialized to 1.
        settings : dict
            A dictionary of settings for the simulation, initialized with default values.
        material_names : list
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
        logger.debug("Start:: initialiser")
        if "Crystal" not in spectroscopy:
            logger.error(f"PowderSenarioTab failed incompatible spectroscopy {spectroscopy}")
        self.refresh_required = True
        self.calculation_required = True
        self.no_calculations_required = 1
        self.spectroscopy = spectroscopy
        self.settings["Scenario type"] = self.spectroscopy
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
        self.material_names = []
        self.p_reflectance = []
        self.s_reflectance = []
        self.p_transmittance = []
        self.s_transmittance = []
        self.p_absorbtance = []
        self.s_absorbtance = []
        self.epsilon = []
        self.raman_spectrum = []
        self.raman_mode_frequencies = np.array([])
        self.raman_mode_intensities = np.array([])
        self.raman_mode_sigmas = np.array([])
        self.layers = []
        self.settings["Laser wavelength nm"] = 532.0
        self.settings["Incident polarisation"] = "p"
        self.settings["Detected polarisation"] = "unpolarised"
        self.settings["Temperature K"] = 298.0
        self.settings["GL point density"] = 20.0
        self.settings["Collection side"] = "superstrate"  # 'superstrate' = backscatter, 'substrate' = forward
        self.settings["Collection angle"] = -1.0          # negative sentinel: default to angle of incidence
        self.settings["Coherent layer summation"] = False
        self.settings["Approximate ES"] = False
        self.settings["Phonon boundary correction"] = "none"  # 'none', 'NAC', or 'slab-environment'
        self.settings["Azimuthal sweep points"] = 36
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
        self.DataBase = MaterialsDataBase(self.settings["Materials database"])
        self.settings["Materials database"] = self.DataBase.get_file_name()
        self.material_names = self.set_material_names()
        # Create the layers - superstrate / dielectric / substrate from the defaults layer settings
        if self.reader is not None:
            self.settings2Layers()
        # Create last tab - CrystalInfraredScenarioTab
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
        label.setToolTip("Set the method for calculating light transmission and reflectanceb;\n Transfer matrix.  This method is fast but can be numerically unstable.\n Scattering matrix. This method is slower but is numerically stable")
        self.form.addRow(label, self.mode_cb)
        #
        # Define the global azimuthal angle widget
        #
        label,layout = self.global_azimuthal_widget()
        self.form.addRow(label, layout)
        #
        # Define the angle of incidence widget
        #
        label, layout = self.angle_of_incidence_widget()
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
        self.form.addRow(self.draw_layer_table())
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
        label,layout = self.partial_incoherence_widget()
        #
        # Smoothing widget
        #
        label,layout = self.smoothing_widget()
        #
        # Crystal Raman specific widgets (Phase 2f + Phase 3a)
        #
        if spectroscopy == "Crystal Raman":
            self._build_raman_widgets()
        #
        # Add a legend option
        #
        self.legend_le = QLineEdit(self)
        self.legend_le.setToolTip("The legend will be used to describe the results in the plot")
        self.legend_le.setText(self.settings["Legend"])
        self.legend_le.textChanged.connect(self.on_legend_le_changed)
        legend_label_text = "Crystal Raman Scenario legend" if spectroscopy == "Crystal Raman" else "Crystal IR Scenario legend"
        label = QLabel(legend_label_text)
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
        logger.debug("Finished:: initialiser")

    def redraw_layer_table(self):
        """Redraw the layer table widget.

        Adds a row for each layer in the gadget.
        The row is drawn using redraw_layer_table_row()
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
            self.redraw_layer_table_row(sequenceNumber,layer,rowCount,firstLayer,lastLayer)
        # Add a 'create new layer' button
        rowCount += 1
        new_layer_cb = self.new_layer_widget()
        new_layer_cb.setStyleSheet("Text-align:left")
        self.layerTable_tw.setRowCount(rowCount)
        self.layerTable_tw.setCellWidget(rowCount-1,0,new_layer_cb)

    def redraw_layer_table_row(self,sequenceNumber,layer,rowCount,firstLayer,lastLayer):
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
        material = layer.get_material()
        material_name = material.get_name()
        layer_button = QPushButton(material_name)
        layer_button.setToolTip("Show the material properties in a new window")
        layer_button.setStyleSheet("Text-align:left")
        layer_button.clicked.connect(lambda x: self.on_layer_button_clicked(x,layer,sequenceNumber))
        self.layerTable_tw.setCellWidget(sequenceNumber,0,layer_button)
        # Handle thickness 
        material_thickness = layer.get_thickness()
        thickness_unit = layer.get_thickness_unit()
        film_thickness_sb = QDoubleSpinBox(self)
        film_thickness_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        film_thickness_sb.setToolTip("Define the thin film thickness in the defined thickness units")
        film_thickness_sb.setRange(0,100000)
        film_thickness_sb.setDecimals(3)
        film_thickness_sb.setSingleStep(0.001)
        film_thickness_sb.setValue(material_thickness)
        film_thickness_sb.valueChanged.connect(lambda x: self.on_film_thickness_sb_changed(x,layer))
        self.layerTable_tw.setCellWidget(sequenceNumber,1,film_thickness_sb)
        # thickness unit
        thickness_unit_cb = QComboBox(self)
        thickness_unit_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thickness_unit_cb.setToolTip("Set the units to be used for thickness; either angs nm, um, mm or cm")
        thickness_unit_cb.addItems( thickness_units )
        index = thickness_unit_cb.findText(thickness_unit, Qt.MatchFixedString)
        thickness_unit_cb.setCurrentIndex(index)
        thickness_unit_cb.activated.connect(lambda x: self.on_thickness_units_cb_activated(x, layer))
        self.layerTable_tw.setCellWidget(sequenceNumber,2,thickness_unit_cb)
        # define hkl
        h_sb = QSpinBox(self)
        h_sb.setToolTip("Define the h dimension of the unique direction")
        h_sb.setRange(-20,20)
        h_sb.setSingleStep(1)
        h_sb.setValue(layer.get_hkl()[0])
        h_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,0,layer))
        h_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        k_sb = QSpinBox(self)
        k_sb.setToolTip("Define the k dimension of the unique direction")
        k_sb.setRange(-20,20)
        k_sb.setSingleStep(1)
        k_sb.setValue(layer.get_hkl()[1])
        k_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,1,layer))
        k_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        l_sb = QSpinBox(self)
        l_sb.setToolTip("Define the l dimension of the unique direction")
        l_sb.setRange(-20,20)
        l_sb.setSingleStep(1)
        l_sb.setValue(layer.get_hkl()[2])
        l_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,2,layer))
        l_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.layerTable_tw.setCellWidget(sequenceNumber,3,h_sb)
        self.layerTable_tw.setCellWidget(sequenceNumber,4,k_sb)
        self.layerTable_tw.setCellWidget(sequenceNumber,5,l_sb)
        # define azimuthal angle
        azimuthal = layer.get_azimuthal()
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
        option_cb.addItems( incoherent_options )
        # We can't use incoherent intensity method with the scattering matrix method
        if self.settings["Mode"] == "Scattering matrix" and layer.get_incoherent_option() == "Incoherent (intensity)":
            layer.set_incoherent_option("Coherent")
        index = option_cb.findText(layer.get_incoherent_option(), Qt.MatchFixedString)
        option_cb.setCurrentIndex(index)
        option_cb.activated.connect(lambda x: self.on_option_cb_activated(x,layer))
        option_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        # Disable the intensity incoherence option if it is using scattering
        if self.settings["Mode"] == "Scattering matrix":
            index = option_cb.findText("Incoherent (intensity)", Qt.MatchFixedString)
            option_cb.model().item(index).setEnabled(False)
        self.layerTable_tw.setCellWidget(sequenceNumber,7,option_cb)
        # Create a toolbar for up down delete
        toolbar = self.create_tool_bar(layer,sequenceNumber,len(self.layers))
        self.layerTable_tw.setCellWidget(sequenceNumber,8,toolbar)
        # Add a Print option if debug is on
        if logger.isEnabledFor(logging.DEBUG):
            printButton = QPushButton("Print")
            printButton.setToolTip("Print the permittivity")
            printButton.clicked.connect(lambda x: self.on_print_button_clicked(x,layer))
            self.layerTable_tw.setCellWidget(sequenceNumber,9,printButton)
        if layer.is_scalar():
            h_sb.setEnabled(False)
            k_sb.setEnabled(False)
            l_sb.setEnabled(False)
            azimuthal_angle_sb.setEnabled(False)
        if firstLayer or lastLayer:
            film_thickness_sb.setEnabled(False)
            thickness_unit_cb.setEnabled(False)
        return

    def draw_layer_table(self):
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
        if logger.isEnabledFor(logging.DEBUG):
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
            self.redraw_layer_table()
        return self.layerTable_tw

    def delete_layer(self,x,layer,layerIndex):
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
        if layerIndex == 0 and layer[new].is_tensor():
            #  Only allow scalar materials as the superstrate
            logger.error("New superstrate material must be a scalar dielectric")
            return
        # Delete the layer
        del self.layers[layerIndex]
        self.generate_layer_settings()
        self.refresh(force=True)
        self.refresh_required=True
        return

    def move_layer_up(self,x,layer,layerIndex):
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
        if layerIndex == 1 and layer.is_tensor():
            #  Only allow scalar materials as the superstrate
            logger.error("New superstrate material must be a scalar dielectric")
            return
        new = layerIndex - 1
        item = self.layers[layerIndex]
        self.layers.pop(layerIndex)
        self.layers.insert(new, item)
        self.generate_layer_settings()
        self.refresh(force=True)
        self.refresh_required=True
        return

    def move_layer_down(self,x,layer,layerIndex):
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
        if layerIndex == 0 and self.layers[new].is_tensor():
            #  Only allow scalar materials as the superstrate
            logger.error("New superstrate material must be a scalar dielectric")
            return
        item = self.layers[layerIndex]
        self.layers.pop(layerIndex)
        self.layers.insert(new, item)
        self.generate_layer_settings()
        self.refresh(force=True)
        self.refresh_required=True
        return

    def create_tool_bar_move_up_button(self,layer,layerIndex,nLayers):
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
        moveUpButton.clicked.connect(lambda x: self.move_layer_up(x,layer,layerIndex))
        moveUpButton.setFixedSize(20,20)
        moveUpButton.setIconSize(QSize(20,20))
        moveUpButton.setStyleSheet("border: none;")
        moveUpButton.setToolTip("Move this layer up the list of layers")
        return moveUpButton

    def create_tool_bar_move_down_button(self,layer,layerIndex,nLayers):
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
        moveDownButton.clicked.connect(lambda x: self.move_layer_down(x,layer,layerIndex))
        moveDownButton.setFixedSize(20,20)
        moveDownButton.setIconSize(QSize(20,20))
        moveDownButton.setStyleSheet("border: none;")
        moveDownButton.setToolTip("Move this layer down the list of layers")
        return moveDownButton

    def create_tool_bar_delete_button(self,layer,layerIndex,nLayers):
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
        deleteButton.clicked.connect(lambda x: self.delete_layer(x,layer,layerIndex))
        deleteButton.setFixedSize(20,20)
        deleteButton.setIconSize(QSize(20,20))
        deleteButton.setStyleSheet("border: none;")
        deleteButton.setToolTip("Delete this layer")
        return deleteButton

    def create_tool_bar(self,layer,layerIndex,nLayers):
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
        moveUpButton   =  self.create_tool_bar_move_up_button(layer,layerIndex,nLayers)
        moveDownButton =  self.create_tool_bar_move_down_button(layer,layerIndex,nLayers)
        deleteButton   =  self.create_tool_bar_delete_button(layer,layerIndex,nLayers)
        nextIndex = layerIndex+1
        # disable any buttons that are irrelevant to the layer
        if layerIndex == 0:
            moveUpButton.setEnabled(False)
        if layerIndex == nLayers-1:
            moveDownButton.setEnabled(False)
        if layerIndex == 0 and self.layers[nextIndex].is_tensor():
            moveDownButton.setEnabled(False)
            deleteButton.setEnabled(False)
        if layerIndex == 1 and layer.is_tensor():
            moveUpButton.setEnabled(False)
        # Add the buttons to the frame and return the frame
        frame_layout.addWidget(moveUpButton)
        frame_layout.addWidget(moveDownButton)
        frame_layout.addWidget(deleteButton)
        return frame
        
    def new_layer_widget(self):
        """Create and return a new layer widget as a QComboBox.

        Returns a QComboBox widget configured for creating a new layer in a graphical user interface.
        The combo box is initialized with a default item for creating a new layer, followed by any existing layer names.

        Parameters
        ----------
        None

        Returns
        -------
        QComboBox A QComboBox widget configured with an initial item for creating a new layer and any additional
        existing material layer names.

        Notes
        -----
        - The first item in the combo box is a placeholder for creating a new layer, indicated by "New layer...". - The
          combo box is connected to the `on_new_layer_cb_activated` method, which should handle the action performed upon
          selecting an item from the combo box.

        """        
        new_layer_cb = QComboBox()
        new_layer_cb.setToolTip("Create a new layer")
        material_names = ["New layer..."]
        material_names += self.material_names
        new_layer_cb.addItems(material_names)
        new_layer_cb.setCurrentIndex(0)
        new_layer_cb.activated.connect(self.on_new_layer_cb_activated)
        return new_layer_cb

    def on_layerTable_itemChanged(self,item):
        """Handle a change to the layer table.

        An empty routine as changes are handled by other routines

        Parameters
        ----------
        item : tem
            The item which has changed

        Returns
        -------
        None

        """
        logger.debug(f"on_layerTable_itemChanged: {item}")
        return

    def on_new_layer_cb_activated(self,index):
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
        new_material_name = self.material_names[index-1]
        if "manual" in new_material_name:
            return
        newMaterial = self.get_material_from_data_base(new_material_name)
        hkl = [0,0,0]
        if newMaterial.is_tensor():
            hkl = [0,0,1]
        new_layer = SingleCrystalLayer(newMaterial,hkl=hkl,azimuthal=0.0,thickness=1.0,thickness_unit="um")
        self.layers.append(new_layer)
        self.generate_layer_settings()
        self.refresh(force=True)
        self.refresh_required = True
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
        material = layer.get_material()
        permittivity_object = material.get_permittivity_object()
        name = material.get_name()
        name = name.replace(" ","_")
        name += "_permittivity.csv"
        logger.info(f"Printing permittivity information to {name}")
        permittivity_object.print(0.0,2000.0,1.0,file=name)
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
        logger.debug(f"on_incoherence_cb_activated {index} {layer.get_name()}")
        option = incoherent_options[index]
        layer.set_incoherent_option(option)
        self.set_no_calculations_required()
        self.generate_layer_settings()
        self.refresh_required = True
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
        logger.debug(f"on_film_thickness_sb_changed {value} {layer.get_name()}")
        layer.set_thickness(value)
        self.generate_layer_settings()
        self.refresh_required = True
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
        logger.debug(f"Start:: on_thickness_units_cb_activated {index} {layer.get_name()}")
        unit = thickness_units[index]
        layer.set_thickness_unit(unit)
        self.generate_layer_settings()
        self.refresh_required = True
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
        logger.debug(f"on_azimuthal_angl_sb_changed {value} {layer.get_name()}")
        layer.set_azimuthal(value)
        layer.change_lab_frame_info()
        self.generate_layer_settings()
        self.refresh_required = True
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
        logger.debug(f"on_hkl_sb_changed {value}")
        hkl = layer.get_hkl()
        hkl[hkorl] = value
        layer.set_hkl(hkl)
        layer.change_lab_frame_info()
        self.generate_layer_settings()
        self.refresh_required = True
        return

    def generate_layer_settings(self):
        """Generate the settings dictionary for every layer.

        Parameters
        ----------
        None

        Returns
        -------
        list
            A list of dictionary settings for each layer.

        """
        logger.debug(f"{self.settings['Legend']} generate_layer_settings")
        self.settings["Layer material names"]  = []
        self.settings["Layer hkls"]            = []
        self.settings["Layer azimuthals"]      = []
        self.settings["Layer thicknesses"]     = []
        self.settings["Layer thickness units"] = []
        self.settings["Layer dielectric flags"] = []
        self.settings["Layer incoherent options"] = []
        for layer in self.layers:
            self.settings["Layer material names"].append(layer.get_material().get_name())
            self.settings["Layer hkls"].append(layer.get_hkl())
            self.settings["Layer azimuthals"].append(layer.get_azimuthal())
            self.settings["Layer thicknesses"].append(layer.get_thickness())
            self.settings["Layer thickness units"].append(layer.get_thickness_unit())
            self.settings["Layer dielectric flags"].append(layer.is_dielectric())
            self.settings["Layer incoherent options"].append(layer.get_incoherent_option())
        return

    def set_material_names(self):
        """Set and append a material name to the list of material names.

        This function retrieves the current list of material names from the database using
        `self.DataBase.get_sheet_names()`, appends 'Dielectric layer' to the list, and returns the updated list.

        Parameters
        ----------
        None

        Returns
        -------
        list
            The list of material names including the appended default 'Dielectric layer'.

        """        
        # Get the list of material names from the database
        material_names = self.DataBase.get_sheet_names()
        material_names.append("Dielectric layer")
        return material_names

    def print_layer_settings(self,message):
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
        This function prints the settings associated with layers, including material names, HKL values (Miller indices
        for the crystallographic planes), azimuthal angles, thicknesses with their units, flags indicating whether the
        layer is dielectric, and options for incoherent scattering, if applicable.

        This function does not return any value. It is primarily used for debugging or informative purposes to quickly
        display the current configuration of layer settings within a class that has `settings` as an attribute storing
        these properties.

        """        
        logger.debug(message)
        logger.debug(f"Layer material names: {self.settings['Layer material names']}")
        logger.debug(f"Layer hkls: {self.settings['Layer hkls']}")
        logger.debug(f"Layer azimuthals: {self.settings['Layer azimuthals']}")
        logger.debug(f"Layer thicknesses: {self.settings['Layer thicknesses']}")
        logger.debug(f"Layer thickness units: {self.settings['Layer thickness units']}")
        logger.debug(f"Layer dielectric flags: {self.settings['Layer dielectric flags']}")
        logger.debug(f"Layer incoherent options: {self.settings['Layer incoherent options']}")

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
        logger.debug(f"{self.settings['Legend']} settings2Layers")
        self.layers = []
        self.material_names = self.set_material_names()
        # Process the settings information and append each layer to the list
        for  name, hkl, azimuthal, thickness, thickness_unit, dielectricFlag, incoherentOption in zip(
                          self.settings["Layer material names"],
                          self.settings["Layer hkls"],
                          self.settings["Layer azimuthals"],
                          self.settings["Layer thicknesses"],
                          self.settings["Layer thickness units"],
                          self.settings["Layer dielectric flags"],
                          self.settings["Layer incoherent options"]):
            if name not in self.material_names:
                logger.error(f"Error material {name} not available {self.material_names}")
                name = "air"
            material = self.get_material_from_data_base(name)
            self.layers.append(SingleCrystalLayer(material,hkl=hkl,azimuthal=azimuthal,
                                     thickness=thickness,thickness_unit=thickness_unit,
                                     incoherentOption=incoherentOption,dielectricFlag=dielectricFlag))
        return

    def get_material_from_data_base(self,name,permittivity=None):
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
            crystalPermittivityObject = self.notebook.settingsTab.get_crystal_permittivity_object()
            material = Materials.External(name,permittivity_object=crystalPermittivityObject,cell=self.cell)
        elif name == "Material defined manually":
            material = Materials.Constant("Material defined manually",permittivity=permittivity)
        else:
            # Get the material from the data base
            # set the units for frequency to Hz for all materials
            material = self.DataBase.get_material(name)
        return material

    def get_dielectric_layer_index(self):
        """Return the index of the dielectric layer in the list of layers.

        - Sorts through the list of layers and reports the index of the first layer that came from the provided DFT
          calculation - If there is no such layer then the value None is returned

        Parameters
        ----------
        None

        Returns
        -------
        int or None
            The index of the dielectric layer in the list.

        """
        for index,layer in enumerate(self.layers):
            if layer.is_dielectric():
                self.dielectricLayer = layer
                return index
        return None

    def on_layer_button_clicked(self,x,layer,layerIndex):
        """Handle a click on the show layer widget.

        A new window is shown with additional information about material in the layer. For materials with a 3x3 tensor
        permittivity the window shows the relationship between the laboratory and crystal coordinates.

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
        show_layer_window = ShowLayerWindow(copy.copy(layer),message=message)
        if show_layer_window.exec():
            # The 'Ok' button was pressed
            # get the new Layer and replace the old one
            self.layers[layerIndex] = show_layer_window.get_layer()
            self.generate_layer_settings()
            self.refresh(force=True)
            self.refresh_required = True
        return

    def global_azimuthal_widget(self):
        """Create a global azimuthal angle widget.

        This widget is shown in the main crystal infrared scenario tab.
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

    def angle_of_incidence_widget(self):
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

    def _build_raman_widgets(self):
        """Build and add Crystal Raman-specific GUI widgets to self.form.

        Called from __init__ only when spectroscopy == 'Crystal Raman'.
        Creates widgets for laser frequency, temperature, GL quadrature points,
        incident/detected polarisation, collection geometry, coherent summation,
        approximate-E_S flag, and phonon boundary correction (Phase 3a).

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # Separator
        label = QLabel("Crystal Raman settings")
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        hbox = QHBoxLayout()
        hbox.addWidget(line)
        hbox.setAlignment(Qt.AlignVCenter)
        self.form.addRow(label, hbox)

        # Laser wavelength
        self.laser_wavelength_sb = QDoubleSpinBox(self)
        self.laser_wavelength_sb.setRange(200.0, 2000.0)
        self.laser_wavelength_sb.setSingleStep(1.0)
        self.laser_wavelength_sb.setDecimals(1)
        self.laser_wavelength_sb.setValue(self.settings["Laser wavelength nm"])
        self.laser_wavelength_sb.valueChanged.connect(self.on_laser_wavelength_sb_changed)
        self.laser_wavelength_sb.setToolTip("Laser wavelength in nm (e.g. 532, 633, 785)")
        label = QLabel("Laser wavelength (nm)")
        label.setToolTip(self.laser_wavelength_sb.toolTip())
        self.form.addRow(label, self.laser_wavelength_sb)

        # Temperature
        self.temperature_sb = QDoubleSpinBox(self)
        self.temperature_sb.setRange(0.0, 2000.0)
        self.temperature_sb.setSingleStep(1.0)
        self.temperature_sb.setDecimals(1)
        self.temperature_sb.setValue(self.settings["Temperature K"])
        self.temperature_sb.valueChanged.connect(self.on_temperature_sb_changed)
        self.temperature_sb.setToolTip("Sample temperature in Kelvin for the Bose-Einstein thermal prefactor")
        label = QLabel("Temperature (K)")
        label.setToolTip(self.temperature_sb.toolTip())
        self.form.addRow(label, self.temperature_sb)

        # GL quadrature point density (points per micron)
        self.gl_density_sb = QDoubleSpinBox(self)
        self.gl_density_sb.setRange(0.001, 1000.0)
        self.gl_density_sb.setSingleStep(1.0)
        self.gl_density_sb.setDecimals(3)
        self.gl_density_sb.setValue(self.settings["GL point density"])
        self.gl_density_sb.valueChanged.connect(self.on_gl_density_sb_changed)
        self.gl_density_sb.setToolTip(
            "Gauss-Legendre quadrature point density (points per µm of layer thickness).\n"
            "20–60 pts/µm is recommended for accuracy.\n"
            "The actual number of points is density × layer thickness, clamped to [5, 1 000 000]."
        )
        label = QLabel("GL quadrature density (pts/µm)")
        label.setToolTip(self.gl_density_sb.toolTip())
        self.form.addRow(label, self.gl_density_sb)

        # Incident polarisation
        self.incident_pol_cb = QComboBox(self)
        self.incident_pol_cb.addItems(["p", "s"])
        idx = self.incident_pol_cb.findText(self.settings["Incident polarisation"], Qt.MatchFixedString)
        if idx >= 0:
            self.incident_pol_cb.setCurrentIndex(idx)
        self.incident_pol_cb.activated.connect(self.on_incident_pol_cb_activated)
        self.incident_pol_cb.setToolTip("Incident laser polarisation (p = in the plane of incidence, s = perpendicular)")
        label = QLabel("Incident polarisation")
        label.setToolTip(self.incident_pol_cb.toolTip())
        self.form.addRow(label, self.incident_pol_cb)

        # Detected polarisation
        self.detected_pol_cb = QComboBox(self)
        self.detected_pol_cb.addItems(["p", "s", "unpolarised"])
        idx = self.detected_pol_cb.findText(self.settings["Detected polarisation"], Qt.MatchFixedString)
        if idx >= 0:
            self.detected_pol_cb.setCurrentIndex(idx)
        self.detected_pol_cb.activated.connect(self.on_detected_pol_cb_activated)
        self.detected_pol_cb.setToolTip("Detected polarisation (p, s, or 'unpolarised' for no analyser — sums |A_p|² + |A_s|²)")
        label = QLabel("Detected polarisation")
        label.setToolTip(self.detected_pol_cb.toolTip())
        self.form.addRow(label, self.detected_pol_cb)

        # Collection side
        self.collection_side_cb = QComboBox(self)
        self.collection_side_cb.addItems(["superstrate", "substrate"])
        idx = self.collection_side_cb.findText(self.settings["Collection side"], Qt.MatchFixedString)
        if idx >= 0:
            self.collection_side_cb.setCurrentIndex(idx)
        self.collection_side_cb.activated.connect(self.on_collection_side_cb_activated)
        self.collection_side_cb.setToolTip("Which side the detector is on:\n'superstrate' = backscattering (default)\n'substrate' = forward scattering (reversed stack for E_S)")
        label = QLabel("Collection side")
        label.setToolTip(self.collection_side_cb.toolTip())
        self.form.addRow(label, self.collection_side_cb)

        # Collection angle
        self.collection_angle_sb = QDoubleSpinBox(self)
        self.collection_angle_sb.setRange(-1.0, 89.9)
        self.collection_angle_sb.setSingleStep(1.0)
        self.collection_angle_sb.setDecimals(1)
        self.collection_angle_sb.setValue(self.settings["Collection angle"])
        self.collection_angle_sb.valueChanged.connect(self.on_collection_angle_sb_changed)
        self.collection_angle_sb.setToolTip("Collection (detector) angle in degrees.\n−1 = use the angle of incidence (default for backscattering)")
        label = QLabel("Collection angle (°, −1 = auto)")
        label.setToolTip(self.collection_angle_sb.toolTip())
        self.form.addRow(label, self.collection_angle_sb)

        # Coherent layer summation
        self.coherent_cb = QCheckBox(self)
        self.coherent_cb.setChecked(self.settings["Coherent layer summation"])
        self.coherent_cb.toggled.connect(self.on_coherent_cb_toggled)
        self.coherent_cb.setToolTip("When checked, sum amplitudes across Raman-active layers before squaring\n(coherent combination). Default: incoherent (sum intensities).")
        label = QLabel("Coherent layer summation")
        label.setToolTip(self.coherent_cb.toolTip())
        self.form.addRow(label, self.coherent_cb)

        # Approximate E_S = E_L
        self.approximate_cb = QCheckBox(self)
        self.approximate_cb.setChecked(self.settings["Approximate ES"])
        self.approximate_cb.toggled.connect(self.on_approximate_cb_toggled)
        self.approximate_cb.setToolTip("When checked, use E_S ≈ E_L (both at the laser frequency).\nDefault: compute E_S separately at ν_S = ν_L − ν_m per mode.")
        label = QLabel("Approximate E_S = E_L")
        label.setToolTip(self.approximate_cb.toolTip())
        self.form.addRow(label, self.approximate_cb)

        # Phonon boundary correction (Phase 3a)
        self.phonon_bc_cb = QComboBox(self)
        self.phonon_bc_cb.addItems(["none", "NAC", "slab-environment"])
        # Handle legacy 'slab' setting from older saved sessions
        bc_text = self.settings["Phonon boundary correction"]
        if bc_text == "slab":
            bc_text = "slab-environment"
            self.settings["Phonon boundary correction"] = bc_text
        idx = self.phonon_bc_cb.findText(bc_text, Qt.MatchFixedString)
        if idx >= 0:
            self.phonon_bc_cb.setCurrentIndex(idx)
        self.phonon_bc_cb.activated.connect(self.on_phonon_bc_cb_activated)
        self.phonon_bc_cb.setToolTip(
            "Phonon boundary correction:\n"
            "'none' — use bulk TO frequencies\n"
            "'NAC' — non-analytic correction using the layer's optical dielectric (eq 111)\n"
            "'slab-environment' — slab depolarisation with screening from surrounding layers (eqs 112–114)"
        )
        label = QLabel("Phonon boundary correction")
        label.setToolTip(self.phonon_bc_cb.toolTip())
        self.form.addRow(label, self.phonon_bc_cb)

        # Separator: azimuthal sweep
        sweep_label = QLabel("Azimuthal sweep")
        sweep_line = QFrame()
        sweep_line.setFrameShape(QFrame.HLine)
        sweep_hbox = QHBoxLayout()
        sweep_hbox.addWidget(sweep_line)
        sweep_hbox.setAlignment(Qt.AlignVCenter)
        self.form.addRow(sweep_label, sweep_hbox)

        # Number of sweep points
        self.sweep_n_points_sb = QSpinBox(self)
        self.sweep_n_points_sb.setRange(4, 360)
        self.sweep_n_points_sb.setSingleStep(4)
        self.sweep_n_points_sb.setValue(self.settings["Azimuthal sweep points"])
        self.sweep_n_points_sb.valueChanged.connect(self.on_sweep_n_points_sb_changed)
        self.sweep_n_points_sb.setToolTip("Number of azimuthal angles (0°–360°) to compute in the sweep")
        sweep_pts_label = QLabel("Sweep points")
        sweep_pts_label.setToolTip(self.sweep_n_points_sb.toolTip())
        self.form.addRow(sweep_pts_label, self.sweep_n_points_sb)

        # Plot button
        self.azimuthal_sweep_btn = QPushButton("Plot azimuthal sweep", self)
        self.azimuthal_sweep_btn.clicked.connect(self.on_azimuthal_sweep_btn_clicked)
        self.azimuthal_sweep_btn.setToolTip(
            "Calculate Raman intensity over 0°–360° azimuthal angles and open a plot window.\n"
            "Multiple windows may be opened to compare settings."
        )
        self.form.addRow("", self.azimuthal_sweep_btn)


    def partial_incoherence_widget(self):
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

    def smoothing_widget(self):
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

    def open_db_button_clicked(self):
        """Open a new materials' database.

        After the database is opened a refresh is requested.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        logger.debug("Start:: open_db_button_clicked")
        self.open_data_base()
        self.refresh(force=True)
        self.refresh_required = True
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
        logger.debug(f"{self.settings['Legend']} on_partially_incoherent_kernel_sb_changed {value}")
        self.refresh_required = True
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
        logger.debug(f"{self.settings['Legend']} on_partially_incoherent_polynomial_sb_changed {value}")
        self.refresh_required = True
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
        logger.debug(f"{self.settings['Legend']} on_partially_incoherent_samples_sb_changed {value}")
        self.refresh_required = True
        self.settings["Partially incoherent samples"] = value
        self.no_calculations_required = value
        return

    def on_percentage_partial_incoherence_sb_changed(self,value):
        """Handle changes in percentage partial incoherence spinner box value.

        This method updates the 'Percentage partial incoherence' setting based on the new value provided by the relevant
        spinner box. It marks that a refresh is required and logs the change using the configured debugger.

        Parameters
        ----------
        value : int or float
            The new value from the percentage partial incoherence spinner box.

        Returns
        -------
        None

        """        
        logger.debug(f"{self.settings['Legend']} on_percentage_partial_incoherence_sb_changed {value}")
        self.refresh_required = True
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
        logger.debug(f"{self.settings['Legend']} on_global_azimuthal_angl_sb_changed {value}")
        self.refresh_required = True
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
        logger.debug(f"{self.settings['Legend']} on_angle_of_incidence_sb_changed {value}")
        self.refresh_required = True
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
        logger.debug(f"{self.settings['Legend']} Start:: refresh, force = {force}")
        if not self.refresh_required and not force :
            logger.debug(f"{self.settings['Legend']} Finished:: refreshing widget aborted {self.refresh_required} {force}")
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
        if self.settings["Materials database"] != self.DataBase.get_file_name():
            self.DataBase = MaterialsDataBase(self.settings["Materials database"])
            self.settings["Materials database"] = self.DataBase.get_file_name()
            self.database_le.setText(self.settings["Materials database"])
            # Update the possible  material names from the database
            self.material_names = self.set_material_names()
        # Generate the layers from the settings
        self.settings2Layers()
        self.generate_layer_settings()
        # Force recalculation
        self.calculation_required = True
        #
        # Now refresh values that need updating
        #
        # Check the mode 
        index = self.mode_cb.findText(self.settings["Mode"], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        # Work out how many calculations are going to be performed
        self.set_no_calculations_required()
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
        # Crystal Raman widgets (Phase 2f)
        if self.spectroscopy == "Crystal Raman":
            self.laser_wavelength_sb.setValue(self.settings["Laser wavelength nm"])
            self.temperature_sb.setValue(self.settings["Temperature K"])
            self.gl_density_sb.setValue(self.settings.get("GL point density", 20.0))
            idx = self.incident_pol_cb.findText(self.settings["Incident polarisation"], Qt.MatchFixedString)
            if idx >= 0:
                self.incident_pol_cb.setCurrentIndex(idx)
            idx = self.detected_pol_cb.findText(self.settings["Detected polarisation"], Qt.MatchFixedString)
            if idx >= 0:
                self.detected_pol_cb.setCurrentIndex(idx)
            idx = self.collection_side_cb.findText(self.settings["Collection side"], Qt.MatchFixedString)
            if idx >= 0:
                self.collection_side_cb.setCurrentIndex(idx)
            self.collection_angle_sb.setValue(self.settings["Collection angle"])
            self.coherent_cb.setChecked(self.settings["Coherent layer summation"])
            self.approximate_cb.setChecked(self.settings["Approximate ES"])
            idx = self.phonon_bc_cb.findText(self.settings["Phonon boundary correction"], Qt.MatchFixedString)
            if idx >= 0:
                self.phonon_bc_cb.setCurrentIndex(idx)
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        # Redraw the layer information widget
        self.redraw_layer_table()
        for _i in range(20):
            QCoreApplication.processEvents()
        self.refresh_required = False
        logger.debug(f"{self.settings['Legend']} Finished:: refresh, force = {force}")
        return

    def set_no_calculations_required(self):
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
            if layer.get_incoherent_option() == "Incoherent (phase averaging)":
                self.number_of_average_incoherent_layers += 1
        # The Raman calculator runs a single pass over phonon modes regardless of
        # the incoherent settings (those only affect IR), so its per-scenario
        # progress-bar quota is just len(vs_cm1).
        if self.spectroscopy == "Crystal Raman":
            self.no_calculations_required = 1
            return
        # First see how many layers are using phase averaging
        if self.number_of_average_incoherent_layers > 0:
            number_of_samples = self.settings["Number of average incoherence samples"]
            self.no_calculations_required = pow(number_of_samples,self.number_of_average_incoherent_layers)
        else:
            self.no_calculations_required = 1
        # Now see if partial incoherence is being used
        if self.settings["Percentage partial incoherence"] > 0:
            self.no_calculations_required *= self.settings["Partially incoherent samples"]
        return

    def on_mode_cb_activated(self, index):
        """Handle the activation of the mode combo box in the GUI.

        Parameters
        ----------
        index : int
            The index of the selected mode in the combo box. This controls what mode the system will use.  Index 0 corresponds to the 'Transfer matrix' mode, while index 1 corresponds to the 'Scattering matrix' mode.


        Returns
        -------
        None

        Notes
        -----
        This function performs several operations based on the selected mode: - Sets the mode in the settings to either
        'Transfer matrix' or 'Scattering matrix'. - If 'Scattering matrix' mode is selected (index 1), it goes through
        all layers and changes any layer with the incoherent option set to 'Incoherent (intensity)' to be 'Coherent'. -
        Invokes set_no_calculations_required to indicate that new calculations are necessary based on the changed
        settings. - Calls generate_layer_settings to update the layer settings according to the new mode. - Forces a
        refresh of the display to reflect any changes. - Sets a flag indicating that a refresh is required.


        """        
        logger.debug(f"{self.settings['Legend']} Start:: on_mode_cb_activated")
        if index == 0:
            self.settings["Mode"] = "Transfer matrix"
        elif index == 1:
            self.settings["Mode"] = "Scattering matrix"
            # If any layers are using intensity incoherence move them to coherent
            for layer in self.layers:
                incoherentOption = layer.get_incoherent_option()
                if incoherentOption == "Incoherent (intensity)":
                    layer.set_incoherent_option("Coherent")
        self.set_no_calculations_required()
        self.generate_layer_settings()
        self.refresh(force=True)
        self.refresh_required = True
        logger.debug(f"{self.settings['Legend']} Mode changed to {self.settings['Mode']}")
        logger.debug(f"{self.settings['Legend']} Finished:: on_mode_cb_activated")
        return

    # ------------------------------------------------------------------
    # Crystal Raman signal handlers (Phase 2f)
    # ------------------------------------------------------------------

    def on_laser_wavelength_sb_changed(self, value):
        """Handle a change in the laser wavelength spin box."""
        self.settings["Laser wavelength nm"] = value
        self.calculation_required = True
        self.refresh_required = True

    def on_temperature_sb_changed(self, value):
        """Handle a change in the temperature spin box."""
        self.settings["Temperature K"] = value
        self.calculation_required = True
        self.refresh_required = True

    def on_gl_density_sb_changed(self, value):
        """Handle a change in the GL quadrature density spin box."""
        self.settings["GL point density"] = value
        self.calculation_required = True
        self.refresh_required = True

    def on_incident_pol_cb_activated(self, index):
        """Handle a change in the incident polarisation combo box."""
        self.settings["Incident polarisation"] = self.incident_pol_cb.currentText()
        self.calculation_required = True
        self.refresh_required = True

    def on_detected_pol_cb_activated(self, index):
        """Handle a change in the detected polarisation combo box."""
        self.settings["Detected polarisation"] = self.detected_pol_cb.currentText()
        self.calculation_required = True
        self.refresh_required = True

    def on_collection_side_cb_activated(self, index):
        """Handle a change in the collection side combo box."""
        self.settings["Collection side"] = self.collection_side_cb.currentText()
        self.calculation_required = True
        self.refresh_required = True

    def on_collection_angle_sb_changed(self, value):
        """Handle a change in the collection angle spin box."""
        self.settings["Collection angle"] = value
        self.calculation_required = True
        self.refresh_required = True

    def on_coherent_cb_toggled(self, checked):
        """Handle a toggle of the coherent layer summation checkbox."""
        self.settings["Coherent layer summation"] = checked
        self.calculation_required = True
        self.refresh_required = True

    def on_approximate_cb_toggled(self, checked):
        """Handle a toggle of the approximate E_S checkbox."""
        self.settings["Approximate ES"] = checked
        self.calculation_required = True
        self.refresh_required = True

    def on_phonon_bc_cb_activated(self, index):
        """Handle a change in the phonon boundary correction combo box."""
        self.settings["Phonon boundary correction"] = self.phonon_bc_cb.currentText()
        self.calculation_required = True
        self.refresh_required = True

    def on_sweep_n_points_sb_changed(self, value):
        """Handle change in azimuthal sweep point count."""
        self.settings["Azimuthal sweep points"] = value

    def on_azimuthal_sweep_btn_clicked(self):
        """Run the azimuthal sweep and open a new plot window."""
        if not hasattr(self, "vs_cm1") or self.vs_cm1 is None or len(self.vs_cm1) == 0:
            QMessageBox.warning(
                self, "No calculation",
                "Please run a Crystal Raman calculation first to define the frequency axis."
            )
            return
        n_points = self.settings["Azimuthal sweep points"]
        psi_values = np.linspace(0.0, 360.0, n_points, endpoint=False)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            sweep_data = self._run_azimuthal_sweep(self.vs_cm1, psi_values)
        finally:
            QApplication.restoreOverrideCursor()
        if sweep_data is None:
            QMessageBox.warning(
                self, "Sweep failed",
                "Azimuthal sweep could not be completed. Check the Crystal Raman settings."
            )
            return
        from PDielec.GUI.AzimuthalSweepWindow import AzimuthalSweepWindow
        window = AzimuthalSweepWindow(
            sweep_data,
            title=f"Azimuthal sweep — {self.settings['Legend']}",
            parent=None,
        )
        window.show()
        if not hasattr(self, "_sweep_windows"):
            self._sweep_windows = []
        self._sweep_windows.append(window)

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
        layers : list
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
        logger.debug(f"{self.settings['Legend']} Start:: partially_incoherent_calculator")
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
            if layer.get_incoherent_option() == "Incoherent (phase averaging)":
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
                layer.set_phase_shift(beta[index])
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
        layers : list
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
        logger.debug(f"{self.settings['Legend']} Start:: partially_incoherent_calculator")
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
            if layer.is_dielectric():
                crystalLayer = layer
        #
        # Loop over the number of samples requred
        #
        crystalDepth = crystalLayer.get_thickness()
        keepCrystalDepth = crystalDepth
        d = crystalDepth
        t = theta
        p  = phi
        q  = psi
        a  = angleOfIncidence
        fractionalIncoherence = self.settings["Percentage partial incoherence"]/100.0
        for _s in range(self.settings["Partially incoherent samples"]):
            crystalDepth = d + d*( -1 + 2*np.random.rand())*fractionalIncoherence
            crystalLayer.set_thickness(crystalDepth)
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
        crystalLayer.set_thickness(keepCrystalDepth)
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
        layers : list
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
        logger.debug(f"{self.settings['Legend']} Entering the coherent_calculator function")
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
        logger.debug(f"{self.settings['Legend']} About to calculate crystal infrared scenario using pool")
        if self.notebook.pool is None:
            self.notebook.start_pool()
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
        logger.debug(f"{self.settings['Legend']} About to extract results for crystal infrared scenario")
        for v,_r,R,_t,T,eps,errors,largest_exponent in results:
            if self.settings["Mode"] == "Transfer matrix" and errors > 0:
                logger.warning(f"Warning exponential overflow occured at frequency {v} {errors} {largest_exponent}")
            p_reflectance.append(R[0]+R[2])
            s_reflectance.append(R[1]+R[3])
            p_transmittance.append(T[0])
            s_transmittance.append(T[1])
            p_absorbtance.append(1.0 - R[0]-R[2]-T[0])
            s_absorbtance.append(1.0 - R[1]-R[3]-T[1])
            epsilon.append(eps)
        logger.debug(f"{self.settings['Legend']} Finished the coherent_calculator function")
        return ( p_reflectance, s_reflectance, p_transmittance, s_transmittance, p_absorbtance, s_absorbtance, epsilon )

    def calculate(self, vs_cm1):
        """Calculate the crystal spectrum for the range of frequencies in vs_cm1.

        Dispatches to the appropriate calculation method based on the spectroscopy
        type set in the settings tab.

        Parameters
        ----------
        vs_cm1 : array_like
            Array of frequencies for which to calculate the crystal spectrum.

        Returns
        -------
        None

        """
        if self.spectroscopy == "Crystal Infrared":
            self._calculate_infrared(vs_cm1)
        elif self.spectroscopy == "Crystal Raman":
            self._calculate_raman(vs_cm1)
        else:
            logger.error(f"{self.settings['Legend']} calculate: unknown spectroscopy type: {self.spectroscopy}")

    def _compute_slab_modes(self, G_total, scaled_tensors, frequencies_cm1, sigmas_cm1,
                            method="NAC", epsilon_above=None, epsilon_below=None):
        """Compute corrected phonon frequencies and Raman tensors for a crystal layer.

        Two methods are supported:

        **NAC** (non-analytic correction, eq 111):
            The screening tensor uses only the layer's own optical dielectric.
            ``S_bg = n̂ n̂^T / (n̂^T ε_i^b n̂)``
            This is the standard NAC with the wave-vector direction replaced by
            the slab normal.

        **slab-environment** (eqs 112–114):
            The screening accounts for the permittivities of the materials above
            and below the active layer.
            ``ε_e^b = ½ n̂^T (ε_above^b + ε_below^b) n̂``
            ``N_bg = [I + (1/ε_e^b) L (ε_i^b − ε_e^b I)]⁻¹``
            ``S_bg = (1/ε_e^b) N_bg L``

        In both cases the corrected dynamical matrix is (eq 109):
            ``D_layer = D_TO + (1/(ε₀ V)) (Z^mw)^T S_bg Z^mw``

        Parameters
        ----------
        G_total : ndarray, shape (3, 3)
            Combined crystal-to-lab rotation (G_psi @ layer.euler).
        scaled_tensors : list of ndarray, each (3, 3)
            Bulk TO Raman tensors in the crystal frame, scaled by sqrt(V_cell).
        frequencies_cm1 : array_like, shape (n_to_modes,)
            Bulk TO phonon frequencies in cm⁻¹ (including acoustic zeros).
        sigmas_cm1 : array_like, shape (n_to_modes,)
            Bulk TO Lorentzian linewidths in cm⁻¹.
        method : {'NAC', 'slab-environment'}
            Which phonon boundary correction to apply.
        epsilon_above : float or None
            Scalar background permittivity of the layer above (superstrate side).
            Required for 'slab-environment'; ignored for 'NAC'.
        epsilon_below : float or None
            Scalar background permittivity of the layer below (substrate side).
            Required for 'slab-environment'; ignored for 'NAC'.

        Returns
        -------
        slab_freqs : ndarray, shape (3N,)
            Corrected phonon frequencies in cm⁻¹ (all 3N modes, including acoustic).
        slab_tensors : list of ndarray, each (3, 3)
            Corrected Raman tensors in crystal frame (3N total).
        slab_sigmas : ndarray, shape (3N,)
            Linewidths inherited from the dominant TO component.

        """
        nAtoms = self.reader.nions
        n_modes = 3 * nAtoms
        volume_au = self.reader.volume * angs2bohr ** 3
        masses_au = np.array(self.reader.masses) * amu

        # Slab normal in lab frame → crystal frame
        n_hat_crystal = G_total.T @ np.array([0.0, 0.0, 1.0])
        norm = np.linalg.norm(n_hat_crystal)
        if norm > 0.0:
            n_hat_crystal /= norm

        # Depolarisation tensor for a thin slab: L = n̂ ⊗ n̂  (eq 106)
        L = np.outer(n_hat_crystal, n_hat_crystal)

        # Optical dielectric tensor ε_i^b (layer material) in crystal frame
        eps_inf = np.array(self.reader.zerof_optical_dielectric, dtype=float)
        if eps_inf.ndim == 1:
            eps_inf = np.diag(eps_inf)

        I3 = np.eye(3)

        if method == "NAC":
            # NAC approximation (eq 111):
            # S_bg^NAC = L / (n̂^T ε_i^b n̂)
            eps_b_normal = n_hat_crystal @ eps_inf @ n_hat_crystal
            S_bg = L / eps_b_normal
        else:
            # Slab-environment approximation (eqs 112–114):
            # ε_e^b = ½ n̂^T (ε_above^b + ε_below^b) n̂  (eq 112)
            eps_e_above = epsilon_above if epsilon_above is not None else 1.0
            eps_e_below = epsilon_below if epsilon_below is not None else 1.0
            epsilon_e = 0.5 * (eps_e_above + eps_e_below)
            # N_bg^slab = [I + (1/ε_e^b) L (ε_i^b − ε_e^b I)]⁻¹  (eq 113)
            N_bg = np.linalg.inv(I3 + (1.0 / epsilon_e) * L @ (eps_inf - epsilon_e * I3))
            # S_bg^slab = (1/ε_e^b) N_bg L  (eq 114)
            S_bg = (1.0 / epsilon_e) * N_bg @ L

        # Mass-weighted Born charges Z'[α, κβ] = Z*[κ, α, β] / √M_κ
        born_charges = np.array(self.reader.born_charges)
        Z_mat = np.zeros((3, n_modes))
        for kappa in range(nAtoms):
            inv_sqrtM = 1.0 / math.sqrt(masses_au[kappa])
            for beta in range(3):
                Z_mat[:, kappa * 3 + beta] = born_charges[kappa, :, beta] * inv_sqrtM

        # Dynamical matrix correction (eq 109):
        # D_layer = D_TO + (1/(ε₀ V)) (Z^mw)^T S_bg Z^mw
        # In Gaussian units: ε₀ → 1/(4π), so factor is 4π/V
        delta_D = (4.0 * np.pi / volume_au) * (Z_mat.T @ S_bg @ Z_mat)
        D_TO = np.array(self.reader.hessian, dtype=float)

        # Diagonalize the corrected dynamical matrix
        eig_val, eig_vec = np.linalg.eigh(D_TO + delta_D)

        # Frequencies in cm⁻¹ (preserve sign for unstable modes)
        slab_freqs = np.array([
            (math.sqrt(abs(ev)) / wavenumber) * (1.0 if ev >= 0.0 else -1.0)
            for ev in eig_val
        ])

        # Build U_TO: rows = TO eigenvectors flattened, shape (n_to_modes, 3N)
        n_to_modes = len(self.reader.mass_weighted_normal_modes)
        U_TO = np.zeros((n_to_modes, n_modes))
        for imode, mode in enumerate(self.reader.mass_weighted_normal_modes):
            col = 0
            for atom in mode:
                U_TO[imode, col:col + 3] = atom
                col += 3

        # Overlap matrix C[n, m] = <u_n^TO | u_m^slab>
        C = U_TO @ eig_vec  # (n_to_modes, 3N)

        # Transform Raman tensors and inherit linewidths from dominant TO component
        sigmas = np.asarray(sigmas_cm1, dtype=float)
        slab_tensors = []
        slab_sigmas = np.zeros(n_modes)
        for p_idx in range(n_modes):
            R_p = np.zeros((3, 3), dtype=float)
            for n_to in range(n_to_modes):
                R_p += C[n_to, p_idx] * np.asarray(scaled_tensors[n_to], dtype=float)
            slab_tensors.append(R_p)
            dominant_to = int(np.argmax(np.abs(C[:, p_idx])))
            slab_sigmas[p_idx] = sigmas[dominant_to] if dominant_to < len(sigmas) else 5.0

        return slab_freqs, slab_tensors, slab_sigmas

    def _build_raman_calculator(self, psi_rad):
        """Build and return a configured LayeredRamanCalculator for the given azimuthal angle.

        Constructs the GTMcore multilayer system, applies global azimuthal rotation,
        builds the list of RamanLayer descriptors (with optional phonon boundary
        correction: NAC or slab-environment), and returns a fully configured
        LayeredRamanCalculator ready to call calculate_mode_intensities().

        Parameters
        ----------
        psi_rad : float
            Global azimuthal rotation angle in radians.

        Returns
        -------
        LayeredRamanCalculator or None
            Configured calculator, or None if prerequisites are missing or invalid.

        """
        raman_tensors = self.reader.get_raman_tensors()
        if raman_tensors is None:
            logger.warning(f"{self.settings['Legend']} _build_raman_calculator: no Raman tensors")
            return None

        theta = 0.0
        phi = 0.0
        angle_of_incidence = np.radians(self.settings["Angle of incidence"])

        # Compute HKL rotation matrices for all tensor layers
        for layer in self.layers:
            if layer.is_tensor():
                hkl = layer.get_hkl()
                if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
                    logger.warning(f"{self.settings['Legend']} _build_raman_calculator: hkl=[0,0,0]")
                    return None
                layer.calculate_euler_matrix()

        # Build GTM multilayer system
        mode = self.settings["Mode"]
        exponent_threshold = self.exponent_threshold
        superstrate = GTM.SemiInfiniteLayer(self.layers[0], exponent_threshold=exponent_threshold)
        substrate   = GTM.SemiInfiniteLayer(self.layers[-1], exponent_threshold=exponent_threshold)
        selected_layers = self.layers[1:-1]
        gtm_layers = []
        for layer in selected_layers:
            incoherent_option = layer.get_incoherent_option()
            gtm_layers.append(gtm_methods[incoherent_option](layer, exponent_threshold=exponent_threshold))
        if mode == "Scattering matrix":
            system = GTM.ScatteringMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtm_layers)
        else:
            system = GTM.TransferMatrixSystem(substrate=substrate, superstrate=superstrate, layers=gtm_layers)

        # Apply global azimuthal rotation to all GTM layers
        system.superstrate.set_euler(theta, phi, psi_rad)
        system.substrate.set_euler(theta, phi, psi_rad)
        for gtm_layer in system.layers:
            gtm_layer.set_euler(theta, phi, psi_rad)

        # Global azimuthal rotation matrix (rotation around z by psi_rad)
        G_psi = np.array([
            [ np.cos(psi_rad), -np.sin(psi_rad), 0.0],
            [ np.sin(psi_rad),  np.cos(psi_rad), 0.0],
            [ 0.0,              0.0,              1.0],
        ])

        # Scale stored Raman tensors by sqrt(volume) to get actual R = sqrt(V) * T_stored
        volume_ang3 = self.reader.volume
        scale = np.sqrt(volume_ang3)
        scaled_tensors = [scale * np.asarray(R, dtype=float) for R in raman_tensors]

        # Phonon frequencies and linewidths from SettingsTab
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        sigmas_cm1      = self.notebook.settingsTab.sigmas_cm1

        # Phonon boundary correction: 'none', 'NAC', or 'slab-environment'
        phonon_bc = self.settings.get("Phonon boundary correction", "none")
        # Handle legacy 'slab' setting
        if phonon_bc == "slab":
            phonon_bc = "slab-environment"
        has_hessian = hasattr(self.reader, "hessian") and self.reader.hessian is not None
        has_born    = len(self.reader.born_charges) > 0
        has_modes   = np.any(self.reader.mass_weighted_normal_modes)
        has_optical = (hasattr(self.reader, "zerof_optical_dielectric")
                       and self.reader.zerof_optical_dielectric is not None)
        can_correct = (phonon_bc in ("NAC", "slab-environment")
                       and has_hessian and has_born and has_modes and has_optical)
        if phonon_bc in ("NAC", "slab-environment") and not can_correct:
            logger.warning(
                f"{self.settings['Legend']} _build_raman_calculator: '{phonon_bc}' requested but "
                "Born charges / hessian / optical dielectric not available — "
                "falling back to bulk TO frequencies"
            )

        # For slab-environment, determine the normal-component permittivity of
        # the layers above and below each dielectric layer.  Use the full layer
        # list (including superstrate [0] and substrate [-1]).
        all_layers = self.layers  # [superstrate, layer1, ..., layerN, substrate]

        def _get_layer_scalar_permittivity(layer):
            """Return a scalar optical permittivity for a layer.

            Uses the material's ``get_optical_permittivity()`` method.
            For scalar materials this returns the scalar value directly.
            For tensor materials the isotropic average (trace/3) is returned.
            """
            mat = layer.get_material()
            try:
                eps = mat.get_optical_permittivity()
            except Exception:
                return 1.0
            if np.isscalar(eps) or (hasattr(eps, "ndim") and eps.ndim == 0):
                return float(np.real(eps))
            # Tensor: return trace/3 as isotropic average
            return float(np.real(np.trace(eps))) / 3.0

        corrected_freqs_shared = None
        corrected_sigmas_shared = None
        if can_correct:
            first_G = None
            first_env = (None, None)
            for i, scl in enumerate(selected_layers):
                if scl.is_dielectric():
                    first_G = G_psi @ scl.euler
                    if phonon_bc == "slab-environment":
                        # Layer above: previous in all_layers (index i in selected = index i+1 in all_layers)
                        above_layer = all_layers[i] if i >= 0 else all_layers[0]
                        # Layer below: next in all_layers
                        below_layer = all_layers[i + 2] if (i + 2) < len(all_layers) else all_layers[-1]
                        first_env = (_get_layer_scalar_permittivity(above_layer),
                                     _get_layer_scalar_permittivity(below_layer))
                    break
            if first_G is not None:
                corrected_freqs_shared, _, corrected_sigmas_shared = self._compute_slab_modes(
                    first_G, scaled_tensors, frequencies_cm1, sigmas_cm1,
                    method=phonon_bc,
                    epsilon_above=first_env[0], epsilon_below=first_env[1],
                )

        # Build RamanLayer descriptors for each Raman-active (dielectric) layer
        raman_layer_list = []
        for sys_idx, scl in enumerate(selected_layers):
            if not scl.is_dielectric():
                continue
            G_total = G_psi @ scl.euler
            if can_correct and corrected_freqs_shared is not None:
                # Determine environment permittivities for this specific layer
                eps_above, eps_below = None, None
                if phonon_bc == "slab-environment":
                    above_layer = all_layers[sys_idx] if sys_idx >= 0 else all_layers[0]
                    below_layer = all_layers[sys_idx + 2] if (sys_idx + 2) < len(all_layers) else all_layers[-1]
                    eps_above = _get_layer_scalar_permittivity(above_layer)
                    eps_below = _get_layer_scalar_permittivity(below_layer)
                _, layer_corrected_tensors, _ = self._compute_slab_modes(
                    G_total, scaled_tensors, frequencies_cm1, sigmas_cm1,
                    method=phonon_bc,
                    epsilon_above=eps_above, epsilon_below=eps_below,
                )
                rl = RamanLayer(
                    layer_index=sys_idx,
                    phonon_frequencies_cm1=corrected_freqs_shared,
                    raman_tensors=layer_corrected_tensors,
                    rotation_matrix=G_total,
                )
            else:
                rl = RamanLayer(
                    layer_index=sys_idx,
                    phonon_frequencies_cm1=frequencies_cm1,
                    raman_tensors=scaled_tensors,
                    rotation_matrix=G_total,
                )
            raman_layer_list.append(rl)

        if not raman_layer_list:
            logger.warning(f"{self.settings['Legend']} _build_raman_calculator: no Raman-active dielectric layers in stack")
            return None

        # Linewidths: use corrected sigmas when phonon correction is active
        linewidths = corrected_sigmas_shared if (can_correct and corrected_sigmas_shared is not None) else sigmas_cm1

        laser_wavelength_nm = self.settings.get("Laser wavelength nm", 532.0)
        laser_freq_cm1    = 1.0e7 / laser_wavelength_nm
        incident_pol      = self.settings.get("Incident polarisation", "p")
        detected_pol      = self.settings.get("Detected polarisation", "unpolarised")
        temperature_K     = self.settings.get("Temperature K", 298.0)
        density           = self.settings.get("GL point density", 20.0)
        # Compute n_gauss from density (pts/µm) and the thickest Raman-active layer
        _min_points = 5
        _max_points = 1_000_000
        max_thick_um = max(
            (system.layers[rl.layer_index].thick * 1.0e6 for rl in raman_layer_list),
            default=1.0,
        )
        n_gauss_raw = int(round(density * max_thick_um))
        if n_gauss_raw < _min_points:
            print(
                f"Crystal Raman: GL density {density} pts/µm × {max_thick_um:.4g} µm = "
                f"{n_gauss_raw} points; using minimum of {_min_points}."
            )
            n_gauss = _min_points
        elif n_gauss_raw > _max_points:
            print(
                f"Crystal Raman: GL density {density} pts/µm × {max_thick_um:.4g} µm = "
                f"{n_gauss_raw} points; using maximum of {_max_points:,}."
            )
            n_gauss = _max_points
        else:
            n_gauss = n_gauss_raw
        collection_side   = self.settings.get("Collection side", "superstrate")
        collection_angle  = self.settings.get("Collection angle", -1.0)
        coherent_layers   = self.settings.get("Coherent layer summation", False)
        approximate_es    = self.settings.get("Approximate ES", False)

        collection_angle_rad = angle_of_incidence if collection_angle < 0.0 else np.radians(collection_angle)

        return LayeredRamanCalculator(
            system=system,
            raman_layers=raman_layer_list,
            laser_frequency_cm1=laser_freq_cm1,
            incident_angle_rad=angle_of_incidence,
            incident_pol=incident_pol,
            detected_pol=detected_pol,
            temperature_K=temperature_K,
            linewidths_cm1=linewidths,
            n_gauss=n_gauss,
            collection_side=collection_side,
            collection_angle_rad=collection_angle_rad,
            coherent_layers=coherent_layers,
            approximate_es=approximate_es,
        )

    def _calculate_raman(self, vs_cm1):
        """Calculate the layered crystal Raman spectrum via GTM field integration.

        Builds a GTMcore multilayer system from the current layer stack, then
        uses LayeredRamanCalculator to evaluate the Raman source-overlap integral
        with Gauss-Legendre quadrature inside each Raman-active (dielectric) layer.

        Results are stored in ``self.raman_spectrum``, ``self.raman_mode_frequencies``,
        ``self.raman_mode_intensities``, and ``self.raman_mode_sigmas``.

        Parameters
        ----------
        vs_cm1 : array_like
            Raman-shift axis in cm⁻¹ at which to evaluate the spectrum.

        Returns
        -------
        None

        """
        logger.debug(f"{self.settings['Legend']} Start:: _calculate_raman")
        if not self.calculation_required:
            logger.debug(f"{self.settings['Legend']} _calculate_raman aborted: calculation_required=False")
            return
        QCoreApplication.processEvents()
        self.vs_cm1 = vs_cm1

        settings = self.notebook.mainTab.settings
        program = settings["Program"]
        filename = self.notebook.mainTab.get_full_file_name()
        if self.reader is None or program == "" or filename == "":
            logger.debug(f"{self.settings['Legend']} _calculate_raman aborted: reader/program/filename not set")
            return

        raman_tensors = self.reader.get_raman_tensors()
        if raman_tensors is None:
            logger.warning(f"{self.settings['Legend']} _calculate_raman: reader has no Raman tensors")
            self.calculation_required = False
            return

        for layer in self.layers:
            if layer.is_tensor():
                hkl = layer.get_hkl()
                if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
                    QMessageBox.about(self, "", f"Unable to calculate surface for scenario {self.settings['Legend']}, hkl=[0,0,0]")
                    return

        psi = np.radians(self.settings["Global azimuthal angle"])
        calculator = self._build_raman_calculator(psi)
        if calculator is None:
            self.calculation_required = False
            return

        ref_layer = calculator.raman_layers[0]
        n_modes = len(ref_layer.phonon_frequencies_cm1)
        n_freqs = len(vs_cm1)
        _count = [0]
        _updated = [0]

        def _progress_callback():
            _count[0] += 1
            target = (n_freqs * _count[0]) // n_modes if n_modes > 0 else n_freqs
            if target > _updated[0]:
                self.notebook.progressbars_update(increment=target - _updated[0])
                _updated[0] = target

        try:
            active_freqs, active_ints, active_sigmas = calculator.calculate_mode_intensities(
                progress_callback=_progress_callback
            )
        finally:
            remaining = n_freqs - _updated[0]
            if remaining > 0:
                self.notebook.progressbars_update(increment=remaining)

        if len(active_freqs) > 0:
            spectrum = lorentzian_broaden(active_freqs, active_ints, active_sigmas, np.asarray(vs_cm1))
        else:
            spectrum = np.zeros(len(vs_cm1))

        self.raman_mode_frequencies = active_freqs
        self.raman_mode_intensities = active_ints
        self.raman_mode_sigmas = active_sigmas
        self.raman_spectrum = spectrum
        self.calculation_required = False
        logger.debug(f"{self.settings['Legend']} Finished:: _calculate_raman")

    def _run_azimuthal_sweep(self, vs_cm1, psi_values):
        """Compute Raman mode intensities and spectra over a range of azimuthal angles.

        Parameters
        ----------
        vs_cm1 : ndarray
            Frequency axis in cm⁻¹.
        psi_values : ndarray
            Azimuthal angles in degrees.

        Returns
        -------
        dict or None
            Keys: 'psi_array', 'mode_frequencies', 'mode_intensities_2d',
            'vs_cm1', 'spectra_2d'. Returns None on failure.

        """
        if self.reader is None:
            logger.warning(f"{self.settings['Legend']} _run_azimuthal_sweep: no reader")
            return None
        if self.reader.get_raman_tensors() is None:
            logger.warning(f"{self.settings['Legend']} _run_azimuthal_sweep: no Raman tensors")
            return None
        for layer in self.layers:
            if layer.is_tensor():
                hkl = layer.get_hkl()
                if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
                    logger.warning(f"{self.settings['Legend']} _run_azimuthal_sweep: hkl=[0,0,0]")
                    return None

        vs = np.asarray(vs_cm1)
        all_intensities = []
        all_spectra = []
        mode_freqs = None

        for psi_deg in psi_values:
            calculator = self._build_raman_calculator(np.radians(psi_deg))
            if calculator is None:
                logger.warning(f"{self.settings['Legend']} _run_azimuthal_sweep: calculator failed at psi={psi_deg:.1f}")
                return None
            freqs, ints, sigmas = calculator.calculate_mode_intensities()
            if mode_freqs is None:
                mode_freqs = freqs
            spectrum = lorentzian_broaden(freqs, ints, sigmas, vs) if len(freqs) > 0 else np.zeros(len(vs))
            n_m = len(mode_freqs) if mode_freqs is not None else 0
            all_intensities.append(ints if len(ints) == n_m else np.zeros(n_m))
            all_spectra.append(spectrum)
            QCoreApplication.processEvents()

        if mode_freqs is None or len(mode_freqs) == 0:
            logger.warning(f"{self.settings['Legend']} _run_azimuthal_sweep: no active modes found")
            return None

        return {
            "psi_array": psi_values,
            "mode_frequencies": mode_freqs,
            "mode_intensities_2d": np.array(all_intensities),
            "vs_cm1": vs,
            "spectra_2d": np.array(all_spectra),
        }

    def _calculate_infrared(self,vs_cm1):
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
        This method updates the object's state by calculating and setting various optical properties (reflectance,
        transmittance, absorbtance) based on the input wavenumbers, layer configurations, and other settings such as the
        angle of incidence. It requires that settings, frequencies, and other parameters be previously defined and
        valid. The method also makes use of different calculator objects depending on the coherence conditions and
        partial incoherence percentage. It handles conditions such as missing program, file reader, or filename settings
        by aborting the calculations and logging appropriate messages through a `debugger`.

        """        
        logger.debug(f"{self.settings['Legend']} Start:: calculate - number of frequencies {len(vs_cm1)}")
        if not self.calculation_required:
            logger.debug(f"{self.settings['Legend']} Finished:: calculate aborted because calculation_required false")
            return
        QCoreApplication.processEvents()
        self.vs_cm1 = vs_cm1
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings["Program"]
        filename = self.notebook.mainTab.get_full_file_name()
        if self.reader is None:
            logger.debug(f"{self.settings['Legend']} Finished:: Calculate aborting - no reader")
            return
        if program == "":
            logger.debug(f"{self.settings['Legend']} Finished:: Calculate aborting - no program")
            return
        if filename == "":
            logger.debug(f"{self.settings['Legend']} Finished:: Calculate aborting - no file")
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
            if layer.is_tensor():
                hkl = layer.get_hkl()
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
        logger.debug(f"{self.settings['Legend']} Finished:: calculate - number of frequencies {len(vs_cm1)}")
        return

    def get_result(self, vs_cm1, plot_type):
        """Return a particular result.

        Parameters
        ----------
        vs_cm1 : list
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
        logger.debug(f"{self.settings['Legend']} Start:: get_result {len(vs_cm1)} {plot_type}")
        self.get_results(vs_cm1)
        logger.debug(f"{self.settings['Legend']} Finished:: get_result {len(vs_cm1)} {plot_type}")
        return{ "Crystal Reflectance (P polarisation)"  : self.p_reflectance,
                "Crystal Reflectance (S polarisation)"  : self.s_reflectance,
                "Crystal Transmittance (P polarisation)": self.p_transmittance,
                "Crystal Transmittance (S polarisation)": self.s_transmittance,
                "Crystal Absorbtance (P polarisation)"  : self.p_absorbtance,
                "Crystal Absorbtance (S polarisation)"  : self.s_absorbtance,
                "Crystal Raman"                         : self.raman_spectrum,
        }.get(plot_type)

    def get_results(self, vs_cm1):
        """Return the results of the crystal infrared calculation.

        If a refresh has been requested it is performed before the calculate method is called.

        Parameters
        ----------
        vs_cm1 : list
            The frequencies in cm-1

        Returns
        -------
        None

        """
        logger.debug(f"{self.settings['Legend']} Start:: get_results {len(vs_cm1)} {self.refresh_required}")
        if len(vs_cm1) > 0 and (self.refresh_required or len(self.vs_cm1) != len(vs_cm1) or self.vs_cm1[0] != vs_cm1[0] or self.vs_cm1[1] != vs_cm1[1]) :
            logger.debug(f"{self.settings['Legend']} get_results recalculating")
            self.refresh()
            self.calculate(vs_cm1)
        else:
            logger.debug(f"{self.settings['Legend']} get_results no need for recalculation")
            #self.notebook.progressbars_update(increment=len(vs_cm1))
        logger.debug(f"{self.settings['Legend']} Finished:: get_results {len(vs_cm1)} {self.refresh_required}")
        return

