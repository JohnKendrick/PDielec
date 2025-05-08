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
"""ViewerTab module."""
import copy
import os
from collections import deque

import numpy as np
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from PDielec import Calculator
from PDielec.Constants import elemental_colours, jmol_elemental_colours, vesta_elemental_colours
from PDielec.GUI.OpenGLWidget import OpenGLWidget
from PDielec.PrimitiveCell import PrimitiveCell

# Need the SuperCell & PrimitiveCell classes
from PDielec.SuperCell import SuperCell

# Import plotting requirements
from PDielec.Utilities import Debug


class ViewerTab(QWidget):
    """A class representing a viewer tab in a graphical user interface, specifically for visualizing molecular structures and vibrational modes.

    This class integrates several widgets to enable the selection and visualization of various 
    parameters relevant to molecular structures and their vibrational modes. It allows for 
    interactive change of settings like atom scaling, super cell configuration, and visualization 
    options such as showing arrows for vibrational modes or animation. It also includes functionalities
    to save the visualization in different formats.

    Parameters
    ----------
    parent : QWidget
        The parent widget to which this widget belongs.
    debug : bool, optional
        Flag to enable debugging messages, by default False.

    Attributes
    ----------
    debug : bool
        If debug mode is enabled.
    refreshRequired : bool
        Indicates if a refresh of the visualization is needed due to settings changes.
    UVW : deque
        A deque containing the vibrational mode vectors.
    settings : dict
        Stores various visualization settings such as selected mode, atom scaling, cell color, etc.
    light_switches : list of bool
        A list indicating the state (on or off) of up to 8 light switches for visualization.
    plot_types : list of str
        Available types of plots for visualization.
    plot_type_index : int
        The index of the currently selected plot type.
    number_of_molecules : int
        The number of molecules in the visualization.
    standard_cell, unit_cell, super_cell, cell_edges, cell_corners : Various
        Attributes related to the molecular/crystalline structure and its visualization.
    element_colours : dict
        A mapping of elements to their colors used in the visualization.
    element_names : list of str
        Names of the elements present in the visualization.
    species : list of str
        The species of the elements/molecules considered in the visualization.
    image_filename : str
        The filename where the visualization can be saved.
    notebook : QWidget
        The parent notebook widget containing this ViewerTab.
    reader : Reader
        The reader object used for reading structural information.
    frequencies_cm1 : list
        List of frequencies of the vibrational modes in cm^-1.

    Methods
    -------
    There are multiple methods for event handling (e.g., on_filename_le_return, on_super_cell_changed_a),
    methods for recalculating visualization parameters (e.g., calculate, calculatePhasePositions),
    and methods for changing visualization settings interactively and programatically (e.g., setColour, refresh).
    Each method typically adjusts the ViewerTab's state or visualization output accordingly based on user interaction or programmatic changes.

    calculate
        Obtains frequencies, normal modes, unit cell of whole molecules and supercell
    calculatePhasePositions
        Calculates the positions of atoms during a vibration
    on_arrow_radius_changed
        Handles the change in arrow radius widget
    on_atom_scaling_changed
        Handles the change in atom scaling widget
    on_bond_radius_changed
        Handles the change in bond radius widget
    on_cell_radius_changed
        Handle the change in unit cell radius widget
    on_coloured_button_clicked
        Handles a click on the coloured button
    on_coloured_element_clicked
        Changes the colour of an element
    on_edit_button_clicked
        Edit the transform matrix in the GUI
    on_filename_button_clicked
        Handles a click on the file name button
    on_filename_le_changed
        Handles a change to the file name
    on_filename_le_return
        Handles a return press in the file name widget
    on_guess_button_clicked
        Guess the transform to a primitive cell
    on_light_switches_cb_activated
        Activates the light switch combobox
    on_maximum_displacement_changed
        Handles a change in the maximum displace allowed in following molecular vibrations
    on_plottype_cb_changed
        Handles a change in the plot type combobox
    on_selected_mode_changed
        Handles a change in the vibrational mode which is being displayed
    on_super_cell_changed_a
        Handles the change in the a supercell parameter
    on_super_cell_changed_b
        Handles the change in the b supercell parameter
    on_super_cell_changed_c
        Handles the change in the c supercell parameter
    on_reset_button_clicked
        Changes the transform back to a unit matrix
    plot
        Plot the 3D view of the molecule
    plot_animation
        Plot an animated 3D view of the vibration being studied
    plot_arrows
        Plot an 3D view of the vibration being studied using arrows
    plot_none
        Plots the molecule only, no vibrations
    refresh
        Refresh the GUI with the latest information
    requestRefresh
        Something has changed in the GUI that will need a refresh
    save_as_cif
        Write out a cif file
    save_mode_as_cif
        Write out a cif file with the mode displacememnt information
    setColour
        Set the colour of the specified element

    """

    def __init__(self, parent, debug=False ):
        """Initialise the viewer tab with default settings, UI components, and signal connections.

        Parameters
        ----------
        parent : QWidget
            The parent widget to which this widget belongs.
        debug : bool, optional
            Flag to enable debugging messages, by default False.

        Notes
        -----
        This method sets up various UI components such as spin boxes for selecting visualization mode,
        super cell size, atom scaling, bond radius, cell radius, arrow radius, maximum displacement,
        and options for plotting type and light switches. It initializes default settings for visualization,
        sets up connections for UI component signals to their respective slot methods, and arranges all UI components
        using layout managers. Also, it sets up OpenGL widget for rendering.

        """        
        super(QWidget, self).__init__(parent)
        self.debugger = Debug(debug,"ViewerTab")
        self.debugger.print("Start:: initialisation")
        self.debug = debug
        self.refreshRequired = True
        self.setWindowTitle("Viewer")
        self.settings = {}
        self.UVW = deque()
        self.settings["Selected mode"] = 4
        self.settings["Atom scaling"] = 0.5
        self.settings["Maximum displacement"]  = 1.0
        self.settings["Bond colour"]           = [  80,  80,  80, 255 ]
        self.settings["Bond radius"]           = 0.1
        self.settings["Cell colour"]           = [ 255,   0,   0, 255 ]
        self.settings["Cell radius"]           = 0.1
        self.settings["Text colour"]           = [ 255, 255,   0, 255 ]
        self.settings["Text size"]             = "24"
        self.settings["Background colour"]     = [ 120, 120, 120, 255 ]
        self.settings["Arrow colour"]          = [   0, 255,   0, 255 ]
        self.settings["Arrow radius"]          = 0.07
        self.settings["Number of phase steps"] = 41
        self.settings["Super Cell"] =  [ 1, 1, 1 ]
        self.settings["Transform"] =  [ [ "1","0","0" ], [ "0","1","0" ], ["0","0","1"] ]
        self.settings["hkl"] =  (0, 0, 1)
        self.settings["uvw"] =  (1, 0, 0)
        self.settings["Element colours"] = None
        self.settings["Element palette"] = "Jmol"
        self.light_switches = [False]*8
        self.light_switches[0] = True
        self.light_switches[1] = True
        self.plot_types = ["Animation","Arrows","No arrows or animation"]
        self.plot_type_index = 1
        self.number_of_molecules = 0
        self.standard_cell = None
        self.super_cell = None
        self.cell_edges = None
        self.cell_corners = None
        # element_colours is a dictionary
        self.element_colours = elemental_colours
        self.element_names = []
        self.species       = []
        self.element_coloured_buttons = []
        self.bond_cell_background_arrow_buttons = []
        self.bond_cell_background_arrow_names = ["Background","Cell","Bonds","Arrows"]
        self.image_filename = ""
        self.transform_tab_entry = None
        self.transformed_cell = None
        self.hkl_tab_entry = None
        self.toggle_names =  ["Show cell labels",
                              "Show cell",
                              "Show orientation",
                              "Show bonds",
                              "Show atoms",
                             ]
        self.settings["Toggle states"] = [ True,
                               True,
                               True,
                               True,
                               True,
                             ]
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        self.frequencies_cm1  = self.notebook.settingsTab.frequencies_cm1
        # Create last tab - ViewerTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The selected mode
        # Mode numbering starts at 1
        #
        self.selected_mode_sb = QSpinBox(self)
        if self.frequencies_cm1 is not None and len(self.frequencies_cm1) > 1:
            self.selected_mode_sb.setRange(1,len(self.frequencies_cm1))
        else:
            self.selected_mode_sb.setRange(1,2000)
        self.selected_mode_sb.setValue(self.settings["Selected mode"])
        self.selected_mode_sb.setToolTip("Set the mode to be visualised")
        self.selected_mode_sb.valueChanged.connect(self.on_selected_mode_changed)
        label = QLabel("Select phonon mode", self)
        label.setToolTip("Set the mode to be visualised")
        form.addRow(label, self.selected_mode_sb)
        # Add frequency of mode
        self.frequency_le = QLineEdit(self)
        self.frequency_le.setEnabled(False)
        self.frequency_le.setText(f"{0.0}")
        label = QLabel("Frequency (cm-1)", self)
        form.addRow(label, self.frequency_le)
        #
        # The transform window popup
        #
        self.transform_tab_entry = self.createTransformTabEntry()
        #
        # The hkl tab entry
        #
        self.hkl_tab_entry = self.createHKLTabEntry()
        #
        # The super-cell widget
        #
        self.super_cell_widget = QWidget(self)
        self.super_cell_widget.setToolTip("Generate a super-cell\nThe super-cell is a multiple of the DFT cell, not the transformed cell")
        self.super_cell_hbox = QHBoxLayout()
        self.super_cell_spinBoxes = []
        super_cell_changed = [self.on_super_cell_changed_a, self.on_super_cell_changed_b, self.on_super_cell_changed_c]
        super_cell_tooltip = ["Size in a", "Size in b", "Size in c"]
        for change_function,tip in zip(super_cell_changed,super_cell_tooltip):
            spinBox = QSpinBox(self)
            spinBox.setToolTip(tip)
            spinBox.setRange(1,20)
            spinBox.setSingleStep(1)
            spinBox.valueChanged.connect(change_function)
            self.super_cell_spinBoxes.append(spinBox)
            self.super_cell_hbox.addWidget(spinBox)
        self.super_cell_widget.setLayout(self.super_cell_hbox)
        #
        # The atom scaling
        #
        self.atom_scaling_sb = QDoubleSpinBox(self)
        self.atom_scaling_sb.setRange(0.01,5.0)
        self.atom_scaling_sb.setSingleStep(0.1)
        self.atom_scaling_sb.setDecimals(2)
        self.atom_scaling_sb.setValue(self.settings["Atom scaling"])
        self.atom_scaling_sb.setToolTip("Scale the covalent radii to determine the size of the atomic spheres")
        self.atom_scaling_sb.valueChanged.connect(self.on_atom_scaling_changed)
        #
        # The bond width
        #
        self.bond_radius_sb = QDoubleSpinBox(self)
        self.bond_radius_sb.setRange(0.01,1.0)
        self.bond_radius_sb.setSingleStep(0.1)
        self.bond_radius_sb.setDecimals(2)
        self.bond_radius_sb.setValue(self.settings["Bond radius"])
        self.bond_radius_sb.setToolTip("Determines the size of the bonds drawn")
        self.bond_radius_sb.valueChanged.connect(self.on_bond_radius_changed)
        #
        # The cell width
        #
        self.cell_radius_sb = QDoubleSpinBox(self)
        self.cell_radius_sb.setRange(0.01,1.0)
        self.cell_radius_sb.setSingleStep(0.1)
        self.cell_radius_sb.setDecimals(2)
        self.cell_radius_sb.setValue(self.settings["Cell radius"])
        self.cell_radius_sb.setToolTip("Determines the radius of the cell outline")
        self.cell_radius_sb.valueChanged.connect(self.on_cell_radius_changed)
        #
        # The arrow radius
        #
        self.arrow_radius_sb = QDoubleSpinBox(self)
        self.arrow_radius_sb.setRange(0.01,1.0)
        self.arrow_radius_sb.setSingleStep(0.1)
        self.arrow_radius_sb.setDecimals(2)
        self.arrow_radius_sb.setValue(self.settings["Arrow radius"])
        self.arrow_radius_sb.setToolTip("Determines the radius of the arrow outline")
        self.arrow_radius_sb.valueChanged.connect(self.on_arrow_radius_changed)
        #
        # The maximum displacement
        #
        self.maximum_displacement_sb = QDoubleSpinBox(self)
        self.maximum_displacement_sb.setRange(0.1,10.0)
        self.maximum_displacement_sb.setSingleStep(0.1)
        self.maximum_displacement_sb.setDecimals(1)
        self.maximum_displacement_sb.setValue(self.settings["Maximum displacement"])
        self.maximum_displacement_sb.setToolTip("Set the size of the maximum displacement")
        self.maximum_displacement_sb.valueChanged.connect(self.on_maximum_displacement_changed)
        #
        # Add a combo box to select which type of plot
        #
        self.plottype_cb = QComboBox(self)
        self.plottype_cb.setToolTip("The plot can either be an animation or the modes can be shown by arrow")
        for choice in self.plot_types:
            self.plottype_cb.addItem(choice)
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        self.plottype_cb.currentIndexChanged.connect(self.on_plottype_cb_changed)
        label = QLabel("Choose the plot type", self)
        label.setToolTip("The plot can either be an animation or the modes can be shown by arrow")
        form.addRow(label, self.plottype_cb)
        #
        # Toggle switches
        #
        self.toggles_cb = QComboBox(self)
        self.toggles_cb.setToolTip("The toggles switches")
        for toggle,state in zip(self.toggle_names,self.settings["Toggle states"]):
            string = f"{toggle} is on" if state else f"{toggle} is off"
            self.toggles_cb.addItem(string)
        self.toggles_cb.activated.connect(self.on_toggles_cb_activated)
        #
        # Light switches
        #
        self.light_switches_cb = QComboBox(self)
        self.light_switches_cb.setToolTip("Toogle the light switches on or off")
        for index,light in enumerate(self.light_switches):
            string = f"switch light {index} off" if light else f"switch light {index} on"
            self.light_switches_cb.addItem(string)
        self.light_switches_cb.activated.connect(self.on_light_switches_cb_activated)
        #
        # Colours list of buttons with element colours
        #
        self.coloured_elements_widget = QWidget(self)
        self.element_coloured_buttons = []
        self.element_coloured_hbox = QHBoxLayout()
        for el in self.species:
            r,g,b,a = self.element_colours[el]
            button = QPushButton(el)
            button.setStyleSheet(f"background-color:rgba( {r}, {g}, {b}, {a});")
            button.clicked.connect(self.on_coloured_element_clicked)
            self.element_coloured_buttons.append(button)
            self.element_coloured_hbox.addWidget(button)
        self.coloured_elements_widget.setLayout(self.element_coloured_hbox)
        self.coloured_elements_widget.setToolTip("Change the colour of the elements")
        #
        # Colours list of buttons with colours
        #
        self.coloured_buttons_widget = QWidget(self)
        self.bond_cell_background_arrow_buttons = []
        hbox = QHBoxLayout()
        button_colours = [  (r,g,b,a) for (r,g,b,a) in [self.settings["Background colour"], self.settings["Cell colour"], self.settings["Bond colour"],self.settings["Arrow colour"]] ]
        self.bond_cell_background_arrow_buttons = []
        for col, name in zip(button_colours, self.bond_cell_background_arrow_names):
            button = QPushButton(name)
            button.setStyleSheet("background-color:rgba( {}, {}, {}, {});".format(*col))
            button.clicked.connect(self.on_coloured_button_clicked)
            self.bond_cell_background_arrow_buttons.append(button)
            hbox.addWidget(button)
        self.coloured_buttons_widget.setToolTip("Change the colours of the background, cell, bonds and arrows")
        self.coloured_buttons_widget.setLayout(hbox)
        #
        # Add a tab widget for the settings
        #
        self.settingsTab = QTabWidget(self)
        self.settingsTab.addTab(self.coloured_elements_widget, "Elements")
        self.settingsTab.addTab(self.coloured_buttons_widget, "Colours")
        self.settingsTab.addTab(self.atom_scaling_sb, "Atom Size")
        self.settingsTab.addTab(self.super_cell_widget, "Super Cell")
        self.settingsTab.addTab(self.transform_tab_entry, "Transform")
        self.settingsTab.addTab(self.hkl_tab_entry, "Surface")
        self.settingsTab.addTab(self.light_switches_cb, "Lighting")
        self.settingsTab.addTab(self.toggles_cb, "Toggles")
        self.settingsTab.addTab(self.maximum_displacement_sb, "Displacement")
        self.settingsTab.addTab(self.bond_radius_sb, "Bond Radius")
        self.settingsTab.addTab(self.cell_radius_sb, "Cell Radius")
        self.settingsTab.addTab(self.arrow_radius_sb, "Arrow Radius")
        label = QLabel("Settings", self)
        form.addRow(label,self.settingsTab)
        #
        # Add buttons to save a png or movie
        #
        hbox = QHBoxLayout()
        self.filename_le = QLineEdit(self)
        self.filename_le.setToolTip("Give a file name in which to save the image.\nIf a return is pressed the filename is scanned to see whatyoe of file has been requested.")
        self.filename_le.setText(self.image_filename)
        self.filename_le.textChanged.connect(self.on_filename_le_changed)
        self.filename_le.returnPressed.connect(self.on_filename_le_return)
        self.png_filename_button = QPushButton("Save as png")
        self.png_filename_button.clicked.connect(self.on_filename_button_clicked)
        self.png_filename_button.setToolTip("Save the image as a png file with arrows")
        self.mp4_filename_button = QPushButton("Save as mp4")
        self.mp4_filename_button.clicked.connect(self.on_filename_button_clicked)
        self.mp4_filename_button.setToolTip("Save the animation as an mp4/avi file")
        self.gif_filename_button = QPushButton("Save as gif")
        self.gif_filename_button.clicked.connect(self.on_filename_button_clicked)
        self.gif_filename_button.setToolTip("Save the animation as a gif file")
        self.cif_filename_button = QPushButton("Save as cif or mcif")
        self.cif_filename_button.clicked.connect(self.on_filename_button_clicked)
        self.cif_filename_button.setToolTip("Save the structure as cif. Save the mode animation as mcif")
        label = QLabel("Image file name", self)
        label.setToolTip("Give a file name in which to save the image.\nIf a return is pressed the filename is scanned to see whatyoe of file has been requested.")
        hbox.addWidget(self.filename_le)
        hbox.addWidget(self.png_filename_button)
        hbox.addWidget(self.mp4_filename_button)
        hbox.addWidget(self.gif_filename_button)
        hbox.addWidget(self.cif_filename_button)
        form.addRow(label, hbox)
        #
        # Add the opengl figure to the bottom
        #
        self.opengl_widget = OpenGLWidget(self, debug=debug)
        self.opengl_widget.setFocusPolicy(Qt.StrongFocus)
        self.opengl_widget.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding))
        form.addRow(self.opengl_widget)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        self.debugger.print("Finished:: initialisation")

    def convert_transform(self):
        """Convert the transform stored as strings in settings to an np.array.

        Parameters
        ----------
        None

        Modifies
        --------
        transform

        """
        self.debugger.print("convert_transform")
        transform = []
        for row in self.settings["Transform"]:
            new_row = [ eval(col) for col in row ]
            transform.append(new_row)
        self.debugger.print("convert_transform", np.array(transform))
        return np.array(transform)
       
    def on_reset_button_clicked(self, item ):
        """Handle a push of the reset transform button.

        When the button is pushed the transform is set to unity.

        Parameters
        ----------
        item : bool
            ignored

        Returns
        -------
        None

        """
        self.debugger.print("on_reset_button_clicked",item)
        self.settings["Transform"] =  [ [ "1","0","0" ], [ "0","1","0" ], ["0","0","1"] ]
        self.debugger.print("on_guess_button_clicked transform",self.settings["Transform"])
        self.refreshRequired = True
        self.refresh()

    def on_guess_button_clicked(self, item):
        """Guess the primitive transform button click.

        The primitive transform is guessed from the spacegroup

        Parameters
        ----------
        item : bool
            ignored

        Returns
        -------
        None

        Modifies
        --------
        Modifies primitive unit cell

        """
        self.debugger.print("on_guess_button_clicked",item)
        transform = self.standard_cell.guess_primitive_transform()
        new = []
        for row in transform:
            new.append( [f"{col:.9f}" for col in row])
        self.settings["Transform"] = new
        self.debugger.print("on_guess_button_clicked transform",self.settings["Transform"])
        self.refreshRequired = True
        self.refresh()

    def on_edit_button_clicked(self, item):
        """Edit the transform.

        The transform window dialog is shown and any edits applied

        Parameters
        ----------
        item : bool
            ignored

        Returns
        -------
        None

        Modifies
        --------
        Modifies primitive unit cell

        """
        self.debugger.print("on_edit_button_clicked",item)
        transform_window = TransformWindow(self.settings["Transform"],
                                           debug=self.debugger.state())
        if transform_window.exec():
            # The 'Ok' button was pressed
            # get the new transform and replace the old one
            self.settings["Transform"] = transform_window.getTransform()
            self.debugger.print("on_edit_button_clicked transform",self.settings["Transform"])
            self.refreshRequired = True
            self.refresh()

    def on_primitive_button_clicked(self,item):
        """Handle a push of the primitive transform button.

        When the button is pushed the primitive transform is loaded from phonopy.yaml

        Parameters
        ----------
        item : bool
            ignored

        Returns
        -------
        None

        Modifies
        --------
        self.settings["Primitive transform"]

        """
        self.debugger.print("on_primitive_button_clicked",item)
        if self.reader is not None and self.reader.primitive_transformation is not None:
            self.settings["Primitive transform"] = self.reader.primitive_transformation
            self.debugger.print("on_primtive_button_clicked transform",self.settings["Primitive transform"])
        # change the primitive transform window popup
        self.refreshRequired = True
        self.refresh()

    def createHKLTabEntry(self):
        """Create a one line entry to allow an entry of the hkl for a surface.

        Parameters
        ----------
        None

        Returns
        -------
        a widget

        """
        self.debugger.print("CreateHKLTabEntry")
        # unique direction (hkl)
        self.h_sb, self.k_sb, self.l_sb = self.hkl_spin_boxes()
        # unique direction (uvw)
        self.u_sb, self.v_sb, self.w_sb = self.uvw_spin_boxes()
        container = QHBoxLayout()
        form = QFormLayout()
        hbox = QHBoxLayout()
        hbox.addWidget(self.h_sb)
        hbox.addWidget(self.k_sb)
        hbox.addWidget(self.l_sb)
        label = QLabel("    hkl:")
        form.addRow(label, hbox)
        container.addLayout(form)
        form = QFormLayout()
        hbox = QHBoxLayout()
        hbox.addWidget(self.u_sb)
        hbox.addWidget(self.v_sb)
        hbox.addWidget(self.w_sb)
        label = QLabel("    uvw:")
        form.addRow(label, hbox)
        container.addLayout(form)
        widget = QWidget(self)
        widget.setToolTip("Define a surface of the transformed cell using (hkl) & [uvw]")
        widget.setLayout(container)
        return widget

    def hkl_spin_boxes(self):
        """Create the spin boxes for hkl.

        Parameters
        ----------
        None

        Returns
        -------
        a set of spin boxes

        """
        infostring = "\nThe surface (hkl) is displayed using the keyboard shortcuts 's' or 'S'.\n's' shows the surface in the plane of the screen.\n'S' shows the surface edge on.\n [uvw] defines a perpendicular to the surface normal.\n(hkl) takes precedence unless it is (000).\nIf hkl is (000) 's' shows [uvw] perpendicular to the screen \nand 'S' shows [uvw] vertical on the screen"
        h,k,l = self.settings["hkl"]
        h_sb = QSpinBox(self)
        h_sb.setToolTip("Define the h dimension of the surface"+infostring)
        h_sb.setRange(-20,20)
        h_sb.setValue(h)
        h_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,0))
        k_sb = QSpinBox(self)
        k_sb.setToolTip("Define the k dimension of the surface"+infostring)
        k_sb.setRange(-20,20)
        k_sb.setValue(k)
        k_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,1))
        l_sb = QSpinBox(self)
        l_sb.setToolTip("Define the l dimension of the surface"+infostring)
        l_sb.setRange(-20,20)
        l_sb.setValue(l)
        l_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,2))
        return h_sb, k_sb, l_sb

    def uvw_spin_boxes(self):
        """Create the spin boxes for uvw.

        Parameters
        ----------
        None

        Returns
        -------
        a set of spin boxes

        """
        infostring = "\nThe surface (hkl) is displayed using the keyboard shortcuts 's' or 'S'.\n's' shows the surface in the plane of the screen.\n'S' shows the surface edge on.\n [uvw] defines a perpendicular to the surface normal.\n(hkl) takes presidence unless it is (000).\nIn this case 's' shows [uvw] perpendicular to the screen and 'S' show [uvw] vertical on the screen\n "
        h,k,l = self.settings["uvw"]
        h_sb = QSpinBox(self)
        h_sb.setToolTip("Define the u dimension of the surface"+infostring)
        h_sb.setRange(-20,20)
        h_sb.setValue(h)
        h_sb.valueChanged.connect(lambda x: self.on_uvw_sb_changed(x,0))
        k_sb = QSpinBox(self)
        k_sb.setToolTip("Define the v dimension of the surface"+infostring)
        k_sb.setRange(-20,20)
        k_sb.setValue(k)
        k_sb.valueChanged.connect(lambda x: self.on_uvw_sb_changed(x,1))
        l_sb = QSpinBox(self)
        l_sb.setToolTip("Define the w dimension of the surface"+infostring)
        l_sb.setRange(-20,20)
        l_sb.setValue(l)
        l_sb.valueChanged.connect(lambda x: self.on_uvw_sb_changed(x,2))
        return h_sb, k_sb, l_sb

    def on_uvw_sb_changed(self,value,index):
        """Handle a change in a value of u,v or, w.

        Parameters
        ----------
        value : int
            The value of u,v or w

        index : integer
          1, 2, or 3 for u, v or w

        Modifies
        --------
        settings["uvw"]

        """
        self.debugger.print("on_uvw_sb_changed", value, index)
        u,v,w = self.settings["uvw"]
        if index == 0:
            u = value
        elif index == 1:
            v = value
        else:
            w = value
        #
        # Do not allow 0,0,0
        #
        if abs(u) + abs(v) + abs(w) == 0:
            u = 1
            v = 0
            w = 0
            self.u_sb.blockSignals(True)
            self.v_sb.blockSignals(True)
            self.w_sb.blockSignals(True)
            self.u_sb.setValue(u)
            self.v_sb.setValue(v)
            self.w_sb.setValue(w)
            self.u_sb.blockSignals(False)
            self.v_sb.blockSignals(False)
            self.w_sb.blockSignals(False)
        self.settings["uvw"] = (u,v,w)
        self.opengl_widget.define_surface_orientations(self.transformed_cell,
                                                       self.settings["hkl"],
                                                       self.settings["uvw"])
        self.debugger.print("on_uvw_sb_changed uvw=", self.settings["uvw"])
        return 
    
    def on_hkl_sb_changed(self,value,index):
        """Handle a change in a value of h,k or, l.

        Parameters
        ----------
        value : int
            The value of h,k or l

        index : integer
          1, 2, or 3 for h, k or l

        Modifies
        --------
        settings["hkl"]

        """
        self.debugger.print("on_hkl_sb_changed", value, index)
        h,k,l = self.settings["hkl"]
        if index == 0:
            h = value
        elif index == 1:
            k = value
        else:
            l = value
        self.settings["hkl"] = (h,k,l)
        self.opengl_widget.define_surface_orientations(self.transformed_cell,
                                                       self.settings["hkl"],
                                                       self.settings["uvw"])
        self.debugger.print("on_hkl_sb_changed hkl=", self.settings["hkl"])
        return 
    

    def createTransformTabEntry(self):
        """Create a one line entry to modify the cell transformation matrix.

        Add buttons to reset, edit and if possible read the transformation matrix

        Parameters
        ----------
        None

        Returns
        -------
        a widget

        """
        self.debugger.print("createTransformTabEntry")
        hbox = QHBoxLayout()
        button = QPushButton("Reset transformation matrix")
        button.setToolTip("Reset the transformation to a unit matrix")
        button.clicked.connect(self.on_reset_button_clicked)
        button.clicked.connect(self.on_reset_button_clicked)
        hbox.addWidget(button)
        if self.reader is not None and self.reader.primitive_transformation is not None:
            # Give a button to load the primitive transformation 
            button = QPushButton("Load transformation from phonopy")
            button.setToolTip("Set the transformation to the primitive transformation matrix in phonopy.yaml")
            button.clicked.connect(self.on_primitive_button_clicked)
            hbox.addWidget(button)
        button = QPushButton("Guess transformation matrix")
        button.setToolTip("Guess the transformation and generate a new cell")
        button.clicked.connect(self.on_guess_button_clicked)
        hbox.addWidget(button)
        button = QPushButton("Edit transformation matrix")
        button.setToolTip("Edit the transformation and generate a new cell")
        button.clicked.connect(self.on_edit_button_clicked)
        hbox.addWidget(button)
        widget = QWidget(self)
        widget.setLayout(hbox)
        return widget

        
    def on_filename_le_return(self):
        """Handle the event triggered by the return key press within the filename input field.

        This method is typically invoked when the return (or enter) key is pressed while inputting a filename. It primarily invokes the filename button click event handler, simulating a click operation on the associated button.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("on filename le return pressed")
        self.on_filename_button_clicked(True)
        return

    def on_filename_le_changed(self,text):
        """Handle changes to the filename line edit.

        This method is triggered when the text in the filename line edit changes. It updates the internal state to reflect the new filename.

        Parameters
        ----------
        text : str
            The new text from the filename line edit.

        Returns
        -------
        None

        """        
        self.debugger.print("on filename le changed", text)
        self.image_filename = text
        return

    def on_filename_button_clicked(self,boolean):
        """Handle the action triggered by clicking the filename button.

        This method checks if the selected filename has a valid extension and exists. If any 
        checks fail, it alerts the user and aborts the operation. Otherwise, it proceeds to 
        save the current state or image of the OpenGL widget in the specified file format.

        Parameters
        ----------
        boolean : any
            This is a dummy parameter.

        Returns
        -------
        None

        Raises
        ------
        QMessageBox
            Shows an alert if the filename is empty or has an invalid extension, indicating that the file name for the image is not valid.

        Notes
        -----
        The valid file extensions are .mp4, .avi, .png, .gif, and .cif. The method modifies the `plot_type_index` 
        attribute of the caller object, saves the data or image into the specified file format, and ensures 
        the application's cursor is restored to its default state after the operation completes.

        """        
        self.debugger.print("on filename button clicked")
        #button = self.sender()
        #text = button.text()
        #
        # Check to see if the filename has been set
        #
        filename = self.image_filename
        root,extension = os.path.splitext(filename)
        if filename == "" or extension not in (".mp4", ".avi", ".png", ".gif", ".cif", ".mcif" ):
            self.debugger.print("Aborting on filename button clicked", filename)
            QMessageBox.about(self,"Image file name", "The file name for the image is not valid "+filename)
            return
        #
        # The filename will be opened in the directory of the output file which was read in mainTab
        #
        filename = os.path.join(self.notebook.mainTab.directory, filename)
        #
        # We need to remember the visualisation settings
        #
        old_plot_type = self.plot_type_index
        QApplication.setOverrideCursor(Qt.WaitCursor)
        if extension == ".png":
            self.opengl_widget.snapshot(filename)
        elif extension in ( ".avi", ".mp4", ".gif" ):
            self.opengl_widget.save_movie(filename)
        elif extension == ".cif":
            self.save_as_cif(filename)
        elif extension == ".mcif":
            self.save_mode_as_cif(filename)
        self.plot_type_index = old_plot_type
        QApplication.restoreOverrideCursor()
        return

    def on_super_cell_changed_a(self,newa):
        """Update the 'a' parameter of the 'Super Cell' setting and refresh the object state.

        This function takes a new value for 'a', updates the 'Super Cell' setting of the object
        accordingly (while keeping 'b' and 'c' parameters unchanged), marks the object as
        requiring refresh, and finally calls the refresh method.

        Parameters
        ----------
        newa : int
            The new value for the 'a' parameter to be set in the 'Super Cell' setting.

        Returns
        -------
        None

        """        
        a,b,c = self.settings["Super Cell"]
        a = newa
        self.settings["Super Cell"] =  [ a, b, c ]
        self.refreshRequired = True
        self.refresh()
        return

    def on_super_cell_changed_b(self,newb):
        """Update the 'b' component of the 'Super Cell' setting and trigger a refresh.

        This function takes a new value for 'b', updates the 'Super Cell' setting of the object
        accordingly (while keeping 'a' and 'c' parameters unchanged), marks the object as
        requiring refresh, and finally calls the refresh method.

        Parameters
        ----------
        newb : int
            The new value for the 'b' component of the 'Super Cell'.

        Returns
        -------
        None

        """        
        a,b,c = self.settings["Super Cell"]
        b = newb
        self.settings["Super Cell"] =  [ a, b, c ]
        self.refreshRequired = True
        self.refresh()
        return

    def on_super_cell_changed_c(self,newc):
        """Update the 'c' component of the 'Super Cell' setting and refresh the object.

        This function takes a new value for 'c', updates the 'Super Cell' setting of the object
        accordingly (while keeping 'a' and 'b' parameters unchanged), marks the object as
        requiring refresh, and finally calls the refresh method.

        Parameters
        ----------
        newc : int
            The new value to set for the 'c' component of the 'Super Cell' setting.

        Returns
        -------
        None

        """        
        a,b,c = self.settings["Super Cell"]
        c = newc
        self.settings["Super Cell"] =  [ a, b, c ]
        self.refreshRequired = True
        self.refresh()
        return

    def on_coloured_element_clicked(self,boolean):
        """Handle a click event on an element with colour functionality.

        A colour dialogue is invoked and the chosen colour is assigned to the element in question
        The settings dictionary is updated with the new colours and a refresh requested

        Parameters
        ----------
        boolean : any
            This parameter is a dummy parameters

        Returns
        -------
        None

        """        
        self.debugger.print("on coloured elements clicked")
        button = self.sender()
        text = button.text()
        colourChooser = QColorDialog()
        if colourChooser.exec_() != QColorDialog.Accepted:
            return
        colour = colourChooser.currentColor()
        rgba = [ colour.red(), colour.green(), colour.blue(), colour.alpha() ]
        self.element_colours[text] = rgba
        #
        # This is a setting based on the element, there is a colour for each element
        #
        self.settings["Element colours"] = [ self.element_colours[el] for el in self.species ]
        self.refreshRequired = True
        self.refresh()
        return

    def on_coloured_button_clicked(self,boolean):
        """Handle the event when a coloured button is clicked.

        This function triggers a color dialog upon the click of a designated button. It captures the selected color and applies it to a specific setting based on the button's label (e.g., Background, Cell, Bonds, Arrows). The settings are modified to include the chosen color, and a refresh is triggered to apply the changes.

        Parameters
        ----------
        boolean : bool
            This parameter is not used in the function.

        Returns
        -------
        None

        """        
        self.debugger.print("on coloured button clicked")
        button = self.sender()
        text = button.text()
        colourChooser = QColorDialog()
        if colourChooser.exec_() != QColorDialog.Accepted:
            return
        colour = colourChooser.currentColor()
        rgba = [ colour.red(), colour.green(), colour.blue(), colour.alpha() ]
        if text == "Background":
            self.settings["Background colour"] = rgba
        elif text == "Cell":
            self.settings["Cell colour"] = rgba
        elif text == "Bonds":
            self.settings["Bond colour"] = rgba
        elif text == "Arrows":
            self.settings["Arrow colour"] = rgba
        self.refreshRequired = True
        self.refresh()
        return

    def on_toggles_cb_activated(self, index):
        """Activate or deactivate the various toglles based on the current state and updates the GUI accordingly.

        Parameters
        ----------
        index : int
            Index of the toggle in `self.toggle_name` list

        Returns
        -------
        None

        """        
        self.debugger.print("on_toggles_cb_activated")
        self.settings["Toggle states"][index] = not self.settings["Toggle states"][index]
        state = self.settings["Toggle states"][index]
        string = f"{self.toggle_names[index]} is on" if state else f"{self.toggle_names[index]} is off"
        self.toggles_cb.setItemText(index,string)
        self.refreshRequired = True
        self.refresh()
        self.plot()
        return

    def on_light_switches_cb_activated(self, index):
        """Activate or deactivate the light switch based on the current state and updates the GUI accordingly.

        Parameters
        ----------
        index : int
            Index of the light switch in `self.light_switches` list that needs to be toggled.

        Returns
        -------
        None

        """        
        self.debugger.print("on_light_switches_cb_activated")
        self.light_switches[index] = not self.light_switches[index]
        string = f"switch light {index} off" if self.light_switches[index] else f"switch light {index} on"
        self.light_switches_cb.setItemText(index,string)
        self.opengl_widget.defineLights()
        self.calculate()
        self.plot()
        return

    def on_maximum_displacement_changed(self,value):
        """Handle changes to the displacement setting.

        This method updates the 'Displacement' setting, recalculates based on the new value, and then replots the relevant data or figures.

        Parameters
        ----------
        value : float
            The new value for the displacement setting.

        Returns
        -------
        None

        """        
        self.debugger.print("on maximum_displacement changed ", value)
        self.settings["Maximum displacement"] = value
        self.calculate()
        self.plot()
        return

    def on_atom_scaling_changed(self,value):
        """Handle changes to atom scaling settings.

        This method updates the atom scaling setting based on a new value,
        recalculates relevant data, and triggers a plot update.

        Parameters
        ----------
        value : float
            The new value for atom scaling.

        Returns
        -------
        None

        """        
        self.debugger.print("on atom_scaling changed ", value)
        self.settings["Atom scaling"] = value
        self.calculate()
        self.plot()
        return

    def on_arrow_radius_changed(self,value):
        """Handle the event where the arrow's radius value is changed.

        This method updates the setting for the arrow's radius, recalculates, and replots based on the new value.

        Parameters
        ----------
        value : float
            The new value for the arrow radius.

        Returns
        -------
        None

        """        
        self.debugger.print("on arrow_radius changed ", value)
        self.settings["Arrow radius"] = value
        self.calculate()
        self.plot()
        return

    def on_cell_radius_changed(self,value):
        """Handle the event when the cell radius setting is changed.

        This method updates the 'Cell radius' setting with the new value, recalculates 
        the necessary data, and then replots the results.

        Parameters
        ----------
        value : float
            The new value for the cell radius.

        Returns
        -------
        None

        """        
        self.debugger.print("on cell_radius changed ", value)
        self.settings["Cell radius"] = value
        self.calculate()
        self.plot()
        return

    def on_bond_radius_changed(self,value):
        """Handle bond radius change events.

        This method updates the 'Bond radius' setting with the new value, recalculates 
        the necessary data, and then replots the results.

        Parameters
        ----------
        value : float
            The new value for the bond radius.

        Returns
        -------
        None

        """        
        self.debugger.print("on_bond_radius_changed")
        self.settings["Bond radius"] = value
        self.calculate()
        self.plot()
        return

    def on_selected_mode_changed(self):
        """Handle the changes when a new vibrational mode is selected in the UI.

        This function updates the application state based on the newly selected vibrational mode. It
        updates the displayed frequency, removes existing arrows from the OpenGL widget, 
        computes maximum displacement for scaling, and plots new arrows with updated settings. 
        It also invokes the calculate and plot functions to reflect changes in the UI.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("on_selected_mode_changed")
        if self.frequencies_cm1 is None or len(self.frequencies_cm1) < 1:
           self.debugger.print("on_selected_mode_changed aborting")
           return
        self.debugger.print("on selected_mode change mode was ", self.settings["Selected mode"])
        self.settings["Selected mode"] = self.selected_mode_sb.value()
        self.debugger.print("on selected_mode change mode is now ", self.settings["Selected mode"])
        self.frequency_le.setText("{:.5f}".format(self.frequencies_cm1[self.settings["Selected mode"]-1]))
        self.opengl_widget.deleteArrows()
        maxR = np.max( np.abs(np.array(self.UVW[self.settings["Selected mode"]-1])))
        if maxR < 1.0E-8:
            maxR = 1.0E-8
            self.plot_type_index = 2
        arrow_scaling = self.settings["Maximum displacement"] / maxR
        for uvw in self.UVW[self.settings["Selected mode"]-1]:
            self.opengl_widget.addArrows( self.settings["Arrow colour"],self.settings["Arrow radius"], uvw, arrow_scaling )
        self.calculate()
        self.plot()
        return

    def on_plottype_cb_changed(self, index):
        """Handle the change in plot type selection.

        This method updates the plot type index based on the user's selection and then
        re-plots the data according to the new plot type.
        Plot types can be: 'Animation', 'Arrows', or 'No arrows or animation'.

        Parameters
        ----------
        index : int
            The index of the newly selected plot type.

        Returns
        -------
        None

        """        
        self.debugger.print("on_plottype_cb_changed")
        self.plot_type_index = index
        self.debugger.print("Plot type index changed to ", self.plot_type_index)
        self.plot()
        return

    def calculate(self):
        """Perform calculations related to the ViewerTab instance, including processing program, file name, calculating frequencies, super cells, normal modes, bonds, center of mass, bounding box, element names, species, covalent radii, and updating the UI with calculated values.

        In addition, calculate the transformed cell if a transformation is given.
        The supercell is calculated from the dft cell.

        Parameters
        ----------
        None 

        Returns
        -------
        None

        """        
        self.debugger.print("Start:: calculate")
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings["Program"]
        filename = self.notebook.mainTab.getFullFileName()
        self.debugger.print("calculate program file name",program, filename)
        self.reader = self.notebook.reader
        if self.reader is None:
            self.debugger.print("Finished:: calculate - reader is None")
            return
        if program == "":
            self.debugger.print("Finished:: calculate - program is blank")
            return
        if filename == "":
            self.debugger.print("Finished:: calculate - filename is blank")
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Assemble the settingsTab settings
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        self.selected_mode_sb.setRange(1,len(self.frequencies_cm1))
        #
        # Define a standard cell as read by the reader
        #
        self.standard_cell = copy.copy(self.reader.get_unit_cell())
        if self.standard_cell is None:
            return
        #
        # Transform the standard cell 
        # this involves folding the cell into itself and applying the transformation
        # In the process there is a reordering of the atoms which is kept track of in the transformed cell
        # NOTE: The transformed cell may have fewer atoms than the standard cell
        #
        cell = copy.copy(self.standard_cell)
        transform = self.convert_transform()
        self.transformed_cell = PrimitiveCell(cell,transformation=transform)
        self.opengl_widget.define_orientations(self.transformed_cell)
        self.opengl_widget.define_surface_orientations(self.transformed_cell,
                                                       self.settings["hkl"],
                                                       self.settings["uvw"])
        self.opengl_widget.set_orientation(self.opengl_widget.orientation)
        #
        # Calculate the whole molecule content of the DFT cell
        #
        tolerance = self.notebook.analysisTab.settings["Bonding tolerance"]
        scale     = self.notebook.analysisTab.settings["Covalent radius scaling"]
        radii     = self.notebook.analysisTab.element_radii
        cell.calculate_molecular_contents( scale=scale,
                                           tolerance=tolerance,
                                           radii=radii)
        atom_masses = cell.get_atomic_masses()
        mass_weighted_normal_modes = self.reader.mass_weighted_normal_modes
        normal_modes = Calculator.normal_modes(atom_masses,mass_weighted_normal_modes)
        # Generate a super cell
        imageSpecifier = self.settings["Super Cell"]
        self.super_cell = SuperCell(cell,imageSpecifier)
        self.normal_modes = self.super_cell.calculateNormalModes(normal_modes)
        self.bonds = self.super_cell.calculateBonds()
        self.snapshot_number = 0
        self.nbonds = len(self.bonds)
        self.XYZ = self.super_cell.calculateXYZ()
        self.natoms = len(self.XYZ)
        self.number_of_modes = len(self.normal_modes)
        # get the cell edges for the bounding box, shifted to the centre of mass origin
        totalMass,centreOfMassXYZ,centreOfMassABC = self.super_cell.calculateCentreOfMass(output=all)
        self.cell_corners,self.cell_edges,self.cell_labels = self.transformed_cell.getBoundingBox(originABC=[0.5,0.5,0.5])
        #
        # self.element_names is a list of element names for each atom
        # self.species is just a unique list of species
        #
        self.element_names = self.super_cell.getElementNames()
        self.species = self.reader.getSpecies()
        covalent_radii = self.notebook.analysisTab.element_radii
        #
        # Overwrite the element colours from the settings entry
        #
        if self.settings["Element colours"] is not None:
            self.debugger.print("refresh - processing colours ",self.species)
            for el,colour in zip(self.species,self.settings["Element colours"]):
                self.element_colours[el] = colour
        #
        # Define the radius and colour for every atom from the element_names
        #
        self.radii = [self.settings["Atom scaling"]*covalent_radii[el] for el in self.element_names ]
        self.colours = [ self.element_colours[el] for el in self.element_names ]
        #
        # reorder the displacement info in the normal modes into U,V and W lists
        # Using deque here rather than a simple list as the memory allocation doesn't have to be contiguous
        #
        self.UVW.clear()
        for displacements in self.normal_modes:
            uvw = deque()
            for i in range(0,len(displacements),3):
                uvw.append( displacements[i:i+3] )
            self.UVW.append(uvw)
        # CalculatePhasePositions stores all the sphere and bond information
        self.calculatePhasePositions()
        # Add the arrows
        self.debugger.print("calculate: Selected mode",self.settings["Selected mode"])
        maxR = np.max( np.abs(np.array(self.UVW[self.settings["Selected mode"]-1])))
        if maxR < 1.0E-8:
            maxR = 1.0E-8
            self.plot_type_index = 2
        arrow_scaling = self.settings["Maximum displacement"] / maxR
        self.opengl_widget.deleteArrows()
        for uvw in self.UVW[self.settings["Selected mode"]-1]:
            self.opengl_widget.addArrows( self.settings["Arrow colour"],self.settings["Arrow radius"], uvw, arrow_scaling )
        self.opengl_widget.setRotationCentre( centreOfMassXYZ )
        self.opengl_widget.setImageSize()
        QApplication.restoreOverrideCursor()
        self.debugger.print("Finished:: calculate")
        return

    def setColour(self, element, colour):
        """Set the colour of a specified element in the interface.

        Parameters
        ----------
        element : str
            The element whose colour is to be set. Can be 'Background', 'Cell', 'Arrow', or any other custom element defined in `element_colours`.
        colour : str
            The new colour to be set for the element. This should be a string describing the colour.

        Returns
        -------
        None

        Examples
        --------
        ```python
        myObject.setColour('Background', '#FFFFFF')  # Sets the background colour to white
        myObject.setColour('Cell', 'red')            # Sets the cell colour to red
        myObject.setColour('customElement', 'blue')  # Sets a custom element's colour to blue
        ```

        """        
        self.debugger.print("setcolour")
        if element in ( "Background" , "background" ):
            self.settings["Background colour"] = colour
        elif element in ( "Cell" , "cell" ):
            self.settings["Cell colour"] = colour
        elif element in ( "Arrow" , "arrow" ):
            self.settings["Arrow colour"] = colour
        else:
            self.element_colours[element] = colour
        self.plot()
        return

    def calculatePhasePositions(self):
        """Calculate and update the phase positions for the normal mode of a molecular structure visualization.

        This function calculates the phase positions for all atoms in a molecular structure over a range of phase steps. It updates the visualization by adding spheres and cylinders to represent atoms and bonds in different phase positions. The phase steps are adjusted to be odd in number, ensuring symmetry around zero phase. The function directly modifies the OpenGL widget used for visualization by deleting existing graphical representations and creating new ones based on the calculated positions.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("calculatePhasePositions")
        # we need the number of phase steps to be odd
        if self.settings["Number of phase steps"]%2 == 0:
            self.settings["Number of phase steps"] += 1
        UVW = np.array( self.UVW[self.settings["Selected mode"]-1] )
        maxR = np.amax(np.abs(UVW))
        if maxR < 1.0E-8:
            maxR = 1.0E-8
            self.plot_type_index = 2
        self.scale_vibrations = self.settings["Maximum displacement"] / maxR
        self.newXYZ       = np.zeros( (self.settings["Number of phase steps"], self.natoms, 3) )
        n2 = int(self.settings["Number of phase steps"]/2)
        delta = 1.0 / float(n2)
        phases = np.arange(-1.0, 1.0+delta-1.0E-10, delta)
        for phase_index,phase in enumerate(phases):
            self.newXYZ[phase_index] = self.XYZ+phase*self.scale_vibrations*UVW
        # end for phase_index
        #
        # Store the results in the opengl widget
        #
        self.opengl_widget.deleteSpheres()
        self.opengl_widget.deleteCylinders()
        self.opengl_widget.createArrays(len(phases))
        self.debugger.print("calculatePhasePositions - adding spheres and cylinders")
        if self.get_toggle_state("Show debug information"):
            self.debugger = Debug(True,"ViewerTab")
            self.debugger.print("Debug is on:: toggle")
            self.debug = True
        else:
            self.debugger.print("Debug is off:: toggle")
            self.debugger = Debug(False,"ViewerTab")
            self.debug = False
        if self.get_toggle_state("Show orientation"):
            self.opengl_widget.show_orientation = True
        else:
            self.opengl_widget.show_orientation = False
        for phase_index in range(len(phases)):
            for col, rad, xyz in zip(self.colours, self.radii, self.newXYZ[phase_index]):
                if self.get_toggle_state("Show atoms"):
                    self.opengl_widget.addSphere(col, rad, xyz, phase=phase_index )
            for p,l in zip(self.cell_corners,self.cell_labels):
                if self.get_toggle_state("Show cell"):
                    self.opengl_widget.addSphere(self.settings["Cell colour"],
                                                 self.settings["Cell radius"],
                                                 p,
                                                 phase=phase_index )
                if self.get_toggle_state("Show cell labels") and self.get_toggle_state("Show cell"):
                    self.opengl_widget.addText(l, 
                                               self.settings["Text colour"], 
                                               self.settings["Text size"], 
                                               p, 
                                               phase=phase_index)
            if self.get_toggle_state("Show bonds"):
                for bond in self.bonds:
                    i,j = bond
                    self.opengl_widget.addCylinder(self.settings["Bond colour"],
                                                   self.settings["Bond radius"],
                                                   self.newXYZ[phase_index,i],
                                                   self.newXYZ[phase_index,j],
                                                   phase=phase_index)
            if self.get_toggle_state("Show cell"):
                for p1,p2 in self.cell_edges:
                    self.opengl_widget.addCylinder(self.settings["Cell colour"], self.settings["Cell radius"], p1, p2, phase=phase_index)
        self.debugger.print("calculatePhasePositions - exiting")
        return

    def save_as_cif(self,filename):
        """Save the crystallographic information file (CIF) for the structure phases.

        The cell is wrtten out

        Parameters
        ----------
        filename : str
            The name of the file where the CIF data will be saved.

        Returns
        -------
        None

        """        
        self.standard_cell.write_cif(filename=filename)

    def save_mode_as_cif(self,filename):
        #
        # Write a single cif file containing all the phase information for displacement along a mode
        # First transform to abc coordinates from xyz
        #
        """Save the crystallographic information file (CIF) for different phases.

        A copy is made of the current unit cell and a cell is wrtten out for each phase of the vibration

        Parameters
        ----------
        filename : str
            The name of the file where the CIF data will be saved.

        Returns
        -------
        None

        Notes
        -----
        This method relies on the `self.standard_cell` object, which should have a method `set_xyz_coordinates` and `write_cif`. It also uses `self.settings` to retrieve the number of phase steps and the selected mode for CIF generation. The function iterates over a calculated range of phases, updates the unit cell coordinates for each phase, and writes the CIF data to the provided filename. The function assumes that `self.newXYZ` is an iterable object containing new XYZ coordinates for each phase. The CIF files are saved with a description that includes the mode and phase.

        """        
        import copy
        unitcell = copy.deepcopy(self.standard_cell)
        n2 = int(self.settings["Number of phase steps"]/2)
        delta = 1.0 / float(n2)
        phases = np.arange(-1.0, 1.0+delta-1.0E-10, delta)
        with open(filename,"w") as fd:
            for phase_index,phase in enumerate(phases):
                description = "mode_"+str(self.settings["Selected mode"]-1)+"_phase_"+str(phase)
                unitcell.set_xyz_coordinates(self.newXYZ[phase_index])
                unitcell.write_cif(description=description,filedescriptor=fd)

    def plot(self):
        """Plot data based on the instance's current plot type.
    
        This method selects a plotting routine based on the value of `self.plot_type_index`.
        It does an early return if `self.reader` is None, indicating there might be no data
        to plot. It logs the method's progress using a debugger's print function.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("Start:: plot")
        if self.reader is None:
            self.debugger.print("Finished:: plot reader is None")
            return
        if self.plot_type_index == 0:
            self.plot_animation()
        elif self.plot_type_index == 1:
            self.plot_arrows()
        else:
            self.plot_none()
        self.debugger.print("Finished:: plot")

    def plot_none(self):
        """Hides arrow visuals, stops any ongoing animation, and updates the Open GL widget.

        This method is intended to reset or clear the current state of the Open GL widget by disabling any arrow visuals, stopping any animations that might be running, and then updating the widget to reflect these changes.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("Start:: plot_animation")
        self.opengl_widget.showArrows(False)
        self.opengl_widget.stopAnimation()
        self.opengl_widget.update()
        self.debugger.print("Finished:: plot_animation")
        return

    def plot_animation(self):
        """Initiate and display an animation of a vibrational mode in an OpenGL widget.

        This method controls the display of a predefined animation within 
        an OpenGL widget. It hides any arrows, triggers a visual update, 
        and starts the animation process. Debug messages are logged at the 
        start and completion of the animation process.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("Start:: plot_animation")
        self.opengl_widget.showArrows(False)
        self.opengl_widget.update()
        self.opengl_widget.startAnimation()
        self.debugger.print("Finished:: plot_animation")
        return

    def plot_arrows(self):
        """Plot arrows showing the vibrational mode of a molecule.

        This method enables the display of arrows on an associated OpenGL widget, halts any ongoing animations, and triggers a refresh of the widget to ensure the changes are visualized.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("Start:: plot_arrows")
        self.opengl_widget.showArrows(True)
        self.opengl_widget.stopAnimation()
        self.opengl_widget.update()
        self.debugger.print("Finished:: plot_arrows")
        return

    def requestRefresh(self):
        """Initiate a refresh request.

        This method sets an internal flag to indicate that a refresh is required.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        self.debugger.print("Start:: requestRefresh")
        self.refreshRequired = True
        self.debugger.print("Finished:: requestRefresh")

    def refresh(self,force=False):
        """Refresh the state of the object, optionally forcing the refresh.

        This method updates the object's state by recalculating and refreshing its components. It checks whether a refresh is necessary or if it has been forced. It iterates through child widgets to block signals, updates settings from a notebook, adjusts visual elements based on these settings, and ultimately recalculates and replots data.

        Parameters
        ----------
        force : bool, optional
            If True, the refresh proceeds regardless of whether it is deemed necessary based on the object state. Default is False.

        Returns
        -------
        None

        Notes
        -----
        - The function begins by checking if a refresh is required or has been forced. If neither condition is met, it aborts the refresh.
        - It temporarily changes the cursor to a waiting cursor to indicate processing.
        - It blocks signals from all child widgets to prevent any unintended interaction during the refresh process.
        - It updates various widget properties and settings based on attributes of the notebook associated with the object.
        - The function finally recalculates and re-plots based on the latest data and settings, then restores the standard cursor and marks the refresh as complete.
        - Debug statements are interspersed throughout for tracking the refresh process's progress.

        """        
        self.debugger.print("Start:: refresh")
        if not self.refreshRequired and not force:
            self.debugger.print("Finished:: refresh aborted",self.refreshRequired,force)
            return
        self.debugger.print("refresh widget",force)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        if "vesta" in self.settings["Element palette"].lower():
            self.element_colours = vesta_elemental_colours
        elif "jmol" in self.settings["Element palette"].lower():
            self.element_colours = jmol_elemental_colours
        else:
            self.element_colours = elemental_colours
        self.calculate()
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        if self.frequencies_cm1 is not None and len(self.frequencies_cm1) > 1:
            self.selected_mode_sb.setRange(1,len(self.frequencies_cm1))
            self.selected_mode_sb.setValue(self.settings["Selected mode"])
            self.frequency_le.setText("{:.5f}".format(self.frequencies_cm1[self.settings["Selected mode"]-1]))
        self.debugger.print("refresh: selected mode is now",self.selected_mode_sb.value())
        for index,(toggle,state) in enumerate(zip(self.toggle_names,self.settings["Toggle states"])):
            string = f"{toggle} is on" if state else f"{toggle} is off"
            self.toggles_cb.setItemText(index,string)
        h,k,l = self.settings["hkl"]
        self.h_sb.setValue(h)
        self.k_sb.setValue(k)
        self.l_sb.setValue(l)
        u,v,w = self.settings["uvw"]
        self.u_sb.setValue(u)
        self.v_sb.setValue(v)
        self.w_sb.setValue(w)
        self.atom_scaling_sb.setValue(self.settings["Atom scaling"])
        self.bond_radius_sb.setValue(self.settings["Bond radius"])
        self.cell_radius_sb.setValue(self.settings["Cell radius"])
        self.maximum_displacement_sb.setValue(self.settings["Maximum displacement"])
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        for index,light in enumerate(self.light_switches):
            string = f"switch light {index} off" if light else f"switch light {index} on"
            self.light_switches_cb.setItemText(index,string)
        #
        # Colours list of buttons with element colours
        #
        self.reader = self.notebook.mainTab.reader
        if self.reader is None:
            self.debugger.print("reader is none ")
            self.element_names = []
            self.species = []
        else:
            if self.standard_cell is None:
                return
            self.element_names = self.standard_cell.element_names
            self.species = self.reader.getSpecies()
            self.debugger.print("refresh - species ",self.species)
            self.settings["Element colours"] = [ self.element_colours[el] for el in self.species ]
        count = self.element_coloured_hbox.count()
        if count == 0:
            self.debugger.print("initialising element colours widget")
            self.element_coloured_buttons = []
            for el in self.species:
                r,g,b,a = self.element_colours[el]
                button = QPushButton(el)
                button.setStyleSheet(f"background-color:rgba( {r}, {g}, {b}, {a});")
                button.clicked.connect(self.on_coloured_element_clicked)
                self.element_coloured_buttons.append(button)
                self.element_coloured_hbox.addWidget(button)
        else:
            self.debugger.print("update element colours widget")
            for el,button in zip(self.species,self.element_coloured_buttons):
                r,g,b,a = self.element_colours[el]
                button.setStyleSheet(f"background-color:rgba( {r}, {g}, {b}, {a});")
        #
        # Colours list of buttons with colours
        #
        button_colours = [  (r,g,b,a) for (r,g,b,a) in [self.settings["Background colour"], self.settings["Cell colour"], self.settings["Bond colour"],self.settings["Arrow colour"]] ]
        for col, button in zip(button_colours, self.bond_cell_background_arrow_buttons):
            button.setStyleSheet("background-color:rgba( {}, {}, {}, {});".format(*col))
        #
        # Filename
        #
        self.filename_le.setText(self.image_filename)
        #
        # Unlock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.plot()
        self.refreshRequired = False
        QApplication.restoreOverrideCursor()
        self.debugger.print("Finished:: refresh")
        return

    def get_toggle_state(self,toggle):
        """Return the state of the given toggle.

        Parameters
        ----------
        toggle : str
            The toggle name

        Returns
        -------
        bool

        """
        self.debugger.print("get_toggle_state: ")
        if toggle not in self.toggle_names:
            return False
        index = self.toggle_names.index(toggle)
        self.debugger.print("get_toggle_state: {index}")
        return self.settings["Toggle states"][index]

class TransformWindow(QDialog):
    """A GUI window for displaying and editing the transform matrix.

    Parameters
    ----------
    transform : 3x3 array of floats
        The current transformation matrix
    message : str
        Message to displayed
    parent : QWidget, optional
        The parent widget of this window. The default is None.
    debug : bool, optional
        Flag to enable or disable debug mode. The default is False.

    Attributes
    ----------
    transform : 3x3 array of strs
        The transform matrix
    message : str
        The message displayed at the top of the window.
    layout : QVBoxLayout
        The main layout of the window.
    buttonBox : QDialogButtonBox
        The dialog box containing Ok and Cancel buttons.

    Methods
    -------
    reDrawTransformTable
        draw the transform table
    getTransform
        Return the transform matrix
    on_transform_table_itemChanged
        Handles a change in the transform table

    """

    def __init__(self, transform, message = "Cell Transform", parent=None, debug=False ):
        """Initialize a Transform Window instance.

        Parameters
        ----------
        transform : a 3x3 list of strs
            The transformation matrix
        message : str, optional
            Custom message to display in the window. The default is an empty string.
        parent : QWidget, optional
            The parent widget of this window. The default is None.
        debug : bool, optional
            Flag to enable or disable debug mode. The default is False.

        Notes
        -----
        This constructor initializes the showTransformationWindow by setting up the UI elements
        including buttons and the layout. It also configures the debug mode according
        to the provided argument.

        """        
        super().__init__(parent)
        self.debugger = Debug(debug,"TransformWindow")
        self.debugger.print("Start:: initialiser")
        # Set up the buttons of the button box
        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        # Initialize the transform as a 3x3 array of strings
        self.transform_as_str = transform
        self.debugger.print("transform:", self.transform_as_str)
        self.message = message
        # The dialog will have a vertical layout
        self.layout = QVBoxLayout(self)
        # Create the transfom table widget
        self.transform_table = QTableWidget()
        self.transform_table.setToolTip("Transformation for a cell")
        self.transform_table.setShowGrid(True)
        self.transform_table.setRowCount(3)
        self.transform_table.setColumnCount(3)
        self.reDrawTransformTable()
        self.transform_table.itemChanged.connect(self.on_transform_table_itemChanged)
        self.layout.addWidget(self.transform_table)
        # Add the button box
        self.layout.addWidget(self.buttonBox)
        self.debugger.print("Finished:: initialiser")

    def on_transform_table_itemChanged(self,item):
        """Handle a change to the transform table.

        The item which has changed is ignored and the whole table is rebuilt

        Parameters
        ----------
        item : the item changed
            The item which has changed

        Returns
        -------
        None

        Modifies
        --------
        self.settings["Transform"]

        """
        self.debugger.print("on_transform_table_itemChanged: ",item)
        if self.transform_table is None:
            return
        for col in range(3):
            for row in range(3):
                self.transform_as_str[row][col] = self.transform_table.item(row,col).text()
        self.debugger.print("New transform", self.transform_as_str)
        return

    def reDrawTransformTable(self):
        """Redraw the 3x3 table representing the cell transformation.

        The table is generated from the strings stored in self.transform_as_str
        If a number is stored it is converted to a string.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Modifies
        --------
        transform_table

        """
        self.debugger.print("reDrawPrimtiveTable: ")
        tw = self.transform_table
        for col in range(3):
            for row in range(3):
                tw.setItem(row,col,QTableWidgetItem( self.transform_as_str[row][col] ) )
        return

    def getTransform(self):
        """Return the primtive transformation matrix.

        The transformation is stored as a 3x3 list of strings.
        The strings are treated as python commands to return a value.:w


        Parameters
        ----------
        None

        Returns
        -------
        3x3 array of floats

        """
        self.debugger.print("getTransform: ")
        return self.transform_as_str
