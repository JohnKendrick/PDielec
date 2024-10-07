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
import os
from collections import deque

import numpy as np
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication,
    QColorDialog,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from PDielec.Constants import elemental_colours
from PDielec.GUI.OpenGLWidget import OpenGLWidget

# Need the SuperCell class
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
    unit_cell, super_cell, cell_edges, cell_corners : Various
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
    on_filename_button_clicked
        Handles a click on the file name button
    on_filename_le_changed
        Handles a change to the file name
    on_filename_le_return
        Handles a return press in the file name widget
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
        global debugger
        debugger = Debug(debug,"ViewerTab")
        debugger.print("Start:: initialisation")
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
        self.settings["Background colour"]     = [ 120, 120, 120, 255 ]
        self.settings["Arrow colour"]          = [   0, 255,   0, 255 ]
        self.settings["Arrow radius"]          = 0.07
        self.settings["Number of phase steps"] = 41
        self.settings["Super Cell"] =  [ 1, 1, 1 ]
        self.light_switches = [False]*8
        self.light_switches[0] = True
        self.light_switches[1] = True
        self.plot_types = ["Animation","Arrows","No arrows or animation"]
        self.plot_type_index = 1
        self.number_of_molecules = 0
        self.unit_cell = None
        self.super_cell = None
        self.cell_edges = None
        self.cell_corners = None
        # element_colours is a dictionary
        self.element_colours = elemental_colours
        self.settings["Element colours"] = None
        self.element_names = []
        self.species       = []
        self.element_coloured_buttons = []
        self.bond_cell_background_arrow_buttons = []
        self.bond_cell_background_arrow_names = ["Background","Cell","Bonds","Arrows"]
        self.image_filename = ""
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
        # The super-cell
        #
        self.super_cell_widget = QWidget(self)
        self.super_cell_widget.setToolTip("Generate a super-cell")
        self.super_cell_hbox = QHBoxLayout()
        self.super_cell_spinBoxes = []
        super_cell_changed = [self.on_super_cell_changed_a, self.on_super_cell_changed_b, self.on_super_cell_changed_c]
        super_cell_tooltip = ["Size in a", "Size in b", "Size in c"]
        for change_function,tip in zip(super_cell_changed,super_cell_tooltip):
            spinBox = QSpinBox(self)
            spinBox.setToolTip(tip)
            spinBox.setRange(1,5)
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
        # Add a comb box to select which type of plot
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
        self.coloured_buttons_widget.setLayout(hbox)
        #
        # Add a tab widget for the settings
        #
        self.settingsTab = QTabWidget(self)
        #self.settingsTab.currentChanged.connect(self.on_settingsTab_currentChanged)
        self.settingsTab.addTab(self.coloured_elements_widget, "Elements")
        self.settingsTab.addTab(self.coloured_buttons_widget, "Colours")
        self.settingsTab.addTab(self.atom_scaling_sb, "Atom Scaling")
        self.settingsTab.addTab(self.super_cell_widget, "Super Cell")
        self.settingsTab.addTab(self.light_switches_cb, "Lighting")
        self.settingsTab.addTab(self.maximum_displacement_sb, "Maximum Displacement")
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
        debugger.print("Finished:: initialisation")

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
        debugger.print("on filename le return pressed")
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
        debugger.print("on filename le changed", text)
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
        debugger.print("on filename button clicked")
        #button = self.sender()
        #text = button.text()
        #
        # Check to see if the filename has been set
        #
        filename = self.image_filename
        root,extension = os.path.splitext(filename)
        if filename == "" or extension not in (".mp4", ".avi", ".png", ".gif", ".cif", ".mcif" ):
            debugger.print("Aborting on filename button clicked", filename)
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
        debugger.print("on coloured elements clicked")
        button = self.sender()
        text = button.text()
        colourChooser = QColorDialog()
        if colourChooser.exec_() != QColorDialog.Accepted:
            return
        colour = colourChooser.currentColor()
        rgba = [ colour.red(), colour.green(), colour.blue(), colour.alpha() ]
        self.element_colours[text] = rgba
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
        debugger.print("on coloured button clicked")
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
        debugger.print("on_light_switches_cb_activated")
        self.light_switches[index] = not self.light_switches[index]
        string = f"switch light {index} off" if self.light_switches[index] else f"switch light {index} on"
        self.light_switches_cb.setItemText(index,string)
        self.opengl_widget.defineLights()
        self.calculate()
        self.plot()
        return

    def on_maximum_displacement_changed(self,value):
        """Handle changes to the maximum displacement setting.

        This method updates the 'Maximum displacement' setting, recalculates based on the new value, and then replots the relevant data or figures.

        Parameters
        ----------
        value : float
            The new value for the maximum displacement setting.

        Returns
        -------
        None

        """        
        debugger.print("on maximum_displacement changed ", value)
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
        debugger.print("on atom_scaling changed ", value)
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
        debugger.print("on arrow_radius changed ", value)
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
        debugger.print("on cell_radius changed ", value)
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
        debugger.print("on_bond_radius_changed")
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
        debugger.print("on_selected_mode_changed")
        if self.frequencies_cm1 is None or len(self.frequencies_cm1) < 1:
           debugger.print("on_selected_mode_changed aborting")
           return
        debugger.print("on selected_mode change mode was ", self.settings["Selected mode"])
        self.settings["Selected mode"] = self.selected_mode_sb.value()
        debugger.print("on selected_mode change mode is now ", self.settings["Selected mode"])
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
        debugger.print("on_plottype_cb_changed")
        self.plot_type_index = index
        debugger.print("Plot type index changed to ", self.plot_type_index)
        self.plot()
        return

    def calculate(self):
        """Perform calculations related to the notebook object, including processing program, file name, calculating frequencies, super cells, normal modes, bonds, center of mass, bounding box, element names, species, covalent radii, and updating the UI with calculated values.

        Parameters
        ----------
        None 

        Returns
        -------
        None

        """        
        debugger.print("Start:: calculate")
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings["Program"]
        filename = self.notebook.mainTab.getFullFileName()
        debugger.print("calculate program file name",program, filename)
        self.reader = self.notebook.reader
        if self.reader is None:
            debugger.print("Finished:: calculate - reader is None")
            return
        if program == "":
            debugger.print("Finished:: calculate - program is blank")
            return
        if filename == "":
            debugger.print("Finished:: calculate - filename is blank")
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Assemble the settingsTab settings
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        self.selected_mode_sb.setRange(1,len(self.frequencies_cm1))
        # Get the cell with whole molecules from the analysis tab
        self.notebook.analysisTab.refresh()
        self.unit_cell = self.notebook.analysisTab.cell_of_molecules
        if self.unit_cell is None:
            return
        # Generate a super cell
        imageSpecifier = self.settings["Super Cell"]
        self.super_cell = SuperCell(self.unit_cell,imageSpecifier)
        # if self.debug:
        #     self.unit_cell.printInfo()
        self.normal_modes = self.super_cell.calculateNormalModes(self.notebook.analysisTab.new_normal_modes)
        self.bonds = self.super_cell.calculateBonds()
        self.snapshot_number = 0
        self.nbonds = len(self.bonds)
        #
        self.XYZ = self.super_cell.calculateXYZ()
        self.natoms = len(self.XYZ)
        self.number_of_modes = len(self.normal_modes)
        # get the cell edges for the bounding box, shifted to the centre of mass origin
        totalMass,centreOfMassXYZ,centreOfMassABC = self.super_cell.calculateCentreOfMass(output=all)
        self.cell_corners,self.cell_edges = self.super_cell.getBoundingBox(centreOfMassABC)
        self.element_names = self.super_cell.getElementNames()
        self.species = self.reader.getSpecies()
        covalent_radii = self.notebook.analysisTab.element_radii
        self.radii = [self.settings["Atom scaling"]*covalent_radii[el] for el in self.element_names ]
        self.colours = [ self.element_colours[el] for el in self.element_names ]
        # reorder the displacement info in the normal modes into U,V and W lists
        # Using deque here rather than a simple list as the memory allocation doesn't have to be contiguous
        self.UVW.clear()
        for displacements in self.normal_modes:
            uvw = deque()
            for i in range(0,len(displacements),3):
                uvw.append( displacements[i:i+3] )
            self.UVW.append(uvw)
        # CalculatePhasePositions stores all the sphere and bond information
        self.calculatePhasePositions()
        # Add the arrows
        debugger.print("calculate: Selected mode",self.settings["Selected mode"])
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
        debugger.print("Finished:: calculate")
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
        debugger.print("setcolour")
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
        debugger.print("calculatePhasePositions")
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
        debugger.print("calculatePhasePositions - adding spheres and cylinders")
        for phase_index in range(len(phases)):
            for col, rad, xyz in zip(self.colours, self.radii, self.newXYZ[phase_index]):
                self.opengl_widget.addSphere(col, rad, xyz, phase=phase_index )
            for p in self.cell_corners:
                self.opengl_widget.addSphere(self.settings["Cell colour"], self.settings["Cell radius"], p, phase=phase_index )
            for bond in self.bonds:
                i,j = bond
                self.opengl_widget.addCylinder(self.settings["Bond colour"], self.settings["Bond radius"], self.newXYZ[phase_index,i], self.newXYZ[phase_index,j], phase=phase_index)
            for p1,p2 in self.cell_edges:
                self.opengl_widget.addCylinder(self.settings["Cell colour"], self.settings["Cell radius"], p1, p2, phase=phase_index)
        debugger.print("calculatePhasePositions - exiting")
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
        self.unit_cell.write_cif(filename=filename)

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
        This method relies on the `self.unit_cell` object, which should have a method `set_xyz_coordinates` and `write_cif`. It also uses `self.settings` to retrieve the number of phase steps and the selected mode for CIF generation. The function iterates over a calculated range of phases, updates the unit cell coordinates for each phase, and writes the CIF data to the provided filename. The function assumes that `self.newXYZ` is an iterable object containing new XYZ coordinates for each phase. The CIF files are saved with a description that includes the mode and phase.

        """        
        import copy
        unitcell = copy.deepcopy(self.unit_cell)
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
        debugger.print("Start:: plot")
        if self.reader is None:
            debugger.print("Finished:: plot reader is None")
            return
        if self.plot_type_index == 0:
            self.plot_animation()
        elif self.plot_type_index == 1:
            self.plot_arrows()
        else:
            self.plot_none()
        debugger.print("Finished:: plot")

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
        debugger.print("Start:: plot_animation")
        self.opengl_widget.showArrows(False)
        self.opengl_widget.stopAnimation()
        self.opengl_widget.update()
        debugger.print("Finished:: plot_animation")
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
        debugger.print("Start:: plot_animation")
        self.opengl_widget.showArrows(False)
        self.opengl_widget.update()
        self.opengl_widget.startAnimation()
        debugger.print("Finished:: plot_animation")
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
        debugger.print("Start:: plot_arrows")
        self.opengl_widget.showArrows(True)
        self.opengl_widget.stopAnimation()
        self.opengl_widget.update()
        debugger.print("Finished:: plot_arrows")
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
        debugger.print("Start:: requestRefresh")
        self.refreshRequired = True
        debugger.print("Finished:: requestRefresh")

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
        debugger.print("Start:: refresh")
        if not self.refreshRequired and not force:
            debugger.print("Finished:: refresh aborted",self.refreshRequired,force)
            return
        debugger.print("refresh widget",force)
        QApplication.setOverrideCursor(Qt.WaitCursor)
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
        debugger.print("refresh: selected mode is now",self.selected_mode_sb.value())
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
            debugger.print("reader is none ")
            self.element_names = []
            self.species = []
        else:
            self.notebook.analysisTab.refresh()
            self.unit_cell = self.notebook.analysisTab.cell_of_molecules
            if self.unit_cell is None:
                return
            self.element_names = self.unit_cell.element_names
            self.species = self.reader.getSpecies()
            debugger.print("refresh - species ",self.species)
            if self.settings["Element colours"] is not None:
                debugger.print("refresh - processing colours ",self.species)
                for el,colour in zip(self.species,self.settings["Element colours"]):
                    self.element_colours[el] = colour
            self.settings["Element colours"] = [ self.element_colours[el] for el in self.species ]
        count = self.element_coloured_hbox.count()
        if count == 0:
            debugger.print("initialising element colours widget")
            self.element_coloured_buttons = []
            for el in self.species:
                r,g,b,a = self.element_colours[el]
                button = QPushButton(el)
                button.setStyleSheet(f"background-color:rgba( {r}, {g}, {b}, {a});")
                button.clicked.connect(self.on_coloured_element_clicked)
                self.element_coloured_buttons.append(button)
                self.element_coloured_hbox.addWidget(button)
        else:
            debugger.print("update element colours widget")
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
        self.calculate()
        self.plot()
        self.refreshRequired = False
        QApplication.restoreOverrideCursor()
        debugger.print("Finished:: refresh")
        return
