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
"""ScenarioTab module."""
# -*- coding: utf8 -*-
import os

from qtpy.QtCore import Qt
from qtpy.QtWidgets import QFileDialog, QHBoxLayout, QLabel, QLineEdit, QPushButton, QWidget
from qtpy.QtWidgets import QComboBox

from PDielec import __file__ as PDielec_init_filename
from PDielec.Materials import MaterialsDataBase
from PDielec.Utilities import Debug


class ScenarioTab(QWidget):
    """A class representing a tab for scenario configurations within a user interface.

    This class provides functionalities for managing and interacting with scenarios.
    It allows users to create, delete, and switch between different scenarios, as well as open and manipulate a materials database.
    :class:`~PDielec.GUI.PowderInfraredScenarioTab`, :class:`~PDielec.GUI.CrystalInfraredScenarioTab`,
    :class:`~PDielec.GUI.PowderRamanScenarioTab` and :class:`~PDielec.GUI.CrystalRamanScenarioTab` inherit from this class.

    Attributes
    ----------
    refresh_required : bool
        Indicates whether a refresh is required.
    noCalculationsRequired : int
        The number of calculations that need to be performed.
    settings : dict
        A dictionary of settings for the scenario.
    notebook : QWidget
        The parent widget, which is expected to be the notebook container for the scenarios.
    scenarioType : type, optional
        The type of the scenario, e.g., Powder Infrared, Crystal Infrared, Powder Raman, Crystal Raman.
    scenarioTypes : list of scenario types
        The list is obtained from the keys of the self.notebook.scenarioTypes dictionary
    vs_cm1 : list
        List containing default values for some settings.
    DataBase : MaterialsDataBase
        An instance of a class for interacting with a materials database.
    openDB_button : QPushButton
        Button to open the materials database.
    openDB_label : QLabel
        Label associated with the openDB_button.
    database_le : QLineEdit
        Line edit showing the path/name of the current materials database.
    database_le_label : QLabel
        Label associated with the database_le.
    legend_le : QLineEdit
        Line edit for specifying a custom legend for the scenario.
    addScenarioButton : QPushButton
        Button to add another scenario.
    deleteScenarioButton : QPushButton
        Button to delete the current scenario.
    switchScenarioCB : QComboBOx
        Dropdown menu to switch between scenario types (e.g., Powder Infrared, Cystal Infrared, Powder Raman or Crystal Raman).

    Methods
    -------
    open_data_base()
        Open the database and set the material names.

    get_no_calculations_required()
        Get the number of spectra that need recalculating from this scenario.

    request_refresh()
        Request a refresh of the scenario

    set_reader(reader)
        Set the reader associated with this scenario

    set_scenario_index(index)
        Set the index for the current scenario and a default legend name based on the index.

    print_settings()
        Print all the settings for this scenario.

    on_legend_le_changed(text)
        Handle a change in the scenario legend.

    add_scenario_buttons()
        Add a set of scenario buttons in an hbox. Return the hbox.

    add_scenario_button_clicked()
        Handle when the "add another scenario" button has been clicked.

    delete_scenario_button_clicked()
        Handle when a delete button has been clicked.

    switchScenarioCBClicked()
        Handle the scenario when a button has been clicked.

    """

    def __init__(self, parent, debug=False):
        """Initialise ScenarioTabs.

        Parameters
        ----------
        parent : QWidget
            The parent widget.
        debug : bool, optional
            Flag to turn debugging messages on or off.

        """
        super(QWidget, self).__init__(parent)
        self.debugger = Debug(debug,"ScenarioTab:")
        self.debugger.print("Start:: initialiser")
        self.refresh_required = True
        self.noCalculationsRequired = 0
        self.settings = {}
        self.notebook = parent
        self.settings["Legend"] = "Unset"
        self.scenarioType = None
        self.settings["Scenario type"] = "Unset"
        self.scenarioTypes = list(self.notebook.scenarioTypes.keys())
        self.vs_cm1 = [0, 0]
        # Deal with the Materials Database here as it is used in all Scenarios
        PDielec_Directory = os.path.dirname(PDielec_init_filename)
        filename  = os.path.join(PDielec_Directory, "MaterialsDataBase.xlsx")
        filename  = os.path.relpath(filename)
        self.settings["Materials database"] = filename
        # Open the database
        self.DataBase = MaterialsDataBase(self.settings["Materials database"],debug=debug)
        # Set up the open database button
        self.openDB_button = QPushButton("Open materials' database")
        self.openDB_button.clicked.connect(self.openDB_button_clicked)
        self.openDB_button.setToolTip("Open a new materials' database (.xlsx file)")
        self.openDB_label = QLabel("Open materials' database")
        self.openDB_label.setToolTip("Open a new materials' database (.xlsx file)")
        # set up the database information line
        self.database_le = QLineEdit(self)
        self.database_le.setToolTip("Provides information about the name of the materials' database")
        self.database_le.setText(self.settings["Materials database"])
        self.database_le.setReadOnly(True)
        self.database_le_label = QLabel("Current materials' database")
        self.database_le_label.setToolTip("Provides information about the name of the materials' database")
        self.debugger.print("Finished:: initialiser")

    def open_data_base(self):
        """Open the database and set the material names.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        selfilter = "Spreadsheet (*.xlsx)"
        filename,myfilter = QFileDialog.getOpenFileName(self,"Open spreadsheet","","Spreadsheet (*.xls);;Spreadsheet (*.xlsx);;All files(*)",selfilter)
        # Process the filename
        if filename == "":
            return
        oldDataBase = self.DataBase
        self.DataBase = MaterialsDataBase(filename,debug=self.debugger.state())
        sheets = self.DataBase.get_sheet_names()
        if not self.DataBase.valid():
            self.DataBase = oldDataBase
            print("Error chosen file is not a materials database",sheets)
            return
        self.settings["Materials database"] = self.DataBase.get_file_name()
        self.database_le.setText(self.settings["Materials database"])
        self.materialNames = sheets
        return

    def get_no_calculations_required(self):
        """Get the number of spectra that need recalculating from this scenario.

        Parameters
        ----------
        None

        Returns
        -------
        int
            The number of spectra requiring recalculation.

        """
        result = self.noCalculationsRequired if self.refresh_required else 0
        self.debugger.print(self.settings["Legend"], "get_no_calculations_required",result)
        return result

    def request_refresh(self):
        """Request a refresh of the scenario.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.debugger.print(self.settings["Legend"], "request_refresh")
        self.refresh_required = True
        return

    def set_reader(self,reader):
        """Set the reader associated with this scenario.

        Parameters
        ----------
        reader : a reader object
            The reader to be used for this scenario

        Returns
        -------
        None

        """
        self.debugger.print(self.settings["Legend"], "set_reader")
        self.refresh_required = True
        self.reader = reader
        return

    def set_scenario_index(self,index):
        """Set the index for the current scenario and a default legend name based on the index.

        Parameters
        ----------
        index : int
            The scenario index.  

        Returns
        -------
        None

        """
        self.debugger.print(self.settings["Legend"], "set_scenario_index",index)
        self.scenarioIndex = index
        text = self.settings["Legend"]
        if text.startswith("Unset") or text.startswith("Scenario ") or text.startswith("Powder scenario ") or text.startswith("Crystal scenario "):
            self.debugger.print(self.settings["Legend"], "set_scenario_index changing scenario legend","Scenario "+str(index+1))
            self.legend_le.setText("Scenario "+str(index + 1))
            self.settings["Legend"] = "Scenario "+str(index + 1)
        return

    def print_settings(self):
        """Print all the settings for this scenario.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.debugger.print(self.settings["Legend"], "print_settings")
        print("#")
        print("# Scenario tab")
        print("#")
        print("tab = self.notebook.scenarios")
        for key in self.settings:
            print(key, self.settings[key])

    def on_legend_le_changed(self,text):
        """Handle a change in the scenario legend.

        Parameters
        ----------
        text : str
            The new legend for this scenario

        Returns
        -------
        None

        """
        self.debugger.print(self.settings["Legend"], "on_legend_le_changed",text)
        # self.refresh_required = True
        self.settings["Legend"] = text
        return

    def  add_scenario_buttons(self):
        """Add a set of scenario buttons in an hbox. Return the hbox.

        Parameters
        ----------
        None

        Returns
        -------
        hbox
            The hbox containing the scenario buttons.

        """
        self.debugger.print(self.settings["Legend"], "add_scenario_buttons start")
        hbox = QHBoxLayout()
        self.addScenarioButton = QPushButton("Add another scenario")
        self.addScenarioButton.setToolTip("Add another scenario the the notebook tabs, \nthe new scenario is added to the end of the current tab list")
        self.addScenarioButton.clicked.connect(self.add_scenario_button_clicked)
        hbox.addWidget(self.addScenarioButton)
        self.deleteScenarioButton = QPushButton("Delete this scenario")
        self.deleteScenarioButton.setToolTip("Delete the current scenario")
        self.deleteScenarioButton.clicked.connect(self.delete_scenario_button_clicked)
        hbox.addWidget(self.deleteScenarioButton)
        self.switchScenarioCB = QComboBox(self)
        self.switchScenarioCB.setToolTip("Switch to a new scenario type")
        self.switchScenarioCB.addItems(self.scenarioTypes)
        index = self.switchScenarioCB.findText(self.scenarioType, Qt.MatchFixedString)
        if index >=0:
            self.switchScenarioCB.setCurrentIndex(index)
        else:
            print("Error in scenarioType", self.scenarioType, self.scenarioTypes)
        self.switchScenarioCB.activated.connect(self.switch_scenario_cb_activated)
        label = QLabel("Switch to a new scenario type:")
        label.setAlignment(Qt.AlignBottom | Qt.AlignRight)
        hbox.addWidget(label)
        hbox.addWidget(self.switchScenarioCB)
        self.debugger.print(self.settings["Legend"], "add_scenario_buttons finish")
        return hbox

    def add_scenario_button_clicked(self):
        """Handle when the "add another scenario" button has been clicked.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # Add another scenario
        self.debugger.print(self.settings["Legend"], "add_scenario_button_clicked")
        self.notebook.add_scenario(copyFromIndex=self.scenarioIndex)
        return

    def delete_scenario_button_clicked(self):
        """Handle when a delete button has been clicked.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.debugger.print(self.settings["Legend"], "delete_scenario_button_clicked")
        self.notebook.delete_scenario(self.scenarioIndex)
        return

    def switch_scenario_cb_activated(self, index):
        """Handle the switch scenario when combo box has been activated.

        Asks the notebook to change the scenario type 

        Parameters
        ----------
        index
            The index in the list of scenarios

        Returns
        -------
        None

        """
        self.debugger.print(self.settings["Legend"],"switch_scenario_cb_activated", index)
        self.scenarioType = self.scenarioTypes[index]
        self.notebook.switch_scenario(self.scenarioIndex, self.scenarioType)
        self.request_refresh()
        self.notebook.plottingTab.refresh(force=True)
        return




