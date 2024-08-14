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

from qtpy.QtWidgets import QFileDialog, QHBoxLayout, QLabel, QLineEdit, QPushButton, QWidget

from PDielec import __file__ as PDielec_init_filename
from PDielec.Materials import MaterialsDataBase
from PDielec.Utilities import Debug


class ScenarioTab(QWidget):
    """A class representing a tab for scenario configurations within a user interface.

    This class provides functionalities for managing and interacting with scenarios.
    It allows users to create, delete, and switch between different scenarios, as well as open and manipulate a materials database.
    :class:`~PDielec.GUI.PowderScenarioTab.PowderScenarioTab` and :class:`~PDielec.GUI.SingleCrystalScenarioTab.SingleCrystalScenarioTab` inherit from this class.

    Attributes
    ----------
    refreshRequired : bool
        Indicates whether a refresh is required.
    noCalculationsRequired : int
        The number of calculations that need to be performed.
    settings : dict
        A dictionary of settings for the scenario.
    notebook : QWidget
        The parent widget, which is expected to be the notebook container for the scenarios.
    scenarioType : type, optional
        The type of the scenario, e.g., powder, crystal.
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
    switchScenarioButton : QPushButton
        Button to switch between scenario types (e.g., powder to crystal).

    Methods
    -------
    openDataBase()
        Open the database and set the material names.

    getNoCalculationsRequired()
        Get the number of spectra that need recalculating from this scenario.

    requestRefresh()
        Request a refresh of the scenario

    set_reader(reader)
        Set the reader associated with this scenario

    setScenarioIndex(index)
        Set the index for the current scenario and a default legend name based on the index.

    print_settings()
        Print all the settings for this scenario.

    on_legend_le_changed(text)
        Handle a change in the scenario legend.

    add_scenario_buttons()
        Add a set of scenario buttons in an hbox. Return the hbox.

    addScenarioButtonClicked()
        Handle when the "add another scenario" button has been clicked.

    deleteScenarioButtonClicked()
        Handle when a delete button has been clicked.

    switchScenarioButtonClicked()
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
        global debugger
        debugger = Debug(debug,"ScenarioTab:")
        debugger.print("Start:: initialiser")
        self.refreshRequired = True
        self.noCalculationsRequired = 0
        self.settings = {}
        self.notebook = parent
        self.settings["Legend"] = "Unset"
        self.scenarioType = None
        self.settings["Scenario type"] = "Unset"
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
        debugger.print("Finished:: initialiser")

    def openDataBase(self):
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
        self.DataBase = MaterialsDataBase(filename,debug=debugger.state())
        sheets = self.DataBase.getSheetNames()
        if not self.DataBase.valid():
            self.DataBase = oldDataBase
            print("Error chosen file is not a materials database",sheets)
            return
        self.settings["Materials database"] = self.DataBase.getFileName()
        self.database_le.setText(self.settings["Materials database"])
        self.materialNames = sheets
        return

    def getNoCalculationsRequired(self):
        """Get the number of spectra that need recalculating from this scenario.

        Parameters
        ----------
        None

        Returns
        -------
        int
            The number of spectra requiring recalculation.

        """
        result = self.noCalculationsRequired if self.refreshRequired else 0
        debugger.print(self.settings["Legend"], "getNoCalculationsRequired",result)
        return result

    def requestRefresh(self):
        """Request a refresh of the scenario.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        debugger.print(self.settings["Legend"], "requestRefresh")
        self.refreshRequired = True
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
        debugger.print(self.settings["Legend"], "set_reader")
        self.refreshRequired = True
        self.reader = reader
        return

    def setScenarioIndex(self,index):
        """Set the index for the current scenario and a default legend name based on the index.

        Parameters
        ----------
        index : int
            The scenario index.  

        Returns
        -------
        None

        """
        debugger.print(self.settings["Legend"], "setScenarioIndex",index)
        self.scenarioIndex = index
        text = self.settings["Legend"]
        if text.startswith("Unset") or text.startswith("Scenario ") or text.startswith("Powder scenario ") or text.startswith("Crystal scenario "):
            debugger.print(self.settings["Legend"], "setScenarioIndex changing scenario legend","Scenario "+str(index+1))
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
        debugger.print(self.settings["Legend"], "print_settings")
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
        debugger.print(self.settings["Legend"], "on_legend_le_changed",text)
        # self.refreshRequired = True
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
        debugger.print(self.settings["Legend"], "add_scenario_buttons start")
        hbox = QHBoxLayout()
        self.addScenarioButton = QPushButton("Add another scenario")
        self.addScenarioButton.setToolTip("Add another scenario the the notebook tabs, \nthe new scenario is added to the end of the current tab list")
        self.addScenarioButton.clicked.connect(self.addScenarioButtonClicked)
        hbox.addWidget(self.addScenarioButton)
        self.deleteScenarioButton = QPushButton("Delete this scenario")
        self.deleteScenarioButton.setToolTip("Delete the current scenario")
        self.deleteScenarioButton.clicked.connect(self.deleteScenarioButtonClicked)
        hbox.addWidget(self.deleteScenarioButton)
        if self.scenarioType == "Powder":
            self.switchScenarioButton = QPushButton("Switch to crystal scenario")
            self.switchScenarioButton.setToolTip("Switch the current scenario to a single crystal scenario")
        else:
            self.switchScenarioButton = QPushButton("Switch to powder scenario")
            self.switchScenarioButton.setToolTip("Switch the current scenario to a powder scenario")
        self.switchScenarioButton.clicked.connect(self.switchScenarioButtonClicked)
        hbox.addWidget(self.switchScenarioButton)
        debugger.print(self.settings["Legend"], "add_scenario_buttons finish")
        return hbox

    def addScenarioButtonClicked(self):
        """Handle when the "add another scenario" button has been clicked.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # Add another scenario
        debugger.print(self.settings["Legend"], "addScenarioButtonClicked")
        self.notebook.addScenario(copyFromIndex=self.scenarioIndex)
        return

    def deleteScenarioButtonClicked(self):
        """Handle when a delete button has been clicked.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        debugger.print(self.settings["Legend"], "deleteScenarioButtonClicked")
        self.notebook.deleteScenario(self.scenarioIndex)
        return

    def switchScenarioButtonClicked(self):
        """Handle the switch scenario when a button has been clicked.

        Asks the notebook to change the scenario type (Powder or Crystal) to the other.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        debugger.print(self.settings["Legend"],"switchScenarioButtonClicked")
        self.notebook.switchScenario(self.scenarioIndex)
        self.requestRefresh()
        self.notebook.plottingTab.refresh(force=True)
        return




