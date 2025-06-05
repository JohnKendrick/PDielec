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
"""NoteBook module."""
import copy
import os
import sys

import psutil
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QApplication, QFileDialog, QMessageBox, QTabWidget, QVBoxLayout, QWidget

from PDielec import Calculator
from PDielec.GUI.AnalysisTab import AnalysisTab
from PDielec.GUI.FitterTab import FitterTab
from PDielec.GUI.MainTab import MainTab
from PDielec.GUI.PlottingTab import PlottingTab
from PDielec.GUI.PowderScenarioTab import PowderScenarioTab
from PDielec.GUI.SettingsTab import SettingsTab
from PDielec.GUI.SingleCrystalScenarioTab import SingleCrystalScenarioTab
from PDielec.GUI.SpreadSheetManager import SpreadSheetManager
from PDielec.GUI.ViewerTab import ViewerTab
from PDielec.Utilities import Debug


class NoteBook(QWidget):
    """A Qt widget that holds tabs for managing different aspects of a notebook-like interface.

    This class manages a complex interface which includes a range of functionalities 
    such as handling various scenarios (e.g., Powder or SingleCrystal), managing computational 
    resources, offering plotting abilities, managing settings, etc. It utilizes multithreading, 
    handles external spreadsheet files for input and output, and provides a user interface 
    for configuring and monitoring the progress of analytical or simulation tasks.

    Parameters
    ----------
    parent : QWidget
        The parent widget.
    program : str
        Description or identifier of the program associated with the notebook.
    filename : str
        Name of the file, if any, associated with the notebook.
    spreadsheet : type
        The spreadsheet object associated with the notebook.
    debug : bool, optional
        Flag to activate debug mode, by default False.
    progressbar : QProgressBar, optional
        A progressbar object to reflect the current progress, by default None.
    scripting : bool, optional
        Flag to indicate if the notebook is used in scripting mode, by default False.
    default_scenario : str, optional
        The type of default scenario to load at initiation, by default 'powder'.
    ncpus : int, optional
        The number of CPUs to use, by default 0 which means autodetect.
    threading : bool, optional
        Flag to enable threading, by default False.

    Attributes
    ----------
    app : QWidget
        The application the notebook is part of.
    reader : type
        Reserved for future use, currently None.
    progressbars : list of QProgressBar
        List of progress bar objects for tracking progress.
    progressbar_status : int
        Current status of the progressbar.
    progressbar_maximum : int
        Maximum value of the progressbar.
    spreadsheet : type
        Spreadsheet object associated with the notebook.
    threading : bool
        State of threading use.
    ncpus : int
        Number of CPUs being used.
    layout : QVBoxLayout
        The layout for adding widgets to the notebook.
    tabOffSet : int
        Offset value for tab indices.
    plottingTab : QWidget
        Tab for managing plotting related functions.
    settingsTab : QWidget
        Tab for managing settings.
    analysisTab : QWidget
        Tab for managing analysis functions.
    viewerTab : QWidget
        Tab for managing 3D viewing options.
    fitterTab : QWidget
        Tab for managing fitting tasks.
    scenarios : list of QWidget
        List containing scenario-specific tabs.
    tabs : QTabWidget
        The QTabWidget managing multiple tabs in the notebook interface.
    debug : bool
        Debug mode state.
    overwriting : bool
        State indicating if overwriting files without prompt is enabled.

    """

    def __init__(self, parent, program, filename, spreadsheet, debug=False, progressbar=None, scripting=False, default_scenario="powder",ncpus=0, threading=False):
        """Initialise the main NoteBook.

        This method initializes the main widget with all necessary components
        such as the settings, plotting, analysis, viewer, and fitter tabs, along
        with the scenarios management logic. It sets up the multiprocessing capabilities,
        reads the provided files, and sets up the debugging interface.

        Parameters
        ----------
        parent : QWidget
            The parent widget which this widget is a part of.
        program : str
            Represents the program used to generate the output
        filename : str
            The path to a file that is relevant for initializing the application.
            This could be a configuration file, a data file, etc.
        spreadsheet : Various possible types
            An object or identifier for a spreadsheet component or functionality
            within the application. The exact type and use can vary based on context.
        debug : bool, optional
            Flag to enable debugging mode. Defaults to False.
        progressbar : QProgressBar or similar, optional
            A progress bar component for the UI. If None, no progress bar is used.
            Defaults to None.
        scripting : bool, optional
            Flag to indicate whether the application is being used in a scripting mode.
            Defaults to False.
        default_scenario : str, optional
            Specifies the default scenario to be loaded at startup. Possible values
            might include 'powder', 'single crystal', etc. Defaults to 'powder'.
        ncpus : int, optional
            The number of CPUs to be used for multiprocessing. If set to 0, the
            application will try to use all physical cores available. Defaults to 0.
        threading : bool, optional
            Flag to enable or disable threading features. Defaults to False.

        Notes
        -----
        This constructor initializes the GUI components of the application,
        sets up multiprocessing/threading as specified, handles file reading
        and application settings, and prepares the different analysis and
        visualization tabs. It is vital for setting up the initial state
        of the application's main window or a major component within a larger
        UI framework.

        """        
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,"NoteBook:")
        debugger.print("Start:: Initialising")
        self.app = parent
        self.reader = None
        self.progressbars=[progressbar]
        if progressbar is None:
            self.progressbars = [ ]
        self.progressbar_status = 0
        self.progressbar_maximum = 0
        self.spreadsheet = None
        self.threading = threading
        if default_scenario == "powder":
            self.currentScenarioTab = PowderScenarioTab
        else:
            self.currentScenarioTab = SingleCrystalScenarioTab
        if ncpus == 0:
            self.ncpus = psutil.cpu_count(logical=False)
        else:
            self.ncpus = ncpus
        self.startPool()
        self.scripting = scripting
        # Overwriting of files is not allowed with a prompt
        # If scripting is used then overwriting is allowed
        self.overwriting = False
        self.debug = debug
        self.layout = QVBoxLayout()
        # The number of tabs before we have scenarios
        self.tabOffSet = 2
        # Set the plotting tab to None in case a scenario tries to read it
        self.plottingTab = None
        self.settingsTab = None
        self.analysisTab = None
        self.viewerTab = None
        self.fitterTab = None
        self.scenarios = None
        #
        # Initialize tab screen
        #
        self.tabs = QTabWidget(self)
        self.tabs.currentChanged.connect(self.on_tabs_currentChanged)
        self.mainTab = MainTab(self, program, filename, spreadsheet, debug=debug)
        self.settingsTab = SettingsTab(self, debug=debug)
        if filename != "" and not self.scripting:
            debugger.print("Refreshing settingsTab in notebook initialisation - filename",filename)
            self.settingsTab.refresh()
        #
        # Open more windows
        #
        debugger.print("Initialising the first scenario")
        self.scenarios = []
        self.scenarios.append( self.currentScenarioTab(self, debug=debug ) )
        self.scenarios[0].setScenarioIndex(0)
        self.scenarios[0].settings["Legend"] = "Scenario 1"
        debugger.print("Finished adding the first scenario")
        #
        # Open the plotting tab
        #
        self.plottingTab = PlottingTab(self, debug=debug)
        if filename != "" and not self.scripting:
            debugger.print("Refreshing plotting because filename is set")
            self.plottingTab.refresh()
        #
        # Open the Analysis tab
        #
        self.analysisTab = AnalysisTab(self, debug=debug)
        if filename != "" and not self.scripting:
            debugger.print("Refreshing analysis because filename is set")
            self.analysisTab.refresh()
        #
        # Open the Viewer tab
        #
        debugger.print("Initialising the viewer tab")
        self.viewerTab = ViewerTab(self, debug=debug)
        #
        # Open the Fitter tab
        #
        debugger.print("Initialising the fitter tab")
        self.fitterTab = FitterTab(self, debug=debug)
        #
        # Add tabs
        #
        debugger.print("Adding all tabs to the notebook")
        self.tabs.addTab(self.mainTab,"Main")
        self.tabs.addTab(self.settingsTab,"Settings")
        for i,tab in enumerate(self.scenarios):
            tab.requestRefresh()
            self.tabs.addTab(tab,"Scenario "+str(i+1))
        self.tabs.addTab(self.plottingTab,"Plotting")
        self.tabs.addTab(self.analysisTab,"Analysis")
        self.tabs.addTab(self.viewerTab,"3D Viewer")
        self.tabs.addTab(self.fitterTab,"Fitter")

        # Add the tab widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        debugger.print("Finished:: Initialising")
        return

    def startPool(self):
        """Initialize a pool of worker processes or threads for computation.

        This method initializes a pool based on the instance's specified number of CPUs and the threading model. It accesses a global debugger variable for potential debugging purposes.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - The `Calculator.get_pool` method is expected to be a class method or static method of the class `Calculator` that initializes a pool of worker processes or threads.
        - The pool is stored in the instance's `pool` attribute.
        - The method uses a global variable `debugger` which should be defined elsewhere in the global scope for debugging purposes.
        - The number of CPUs (`ncpus`) and the threading model (`threading`) are not parameters of this method, but are expected to be attributes of the instance (`self`).

        """        
        global debugger
        self.pool = Calculator.get_pool(self.ncpus,self.threading, debugger = debugger)
        return

    def requestRefresh(self):
        """Request a refresh operation.

        This method toggles the refresh requirement state to true, indicating that a refresh is needed.

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

    def addScenario(self,scenarioType=None,copyFromIndex=-2):
        """Add a new scenario tab.

        If a copy is requested then copy the appropriate information into the new tab.

        Parameters
        ----------
        scenarioType : string
             scenarioType can be one of 'Powder' or 'SingleCrystal'
        copyFromIndex : int
             if copyFromIndex is not -2 then use the index to determine the scenario type
             if copyFromIndex is -2 and the scenarioType has not been specified then just use the last scenario type there is
             Otherwise just find the last scenario type to copy

        Returns
        -------
        None

        """
        debugger.print("Start:: addScenario for scenarioType", scenarioType,copyFromIndex)
        if copyFromIndex != -2:
            # If the copyFromIndex is not -2 then we override the scenarioType
            last = self.scenarios[copyFromIndex]
            scenarioType = last.scenarioType
            debugger.print("scenario type has been set from copyFromIndex",scenarioType)
        elif scenarioType is None:
            # The default behaviour with no parameters in the call, use the last scenario in the list
            last = self.scenarios[-1]
            scenarioType = last.scenarioType
            debugger.print("scenario type has been set from the last scenario",scenarioType)
        else:
            # copyFromIndex is default so we find the last scenario of scenarioType in the list
            last = None
            for scenario in self.scenarios:
                if scenarioType == scenario.scenarioType:
                    last = scenario
            # end for
        # Create a new scenario
        if scenarioType == "Powder":
            self.currentScenarioTab = PowderScenarioTab
        else:
            self.currentScenarioTab = SingleCrystalScenarioTab
        # Add the scenario to the end of the list
        debugger.print("Appending the new scenario")
        self.scenarios.append(self.currentScenarioTab(self, self.debug))
        # If we have found a previous scenario of the same time set the settings to it
        debugger.print("Checking the value of last",last)
        if last is not None:
            debugger.print("Copying settings from old to new scenario")
            self.scenarios[-1].settings = copy.deepcopy(last.settings)
        self.scenarios[-1].requestRefresh()
        self.scenarios[-1].refresh()
        n = len(self.scenarios)
        self.tabs.insertTab(self.tabOffSet+n-1,self.scenarios[-1],"Scenario "+str(n))
        self.tabs.setCurrentIndex(self.tabOffSet+n-1)
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
            self.tabs.setTabText(self.tabOffSet+i,"Scenario "+str(i+1))
        debugger.print("Finished:: addScenario for scenarioType", scenarioType,copyFromIndex)
        return

    def print_settings(self, filename=None):
        # Print the settings of all the settings that have been used to a file settings.py
        """Print the current program settings to a file.

        Each tab in the notebook is processed in turn, including all the scenarios.
        This allows a script to be written which can be used to recall a configuration of the program.
        Some settings need careful handling, for example the sigmas_cm1 values for the settingsTab

        Parameters
        ----------
        filename : str, optional
            The name of the file to save the settings. If not specified, a save file dialog will be shown.

        Returns
        -------
        None

        """        
        debugger.print("Start:: print_settings, filename=",filename)
        qf = QFileDialog()
        qf.setWindowTitle("Save the program settings to a file")
        debugger.print("print_settings, directory=",self.mainTab.directory)
        qf.setDirectory(self.mainTab.directory)
        if filename is None:
            filename,selection = qf.getSaveFileName()
        if filename == "":
            debugger.print("Start:: print_settings, filename is blank")
            return
        print("Current settings will be saved to "+filename)
        with open(filename,"w") as fd:
            # Handle the special case of the first scenario
            print("#",file=fd)
            print("# Handle the special case of the first scenario",file=fd)
            print("#",file=fd)
            print('self.notebook.switchScenario(0,scenarioType="'+self.scenarios[0].scenarioType+'")',file=fd )
            print("#",file=fd)
            # Print settings of mainTab
            self.print_tab_settings(self.mainTab, "mainTab",fd)
            # Print settings of settingsTab
            self.print_tab_settings(self.settingsTab, "settingsTab",fd)
            print("tab.sigmas_cm1 =",self.settingsTab.sigmas_cm1,file=fd)
            # Print settings of all scenarios
            for i,tab in enumerate(self.scenarios):
                if i == 0:
                    self.print_tab_settings(tab, f"scenarios[{i}]", fd, new_scenario = False)
                else:
                    self.print_tab_settings(tab, f"scenarios[{i}]", fd, new_scenario = True)
            # print analysis tab settings
            self.print_tab_settings(self.analysisTab, "analysisTab",fd)
            # print viewer tab settings
            self.print_tab_settings(self.viewerTab, "viewerTab",fd)
            # print fitter tab settings
            self.print_tab_settings(self.fitterTab, "fitterTab",fd)
            # print plotting tab settings
            self.print_tab_settings(self.plottingTab, "plottingTab",fd)
        debugger.print("Finished:: print_settings, filename=",filename)
        return

    def print_tab_settings(self,tab,title,fd,new_scenario = False):
        """Print the configuration settings of a specified tab to a file descriptor.

        Parameters
        ----------
        tab : Object
            An object representing a tab that contains settings and other attributes.
        title : str
            The title of the tab.
        fd : file descriptor
            An open file descriptor where the settings should be printed.
        new_scenario : bool, optional
            A flag indicating whether this is a new scenario, by default False.

        Returns
        -------
        None

        Notes
        -----
        This function iterates through the settings of the provided tab object and
        prints each setting to the provided file descriptor. If `new_scenario` is
        True, it adds a scenario definition to the file. It handles special cases
        for certain settings like 'Optical permittivity' and 'Mass definition',
        and formats string values and others appropriately for printing.

        """        
        debugger.print("Start:: print_tab_settings")
        print("#",file=fd)
        print("#",file=fd)
        if new_scenario:
            print('self.notebook.addScenario(scenarioType="'+tab.scenarioType+'")',file=fd )
        print("tab = self.notebook."+title,file=fd)
        for item in tab.settings:
            if item == "Optical permittivity" and not tab.settings["Optical permittivity edited"]:
                    pass
            elif item == "Mass definition":
                print("tab.settings['"+item+f"'] = '{tab.settings[item]}'",file=fd)
                # Check to see if the mass_definition is gui, if so set all the masses
                if tab.settings[item] == "gui":
                    for c in tab.masses_dictionary:
                        print("tab.masses_dictionary['"+c+"'] = ",tab.masses_dictionary[c],file=fd)
            else:
                value = tab.settings[item]
                if "str" in str(type(value)):
                    print("tab.settings['"+item+f"'] = '{tab.settings[item]}'",file=fd)
                else:
                    print("tab.settings['"+item+"'] = ", tab.settings[item],file=fd)
        debugger.print("Finished:: print_tab_settings")

    def deleteAllScenarios(self):
        """Delete all scenarios except the first one.

        This method sequentially deletes each scenario from the end of the collection until only one scenario is left. It also removes the corresponding tabs.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("Start:: deleteAllScenarios")
        # Don't delete the last scenario
        index = len(self.scenarios)-1
        while len(self.scenarios) > 1:
            self.tabs.removeTab(self.tabOffSet+index)
            del self.scenarios[index]
            index -= 1
        debugger.print("Finished:: deleteAllScenarios")
        return

    def deleteScenario(self,index):
        """Delete a scenario from the scenarios list and update tabs accordingly.

        Parameters
        ----------
        index : int
            The index of the scenario to be deleted.

        Returns
        -------
        None

        Notes
        -----
        This function will remove a scenario from the scenarios list at the specified index
        and update the tabs in the UI to reflect this change. It ensures that there is at least
        one scenario remaining. If the scenario to be deleted is the first one, the selection
        moves to the next available scenario. Otherwise, it selects the previous scenario.

        """        
        debugger.print("Start:: deleteScenario",index)
        # Don't delete the last scenario
        if len(self.scenarios) > 1:
            self.tabs.removeTab(self.tabOffSet+index)
            del self.scenarios[index]
            for i,scenario in enumerate(self.scenarios):
                scenario.setScenarioIndex(i)
                self.tabs.setTabText(self.tabOffSet+i,"Scenario "+str(i+1))
            if index-1 < 0:
                index += 1
            self.tabs.setCurrentIndex(self.tabOffSet+index-1)
        debugger.print("Finished:: deleteScenario",index)
        return

    def switchScenario(self,index,scenarioType=None):
        """Switch the scenario tab based on the scenario type.

        Parameters
        ----------
        index : int
            The index of the scenario to switch to.
        scenarioType : str, optional
            The type of scenario to switch to. Can be 'Powder' or None. If None, the
            scenario type is determined by the current scenario's type.
            If the scenarioType is anything else a 'SingleCrystal' scenario is switched to.

        Returns
        -------
        None

        Notes
        -----
        - This method switches the current scenario to a new scenario based on the
          given scenario type or the current scenario's type if no scenario type is
          provided.
        - The method updates the current scenario tab and refreshes the UI to reflect
          the new scenario.
        - If `scenarioType` is 'Powder', switch to a Powder scenario tab. Otherwise,
          switch to a Single Crystal scenario tab.
        - This function also updates all scenario tabs' names based on their index
          and requests a refresh of the currently selected scenario.

        See Also
        --------
        SingleCrystalScenarioTab, PowderScenarioTab : Classes representing different types of scenario tabs.

        """        
        debugger.print("Start:: switch for scenario", index+1)
        # Replace the scenario with the other scenario type
        scenario = self.scenarios[index]
        debugger.print("Current scenario type", scenario.scenarioType, scenarioType)
        #
        # If scenarioType is specified in the call then force that type
        # Otherwise switch type
        #
        if scenarioType is None:
            if scenario.scenarioType == "Powder":
                self.currentScenarioTab = SingleCrystalScenarioTab
            else:
                self.currentScenarioTab = PowderScenarioTab
            # end if
        elif scenarioType == "Powder":
            self.currentScenarioTab = PowderScenarioTab
        else:
            self.currentScenarioTab = SingleCrystalScenarioTab
            # end if
        #end if
        self.scenarios[index] =  self.currentScenarioTab(self, self.debug)
        scenario = self.scenarios[index]
        debugger.print("Current scenario type now", scenario.scenarioType)
        self.tabs.removeTab(self.tabOffSet+index)
        self.tabs.insertTab(self.tabOffSet+index,scenario,"Scenario "+str(index+1) )
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
            self.tabs.setTabText(self.tabOffSet+i,"Scenario "+str(i+1))
        self.scenarios[index].requestRefresh()
        if not self.scripting:
            self.scenarios[index].refresh()
        self.tabs.setCurrentIndex(self.tabOffSet+index)
        debugger.print("Finished:: switch for scenario", index+1)
        return

    def refresh(self,force=False):
        """Refresh the current state, optionally forcing a refresh regardless of scripting constraints.

        Parameters
        ----------
        force : bool, optional
            Force a refresh even if scripting is active, by default False.

        Returns
        -------
        None

        Notes
        -----
        This method initiates a refresh process on various components such as the main settings, scenarios, and several tabs including plotting, analysis, viewer, and fitter. It adjusts the active tab based on the current number of scenarios. If 'force' is set to True, the refresh process is executed disregarding any active scripting conditions.

        """        
        debugger.print("Started:: newrefresh",force)
        if not force and self.scripting:
            debugger.print("Finished:: newrefresh Notebook aborting refresh because of scripting")
            return
        ntabs = 2 + len(self.scenarios) + 4
        # Do a refresh on the mainTab and the settingsTab
        # Add all the scenarios
        # This should caused anything that needs reading in to be read and processed
        self.mainTab.refresh(force=force)
        self.settingsTab.refresh(force=force)
        for tab in self.scenarios:
            tab.refresh(force=force)
        # Request refreshes on everything else
        self.plottingTab.requestRefresh()
        self.analysisTab.requestRefresh()
        self.viewerTab.requestRefresh()
        self.fitterTab.requestRefresh()
        # In a script we do not change the tab index, but we need the analysis tab and the plotter tab to be refreshed
        # So do it here, leave the GUI after a script showing the plotter tab
        self.tabs.setCurrentIndex(ntabs-3)
        self.tabs.setCurrentIndex(ntabs-4)
        debugger.print("Finished:: newrefresh",force)

    def writeSpreadsheet(self):
        """Write data to an Excel spreadsheet.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This function assumes that the Excel spreadsheet is an attribute of the object this method belongs to. It attempts to write data to various tabs within the spreadsheet, namely 'mainTab', 'settingsTab', 'analysisTab', and 'plottingTab'. The method opens the spreadsheet, writes data to these tabs if the spreadsheet is not `None`, and then closes the spreadsheet. 

        """        
        debugger.print("Start:: Write spreadsheet")
        self.open_excel_spreadsheet()
        if self.spreadsheet is not None:
            self.mainTab.writeSpreadsheet()
            self.settingsTab.writeSpreadsheet()
            self.analysisTab.writeSpreadsheet()
            self.plottingTab.writeSpreadsheet()
            self.spreadsheet.close()
        debugger.print("Finished:: Write spreadsheet")

    def open_excel_spreadsheet(self):
        """Open an Excel spreadsheet based on the filename set in settings.

        This method tries to open an Excel spreadsheet file (.xlsx) whose name is provided in the 'Excel file name' setting of the 'mainTab' attribute. It checks for the validity of the filename (i.e., whether it ends in '.xlsx') and existence in the specified directory. On failure, it alerts the user accordingly.

        Parameters
        ----------
        None

        Returns
        -------
        None


        Raises
        ------
        QMessageBox
            - If the spreadsheet name is not valid (does not end in .xlsx).
            - If the spreadsheet name is empty.

        """        
        debugger.print("Start:: open_spreadsheet clicked")
        if len(self.mainTab.settings["Excel file name"]) > 5 and self.mainTab.settings["Excel file name"][-5:] == ".xlsx":
            self.directory = self.mainTab.directory
            # open the file name with the directory of the output file name
            self.openSpreadSheet(os.path.join(self.directory,self.mainTab.settings["Excel file name"]))
        elif len(self.mainTab.settings["Excel file name"]) > 1 and self.mainTab.settings["Excel file name"][-5:] != ".xlsx":
            # The file isn't valid so tell the user there is a problem
            debugger.print("open_spreadsheet spreadsheet name is not valid",self.mainTab.settings["Excel file name"])
            QMessageBox.about(self,"Spreadsheet name","File name of spreadsheet must end in  .xlsx")
        else:
            debugger.print("open_spreadsheet spreadsheet name is empty")
        debugger.print("Finished:: open_spreadsheet clicked")
        return

    def openSpreadSheet(self,filename):
        """Open or create a spreadsheet file.

        This function checks whether the specified spreadsheet (.xlsx) file exists. If the file exists and overwriting is allowed or confirmed by the user, it opens and overwrites the file. If the file does not exist, it creates a new spreadsheet file. The function also closes any previously opened spreadsheet before attempting to open or create a new one.

        Parameters
        ----------
        filename : str
            The name of the spreadsheet file to open or create.

        Returns
        -------
        None


        Notes
        -----
        - The function relies on the 'SpreadSheetManager' class for handling spreadsheet operations.
        - The function raises no exceptions, but will print a message if the provided filename does not have a '.xlsx' extension.

        """        
        debugger.print("Start:: openSpreadSheet", filename)
        if self.spreadsheet is not None:
            self.spreadsheet.close()
        if filename[-5:] == ".xlsx":
            if os.path.exists(filename):
                debugger.print("Spreadsheet file already exists",self.directory)
                if self.overwriting:
                    debugger.print("Overwriting existing spreadsheet anyway",filename)
                    self.spreadsheet = SpreadSheetManager(filename)
                else:
                    answer = QMessageBox.question(self,"","Spreadsheet already exists.  Continue?", QMessageBox.Yes | QMessageBox.No)
                    if answer == QMessageBox.Yes:
                        debugger.print("Overwriting existing spreadsheet",filename)
                        self.spreadsheet = SpreadSheetManager(filename)
            else:
                debugger.print("Creating a new spreadsheet",filename)
                self.spreadsheet = SpreadSheetManager(filename)
        else:
           print("spreadsheet name not valid", filename)
        debugger.print("Finished:: openSpreadSheet", filename)
        return

    def on_tabs_currentChanged(self, tabindex):
        """Handle tab change events and refresh content accordingly.

        This function responds to changes in the current tab index within a tabbed interface. It refreshes the content of the new current tab based on the index of the tab. This includes refreshing content in tabs corresponding to settings, plotting, analysis, viewing, fitting, or specific scenarios.

        Parameters
        ----------
        tabindex : int
            The index of the newly selected tab.

        Returns
        -------
        None

        Notes
        -----
        - The function first checks if scripting is currently active; if so, it exits early without refreshing to avoid conflicts.
        - The `ntabs` variable calculates the total number of tabs dynamically based on the number of scenarios present.
        - The function determines which tab has been selected based on the `tabindex` and calls the appropriate refresh function for the content of that tab.
        - For predefined tabs (such as settings, plotting, analysis, viewing, and fitting tabs), direct refresh calls are made.
        - For scenario-specific tabs, which are dynamically added based on the number of scenarios, the function calculates the appropriate scenario index and triggers a refresh for the selected scenario.

        """        
        debugger.print("Start:: on_tabs_currentChanged", tabindex)
        # 
        # If scripting do not refresh tabs
        #
        if self.scripting:
            debugger.print("Finished:: Exiting on_tabs_currentChanged without refreshing")
            return
        #       Number of tabs
        ntabs = 2 + len(self.scenarios) + 4
        debugger.print("Number of tabs",ntabs)
        if tabindex == ntabs-1:
            # fitter tab
            debugger.print("Calling fitterTab refresh")
            self.fitterTab.refresh()
        elif tabindex == ntabs-2:
            # viewer tab
            debugger.print("Calling viewerTab refresh")
            self.viewerTab.refresh()
        elif tabindex == ntabs-3:
            # analysis tab
            debugger.print("Calling analysisTab refresh")
            self.analysisTab.refresh()
        elif tabindex == ntabs-4:
            # plottings tab
            debugger.print("Calling plottingTab refresh")
            self.plottingTab.refresh()
        elif tabindex == 1:
            # settings tab
            debugger.print("Calling settingsTab refresh")
            self.settingsTab.refresh()
        else :
            #  Refresh scenario tabs
            scenarioTabIndex = tabindex-2
            if scenarioTabIndex >= 0 and scenarioTabIndex < len(self.scenarios):
                self.scenarios[scenarioTabIndex].refresh()
        debugger.print("Exiting on_tabs_currentChanged()")
        debugger.print("Finished:: on_tabs_currentChanged", tabindex)
        return

    def keyPressEvent(self, e):
        """Handle key press events for the application.

        Parameters
        ----------
        e : QKeyEvent
            An event parameter containing details of the key that was pressed.

        Returns
        -------
        None

        Notes
        -----
        This function checks for specific key combinations (Control + S, and Control + C) and performs actions accordingly:
        - Control + S: Calls the `print_settings` method.
        - Control + C: Prints a message and exits the program.

        """        
        debugger.print("Start:: keyPressEvent")
        if (e.key() == Qt.Key_S)  and QApplication.keyboardModifiers() and Qt.ControlModifier:
            print("Control S has been pressed")
            self.print_settings()
        elif (e.key() == Qt.Key_C)  and QApplication.keyboardModifiers() and Qt.ControlModifier:
            print("Control C has been pressed")
            print("The program will close down")
            sys.exit()
        debugger.print("Finished:: keyPressEvent")
        return

    def progressbars_set_maximum( self, maximum ):
        """Set the maximum value for all progress bars in an object and reset their current status.

        Parameters
        ----------
        maximum : int
            The maximum value to set for each progress bar.

        Returns
        -------
        None
    
        Notes
        -----
        This method sets the maximum value of all progress bars stored in the object's `progressbars` attribute. It also resets the progress to 0. 

        """        
        debugger.print("Start:: progressbars_set_maximum",maximum)
        self.progressbar_status = 0
        self.progressbar_maximum = maximum
        for bar in self.progressbars:
            bar.setMaximum(maximum)
            bar.setValue(self.progressbar_status)
        debugger.print("Finished:: progressbars_set_maximum",maximum)
        return

    def progressbars_update( self, increment=1 ):
        """Update the progress bars status by a specified increment.

        Parameters
        ----------
        increment : int, optional
            The value by which to increment the progress bar status, by default 1.

        Returns
        -------
        None

        This method increments the progress bar status stored in `progressbar_status` by the specified `increment` amount. It then sets this updated value as the new value for all progress bars stored in the `progressbars` list attribute of the instance.

        """        
        self.progressbar_status += increment
        for bar in self.progressbars:
            bar.setValue(self.progressbar_status)
        return

    def progressbars_add( self, bar ):
        """Add a progress bar to the list of progress bars.

        Parameters
        ----------
        bar : object
            The progress bar object to be added.

        Returns
        -------
        None

        Notes
        -----
        After adding the new progress bar to the list, this method updates the maximum value 
        of all progress bars by calling `self.progressbars_set_maximum` with the current 
        maximum value defined in `self.progressbar_maximum`.

        """        
        self.progressbars.append(bar)
        self.progressbars_set_maximum(self.progressbar_maximum)
        return


