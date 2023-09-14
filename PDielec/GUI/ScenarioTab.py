# -*- coding: utf8 -*-
from PyQt5.QtWidgets   import  QWidget, QHBoxLayout, QPushButton, QLabel, QLineEdit
from PyQt5.QtWidgets   import  QFileDialog
from PyQt5.QtCore      import  Qt
from PDielec.Utilities import  Debug
from PDielec           import __file__  as PDielec_init_filename
from PDielec.Materials import MaterialsDataBase
from PDielec.Materials import Material
import os

class ScenarioTab(QWidget):
    def __init__(self, parent, debug=False):
        """This the base class for all ScenarioTabs"""
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'ScenarioTab:')
        debugger.print('Start:: initialiser')
        self.refreshRequired = True
        self.noCalculationsRequired = 0
        self.settings = {}
        self.notebook = parent
        self.settings['Legend'] = 'Unset'
        self.scenarioType = None
        self.settings['Scenario type'] = 'Unset'
        self.vs_cm1 = [0, 0]
        # Deal with the Materials Database here as it is used in all Scenarios
        PDielec_Directory = os.path.dirname(PDielec_init_filename)
        self.settings['Materials database'] = os.path.join(PDielec_Directory, 'MaterialsDataBase.xlsx')
        # Open the database
        self.DataBase = MaterialsDataBase(self.settings['Materials database'],debug=debug)
        # Set up the open database button
        self.openDB_button = QPushButton("Open materials' database")
        self.openDB_button.clicked.connect(self.openDB_button_clicked)
        self.openDB_button.setToolTip("Open a new materials' database (.xlsx file)")
        self.openDB_label = QLabel("Open materials' database")
        self.openDB_label.setToolTip("Open a new materials' database (.xlsx file)")
        # set up the database information line
        self.database_le = QLineEdit(self)
        self.database_le.setToolTip("Provides information about the name of the materials' database")
        self.database_le.setText(self.settings['Materials database'])
        self.database_le.setReadOnly(True)
        self.database_le_label = QLabel("Current materials' database")
        self.database_le_label.setToolTip("Provides information about the name of the materials' database")
        debugger.print('Finished:: initialiser')

    def openDataBase(self):
        '''Open the database and set the material names'''
        selfilter = 'Spreadsheet (*.xlsx)'
        filename,myfilter = QFileDialog.getOpenFileName(self,'Open spreadsheet','','Spreadsheet (*.xls);;Spreadsheet (*.xlsx);;All files(*)',selfilter)
        # Process the filename
        if filename == '':
            return
        oldDataBase = self.DataBase
        self.DataBase = MaterialsDataBase(self.settings['Materials database'],debug=debugger.state())
        sheets = self.DataBase.getSheetNames()
        if not self.DataBase.valid():
            self.DataBase = oldDataBase
            print('Error chosen file is not a materials database',sheets)
            return
        self.settings['Materials database'] = filename
        self.database_le.setText(self.settings['Materials database'])
        self.materialNames = sheets
        return

    def getNoCalculationsRequired(self):
        """Get the number of spectra that need recalculating from this scenario"""
        if self.refreshRequired:
            result = self.noCalculationsRequired
        else:
            result = 0
        debugger.print(self.settings['Legend'], 'getNoCalculationsRequired',result)
        return result

    def requestRefresh(self):
        """Request a refresh of the scenario"""
        debugger.print(self.settings['Legend'], 'requestRefresh')
        self.refreshRequired = True
        return

    def set_reader(self,reader):
        '''Set the reader associated with this scenario'''
        debugger.print(self.settings['Legend'], 'set_reader')
        self.refreshRequired = True
        self.reader = reader
        return

    def setScenarioIndex(self,index):
        '''Set the index for the current scenario.  Also set a default legend name based on the index'''
        debugger.print(self.settings['Legend'], 'setScenarioIndex',index)
        self.scenarioIndex = index
        text = self.settings['Legend']
        if text.startswith('Unset') or text.startswith('Scenario ') or text.startswith('Powder scenario ') or text.startswith('Crystal scenario '):
            debugger.print(self.settings['Legend'], 'setScenarioIndex changing scenario legend','Scenario '+str(index+1))
            self.legend_le.setText('Scenario '+str(index + 1))
            self.settings['Legend'] = 'Scenario '+str(index + 1)
        return

    def print_settings(self):
        ''' Print all the settings for this scenario'''
        debugger.print(self.settings['Legend'], 'print_settings')
        print('#')
        print('# Scenario tab')
        print('#')
        print('tab = self.notebook.scenarios')
        for key in self.settings:
            print(key, self.settings[key])

    def on_legend_le_changed(self,text):
        ''' Handle a change in the scenario legend'''
        debugger.print(self.settings['Legend'], 'on_legend_le_changed',text)
        # self.refreshRequired = True
        self.settings['Legend'] = text
        return

    def  add_scenario_buttons(self):
        """Add a set of scenario buttons in an hbox.  Return the hbox"""
        debugger.print(self.settings['Legend'], 'add_scenario_buttons start')
        hbox = QHBoxLayout()
        self.addScenarioButton = QPushButton('Add another scenario')
        self.addScenarioButton.setToolTip('Add another scenario the the notebook tabs, \nthe new scenario is added to the end of the current tab list')
        self.addScenarioButton.clicked.connect(self.addScenarioButtonClicked)
        hbox.addWidget(self.addScenarioButton)
        self.deleteScenarioButton = QPushButton('Delete this scenario')
        self.deleteScenarioButton.setToolTip('Delete the current scenario')
        self.deleteScenarioButton.clicked.connect(self.deleteScenarioButtonClicked)
        hbox.addWidget(self.deleteScenarioButton)
        if self.scenarioType == 'Powder':
            self.switchScenarioButton = QPushButton('Switch to crystal scenario')
            self.switchScenarioButton.setToolTip('Switch the current scenario to a single crystal scenario')
        else:
            self.switchScenarioButton = QPushButton('Switch to powder scenario')
            self.switchScenarioButton.setToolTip('Switch the current scenario to a powder scenario')
        self.switchScenarioButton.clicked.connect(self.switchScenarioButtonClicked)
        hbox.addWidget(self.switchScenarioButton)
        debugger.print(self.settings['Legend'], 'add_scenario_buttons finish')
        return hbox

    def addScenarioButtonClicked(self):
        ''' Handle add another scenario button has been clicked'''
        # Add another scenario
        debugger.print(self.settings['Legend'], 'addScenarioButtonClicked')
        self.notebook.addScenario(copyFromIndex=self.scenarioIndex)
        return

    def deleteScenarioButtonClicked(self):
        ''' Handle switch a delete button has been clicked'''
        debugger.print(self.settings['Legend'], 'deleteScenarioButtonClicked')
        self.notebook.deleteScenario(self.scenarioIndex)
        return

    def switchScenarioButtonClicked(self):
        ''' Handle switch a scenario button has been clicked'''
        debugger.print(self.settings['Legend'],'switchScenarioButtonClicked')
        self.notebook.switchScenario(self.scenarioIndex)
        return




