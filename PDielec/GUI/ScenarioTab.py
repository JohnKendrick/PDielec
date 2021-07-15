# -*- coding: utf8 -*-
from PyQt5.QtWidgets   import  QWidget, QHBoxLayout, QPushButton
from PyQt5.QtCore      import  Qt
from PDielec.Utilities import  Debug

class ScenarioTab(QWidget):
    def __init__(self, parent, debug=False):
        """This the base class for all ScenarioTabs"""
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'ScenarioTab:')
        debugger.print('In the initialiser')
        self.requireRefresh = True
        self.settings = {}
        self.notebook = parent
        self.settings['Legend'] = 'Unset'
        self.scenarioType = None
        self.settings['Scenario type'] = 'Unset'
        self.vs_cm1 = [0, 0]

    def set_reader(self,reader):
        self.requireRefresh = True
        self.reader = reader

    def setScenarioIndex(self,index):
        self.scenarioIndex = index
        text = self.legend_le.text()
        if text.startswith('Unset') or text.startswith('Scenario ') or text.startswith('Powder scenario ') or text.startswith('Crystal scenario '):
            self.legend_le.setText('Scenario '+str(index + 1))
        return

    def print_settings(self):
        print('#')
        print('# Scenario tab')
        print('#')
        print('tab = self.notebook.scenarios')
        for key in self.settings:
            print(key, self.settings[key])

    def on_legend_le_changed(self,text):
        debugger.print('on legend change', text)
        self.requireRefresh = True
        self.settings['Legend'] = text

    def  add_scenario_buttons(self):
        """Add a set of scenario buttons in an hbox.  Return the hbox"""
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
        return hbox

    def addScenarioButtonClicked(self):
        # Add another scenario
        debugger.print('addScenarioButtonClicked')
        self.notebook.addScenario(copyFromIndex=self.scenarioIndex)

    def deleteScenarioButtonClicked(self):
        # Delete a scenario
        debugger.print('deleteScenarioButtonClicked')
        self.notebook.deleteScenario(self.scenarioIndex)

    def switchScenarioButtonClicked(self):
        # Switch a scenario
        debugger.print('SwitchScenarioButtonClicked')
        self.notebook.switchScenario(self.scenarioIndex)




