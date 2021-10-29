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
        debugger.print('Start:: initialiser')
        self.refreshRequired = True
        self.settings = {}
        self.notebook = parent
        self.settings['Legend'] = 'Unset'
        self.scenarioType = None
        self.settings['Scenario type'] = 'Unset'
        self.vs_cm1 = [0, 0]
        debugger.print('Finished:: initialiser')

    def requestRefresh(self):
        debugger.print('Start:: requestRefresh')
        self.refreshRequired = True
        debugger.print('Finished:: requestRefresh')
        return

    def set_reader(self,reader):
        debugger.print('Start:: set_reader')
        self.refreshRequired = True
        self.reader = reader
        debugger.print('Finished:: set_reader')
        return

    def setScenarioIndex(self,index):
        debugger.print('Start:: setScenarioIndex',index)
        self.scenarioIndex = index
        text = self.settings['Legend']
        if text.startswith('Unset') or text.startswith('Scenario ') or text.startswith('Powder scenario ') or text.startswith('Crystal scenario '):
            debugger.print('setScenarioIndex changing scenario legend','Scenario+str(index+1)')
            self.legend_le.setText('Scenario '+str(index + 1))
            self.settings['Legend'] = 'Scenario '+str(index + 1)
        debugger.print('Finished:: setScenarioIndex',index)
        return

    def print_settings(self):
        debugger.print('Start:: print_settings')
        print('#')
        print('# Scenario tab')
        print('#')
        print('tab = self.notebook.scenarios')
        for key in self.settings:
            print(key, self.settings[key])
        debugger.print('Finished:: print_settings')

    def on_legend_le_changed(self,text):
        debugger.print('on legend change', text)
        self.refreshRequired = True
        self.settings['Legend'] = text
        return

    def  add_scenario_buttons(self):
        """Add a set of scenario buttons in an hbox.  Return the hbox"""
        debugger.print('Start:: add_scenario_buttons')
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
        debugger.print('Finished:: add_scenario_buttons')
        return hbox

    def addScenarioButtonClicked(self):
        # Add another scenario
        debugger.print('addScenarioButtonClicked')
        self.notebook.addScenario(copyFromIndex=self.scenarioIndex)
        return

    def deleteScenarioButtonClicked(self):
        # Delete a scenario
        debugger.print('deleteScenarioButtonClicked')
        self.notebook.deleteScenario(self.scenarioIndex)
        return

    def switchScenarioButtonClicked(self):
        # Switch a scenario
        debugger.print('SwitchScenarioButtonClicked')
        self.notebook.switchScenario(self.scenarioIndex)
        return




