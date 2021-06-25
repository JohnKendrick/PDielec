# -*- coding: utf8 -*-
from PyQt5.QtWidgets   import  QWidget
from PyQt5.QtCore      import  Qt
from PDielec.Utilities import  Debug

class ScenarioTab(QWidget):
    def __init__(self, parent, debug=False):
        """This the base class for all ScenarioTabs"""
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'ScenarioTab:')
        self.dirty = True
        self.settings = {}
        self.notebook = parent
        self.notebook.plottingCalculationRequired = True
        self.notebook.fittingCalculationRequired = True

    def set_reader(self,reader):
        self.dirty = True
        self.notebook.plottingCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.reader = reader

    def setScenarioIndex(self,index):
        self.scenarioIndex = index
        text = self.legend_le.text()
        if text == 'Scenario legend':
            self.legend_le.setText('Scenario '+str(index + 1))
        return

    def print_settings(self):
        print('#')
        print('# Scenario tab')
        print('#')
        print('tab = self.notebook.scenarios')
        for key in self.settings:
            print(key, self.settings[key])

