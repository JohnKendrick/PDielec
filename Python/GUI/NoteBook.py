import sys
import os.path
import copy
from PyQt5.QtWidgets            import  QPushButton, QWidget, QTabWidget
from PyQt5.QtWidgets            import  QVBoxLayout
from Python.GUI.MainTab         import  MainTab
from Python.GUI.SettingsTab     import  SettingsTab
from Python.GUI.ScenarioTab     import  ScenarioTab
from Python.GUI.PlottingTab     import  PlottingTab
from Python.GUI.AnalysisTab     import  AnalysisTab
from Python.GUI.ViewerTab       import  ViewerTab
from Python.Utilities           import  Debug
 
class NoteBook(QWidget):        
 
    def __init__(self, parent, program, filename, debug=False):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'NoteBook:')
        self.reader = None
        self.newCalculationRequired = True
        self.debug = debug
        self.layout = QVBoxLayout()
        # The number of tabs before we have scenarios
        self.tabOffSet = 2
        # Initialize tab screen
        self.tabs = QTabWidget(self)
        self.tabs.currentChanged.connect(self.on_tabs_currentChanged)
        self.mainTab = MainTab(self, program, filename, debug=debug)
        self.settingsTab = SettingsTab(self, debug=debug)
        if filename != '':
            self.settingsTab.refresh()
            self.settingsTab.calculateButtonClicked()
        # Open more windows
        self.scenarios = []
        self.scenarios.append( ScenarioTab(self, debug=debug ) )
        self.scenarios[0].setScenarioIndex(0)
        # Open the plotting tab
        self.plottingTab = PlottingTab(self, debug=debug)
        if filename != '':
            self.plottingTab.refresh()
        #
        # Open the Analysis tab
        self.analysisTab = AnalysisTab(self, debug=debug)
        if filename != '':
            self.analysisTab.refresh()
        #
        # Open the Viewer tab
        self.viewerTab = ViewerTab(self, debug=debug)
        if filename != '':
            self.viewerTab.refresh()
        #
        #self.tabs.resize(300,200) 
        # Add tabs
        self.tabs.addTab(self.mainTab,'Main')
        self.tabs.addTab(self.settingsTab,'Settings')
        for i,tab in enumerate(self.scenarios):
            self.tabs.addTab(tab,'Scenario '+str(i+1))
        self.tabs.addTab(self.plottingTab,'Plotting')
        self.tabs.addTab(self.analysisTab,'Analysis')
        self.tabs.addTab(self.viewerTab,'3D Viewer')

        # Add the tab widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        return

    def addScenario(self,index):
        self.scenarios.append( ScenarioTab(self, self.debug) )
        # debugger.print('Settings for scenario', index)
        # self.scenarios[index].print_settings()
        self.scenarios[-1].settings = copy.deepcopy(self.scenarios[index].settings)
        # debugger.print('Settings for new scenario')
        # self.scenarios[-1].print_settings()
        self.scenarios[-1].refresh()
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
        n = len(self.scenarios)
        self.tabs.insertTab(self.tabOffSet+n-1,self.scenarios[-1],'Scenario '+str(n))
        self.tabs.setCurrentIndex(self.tabOffSet+n-1)
        return

    def deleteScenario(self,index):
        # Don't delete the last scenario
        if len(self.scenarios) > 1:
            self.tabs.removeTab(self.tabOffSet+index)
            del self.scenarios[index]
            for i,scenario in enumerate(self.scenarios):
                scenario.setScenarioIndex(i)
                self.tabs.setTabText(self.tabOffSet+i,'Scenario '+str(i+1))
        return

    def on_tabs_currentChanged(self, tabindex):
        debugger.print('Tab index changed', tabindex)
        #       Number of tabs
        ntabs = 2 + len(self.scenarios) + 2
        self.newCalculationRequired = True
        if tabindex == ntabs-1:
            # analysis tab
            self.viewerTab.refresh()
        if tabindex == ntabs-2:
            # analysis tab
            self.analysisTab.refresh()
        elif tabindex == ntabs-3:
            # plottings tab
            self.plottingTab.refresh()
