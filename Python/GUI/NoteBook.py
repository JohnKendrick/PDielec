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
 
    def __init__(self, parent, program, filename, spreadsheet, debug=False, progressbar=False, script=False):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'NoteBook:')
        self.reader = None
        self.progressbar=progressbar
        self.spreadsheet = None
        self.newPlottingCalculationRequired = True
        self.newAnalysisCalculationRequired = True
        self.debug = debug
        self.script = script
        self.layout = QVBoxLayout()
        # The number of tabs before we have scenarios
        self.tabOffSet = 2
        #
        # Initialize tab screen
        #
        self.tabs = QTabWidget(self)
        self.tabs.currentChanged.connect(self.on_tabs_currentChanged)
        self.mainTab = MainTab(self, program, filename, spreadsheet, debug=debug)
        self.settingsTab = SettingsTab(self, debug=debug)
        if filename != '' and not self.script:
            debugger.print('Refreshing settingsTab in notebook initialisation - filename',filename)
            self.settingsTab.refresh()
        #
        # Open more windows
        #
        self.scenarios = []
        self.scenarios.append( ScenarioTab(self, debug=debug ) )
        self.scenarios[0].setScenarioIndex(0)
        #
        # Open the plotting tab
        #
        self.plottingTab = PlottingTab(self, debug=debug)
        if filename != '' and not self.script:
            debugger.print('Refreshing plotting because filename is set')
            self.plottingTab.refresh()
        #
        # Open the Analysis tab
        #
        self.analysisTab = AnalysisTab(self, debug=debug)
        if filename != '' and not self.script:
            debugger.print('Refreshing analysis because filename is set')
            self.analysisTab.refresh()
        #
        # Open the Viewer tab
        #
        self.viewerTab = ViewerTab(self, debug=debug)
        #
        # Add tabs
        #
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

    def addScenario(self,copyFromIndex=-2):
        debugger.print('Settings for scenario', copyFromIndex)
        self.newPlottingCalculationRequired = True
        self.scenarios.append( ScenarioTab(self, self.debug) )
        self.scenarios[-1].settings = copy.deepcopy(self.scenarios[copyFromIndex].settings)
        # debugger.print('Settings for new scenario')
        # self.scenarios[-1].print_settings()
        self.scenarios[-1].refresh(force=True)
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
        n = len(self.scenarios)
        self.tabs.insertTab(self.tabOffSet+n-1,self.scenarios[-1],'Scenario '+str(n))
        self.tabs.setCurrentIndex(self.tabOffSet+n-1)
        return

    def deleteScenario(self,index):
        # Don't delete the last scenario
        if len(self.scenarios) > 1:
            self.newPlottingCalculationRequired = True
            self.tabs.removeTab(self.tabOffSet+index)
            del self.scenarios[index]
            for i,scenario in enumerate(self.scenarios):
                scenario.setScenarioIndex(i)
                self.tabs.setTabText(self.tabOffSet+i,'Scenario '+str(i+1))
        return

    def refresh(self,force=False):
        if self.script:
            debugger.print('Notebook aborting refresh because of scripting')
            return
        debugger.print('Notebook refresh changed',force)
        ntabs = 2 + len(self.scenarios) + 3
        self.mainTab.refresh(force=force)
        self.settingsTab.refresh(force=force)
        for tab in self.scenarios:
            tab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-3)
        self.plottingTab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-3)
        self.analysisTab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-2)
        self.viewerTab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-1)

    def on_tabs_currentChanged(self, tabindex):
        debugger.print('Tab index changed', tabindex)
        if self.script:
            return
        #       Number of tabs
        ntabs = 2 + len(self.scenarios) + 3
        if tabindex == ntabs-1:
            # viewer tab
            self.viewerTab.refresh()
        elif tabindex == ntabs-2:
            # analysis tab
            self.analysisTab.refresh()
        elif tabindex == ntabs-3:
            # plottings tab
            self.plottingTab.refresh()
