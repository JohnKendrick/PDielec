import os.path
import sys
import copy
from PyQt5.QtWidgets            import  QWidget, QTabWidget
from PyQt5.QtWidgets            import  QVBoxLayout
from PyQt5.QtWidgets            import  QApplication
from PyQt5.QtWidgets            import  QFileDialog
from PyQt5.QtCore               import  Qt
from PDielec.GUI.MainTab         import  MainTab
from PDielec.GUI.SettingsTab     import  SettingsTab
from PDielec.GUI.ScenarioTab     import  ScenarioTab
from PDielec.GUI.PlottingTab     import  PlottingTab
from PDielec.GUI.AnalysisTab     import  AnalysisTab
from PDielec.GUI.ViewerTab       import  ViewerTab
from PDielec.GUI.FitterTab       import  FitterTab
from PDielec.Utilities           import  Debug

class NoteBook(QWidget):

    def __init__(self, parent, program, filename, spreadsheet, debug=False, progressbar=False, scripting=False, ncpus=0, threading=False):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'NoteBook:')
        self.reader = None
        self.progressbar=progressbar
        self.spreadsheet = None
        self.threading = threading
        self.ncpus = ncpus
        self.scripting = scripting
        self.debug = debug
        self.plottingCalculationRequired = True
        self.analysisCalculationRequired = True
        self.visualerCalculationRequired = True
        self.fittingCalculationRequired = True
        self.old_tab_index = None
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
        if filename != '' and not self.scripting:
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
        if filename != '' and not self.scripting:
            debugger.print('Refreshing plotting because filename is set')
            self.plottingTab.refresh()
        #
        # Open the Analysis tab
        #
        self.analysisTab = AnalysisTab(self, debug=debug)
        if filename != '' and not self.scripting:
            debugger.print('Refreshing analysis because filename is set')
            self.analysisTab.refresh()
        #
        # Open the Viewer tab
        #
        self.viewerTab = ViewerTab(self, debug=debug)
        #
        # Open the Fitter tab
        #
        self.fitterTab = FitterTab(self, debug=debug)
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
        self.tabs.addTab(self.fitterTab,'Fitter')

        # Add the tab widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        return

    def addScenario(self,copyFromIndex=-2):
        debugger.print('Settings for scenario', copyFromIndex)
        self.plottingCalculationRequired = True
        self.scenarios.append( ScenarioTab(self, self.debug) )
        self.scenarios[-1].settings = copy.deepcopy(self.scenarios[copyFromIndex].settings)
        # debugger.print('Settings for new scenario')
        self.scenarios[-1].refresh(force=True)
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
        n = len(self.scenarios)
        self.tabs.insertTab(self.tabOffSet+n-1,self.scenarios[-1],'Scenario '+str(n))
        self.tabs.setCurrentIndex(self.tabOffSet+n-1)
        return

    def print_settings(self):
        # Print the settings of all the settings that have been used to a file settings.py
        qf = QFileDialog()
        qf.setWindowTitle('Save the program settings to a file')
        filename,selection = qf.getSaveFileName()
        if filename == '':
            return
        print('Current settings will be saved to '+filename)
        fd = open(filename,'w')
        ntabs = 2+ len(self.scenarios) + 4
        self.print_tab_settings(self.mainTab, 'mainTab',fd)
        print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.settingsTab, 'settingsTab',fd)
        print('tab.sigmas_cm1 =',self.settingsTab.sigmas_cm1,file=fd)
        print('tab.refresh(force=True)',file=fd)
        requireNewScenario = False
        for i,tab in enumerate(self.scenarios):
            self.print_tab_settings(tab, 'scenarios[{}]'.format(i), fd, new_scenario = requireNewScenario)
            print('tab.refresh(force=True)',file=fd)
            requireNewScenario = True
        self.print_tab_settings(self.plottingTab, 'plottingTab',fd)
        print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.analysisTab, 'analysisTab',fd)
        print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.viewerTab, 'viewerTab',fd)
        print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.fitterTab, 'fitterTab',fd)
        print('tab.refresh(force=True)',file=fd)
        fd.close()
        return

    def print_tab_settings(self,tab,title,fd,new_scenario = False):
        print('#',file=fd)
        print('#',file=fd)
        if new_scenario:
            print('self.notebook.addScenario()',file=fd)
        print('tab = self.notebook.'+title,file=fd)
        for item in tab.settings:
            value = tab.settings[item]
            if 'str' in str(type(value)):
                print('tab.settings[\''+item+'\'] = \'{}\''.format(tab.settings[item]),file=fd)
            else:
                print('tab.settings[\''+item+'\'] = ', tab.settings[item],file=fd)

    def deleteScenario(self,index):
        # Don't delete the last scenario
        if len(self.scenarios) > 1:
            self.plottingCalculationRequired = True
            self.tabs.removeTab(self.tabOffSet+index)
            del self.scenarios[index]
            for i,scenario in enumerate(self.scenarios):
                scenario.setScenarioIndex(i)
                self.tabs.setTabText(self.tabOffSet+i,'Scenario '+str(i+1))
        return

    def refresh(self,force=False):
        if self.scripting:
            debugger.print('Notebook aborting refresh because of scripting')
            return
        debugger.print('Notebook refresh changed',force)
        ntabs = 2 + len(self.scenarios) + 4
        self.mainTab.refresh(force=force)
        self.settingsTab.refresh(force=force)
        for tab in self.scenarios:
            tab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-5)
        self.plottingTab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-4)
        self.analysisTab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-3)
        self.viewerTab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-2)
        self.fitterTab.refresh(force=force)
        self.tabs.setCurrentIndex(ntabs-1)

    def write_spreadsheet(self):
        debugger.print('Write spreadsheet')
        self.mainTab.write_spreadsheet()
        self.settingsTab.write_spreadsheet()
        self.analysisTab.write_spreadsheet()
        self.plottingTab.write_spreadsheet()

    def on_tabs_currentChanged(self, tabindex):
        debugger.print('Tab index changed', tabindex)
        if self.scripting:
            return
        # See if we have to up date a tab we have left
        if self.old_tab_index is not None:
            if self.old_tab_index == 0:
                self.mainTab.refresh()
            elif self.old_tab_index == 1:
                self.settingsTab.refresh()
        # end if
        #       Number of tabs
        ntabs = 2 + len(self.scenarios) + 4
        if tabindex == ntabs-1:
            # fitter tab
            self.fitterTab.refresh()
        elif tabindex == ntabs-2:
            # viewer tab
            self.viewerTab.refresh()
        elif tabindex == ntabs-3:
            # analysis tab
            self.analysisTab.refresh()
        elif tabindex == ntabs-4:
            # plottings tab
            self.plottingTab.refresh()
        self.old_tab_index = tabindex

    def keyPressEvent(self, e):
        if (e.key() == Qt.Key_S)  and QApplication.keyboardModifiers() and Qt.ControlModifier:
            print('Control S has been pressed')
            self.print_settings()
        elif (e.key() == Qt.Key_C)  and QApplication.keyboardModifiers() and Qt.ControlModifier:
            print('Control C has been pressed')
            print('The program will close down')
            sys.exit()
