import sys
import copy
import psutil
from PyQt5.QtWidgets                      import  QWidget, QTabWidget
from PyQt5.QtWidgets                      import  QVBoxLayout
from PyQt5.QtWidgets                      import  QApplication
from PyQt5.QtWidgets                      import  QFileDialog
from PyQt5.QtCore                         import  Qt
from PDielec.GUI.MainTab                  import  MainTab
from PDielec.GUI.SettingsTab              import  SettingsTab
from PDielec.GUI.PowderScenarioTab        import  PowderScenarioTab
from PDielec.GUI.SingleCrystalScenarioTab import  SingleCrystalScenarioTab
from PDielec.GUI.PlottingTab              import  PlottingTab
from PDielec.GUI.AnalysisTab              import  AnalysisTab
from PDielec.GUI.ViewerTab                import  ViewerTab
from PDielec.GUI.FitterTab                import  FitterTab
from PDielec.Utilities                    import  Debug

class NoteBook(QWidget):

    def __init__(self, parent, program, filename, spreadsheet, debug=False, progressbar=None, scripting=False, default_scenario='powder',ncpus=0, threading=False):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'NoteBook:')
        self.reader = None
        self.progressbars=[progressbar]
        if progressbar is None:
            self.progressbars = [ ]
        self.progressbar_status = 0
        self.progressbar_maximum = 0
        self.spreadsheet = None
        self.threading = threading
        if default_scenario == 'powder':
            self.currentScenarioTab = PowderScenarioTab
        else:
            self.currentScenarioTab = SingleCrystalScenarioTab
        if ncpus == 0:
            self.ncpus = psutil.cpu_count(logical=False)
        else:
            self.ncpus = ncpus
        #
        # Set threading for mkl
        #
        try:
            import mkl
            mkl.set_num_threads(self.ncpus)
        except:
            pass
        self.scripting = scripting
        self.debug = debug
        #jk self.old_tab_index = None
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
        if filename != '' and not self.scripting:
            debugger.print('Refreshing settingsTab in notebook initialisation - filename',filename)
            self.settingsTab.refresh()
        #
        # Open more windows
        #
        debugger.print('Initialising the first scenario')
        self.scenarios = []
        self.scenarios.append( self.currentScenarioTab(self, debug=debug ) )
        self.scenarios[0].setScenarioIndex(0)
        self.scenarios[0].settings['Legend'] = 'Scenario 1'
        debugger.print('Finished adding the first scenario')
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
        debugger.print('Initialising the viewer tab')
        self.viewerTab = ViewerTab(self, debug=debug)
        #
        # Open the Fitter tab
        #
        debugger.print('Initialising the fitter tab')
        self.fitterTab = FitterTab(self, debug=debug)
        #
        # Add tabs
        #
        debugger.print('Adding all tabs to the notebook')
        self.tabs.addTab(self.mainTab,'Main')
        self.tabs.addTab(self.settingsTab,'Settings')
        for i,tab in enumerate(self.scenarios):
            tab.refresh(force=True)
            self.tabs.addTab(tab,'Scenario '+str(i+1))
        self.tabs.addTab(self.plottingTab,'Plotting')
        self.tabs.addTab(self.analysisTab,'Analysis')
        self.tabs.addTab(self.viewerTab,'3D Viewer')
        self.tabs.addTab(self.fitterTab,'Fitter')

        # Add the tab widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        return

    def setRefreshRequest(self):
        self.requireRefresh = True
        debugger.print('setRefreshRequest')
        return

    def addScenario(self,scenarioType=None,copyFromIndex=-2):
        """Add Scenario is used by the script to add a new scenario"""
        debugger.print('addScenario for scenarioType', scenarioType,copyFromIndex)
        if copyFromIndex != -2:
            # If the copyFromIndex is not -2 then we override the scenarioType
            last = self.scenarios[copyFromIndex]
            scenarioType = last.scenarioType
        elif scenarioType == None:
            # The default behaviour with no parameters in the call, use the last scenario in the list
            last = self.scenarios[-1]
            scenarioType = last.scenarioType
        else:
            # copyFromIndex is default so we find the last scenario of scenarioType in the list
            last = None
            for scenario in self.scenarios:
                if scenarioType == scenario.scenarioType:
                    last = scenario
            # end for
        # Create a new scenario
        if scenarioType == 'Powder':
            self.currentScenarioTab = PowderScenarioTab
        else:
            self.currentScenarioTab = SingleCrystalScenarioTab
        # Add the scenario to the end of the list
        self.scenarios.append(self.currentScenarioTab(self, self.debug))
        # If we have found a previous scenario of the same time set the settings to it
        if last is not None:
            self.scenarios[-1].settings = copy.deepcopy(last.settings)
            self.scenarios[-1].refresh(force=True)
        n = len(self.scenarios)
        self.tabs.insertTab(self.tabOffSet+n-1,self.scenarios[-1],'Scenario '+str(n))
        self.tabs.setCurrentIndex(self.tabOffSet+n-1)
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
            self.tabs.setTabText(self.tabOffSet+i,'Scenario '+str(i+1))
        return

    def print_settings(self, filename=None):
        # Print the settings of all the settings that have been used to a file settings.py
        debugger.print('print_settings, filename=',filename)
        qf = QFileDialog()
        qf.setWindowTitle('Save the program settings to a file')
        debugger.print('print_settings, directory=',self.mainTab.directory)
        qf.setDirectory(self.mainTab.directory)
        if filename == None:
            filename,selection = qf.getSaveFileName()
        if filename == '':
            return
        print('Current settings will be saved to '+filename)
        fd = open(filename,'w')
        # Handle the special case of the first scenario
        print('#',file=fd)
        print('# Handle the special case of the first scenario',file=fd)
        print('#',file=fd)
        print('self.notebook.switchScenario(0,scenarioType=\"'+self.scenarios[0].scenarioType+'\")',file=fd )
        print('#',file=fd)
        # Print settings of mainTab
        self.print_tab_settings(self.mainTab, 'mainTab',fd)
        print('tab.refresh(force=True)',file=fd)
        # Print settings of settingsTab
        self.print_tab_settings(self.settingsTab, 'settingsTab',fd)
        print('tab.sigmas_cm1 =',self.settingsTab.sigmas_cm1,file=fd)
        print('tab.refresh(force=True)',file=fd)
        # Print settings of scenarios
        for i,tab in enumerate(self.scenarios):
            if i == 0:
                self.print_tab_settings(tab, 'scenarios[{}]'.format(i), fd, new_scenario = False)
            else:
                self.print_tab_settings(tab, 'scenarios[{}]'.format(i), fd, new_scenario = True)
            print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.analysisTab, 'analysisTab',fd)
        print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.viewerTab, 'viewerTab',fd)
        print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.fitterTab, 'fitterTab',fd)
        print('tab.refresh(force=True)',file=fd)
        self.print_tab_settings(self.plottingTab, 'plottingTab',fd)
        print('tab.refresh(force=True)',file=fd)
        fd.close()
        return

    def print_tab_settings(self,tab,title,fd,new_scenario = False):
        print('#',file=fd)
        print('#',file=fd)
        if new_scenario:
            print('self.notebook.addScenario(scenarioType=\"'+self.scenarios[0].scenarioType+'\")',file=fd )
        print('tab = self.notebook.'+title,file=fd)
        for item in tab.settings:
            if item == 'Optical permittivity' and not tab.settings['Optical permittivity edited']:
                    pass
            else:
                value = tab.settings[item]
                if 'str' in str(type(value)):
                    print('tab.settings[\''+item+'\'] = \'{}\''.format(tab.settings[item]),file=fd)
                else:
                    print('tab.settings[\''+item+'\'] = ', tab.settings[item],file=fd)

    def deleteScenario(self,index):
        # Don't delete the last scenario
        if len(self.scenarios) > 1:
            self.tabs.removeTab(self.tabOffSet+index)
            del self.scenarios[index]
            for i,scenario in enumerate(self.scenarios):
                scenario.setScenarioIndex(i)
                self.tabs.setTabText(self.tabOffSet+i,'Scenario '+str(i+1))
            if index-1 < 0:
                index += 1
            self.tabs.setCurrentIndex(self.tabOffSet+index-1)
        return

    def switchScenario(self,index,scenarioType=None):
        debugger.print('switch for scenario', index+1)
        # Replace the scenario with the other scenario type
        scenario = self.scenarios[index]
        debugger.print('Current scenario type', scenario.scenarioType, scenarioType)
        #
        # If scenarioType is specified in the call then force that type
        # Otherwise switch type
        #
        if scenarioType == None:
            if scenario.scenarioType == 'Powder':
                self.currentScenarioTab = SingleCrystalScenarioTab
            else:
                self.currentScenarioTab = PowderScenarioTab
            # end if
        else:
            if scenarioType == 'Powder':
                self.currentScenarioTab = PowderScenarioTab
            else:
                self.currentScenarioTab = SingleCrystalScenarioTab
            # end if
        #end if
        self.scenarios[index] =  self.currentScenarioTab(self, self.debug)
        scenario = self.scenarios[index]
        debugger.print('Current scenario type now', scenario.scenarioType)
        self.tabs.removeTab(self.tabOffSet+index)
        self.tabs.insertTab(self.tabOffSet+index,scenario,'Scenario '+str(index+1) )
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
            self.tabs.setTabText(self.tabOffSet+i,'Scenario '+str(i+1))
        self.tabs.setCurrentIndex(self.tabOffSet+index)
        self.scenarios[index].refresh(force=True)
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

    def writeSpreadsheet(self):
        debugger.print('Write spreadsheet')
        self.mainTab.writeSpreadsheet()
        self.settingsTab.writeSpreadsheet()
        self.analysisTab.writeSpreadsheet()
        self.plottingTab.writeSpreadsheet()

    def on_tabs_currentChanged(self, tabindex):
        debugger.print('Tab index changed', tabindex)
        debugger.print('on_tabs_currentChanged requireRefresh=',self.scenarios[0].requireRefresh)
        #jk if tabindex == 3:
        #jk     import pdb
        #jk     from PyQt5.QtCore import pyqtRemoveInputHook
        #jk     pyqtRemoveInputHook()
        #jk     pdb.set_trace()
        # 
        # If scripting do not refresh tabs
        #
        if self.scripting:
            return
        # See if we have to update a tab we have left
        #jk if self.old_tab_index is not None:
        #jk     if self.old_tab_index == 0:
        #jk         self.mainTab.refresh()
        #jk     elif self.old_tab_index == 1:
        #jk         self.settingsTab.refresh()
        #jk # end if
        #       Number of tabs
        ntabs = 2 + len(self.scenarios) + 4
        if tabindex == ntabs-1:
            # fitter tab
            debugger.print('Calling fitterTab refresh')
            self.fitterTab.refresh()
        elif tabindex == ntabs-2:
            # viewer tab
            debugger.print('Calling viewerTab refresh')
            self.viewerTab.refresh()
        elif tabindex == ntabs-3:
            # analysis tab
            debugger.print('Calling analysisTab refresh')
            self.analysisTab.refresh()
        elif tabindex == ntabs-4:
            # plottings tab
            debugger.print('Calling plottingTab refresh')
            debugger.print('requireRefresh value before call',self.scenarios[0].requireRefresh)
            self.plottingTab.refresh()
        #jk self.old_tab_index = tabindex

    def keyPressEvent(self, e):
        if (e.key() == Qt.Key_S)  and QApplication.keyboardModifiers() and Qt.ControlModifier:
            print('Control S has been pressed')
            self.print_settings()
        elif (e.key() == Qt.Key_C)  and QApplication.keyboardModifiers() and Qt.ControlModifier:
            print('Control C has been pressed')
            print('The program will close down')
            sys.exit()

    def progressbars_set_maximum( self, maximum ):
        self.progressbar_status = 0
        self.progressbar_maximum = maximum
        for bar in self.progressbars:
            bar.setMaximum(maximum)
            bar.setValue(self.progressbar_status)
        return

    def progressbars_update( self, increment=1 ):
        self.progressbar_status += increment
        for bar in self.progressbars:
            bar.setValue(self.progressbar_status)
        return

    def progressbars_add( self, bar ):
        self.progressbars.append(bar)
        bar.setMaximum(self.progressbar_maximum)
        bar.setValue(self.progressbar_status)
        return


