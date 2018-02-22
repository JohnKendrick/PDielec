import sys
import os.path
import copy
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtWidgets import  QPushButton, QWidget, QAction, QTabWidget
from PyQt5.QtWidgets import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets import  QFileDialog, QCheckBox
from PyQt5.QtWidgets import  QVBoxLayout, QHBoxLayout
from PyQt5.QtGui     import  QIcon
from PyQt5.QtCore    import  pyqtSlot
from PyQt5.QtCore    import  Qt
from MainTab         import  MainTab
from SettingsTab     import  SettingsTab
from ScenarioTab     import  ScenarioTab
from PlottingTab     import  PlottingTab
 
class NoteBook(QWidget):        
 
    def __init__(self, parent, program, filename):   
        super(QWidget, self).__init__(parent)
        self.reader = None
        self.layout = QVBoxLayout()
        # The number of tabs before we have scenarios
        self.tabOffSet = 2
        # Initialize tab screen
        self.tabs = QTabWidget(self)
        self.mainTab = MainTab(self, program, filename)
        self.settingsTab = SettingsTab(self)
        if filename != "":
            self.settingsTab.refresh()
            self.settingsTab.calculateButtonClicked()
        # Open more windows
        self.scenarios = []
        self.scenarios.append( ScenarioTab(self, self.mainTab, self.settingsTab ) )
        self.scenarios[0].setScenarioIndex(0)
        # Open the plotting tab
        self.plottingTab = PlottingTab(self)
        #self.tabs.resize(300,200) 
 
        # Add tabs
        self.tabs.addTab(self.mainTab,"Main")
        self.tabs.addTab(self.settingsTab,"Settings")
        for i,tab in enumerate(self.scenarios):
            self.tabs.addTab(tab,"Scenario "+str(i+1))
        # tell the main tab about the scenarios
        self.mainTab.scenarios = self.scenarios
        self.tabs.addTab(self.plottingTab,"Plotting")

        # Add the tab widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        return

    def addScenario(self,index):
        self.scenarios.append( ScenarioTab(self, self.mainTab, self.settingsTab) )
        self.scenarios[-1].settings = copy.deepcopy(self.scenarios[index].settings)
        self.scenarios[-1].refresh()
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
        n = len(self.scenarios)
        self.tabs.insertTab(self.tabOffSet+n-1,self.scenarios[-1],"Scenario "+str(n))
        self.tabs.setCurrentIndex(self.tabOffSet+n-1)
        return

    def deleteScenario(self,index):
        self.tabs.removeTab(self.tabOffSet+index)
        del self.scenarios[index]
        for i,scenario in enumerate(self.scenarios):
            scenario.setScenarioIndex(i)
            self.tabs.setTabText(self.tabOffSet+i,"Scenario "+str(i))
        return
