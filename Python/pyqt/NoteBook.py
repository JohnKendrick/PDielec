import sys
import os.path
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
 
class NoteBook(QWidget):        
 
    def __init__(self, parent):   
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout()
 
        # Initialize tab screen
        self.tabs = QTabWidget(self)
        self.tab1 = MainTab(self)
        self.tab2 = SettingsTab(self, self.tab1)
        # tell the main window the tab for the settings
        self.tab1.settingsTab = self.tab2
        # tell the main window what its controlling notebook is
        self.tab1.notebook = self
        # Open more windows
        self.tab3 = QWidget(self)
        self.tabs.resize(300,200) 
 
        # Add tabs
        self.tabs.addTab(self.tab1,"Main")
        self.tabs.addTab(self.tab2,"Settings")
        self.tabs.addTab(self.tab3,"Scenario 1")
 
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        return
