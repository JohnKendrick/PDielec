import sys
import os.path
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtWidgets import  QPushButton, QWidget, QAction, QTabWidget
from PyQt5.QtWidgets import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets import  QFileDialog, QCheckBox
from PyQt5.QtWidgets import  QVBoxLayout, QHBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtCore import Qt
 
class NoteBook(QWidget):        
 
    def __init__(self, parent):   
        super(QWidget, self).__init__(parent)
        self.settings = {}
        self.settings["eckart"] = True
        self.settings["born"] = False
        self.settings["program"] = "Castep"
        self.settings["filename"] = ""
        self.layout = QVBoxLayout()
 
        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tab1 = MainTab()	
#        self.tab2 = QWidget()
#        self.tab3 = QWidget()
#        self.tabs.resize(300,200) 
 
        # Add tabs
        self.tabs.addTab(self.tab1,"Main")
#        self.tabs.addTab(self.tab2,"Settings")
#        self.tabs.addTab(self.tab3,"Scenario 1")
 
        self.setLayout(self.layout)
        return
