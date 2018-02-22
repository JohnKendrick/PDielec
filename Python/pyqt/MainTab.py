import sys
import os.path
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtWidgets import  QPushButton, QWidget, QAction, QTabWidget
from PyQt5.QtWidgets import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets import  QFileDialog, QCheckBox
from PyQt5.QtWidgets import  QFormLayout
from PyQt5.QtWidgets import  QVBoxLayout, QHBoxLayout, QMessageBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtCore import Qt
from Python.Utilities import get_reader
import numpy as np
 
class MainTab(QWidget):        
 
    def __init__(self, parent, program, filename):   
        super(QWidget, self).__init__(parent)
        self.settings = {}
        self.settings["program"] = "Castep"
        self.settings["filename"] = filename
        if program != "":
            self.settings["program"] = program.capitalize()
        self.notebook = parent
        self.reader = None
#        self.layout = QVBoxLayout()
 
        # Create first tab - MAIN
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The program combobox
        #
        self.program_cb = QComboBox(self)
        self.program_cb.setToolTip("Choose QM/MM program")
        self.program_cb.addItem('Abinit')
        self.program_cb.addItem('Castep')
        self.program_cb.addItem('Crystal')
        self.program_cb.addItem('Experiment')
        self.program_cb.addItem('Gulp')
        self.program_cb.addItem('Phonopy')
        self.program_cb.addItem('Quantum Espresso')
        self.program_cb.addItem('Vasp')
        index = self.program_cb.findText(self.settings["program"], Qt.MatchFixedString)
        if index >=0:
            self.program_cb.setCurrentIndex(index)
        else:
            self.settings["program"] = "Castep"
            index = 1
            self.program_cb.setCurrentIndex(index)
        self.program_cb.currentIndexChanged.connect(self.on_program_cb_changed)
        form.addRow(QLabel("QM/MM Program:"), self.program_cb)
        #
        # The file selector
        #
        self.file_le = QLineEdit(self)
        self.file_le.setToolTip("Choose output file for analysis (press return for a file chooser)")
        self.file_le.setText(filename)
        self.file_le.returnPressed.connect(self.on_file_le_return)
        self.file_le.textChanged.connect(self.on_file_le_changed)
        form.addRow(QLabel("Output file name:"), self.file_le)
        #
        # Final button
        #
        self.pushButton1 = QPushButton("Read Output File")
        self.pushButton1.setToolTip("Read the output file specified and list the phonon frequencies found")
        self.pushButton1.clicked.connect(self.pushButton1Clicked)
        form.addRow(self.pushButton1)
        vbox.addLayout(form)
        # output window
        self.listw_l = QLabel("Frequencies from "+self.settings["filename"]+":", self)
        vbox.addWidget(self.listw_l)
        self.listw = QListWidget(self)
        vbox.addWidget(self.listw)
        # finalise the layout
        self.setLayout(vbox)
        # If the filename was given the force it to be read and processed
        if filename != "":
            self.pushButton1Clicked()

    def pushButton1Clicked(self):
        print("Reading output file ", self.settings["filename"])
        if not os.path.isfile(self.settings["filename"]):
            QMessageBox.about(self,"Processing output file","The filename for the output file to be processed is not correct: "+self.settings["filename"])
            return
        self.listw_l.setText("Frequencies from "+self.settings["filename"]+":")
        self.reader = get_reader(self.settings["program"],self.settings["filename"])
        # tell the notebook that we have read the info and we have a reader
        self.notebook.reader = self.reader
        self.reader.read_output()
        self.reader.print_info()
        frequencies_cm1 = np.sort(self.reader.frequencies)
        for f in frequencies_cm1:
            self.listw.addItem("{0:.3f}".format(f))
        # tell the settings tab to update the widgets that depend on the contents of the reader
        if hasattr(self.notebook, 'settingsTab'):
            self.notebook.settingsTab.refresh()
        if hasattr(self.notebook, 'scenarios'):
            for tab in self.notebook.scenarios:
                tab.refresh()
        

    def on_born_changed(self):
        print("on born change ", self.born_cb.isChecked())
        self.settings["neutral"] = self.born_cb.isChecked()
        print("on born change ", self.settings["neutral"])

    def on_eckart_changed(self):
        print("on eckart change ", self.eckart_cb.isChecked())
        self.settings["eckart"] = self.eckart_cb.isChecked()
        print("on eckart change ", self.settings["eckart"])

    def on_file_le_return(self):
        print("on file return ", self.file_le.text())
        # Does the file exist?
        self.settings["filename"] = self.file_le.text()
        if os.path.isfile(self.settings["filename"]):
            self.settings["filename"] = os.path.relpath(self.settings["filename"])
            # The file exists to treat it as though the button has been pressed
            self.pushButton1Clicked()
            return
        # The file doesn't exist so open a file chooser
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.settings["filename"], _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        self.settings["filename"] = os.path.relpath(self.settings["filename"])
        if self.settings["filename"]:
            print(self.settings["filename"])
        self.file_le.setText(self.settings["filename"])
 
    def on_file_le_changed(self, text):
        print("on file changed", text)
        self.settings["filename"] = text
 
    def on_program_cb_changed(self,index):
        print("on program combobox changed", index)
        print("on program combobox changed", self.program_cb.currentText())
        self.settings["program"] = self.program_cb.currentText()
