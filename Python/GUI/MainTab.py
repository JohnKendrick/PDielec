import sys
import os.path
from PyQt5.QtWidgets import  QWidget
from PyQt5.QtWidgets import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets import  QFileDialog, QPushButton
from PyQt5.QtWidgets import  QFormLayout
from PyQt5.QtWidgets import  QVBoxLayout, QHBoxLayout, QMessageBox
from PyQt5.QtCore import Qt
from Python.Utilities import get_reader, Debug
import numpy as np
 
class MainTab(QWidget):        
 
    def __init__(self, parent, program, filename, debug=False):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug, 'MainTab:')
        self.debug = debug
        self.settings = {}
        self.settings['program'] = 'Castep'
        self.settings['filename'] = filename
        self.settings['qmprogram'] = ''
        if program != '':
            self.settings['program'] = program.lower()
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
        self.program_cb.setToolTip('Choose QM/MM program')
        self.program_cb.addItem('Abinit')
        self.program_cb.addItem('Castep')
        self.program_cb.addItem('Crystal')
        self.program_cb.addItem('Vasp')
        self.program_cb.addItem('Gulp')
        self.program_cb.addItem('Phonopy - VASP')
#        self.program_cb.addItem('Phonopy - QE') 
#        self.program_cb.addItem('Phonopy - Crystal')
#        self.program_cb.addItem('Quantum Espresso')
        self.program_cb.addItem('Experiment')
        prtext = self.settings['program'].capitalize()
        qmtext = self.settings['qmprogram'].capitalize()
        if qmtext == 'Qe':
            qmtext = 'QE'
        if prtext == 'Phonopy':
            if qmtext != '':
                prtext = prtext + ' - ' + qmtext
        elif prtext == 'QE':
            prtext = 'Quantum Espresso'
        index = self.program_cb.findText(prtext, Qt.MatchFixedString)
        if index >=0:
            self.program_cb.setCurrentIndex(index)
        else:
            self.settings['program'] = 'Castep'
            index = 1
            self.program_cb.setCurrentIndex(index)
        self.program_cb.currentIndexChanged.connect(self.on_program_cb_changed)
        label = QLabel('QM/MM Program ')
        label.setToolTip('Choose QM/MM program')
        form.addRow(label, self.program_cb)
        #
        # The file selector
        #
        self.file_le = QLineEdit(self)
        self.file_le.setToolTip('Choose output file for analysis (press return for a file chooser)')
        self.file_le.setText(filename)
        self.file_le.returnPressed.connect(self.on_file_le_return)
        self.file_le.textChanged.connect(self.on_file_le_changed)
        label = QLabel('Output file name:')
        label.setToolTip('Choose output file for analysis (press return for a file chooser)')
        form.addRow(label, self.file_le)
        #
        # Final button
        #
        #self.pushButton1 = QPushButton('Read Output File')
        #self.pushButton1.setToolTip('Read the output file specified and list the phonon frequencies found')
        #self.pushButton1.clicked.connect(self.pushButton1Clicked)
        #form.addRow(self.pushButton1)
        vbox.addLayout(form)
        # output window
        self.listw_l = QLabel('Frequencies from '+self.settings['filename']+':', self)
        vbox.addWidget(self.listw_l)
        self.listw = QListWidget(self)
        vbox.addWidget(self.listw)
        # finalise the layout
        self.setLayout(vbox)
        # If the filename was given the force it to be read and processed
        if filename != '':
            self.pushButton1Clicked()

    def pushButton1Clicked(self):
        debugger.print('Reading output file ', self.settings['filename'])
        if not os.path.isfile(self.settings['filename']):
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['filename'])
            return
        self.listw_l.setText('Frequencies from '+self.settings['filename']+':')
        self.reader = get_reader(self.settings['program'],[ self.settings['filename'] ], self.settings['qmprogram'] )
        if self.reader is None:
            print('Error in reading files - program  is ',self.settings['program'])
            print('Error in reading files - filename is ',self.settings['filename'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['filename'])
            return
        #self.reader.debug = self.debug
        # tell the notebook that we have read the info and we have a reader
        self.notebook.reader = self.reader
        try:
            self.reader.read_output()
        except:
            print('Error in reading output files - program  is ',self.settings['program'])
            print('Error in reading output files - filename is ',self.settings['filename'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['filename'])
            return
        if self.debug:
            self.reader.print_info()
        frequencies_cm1 = np.sort(self.reader.frequencies)
        for f in frequencies_cm1:
            self.listw.addItem('{0:.3f}'.format(f))
        # tell the settings tab to update the widgets that depend on the contents of the reader
        debugger.print('Main tab: processing a return')
        if hasattr(self.notebook, 'settingsTab'):
            self.notebook.settingsTab.refresh()
        if hasattr(self.notebook, 'scenarios'):
            debugger.print('notebook has {} scenarios'.format(len(self.notebook.scenarios)))
            for tab in self.notebook.scenarios:
                tab.refresh()
        else:
            debugger.print('notebook has no scenarios yet')
        if hasattr(self.notebook, 'plottingTab'):
            debugger.print('notebook has a plottingTab')
            self.notebook.plottingTab.refresh()
        else:
            debugger.print('notebook has no plottab yet')

    def on_born_changed(self):
        debugger.print('on born change ', self.born_cb.isChecked())
        self.settings['neutral'] = self.born_cb.isChecked()
        debugger.print('on born change ', self.settings['neutral'])

    def on_eckart_changed(self):
        debugger.print('on eckart change ', self.eckart_cb.isChecked())
        self.settings['eckart'] = self.eckart_cb.isChecked()
        debugger.print('on eckart change ', self.settings['eckart'])

    def on_file_le_return(self):
        debugger.print('on file return ', self.file_le.text())
        # Does the file exist?
        self.settings['filename'] = self.file_le.text()
        if os.path.isfile(self.settings['filename']):
            # self.settings['filename'] = os.path.relpath(self.settings['filename'])
            # The file exists to treat it as though the button has been pressed
            self.pushButton1Clicked()
            return
        # The file doesn't exist so open a file chooser
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.settings['filename'], _ = QFileDialog.getOpenFileName(self,'QFileDialog.getOpenFileName()', '','All Files (*)', options=options)
        try:
            #self.settings['filename'] = os.path.relpath(self.settings['filename'])
            pass
        except:
            self.settings['filename'] = ''
        if self.settings['filename']:
            debugger.print('new file name', self.settings['filename'])
        self.file_le.setText(self.settings['filename'])
        if os.path.isfile(self.settings['filename']):
            #self.settings['filename'] = os.path.relpath(self.settings['filename'])
            # The file exists to treat it as though the button has been pressed
            self.pushButton1Clicked()
            return
 
    def on_file_le_changed(self, text):
        debugger.print('on file changed', text)
        self.settings['filename'] = text
 
    def on_program_cb_changed(self,index):
        debugger.print('on program combobox changed', index)
        debugger.print('on program combobox changed', self.program_cb.currentText())
        text = self.program_cb.currentText()
        self.settings['program'] = text
        self.settings['qmprogram'] = ''
        if text == 'Phonopy - VASP':
            self.settings['program']   = 'phonopy'
            self.settings['qmprogram'] = 'vasp'
        elif text == 'Phonopy - QE':
            self.settings['program']   = 'phonopy'
            self.settings['qmprogram'] = 'qe'
        elif text == 'Phonopy - Crystal':
            self.settings['program']   = 'phonopy'
            self.settings['qmprogram'] = 'crystal'
        elif text == 'Quantum Espresso':
            self.settings['program']   = 'qe'
        debugger.print('Program is now  ', self.settings['program'])
        debugger.print('QMProgram is now', self.settings['qmprogram'])
