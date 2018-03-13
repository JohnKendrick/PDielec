import sys
import os.path
from PyQt5.QtWidgets  import  QWidget
from PyQt5.QtWidgets  import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QFileDialog, QPushButton
from PyQt5.QtWidgets  import  QFormLayout
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QMessageBox
from PyQt5.QtCore     import  Qt, QCoreApplication
from Python.Utilities import  get_reader, Debug
from Python.GUI.SpreadSheetManager import SpreadSheetManager
import numpy as np
 
class MainTab(QWidget):        
 
    def __init__(self, parent, program, filename, excelfile, debug=False):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug, 'MainTab')
        self.debug = debug
        self.settings = {}
        self.settings['Program'] = 'Castep'
        if program != '':
            self.settings['Program'] = program.lower()
        self.settings['Output file name'] = filename
        self.settings['Excel file name'] = excelfile
        self.settings['QM program'] = ''
        self.notebook = parent
        self.notebook.plottingCalculationRequired = True
        self.notebook.analysisCalculationRequired = True
        self.reader = None
        self.frequencies_cm1 = None
        self.dirty = True
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
        self.program_cb.addItem('Phonopy - VASP')
        self.program_cb.addItem('Quantum Espresso')
        self.program_cb.addItem('Gulp')
#        self.program_cb.addItem('Phonopy - QE') 
#        self.program_cb.addItem('Phonopy - Crystal')
#        self.program_cb.addItem('Experiment')
        prtext = self.settings['Program'].capitalize()
        qmtext = self.settings['QM program'].capitalize()
        if prtext == 'Qe':
            prtext = 'Quantum Espresso'
        elif prtext == 'Vasp':
            prtext = 'VASP'
            self.settings['QM program'] = 'vasp'
        elif prtext == 'Phonopy':
            if qmtext != '':
                prtext = prtext + ' - ' + qmtext
        index = self.program_cb.findText(prtext, Qt.MatchFixedString)
        if index >=0:
            self.program_cb.setCurrentIndex(index)
        elif self.settings['Program'] == 'phonopy':
            index = 4
            self.settings['QM program'] = 'vasp'
            self.program_cb.setCurrentIndex(index)
        else:
            self.settings['Program'] = 'castep'
            index = 1
            self.program_cb.setCurrentIndex(index)
        self.program_cb.activated.connect(self.on_program_cb_activated)
        label = QLabel('QM/MM Program ')
        label.setToolTip('Choose QM/MM program')
        form.addRow(label, self.program_cb)
        #
        # The file selector
        #
        self.file_le = QLineEdit(self)
        self.file_le.setToolTip('Choose output file for analysis (press return for a file chooser)')
        self.file_le.setText(self.settings['Output file name'])
        self.file_le.returnPressed.connect(self.on_file_le_return)
        self.file_le.textChanged.connect(self.on_file_le_changed)
        label = QLabel('Analyse this output file')
        label.setToolTip('Choose output file for analysis (press return for a file chooser)')
        form.addRow(label, self.file_le)
        #
        # The store results
        #
        self.resultsfile_le = QLineEdit(self)
        self.resultsfile_le.setToolTip('Provide the name of an .xlsx file if results are to be stored in a spreadsheet')
        self.resultsfile_le.setText(self.settings['Excel file name'])
        self.resultsfile_le.returnPressed.connect(self.on_resultsfile_le_return)
        self.resultsfile_le.textChanged.connect(self.on_resultsfile_le_changed)
        label = QLabel('Excel spread sheet')
        label.setToolTip('Provide the name of an .xlsx file if results are to be stored in a spreadsheet')
        form.addRow(label, self.resultsfile_le)
        # add form layout 
        vbox.addLayout(form)
        # output window
        self.listw_l = QLabel('Frequencies from '+self.settings['Output file name'], self)
        vbox.addWidget(self.listw_l)
        self.listw = QListWidget(self)
        vbox.addWidget(self.listw)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        # If the filename was given then force it to be read and processed
        if filename != '':
            self.read_output_file()
        QCoreApplication.processEvents()
        # If there is a spreadsheet the write it
        if self.settings['Excel file name'] != '':
            self.openSpreadSheet(self.settings['Excel file name'])
        if self.notebook.spreadsheet is not None:
            self.write_spreadsheet()
        QCoreApplication.processEvents()

    def write_spreadsheet(self):
        sp = self.notebook.spreadsheet
        if sp is None:
            return
        debugger.print('Reading output file ', self.settings['Output file name'])
        sp.selectWorkSheet('Main')
        debugger.print('write_spreadsheet',self.settings)
        sp.writeNextRow( ['Main Tab Settings'], col=1 )
        for item in sorted(self.settings):
            sp.writeNextRow([item,self.settings[item]], col=1, check=1)
        sp.writeNextRow( [''], col=1 )
        sp.writeNextRow( ['Frequencies (cm1) as read from the output file'], col=1 )
        for ifreq in enumerate(self.frequencies_cm1):
            sp.writeNextRow(ifreq, col=1, check=1)
            
    def read_output_file(self):
        if not os.path.isfile(self.settings['Output file name']):
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['Output file name'])
            return
        self.listw_l.clear()
        self.listw_l.setText('Frequencies from '+self.settings['Output file name'])
        self.reader = get_reader(self.settings['Program'],[ self.settings['Output file name'] ], self.settings['QM program'] )
        if self.reader is None:
            print('Error in reading files - program  is ',self.settings['Program'])
            print('Error in reading files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['Output file name'])
            return
        #switch on debugging in the reader
        self.reader.debug = self.debug
        try:
            self.reader.read_output()
        except:
            print('Error in reading output files - program  is ',self.settings['Program'])
            print('Error in reading output files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['Output file name'])
            return
        if len(self.reader.unit_cells) == 0:
            print('Error in reading output files - program  is ',self.settings['Program'])
            print('Error in reading output files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['Output file name'])
            return
        # tell the notebook that we have read the info and we have a reader
        self.notebook.reader = self.reader
        self.notebook.plottingCalculationRequired = True
        self.notebook.analysisCalculationRequired = True
        if self.debug:
            self.reader.print_info()
        self.frequencies_cm1 = np.sort(self.reader.frequencies)
        for f in self.frequencies_cm1:
            self.listw.addItem('{0:.3f}'.format(f))
        # tell the settings tab to update the widgets that depend on the contents of the reader
        debugger.print('processing a return')
        if hasattr(self.notebook, 'settingsTab'):
            debugger.print('about to refresh settings')
            self.notebook.settingsTab.refresh()
        # Update any scenarios
        if hasattr(self.notebook, 'scenarios'):
            debugger.print('about to refresh scenarios')
            debugger.print('notebook has {} scenarios'.format(len(self.notebook.scenarios)))
            for tab in self.notebook.scenarios:
                tab.refresh()
        else:
            debugger.print('notebook has no scenarios yet')
        # Update the plotting tab
        if hasattr(self.notebook, 'plottingTab'):
            debugger.print('about to refresh plottingtab')
            self.notebook.plottingTab.refresh()
        else:
            debugger.print('notebook has no plotting tab yet')
        # Update the analysis tab
        if hasattr(self.notebook, 'analysisTab'):
            debugger.print('about to refresh analysisTab')
            self.notebook.analysisTab.refresh()
        else:
            debugger.print('notebook has no analysis tab yet')
        # Update the viewer tab
        if hasattr(self.notebook, 'viewerTab'):
            debugger.print('about to refresh viewerTab')
            self.notebook.viewerTab.refresh()
        else:
            debugger.print('notebook has no viewer tab yet')

    def on_resultsfile_le_return(self):
        debugger.print('on resultsfile return ', self.resultsfile_le.text())
        # Does the file end with .xlsx?
        text = self.resultsfile_le.text()
        self.settings['Excel file name'] = text
        if text[-5:] == '.xlsx':
            self.openSpreadSheet(text)
            # The file is valid so just return
        else:
            # The file doesn't exist so open a file chooser
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            text, _ = QFileDialog.getOpenFileName(self,'QFileDialog.getOpenFileName()', '','All Files (*.xlsx)', options=options)
            if text != '':
                debugger.print('new resultsfile name', text)
                self.resultsfile_le.setText(text)
                self.openSpreadSheet(text)
        if self.notebook.spreadsheet is not None:
            if os.path.isfile(self.settings['Output file name']):
                debugger.print('on resultsfile about to press the button ')
                # The file exists to treat it as though the button has been pressed
                self.pushButton1Clicked()
        debugger.print('on resultsfile return ')
 
    def on_resultsfile_le_changed(self, text):
        debugger.print('on resultsfile changed', text)
        self.openSpreadSheet(text)
        self.dirty = True

    def openSpreadSheet(self,text):
        if self.notebook.spreadsheet is not None:
            self.notebook.spreadsheet.close()
        if text[-5:] == '.xlsx':
            self.notebook.spreadsheet = SpreadSheetManager(text)
            self.settings['Excel file name'] = text
        else:
            print('spreadsheet name not valid', text)

    def on_file_le_return(self):
        debugger.print('on file return ', self.file_le.text())
        # Does the file exist?
        self.settings['Output file name'] = self.file_le.text()
        if os.path.isfile(self.settings['Output file name']):
            # The file exists to treat it as though the button has been pressed
            self.read_output_file()
            return
        # The file doesn't exist so open a file chooser
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.settings['Output file name'], _ = QFileDialog.getOpenFileName(self,'QFileDialog.getOpenFileName()', '','All Files (*)', options=options)
        if self.settings['Output file name']:
            debugger.print('new file name', self.settings['Output file name'])
        self.file_le.setText(self.settings['Output file name'])
        if os.path.isfile(self.settings['Output file name']):
            #self.settings['Output file name'] = os.path.relpath(self.settings['Output file name'])
            # The file exists to treat it as though the button has been pressed
            self.read_output_file()
            return
 
    def on_file_le_changed(self, text):
        debugger.print('on file changed', text)
        self.settings['Output file name'] = text
        self.dirty = True
 
    def on_program_cb_activated(self,index):
        debugger.print('on program combobox activated', index)
        debugger.print('on program combobox activated', self.program_cb.currentText())
        text = self.program_cb.currentText()
        self.settings['Program'] = text
        self.settings['QM program'] = ''
        if text == 'Phonopy - VASP':
            self.settings['Program']   = 'phonopy'
            self.settings['QM program'] = 'vasp'
        elif text == 'Phonopy - QE':
            self.settings['Program']   = 'phonopy'
            self.settings['QM program'] = 'qe'
        elif text == 'Phonopy - Crystal':
            self.settings['Program']   = 'phonopy'
            self.settings['QM program'] = 'crystal'
        elif text == 'Quantum Espresso':
            self.settings['Program']   = 'qe'
        debugger.print('Program is now  ', self.settings['Program'])
        debugger.print('QM program is now', self.settings['QM program'])
        self.dirty = True

    def refresh(self,force=False):
        debugger.print('refresh', force)
        if not self.dirty and not force:
            return
        self.notebook.newPlottingCalculationRequired = True
        self.notebook.newAnalysisCalculationRequired = True
        prtext = self.settings['Program'].capitalize()
        qmtext = self.settings['QM program'].capitalize()
        if prtext == 'Qe':
            prtext = 'Quantum Espresso'
        elif prtext == 'Vasp':
            prtext = 'VASP'
            self.settings['QM program'] = 'vasp'
        elif prtext == 'Phonopy':
            if qmtext != '':
                prtext = prtext + ' - ' + qmtext
        index = self.program_cb.findText(prtext, Qt.MatchFixedString)
        if index >=0:
            self.program_cb.setCurrentIndex(index)
        elif self.settings['Program'] == 'phonopy':
            index = 4
            self.settings['QM program'] = 'vasp'
            self.program_cb.setCurrentIndex(index)
        else:
            self.settings['Program'] = 'castep'
            index = 1
            self.program_cb.setCurrentIndex(index)
        self.file_le.setText(self.settings['Output file name'])
        self.resultsfile_le.setText(self.settings['Excel file name'])
        self.dirty = False
        return
