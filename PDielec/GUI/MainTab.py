import sys
import os.path
from PyQt5.QtWidgets  import  QWidget
from PyQt5.QtWidgets  import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QFileDialog, QPushButton, QCheckBox
from PyQt5.QtWidgets  import  QFormLayout, QApplication
from PyQt5.QtWidgets  import  QVBoxLayout, QMessageBox
from PyQt5.QtCore     import  Qt, QCoreApplication
from PDielec.Utilities import  get_reader, Debug
from PDielec.GUI.SpreadSheetManager import SpreadSheetManager
import numpy as np

class MainTab(QWidget):

    def __init__(self, parent, program, filename, excelfile, debug=False):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug, 'MainTab')
        debugger.print('Initialising ',program, filename,excelfile)
        self.debug = debug
        self.settings = {}
        if program != '':
            self.settings['Program'] = program.lower()
        else:
            self.settings['Program'] = 'castep'
        self.settings['Output file name'] = filename
        self.directory = os.path.dirname(self.settings['Output file name'])
        self.settings['Excel file name'] = excelfile
        self.settings['QM program'] = ''
        self.settings['Hessian symmetrisation'] = 'symm'
        self.notebook = parent
        self.notebook.plottingCalculationRequired = True
        self.notebook.analysisCalculationRequired = True
        self.notebook.visualerCalculationRequired = True
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
        # The method for symmetrising the hessian
        #
        self.hessian_symmetry_cb = QCheckBox(self)
        self.hessian_symmetry_cb.setToolTip('The Crystal program uses a different method for symmetrising the hessian.  Check this flag if you want to use the same method as used by Crystal14')
        self.hessian_symmetry_cb.setText('')
        self.hessian_symmetry_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings['Program'] == 'crystal':
            self.settings['Hessian symmetrisation'] = 'crystal'
            self.hessian_symmetry_cb.setCheckState(Qt.Checked)
            self.hessian_symmetry_cb.setEnabled(True)
        else:
            self.hessian_symmetry_cb.setCheckState(Qt.Unchecked)
            self.hessian_symmetry_cb.setEnabled(False)
        self.hessian_symmetry_cb.stateChanged.connect(self.on_hessian_symmetry_changed)
        label = QLabel('Symmetrise the hessian for Crystal14')
        label.setToolTip('The Crystal program uses a different method for symmetrising the hessian.  Check this flag if you want to use the same method as used by Crystal14')
        form.addRow(label, self.hessian_symmetry_cb)
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
        #
        # Calculation button
        #
        self.calculation_button = QPushButton('Read the output file and start the calculation')
        self.calculation_button.setToolTip('Pressing this button initiates reading of the output file and processing of all the options')
        self.calculation_button.clicked.connect(self.on_calculation_button_clicked)
        #label = QLabel('Read and calculate')
        #label.setToolTip('Pressing this button initiates reading of the output file and processing of all the options')
        form.addRow(self.calculation_button)
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
            debugger.print('Reading output file in maintab initialisation')
            self.on_calculation_button_clicked()
        QCoreApplication.processEvents()

    def on_calculation_button_clicked(self):
        debugger.print('Calculation button clicked')
        if len(self.settings['Excel file name']) > 5 and self.settings['Excel file name'][-5:] == '.xlsx':
            self.directory = os.path.dirname(self.settings['Output file name'])
            # open the file name with the directory of the output file name
            self.openSpreadSheet(os.path.join(self.directory,self.settings['Excel file name']))
        elif len(self.settings['Excel file name']) > 1 and self.settings['Excel file name'][-5:] != '.xlsx':
            # The file isn't valid so tell the user there is a problem
            QMessageBox.about(self,'Spreadsheet name','File name of spreadsheet must end in  .xlsx')
            return
        #
        # Read the output file
        #
        self.read_output_file()
        self.dirty = False
        # If there is a spread sheet write the spread sheet
        # if self.notebook.spreadsheet is not None:
        #     self.write_spreadsheet()

    def write_spreadsheet(self):
        sp = self.notebook.spreadsheet
        if sp is None:
            return
        debugger.print('Reading output file ', self.settings['Output file name'])
        sp.selectWorkSheet('Main')
        sp.delete()
        debugger.print('write_spreadsheet',self.settings)
        sp.writeNextRow( ['Main Tab Settings'], col=1 )
        for item in sorted(self.settings):
            sp.writeNextRow([item,self.settings[item]], col=1, check=1)
        sp.writeNextRow( [''], col=1 )
        sp.writeNextRow( ['Frequencies (cm1) as read from the output file'], col=1 )
        for ifreq in enumerate(self.frequencies_cm1):
            sp.writeNextRow(ifreq, col=1, check=1)

    def read_output_file(self):
        if self.settings['Output file name'] == '':
            return
        if not os.path.isfile(self.settings['Output file name']):
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['Output file name'])
            return
        debugger.print('Read output file - clear list widget')
        self.listw.clear()
        self.listw_l.setText('Frequencies from '+self.settings['Output file name'])
        self.reader = get_reader(self.settings['Program'],[ self.settings['Output file name'] ], self.settings['QM program'] )
        if self.reader is None:
            print('Error in reading files - program  is ',self.settings['Program'])
            print('Error in reading files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','A reader has not be created for this filename: '+self.settings['Output file name'])
            return
        #switch on debugging in the reader
        #self.reader.debug = self.debug
        self.reader.hessian_symmetrisation = self.settings['Hessian symmetrisation']
        try:
            self.reader.read_output()
        except:
            print('Error in reading output files - program  is ',self.settings['Program'])
            print('Error in reading output files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','Error on reading the output file using read_output(): '+self.settings['Output file name'])
            return
        if len(self.reader.unit_cells) == 0:
            print('Error in reading output files - program  is ',self.settings['Program'])
            print('Error in reading output files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The output file has no unit cells in it: '+self.settings['Output file name'])
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Update the checkbox
        if self.settings['Hessian symmetrisation'] == 'crystal':
            self.hessian_symmetry_cb.setCheckState(Qt.Checked)
        else:
            self.hessian_symmetry_cb.setCheckState(Qt.Unchecked)
        if self.settings['Program'] == 'crystal':
            self.hessian_symmetry_cb.setEnabled(True)
        else:
            self.hessian_symmetry_cb.setEnabled(False)
        # tell the notebook that we have read the info and we have a reader
        self.notebook.reader = self.reader
        self.notebook.plottingCalculationRequired = True
        self.notebook.analysisCalculationRequired = True
        self.notebook.visualerCalculationRequired = True
        self.directory = os.path.dirname(self.settings['Output file name'])
        if self.debug:
            self.reader.print_info()
        self.frequencies_cm1 = np.sort(self.reader.frequencies)
        for f in self.frequencies_cm1:
            self.listw.addItem('{0:.3f}'.format(f))
        # tell the settings tab to update the widgets that depend on the contents of the reader
        debugger.print('processing a return')
        if hasattr(self.notebook, 'settingsTab'):
            debugger.print('about to refresh settings')
            self.notebook.settingsTab.refresh(True)
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
        # Update the fitter tab
        if hasattr(self.notebook, 'fitterTab'):
            debugger.print('about to refresh fitterTab')
            self.notebook.fitterTab.refresh()
        else:
            debugger.print('notebook has no fitter tab yet')
        QApplication.restoreOverrideCursor()

    def on_hessian_symmetry_changed(self):
        debugger.print('on hessian symmetry changed')
        if self.hessian_symmetry_cb.isChecked():
            self.settings['Hessian symmetrisation'] = 'crystal'
        else:
            self.settings['Hessian symmetrisation'] = 'symm'
        self.dirty = True

    def on_resultsfile_le_return(self):
        debugger.print('on resultsfile return ', self.resultsfile_le.text())

    def on_resultsfile_le_changed(self, text):
        debugger.print('on resultsfile changed', text)
        text = self.resultsfile_le.text()
        self.settings['Excel file name'] = text
        self.dirty = True

    def openSpreadSheet(self,text):
        if self.notebook.spreadsheet is not None:
            self.notebook.spreadsheet.close()
        if text[-5:] == '.xlsx':
            self.notebook.spreadsheet = SpreadSheetManager(text)
            self.dirty = True
        else:
           print('spreadsheet name not valid', text)
           exit()

    def on_file_le_return(self):
        debugger.print('on file return ', self.file_le.text())
        # Does the file exist?
        self.settings['Output file name'] = self.file_le.text()
        # Open a file chooser
        #options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getOpenFileName(self,'Open MM/QM Output file','','Castep (*.castep);;Abinit (*.out);;Gulp (*.gout);;VASP (OUTCAR*);; QE (*.dynG);; Crystal 14 (*.out);; Phonopy (OUTCAR*);; All Files (*)')
        if filename != '':
            self.settings['Output file name'] = filename
            self.file_le.setText(self.settings['Output file name'])
            debugger.print('new file name', self.settings['Output file name'])
            self.directory = os.path.dirname(self.settings['Output file name'])
            self.dirty = True
        return

    def on_file_le_changed(self, text):
        debugger.print('on file changed', text)
        self.settings['Output file name'] = text
        self.dirty = True

    def on_program_cb_activated(self,index):
        debugger.print('on program combobox activated', index)
        debugger.print('on program combobox activated', self.program_cb.currentText())
        text = self.program_cb.currentText()
        self.settings['Program'] = text.lower()
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
        if self.settings['Program'] == 'crystal':
            self.hessian_symmetry_cb.setEnabled(True)
            self.hessian_symmetry_cb.setCheckState(Qt.Checked)
        else:
            self.hessian_symmetry_cb.setEnabled(True)
            self.hessian_symmetry_cb.setCheckState(Qt.Unchecked)
            self.hessian_symmetry_cb.setEnabled(False)
        self.dirty = True

    def refresh(self,force=False):
        debugger.print('refresh', force)
        if not self.dirty and not force:
            return
        self.notebook.newPlottingCalculationRequired = True
        self.notebook.newAnalysisCalculationRequired = True
        self.notebook.newVisualerCalculationRequired = True
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
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
        if self.settings['Hessian symmetrisation'] == 'crystal':
            self.hessian_symmetry_cb.setCheckState(Qt.Checked)
        else:
            self.hessian_symmetry_cb.setCheckState(Qt.Unchecked)
        if self.settings['Program'] == 'crystal':
            self.hessian_symmetry_cb.setEnabled(True)
        else:
            self.hessian_symmetry_cb.setEnabled(False)
        self.file_le.setText(self.settings['Output file name'])
        self.resultsfile_le.setText(self.settings['Excel file name'])
        self.on_calculation_button_clicked()
        self.dirty = False
        #
        # UnBlock signals
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        return
