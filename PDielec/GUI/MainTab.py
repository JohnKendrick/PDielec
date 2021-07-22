import os.path
from PyQt5.QtWidgets                import  QWidget
from PyQt5.QtWidgets                import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets                import  QFileDialog, QPushButton, QCheckBox
from PyQt5.QtWidgets                import  QFormLayout
from PyQt5.QtWidgets                import  QHBoxLayout, QVBoxLayout, QMessageBox
from PyQt5.QtCore                   import  Qt, QCoreApplication
from PDielec.Utilities              import  pdgui_get_reader, Debug
from PDielec.GUI.SpreadSheetManager import  SpreadSheetManager
import numpy as np

class MainTab(QWidget):

    def __init__(self, parent, program, filename, excelfile, debug=False):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug, 'MainTab')
        debugger.print('Initialising ',program, filename, excelfile)
        self.debug = debug
        self.settings = {}
        if program != '':
            self.settings['Program'] = program.lower()
        else:
            self.settings['Program'] = 'castep'
        self.settings['Output file name'] = filename
        self.directory = os.path.dirname(self.settings['Output file name'])
        self.settings['Excel file name'] = excelfile
        self.settings['Script file name'] = ''
        self.settings['QM program'] = ''
        self.settings['Hessian symmetrisation'] = 'symm'
        self.notebook = parent
        self.reader = None
        self.frequencies_cm1 = None
        self.refreshRequired = True
        self.calculationRequired = True
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
        self.program_cb.addItem('Experiment')
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
        # Store results
        #
        self.resultsfile_le = QLineEdit(self)
        self.resultsfile_le.setToolTip('Provide the name of an .xlsx file if results are to be stored in a spreadsheet.\nIf specified, the file will be written when the program exits.\nThe default directory is that of the QM/MM output file which has been read in.')
        self.resultsfile_le.setText(self.settings['Excel file name'])
        self.resultsfile_le.returnPressed.connect(self.on_excel_button_clicked)
        self.resultsfile_le.textChanged.connect(self.on_resultsfile_le_changed)
        label = QLabel('Excel spread sheet')
        label.setToolTip('Provide the name of an .xlsx file if results are to be stored in a spreadsheet.\nIf specified the file will be written when the program exits.\nThe default directory is that of the QM/MM output file which has been read in.')
        excel_button = QPushButton('Save results  ')
        excel_button.setToolTip('Save the results of calculation to the excel spreadsheet specified')
        excel_button.clicked.connect(self.on_excel_button_clicked)
        hbox = QHBoxLayout()
        hbox.addWidget(self.resultsfile_le)
        hbox.addWidget(excel_button)
        form.addRow(label, hbox)
        #
        # Store a script
        #
        self.scriptsfile_le = QLineEdit(self)
        self.scriptsfile_le.setToolTip('Provide the name of a python script file to save the program settings to when the \"Save settings\" button is pressed.\nThe directory where the file is saved is the same as the directory containing the QM/MM output file read in.')
        self.scriptsfile_le.setText(self.settings['Script file name'])
        self.scriptsfile_le.returnPressed.connect(self.on_script_button_clicked)
        self.scriptsfile_le.textChanged.connect(self.on_scriptsfile_le_changed)
        label = QLabel('Script filename')
        label.setToolTip('Provide the name of a python script file to save the program settings to when the \"Save settings\" button is pressed.\nThe directory where the file is saved is the same as the directory containing the QM/MM output file read in.')
        script_button = QPushButton('Save settings')
        script_button.setToolTip('Save the setttings of calculation to a python file.\nThe directory where the file is saved is the same as the directory containing the QM/MM output file read in.')
        script_button.clicked.connect(self.on_script_button_clicked)
        hbox = QHBoxLayout()
        hbox.addWidget(self.scriptsfile_le)
        hbox.addWidget(script_button)
        form.addRow(label, hbox)
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
        # output window for unit cell
        self.cell_window_l = QLabel('Unit-cell (Angstrom) from '+self.settings['Output file name'], self)
        vbox.addWidget(self.cell_window_l)
        self.cell_window_w = QListWidget(self)
        fm = self.cell_window_w.fontMetrics()
        h = fm.ascent() + fm.descent()
        self.cell_window_w.setMaximumHeight(6*h)
        vbox.addWidget(self.cell_window_w)
        # output window for frequencies
        self.frequencies_window_l = QLabel('Frequencies from '+self.settings['Output file name'], self)
        vbox.addWidget(self.frequencies_window_l)
        self.frequencies_window = QListWidget(self)
        vbox.addWidget(self.frequencies_window)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        self.open_excel_spreadsheet()
        # If the filename was given then force it to be read and processed
        if filename != '':
            debugger.print('Reading output file in maintab initialisation')
            self.on_calculation_button_clicked()
        QCoreApplication.processEvents()

    def on_script_button_clicked(self):
        debugger.print('on_script_button clicked')
        self.directory = os.path.dirname(self.settings['Output file name'])
        debugger.print('on_script_button clicked, directory=',self.directory)
        self.notebook.print_settings(filename=os.path.join(self.directory,self.settings['Script file name']))

    def open_excel_spreadsheet(self):
        debugger.print('open_spreadsheet clicked')
        if len(self.settings['Excel file name']) > 5 and self.settings['Excel file name'][-5:] == '.xlsx':
            self.directory = os.path.dirname(self.settings['Output file name'])
            # open the file name with the directory of the output file name
            self.openSpreadSheet(os.path.join(self.directory,self.settings['Excel file name']))
        elif len(self.settings['Excel file name']) > 1 and self.settings['Excel file name'][-5:] != '.xlsx':
            # The file isn't valid so tell the user there is a problem
            QMessageBox.about(self,'Spreadsheet name','File name of spreadsheet must end in  .xlsx')
        return

    def on_excel_button_clicked(self):
        debugger.print('on_excel_button clicked')
        self.open_excel_spreadsheet()
        if self.notebook.spreadsheet is not None:
            self.notebook.writeSpreadsheet()
        return

    def on_calculation_button_clicked(self):
        debugger.print('Calculation button clicked')
        #
        # Read the output file
        #
        self.read_output_file()
        self.calculationRequired = False

    def writeSpreadsheet(self):
        sp = self.notebook.spreadsheet
        if sp is None:
            return
        debugger.print('Reading output file ', self.settings['Output file name'])
        sp.selectWorkSheet('Main')
        sp.delete()
        debugger.print('writeSpreadsheet',self.settings)
        sp.writeNextRow( ['Main Tab Settings'], col=1 )
        for item in sorted(self.settings):
            sp.writeNextRow([item,self.settings[item]], col=1, check=1)
        sp.writeNextRow( [''], col=1 )
        if self.frequencies_cm1 is not None:
            sp.writeNextRow( ['Frequencies (cm1) as read from the output file'], col=1 )
            for ifreq in enumerate(self.frequencies_cm1):
                sp.writeNextRow(ifreq, col=1, check=1)

    def read_output_file(self):
        if self.settings['Output file name'] == '':
            return
        if not os.path.isfile(self.settings['Output file name']):
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+self.settings['Output file name'])
            return
        debugger.print('Read output file - clear list widgets')
        self.cell_window_w.clear()
        self.cell_window_l.setText('Unit-cell (Angstrom) from '+self.settings['Output file name'])
        self.frequencies_window.clear()
        self.frequencies_window_l.setText('Frequencies from '+self.settings['Output file name'])
        self.reader = pdgui_get_reader(self.settings['Program'],[ self.settings['Output file name'] ], self.settings['QM program'] )
        if self.reader is None:
            print('Error in reading files - program  is ',self.settings['Program'])
            print('Error in reading files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','A reader has not be created for this filename: '+self.settings['Output file name'])
            return
        #switch on debugging in the reader
        #self.reader.debug = self.debug
        self.reader.hessian_symmetrisation = self.settings['Hessian symmetrisation']
        if self.debug:
            self.reader.read_output()
        else:
            try:
                self.reader.read_output()
            except:
                print('Error in reading output files - program  is ',self.settings['Program'])
                print('Error in reading output files - filename is ',self.settings['Output file name'])
                print('Need to choose the file and program properly')
                QMessageBox.about(self,'Processing output file','Error on reading the output file using read_output(): '+self.settings['Output file name'])
                return
            # end try
        # end if debug
        if len(self.reader.unit_cells) == 0:
            print('Error in reading output files - program  is ',self.settings['Program'])
            print('Error in reading output files - filename is ',self.settings['Output file name'])
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The output file has no unit cells in it: '+self.settings['Output file name'])
            return
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
        self.directory = os.path.dirname(self.settings['Output file name'])
        if self.debug:
            self.reader.print_info()
        cell = self.reader.unit_cells[-1]
        a = cell.lattice[0]
        b = cell.lattice[1]
        c = cell.lattice[2]
        self.cell_window_w.addItem('                                  ')
        self.cell_window_w.addItem(f'a = {a[0]: 4.5f}   {a[1]:4.5f}   {a[2]: 4.5f}')
        self.cell_window_w.addItem(f'b = {b[0]: 4.5f}   {b[1]:4.5f}   {b[2]: 4.5f}')
        self.cell_window_w.addItem(f'c = {c[0]: 4.5f}   {c[1]:4.5f}   {c[2]: 4.5f}')
        self.cell_window_w.addItem('                                  ')
        self.frequencies_cm1 = np.sort(self.reader.frequencies)
        for f in self.frequencies_cm1:
            self.frequencies_window.addItem('{0:.3f}'.format(f))
        # tell the settings tab to update the widgets that depend on the contents of the reader
        debugger.print('processing a return')
        if self.notebook.settingsTab is not None:
            debugger.print('about to refresh settings')
            # There is a subtle problem with the sigmas, they need to be reset on reading a new file
            self.notebook.settingsTab.sigmas_cm1 = []
            self.notebook.settingsTab.refresh(force=True)
        # Update any scenarios
        if self.notebook.scenarios is not None:
            debugger.print('about to refresh scenarios')
            debugger.print('notebook has {} scenarios'.format(len(self.notebook.scenarios)))
            for tab in self.notebook.scenarios:
                tab.refresh(force=True)
        else:
            debugger.print('notebook has no scenarios yet')
        # Update the plotting tab
        if self.notebook.plottingTab is not None:
            debugger.print('about to refresh plottingtab')
            self.notebook.plottingTab.refresh(force=True)
        else:
            debugger.print('notebook has no plotting tab yet')
        # Update the analysis tab
        if self.notebook.analysisTab is not None:
            debugger.print('about to refresh analysisTab')
            self.notebook.analysisTab.refresh(force=True)
        else:
            debugger.print('notebook has no analysis tab yet')
        # Update the viewer tab
        if self.notebook.viewerTab is not None:
            debugger.print('about to refresh viewerTab')
            self.notebook.viewerTab.refresh(force=True)
        else:
            debugger.print('notebook has no viewer tab yet')
        # Update the fitter tab
        if self.notebook.fitterTab is not None:
            # There is a subtle problem with the modes that are to be fitted, they need resetting
            self.notebook.fitterTab.modes_fitted = []
            debugger.print('about to refresh fitterTab')
            self.notebook.fitterTab.refresh(force=True)
        else:
            debugger.print('notebook has no fitter tab yet')

    def on_hessian_symmetry_changed(self):
        debugger.print('on hessian symmetry changed')
        if self.hessian_symmetry_cb.isChecked():
            self.settings['Hessian symmetrisation'] = 'crystal'
        else:
            self.settings['Hessian symmetrisation'] = 'symm'
        self.calculationRequired = True

    def on_scriptsfile_le_changed(self, text):
        debugger.print('on scriptsfile changed', text)
        text = self.scriptsfile_le.text()
        self.settings['Script file name'] = text
        self.calculationRequired = True

    def on_resultsfile_le_changed(self, text):
        debugger.print('on resultsfile changed', text)
        text = self.resultsfile_le.text()
        self.settings['Excel file name'] = text
        self.calculationRequired = True

    def openSpreadSheet(self,text):
        if self.notebook.spreadsheet is not None:
            self.notebook.spreadsheet.close()
        if text[-5:] == '.xlsx':
            self.notebook.spreadsheet = SpreadSheetManager(text)
            self.calculationRequired = True
        else:
           print('spreadsheet name not valid', text)

    def on_file_le_return(self):
        debugger.print('on file return ', self.file_le.text())
        # Does the file exist?
        self.settings['Output file name'] = self.file_le.text()
        # Open a file chooser
        #options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getOpenFileName(self,'Open MM/QM Output file','','Castep (*.castep);;Abinit (*.out);;Gulp (*.gout);;VASP (OUTCAR*);; QE (*.dynG);; Crystal 14 (*.out);; Phonopy (OUTCAR*);; Experiment (*.exp);; All Files (*)')
        if filename != '':
            self.settings['Output file name'] = filename
            self.file_le.setText(self.settings['Output file name'])
            debugger.print('new file name', self.settings['Output file name'])
            self.directory = os.path.dirname(self.settings['Output file name'])
            self.calculationRequired = True
        return

    def on_file_le_changed(self, text):
        debugger.print('on file changed', text)
        self.settings['Output file name'] = text
        self.calculationRequired = True

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
        self.calculationRequired = True

    def refresh(self,force=False):
        debugger.print('refresh', force)
        if not self.refreshRequired and not force:
            return
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
        # JK self.on_calculation_button_clicked()
        self.open_excel_spreadsheet()
        self.refreshRequired = False
        #
        # UnBlock signals
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        return
