import os.path
from PyQt5.QtWidgets                import  QWidget
from PyQt5.QtWidgets                import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets                import  QFileDialog, QPushButton, QCheckBox
from PyQt5.QtWidgets                import  QFormLayout
from PyQt5.QtWidgets                import  QHBoxLayout, QVBoxLayout, QMessageBox
from PyQt5.QtCore                   import  Qt, QCoreApplication, QSize
from PDielec.Utilities              import  find_program_from_name, pdgui_get_reader, Debug
from PDielec.GUI.SpreadSheetManager import  SpreadSheetManager
import numpy as np

class MainTab(QWidget):

    def __init__(self, parent, program, filename, excelfile, debug=False):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug, 'MainTab')
        debugger.print('Start:: initialising ',program, filename, excelfile)
        self.debug = debug
        self.settings = {}
        if program != '':
            self.settings['Program'] = program.lower()
        else:
            self.settings['Program'] = 'castep'
        self.directory = os.path.dirname(filename)
        self.settings['Output file name'] = os.path.basename(filename)
        self.settings['Excel file name'] = excelfile
        self.settings['Script file name'] = ''
        self.settings['QM program'] = ''
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
        self.program_cb.addItem('Experiment')
        self.program_cb.addItem('Gulp')
        self.program_cb.addItem('Phonopy - VASP')
        self.program_cb.addItem('Quantum Espresso')
        self.program_cb.addItem('Vasp')
        self.program_cb.addItem('PDGui')
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
        self.file_le.setToolTip('Enter output file or script for processing (press return)')
        self.file_le.setText(self.settings['Output file name'])
        self.file_le.returnPressed.connect(self.on_file_le_return)
        self.file_le.textChanged.connect(self.on_file_le_changed)
        label = QLabel('Analyse this output file')
        label.setToolTip('Enter output file for analysis or choose a script to process (press return)')
        file_button = QPushButton(' File manager ')
        file_button.setToolTip('Open a file manager to choose a file for analysis or choose a script to process')
        file_button.clicked.connect(self.on_file_button_clicked)
        file_button.resize(file_button.sizeHint())
        file_button_size = file_button.size()
        hbox = QHBoxLayout()
        hbox.addWidget(self.file_le)
        hbox.addWidget(file_button)
        form.addRow(label, hbox)
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
        excel_button = QPushButton(' Save results ')
        excel_button.setToolTip('Save the results of calculation to the excel spreadsheet specified')
        excel_button.clicked.connect(self.on_excel_button_clicked)
        excel_button.resize(excel_button.sizeHint())
        excel_button_size = excel_button.size()
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
        script_button = QPushButton(' Save settings ')
        script_button.setToolTip('Save the setttings of calculation to a python file.\nThe directory where the file is saved is the same as the directory containing the QM/MM output file read in.')
        script_button.clicked.connect(self.on_script_button_clicked)
        script_button.resize(script_button.sizeHint())
        script_button_size = script_button.size()
        maximum_button_width = max(script_button_size.width(), excel_button_size.width(), file_button_size.width())
        maximum_button_height = max(script_button_size.height(), excel_button_size.height(), file_button_size.height())
        new_button_size = QSize(maximum_button_width, maximum_button_height)
        script_button.setFixedSize(new_button_size)
        excel_button.setFixedSize(new_button_size)
        file_button.setFixedSize(new_button_size)
        hbox = QHBoxLayout()
        hbox.addWidget(self.scriptsfile_le)
        hbox.addWidget(script_button)
        form.addRow(label, hbox)
        #jk #
        #jk # Calculation button
        #jk #
        #jk self.calculation_button = QPushButton('Read the output file and start the calculation')
        #jk self.calculation_button.setToolTip('Pressing this button initiates reading of the output file and processing of all the options')
        #jk self.calculation_button.clicked.connect(self.on_calculation_button_clicked)
        #jk #label = QLabel('Read and calculate')
        #jk #label.setToolTip('Pressing this button initiates reading of the output file and processing of all the options')
        #jk form.addRow(self.calculation_button)
        #jk # add form layout
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
        # If the filename was given then force it to be read and processed
        if filename != '':
            debugger.print('Reading output file in maintab initialisation')
            self.on_calculation_button_clicked()
        QCoreApplication.processEvents()
        debugger.print('Finished:: initialising ')
        return

    def on_script_button_clicked(self):
        debugger.print('Start:: on_script_button clicked')
        debugger.print('on_script_button clicked, directory=',self.directory)
        filename=os.path.join(self.directory,self.settings['Script file name'])
        if os.path.exists(filename):
            debugger.print('Script file already exists',self.directory)
            if self.notebook.overwriting:
                debugger.print('Overwriting existing script anyway',filename)
                self.notebook.print_settings(filename=filename)
            else:
                answer = QMessageBox.question(self,'','Script file already exists.  Continue?', QMessageBox.Yes | QMessageBox.No)
                if answer == QMessageBox.Yes:
                    debugger.print('Overwriting existing script',filename)
                    self.notebook.print_settings(filename=filename)
        else:
            debugger.print('Creating a new script',filename)
            self.notebook.print_settings(filename=filename)
        debugger.print('Finished:: on_script_button clicked')
        return

    def on_excel_button_clicked(self):
        debugger.print('Start:: on_excel_button clicked')
        self.notebook.writeSpreadsheet()
        debugger.print('Finished:: on_excel_button clicked')
        return

    def on_calculation_button_clicked(self):
        debugger.print('Start:: on_calculation_button_clicked')
        #
        # Read the output file
        #
        self.read_output_file()
        if self.notebook.settingsTab is not None:
            self.notebook.settingsTab.refresh(force=True)
        self.calculationRequired = False
        debugger.print('Finished:: on_calculation_button_clicked')

    def writeSpreadsheet(self):
        debugger.print('Start:: writeSpreadsheet')
        sp = self.notebook.spreadsheet
        if sp is None:
            debugger.print('Finished:: writeSpreadsheet sp is None')
            return
        sp.selectWorkSheet('Main')
        sp.delete()
        debugger.print('writeSpreadsheet',self.settings)
        sp.writeNextRow( ['Main Tab Settings'], col=1 )
        sp.writeNextRow( ['Directory',self.directory], col=1 )
        for item in sorted(self.settings):
            sp.writeNextRow([item,self.settings[item]], col=1, check=1)
        sp.writeNextRow( [''], col=1 )
        if self.frequencies_cm1 is not None:
            sp.writeNextRow( ['Frequencies (cm1) as read from the output file'], col=1 )
            for ifreq in enumerate(self.frequencies_cm1):
                sp.writeNextRow(ifreq, col=1, check=1)
        debugger.print('Finished:: writeSpreadsheet')
        return

    def read_output_file(self):
        debugger.print('Start:: read_output_file')
        if self.settings['Output file name'] == '':
            debugger.print('Finished:: read_output_file output file is blank')
            return
        filename = os.path.join(self.directory,self.settings['Output file name'])
        if not os.path.isfile(filename):
            QMessageBox.about(self,'Processing output file','The filename for the output file to be processed is not correct: '+filename)
            debugger.print('Finished:: read_output_file output file does not exist')
            return
        debugger.print('Read output file - clear list widgets')
        self.cell_window_w.clear()
        self.cell_window_l.setText('Unit-cell (Angstrom) from '+self.settings['Output file name'])
        self.frequencies_window.clear()
        self.frequencies_window_l.setText('Frequencies from '+self.settings['Output file name'])
        self.reader = pdgui_get_reader(self.settings['Program'],[ filename ], self.settings['QM program'] )
        if self.reader is None:
            print('Error in reading files - program  is ',self.settings['Program'])
            print('Error in reading files - filename is ',filename)
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','A reader has not be created for this filename: '+filename)
            debugger.print('Finished:: read_output_file reader is none')
            return
        #switch on debugging in the reader
        self.reader.debug = self.debug
        if self.debug:
            self.reader.read_output()
        else:
            try:
                self.reader.read_output()
            except:
                print('Error in reading output files - program  is ',self.settings['Program'])
                print('Error in reading output files - filename is ',filename)
                print('Need to choose the file and program properly')
                QMessageBox.about(self,'Processing output file','Error on reading the output file using read_output(): '+filename)
                debugger.print('Finished:: read_output_file error on reading')
                return
            # end try
        # end if debug
        if len(self.reader.unit_cells) == 0:
            print('Error in reading output files - program  is ',self.settings['Program'])
            print('Error in reading output files - filename is ',filename)
            print('Need to choose the file and program properly')
            QMessageBox.about(self,'Processing output file','The output file has no unit cells in it: '+filename)
            debugger.print('Finished:: read_output_file output file has no unit cell')
            return
        # tell the notebook that we have read the info and we have a reader
        self.notebook.reader = self.reader
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
        debugger.print('processing a return in reading the output file')
        # Update any scenarios
        if self.notebook.scenarios is not None:
            debugger.print('about to refresh scenarios')
            debugger.print('notebook has {} scenarios'.format(len(self.notebook.scenarios)))
            for tab in self.notebook.scenarios:
                tab.requestRefresh()
        else:
            debugger.print('notebook has no scenarios yet')
        # Update the plotting tab
        if self.notebook.plottingTab is not None:
            debugger.print('about to refresh plottingtab')
            self.notebook.plottingTab.requestRefresh()
        else:
            debugger.print('notebook has no plotting tab yet')
        # Update the analysis tab
        if self.notebook.analysisTab is not None:
            debugger.print('about to refresh analysisTab')
            self.notebook.analysisTab.requestRefresh()
        else:
            debugger.print('notebook has no analysis tab yet')
        # Update the viewer tab
        if self.notebook.viewerTab is not None:
            debugger.print('about to refresh viewerTab')
            self.notebook.viewerTab.requestRefresh()
        else:
            debugger.print('notebook has no viewer tab yet')
        # Update the fitter tab
        if self.notebook.fitterTab is not None:
            # There is a subtle problem with the modes that are to be fitted, they need resetting
            self.notebook.fitterTab.modes_fitted = []
            debugger.print('about to refresh fitterTab')
            self.notebook.fitterTab.requestRefresh()
        else:
            debugger.print('notebook has no fitter tab yet')
        debugger.print('Finished:: read_output_file')

    def on_scriptsfile_le_changed(self, text):
        debugger.print('Start:: on_scriptsfile_changed', text)
        text = self.scriptsfile_le.text()
        self.settings['Script file name'] = text
        debugger.print('Finished:: on_scriptsfile_changed', text)

    def on_resultsfile_le_changed(self, text):
        debugger.print('Start:: on_resultsfile_changed', text)
        text = self.resultsfile_le.text()
        self.settings['Excel file name'] = text
        debugger.print('Finished:: on_resultsfile_changed', text)

    def on_file_button_clicked(self):
        debugger.print('Start:: on_file_button_clicked ', self.file_le.text())
        # Open a file chooser
        #options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        if self.settings['Program'] == 'castep':
            selfilter = 'Castep (*.castep)'
        elif self.settings['Program'] == 'abinit':
            selfilter = 'Abinit (*.out)'
        elif self.settings['Program'] == 'gulp':
            selfilter = 'Gulp (*.gout)'
        elif self.settings['Program'] == 'vasp':
            selfilter = 'VASP (OUTCAR*)'
        elif self.settings['Program'] == 'qe':
            selfilter = 'QE (*.dynG)'
        elif self.settings['Program'] == 'crystal':
            selfilter = 'Crystal 14 (*.out)'
        elif self.settings['Program'] == 'phonopy':
            selfilter = 'Phonopy (OUTCAR*)'
        elif self.settings['Program'] == 'experiment':
            selfilter = 'Experiment (*.exp)'
        elif self.settings['Program'] == 'pdgui':
            selfilter = 'PDGui (*.py)'
        else:
            selfilter = 'All files (*)'
        filename,myfilter = QFileDialog.getOpenFileName(self,'Open MM/QM Output file','','Abinit (*.out);;Castep (*.castep);;Crystal 14 (*.out);;Experiment (*.exp);;Gulp (*.gout);;Phonopy (OUTCAR*);;QE (*.dynG);;VASP (OUTCAR*);;PDGui (*.py);;All files(*)',selfilter)
        # Process the filename
        if filename != '':
            program = find_program_from_name(filename)
            if program == '':
                debugger.print('Program not found from filename',filename)
                debugger.print('Finished:: on_file_le_return ')
                return
            self.settings['Program'] = program
            self.directory = os.path.dirname(filename)
            self.notebook.app.setMyWindowTitle(self.directory)
            self.settings['Output file name'] = os.path.basename(filename)
            self.file_le.setText(self.settings['Output file name'])
            debugger.print('new file name', self.directory, self.settings['Output file name'])
            self.notebook.deleteAllScenarios()
            if self.settings['Program'] == 'pdgui':
                #
                # If we have changed the file then reset some of the settingsTab settings
                # before we read in the script.
                #
                if self.notebook.settingsTab is not None:
                    self.notebook.settingsTab.settings['Optical permittivity'] = None
                    self.notebook.settingsTab.settings['Mass definition'] = 'average'
                    self.notebook.settingsTab.sigmas_cm1 = []
                    self.notebook.settingsTab.requestRefresh()
                self.notebook.app.readScript(filename)
            else:
                #
                # If we have changed the file then reset some of the settingsTab settings
                #
                if self.notebook.settingsTab is not None:
                    self.notebook.settingsTab.settings['Optical permittivity'] = None
                    self.notebook.settingsTab.settings['Mass definition'] = 'average'
                    self.notebook.settingsTab.sigmas_cm1 = []
                    self.notebook.settingsTab.requestRefresh()
                self.refreshRequired = True
                self.calculationRequired = True
                self.refresh()
        debugger.print('Finished:: on_file_button_clicked ', self.file_le.text())
        return

    def getFullFileName(self):
        return os.path.join(self.directory,self.settings['Output file name'])

    def getRelativeFileName(self):
        return self.settings['Output file name']

    def on_file_le_return(self):
        debugger.print('Start:: on_file_le_return ', self.file_le.text())
        filename = self.file_le.text()
        if filename != '':
            program = find_program_from_name(filename)
            if program == '':
                debugger.print('Finished:: on_file_le_return ')
                return
            self.settings['Program'] = program
            self.directory = os.path.dirname(os.path.abspath(filename))
            self.notebook.app.setMyWindowTitle(self.directory)
            self.settings['Output file name'] = os.path.basename(filename)
            debugger.print('new file name', self.settings['Output file name'])
            self.notebook.deleteAllScenarios()
            self.settings['Program'] = find_program_from_name(self.settings['Output file name'])
            if self.settings['Program'] == 'pdgui':
                #
                # If we are reading a script the reset the defaults before reading the script
                #
                if self.notebook.settingsTab is not None:
                    self.notebook.settingsTab.settings['Optical permittivity'] = None
                    self.notebook.settingsTab.settings['Mass definition'] = 'average'
                    self.notebook.settingsTab.sigmas_cm1 = []
                self.notebook.app.readScript(filename)
            else:
                #
                # If we have changed the file then reset some of the settingsTab settings
                # this is only done if it is not a script
                #
                if self.notebook.settingsTab is not None:
                    self.notebook.settingsTab.settings['Optical permittivity'] = None
                    self.notebook.settingsTab.settings['Mass definition'] = 'average'
                    self.notebook.settingsTab.sigmas_cm1 = []
                self.refreshRequired = True
                self.calculationRequired = True
                self.refresh()
        debugger.print('Finished:: on_file_le_return ', self.file_le.text())
        return

    def on_file_le_changed(self, text):
        debugger.print('Start:: on_file_changed', text)
        self.settings['Output file name'] = text
        self.calculationRequired = True
        debugger.print('Finished:: on_file_changed', text)

    def on_program_cb_activated(self,index):
        debugger.print('Start:: on_program_combobox_activated', index)
        debugger.print('on_program_combobox_activated', self.program_cb.currentText())
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
        self.calculationRequired = True
        debugger.print('Finished:: on_program_combobox_activated', index)

    def requestRefresh(self):
        debugger.print('Start:: requestRefresh')
        self.refreshRequired = True
        debugger.print('Finished:: requestRefresh')

    def refresh(self,force=False):
        debugger.print('Start:: refresh', force)
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
        self.file_le.setText(self.settings['Output file name'])
        self.resultsfile_le.setText(self.settings['Excel file name'])
        if self.calculationRequired:
            debugger.print('on_calculation_button_clicked called from MainTab.refresh()')
            self.on_calculation_button_clicked()
        self.refreshRequired = False
        #
        # UnBlock signals
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        debugger.print('Finished:: refresh', force)
        return
