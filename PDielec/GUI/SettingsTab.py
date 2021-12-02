# -*- coding: utf8 -*-
import numpy as np
import PDielec.Calculator as Calculator
from PyQt5.QtWidgets            import  QWidget, QApplication
from PyQt5.QtWidgets            import  QComboBox, QLabel
from PyQt5.QtWidgets            import  QCheckBox
from PyQt5.QtGui                import  QFont
from PyQt5.QtWidgets            import  QVBoxLayout, QFormLayout
from PyQt5.QtWidgets            import  QDoubleSpinBox, QTableWidget, QTableWidgetItem
from PyQt5.QtWidgets            import  QSizePolicy
from PyQt5.QtCore               import  Qt, QSize, QCoreApplication
from PDielec.Constants          import  wavenumber, amu, PI, angstrom
from PDielec.Constants          import  average_masses, isotope_masses
from PDielec.Utilities          import  Debug
from PDielec.DielectricFunction import  DielectricFunction
from functools                  import  partial

class FixedQTableWidget(QTableWidget):
    def __init__(self, *args, parent=None, rows=None, columns=None):
        self.columns = columns
        self.rows = rows
        super(QTableWidget, self).__init__(*args)
        # fnt = QFont()
        # fnt.setPointSize(12)
        # fnt.setFamily("Arial")
        # self.setFont(fnt)

    def sizeHint(self):
        width = 0
        if self.columns == None:
            columns = self.columnCount()
        else:
            columns = self.columns
        for i in range(columns):
            width += self.columnWidth(i)
        width += self.verticalHeader().sizeHint().width()
        width += self.verticalScrollBar().sizeHint().width()
        width += self.frameWidth()*2
        height = 0
        if self.rows == None:
            rows = self.rowCount()
        else:
            rows = self.rows
        for i in range(rows):
            height += self.rowHeight(i)
        height += self.verticalHeader().sizeHint().width()
        height += self.horizontalScrollBar().sizeHint().height()
        height += self.frameWidth()*2
        return QSize(width,height)

class SettingsTab(QWidget):
    def __init__(self, parent, debug=False):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug, 'SettingsTab:')
        debugger.print('Start:: initialisation')
        self.notebook = parent
        self.refreshRequired = True
        self.calculationRequired = False
        self.settings = {}
        self.settings['Eckart flag'] = True
        self.settings['Neutral Born charges'] = False
        self.settings['Sigma value'] = 5
        self.mass_definition_options = ['average','program','isotope','gui']
        self.settings['Mass definition'] = 'average'
        self.settings['Optical permittivity'] = None
        self.settings['Optical permittivity edited'] = False
        self.masses_dictionary = {}
        self.modes_selected = []
        self.frequencies_cm1 = []
        self.intensities = []
        self.sigmas_cm1 = []
        self.oscillator_strengths = []
        self.mass_weighted_normal_modes = None
        self.CrystalPermittivity = None
        self.vs_cm1 = []
        self.crystal_permittivity = []
        # get the reader from the main tab
        self.reader = self.notebook.reader
        # Create second tab - SettingsTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The eckart checkbox
        #
        self.eckart_cb = QCheckBox(self)
        self.eckart_cb.setToolTip('Applying Eckart conditions ensures three zero translation mode)')
        self.eckart_cb.setText('')
        self.eckart_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings['Eckart flag']:
            self.eckart_cb.setCheckState(Qt.Checked)
        else:
            self.eckart_cb.setCheckState(Qt.Unchecked)
        self.eckart_cb.stateChanged.connect(self.on_eckart_changed)
        if self.settings['Eckart flag']:
            self.eckart_cb.setCheckState(Qt.Checked)
        else:
            self.eckart_cb.setCheckState(Qt.Unchecked)
        self.eckart_cb.stateChanged.connect(self.on_eckart_changed)
        form.addRow(QLabel('Apply Eckart conditions?', self), self.eckart_cb)
        #
        # Add the Born neutral condition
        #
        self.born_cb = QCheckBox(self)
        self.born_cb.setToolTip('Applying Born charge neutrality ensures unit cell has zero charge')
        self.born_cb.setText('')
        self.born_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings['Neutral Born charges']:
            self.born_cb.setCheckState(Qt.Checked)
        else:
            self.born_cb.setCheckState(Qt.Unchecked)
        self.born_cb.stateChanged.connect(self.on_born_changed)
        form.addRow(QLabel('Apply Born charge neutrality?',self),self.born_cb)
        #
        # The mass definition combobox
        #
        self.mass_cb = QComboBox(self)
        self.mass_cb.setToolTip('The atomic masses used to calculate frequencies and intensities can be give here')
        self.mass_cb.addItem('Average natural abundance')
        self.mass_cb.addItem('Mass taken from QM/MM program')
        self.mass_cb.addItem('Most common isotope mass')
        self.mass_cb.addItem('Masses set individually')
        self.mass_cb.model().item(3).setEnabled(False)
        # set default to average natural abundance
        self.mass_cb.activated.connect(self.on_mass_cb_activated)
        self.current_mass_definition_index = self.mass_definition_options.index(self.settings['Mass definition'])
        self.mass_cb.setCurrentIndex(0)
        form.addRow(QLabel('Atomic mass defintion:', self), self.mass_cb)
        # Create Table containing the masses - block signals until the table is loaded
        self.element_masses_tw = FixedQTableWidget(parent=self)
        self.element_masses_tw.setToolTip('Individual element masses can be modified here')
        self.element_masses_tw.itemClicked.connect(self.on_element_masses_tw_itemClicked)
        self.element_masses_tw.itemChanged.connect(self.on_element_masses_tw_itemChanged)
        self.element_masses_tw.setRowCount(1)
        self.element_masses_tw.blockSignals(True)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        self.element_masses_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel('Atomic masses', self), self.element_masses_tw)
        #
        # The lorentzian width - sigma
        #
        self.sigma_sb = QDoubleSpinBox(self)
        self.sigma_sb.setRange(0.1,100.0)
        self.sigma_sb.setSingleStep(0.1)
        self.sigma_sb.setDecimals(1)
        self.sigma_sb.setValue(self.settings['Sigma value'])
        self.sigma_sb.setToolTip('Set the default Lorentzian width factor (σ) in cm-1')
        self.sigma_sb.valueChanged.connect(self.on_sigma_changed)
        label = QLabel('Lorentzian width (σ):', self)
        label.setToolTip('Set the default Lorentzian width factor (σ) in cm-1')
        form.addRow(label, self.sigma_sb)
        #
        # Create the Optical permittivity table widget and block signals until a click on the widget
        #
        self.optical_tw = FixedQTableWidget(3,3,parent=self)
        self.optical_tw.setToolTip('The optical permittivity is taken from the calculation where this is possible.  If it is not availble suitbale values should be provided here')
        # Set the header names
        self.optical_tw.setHorizontalHeaderLabels(['x','y','z'])
        self.optical_tw.setVerticalHeaderLabels  (['x','y','z'])
        self.optical_tw.itemClicked.connect(self.on_optical_tw_itemClicked)
        self.optical_tw.itemChanged.connect(self.on_optical_tw_itemChanged)
        # Block the widget until the optical permittivity is loaded
        self.optical_tw.blockSignals(True)
        self.optical_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel('Optical permittivity:', self), self.optical_tw)
        #
        vbox.addLayout(form)
        # output window
        # Create Table containing the IR active modes
        self.output_tw = FixedQTableWidget(parent=self)
        self.output_tw.setToolTip('Output showing the frequencies and strengths of the IR active modes only')
        self.output_tw.itemChanged.connect(self.on_output_tw_itemChanged)
        self.output_tw.setRowCount(1)
        self.output_tw.blockSignals(True)
        vbox.addWidget(self.output_tw)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        debugger.print('Finished:: initialisation')

    def setElementMass(self,element,mass):
        debugger.print('Start::  setElementMass',element,mass)
        self.settings['Mass definition'] = 'gui'
        self.masses_dict[element] = mass
        self.mass_cb.setCurrentIndex(3)
        self.set_masses_tw()
        self.refreshRequired = True
        self.refresh()
        debugger.print('Finished::  setElementMass',element,mass)

    def createIntensityTable(self):
        debugger.print('Start:: createIntensityTable')
        self.reader = self.notebook.reader
        # Only calculate if the reader is set
        if self.reader is None:
            debugger.print('createIntensityTable aborting as now reader available')
            return
        if self.settings['Neutral Born charges']:
            self.reader.neutralise_born_charges()
        else:
            self.reader.reset_born_charges()
        self.reader.eckart = self.settings['Eckart flag']
        mass_dictionary = []
        self.reader.reset_masses()
        if self.settings['Mass definition'] == 'average':
            self.reader.change_masses(average_masses, mass_dictionary)
        elif self.settings['Mass definition'] == 'program':
            pass
        elif self.settings['Mass definition'] == 'isotope':
            self.reader.change_masses(isotope_masses, mass_dictionary)
        elif self.settings['Mass definition'] == 'gui':
            self.reader.change_masses(self.masses_dictionary, mass_dictionary)
        else:
            print('Error unkown mass definition', self.settings['Mass definition'] )
        QCoreApplication.processEvents()
        self.mass_weighted_normal_modes = self.reader.calculate_mass_weighted_normal_modes()
        # convert cm-1 to au
        self.frequencies_cm1 = self.reader.frequencies
        frequencies = np.array(self.frequencies_cm1) * wavenumber
        if len(self.sigmas_cm1) == 0:
            self.sigmas_cm1 = [ self.settings['Sigma value'] for i in self.frequencies_cm1 ]
        born_charges = np.array(self.reader.born_charges)
        if self.reader.type == 'Experimental output':
            self.oscillator_strengths = self.reader.oscillator_strengths
        else:
            #
            # calculate normal modes in xyz coordinate space
            masses = np.array(self.reader.masses) * amu
            normal_modes = Calculator.normal_modes(masses, self.mass_weighted_normal_modes)
            # from the normal modes and the born charges calculate the oscillator strengths of each mode
            self.oscillator_strengths = Calculator.oscillator_strengths(normal_modes, born_charges)
        # calculate the intensities from the trace of the oscillator strengths
        self.intensities = Calculator.infrared_intensities(self.oscillator_strengths)
        # Decide which modes to select
        self.modes_selected = []
        mode_list = []
        for index,(f,intensity) in enumerate(zip(self.frequencies_cm1,self.intensities)):
            if f > 10.0 and intensity > 1.0E-6:
                self.modes_selected.append(True)
                mode_list.append(index)
            else:
                self.modes_selected.append(False)
        # Calculate the ionic contribution to the permittivity
        frequencies_au = wavenumber*np.array(self.frequencies_cm1)
        volume_au = self.reader.volume*angstrom*angstrom*angstrom
        self.epsilon_ionic = Calculator.ionic_permittivity(mode_list, self.oscillator_strengths, frequencies_au, volume_au )
        # Make an np array of epsilon infinity
        epsilon_inf = np.array(self.settings['Optical permittivity'])
        #
        # If the reader already has a Dielectric Constant then use this
        #
        if self.reader.CrystalPermittivity:
            self.CrystalPermittivity = self.reader.CrystalPermittivity
        else:
            # Otherwise create one from the data that has been read in
            drude = False
            drude_plasma_au = 0
            drude_sigma_au = 0
            sigmas_au = np.array(self.sigmas_cm1)*wavenumber
            self.CrystalPermittivity = DielectricFunction('dft',parameters=(
                                         mode_list, frequencies_au, sigmas_au, self.oscillator_strengths,
                                         volume_au, drude, drude_plasma_au, drude_sigma_au) )
        # Add the optical permittivity to the dielctric function
        self.CrystalPermittivity.setEpsilonInfinity(epsilon_inf)
        #
        # Prepare to finish
        #
        self.output_tw.setRowCount(len(self.sigmas_cm1))
        self.output_tw.setColumnCount(5)
        self.output_tw.setHorizontalHeaderLabels(['   Sigma   \n(cm-1)', ' Frequency \n(cm-1)', '  Intensity  \n(Debye2/Å2/amu)', 'Integrated Molar Absorption\n(L/mole/cm2)', 'Absorption maximum\n(L/mole/cm)'])
        QCoreApplication.processEvents()
        self.redraw_output_tw()
        QCoreApplication.processEvents()
        # if self.notebook.spreadsheet is not None:
        #     self.writeSpreadsheet()
        QCoreApplication.processEvents()
        debugger.print('Finished:: createIntensityTable')
        return

    def requestRefresh(self):
        debugger.print('Start:: requestRefresh')
        self.refreshRequired = True
        # self.refresh()
        debugger.print('Finished:: requestRefresh')
        return

    def writeSpreadsheet(self):
        debugger.print('Start:: writeSpreadsheet')
        sp = self.notebook.spreadsheet
        if sp is None:
            debugger.print('Finished:: writeSpreadsheet - Aborting write of spreadsheet')
            return
        debugger.print('Writing of spreadsheet')
        sp.selectWorkSheet('Settings')
        sp.delete
        sp.writeNextRow(['Settings and calculations of frequencies and absorption'], row=0, col=1)
        for item in sorted(self.settings):
            if item == 'Optical permittivity' and self.settings[item] is not None:
               sp.writeNextRow([item], col=1)
               sp.writeNextRow(self.settings[item][0], col=2, check=1)
               sp.writeNextRow(self.settings[item][1], col=2, check=1)
               sp.writeNextRow(self.settings[item][2], col=2, check=1)
               eps = self.epsilon_ionic.tolist()
               sp.writeNextRow(['Ionic Contribution to permittivity at zero frequency'], col=1)
               sp.writeNextRow(eps[0], col=2, check=1)
               sp.writeNextRow(eps[1], col=2, check=1)
               sp.writeNextRow(eps[2], col=2, check=1)
               eps = self.epsilon_ionic + np.array(self.settings[item])
               eps = eps.tolist()
               sp.writeNextRow(['Total Contribution to permittivity at zero frequency'], col=1)
               sp.writeNextRow(eps[0], col=2, check=1)
               sp.writeNextRow(eps[1], col=2, check=1)
               sp.writeNextRow(eps[2], col=2, check=1)
            elif item == 'Mass definition':
               sp.writeNextRow([item,self.settings[item]], col=1, check=1)
               for c in self.masses_dictionary:
                   sp.writeNextRow([c,self.masses_dictionary[c]],col=3,check=1)
            else:
               sp.writeNextRow([item,self.settings[item]], col=1, check=1)
        sp.writeNextRow([''], col=1)
        sp.writeNextRow(['Mode','Include?','Sigma(cm-1)','Frequency(cm-1)','Intensity(Debye2/Angs2/amu','Integrated Molar Absorption(/L/mole/cm2','Absorption maximum (L/mole/cm)'], col=1)
        for mode,(f,intensity,sigma,selected) in enumerate(zip(self.frequencies_cm1, self.intensities, self.sigmas_cm1,self.modes_selected)):
            yn = 'No'
            if selected:
                yn = 'Yes'
            sp.writeNextRow([mode, yn, sigma, f, intensity, 4225.6*intensity, 2*4225.6*intensity/sigma/PI], col=1)
        debugger.print('Finished:: writeSpreadsheet')

    def redraw_output_tw(self):
        debugger.print('Start:: redraw_output_tw')
        # If the frequencies haven't been set yet just don't try to do anything
        self.output_tw.blockSignals(True)
        for i,(f,sigma,intensity) in enumerate(zip(self.frequencies_cm1, self.sigmas_cm1, self.intensities)):
            # Sigma and check / unchecked column
            items = []
            itemFlags = []
            item = QTableWidgetItem('{0:.1f}'.format(sigma))
            if self.modes_selected[i]:
                item.setCheckState(Qt.Checked)
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable )
                otherFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled
            else:
                #itemFlags.append( item.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled & ~Qt.ItemIsSelectable & ~Qt.ItemIsEditable )
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled )
                item.setCheckState(Qt.Unchecked)
                otherFlags = item.flags() & Qt.NoItemFlags
            items.append(item)
            # Frequency column cm-1
            items.append(QTableWidgetItem('{0:.4f}'.format(f) ) )
            itemFlags.append( otherFlags )
            # Intensity column Debye2/Angs2/amu
            items.append(QTableWidgetItem('{0:.4f}'.format(intensity) ) )
            itemFlags.append( otherFlags )
            # Integrated molar absorption L/mole/cm/cm
            items.append(QTableWidgetItem('{0:.2f}'.format(intensity*4225.6) ) )
            itemFlags.append( otherFlags )
            # Maximum extinction L/mole/cm
            items.append(QTableWidgetItem('{0:.2f}'.format(2*intensity*4225.6/self.sigmas_cm1[i]/PI) ) )
            itemFlags.append( otherFlags )
            for j,(item,flag) in enumerate(zip(items,itemFlags)):
                item.setFlags(flag)
                item.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                self.output_tw.setItem(i, j, item )
        # Release the block on signals for the frequency output table
        self.output_tw.resizeColumnsToContents()
        self.output_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print('Finished:: redraw_output_tw')

    def on_sigma_changed(self):
        debugger.print('Start:: redraw_output_tw')
        self.settings['Sigma value'] = self.sigma_sb.value()
        self.sigmas_cm1 = [ self.settings['Sigma value'] for i in self.frequencies_cm1 ]
        self.redraw_output_tw()
        debugger.print('on sigma change ', self.settings['Sigma value'])
        self.refreshRequired = True
        self.refresh()
        debugger.print('Finished:: on_sigma_changed')

    def on_mass_cb_activated(self,index):
        debugger.print('Start:: on_mass_combobox_activated', self.mass_cb.currentText())
        self.settings['Mass definition'] = self.mass_definition_options[index]
        self.current_mass_definition_index = index
        if index < 3:
            self.mass_cb.model().item(3).setEnabled(False)
        # Modify the element masses
        self.set_masses_tw()
        self.refreshRequired = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print('Finished:: on_mass_combobox_activated', self.mass_cb.currentText())

    def set_masses_tw(self):
        debugger.print('Start:: set_masses_tw')
        if self.reader:
            self.element_masses_tw.blockSignals(True)
            species = self.reader.getSpecies()
            # set the initial dictionary according to the mass_definition
            masses = []
            if self.settings['Mass definition'] == 'average':
                self.mass_cb.setCurrentIndex(0)
                for element in species:
                    mass = average_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings['Mass definition'] == 'program':
                self.mass_cb.setCurrentIndex(1)
                self.reader.reset_masses()
                masses = self.reader.masses_per_type
                for mass,element in zip(masses,species):
                    self.masses_dictionary[element] = mass
            elif self.settings['Mass definition'] == 'isotope':
                self.mass_cb.setCurrentIndex(2)
                for element in species:
                    mass = isotope_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings['Mass definition'] == 'gui':
                for element in species:
                    mass = self.masses_dictionary[element]
                    masses.append(mass)
            else:
                 debugger.print('Error mass_definition not recognised', self.settings['Mass definition'])
            self.element_masses_tw.setColumnCount(len(masses))
            self.element_masses_tw.setHorizontalHeaderLabels(species)
            self.element_masses_tw.setVerticalHeaderLabels([''])
            debugger.print('masses_dictionary',self.masses_dictionary)
            debugger.print('masses',masses)
            # set masses of the elements in the table widget according to the mass definition
            for i,(mass,element) in enumerate(zip(masses,species)):
                debugger.print('set_masses_tw', self.settings['Mass definition'],i,mass,element)
                qw = QTableWidgetItem()
                if self.settings['Mass definition'] == 'program':
                    self.element_masses_tw.blockSignals(True)
                    qw.setText('{0:.6f}'.format(mass))
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    self.element_masses_tw.setItem(0,i, qw )
                elif self.settings['Mass definition'] == 'average':
                    self.element_masses_tw.blockSignals(True)
                    qw.setText('{0:.6f}'.format(average_masses[element]))
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    debugger.print('average',average_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings['Mass definition'] == 'isotope':
                    self.element_masses_tw.blockSignals(True)
                    qw.setText('{0:.6f}'.format(isotope_masses[element]))
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    debugger.print('isotope',isotope_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings['Mass definition'] == 'gui':
                    self.element_masses_tw.blockSignals(True)
                    qw.setText('{0:.6f}'.format(self.masses_dictionary[element]))
                    qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                    debugger.print('gui',self.masses_dictionary[element])
                    self.element_masses_tw.setItem(0,i, qw )
                else:
                    print('Mass definition not processed', self.settings['Mass definition'])
            # unblock the table signals
            self.element_masses_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print('Finished:: set_masses_tw')

    def on_output_tw_itemChanged(self, item):
        self.output_tw.blockSignals(True)
        debugger.print('Start:: on_output_tw_itemChanged)', item.row(), item.column() )
        col = item.column()
        row = item.row()
        if col == 0:
            # If this is the first column alter the check status but reset the sigma value
            if item.checkState() == Qt.Checked:
                self.modes_selected[row] = True
            else:
                 self.modes_selected[row] = False
            new_value = float(item.text())
            if new_value != self.sigmas_cm1:
                self.sigmas_cm1[row] = new_value
                self.redraw_output_tw()
        elif col == 1:
            self.redraw_output_tw()
        else:
            self.redraw_output_tw()
        QCoreApplication.processEvents()
        debugger.print('Finished:: on_output_tw_itemChanged')

    def on_element_masses_tw_itemClicked(self, item):
        debugger.print('Start:: on_element_masses_tw_itemClicked)', item.row(),item.column() )
        self.element_masses_tw.blockSignals(False)
        debugger.print('Finished:: on_element_masses_tw_itemClicked)' )

    def on_element_masses_tw_itemChanged(self, item):
        debugger.print('Start:: on_element_masses_tw_itemChanged)', item.row(), item.column() )
        elements = self.reader.getSpecies()
        col = item.column()
        self.mass_cb.model().item(3).setEnabled(True)
        self.settings['Mass definition'] = 'gui'
        self.mass_cb.setCurrentIndex(3)
        self.masses_dictionary[elements[col]] = float(item.text())
        debugger.print('masses_dictionary', self.masses_dictionary)
        self.refreshRequired = True
        self.refresh()
        debugger.print('Finished:: on_element_masses_tw_itemChanged)' )

    def on_optical_tw_itemChanged(self, item):
        debugger.print('Start::on_optical_itemChanged)', item.row(), item.column(), item.text() )
        self.settings['Optical permittivity'][item.row()][item.column()] = float(item.text())
        self.settings['Optical permittivity'][item.column()][item.row()] = float(item.text())
        self.settings['Optical permittivity edited'] = True
        self.refresh_optical_permittivity_tw()
        self.refreshRequired = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print('Finished::on_optical_itemChanged)')
        return

    def on_optical_tw_itemClicked(self, item):
        debugger.print('Start:: on_optical_itemClicked)', item.row(), item.column() )
        self.optical_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print('Finished:: on_optical_itemClicked)')
        return

    def refresh(self, force=False):
        debugger.print('Start:: refresh', force )
        if not self.reader and self.notebook.reader:
            self.refreshRequired = True
        if not self.refreshRequired:
            debugger.print('Finished:: refresh not required',force)
            return
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        self.reader = self.notebook.reader
        if self.reader:
            # Masses
            self.set_masses_tw()
            # Optical dielectric
            if self.settings['Optical permittivity'] is None:
                self.set_optical_permittivity_tw()
            else:
                self.refresh_optical_permittivity_tw()
        self.sigma_sb.setValue(self.settings['Sigma value'])
        if self.settings['Eckart flag']:
            self.eckart_cb.setCheckState(Qt.Checked)
        else:
            self.eckart_cb.setCheckState(Qt.Unchecked)
        if self.settings['Neutral Born charges']:
            self.born_cb.setCheckState(Qt.Checked)
        else:
            self.born_cb.setCheckState(Qt.Unchecked)
        #
        # Unlock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.createIntensityTable()
        self.refreshRequired = False
        # Refresh the scenarios that depend on this dielectric
        if self.notebook.scenarios is not None:
            for scenario in self.notebook.scenarios:
                scenario.requestRefresh()
        if self.notebook.plottingTab is not None:
            self.notebook.plottingTab.requestRefresh()
        if self.notebook.analysisTab is not None:
            self.notebook.analysisTab.requestRefresh()
        if self.notebook.fitterTab is not None:
            self.notebook.fitterTab.requestRefresh()
        QCoreApplication.processEvents()
        self.refreshRequired = False
        self.calculationRequired = True
        debugger.print('Finished:: refresh', force )
        return

    def refresh_optical_permittivity_tw(self):
        debugger.print('Start:: refresh_optical_permittivity')
        optical = self.settings['Optical permittivity']
        self.optical_tw.blockSignals(True)
        for i,row in enumerate(optical):
            for j, value in enumerate(row):
                qw = QTableWidgetItem('{0:.4f}'.format(value) )
                qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                self.optical_tw.setItem(i,j,qw)
        self.optical_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print('Finished:: refresh_optical_permittivity')
        return

    def set_optical_permittivity_tw(self):
        debugger.print('Start:: set_optical_permittivity_tw')
        self.settings['Optical permittivity'] = self.reader.zerof_optical_dielectric
        self.refresh_optical_permittivity_tw()
        self.refreshRequired = True
        QCoreApplication.processEvents()
        debugger.print('Finished:: set_optical_permittivity_tw')
        return

    def on_born_changed(self):
        debugger.print('Start:: on_born_change ', self.born_cb.isChecked())
        self.settings['Neutral Born charges'] = self.born_cb.isChecked()
        debugger.print('on born change ', self.settings['Neutral Born charges'])
        self.refreshRequired = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print('Finished:: on_born_change ', self.born_cb.isChecked())
        return

    def on_eckart_changed(self):
        debugger.print('Start:: on_eckart_change ', self.eckart_cb.isChecked())
        self.settings['Eckart flag'] = self.eckart_cb.isChecked()
        debugger.print('on eckart change ', self.settings['Eckart flag'])
        self.refreshRequired = True
        self.refresh()
        QCoreApplication.processEvents()
        debugger.print('Finished:: on_eckart_change ', self.eckart_cb.isChecked())
        return

    def calculate(self,vs_cm1):
        """Calculate the permittivity of the crystal over the range of frequencies in vs_cm1"""
        debugger.print('Start:: calculate ')
        if len(vs_cm1) == 0:
            debugger.print('Finished:: calculate aborted vs_cm1 has not been set')
            return
        self.vs_cm1 = vs_cm1.copy()
        self.CrystalPermittivity.setUnits('cm-1')
        dielectricFunction = self.CrystalPermittivity.function()
        # 
        self.crystal_permittivity = []
        debugger.print('About to calculate settings crystal dielectric using pool')
        # Loop over the frequencies and calculate the crystal dielectric for each frequency
        for dielecv in self.notebook.pool.imap(dielectricFunction, vs_cm1, chunksize=20):
            self.crystal_permittivity.append(dielecv)
            self.notebook.progressbars_update()
        QCoreApplication.processEvents()
        self.calculationRequired = False
        debugger.print('Finished:: calculate ')
        return

    def get_crystal_permittivity(self,vs_cm1):
        """Return the crystal permittivity"""
        debugger.print('Start:: get_crystal_permittivity', self.refreshRequired)
        if self.calculationRequired or self.refreshRequired or  ( len(self.vs_cm1) != len(vs_cm1) ) or ( self.vs_cm1[0] != vs_cm1[0] ) or ( self.vs_cm1[1] != vs_cm1[1] ) :
            debugger.print('get_crystal_permittivity refreshing and recalculating' )
            self.refresh()
            self.calculate(vs_cm1)
        debugger.print('Finished:: get_crystal_permittivity', self.refreshRequired)
        return self.crystal_permittivity

