# -*- coding: utf8 -*-
import sys
import os.path
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QCheckBox
from PyQt5.QtWidgets  import  QVBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox, QTableWidget, QTableWidgetItem
from PyQt5.QtWidgets  import  QSizePolicy
from PyQt5.QtCore     import  Qt, QSize
from Python.Utilities import  get_reader
from Python.Constants import  wavenumber, amu, PI
from Python.Constants import  average_masses, isotope_masses
from Python.Utilities import  Debug

class FixedQTableWidget(QTableWidget):
    def __init__(self, *args, parent=None):   
        super(QTableWidget, self).__init__(*args)

    def sizeHint(self):
        width = 0
        for i in range(self.columnCount()):
            width += self.columnWidth(i)
        width += self.verticalHeader().sizeHint().width()
        width += self.verticalScrollBar().sizeHint().width()
        width += self.frameWidth()*2
        height = 0
        for i in range(self.rowCount()):
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
        self.notebook = parent
        self.settings = {}
        self.settings['eckart'] = True
        self.settings['neutral'] = False
        self.settings['hessian_symmetrisation'] = 'symm'
        self.settings['sigma'] = 5
        self.mass_definition_options = ['average','program','isotope','gui']
        self.settings['mass_definition'] = 'average'
        self.masses_dictionary = {}
        self.modes_selected = []
        self.frequencies_cm1 = []
        self.intensities = []
        self.sigmas_cm1 = []
        self.oscillator_strengths = []
        self.mass_weighted_normal_modes = None
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
        if self.settings['eckart']:
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
        if self.settings['neutral']:
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
        # set default to average natural abundance
        self.mass_cb.currentIndexChanged.connect(self.on_mass_cb_changed)
        self.mass_cb.setCurrentIndex(0)
        form.addRow(QLabel('Atomic mass defintion:', self), self.mass_cb)
        # Create Table containing the masses - block signals until the table is loaded
        self.element_masses_tw = FixedQTableWidget(self)
        self.element_masses_tw.setToolTip('Individual element masses can be modified here')
        self.element_masses_tw.itemClicked.connect(self.on_element_masses_tw_itemClicked)
        self.element_masses_tw.itemChanged.connect(self.on_element_masses_tw_itemChanged)
        self.element_masses_tw.setRowCount(1)
        self.element_masses_tw.blockSignals(True)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        self.element_masses_tw.setSizePolicy(sizePolicy)
        debugger.print('element_masses_tw sizeHint',self.element_masses_tw.sizeHint())
        form.addRow(QLabel('Atomic masses', self), self.element_masses_tw)
        #
        # The lorentzian width - sigma
        #
        self.sigma_sb = QSpinBox(self)
        self.sigma_sb.setRange(1,100)
        self.sigma_sb.setValue(self.settings['sigma'])
        self.sigma_sb.setToolTip('Set the default Lorentzian width factor (σ) in cm-1')
        self.sigma_sb.valueChanged.connect(self.on_sigma_changed)
        label = QLabel('Lorentzian width (σ):', self)
        label.setToolTip('Set the default Lorentzian width factor (σ) in cm-1')
        form.addRow(label, self.sigma_sb)
        #
        # Create the Optical permittivity table widget and block signals until a click on the widget
        #
        self.optical_tw = FixedQTableWidget(3,3)
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
        # Final button - need for this has been removed refreshes are done after every change
        #
        #self.calculateButton = QPushButton('Process and calculate intensities')
        #self.calculateButton.setToolTip('Process and calculate intensities using the MM/QM phonon output')
        #self.calculateButton.clicked.connect(self.calculateButtonClicked)
        #form.addRow(self.calculateButton)
        vbox.addLayout(form)
        # output window
        # Create Table containing the IR active modes
        self.output_tw = FixedQTableWidget(self)
        self.output_tw.setToolTip('Output showing the frequencies and strengths of the IR active modes only')
        self.output_tw.itemClicked.connect(self.on_output_tw_itemClicked)
        self.output_tw.itemChanged.connect(self.on_output_tw_itemChanged)
        self.output_tw.setRowCount(1)
        self.output_tw.blockSignals(True)
        vbox.addWidget(self.output_tw)
        # finalise the layout
        self.setLayout(vbox)

    def calculateButtonClicked(self):
        debugger.print('Button 1 pressed')
        self.reader = self.notebook.reader
        # Only calculate if the reader is set
        if self.reader is None:
            return
        self.reader.hessian_symmetrisation = self.settings['hessian_symmetrisation']
        # Do I still need to re-read the output file ?????
        #self.reader.read_output()
        if self.settings['neutral']:
            self.reader.neutralise_born_charges()
        self.reader.eckart = self.settings['eckart']
        mass_dictionary = []
        self.reader.reset_masses()
        if self.settings['mass_definition'] == 'average':
            self.reader.change_masses(average_masses, mass_dictionary)
        elif self.settings['mass_definition'] == 'isotope':
            self.reader.change_masses(isotope_masses, mass_dictionary)
        elif self.settings['mass_definition'] == 'gui':
            self.reader.change_masses(self.masses_dictionary, mass_dictionary)
        elif self.settings['mass_definition'] == 'program':
            pass
        else:
            print('Error unkown mass definition', self.settings['mass_definition'] )
        self.mass_weighted_normal_modes = self.reader.calculate_mass_weighted_normal_modes()
        # convert cm-1 to au
        self.frequencies_cm1 = self.reader.frequencies
        frequencies = np.array(self.frequencies_cm1) * wavenumber
        self.sigmas_cm1 = [ self.settings['sigma'] for i in self.frequencies_cm1 ]
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
        for f,intensity in zip(self.frequencies_cm1,self.intensities):
            if f > 10.0 and intensity > 1.0E-6:
                self.modes_selected.append(True)
            else:
                self.modes_selected.append(False)
        self.output_tw.setRowCount(len(self.sigmas_cm1))
        self.output_tw.setColumnCount(5)
        self.output_tw.setHorizontalHeaderLabels(['   Sigma   \n(cm-1)', ' Frequency \n(cm-1)', '  Intensity  \n(Debye2/Å2/amu)', 'Integrated Molar Absorption\n(L/mole/cm2)', 'Absorption maximum\n(L/mole/cm)'])
        self.redraw_output_tw()

    def redraw_output_tw(self):
        # If the frequencies haven't been set yet just don't try to do anything
        if not self.frequencies_cm1:
            return
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
                item.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
                self.output_tw.setItem(i, j, item )
        # Release the block on signals for the frequency output table
        self.output_tw.resizeColumnsToContents()
        self.output_tw.blockSignals(False)

    def on_sigma_changed(self):
        self.settings['sigma'] = self.sigma_sb.value()
        self.sigmas_cm1 = [ self.settings['sigma'] for i in self.frequencies_cm1 ]
        self.redraw_output_tw()
        debugger.print('on sigma change ', self.settings['sigma'])

    def on_mass_cb_changed(self,index):
        debugger.print('on mass combobox changed', index)
        debugger.print('on mass combobox changed', self.mass_cb.currentText())
        self.settings['mass_definition'] = self.mass_definition_options[index]
        # Modify the element masses
        self.set_masses_tw()
        self.calculateButtonClicked()

    def set_masses_tw(self):
        debugger.print('set masses')
        if self.reader:
            self.element_masses_tw.blockSignals(True)
            species = self.reader.species
            if self.settings['mass_definition'] == 'gui':
                self.settings['mass_definition'] == 'average'
                self.mass_cb.setCurrentIndex(0)
            # set the initial dictionary according to the mass_definition
            masses = []
            if self.settings['mass_definition'] == 'average':
                for element in species:
                    mass = average_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings['mass_definition'] == 'isotope':
                for element in species:
                    mass = isotope_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings['mass_definition'] == 'program':
                 self.reader.reset_masses()
                 masses = self.reader.masses_per_type
                 for mass,element in zip(masses,species):
                    self.masses_dictionary[element] = mass
            elif self.settings['mass_definition'] == 'gui':
                 pass
            else:
                 debugger.print('Error mass_definition not recognised', self.settings['mass_definition'])
            self.element_masses_tw.setColumnCount(len(masses))
            self.element_masses_tw.setHorizontalHeaderLabels(species)
            self.element_masses_tw.setVerticalHeaderLabels([''])
            debugger.print('masses_dictionary',self.masses_dictionary)
            debugger.print('masses',masses)
            # set masses of the elements in the table widget according to the mass definition
            for i,(mass,element) in enumerate(zip(masses,species)):
                debugger.print('set_masses_tw', self.settings['mass_definition'],i,mass,element)
                qw = QTableWidgetItem()
                if self.settings['mass_definition'] == 'program':
                    self.element_masses_tw.blockSignals(True)
                    qw.setText('{0:.6f}'.format(mass))
                    qw.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
                    self.element_masses_tw.setItem(0,i, qw )
                elif self.settings['mass_definition'] == 'average':
                    self.element_masses_tw.blockSignals(True)
                    qw.setText('{0:.6f}'.format(average_masses[element]))
                    qw.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
                    debugger.print('average',average_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings['mass_definition'] == 'isotope':
                    self.element_masses_tw.blockSignals(True)
                    qw.setText('{0:.6f}'.format(isotope_masses[element]))
                    qw.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
                    debugger.print('isotope',isotope_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                else:
                    print('mass_definition not processed', self.settings['mass_definition'])
            # unblock the table signals
            self.element_masses_tw.blockSignals(False)

    def on_output_tw_itemClicked(self, item):
        debugger.print('on_output_tw_itemClicked)', item.row(), item.column() )
        #col = item.column()
        #row = item.row()
        # Handle the case that an unchecked (greyed out) sigma is clicked
        #self.output_tw.blockSignals(True)
        #if col == 0:
        #    if not self.modes_selected[row]:
        #        debugger.print('on_output_tw_item_clicked changing mode selection', row)
        #        self.modes_selected[row] = True
        #        debugger.print('on_output_tw_item_clicked mode selected is', self.modes_selected[row])
        #        item.setFlags(item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable )
        #        item.setCheckState(Qt.Checked)
        #self.output_tw.blockSignals(False)

    def on_output_tw_itemChanged(self, item):
        self.output_tw.blockSignals(True)
        debugger.print('on_output_tw_itemChanged)', item.row(), item.column() )
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

    def on_element_masses_tw_itemClicked(self, item):
        debugger.print('on_element_masses_tw_itemClicked)', item.row(),item.column() )
        self.element_masses_tw.blockSignals(False)

    def on_element_masses_tw_itemChanged(self, item):
        debugger.print('on_element_masses_tw_itemChanged)', item.row(), item.column() )
        col = item.column()
        self.settings['mass_definition'] = 'gui'
        self.masses_dictionary[self.reader.species[col]] = float(item.text())
        self.calculateButtonClicked()
        debugger.print('masses_dictionary', self.masses_dictionary)

    def on_optical_tw_itemChanged(self, item):
        debugger.print('on_optical_itemChanged)', item.row(), item.column() )
        self.settings['optical_permittivity'][item.row()][item.column()] = float(item.text())
        self.settings['optical_permittivity'][item.column()][item.row()] = float(item.text())
        debugger.print('optical permittivity')
        debugger.print(self.settings['optical_permittivity'])

    def on_optical_tw_itemClicked(self, item):
        debugger.print('on_optical_itemClicked)', item.row(), item.column() )
        self.optical_tw.blockSignals(False)

    def refresh(self):
        # Refresh the widgets that depend on the reader
        self.reader = self.notebook.reader
        if self.reader:
            # Masses
            self.set_masses_tw()
            debugger.print('element_masses_tw sizeHint',self.element_masses_tw.sizeHint())
            debugger.print('element_masses_tw size',self.element_masses_tw.size())
            # Optical dielectric
            self.set_optical_permittivity_tw()
            debugger.print('optical_tw sizeHint',self.optical_tw.sizeHint())
            debugger.print('optical_tw size',self.optical_tw.size())
        self.calculateButtonClicked()

    def set_optical_permittivity_tw(self):
        optical = self.reader.zerof_optical_dielectric
        self.settings['optical_permittivity'] = optical
        self.optical_tw.blockSignals(True)
        for i,row in enumerate(optical):
            for j, value in enumerate(row):
                qw = QTableWidgetItem('{0:.4f}'.format(value) )
                qw.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
                self.optical_tw.setItem(i,j,qw)
        self.optical_tw.blockSignals(False)

    def on_born_changed(self):
        debugger.print('on born change ', self.born_cb.isChecked())
        self.settings['neutral'] = self.born_cb.isChecked()
        debugger.print('on born change ', self.settings['neutral'])
        self.calculateButtonClicked()

    def on_eckart_changed(self):
        debugger.print('on eckart change ', self.eckart_cb.isChecked())
        self.settings['eckart'] = self.eckart_cb.isChecked()
        debugger.print('on eckart change ', self.settings['eckart'])
        self.calculateButtonClicked()

