import sys
import os.path
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QMainWindow, QApplication
from PyQt5.QtWidgets  import  QPushButton, QWidget, QAction, QTabWidget
from PyQt5.QtWidgets  import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QFileDialog, QCheckBox
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QMessageBox, QGroupBox, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox, QTableWidget, QTableWidgetItem
from PyQt5.QtWidgets  import  QSizePolicy, QSpacerItem
from PyQt5.QtGui      import  QIcon
from PyQt5.QtCore     import  pyqtSlot
from PyQt5.QtCore     import  Qt, QSize
from Python.Utilities import  get_reader
from Python.Constants import  support_matrix_db, wavenumber, amu
from Python.Constants import  average_masses, isotope_masses

class FixedQTableWidget(QTableWidget):
    def __init__(self, *args, parent=None):   
        super(QTableWidget, self).__init__(*args)

    def sizeHint(self):
        width = 0
        for i in range(self.columnCount()):
            print("column width",i,self.columnWidth(i))
            width += self.columnWidth(i)
        width += self.verticalHeader().sizeHint().width()
        width += self.verticalScrollBar().sizeHint().width()
        width += self.frameWidth()*2
        height = 0
        for i in range(self.rowCount()):
            print("row height",i,self.rowHeight(i))
            height += self.rowHeight(i)
        height += self.verticalHeader().sizeHint().width()
        height += self.horizontalScrollBar().sizeHint().height()
        height += self.frameWidth()*2
        #jk self.resize(QSize(width,height))
        return QSize(width,height)
 
class SettingsTab(QWidget):
    def __init__(self, parent, MainTab):   
        super(QWidget, self).__init__(parent)
        self.settings = {}
        self.settings["eckart"] = True
        self.settings["neutral"] = False
        self.settings["hessian_symmetrisation"] = "symm"
        self.settings["sigma"] = 5
        self.settings["vmin"] = 0
        self.settings["vmax"] = 400
        self.settings["spreadsheet"] = ""
        matrix = "ptfe"
        self.settings["matrix"] = matrix
        self.settings["matrix_density"] = support_matrix_db[matrix][0]
        self.settings["matrix_permittivity"] = support_matrix_db[matrix][1]
        self.mass_definition_options = ["average","program","isotope","gui"]
        self.settings["mass_definition"] = "average"
        self.settings["masses_dictionary"] = {}
        self.reader = MainTab.reader
        # Create second tab - SettingsTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The eckart checkbox
        #
        self.eckart_cb = QCheckBox(self)
        self.eckart_cb.setToolTip("Applying Eckart conditions ensures three zero translation mode)")
        self.eckart_cb.setText("")
        self.eckart_cb.setLayoutDirection(Qt.RightToLeft)
        self.eckart_cb.stateChanged.connect(self.on_eckart_changed)
        if self.settings["eckart"]:
            self.eckart_cb.setCheckState(Qt.Checked)
        else:
            self.eckart_cb.setCheckState(Qt.Unchecked)
        form.addRow(QLabel("Apply Eckart conditions?", self), self.eckart_cb)
        #
        # Add the Born neutral condition
        #
        self.born_cb = QCheckBox(self)
        self.born_cb.setToolTip("Applying Born charge neutrality ensures unit cell has zero charge")
        self.born_cb.setText("")
        self.born_cb.setLayoutDirection(Qt.RightToLeft)
        self.born_cb.stateChanged.connect(self.on_born_changed)
        if self.settings["neutral"]:
            self.born_cb.setCheckState(Qt.Checked)
        else:
            self.born_cb.setCheckState(Qt.Unchecked)
        form.addRow(QLabel("Apply Born charge neutrality?",self),self.born_cb)
        #
        # The mass definition combobox
        #
        self.mass_cb = QComboBox(self)
        self.mass_cb.setToolTip("The atomic masses used to calculate frequencies and intensities can be give here")
        self.mass_cb.addItem('Average natural abundance')
        self.mass_cb.addItem('Mass taken from QM/MM program')
        self.mass_cb.addItem('Most common isotope mass')
        # set default to average natural abundance
        self.mass_cb.setCurrentIndex(0)
        self.mass_cb.currentIndexChanged.connect(self.on_mass_cb_changed)
        form.addRow(QLabel("Atomic mass defintion:", self), self.mass_cb)
        # Create Table containing the masses - block signals until the table is loaded
        self.element_masses_tw = FixedQTableWidget(self)
        self.element_masses_tw.setToolTip("Individual element masses can be modified here")
        self.element_masses_tw.itemClicked.connect(self.on_element_masses_tw_itemClicked)
        self.element_masses_tw.itemChanged.connect(self.on_element_masses_tw_itemChanged)
        self.element_masses_tw.setRowCount(1)
        self.element_masses_tw.blockSignals(True)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        self.element_masses_tw.setSizePolicy(sizePolicy)
        print("element_masses_tw sizeHint",self.element_masses_tw.sizeHint())
        form.addRow(QLabel("Atomic masses", self), self.element_masses_tw)
        #form.addItem(QSpacerItem(20, 40, QSizePolicy.Expanding, QSizePolicy.Expanding) )
        #
        # Support matrix
        #
        self.matrix_cb = QComboBox(self)
        self.matrix_cb.setToolTip("Define the permittivity and density of the support matrix")
        for matrix in support_matrix_db:
            print("matrix:",matrix, support_matrix_db[matrix])
            self.matrix_cb.addItem(matrix)
        index = self.matrix_cb.findText(self.settings["matrix"], Qt.MatchFixedString)
        if index >=0:
            self.matrix_cb.setCurrentIndex(index)
        else:
            print("support matrix index was not 0")
        self.matrix_cb.currentIndexChanged.connect(self.on_matrix_cb_changed)
        form.addRow(QLabel("Support matrix:",self), self.matrix_cb)
        self.density_le = QLineEdit(self) 
        self.density_le.setToolTip("Define the support matrix density")
        self.density_le.setText(str(self.settings["matrix_density"]))
        self.density_le.textChanged.connect(self.on_density_le_changed)
        form.addRow(QLabel("Support density", self), self.density_le)
        #
        # Support matrix permittivity
        #
        self.permittivity_le = QLineEdit(self) 
        self.permittivity_le.setToolTip("Define the support matrix permittivity")
        self.permittivity_le.setText(str(self.settings["matrix_permittivity"]))
        self.permittivity_le.textChanged.connect(self.on_permittivity_le_changed)
        form.addRow(QLabel("Support permittivity", self), self.permittivity_le)
        #
        # The lorentzian width - sigma
        #
        self.sigma_sb = QSpinBox(self)
        self.sigma_sb.setRange(1,100)
        self.sigma_sb.setValue(self.settings["sigma"])
        self.sigma_sb.setToolTip("Set the default Lorentzian width factor (sigma) in cm-1")
        self.sigma_sb.valueChanged.connect(self.on_sigma_changed)
        form.addRow(QLabel("Lorentzian width (sigma):", self), self.sigma_sb)
        #
        # The minimum frequency
        #
        self.vmin_sb = QSpinBox(self)
        self.vmin_sb.setRange(-200,9000)
        self.vmin_sb.setValue(self.settings["vmin"])
        self.vmin_sb.setToolTip("Set the minimum frequency to be considered)")
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        form.addRow(QLabel("Minimum frequency:", self), self.vmin_sb)
        #
        # The maximum frequency
        #
        self.vmax_sb = QSpinBox(self)
        self.vmax_sb.setRange(-200,9000)
        self.vmax_sb.setValue(self.settings["vmax"])
        self.vmax_sb.setToolTip("Set the maximum frequency to be considered)")
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        form.addRow(QLabel("Maximum frequency:", self), self.vmax_sb)
        #
        # Create the Optical permittivity table widget and block signals until a click on the widget
        #
        self.optical_tw = FixedQTableWidget(3,3)
        self.optical_tw.setToolTip("The optical permittivity is taken from the calculation where this is possible.  If it is not availble suitbale values should be provided here")
        # Set the header names
        self.optical_tw.setHorizontalHeaderLabels(["x","y","z"])
        self.optical_tw.setVerticalHeaderLabels  (["x","y","z"])
        self.optical_tw.itemClicked.connect(self.on_optical_tw_itemClicked)
        self.optical_tw.itemChanged.connect(self.on_optical_tw_itemChanged)
        self.optical_tw.blockSignals(True)
        self.optical_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel("Optical permittivity:", self), self.optical_tw)
        # 
        # Store results in a file?
        #
        self.file_store_le = QLineEdit(self) 
        self.file_store_le.setToolTip("Store the results in a .csv or .xlsx file")
        self.file_store_le.setText(str(self.settings["spreadsheet"]))
        self.file_store_le.textChanged.connect(self.on_file_store_le_changed)
        form.addRow(QLabel("Output spreadsheet", self), self.file_store_le)
        #
        # Final button
        #
        self.pushButton1 = QPushButton("Process and calculate intensities")
        self.pushButton1.setToolTip("Process and calculate intensities using the MM/QM phonon output")
        self.pushButton1.clicked.connect(self.pushButton1Clicked)
        form.addRow(self.pushButton1)
        vbox.addLayout(form)
        # output window
        # Create Table containing the IR active modes
        self.output_tw = FixedQTableWidget(self)
        self.output_tw.setToolTip("Output showing the frequencies and strengths of the IR active modes only")
        self.output_tw.itemClicked.connect(self.on_output_tw_itemClicked)
        self.output_tw.itemChanged.connect(self.on_output_tw_itemChanged)
        self.output_tw.setRowCount(1)
        self.output_tw.blockSignals(True)
        vbox.addWidget(self.output_tw)
        # finalise the layout
        self.setLayout(vbox)

    def pushButton1Clicked(self):
        print("Button 1 pressed")
        self.reader.read_output()
        self.reader.neutral = self.settings["neutral"]
        self.reader.eckart = self.settings["eckart"]
        self.reader.hessian_symmetrisation = self.settings["hessian_symmetrisation"]
        mass_dictionary = []
        self.reader.reset_masses()
        if self.settings["mass_definition"] == "average":
            self.reader.change_masses(average_masses, mass_dictionary)
        elif self.settings["mass_definition"] == "isotope":
            self.reader.change_masses(isotope_masses, mass_dictionary)
        elif self.settings["mass_definition"] == "gui":
            self.reader.change_masses(self.settings["masses_dictionary"], mass_dictionary)
        elif self.settings["mass_definition"] == "program":
            pass
        else:
            print("Error unkown mass definition", self.settings["mass_definition"] )
        mass_weighted_normal_modes = self.reader.calculate_mass_weighted_normal_modes()
        # convert sigmas to wavenumbers
        self.frequencies_cm1 = self.reader.frequencies
        self.sigmas = [ self.settings["sigma"] for i in self.frequencies_cm1 ]
        sigmas_cm1 = np.array(self.sigmas)
        sigmas = np.array(self.sigmas)*wavenumber
        born_charges = np.array(self.reader.born_charges)
        if self.reader.type == 'Experimental output':
            oscillator_strengths = self.reader.oscillator_strengths
        else:
            #
            # calculate normal modes in xyz coordinate space
            masses = np.array(self.reader.masses) * amu
            normal_modes = Calculator.normal_modes(masses, mass_weighted_normal_modes)
            # from the normal modes and the born charges calculate the oscillator strengths of each mode
            oscillator_strengths = Calculator.oscillator_strengths(normal_modes, born_charges)
        # calculate the intensities from the trace of the oscillator strengths
        intensities = Calculator.infrared_intensities(oscillator_strengths)
 
        self.output_tw.setRowCount(len(self.sigmas))
        self.output_tw.setColumnCount(4)
        self.output_tw.setHorizontalHeaderLabels(['Include?', 'Frequency', 'sigma', 'Intensity'])
        for i,(f,sigma,intensity) in enumerate(zip(self.frequencies_cm1, self.sigmas, intensities)):
            self.output_tw.setItem(i, 0, QTableWidgetItem(" ") )
            self.output_tw.setItem(i, 1, QTableWidgetItem(str(f) ) )
            self.output_tw.setItem(i, 2, QTableWidgetItem(str(sigma) ) )
            self.output_tw.setItem(i, 3, QTableWidgetItem(str(intensity) ) )

    def on_file_store_le_changed(self,text):
        self.settings["spreadsheet"] = text
        print("on file_store_le change ", self.settings["spreadsheet"])

    def on_sigma_changed(self):
        self.settings["sigma"] = self.sigma_sb.value()
        print("on sigma change ", self.settings["sigma"])

    def on_vmin_changed(self):
        self.settings["vmin"] = self.vmin_sb.value()
        print("on vmin change ", self.settings["vmin"])

    def on_vmax_changed(self):
        self.settings["vmax"] = self.vmin_sb.value()
        print("on vmax change ", self.settings["vmax"])

    def on_mass_cb_changed(self,index):
        print("on mass combobox changed", index)
        print("on mass combobox changed", self.mass_cb.currentText())
        self.settings["mass_definition"] = self.mass_definition_options[index]
        # Modify the element masses
        self.set_masses_tw()

    def set_masses_tw(self):
        print("set masses")
        if self.reader:
            self.element_masses_tw.blockSignals(True)
            species = self.reader.species
            if self.settings["mass_definition"] == "gui":
                self.settings["mass_definition"] == "average"
                self.mass_cb.setCurrentIndex(0)
            # set the initial dictionary according to the mass_definition
            masses = []
            self.masses_dictionary = {}
            if self.settings["mass_definition"] == "average":
                for element in species:
                    mass = average_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings["mass_definition"] == "isotope":
                for element in species:
                    mass = isotope_masses[element]
                    masses.append(mass)
                    self.masses_dictionary[element] = mass
            elif self.settings["mass_definition"] == "program":
                 self.reader.reset_masses()
                 masses = self.reader.masses_per_type
                 for mass,element in zip(masses,species):
                    self.masses_dictionary[element] = mass
            elif self.settings["mass_definition"] == "gui":
                 pass
            else:
                 print("Error mass_definition not recognised", self.settings["mass_definition"])
            self.element_masses_tw.setColumnCount(len(masses))
            self.element_masses_tw.setHorizontalHeaderLabels(species)
            self.element_masses_tw.setVerticalHeaderLabels([""])
            # set masses of the elements in the table widget according to the mass definition
            for i,(mass,element) in enumerate(zip(masses,species)):
                print("set_masses_tw", self.settings["mass_definition"],i,mass,element)
                qw = QTableWidgetItem()
                if self.settings["mass_definition"] == "program":
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(str(mass))
                    self.element_masses_tw.setItem(0,i, qw )
                elif self.settings["mass_definition"] == "average":
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(str(average_masses[element]))
                    print("average",average_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                elif  self.settings["mass_definition"] == "isotope":
                    self.element_masses_tw.blockSignals(True)
                    qw.setText(str(isotope_masses[element]))
                    print("isotope",isotope_masses[element])
                    self.element_masses_tw.setItem(0,i, qw )
                else:
                    print("mass_definition not processed", self.settings["mass_definition"])
            # unblock the table signals
            self.element_masses_tw.blockSignals(False)

    def on_matrix_cb_changed(self,index):
        print("on matrix combobox changed", index)
        print("on matrix combobox changed", self.matrix_cb.currentText())
        matrix = self.matrix_cb.currentText()
        self.settings["matrix"] = matrix
        self.settings["matrix_density"] = support_matrix_db[matrix][0]
        self.settings["matrix_permittivity"] = support_matrix_db[matrix][1]
        self.density_le.setText(str(self.settings["matrix_density"]))
        self.permittivity_le.setText(str(self.settings["matrix_permittivity"]))

    def on_density_le_changed(self,text):
        self.settings["density"] = float(text)
        print("on density line edit changed", text)

    def on_permittivity_le_changed(self,text):
        self.settings["permittivity"] = float(text)
        print("on density line edit changed", text)

    def set_reader(self,reader):
        self.reader = reader

    def on_output_tw_itemClicked(self, item):
        print("on_output_tw_itemClicked)", item)
        print("on_output_tw_itemClicked)", item.row())
        print("on_output_tw_itemClicked)", item.column())
        self.output_tw.blockSignals(False)

    def on_output_tw_itemChanged(self, item):
        print("on_output_tw_itemChanged)", item)
        print("on_output_tw_itemChanged)", item.row())
        print("on_output_tw_itemChanged)", item.column())
        col = item.column()
        row = item.row()
        new_value = float(item.text())
        print("output new value", self.settings["sigma"], new_value)

    def on_element_masses_tw_itemClicked(self, item):
        print("on_element_masses_tw_itemClicked)", item)
        print("on_element_masses_tw_itemClicked)", item.row())
        print("on_element_masses_tw_itemClicked)", item.column())
        self.element_masses_tw.blockSignals(False)

    def on_element_masses_tw_itemChanged(self, item):
        print("on_element_masses_tw_itemChanged)", item)
        print("on_element_masses_tw_itemChanged)", item.row())
        print("on_element_masses_tw_itemChanged)", item.column())
        col = item.column()
        self.settings["mass_definition"] = "gui"
        self.settings["masses_dictionary"][self.reader.species[col]] = float(item.text())
        print("masses_dictionary", self.settings["masses_dictionary"])

    def on_optical_tw_itemChanged(self, item):
        print("on_optical_itemChanged)", item)
        print("on_optical_itemChanged)", item.row())
        print("on_optical_itemChanged)", item.column())
        self.settings["optical_permittivity"][item.row()][item.column()] = float(item.text())
        self.settings["optical_permittivity"][item.column()][item.row()] = float(item.text())
        self.set_optical_permittivity_tw()
        print("optical permittivity")
        print(self.settings["optical_permittivity"])

    def on_optical_tw_itemClicked(self, item):
        print("on_optical_itemClicked)", item)
        print("on_optical_itemClicked)", item.row())
        print("on_optical_itemClicked)", item.column())
        self.optical_tw.blockSignals(False)

    def refresh(self):
        # Refresh the widgets that depend on the reader
        if self.reader:
            # Masses
            self.set_masses_tw()
            print("element_masses_tw sizeHint",self.element_masses_tw.sizeHint())
            print("element_masses_tw size",self.element_masses_tw.size())
            # Optical dielectric
            self.set_optical_permittivity_tw()
            print("optical_tw sizeHint",self.optical_tw.sizeHint())
            print("optical_tw size",self.optical_tw.size())

    def set_optical_permittivity_tw(self):
        optical = self.reader.zerof_optical_dielectric
        self.settings["optical_permittivity"] = optical
        self.optical_tw.blockSignals(True)
        for i,row in enumerate(optical):
            for j, value in enumerate(row):
                self.optical_tw.setItem(i,j, QTableWidgetItem(str(value) ))
        self.optical_tw.blockSignals(False)

    def on_born_changed(self):
        print("on born change ", self.born_cb.isChecked())
        self.settings["neutral"] = self.born_cb.isChecked()
        print("on born change ", self.settings["neutral"])

    def on_eckart_changed(self):
        print("on eckart change ", self.eckart_cb.isChecked())
        self.settings["eckart"] = self.eckart_cb.isChecked()
        print("on eckart change ", self.settings["eckart"])

