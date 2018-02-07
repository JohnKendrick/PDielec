import sys
import os.path
from PyQt5.QtWidgets  import  QMainWindow, QApplication
from PyQt5.QtWidgets  import  QPushButton, QWidget, QAction, QTabWidget
from PyQt5.QtWidgets  import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QFileDialog, QCheckBox
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QMessageBox, QGroupBox, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox, QTableWidget, QTableWidgetItem
from PyQt5.QtGui      import  QIcon
from PyQt5.QtCore     import  pyqtSlot
from PyQt5.QtCore     import  Qt
from Python.Utilities import  get_reader
from Python.Constants import  support_matrix_db
 
class SettingsTab(QWidget):
 
    def __init__(self, parent, MainTab):   
        super(QWidget, self).__init__(parent)
        self.settings = {}
        self.settings["vmin"] = 0
        self.settings["vmax"] = 400
        self.settings["mass"] = "average"
        matrix = "ptfe"
        self.settings["matrix"] = matrix
        self.settings["matrix_density"] = support_matrix_db[matrix][0]
        self.settings["matrix_permittivity"] = support_matrix_db[matrix][1]
        self.reader = MainTab.reader
        self.layout = QVBoxLayout()
        # Create second tab - SettingsTab
        vbox = QVBoxLayout(self)
        # The mass definition combobox
        box = QFormLayout()
        #self.mass_l  = QLabel("Atomic Mass Definition:",self)
        #self.mass_l.setToolTip("The atomic masses used to calculate frequencies and intensities can be give here")
        self.mass_cb = QComboBox(self)
        self.mass_cb.setToolTip("The atomic masses used to calculate frequencies and intensities can be give here")
        self.mass_cb.addItem('Average natural abundance')
        self.mass_cb.addItem('Mass taken from QM/MM program')
        self.mass_cb.addItem('Most common isotope mass')
        index = self.mass_cb.findText(self.settings["mass"], Qt.MatchFixedString)
        if index >=0:
            self.mass_cb.setCurrentIndex(index)
        else:
            print("mass index was not 0")
        self.mass_cb.currentIndexChanged.connect(self.on_mass_cb_changed)
        #hbox.addWidget(self.mass_l)
        #hbox.addWidget(self.mass_cb)
        #vbox.addLayout(hbox)
        box.addRow(QLabel("Atomic mass defintion:", self), self.mass_cb)
        # Create Table containing the masses
        self.element_masses_tw = QTableWidget(self)
        self.element_masses_tw.setToolTip("Individual element masses can be modified here")
        self.element_masses_tw.itemClicked.connect(self.on_element_masses_tw_itemClicked)
        self.element_masses_tw.setRowCount(1)
        if self.reader:
            masses = self.reader.masses_per_type
            self.element_masses_tw.setColumnCount(len(masses))
            for i,mass in enumerate(masses):
                self.element_masses_tw.setItem(0,i, str(mass) )
        box.addRow(QLabel("Atomic masses", self), self.element_masses_tw)
        vbox.addLayout(box)
        # Support matrix
        hbox = QHBoxLayout()
        self.matrix_l  = QLabel("Support matrix:",self)
        self.matrix_l.setToolTip("Define the permittivity and density of the support matrix")
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
        hbox.addWidget(self.matrix_l)
        hbox.addWidget(self.matrix_cb)
        vbox.addLayout(hbox)
        # Support matrix density
        hbox = QHBoxLayout()
        self.density_l = QLabel("Support density:", self)
        self.density_l.setToolTip("Define the support matrix density")
        self.density_le = QLineEdit(self) 
        self.density_le.setToolTip("Define the support matrix density")
        self.density_le.setText(str(self.settings["matrix_density"]))
        self.density_le.textChanged.connect(self.on_density_le_changed)
        hbox.addWidget(self.density_l)
        hbox.addWidget(self.density_le)
        vbox.addLayout(hbox)
        # Support matrix permittivity
        hbox = QHBoxLayout()
        self.permittivity_l = QLabel("Support permittivity:", self)
        self.permittivity_l.setToolTip("Define the support matrix permittivity")
        self.permittivity_le = QLineEdit(self) 
        self.permittivity_le.setToolTip("Define the support matrix permittivity")
        self.permittivity_le.setText(str(self.settings["matrix_permittivity"]))
        self.permittivity_le.textChanged.connect(self.on_permittivity_le_changed)
        hbox.addWidget(self.permittivity_l)
        hbox.addWidget(self.permittivity_le)
        vbox.addLayout(hbox)

        # The minimum frequency
        hbox = QHBoxLayout()
        self.vmin_l = QLabel("Minimum frequency:", self)
        self.vmin_l.setToolTip("Set the minimum frequency to be considered)")
        self.vmin_sb = QSpinBox(self)
        self.vmin_sb.setRange(-200,9000)
        self.vmin_sb.setValue(self.settings["vmin"])
        self.vmin_sb.setToolTip("Set the minimum frequency to be considered)")
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        hbox.addWidget(self.vmin_l)
        hbox.addWidget(self.vmin_sb)
        vbox.addLayout(hbox)
        # The maximum frequency
        hbox = QHBoxLayout()
        self.vmax_l = QLabel("Maximum frequency:", self)
        self.vmax_l.setToolTip("Set the maximum frequency to be considered)")
        self.vmax_sb = QSpinBox(self)
        self.vmax_sb.setRange(-200,9000)
        self.vmax_sb.setValue(self.settings["vmax"])
        self.vmax_sb.setToolTip("Set the maximum frequency to be considered)")
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        hbox.addWidget(self.vmax_l)
        hbox.addWidget(self.vmax_sb)
        vbox.addLayout(hbox)
        # Optical permittivity
        hbox = QHBoxLayout()
        self.optical_l = QLabel("Optical permittivity")
        self.optical_l.setToolTip("The optical permittivity is taken from the calculation where this is possible.  If it is not availble suitbale values should be provided here")
        self.optical_tw = QTableWidget(3,3)
        self.optical_tw.itemClicked.connect(self.on_optical_tw_itemClicked)
        if self.reader:
            optical = self.reader.zerof_optical_dielectric
            for i,row in enumerate(optical):
                for j, value in enumerate(row):
                    self.element_masses_tw.setItem(i,j, str(value) )
        # Support matrix
        hbox.addWidget(self.optical_l)
        hbox.addWidget(self.optical_tw)
        vbox.addLayout(hbox)
        # Final button
        self.pushButton1 = QPushButton("Read Output File")
        self.pushButton1.setToolTip("Read the output file specified and list the phonon frequencies found")
        self.pushButton1.clicked.connect(self.pushButton1Clicked)
        vbox.addWidget(self.pushButton1)
        # output window
        self.listw_l = QLabel("Frequencies from:", self)
        vbox.addWidget(self.listw_l)
        self.listw = QListWidget(self)
        vbox.addWidget(self.listw)
        # finalise the layout
        self.setLayout(vbox)

    def pushButton1Clicked(self):
        print("Button 1 pressed")
        self.reader.neutral = self.settings["neutral"]
        self.reader.eckart = self.settings["eckart"]
        self.reader.hessian_symmetrisation = self.settings["hessian_symmetrisation"]
        self.reader.read_output()
        self.reader.print_info()
        frequencies_cm1 = self.reader.frequencies
        for f in frequencies_cm1:
            self.listw.addItem(str(f))

    def on_vmin_changed(self):
        self.settings["vmin"] = self.vmin_sb.value()
        print("on vmin change ", self.settings["vmin"])

    def on_vmax_changed(self):
        self.settings["vmax"] = self.vmin_sb.value()
        print("on vmax change ", self.settings["vmax"])

    def on_mass_cb_changed(self,index):
        print("on mass combobox changed", index)
        print("on mass combobox changed", self.mass_cb.currentText())
        options = ["average","program","isotope"]
        self.settings["mass"] = options[index]
        # Modify the element masses
        print("on_mass_change",self.reader)
        if self.reader:
            masses = self.reader.masses_per_type
            self.element_masses_tw.setColumnCount(len(masses))
            for i,mass in enumerate(masses):set_reader
                qw = QTableWidgetItem()
                qw.setText(str(mass))
                self.element_masses_tw.setItem(0,i, qw )

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

    def on_element_masses_tw_itemClicked(self, item):
        print("on_element_masses_tw_itemClicked)", item)
        print("on_element_masses_tw_itemClicked)", item.row())
        print("on_element_masses_tw_itemClicked)", item.column())

    def on_optical_tw_itemClicked(self, item):
        print("on_optical_tw)", item)
        print("on_optical)", item.row())
        print("on_optical)", item.column())
