import os.path
import os
import numpy as np
from multiprocess import Array
import PDielec.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QProgressBar, QApplication
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox,QDoubleSpinBox
from PyQt5.QtWidgets  import  QSizePolicy
from PyQt5.QtCore     import  QCoreApplication, Qt
from PDielec.Constants import  wavenumber, PI, avogadro_si, angstrom
import ctypes
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PDielec.Utilities import Debug

def set_affinity_on_worker():
    '''When a new worker process is created, the affinity is set to all CPUs'''
    #JK print('I'm the process %d, setting affinity to all CPUs.' % os.getpid())
    #JK Commented out for the time being
    #JK os.system('taskset -p 0xff %d > /dev/null' % os.getpid())


class SingleCrystalTab(QWidget):
    def __init__(self, parent, debug=False ):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'SingleCrystalTab')
        self.dirty = True
        self.settings = {}
        self.subplot = None
        self.setWindowTitle('SingleCrystal')
        self.settings['Minimum frequency'] = 0
        self.settings['Maximum frequency'] = 400
        self.settings['Frequency increment'] = 0.2
        self.molar_definitions = ['Unit cells','Atoms','Molecules']
        self.settings['Molar definition'] = 'Unit cells'
        self.settings['Number of atoms'] = 1
        self.settings['Plot title'] = 'Plot Title'
        self.settings['Unique direction - h'] = 0
        self.settings['Unique direction - k'] = 0
        self.settings['Unique direction - l'] = 1
        self.settings['Rotation about Z'] = 0.0
        self.settings['Angle of incidence'] = 90.0
        self.settings['Polarisation angle'] = 0.0
        self.settings['Superstrate depth'] = 100.0
        self.settings['Superstrate dielectric'] = 1.0
        self.settings['Substrate depth'] = 10.0
        self.settings['Substrate dielectric'] = 1.0
        self.settings['Crystal depth'] = 1.0
        self.xaxes = []
        self.directions = []
        self.frequency_units = None
        self.molar_cb_current_index = 0
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Create last tab - SingleCrystalTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The minimum frequency
        #
        self.vmin_sb = QSpinBox(self)
        self.vmin_sb.setRange(0,9000)
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmin_sb.setToolTip('Set the minimum frequency to be considered)')
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        label = QLabel('Minimum frequency:', self)
        label.setToolTip('Set the minimum frequency to be considered)')
        form.addRow(label, self.vmin_sb)
        #
        # The maximum frequency
        #
        self.vmax_sb = QSpinBox(self)
        self.vmax_sb.setRange(0,9000)
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.vmax_sb.setToolTip('Set the maximum frequency to be considered)')
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        label = QLabel('Maximum frequency:', self)
        label.setToolTip('Set the maximum frequency to be considered)')
        form.addRow(label, self.vmax_sb)
        #
        # Choose a suitable increment
        #
        self.vinc_sb = QDoubleSpinBox(self)
        self.vinc_sb.setRange(0.0001,5.0)
        self.vinc_sb.setSingleStep(0.1)
        self.vinc_sb.setDecimals(4)
        self.vinc_sb.setToolTip('Choose an increment for the frequency when plotting')
        self.vinc_sb.setValue(self.settings['Frequency increment'])
        self.vinc_sb.valueChanged.connect(self.on_vinc_changed)
        label = QLabel('Frequency increment', self)
        label.setToolTip('Choose an increment for the frequency when plotting')
        form.addRow(label, self.vinc_sb)
        #
        # Define molar quantity
        #
        self.molar_cb = QComboBox(self)
        self.molar_cb.setToolTip('Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given')
        self.molar_cb.addItems(self.molar_definitions)
        try:
            self.molar_cb_current_index = self.molar_definitions.index(self.settings["Molar definition"])
        except:
            self.molar_cb_current_index = 0
            self.settings["Molar definition"] = self.molar_definitions[self.molar_cb_current_index]
        self.molar_cb.setCurrentIndex(self.molar_cb_current_index)
        self.molar_cb.activated.connect(self.on_molar_cb_activated)
        label = QLabel('Molar definition', self)
        label.setToolTip('Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given')
        form.addRow(label, self.molar_cb)
        #
        # Number of atoms in a molecule
        #
        self.natoms_sb = QSpinBox(self)
        self.natoms_sb.setToolTip('Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed')
        self.natoms_sb.setRange(1,500)
        self.natoms_sb.setValue(self.settings['Number of atoms'])
        self.natoms_sb.valueChanged.connect(self.on_natoms_changed)
        self.natoms_sb.setEnabled(False)
        label = QLabel('Number of atoms per molecule', self)
        label.setToolTip('Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed')
        form.addRow(label, self.natoms_sb)
        #
        # Set the plot title
        #
        self.title_le = QLineEdit(self)
        self.title_le.setToolTip('Set the plot title')
        self.title_le.setText(self.settings['Plot title'])
        self.title_le.textChanged.connect(self.on_title_changed)
        label = QLabel('Plot title', self)
        label.setToolTip('Set the plot title')
        form.addRow(label, self.title_le)
        #
        # Set the x-axis frequency units
        #
        self.funits_cb = QComboBox(self)
        self.funits_cb.setToolTip('Set the frequency units for the x-axis')
        self.funits_cb.addItems( ['wavenumber','THz'] )
        self.frequency_units = 'wavenumber'
        self.funits_cb.activated.connect(self.on_funits_cb_activated)
        label = QLabel('Frequency units for the x-axis', self)
        label.setToolTip('Set the frequency units for the x-axis')
        form.addRow(label, self.funits_cb)
        #
        # Define the slab surface in crystal coordinates
        #
        self.h_sb = QSpinBox(self)
        self.h_sb.setToolTip('Define the h dimension of the unique direction')
        self.h_sb.setRange(-20,20)
        self.h_sb.setSingleStep(1)
        self.h_sb.setValue(self.settings['Unique direction - h'])
        self.h_sb.valueChanged.connect(self.on_h_sb_changed)
        self.k_sb = QSpinBox(self)
        self.k_sb.setToolTip('Define the k dimension of the unique direction')
        self.k_sb.setRange(-20,20)
        self.k_sb.setSingleStep(1)
        self.k_sb.setValue(self.settings['Unique direction - k'])
        self.k_sb.valueChanged.connect(self.on_k_sb_changed)
        self.l_sb = QSpinBox(self)
        self.l_sb.setToolTip('Define the l dimension of the unique direction')
        self.l_sb.setRange(-20,20)
        self.l_sb.setSingleStep(1)
        self.l_sb.setValue(self.settings['Unique direction - l'])
        self.l_sb.valueChanged.connect(self.on_l_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.h_sb)
        hbox.addWidget(self.k_sb)
        hbox.addWidget(self.l_sb)
        self.hkl_label = QLabel('Crystal surface (hkl)',self)
        self.hkl_label.setToolTip('Define the crystal surface (hkl). Defines the unique direction in crystallographic units.')
        form.addRow(self.hkl_label, hbox)
        #
        # Define the rotation angle of the slab and the angle of incidence
        #
        self.rotation_about_z_sb = QSpinBox(self)
        self.rotation_about_z_sb.setToolTip('Define the slab rotation angle around the Z-laboratory frame')
        self.rotation_about_z_sb.setRange(0,360)
        self.rotation_about_z_sb.setSingleStep(10)
        self.rotation_about_z_sb.setValue(self.settings['Rotation around Z'])
        self.rotation_about_z_sb.valueChanged.connect(self.on_rotation_about_z_sb_changed)
        self.angle_of_incidence_sb = QSpinBox(self)
        self.angle_of_incidence_sb.setToolTip('Define the angle of incidence')
        self.angle_of_incidence_sb.setRange(0,180)
        self.angle_of_incidence_sb.setSingleStep(10)
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        self.angle_of_incidence_sb.valueChanged.connect(self.on_angle_of_incidence_sb_changed)
        self.polarisation_angle_sb = QSpinBox(self)
        self.polarisation_angle_sb.setToolTip('Define the polarisation angle')
        self.polarisation_angle_sb.setRange(0,90)
        self.polarisation_angle_sb.setSingleStep(10)
        self.polarisation_angle_sb.setValue(self.settings['Polarisation angle'])
        self.polarisation_angle_sb.valueChanged.connect(self.on_polarisation_angle_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.rotation_about_z_sb)
        hbox.addWidget(self.angle of incidence_sb)
        hbox.addWidget(self.polarisation_angle_sb)
        self.angles_label = QLabel('Angles of surface orientation, incidence and polarisation',self)
        self.angles_label.setToolTip('Define the angles of rotation about the Z-laboratory frame, the angle of incidence and the angle of polarisation relative.')
        form.addRow(self.angles_label, hbox)
        #
        # Define the superstrate
        #
        self.superstrate_depth_sb = QSpinBox(self)
        self.superstrate_depth_sb.setToolTip('Define the superstrate depth')
        self.superstrate_depth_sb.setRange(0,1000)
        self.superstrate_depth_sb.setSingleStep(10)
        self.superstrate_depth_sb.setValue(self.settings['Superstrate depth'])
        self.superstrate_depth_sb.valueChanged.connect(self.on_superstrate_depth_sb_change)
        self.superstrate_depth_sb.setToolTip('Define the superstrate dielectric')
        self.superstrate_depth_sb.setRange(1,1000)
        self.superstrate_depth_sb.setSingleStep(0.1)
        self.superstrate_depth_sb.setValue(self.settings['Superstrate dielectric'])
        self.superstrate_depth_sb.valueChanged.connect(self.on_superstrate_dielectric_sb_change)
        hbox = QHBoxLayout()
        hbox.addWidget(self.superstrate_depth_sb_changed)
        hbox.addWidget(self.angle of incidence_sb)
        self.angles_label = QLabel('Superstrate depth(mm) and dielectric',self)
        self.angles_label.setToolTip('Define the properties of the superstrate; depth in millimetres and dielectric constant.')
        form.addRow(self.angles_label, hbox)
        #
        #
        # Define the substrate
        #
        self.substrate_depth_sb = QSpinBox(self)
        self.substrate_depth_sb.setToolTip('Define the substrate depth')
        self.substrate_depth_sb.setRange(0,1000)
        self.substrate_depth_sb.setSingleStep(10)
        self.substrate_depth_sb.setValue(self.settings['Substrate depth'])
        self.substrate_depth_sb.valueChanged.connect(self.on_substrate_depth_sb_change)
        self.substrate_depth_sb.setToolTip('Define the substrate dielectric')
        self.substrate_depth_sb.setRange(1,1000)
        self.substrate_depth_sb.setSingleStep(0.1)
        self.substrate_depth_sb.setValue(self.settings['Substrate dielectric'])
        self.substrate_depth_sb.valueChanged.connect(self.on_substrate_dielectric_sb_change)
        hbox = QHBoxLayout()
        hbox.addWidget(self.substrate_depth_sb_changed)
        hbox.addWidget(self.angle of incidence_sb)
        self.angles_label = QLabel('Substrate depth(mm) and dielectric',self)
        self.angles_label.setToolTip('Define the properties of the substrate; depth in millimetres and dielectric constant.')
        form.addRow(self.angles_label, hbox)
        #
        #
        # Define the crystal depth
        #
        self.crystal_depth_sb = QSpinBox(self)
        self.crystal_depth_sb.setToolTip('Define the crystal depth')
        self.crystal_depth_sb.setRange(0,1000)
        self.crystal_depth_sb.setSingleStep(1)
        self.crystal_depth_sb.setValue(self.settings['Crystal depth'])
        self.crystal_depth_sb.valueChanged.connect(self.on_crystal_depth_sb_change)
        hbox = QHBoxLayout()
        hbox.addWidget(self.substrate_depth_sb_changed)
        hbox.addWidget(self.angle of incidence_sb)
        self.angles_label = QLabel('Crystal depth(mm)',self)
        self.angles_label.setToolTip('Define the depth of the crystal in millimetres.')
        form.addRow(self.angles_label, hbox)
        #
        #
        # Final button
        #
        hbox = QHBoxLayout()
        self.molarAbsorptionButton = QPushButton('Plot molar absorption')
        self.molarAbsorptionButton.setToolTip('Plot the molar absorption')
        self.molarAbsorptionButton.clicked.connect(self.molarAbsorptionButtonClicked)
        hbox.addWidget(self.molarAbsorptionButton)
        self.absorptionButton = QPushButton('Plot absorption')
        self.absorptionButton.setToolTip('Plot the absorption')
        self.absorptionButton.clicked.connect(self.absorptionButtonClicked)
        hbox.addWidget(self.absorptionButton)
        self.reflectivityButton = QPushButton('Plot relectivity')
        self.reflectivityButton.setToolTip('Plot the reflectivity')
        self.reflectivityButton.clicked.connect(self.reflectivityButtonCicked)
        hbox.addWidget(self.reflectivityButton)
        form.addRow(hbox)
        # Add a progress bar
        self.progressbar = QProgressBar(self)
        self.progressbar.setToolTip('Show the progress of any calculations')
        self.progressbar.setMinimum(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setMinimum(0)
        label = QLabel('Calculation progress', self)
        label.setToolTip('Show the progress of any calculations')
        form.addRow(label,self.progressbar)
        # Add the matplotlib figure to the bottom
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding))
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        form.addRow(self.canvas)
        form.addRow(self.toolbar)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()

    def self.on_superstrate_depth_sb_changed(self,value)
        debugger.print('on_superstrate_depth_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Superstrate depth'] = value

    def self.on_superstrate_dielectric_sb_changed(self,value)
        debugger.print('on_superstrate_dielectric_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Superstrate dielectric'] = value

    def self.on_substrate_depth_sb_changed(self,value)
        debugger.print('on_substrate_depth_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Substrate depth'] = value

    def self.on_substrate_dielectric_sb_changed(self,value)
        debugger.print('on_substrate_dielectric_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Substrate dielectric'] = value

    def self.on_rotation_about_z_sb_changed(self,value)
        debugger.print('on_rotation_about_z_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Rotation about Z'] = value

    def self.on_angle_of_incidence_sb_changed(self,value)
        debugger.print('on_angle_of_incidence_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Angle of incidence'] = value

    def self.on_polarisation_angle_sb_changed(self,value)
        debugger.print('on_polarisation_angle_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Polarisation angle'] = value

    def on_h_sb_changed(self,value):
        debugger.print('on_h_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Unique direction - h'] = value

    def on_k_sb_changed(self,value):
        debugger.print('on_k_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Unique direction - k'] = value

    def on_l_sb_changed(self,value):
        debugger.print('on_l_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Unique direction - l'] = value


    def molarAbsorptionButtonClicked(self):
        debugger.print('molarAbsorptionButtonClicked pressed')
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        if not self.notebook.singleCrystalCalculationRequired:
            self.plot(self.xaxes, self.molarAbsorptionCoefficients, r'Molar Absorption Coefficient $\mathdefault{(L mole^{-1} cm^{-1})}$')

    def absorptionButtonClicked(self):
        debugger.print('absorptionButtonClicked pressed')
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        if not self.notebook.singleCrystalCalculationRequired:
            self.plot(self.xaxes, self.absorptionCoefficients, r'Absorption Coefficient $\mathdefault{(cm^{-1})}$')

    def reflectivityButtonClicked(self):
        debugger.print('reflectivityButtonClicked pressed')
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        if not self.notebook.singleCrystalCalculationRequired:
            self.plot(self.xaxes, self.reflectivities, 'Reflectivities')

    def on_title_changed(self,text):
        self.settings['Plot title'] = text
        if self.subplot is not None:
            self.subplot.set_title(self.settings['Plot title'])
            self.canvas.draw_idle()
        debugger.print('on title change ', self.settings['Plot title'])

    def on_vinc_changed(self,value):
        self.settings['Frequency increment'] = value
        self.notebook.singleCrystalCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        debugger.print('on vinc change ', self.settings['Frequency increment'])

    def on_vmin_changed(self):
        self.settings['Minimum frequency'] = self.vmin_sb.value()
        self.notebook.singleCrystalCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        debugger.print('on vmin change ', self.settings['Minimum frequency'])

    def on_vmax_changed(self):
        self.settings['Maximum frequency'] = self.vmax_sb.value()
        self.notebook.singleCrystalCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        debugger.print('on vmax change ', self.settings['Maximum frequency'])

    def refresh(self,force=False):
        if not self.dirty and not force and not self.notebook.singleCrystalCalculationRequired:
            debugger.print('refreshing widget aborted', self.dirty,force,self.notebook.singleCrystalCalculationRequired)
            return
        debugger.print('refreshing widget', force)
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Now refresh values
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.vinc_sb.setValue(self.settings['Frequency increment'])
        try:
            self.molar_cb_current_index = self.molar_definitions.index(self.settings["Molar definition"])
        except:
            self.molar_cb_current_index = 0
            self.settings["Molar definition"] = self.molar_definitions[self.molar_cb_current_index]
        self.molar_cb.setCurrentIndex(self.molar_cb_current_index)
        self.natoms_sb.setValue(self.settings['Number of atoms'])
        self.title_le.setText(self.settings['Plot title'])
        # Refresh the widgets that depend on the reader
        self.reader = self.notebook.reader
        if self.reader is not None:
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
        # Flag a recalculation will be required
        self.notebook.singleCrystalCalculationRequired = True
        # Reset the progress bar
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        self.molarAbsorptionButtonClicked()
        self.dirty = False
        self.notebook.singleCrystalCalculationRequired = False
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        return

    def on_natoms_changed(self, value):
        self.notebook.singleCrystalCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        self.settings['Number of atoms'] = value
        debugger.print('on natoms changed ', self.settings['Number of atoms'])

    def on_funits_cb_activated(self, index):
        if index == 0:
            self.frequency_units = 'wavenumber'
        else:
            self.frequency_units = 'THz'
        self.replot()
        debugger.print('Frequency units changed to ', self.frequency_units)

    def on_molar_cb_activated(self, index):
        self.molar_cb_current_index = index
        self.settings['Molar definition'] = self.molar_definitions[index]
        self.notebook.singleCrystalCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        if self.settings['Molar definition'] == 'Molecules':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 * self.settings['Number of atoms'] / self.reader.nions)
            self.natoms_sb.setEnabled(True)
        elif self.settings['Molar definition'] == 'Unit cells':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
            self.natoms_sb.setEnabled(False)
        elif self.settings['Molar definition'] == 'Atoms':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 / self.reader.nions)
            self.natoms_sb.setEnabled(False)
        debugger.print('The concentration has been set', self.settings['Molar definition'], self.settings['concentration'])

    def calculate(self):
        debugger.print('Calculate')
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['Program']
        filename = settings['Output file name']
        reader = self.notebook.mainTab.reader
        if reader is None:
            return
        if program == '':
            return
        if filename == '':
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Assemble the settingsTab settings
        self.notebook.singleCrystalCalculationRequired = False
        settings = self.notebook.settingsTab.settings
        epsilon_inf = np.array(settings['Optical permittivity'])
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        frequencies = np.array(frequencies_cm1) * wavenumber
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = reader.volume*angstrom*angstrom*angstrom
        vmin = self.settings['Minimum frequency']
        vmax = self.settings['Maximum frequency']
        vinc = self.settings['Frequency increment']
        self.scenarios = self.notebook.scenarios
        mode_list = []
        drude_sigma = 0
        drude_plasma = 0
        drude = False
        for mode_index,selected in enumerate(modes_selected):
            if selected:
                mode_list.append(mode_index)
        #
        # Calculate the ionic permittivity at zero frequency
        #
        cell = reader.unit_cells[-1]
        self.legends = []
        # Set up parameters for the parallel execution
        call_parameters = []
        for v in np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc)):
            vau = v * wavenumber
            call_parameters.append( (v, vau, mode_list, frequencies, sigmas, oscillator_strengths,
                                     volume, epsilon_inf, drude, drude_plasma, drude_sigma) )

        # Initialise the progress bare
        maximum_progress = len(call_parameters) * (1 + len(self.scenarios))
        progress = 0
        self.progressbar.setMaximum(maximum_progress)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setMaximum(maximum_progress)
        QCoreApplication.processEvents()
        number_of_processors = self.notebook.ncpus
        # Switch off mkl threading
        try:
            import mkl
            mkl.set_num_threads(1)
        except:
            pass
        #jk start = time.time()
        if self.notebook.threading:
            from multiprocess.dummy import Pool
            pool = Pool(number_of_processors)
        else:
            from multiprocess import Pool
            pool = Pool(number_of_processors, initializer=set_affinity_on_worker, maxtasksperchild=10)
        # start parallel execution of dielectric constant
        dielecv_results = []
        for result in pool.imap(Calculator.parallel_dielectric, call_parameters,chunksize=40):
            dielecv_results.append(result)
            progress += 1
            self.progressbar.setValue(progress)
            if self.notebook.progressbar is not None:
                self.notebook.progressbar.setValue(progress)
        pool.close()
        pool.join()
        QCoreApplication.processEvents()
        #jk print('Dielec calculation duration ', time.time()-start)
        nplots = len(dielecv_results)
        # Prepare variables for setting up parameters for parallel call
        concentration = self.settings['concentration']
        hkl = [ self.settings['Unique direction - h'] , self.settings['Unique direction - k'], self.settings['Unique direction - l'] ]
        rotz                  = self.settings['Rotation about Z']
        angleOfIncidence      = self.settings['Angle of incidence']
        polarisationAngle     = self.settings['Polarisation angle']
        superstrateDepth      = self.settings['Superstrate depth']
        superstrateDielectric = self.settings['Superstrate dielectric']
        substrateDepth        = self.settings['Substrate depth']
        substrateDielecitric  = self.settings['Substrate dielectric']
        crystalDepth          = self.settings['Crystal depth']
        # Prepare parallel processing options
        if self.notebook.threading:
            from multiprocess.dummy import Pool
            pool = Pool(number_of_processors)
        else:
            from multiprocess import Pool
            pool = Pool(number_of_processors, initializer=set_affinity_on_worker, maxtasksperchild=10)
        #
        # Prepare parameters for a parallel call to the layered absorption / reflection
        #
        nplot = 0
        # Assemble all the parameters we need for parallel execution
        for v,vau,dielecv in dielecv_results:
            data = ''
            call_parameters.append( (v,vau,dielecv,
                                     rotz                  ,
                                     angleOfIncidence      ,
                                     polarisationAngle     ,
                                     superstrateDepth      ,
                                     superstrateDielectric ,
                                     substrateDepth        ,
                                     substrateDielecitric  ,
                                     crystalDepth
                                     ) )
            nplot += 1
        results = []
        for result in pool.imap(Calculator.solve_singleCrystalEquations, call_parameters, chunksize=40):
            results.append(result)
            progress += 1
            self.progressbar.setValue(progress)
            if self.notebook.progressbar is not None:
                self.notebook.progressbar.setValue(progress)
        QCoreApplication.processEvents()
        # Initialise plotting variables
        self.xaxes                       = []
        self.reflectivities              = []
        self.absorptionCoefficients      = []
        self.molarAbsorptionCoefficients = []
        for v,nplot,absorption_coefficient,molar_absorption_coefficient,reflectance in results:
             self.xaxis.append(v)
             self.reflectivities.append(np.real(trace))
             self.absorptionCoefficient.append(absorption_coefficient)
             self.molarAbsorptionCoefficient.append(molar_absorption_coefficient)
        # Close parallel processing down
        pool.close()
        pool.join()
        #
        # Switch on mkl threading
        #
        try:
            import mkl
            mkl.set_num_threads(number_of_processors)
        except:
            pass
        #jk print('Dielec calculation duration ', time.time()-start)
        #if self.notebook.spreadsheet is not None:
        #    self.write_spreadsheet()
        self.dirty = False
        QApplication.restoreOverrideCursor()
        QCoreApplication.processEvents()

    def write_spreadsheet(self):
        if self.notebook.spreadsheet is None:
            return
        sp = self.notebook.spreadsheet
        sp.selectWorkSheet('Single Crystal')
        sp.delete()
        sp.writeNextRow(['Settings for the single crystal calculation of absorption and reflection'],col=1)
        sp.writeNextRow([''],col=1)
        sp.writeNextRow([ 'Minimum frequency',            self.settings['Minimum frequency'] ])
        sp.writeNextRow([ 'Maximum frequency',            self.settings['Maximum frequency'] ])
        sp.writeNextRow([ 'Frequency increment',          self.settings['Frequency increment'] ])
        sp.writeNextRow([ 'Molar definition',             self.settings['Molar definition'] ])
        sp.writeNextRow([ 'Number of atoms in a molecule',self.settings['Number of atoms'] ])
        sp.writeNextRow([ 'Plot title',                   self.settings['Plot title'] ])
        sp.writeNextRow([ 'Surface definition (h)',       self.settings['Unique direction - h'] ])
        sp.writeNextRow([ 'Surface definition (k)',       self.settings['Unique direction - k'] ])
        sp.writeNextRow([ 'Surface definition (.)',       self.settings['Unique direction - l'] ])
        sp.writeNextRow([ 'Rotation angle about Z',       self.settings['Rotation about Z'] ])
        sp.writeNextRow([ 'Angle of incidence',           self.settings['Angle of incidence'] ])
        sp.writeNextRow([ 'Polarisation angle',           self.settings['Polarisation angle'] ])
        sp.writeNextRow([ 'Superstrate depth(mm)',        self.settings['Superstrate depth'] ])
        sp.writeNextRow([ 'Superstrate dielectric',       self.settings['Superstrate dielectric'] ])
        sp.writeNextRow([ 'Substrate depth(mm)',          self.settings['Substrate depth'] ])
        sp.writeNextRow([ 'Substrate dielectric',         self.settings['Substrate dielectric'] ])
        sp.writeNextRow([ 'Crystal depth(mm)',            self.settings['Crystal depth'] ])
        # Now deal with Molar absorption, absorption, reflectivities
        self.write_results(sp, 'Crystal Molar Absorption',self.xaxes, self.molarAbsorptionCoefficients)
        self.write_results(sp, 'Crystal Absorption',      self.xaxes, self.absorptionCoefficients)
        self.write_results(sp, 'Crystal Reflectance',     self.xaxes, self.reflectivities)

    def write_results(self, sp, name, vss, yss):
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ['frequencies (cm-1)']
        #for isc,ys in enumerate(yss):
        #    headers.append('Scenario'+str(isc))
        headers.extend(self.legends)
        sp.writeNextRow(headers,row=0, col=1)
        for iv,v in enumerate(vss[0]):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)

    def plot(self,xs,ys,ylabel):
        # import matplotlib.pyplot as pl
        # mp.use('Qt5Agg')
        self.subplot = None
        self.remember_xs = xs
        self.remember_ys = ys
        self.remember_ylabel = ylabel
        self.figure.clf()
        if self.frequency_units == 'wavenumber':
            xlabel = r'Frequency $\mathdefault{(cm^{-1})}}$'
            scale = 1.0
        else:
            xlabel = r'THz'
            scale = 0.02998
        self.subplot = self.figure.add_subplot(111)
        for scenario,x,y in zip(self.scenarios,xs,ys):
            x = np.array(x)
            legend = scenario.settings['Legend']
            line, = self.subplot.plot(scale*x,y,lw=2, label=legend )
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        self.subplot.legend(loc='best')
        self.subplot.set_title(self.settings['Plot title'])
        self.canvas.draw_idle()

    def replot(self):
        if self.subplot is not None:
            self.plot(self.remember_xs, self.remember_ys, self.remember_ylabel)


