import os.path
import os
import numpy as np
import PDielec.Calculator         as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QProgressBar, QApplication
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox,QDoubleSpinBox
from PyQt5.QtWidgets  import  QSizePolicy
from PyQt5.QtCore     import  QCoreApplication, Qt
from PDielec.DielectricFunction import DielectricFunction
from PDielec.Constants import  wavenumber, PI, avogadro_si, angstrom, speed_light_si
import ctypes
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PDielec.Utilities import Debug
import PDielec.GTMcore as GTM

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
        self.settings['Minimum frequency'] = 100
        self.settings['Maximum frequency'] = 400
        self.settings['Frequency increment'] = 0.2
        self.molar_definitions = ['Unit cells','Atoms','Molecules']
        self.settings['Molar definition'] = 'Unit cells'
        self.settings['Number of atoms'] = 1
        self.settings['Plot title'] = 'Plot Title'
        self.settings['Unique direction - h'] = 0
        self.settings['Unique direction - k'] = 0
        self.settings['Unique direction - l'] = 1
        self.settings['Azimuthal angle'] = 0.0
        self.settings['Angle of incidence'] = 45.0
        self.settings['Superstrate dielectric'] = 1.0
        self.settings['Substrate dielectric'] = 1.0
        self.settings['Superstrate depth'] = 999.0
        self.settings['Substrate depth'] = 999.0
        self.settings['Crystal depth'] = 10.0
        self.xaxis = []
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
        self.vmin_sb.setRange(1,9000)
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
        #jk #
        #jk # Define molar quantity
        #jk #
        #jk self.molar_cb = QComboBox(self)
        #jk self.molar_cb.setToolTip('Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given')
        #jk self.molar_cb.addItems(self.molar_definitions)
        #jk try:
        #jk     self.molar_cb_current_index = self.molar_definitions.index(self.settings["Molar definition"])
        #jk except:
        #jk     self.molar_cb_current_index = 0
        #jk     self.settings["Molar definition"] = self.molar_definitions[self.molar_cb_current_index]
        #jk self.molar_cb.setCurrentIndex(self.molar_cb_current_index)
        #jk self.molar_cb.activated.connect(self.on_molar_cb_activated)
        #jk label = QLabel('Molar definition', self)
        #jk label.setToolTip('Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given')
        #jk form.addRow(label, self.molar_cb)
        #jk #
        #jk # Number of atoms in a molecule
        #jk #
        #jk self.natoms_sb = QSpinBox(self)
        #jk self.natoms_sb.setToolTip('Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed')
        #jk self.natoms_sb.setRange(1,500)
        #jk self.natoms_sb.setValue(self.settings['Number of atoms'])
        #jk self.natoms_sb.valueChanged.connect(self.on_natoms_changed)
        #jk self.natoms_sb.setEnabled(False)
        #jk label = QLabel('Number of atoms per molecule', self)
        #jk label.setToolTip('Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed')
        #jk form.addRow(label, self.natoms_sb)
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
        self.azimuthal_angle_sb = QDoubleSpinBox(self)
        self.azimuthal_angle_sb.setToolTip('Define the slab azimuthal angle')
        self.azimuthal_angle_sb.setRange(0,360)
        self.azimuthal_angle_sb.setSingleStep(10)
        self.azimuthal_angle_sb.setValue(self.settings['Azimuthal angle'])
        self.azimuthal_angle_sb.valueChanged.connect(self.on_azimuthal_angle_sb_changed)
        self.angle_of_incidence_sb = QDoubleSpinBox(self)
        self.angle_of_incidence_sb.setToolTip('Define the angle of incidence')
        self.angle_of_incidence_sb.setRange(0,180)
        self.angle_of_incidence_sb.setSingleStep(10)
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        self.angle_of_incidence_sb.valueChanged.connect(self.on_angle_of_incidence_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.azimuthal_angle_sb)
        hbox.addWidget(self.angle_of_incidence_sb)
        self.angles_label = QLabel('Angles of surface orientation and incidence',self)
        self.angles_label.setToolTip('Define the azimuthal angle and the angle of incidence.')
        form.addRow(self.angles_label, hbox)
        #
        # Define the superstrate and substrate dielectrics
        #
        self.superstrate_dielectric_sb = QDoubleSpinBox(self)
        self.superstrate_dielectric_sb.setToolTip('Define the superstrate dielectric')
        self.superstrate_dielectric_sb.setRange(1,1000)
        self.superstrate_dielectric_sb.setSingleStep(0.1)
        self.superstrate_dielectric_sb.setValue(self.settings['Superstrate dielectric'])
        self.superstrate_dielectric_sb.valueChanged.connect(self.on_superstrate_dielectric_sb_changed)
        self.substrate_dielectric_sb = QDoubleSpinBox(self)
        self.substrate_dielectric_sb.setToolTip('Define the substrate dielectric')
        self.substrate_dielectric_sb.setRange(1,1000)
        self.substrate_dielectric_sb.setSingleStep(0.1)
        self.substrate_dielectric_sb.setValue(self.settings['Substrate dielectric'])
        self.substrate_dielectric_sb.valueChanged.connect(self.on_substrate_dielectric_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.superstrate_dielectric_sb)
        hbox.addWidget(self.substrate_dielectric_sb)
        self.superstrate_label = QLabel('Super- and substrate dielectrics',self)
        self.superstrate_label.setToolTip('Define the dielectric constants of the superstrate and substrates')
        form.addRow(self.superstrate_label, hbox)
        #
        #
        # Define the crystal depth
        #
        self.crystal_depth_sb = QDoubleSpinBox(self)
        self.crystal_depth_sb.setToolTip('Define the crystal depth')
        self.crystal_depth_sb.setRange(0,10000)
        self.crystal_depth_sb.setSingleStep(1)
        self.crystal_depth_sb.setValue(self.settings['Crystal depth'])
        self.crystal_depth_sb.valueChanged.connect(self.on_crystal_depth_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.crystal_depth_sb)
        self.crystal_depth_label = QLabel('Crystal depth(mu)',self)
        self.crystal_depth_label.setToolTip('Define the depth of the crystal in millimetres.')
        form.addRow(self.crystal_depth_label, hbox)
        #
        #
        # Final button
        #
        hbox = QHBoxLayout()
        self.plotButton = QPushButton('Plot transmission and reflection')
        self.plotButton.setToolTip('Plot the molar absorption')
        self.plotButton.clicked.connect(self.plotButtonClicked)
        hbox.addWidget(self.plotButton)
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

    def on_crystal_depth_sb_changed(self,value):
        debugger.print('on_crystal_depth_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Crystal depth'] = value

    def on_superstrate_dielectric_sb_changed(self,value):
        debugger.print('on_superstrate_dielectric_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Superstrate dielectric'] = value

    def on_substrate_dielectric_sb_changed(self,value):
        debugger.print('on_substrate_dielectric_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Substrate dielectric'] = value

    def on_azimuthal_angle_sb_changed(self,value):
        debugger.print('on_azimuthal_angl_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Rotation about Z'] = value

    def on_angle_of_incidence_sb_changed(self,value):
        debugger.print('on_angle_of_incidence_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Angle of incidence'] = value

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


    def plotButtonClicked(self):
        debugger.print('plotButtonClicked pressed')
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        if not self.notebook.singleCrystalCalculationRequired:
            yaxes = [self.p_transmission, self.s_transmission, self.p_reflectivity, self.s_reflectivity, ]
            legends = [ r'$T_p$', r'$T_s$', r'$R_p$', r'$R_s$']
            ylabel = 'Arbitrary units'
            self.plot(self.xaxis, yaxes, ylabel, legends)

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
        self.plotButtonClicked()
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
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        frequencies = np.array(frequencies_cm1) * wavenumber
        # Get the last unit cell in the reader
        cell = reader.unit_cells[-1]
        # Calculate frequency range
        vmin = self.settings['Minimum frequency']
        vmax = self.settings['Maximum frequency']
        vinc = self.settings['Frequency increment']
        vs = np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc))
        # Set up the layered system using GTM calls
        S = GTM.System()
        # The dielectric variables are functions of frequency
        superstrateDielectric = self.settings['Superstrate dielectric']
        substrateDielectric   = self.settings['Substrate dielectric']
        superstrateDielectricFunction = DielectricFunction(epsType='constant',units='hz',parameters=superstrateDielectric).function()
        substrateDielectricFunction   = DielectricFunction(epsType='constant',units='hz',parameters=substrateDielectric).function()
        # The crystal dielectric has already been defined in the SettingsTab
        # Make sure the system knows that frequency will be supplied using Hz
        self.notebook.settingsTab.CrystalPermittivity.setUnits('hz')
        crystalPermittivityFunction     = self.notebook.settingsTab.CrystalPermittivity.function()
        # Create 3 layers, thickness is converted from microns to metres
        superstrateDepth = self.settings['Superstrate depth']
        superstrate      = GTM.Layer(thickness=superstrateDepth*1e-6,epsilon1=superstrateDielectricFunction)
        substrateDepth   = self.settings['Substrate depth']
        substrate        = GTM.Layer(thickness=substrateDepth*1e-6,  epsilon1=substrateDielectricFunction)
        crystalDepth     = self.settings['Crystal depth']
        crystal          = GTM.Layer(thickness=crystalDepth*1e-6,    epsilon=crystalPermittivityFunction)
        # Creat the system with the layers 
        system = GTM.System(substrate=substrate, superstrate=superstrate, layers=[crystal])
        # Prepare variables for setting up parameters for parallel call
        concentration = self.settings['concentration']
        hkl = [ self.settings['Unique direction - h'] , self.settings['Unique direction - k'], self.settings['Unique direction - l'] ]
        # convert normal to plane to a direction in xyz coordinates
        normal_to_plane_xyz = cell.convert_hkl_to_xyz(hkl)
        normal_to_plane_xyz /=  np.linalg.norm(normal_to_plane_xyz)
        # Rotate this normal so it becomes the laboratory Z- direction
        labZ = np.array( [0,0,1] )
        # Define the euler angles in radians; theta rotates z to Z and psi is the azimuthal angle
        euler_theta           = np.arccos(np.dot(labZ,normal_to_plane_xyz))
        euler_phi             = 0.0
        euler_psi             = np.pi / 180.0 * self.settings['Azimuthal angle']
        # Tricky setting for 90 degrees
        angle = self.settings['Angle of incidence']
        if angle >= 89.99:
            angle = 89.99
        angleOfIncidence      = np.pi / 180.0 * angle
        # Rotate the dielectric constants to the laboratory frame
        system.substrate.set_euler(theta=euler_theta, phi=euler_phi, psi=euler_psi)
        system.superstrate.set_euler(theta=euler_theta, phi=euler_phi, psi=euler_psi)
        for layer in system.layers:
            layer.set_euler(theta=euler_theta, phi=euler_phi, psi=euler_psi)
        # Prepare parallel processing options
        number_of_processors = self.notebook.ncpus
        if self.notebook.threading:
            from multiprocess.dummy import Pool
            pool = Pool(number_of_processors)
        else:
            from multiprocess import Pool
            pool = Pool(number_of_processors, initializer=set_affinity_on_worker, maxtasksperchild=10)
        #
        # Prepare parameters for a parallel call to the layered absorption / reflection
        #
        call_parameters = []
        # Assemble all the parameters we need for parallel execution
        # About to call
        results = []
        progress = 0
        for v in vs:
            data = ''
            call_parameters.append( (v, angleOfIncidence, system) )
            # results.append( Calculator.solve_single_crystal_equations( (v, angleOfIncidence, system) ) )
        for result in pool.imap(Calculator.solve_single_crystal_equations, call_parameters, chunksize=40):
            results.append(result)
            progress += 1
            self.progressbar.setValue(progress)
            if self.notebook.progressbar is not None:
                self.notebook.progressbar.setValue(progress)
        QCoreApplication.processEvents()
        # Initialise plotting variables
        self.xaxis          = []
        self.p_reflectivity = []
        self.s_reflectivity = []
        self.p_transmission = []
        self.s_transmission = []
        for v,r,R,t,T in results:
            self.xaxis.append(v)
            self.p_reflectivity.append(R[0])
            self.s_reflectivity.append(R[1])
            self.p_transmission.append(T[0])
            self.s_transmission.append(T[1])
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
        sp.writeNextRow([ 'Minimum frequency',            self.settings['Minimum frequency'] ],col=1)
        sp.writeNextRow([ 'Maximum frequency',            self.settings['Maximum frequency'] ],col=1)
        sp.writeNextRow([ 'Frequency increment',          self.settings['Frequency increment'] ],col=1)
        sp.writeNextRow([ 'Molar definition',             self.settings['Molar definition'] ],col=1)
        sp.writeNextRow([ 'Number of atoms in a molecule',self.settings['Number of atoms'] ],col=1)
        sp.writeNextRow([ 'Plot title',                   self.settings['Plot title'] ],col=1)
        sp.writeNextRow([ 'Surface definition (h)',       self.settings['Unique direction - h'] ],col=1)
        sp.writeNextRow([ 'Surface definition (k)',       self.settings['Unique direction - k'] ],col=1)
        sp.writeNextRow([ 'Surface definition (l)',       self.settings['Unique direction - l'] ],col=1)
        sp.writeNextRow([ 'Azimuthal angle',              self.settings['Azimuthal angle'] ],col=1)
        sp.writeNextRow([ 'Angle of incidence',           self.settings['Angle of incidence'] ],col=1)
        sp.writeNextRow([ 'Superstrate dielectric',       self.settings['Superstrate dielectric'] ],col=1)
        sp.writeNextRow([ 'Substrate dielectric',         self.settings['Substrate dielectric'] ],col=1)
        sp.writeNextRow([ 'Crystal depth(mu)',            self.settings['Crystal depth'] ],col=1)
        # Now deal with Molar absorption, absorption, reflectivities
        headings = ['R_p', 'R_s', 'T_p', 'T_s']
        self.write_results(sp, 'Crystal R&T',     self.xaxis, [self.p_reflectivity, self.s_reflectivity, self.p_transmission, self.s_transmission], headings)

    def write_results(self, sp, name, vs, yss, headings):
        """ 
        sp        is the spreadsheet object
        name      is the worksheet name used for writing
        vs        an np.array of the frequencies
        yss       a list of np.arrays of the reflections and transmissions ] 
        headings  the heading names for the yss
        """
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ['frequencies (cm-1)']
        headers.extend(headings)
        sp.writeNextRow(headers,row=0, col=1)
        for iv,v in enumerate(vs):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)

    def plot(self,x,ys,ylabel,legends):
        # import matplotlib.pyplot as pl
        # mp.use('Qt5Agg')
        self.subplot = None
        self.remember_x  = x
        self.remember_ys = ys
        self.remember_ylabel = ylabel
        self.remember_legends = legends
        self.figure.clf()
        if self.frequency_units == 'wavenumber':
            xlabel = r'Frequency $\mathdefault{(cm^{-1})}}$'
            scale = 1.0
        else:
            xlabel = r'THz'
            scale = 0.02998
        x = np.array(x)
        self.subplot = self.figure.add_subplot(111)
        for y,legend in zip(ys,legends):
            line, = self.subplot.plot(scale*x,y,lw=2, label=legend )
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        self.subplot.legend(loc='best')
        self.subplot.set_title(self.settings['Plot title'])
        self.canvas.draw_idle()

    def replot(self):
        if self.subplot is not None:
            self.plot(self.remember_x, self.remember_ys, self.remember_ylabel, self.remember_legends)


