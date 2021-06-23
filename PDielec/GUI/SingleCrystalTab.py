import os.path
import os
import numpy as np
import copy
import PDielec.Calculator         as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit, QListWidget
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
        self.settings['Unique direction - h'] = 0
        self.settings['Unique direction - k'] = 0
        self.settings['Unique direction - l'] = 1
        self.settings['Azimuthal angle'] = 0.0
        self.settings['Angle of incidence'] = 0.0
        self.settings['Superstrate dielectric'] = 1.0
        self.settings['Substrate dielectric'] = 1.0
        self.settings['Superstrate depth'] = 999.0
        self.settings['Substrate depth'] = 999.0
        self.settings['Film thickness'] = 10.0
        self.settings['Mode'] = 'Thick slab'
        self.settings['Plot title'] = 'Single crystal '+self.settings['Mode']
        self.settings['Use default plot title'] = True
        self.settings['Frequency units'] = 'wavenumber'
        self.p_reflectivity = []
        self.s_reflectivity = []
        self.p_transmission = []
        self.s_transmission = []
        self.epsilon = []
        self.plotLastButtonPressed = self.plotReflectanceButtonClicked
        self.frequencies_cm1 = []
        self.directions = []
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Get the last unit cell in the reader
        if self.reader is not None:
            self.cell = self.reader.unit_cells[-1]
        # Create last tab - SingleCrystalTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # Chose mode of operation
        #
        self.mode_cb = QComboBox(self)
        self.mode_cb.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two media; the incident and the crystal),\n Coherent thin film and incoherent thin film assumes there are three media; the incident, the crystal and the substrate')
        self.mode_cb.addItems( ['Thick slab','Coherent thin film','Incoherent thick film'] )
        self.settings['Mode'] = 'Thick slab'
        self.mode_cb.activated.connect(self.on_mode_cb_activated)
        label = QLabel('Single crystal mode', self)
        label.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two media; the incident and the crystal),\n Coherent thin film iand incoherent thin film assumes there are three media; the incident, the crystal and the substrate')
        form.addRow(label, self.mode_cb)
        #
        # The minimum frequency
        #
        self.vmin_sb = QSpinBox(self)
        self.vmin_sb.setRange(1,9000)
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmin_sb.setToolTip('Set the minimum frequency to be considered)')
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        # form.addRow(label, self.vmin_sb)
        #
        # The maximum frequency
        #
        self.vmax_sb = QSpinBox(self)
        self.vmax_sb.setRange(0,9000)
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.vmax_sb.setToolTip('Set the maximum frequency to be considered)')
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        # form.addRow(label, self.vmax_sb)
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
        #
        label = QLabel('Frequency min, max and increment', self)
        label.setToolTip('Choose minimum, maximum and increment frequencies when plotting')
        #
        hbox = QHBoxLayout()
        hbox.addWidget(self.vmin_sb)
        hbox.addWidget(self.vmax_sb)
        hbox.addWidget(self.vinc_sb)
        form.addRow(label, hbox)
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
        self.settings['Frequency units'] = 'wavenumber'
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
        hkl_label = QLabel('Crystal surface (hkl)',self)
        hkl_label.setToolTip('Define the crystal surface (hkl). Defines the unique direction in crystallographic units.')
        form.addRow(hkl_label, hbox)
        #
        # Define the rotation angle of the slab and the angle of incidence
        #
        self.azimuthal_angle_sb = QDoubleSpinBox(self)
        self.azimuthal_angle_sb.setToolTip('Define the slab azimuthal angle')
        self.azimuthal_angle_sb.setRange(-180,360)
        self.azimuthal_angle_sb.setSingleStep(10)
        self.azimuthal_angle_sb.setValue(self.settings['Azimuthal angle'])
        self.azimuthal_angle_sb.valueChanged.connect(self.on_azimuthal_angle_sb_changed)
        self.angle_of_incidence_sb = QDoubleSpinBox(self)
        self.angle_of_incidence_sb.setToolTip('Define the angle of incidence, (normal incidence is 0 degrees)')
        self.angle_of_incidence_sb.setRange(0,90)
        self.angle_of_incidence_sb.setSingleStep(5)
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        self.angle_of_incidence_sb.valueChanged.connect(self.on_angle_of_incidence_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.azimuthal_angle_sb)
        hbox.addWidget(self.angle_of_incidence_sb)
        angles_label = QLabel('Angles of surface orientation and incidence',self)
        angles_label.setToolTip('Define the azimuthal angle and the angle of incidence, (normal incidence is 0 degrees).')
        form.addRow(angles_label, hbox)
        #
        # Provide information on the lab frame
        #
        labframe_l = QLabel('Lab frame information', self)
        labframe_l.setToolTip('The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y')
        self.labframe_w = QListWidget(self)
        fm = self.labframe_w.fontMetrics()
        h = fm.ascent() + fm.descent()
        self.labframe_w.setMaximumHeight(6*h)
        self.labframe_w.setToolTip('The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y')
        form.addRow(labframe_l,self.labframe_w)
        #
        # Define the superstrate and substrate dielectrics
        #
        self.superstrate_dielectric_sb = QDoubleSpinBox(self)
        self.superstrate_dielectric_sb.setToolTip('Define the incident medium permittivity')
        self.superstrate_dielectric_sb.setRange(1,1000)
        self.superstrate_dielectric_sb.setSingleStep(0.1)
        self.superstrate_dielectric_sb.setValue(self.settings['Superstrate dielectric'])
        self.superstrate_dielectric_sb.valueChanged.connect(self.on_superstrate_dielectric_sb_changed)
        self.substrate_dielectric_sb = QDoubleSpinBox(self)
        self.substrate_dielectric_sb.setToolTip('Define the substrate permittivity')
        self.substrate_dielectric_sb.setRange(1,1000)
        self.substrate_dielectric_sb.setSingleStep(0.1)
        self.substrate_dielectric_sb.setValue(self.settings['Substrate dielectric'])
        self.substrate_dielectric_sb.valueChanged.connect(self.on_substrate_dielectric_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.superstrate_dielectric_sb)
        hbox.addWidget(self.substrate_dielectric_sb)
        superstrate_label = QLabel('Incident medium and substrate permittivities',self)
        superstrate_label.setToolTip('Define the permttivites of the incident medium and substrate')
        form.addRow(superstrate_label, hbox)
        #
        #
        # Define the Film thickness
        #
        self.film_thickness_sb = QDoubleSpinBox(self)
        self.film_thickness_sb.setToolTip('Define the Film thickness')
        self.film_thickness_sb.setRange(0,10000)
        self.film_thickness_sb.setSingleStep(1)
        self.film_thickness_sb.setValue(self.settings['Film thickness'])
        self.film_thickness_sb.valueChanged.connect(self.on_film_thickness_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.film_thickness_sb)
        film_thickness_label = QLabel('Film thickness(mu)',self)
        film_thickness_label.setToolTip('Define the depth of the crystal in millimetres.')
        form.addRow(film_thickness_label, hbox)
        #
        #
        # Final button
        #
        hbox = QHBoxLayout()
        self.plotTransmissionButton = QPushButton('Plot T')
        self.plotTransmissionButton.setToolTip('Plot transmittance')
        self.plotTransmissionButton.clicked.connect(self.plotTransmissionButtonClicked)
        self.plotReflectanceButton = QPushButton('Plot R')
        self.plotReflectanceButton.setToolTip('Plot reflectance')
        self.plotReflectanceButton.clicked.connect(self.plotReflectanceButtonClicked)
        self.plotTplusRButton = QPushButton('Plot T+R')
        self.plotTplusRButton.setToolTip('Plot sum of transmittance and reflectance')
        self.plotTplusRButton.clicked.connect(self.plotTplusRButtonClicked)
        self.plotRealEButton = QPushButton('Plot real eps')
        self.plotRealEButton.setToolTip('Plot real permittivity')
        self.plotRealEButton.clicked.connect(self.plotRealEButtonClicked)
        self.plotImagEButton = QPushButton('Plot imag eps')
        self.plotImagEButton.setToolTip('Plot imaginary permittivity')
        self.plotImagEButton.clicked.connect(self.plotImagEButtonClicked)
        hbox.addWidget(self.plotTransmissionButton)
        hbox.addWidget(self.plotReflectanceButton)
        hbox.addWidget(self.plotTplusRButton)
        hbox.addWidget(self.plotRealEButton)
        hbox.addWidget(self.plotImagEButton)
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

    def on_film_thickness_sb_changed(self,value):
        debugger.print('on_film_thickness_sb', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Film thickness'] = value

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
        self.settings['Azimuthal angle'] = value
        self.calculate_euler_angles()

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
        self.calculate_euler_angles()

    def on_k_sb_changed(self,value):
        debugger.print('on_k_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Unique direction - k'] = value
        self.calculate_euler_angles()

    def on_l_sb_changed(self,value):
        debugger.print('on_l_sb_changed', value)
        self.dirty = True
        self.notebook.singleCrystalCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        self.settings['Unique direction - l'] = value
        self.calculate_euler_angles()

    def plotTransmissionButtonClicked(self):
        debugger.print('plotTransmissionButtonClicked pressed')
        self.plotLastButtonPressed = self.plotTransmissionButtonClicked
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        yaxes = [self.p_transmission, self.s_transmission, ]
        legends = [ r'$T_p$', r'$T_s$']
        ylabel = 'Fraction of incident light'
        self.plot(self.frequencies_cm1, yaxes, ylabel, legends, 'Transmitance')

    def plotReflectanceButtonClicked(self):
        debugger.print('plotReflectanceButtonClicked pressed')
        self.plotLastButtonPressed = self.plotReflectanceButtonClicked
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        yaxes = [self.p_reflectivity, self.s_reflectivity, ]
        legends = [ r'$R_p$', r'$R_s$']
        ylabel = 'Fraction of incident light'
        self.plot(self.frequencies_cm1, yaxes, ylabel, legends, 'Reflectance')

    def plotTplusRButtonClicked(self):
        debugger.print('plotTplusRButtonClicked pressed')
        self.plotLastButtonPressed = self.plotTplusRButtonClicked
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        tr_p = np.array( self.p_transmission ) + np.array( self.p_reflectivity )
        tr_s = np.array( self.s_transmission ) + np.array( self.s_reflectivity )
        yaxes = [tr_p,tr_s, ]
        legends = [ r'$T_p+R_p$', r'$T_s+R_s$']
        ylabel = 'Fraction of incident light'
        self.plot(self.frequencies_cm1, yaxes, ylabel, legends, 'Transmitance+Reflectance')

    def plotRealEButtonClicked(self):
        debugger.print('plotRealEButtonClicked pressed')
        self.plotLastButtonPressed = self.plotRealEButtonClicked
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        eps = np.real( self.epsilon )
        yaxes = [eps[:,0,0], eps[:,1,1], eps[:,2,2], eps[:,0,1], eps[:,0,2], eps[:,1,2] ]
        legends = [ r'$\epsilon_{xx}$', r'$\epsilon_{yy}$', r'$\epsilon_{zz}$', r'$\epsilon_{xy}$',r'$\epsilon_{xz}$',r'$\epsilon_{yz}$' ]
        ylabel = 'Real permittivity'
        self.plot(self.frequencies_cm1, yaxes, ylabel, legends, 'Real permittivity')

    def plotImagEButtonClicked(self):
        debugger.print('plotImagEButtonClicked pressed')
        self.plotLastButtonPressed = self.plotImagEButtonClicked
        if self.notebook.singleCrystalCalculationRequired:
            self.calculate()
        eps = np.imag( self.epsilon )
        yaxes = [eps[:,0,0], eps[:,1,1], eps[:,2,2], eps[:,0,1], eps[:,0,2], eps[:,1,2] ]
        legends = [ r'$\epsilon_{xx}$', r'$\epsilon_{yy}$', r'$\epsilon_{zz}$', r'$\epsilon_{xy}$',r'$\epsilon_{xz}$',r'$\epsilon_{yz}$' ]
        ylabel = 'Imaginary permittivity'
        self.plot(self.frequencies_cm1, yaxes, ylabel, legends, 'Imaginary permittivity')

    def on_title_changed(self,text):
        self.settings['Plot title'] = text
        if len(text) <= 0:
            self.settings['Use default plot title'] == True
        else:
            self.settings['Use default plot title'] == False
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
        self.reader = self.notebook.mainTab.reader
        if self.reader is not None:
            self.cell = self.reader.unit_cells[-1]
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.vinc_sb.setValue(self.settings['Frequency increment'])
        self.title_le.setText(self.settings['Plot title'])
        self.h_sb.setValue(self.settings['Unique direction - h'])
        self.k_sb.setValue(self.settings['Unique direction - k'])
        self.l_sb.setValue(self.settings['Unique direction - l'])
        self.azimuthal_angle_sb.setValue(self.settings['Azimuthal angle'])
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        self.superstrate_dielectric_sb.setValue(self.settings['Superstrate dielectric'])
        self.substrate_dielectric_sb.setValue(self.settings['Substrate dielectric'])
        self.film_thickness_sb.setValue(self.settings['Film thickness'])
        index = self.mode_cb.findText(self.settings['Mode'], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        index = self.funits_cb.findText(self.settings['Frequency units'], Qt.MatchFixedString)
        self.funits_cb.setCurrentIndex(index)
        # Refresh the widgets that depend on the reader
        self.reader = self.notebook.reader
        self.calculate_euler_angles()
        # Flag a recalculation will be required
        self.notebook.singleCrystalCalculationRequired = True
        # Reset the progress bar
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        self.plotLastButtonPressed()
        self.dirty = False
        self.notebook.singleCrystalCalculationRequired = False
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        return

    def on_mode_cb_activated(self, index):
        if index == 0:
            self.settings['Mode'] = 'Thick slab'
        elif index == 1:
            self.settings['Mode'] = 'Coherent thin film'
        else:
            self.settings['Mode'] = 'Incoherent thin film'
        if self.settings['Use default plot title'] == True:
            self.settings['Plot title'] = 'Single crystal '+self.settings['Mode']
        self.notebook.singleCrystalCalculationRequired = True
        self.dirty = True
        self.refresh()
        debugger.print('Mode changed to ', self.settings['Mode'])

    def on_funits_cb_activated(self, index):
        if index == 0:
            self.settings['Frequency units'] = 'wavenumber'
        else:
            self.settings['Frequency units'] = 'THz'
        self.replot()
        debugger.print('Frequency units changed to ', self.settings['Frequency units'])

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
        if self.reader is None:
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
        crystalDepth     = self.settings['Film thickness']
        crystal          = GTM.Layer(thickness=crystalDepth*1e-6,    epsilon=crystalPermittivityFunction)
        # Creat the system with the layers 
        if self.settings['Mode'] == 'Thick slab':
            system = GTM.System(substrate=crystal, superstrate=superstrate, layers=[])
        elif self.settings['Mode'] == 'Coherent thin film':
            system = GTM.System(substrate=substrate, superstrate=superstrate, layers=[crystal])
        else:
            system = GTM.System(substrate=substrate, superstrate=superstrate, layers=[crystal])
        # Determine the euler angles
        theta,phi,psi = self.calculate_euler_angles()
        # Set the angle of incidence in radians
        angle = self.settings['Angle of incidence']
        angleOfIncidence      = np.pi / 180.0 * angle
        # Rotate the dielectric constants to the laboratory frame
        system.substrate.set_euler(theta, phi, psi)
        system.superstrate.set_euler(theta, phi, psi)
        for layer in system.layers:
            layer.set_euler(theta, phi, psi)
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
        results = []
        # Assemble all the parameters we need for parallel execution
        # About to call
        self.progressbar.setMaximum(len(vs))
        progress = 0
        for v in vs:
            #jk system = copy.deepcopy(system)
            call_parameters.append( (v, angleOfIncidence, system) )
            # results.append( Calculator.solve_single_crystal_equations( (v, angleOfIncidence, system) ) )
        for result in pool.map(Calculator.solve_single_crystal_equations, call_parameters, chunksize=40):
            results.append(result)
            progress += 1
            self.progressbar.setValue(progress)
            if self.notebook.progressbar is not None:
                self.notebook.progressbar.setValue(progress)
        QCoreApplication.processEvents()
        # Initialise plotting variables
        self.frequencies_cm1= []
        self.p_reflectivity = []
        self.s_reflectivity = []
        self.p_transmission = []
        self.s_transmission = []
        self.epsilon = []
        for v,r,R,t,T,epsilon in results:
            self.frequencies_cm1.append(v)
            self.p_reflectivity.append(R[0]+R[2])
            self.s_reflectivity.append(R[1]+R[3])
            self.p_transmission.append(T[0])
            self.s_transmission.append(T[1])
            self.epsilon.append(epsilon)
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
        sp.writeNextRow([ 'Single crystal mode',          self.settings['Mode'] ],col=1)
        sp.writeNextRow([ 'Minimum frequency',            self.settings['Minimum frequency'] ],col=1)
        sp.writeNextRow([ 'Maximum frequency',            self.settings['Maximum frequency'] ],col=1)
        sp.writeNextRow([ 'Frequency increment',          self.settings['Frequency increment'] ],col=1)
        sp.writeNextRow([ 'Plot title',                   self.settings['Plot title'] ],col=1)
        sp.writeNextRow([ 'Surface definition (h)',       self.settings['Unique direction - h'] ],col=1)
        sp.writeNextRow([ 'Surface definition (k)',       self.settings['Unique direction - k'] ],col=1)
        sp.writeNextRow([ 'Surface definition (l)',       self.settings['Unique direction - l'] ],col=1)
        sp.writeNextRow([ 'Azimuthal angle',              self.settings['Azimuthal angle'] ],col=1)
        sp.writeNextRow([ 'Angle of incidence',           self.settings['Angle of incidence'] ],col=1)
        sp.writeNextRow([ 'Superstrate dielectric',       self.settings['Superstrate dielectric'] ],col=1)
        sp.writeNextRow([ 'Substrate dielectric',         self.settings['Substrate dielectric'] ],col=1)
        sp.writeNextRow([ 'Film thickness(mu)',           self.settings['Film thickness'] ],col=1)
        headings = ['R_p', 'R_s', 'T_p', 'T_s']
        self.write_results(sp, 'Crystal R&T',     self.frequencies_cm1, [self.p_reflectivity, self.s_reflectivity, self.p_transmission, self.s_transmission], headings)

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

    def plot(self,x,ys,ylabel,legends, subtitle):
        # import matplotlib.pyplot as pl
        # mp.use('Qt5Agg')
        self.subplot = None
        self.remember_x  = x
        self.remember_ys = ys
        self.remember_ylabel = ylabel
        self.remember_legends = legends
        self.remember_subtitle = subtitle
        self.figure.clf()
        if self.settings['Frequency units'] == 'wavenumber':
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
        self.subplot.set_title(self.settings['Plot title']+' ('+subtitle+')')
        self.canvas.draw_idle()

    def replot(self):
        if self.subplot is not None:
            self.plot(self.remember_x, self.remember_ys, self.remember_ylabel, self.remember_legends, self.remember_subtitle)

    def calculate_euler_angles(self):
        '''Calculate the Euler angles for the crystal to lab transformation'''
        # Get plane specification
        hkl = [ self.settings['Unique direction - h'] , self.settings['Unique direction - k'], self.settings['Unique direction - l'] ]
        sum2 = hkl[0]*hkl[0] + hkl[1]*hkl[1] + hkl[2]*hkl[2]
        if sum2 < 1:
            return 0,0,0
        x = 0
        y = 1
        z = 2
        # convert normal to plane to a direction in xyz coordinates
        planez = self.cell.convert_hkl_to_xyz(hkl)
        planez /=  np.linalg.norm(planez)
        plane = np.zeros( (3,3) )
        lab   = np.identity(3)
        plane[z] = planez
        if plane[z][2] < 0.99999999 and plane[z][2] > -0.99999999:
            plane[x] = np.cross(plane[z], lab[z])
            plane[y] = np.cross(plane[z], plane[x])
            plane[x] /= np.linalg.norm(plane[y])
            plane[y] /= np.linalg.norm(plane[y])
        else:
            plane[x] = lab[x]
            plane[y] = lab[y]
            plane[z] = lab[z]
        # Calculate the rotation matrix which transforms us to a unit matrix
        rotation = np.linalg.pinv(plane)
        # Because the formula used for the Euler transform is based on an active transform
        # We calculate the angles using the transpose of the rotation matrix
        rotation = rotation.T
        # Negative sign for angles because of the passive / active problem
        alpha = -np.arctan2(rotation[z][0],-rotation[z][1])
        beta  = -np.arccos(rotation[z][2])
        gamma = -np.arctan2(rotation[x][2],rotation[y][2])
        gamma = -self.settings['Azimuthal angle'] * np.pi / 180.0
        #
        # Some confusion here as to the role of the euler angles
        #
        psi             = alpha
        theta           = beta
        phi             = gamma
        normal_to_plane_lab = Calculator.euler_rotation(plane[z], theta, phi, psi)
        if normal_to_plane_lab[2] < 0.9999 and normal_to_plane_lab[2] > -0.9999:
            print('Error in Euler rotations - surface normal is not along Z-axis', normal_to_plane_lab)
            exit()
        a = Calculator.euler_rotation(self.cell.lattice[0], theta, phi, psi)
        b = Calculator.euler_rotation(self.cell.lattice[1], theta, phi, psi)
        c = Calculator.euler_rotation(self.cell.lattice[2], theta, phi, psi)
        self.labframe_w.clear()
        self.labframe_w.addItem('a-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(a[0],a[1],a[2]) )
        self.labframe_w.addItem('b-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(b[0],b[1],b[2]) )
        self.labframe_w.addItem('c-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(c[0],c[1],c[2]) )
        a = a / np.linalg.norm(a)
        b = b / np.linalg.norm(b)
        c = c / np.linalg.norm(c)
        # print('Projection of a,b,c onto the lab Y-axis (s-pol)', a[1],b[1],c[1])
        # print('Projection of a,b,c onto the lab X-axis (p-pol)', a[0],b[0],c[0])
        return (theta, phi, psi)

