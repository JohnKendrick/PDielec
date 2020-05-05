import sys
import os.path
import os
import numpy as np
osname = os.name
from multiprocessing import cpu_count, Array
import PDielec.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QProgressBar, QApplication
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox,QDoubleSpinBox
from PyQt5.QtWidgets  import  QSizePolicy
from PyQt5.QtCore     import  QCoreApplication, Qt
from PDielec.Constants import  wavenumber, amu, PI, avogadro_si, angstrom
from PDielec.Constants import  average_masses, isotope_masses
import ctypes
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PDielec.Utilities import Debug
import time

def set_affinity_on_worker():
    '''When a new worker process is created, the affinity is set to all CPUs'''
    #JK print('I'm the process %d, setting affinity to all CPUs.' % os.getpid())
    #JK Commented out for the time being
    #JK os.system('taskset -p 0xff %d > /dev/null' % os.getpid())


class PlottingTab(QWidget):
    def __init__(self, parent, debug=False ):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'PlottingTab')
        self.dirty = True
        self.settings = {}
        self.subplot = None
        self.setWindowTitle('Plotting')
        self.settings['Minimum frequency'] = 0
        self.settings['Maximum frequency'] = 400
        self.settings['Frequency increment'] = 0.2
        self.molar_definitions = ['Unit cells','Atoms','Molecules']
        self.settings['Molar definition'] = 'Unit cells'
        self.settings['Number of atoms'] = 1
        self.settings['Plot title'] = 'Plot Title'
        self.xaxes = []
        self.directions = []
        self.depolarisations = []
        self.frequency_units = None
        self.molar_cb_current_index = 0
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Create the parallel environment
        if self.notebook.threading:
            from multiprocessing.dummy import Pool
        else:
            from multiprocessing import Pool
        # Create last tab - PlottingTab
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
        # Final button
        #
        hbox = QHBoxLayout()
        self.molarAbsorptionButton = QPushButton('Plot molar absorption')
        self.molarAbsorptionButton.setToolTip('Plot the molar absorption using the scenarios')
        self.molarAbsorptionButton.clicked.connect(self.molarAbsorptionButtonClicked)
        hbox.addWidget(self.molarAbsorptionButton)
        self.absorptionButton = QPushButton('Plot absorption')
        self.absorptionButton.setToolTip('Plot the absorption using the scenarios')
        self.absorptionButton.clicked.connect(self.absorptionButtonClicked)
        hbox.addWidget(self.absorptionButton)
        self.realPermButton = QPushButton('Plot real permittivity')
        self.realPermButton.setToolTip('Plot the real permittivity using the scenarios')
        self.realPermButton.clicked.connect(self.realPermButtonClicked)
        hbox.addWidget(self.realPermButton)
        self.imagPermButton = QPushButton('Plot imaginary permittivity')
        self.imagPermButton.setToolTip('Plot the imaginary permittivity using the scenarios')
        self.imagPermButton.clicked.connect(self.imagPermButtonClicked)
        hbox.addWidget(self.imagPermButton)
        self.atrButton = QPushButton('Plot ATR')
        self.atrButton.setToolTip('Plot Reflectance ATR')
        self.atrButton.clicked.connect(self.atrButtonClicked)
        hbox.addWidget(self.atrButton)
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

    def molarAbsorptionButtonClicked(self):
        debugger.print('molarAbsorptionButtonClicked pressed')
        if self.notebook.plottingCalculationRequired:
            self.calculate()
        if not self.notebook.plottingCalculationRequired:
            self.plot(self.xaxes, self.molarAbsorptionCoefficients, r'Molar Absorption Coefficient $\mathdefault{(L mole^{-1} cm^{-1})}$')

    def absorptionButtonClicked(self):
        debugger.print('absorptionButtonClicked pressed')
        if self.notebook.plottingCalculationRequired:
            self.calculate()
        if not self.notebook.plottingCalculationRequired:
            self.plot(self.xaxes, self.absorptionCoefficients, r'Absorption Coefficient $\mathdefault{(cm^{-1})}$')

    def realPermButtonClicked(self):
        debugger.print('realPermButtonClicked pressed')
        if self.notebook.plottingCalculationRequired:
            self.calculate()
        if not self.notebook.plottingCalculationRequired:
            self.plot(self.xaxes, self.realPermittivities, 'Real Component of Permittivity')

    def imagPermButtonClicked(self):
        debugger.print('imagPermButtonClicked pressed')
        if self.notebook.plottingCalculationRequired:
            self.calculate()
        if not self.notebook.plottingCalculationRequired:
            self.plot(self.xaxes, self.imagPermittivities, 'Imaginary Component of Permittivity')

    def atrButtonClicked(self):
        debugger.print('atrButtonClicked pressed')
        if self.notebook.plottingCalculationRequired:
            self.calculate()
        if not self.notebook.plottingCalculationRequired:
            self.plot(self.xaxes, self.sp_atrs, 'ATR absorption')

    def on_title_changed(self,text):
        self.settings['Plot title'] = text
        if self.subplot is not None:
            self.subplot.set_title(self.settings['Plot title'])
            self.canvas.draw_idle()
        debugger.print('on title change ', self.settings['Plot title'])

    def on_vinc_changed(self,value):
        self.settings['Frequency increment'] = value
        self.notebook.plottingCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        debugger.print('on vinc change ', self.settings['Frequency increment'])

    def on_vmin_changed(self):
        self.settings['Minimum frequency'] = self.vmin_sb.value()
        self.notebook.plottingCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        debugger.print('on vmin change ', self.settings['Minimum frequency'])

    def on_vmax_changed(self):
        self.settings['Maximum frequency'] = self.vmax_sb.value()
        self.notebook.plottingCalculationRequired = True
        self.dirty = True
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        QCoreApplication.processEvents()
        debugger.print('on vmax change ', self.settings['Maximum frequency'])

    def refresh(self,force=False):
        if not self.dirty and not force and not self.notebook.plottingCalculationRequired:
            debugger.print('refreshing widget aborted', self.dirty,force,self.notebook.plottingCalculationRequired)
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
        self.notebook.plottingCalculationRequired = True
        # Reset the progress bar
        self.progressbar.setValue(0)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setValue(0)
        self.molarAbsorptionButtonClicked()
        self.dirty = False
        self.notebook.plottingCalculationRequired = False
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        return

    def on_natoms_changed(self, value):
        self.notebook.plottingCalculationRequired = True
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
        self.notebook.plottingCalculationRequired = True
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
        debugger.print('calculate')
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
        self.notebook.plottingCalculationRequired = False
        settings = self.notebook.settingsTab.settings
        eckart = settings['Eckart flag']
        neutral = settings['Neutral Born charges']
        epsilon_inf = np.array(settings['Optical permittivity'])
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        frequencies = np.array(frequencies_cm1) * wavenumber
        intensities = self.notebook.settingsTab.intensities
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = reader.volume*angstrom*angstrom*angstrom
        vmin = self.settings['Minimum frequency']
        vmax = self.settings['Maximum frequency']
        vinc = self.settings['Frequency increment']
        calling_parameters = []
        self.scenarios = self.notebook.scenarios
        mode_list = []
        drude_sigma = 0
        drude_plasma = 0
        drude = False
        for mode_index,selected in enumerate(modes_selected):
            if selected:
                mode_list.append(mode_index)
        #debugger.print('mode_list', mode_list)
        # Calculate the ionic permittivity at zero frequency
        epsilon_ionic = self.notebook.settingsTab.epsilon_ionic
        epsilon_total = epsilon_inf + epsilon_ionic
        cell = reader.unit_cells[-1]
        self.directions = []
        self.depolarisations = []
        # Process the shapes and the unique directions.
        # and calculate the depolarisation matrices
        self.legends = []
        for scenario in self.scenarios:
            self.legends.append(scenario.settings['Legend'])
            #debugger.print('Scenario ',scenario.scenarioIndex)
            shape = scenario.settings['Particle shape']
            #debugger.print('Particle shape ',shape)
            hkl = [scenario.settings['Unique direction - h'], scenario.settings['Unique direction - k'], scenario.settings['Unique direction - l']]
            aoverb = scenario.settings['Ellipsoid a/b']
            if shape == 'Ellipsoid':
                direction = cell.convert_abc_to_xyz(hkl)
                depolarisation = Calculator.initialise_ellipsoid_depolarisation_matrix(direction,aoverb)
            elif shape == 'Plate':
                direction = cell.convert_hkl_to_xyz(hkl)
                depolarisation = Calculator.initialise_plate_depolarisation_matrix(direction)
            elif shape == 'Needle':
                direction = cell.convert_abc_to_xyz(hkl)
                depolarisation = Calculator.initialise_needle_depolarisation_matrix(direction)
            else:
                depolarisation = Calculator.initialise_sphere_depolarisation_matrix()
                direction = np.array( [] )
            direction = direction / np.linalg.norm(direction)
            self.directions.append(direction)
            #debugger.print('direction',direction)
            self.depolarisations.append(depolarisation)
        # Set up the parallel processing requirements before looping over the frequencies
        call_parameters = []
        for v in np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc)):
            vau = v * wavenumber
            call_parameters.append( (v, vau, mode_list, frequencies, sigmas, oscillator_strengths,
                                     volume, epsilon_inf, drude, drude_plasma, drude_sigma) )

        maximum_progress = len(call_parameters) * (1 + len(self.scenarios))
        progress = 0
        self.progressbar.setMaximum(maximum_progress)
        if self.notebook.progressbar is not None:
            self.notebook.progressbar.setMaximum(maximum_progress)
        QCoreApplication.processEvents()
        # Lets find ouit how many processors to use
        if self.notebook.ncpus == 0:
            number_of_processors = cpu_count()
        else:
            number_of_processors = self.notebook.ncpus
        #jk start = time.time()
        if self.notebook.threading:
            from multiprocessing.dummy import Pool
            pool = Pool(number_of_processors)
        else:
            from multiprocessing import Pool
            pool = Pool(number_of_processors, initializer=set_affinity_on_worker, maxtasksperchild=10)
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
        # Allocate space for the shared memory, we need twice as much as we have a complex data type
        shared_array_base = Array(ctypes.c_double, 2*nplots*9)
        previous_solution_shared = np.ctypeslib.as_array(shared_array_base.get_obj())
        # Convert the space allocated to complex
        previous_solution_shared.dtype = np.complex128
        # Reshape the array and fill everything with zero's
        previous_solution_shared = previous_solution_shared.reshape(nplots, 3,3)
        previous_solution_shared.fill(0.0+0.0j)
        # Prepare parallel call parameters for the loop over frequencies, methods, volume fractions
        concentration = self.settings['concentration']
        # Assemble each scenario settings
        self.xaxes = []
        self.realPermittivities = []
        self.imagPermittivities = []
        self.absorptionCoefficients = []
        self.molarAbsorptionCoefficients = []
        self.sp_atrs = []
        #jk start = time.time()
        if self.notebook.threading:
            from multiprocessing.dummy import Pool
            pool = Pool(number_of_processors)
        else:
            from multiprocessing import Pool
            pool = Pool(number_of_processors, initializer=set_affinity_on_worker, maxtasksperchild=10)
        for i,(scenario,L) in enumerate(zip(self.scenarios,self.depolarisations)):
            #debugger.print('Scenario ',i,L)
            matrix = scenario.settings['Matrix']
            method = scenario.settings['Effective medium method'].lower()
            matrix_density = scenario.settings['Matrix density']
            matrix_permittivity = np.identity(3) * scenario.settings['Matrix permittivity']
            mass_fraction = scenario.settings['Mass fraction']
            volume_fraction = scenario.settings['Volume fraction']
            particle_size_mu = scenario.settings['Particle size(mu)']
            particle_sigma_mu = scenario.settings['Particle size distribution sigma(mu)']
            shape = scenario.settings['Particle shape'].lower()
            hkl = [scenario.settings['Unique direction - h'], scenario.settings['Unique direction - k'], scenario.settings['Unique direction - l']]
            atr_refractive_index = scenario.settings['ATR material refractive index']
            atr_theta = scenario.settings['ATR theta']
            atr_spolfraction = scenario.settings['ATR S polarisation fraction']
            aoverb = scenario.settings['Ellipsoid a/b']
            bubble_vf = scenario.settings['Bubble volume fraction']
            bubble_radius = scenario.settings['Bubble radius']
            vf_type = ''
            call_parameters = []
            nplot = 0
            for v,vau,dielecv in dielecv_results:
                # convert the size to a dimensionless number which is 2*pi*size/lambda
                lambda_mu = 1.0E4 / (v + 1.0e-12)
                if particle_size_mu < 1.0e-12:
                    particle_size_mu = 1.0e-12
                size = 2.0*PI*particle_size_mu / lambda_mu
                data = ''
                call_parameters.append( (v,vau,dielecv,method,volume_fraction,vf_type,particle_size_mu,particle_sigma_mu,size,nplot,
                                         matrix_permittivity,shape,data,L,concentration,atr_refractive_index,atr_theta,atr_spolfraction,bubble_vf,bubble_radius,previous_solution_shared) )
                nplot += 1
            results = []
            for result in pool.imap(Calculator.solve_effective_medium_equations, call_parameters, chunksize=40):
                results.append(result)
                progress += 1
                self.progressbar.setValue(progress)
                if self.notebook.progressbar is not None:
                    self.notebook.progressbar.setValue(progress)
            QCoreApplication.processEvents()
            xaxis = []
            realPermittivity = []
            imagPermittivity = []
            absorptionCoefficient = []
            molarAbsorptionCoefficient = []
            sp_atr = []
            for v,nplot,method,vf_type,size_mu,size_sigma,shape,data,trace, absorption_coefficient,molar_absorption_coefficient,spatr in results:
                 xaxis.append(v)
                 realPermittivity.append(np.real(trace))
                 imagPermittivity.append(np.imag(trace))
                 absorptionCoefficient.append(absorption_coefficient)
                 molarAbsorptionCoefficient.append(molar_absorption_coefficient)
                 sp_atr.append(spatr)
            self.xaxes.append(xaxis)
            self.realPermittivities.append(realPermittivity)
            self.imagPermittivities.append(imagPermittivity)
            self.absorptionCoefficients.append(absorptionCoefficient)
            self.molarAbsorptionCoefficients.append(molarAbsorptionCoefficient)
            self.sp_atrs.append(sp_atr)
        pool.close()
        pool.join()
        #jk print('Dielec calculation duration ', time.time()-start)
        #if self.notebook.spreadsheet is not None:
        #    self.write_spreadsheet()
        self.dirty = False
        QApplication.restoreOverrideCursor()
        QCoreApplication.processEvents()

    def write_spreadsheet(self):
        if self.notebook.spreadsheet is None:
            return
        # Deal with Scenarios first
        sp = self.notebook.spreadsheet
        sp.selectWorkSheet('Scenarios')
        sp.delete()
        sp.writeNextRow(['A list of the scenarios used the calculation of the effective medium'],col=1)
        for index,(scenario,direction,depolarisation) in enumerate(zip(self.scenarios,self.directions,self.depolarisations)):
            sp.writeNextRow([''],col=1)
            sp.writeNextRow(['Scenario '+str(index)],col=1,check=1)
            settings = scenario.settings
            for key in sorted(settings,key=str.lower):
                sp.writeNextRow([key, settings[key]],col=1,check=1)
            sp.writeNextRow(['Normalised unique direction']+direction.tolist(), col=1,check=1)
            sp.writeNextRow(['Depolarisation matrix'], col=1,check=1)
            sp.writeNextRow(depolarisation[0].tolist(), col=2, check=1)
            sp.writeNextRow(depolarisation[1].tolist(), col=2, check=1)
            sp.writeNextRow(depolarisation[2].tolist(), col=2, check=1)
        # Now deal with Molar absorption, absorption, real and imaginary permittivity
        self.write_results(sp, 'Molar Absorption',self.xaxes, self.molarAbsorptionCoefficients)
        self.write_results(sp, 'Absorption',self.xaxes, self.absorptionCoefficients)
        self.write_results(sp, 'Real Permittivity',self.xaxes, self.realPermittivities)
        self.write_results(sp, 'Imaginary Permittivity',self.xaxes, self.imagPermittivities)
        self.write_results(sp, 'ATR Reflectance',self.xaxes, self.sp_atrs)

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


