import sys
import os.path
import numpy as np
from multiprocessing import Pool, cpu_count, Array
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QFileDialog, QFrame
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox
from PyQt5.QtCore     import  Qt
from Python.Constants import  wavenumber, amu, PI, avogadro_si
from Python.Constants import  average_masses, isotope_masses
import ctypes
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

def set_affinity_on_worker():
    """When a new worker process is created, the affinity is set to all CPUs"""
    #JK print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    #JK Commented out for the time being
    #JK os.system("taskset -p 0xff %d > /dev/null" % os.getpid())


class PlottingTab(QWidget):
    def __init__(self, parent ):   
        super(QWidget, self).__init__(parent)
        self.settings = {}
        self.subplot = None
        self.setWindowTitle("Plotting")
        self.settings["vmin"] = 0
        self.settings["vmax"] = 400
        self.settings["vinc"] = 0.2
        self.settings["spreadsheet"] = ""
        self.molar_definitions = ["Unit cells","Atoms","Molecules"]
        self.settings["molar_definition"] = "Unit cells"
        self.settings["natoms"] = 1
        self.settings["title"] = "Plot Title"
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Create last tab - PlottingTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The minimum frequency
        #
        self.vmin_sb = QSpinBox(self)
        self.vmin_sb.setRange(0,9000)
        self.vmin_sb.setValue(self.settings["vmin"])
        self.vmin_sb.setToolTip("Set the minimum frequency to be considered)")
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        label = QLabel("Minimum frequency:", self)
        label.setToolTip("Set the minimum frequency to be considered)")
        form.addRow(label, self.vmin_sb)
        #
        # The maximum frequency
        #
        self.vmax_sb = QSpinBox(self)
        self.vmax_sb.setRange(0,9000)
        self.vmax_sb.setValue(self.settings["vmax"])
        self.vmax_sb.setToolTip("Set the maximum frequency to be considered)")
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        label = QLabel("Maximum frequency:", self)
        label.setToolTip("Set the maximum frequency to be considered)")
        form.addRow(label, self.vmax_sb)
        # 
        # Choose a suitable increment
        #
        self.vinc_le = QLineEdit(self) 
        self.vinc_le.setToolTip("Choose an increment for the frequency when plotting")
        self.vinc_le.setText("{0:.2f}".format(self.settings["vinc"]))
        self.vinc_le.textChanged.connect(self.on_vinc_changed)
        label = QLabel("Frequency increment", self)
        label.setToolTip("Choose an increment for the frequency when plotting")
        form.addRow(label, self.vinc_le)
        # 
        # Define molar quantity
        #
        self.molar_cb = QComboBox(self) 
        self.molar_cb.blockSignals(True)
        self.molar_cb.setToolTip("Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given")
        self.molar_cb.currentIndexChanged.connect(self.on_molar_cb_changed)
        for definition in self.molar_definitions:
            self.molar_cb.addItem(definition)
        self.molar_cb.setCurrentIndex(0)
        label = QLabel("Molar definition", self)
        label.setToolTip("Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given")
        form.addRow(label, self.molar_cb)
        self.molar_cb.blockSignals(False)
        #
        # Number of atoms in a molecule
        # 
        self.natoms_sb = QSpinBox(self)
        self.natoms_sb.setToolTip("Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed")
        self.natoms_sb.setRange(1,500)
        self.natoms_sb.setValue(self.settings["natoms"])
        self.natoms_sb.valueChanged.connect(self.on_natoms_changed)
        self.natoms_sb.setEnabled(False)
        label = QLabel("Number of atoms per molecule", self)
        label.setToolTip("Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed")
        form.addRow(label, self.natoms_sb)
        # 
        # Store results in a file?
        #
        self.file_store_le = QLineEdit(self) 
        self.file_store_le.setToolTip("Store the results in a .csv or .xlsx file")
        self.file_store_le.setText(self.settings["spreadsheet"])
        self.file_store_le.textChanged.connect(self.on_file_store_le_changed)
        form.addRow(QLabel("Output spreadsheet", self), self.file_store_le)
        # 
        # Set the plot title         
        #
        self.title_le = QLineEdit(self) 
        self.title_le.setToolTip("Set the plot title")
        self.title_le.setText(self.settings["title"])
        self.title_le.textChanged.connect(self.on_title_changed)
        label = QLabel("Plot title", self)
        label.setToolTip("Set the plot title")
        form.addRow(label, self.title_le)
        #
        # Final button
        #
        hbox = QHBoxLayout()
        self.molarAbsorptionButton = QPushButton("Plot the molar absorption")
        self.molarAbsorptionButton.setToolTip("Plot the molar absorption using the scenarios")
        self.molarAbsorptionButton.clicked.connect(self.molarAbsorptionButtonClicked)
        hbox.addWidget(self.molarAbsorptionButton)
        self.absorptionButton = QPushButton("Plot the absorption")
        self.absorptionButton.setToolTip("Plot the absorption using the scenarios")
        self.absorptionButton.clicked.connect(self.absorptionButtonClicked)
        hbox.addWidget(self.absorptionButton)
        self.realPermButton = QPushButton("Plot the real permittivity")
        self.realPermButton.setToolTip("Plot the real permittivity using the scenarios")
        self.realPermButton.clicked.connect(self.realPermButtonClicked)
        hbox.addWidget(self.realPermButton)
        self.imagPermButton = QPushButton("Plot the imaginary permittivity")
        self.imagPermButton.setToolTip("Plot the imaginary permittivity using the scenarios")
        self.imagPermButton.clicked.connect(self.imagPermButtonClicked)
        hbox.addWidget(self.imagPermButton)
        form.addRow(hbox)
        self.info_la = QLabel("Calculation status", self)
        self.info_la.setFrameShape(QFrame.Panel)
        self.info_la.setFrameShadow(QFrame.Sunken)
        self.info_la.setLineWidth(3)
        label = QLabel("Calculation status", self)
        form.addRow(label,self.info_la)
        # Add the matplotlib figure to the bottom 
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        form.addRow(self.canvas)
        form.addRow(self.toolbar)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)

    def molarAbsorptionButtonClicked(self):
        print("molarAbsorptionButtonClicked pressed")
        if self.notebook.newCalculationRequired:
            self.calculate()
        self.plot(self.xaxes, self.molarAbsorptionCoefficients, "molar absorption coefficient")

    def absorptionButtonClicked(self):
        print("absorptionButtonClicked pressed")
        if self.notebook.newCalculationRequired:
            self.calculate()
        self.plot(self.xaxes, self.absorptionCoefficients, "absorption coefficient")

    def realPermButtonClicked(self):
        print("realPermButtonClicked pressed")
        if self.notebook.newCalculationRequired:
            self.calculate()
        self.plot(self.xaxes, self.realPermittivities, "real permittivity")

    def imagPermButtonClicked(self):
        print("imagPermButtonClicked pressed")
        if self.notebook.newCalculationRequired:
            self.calculate()
        self.plot(self.xaxes, self.imagPermittivities, "imaginary permittivity")

    def on_file_store_le_changed(self,text):
        self.settings["spreadsheet"] = text
        print("on file_store_le change ", self.settings["spreadsheet"])

    def on_title_changed(self,text):
        self.settings["title"] = text
        if self.subplot is not None:
            self.subplot.set_title(self.settings["title"])
            self.canvas.draw_idle()
        print("on title change ", self.settings["title"])

    def on_vinc_changed(self,text):
        self.settings["vinc"] = float(text)
        self.notebook.newCalculationRequired = True
        print("on vinc change ", self.settings["vinc"])

    def on_vmin_changed(self):
        self.settings["vmin"] = self.vmin_sb.value()
        self.notebook.newCalculationRequired = True
        print("on vmin change ", self.settings["vmin"])

    def on_vmax_changed(self):
        self.settings["vmax"] = self.vmax_sb.value()
        self.notebook.newCalculationRequired = True
        print("on vmax change ", self.settings["vmax"])

    def refresh(self):
        print("Plotting tab - refreshing widget")
        self.info_la.setText("Refreshing screen")
        # Refresh the widgets that depend on the reader
        self.reader = self.notebook.reader
        if self.reader is not None:
            self.settings["concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
        # Flag a recalculation will be required
        self.notebook.newCalculationRequired = True
        self.info_la.setText("Finished refreshing screen")
        return

    def on_natoms_changed(self, value):
        self.notebook.newCalculationRequired = True
        self.settings["natoms"] = value
        print("on natoms changed ", self.settings["natoms"])

    def on_molar_cb_changed(self, index):
        self.settings["molar_definition"] = self.molar_definitions[index]
        self.notebook.newCalculationRequired = True
        if self.settings["molar_definition"] == "Molecules":
            self.settings["concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 * self.settings["natoms"] / self.reader.nions)
            self.natoms_sb.setEnabled(True)
        elif self.settings["molar_definition"] == "Unit cells":
            self.settings["concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
            self.natoms_sb.setEnabled(False)
        elif self.settings["molar_definition"] == "Atoms":
            self.settings["concentration"] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 / self.reader.nions)
            self.natoms_sb.setEnabled(False)
        print("The concentration has been set", self.settings["molar_definition"], self.settings["concentration"])

    def calculate(self):
        print("calculate")
        self.info_la.setText("Starting calculation")
        self.notebook.newCalculationRequired = False
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings["program"]
        filename = settings["filename"]
        reader = self.notebook.mainTab.reader
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        eckart = settings["eckart"]
        neutral = settings["neutral"]
        hessian_symm = settings["hessian_symmetrisation"]
        epsilon_inf = np.array(settings["optical_permittivity"])
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        frequencies = np.array(frequencies_cm1) * wavenumber
        intensities = self.notebook.settingsTab.intensities
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = reader.volume
        vmin = self.settings["vmin"]
        vmax = self.settings["vmax"]
        vinc = self.settings["vinc"]
        calling_parameters = []
        scenarios = self.notebook.scenarios
        mode_list = []
        drude_sigma = 0
        drude_plasma = 0
        drude = False
        for mode_index,selected in enumerate(modes_selected):
            if selected:
                mode_list.append(mode_index)
        print("mode_list", mode_list)
        print("length mode_list",len(mode_list))
        print("length frequencies",len(frequencies_cm1))
        # Calculate the ionic permittivity at zero frequency
        epsilon_ionic = Calculator.ionic_permittivity(mode_list, oscillator_strengths, frequencies, volume )
        print('epsilon_ionic',epsilon_ionic)
        epsilon_total = epsilon_inf + epsilon_ionic
        print('epsilon_total',epsilon_ionic)
        cell = reader.unit_cells[-1]
        directions = []
        depolarisations = []
        # Process the shapes and the unique directions.
        # and calculate the depolarisation matrices
        for scenario in scenarios:
            print("Scenario ",scenario.scenarioIndex)
            shape = scenario.settings["shape"]
            print("shape ",shape)
            hkl = [scenario.settings["h"], scenario.settings["k"], scenario.settings["l"]] 
            print("shape ",shape)
            aoverb = scenario.settings["aoverb"]
            if shape == "Ellipsoid":
                direction = cell.convert_abc_to_xyz(hkl)
                depolarisation = Calculator.initialise_ellipsoid_depolarisation_matrix(direction,aoverb)
            elif shape == "Plate":
                direction = cell.convert_hkl_to_xyz(hkl)
                depolarisation = Calculator.initialise_plate_depolarisation_matrix(direction)
            elif shape == "Needle":
                direction = cell.convert_abc_to_xyz(hkl)
                depolarisation = Calculator.initialise_needle_depolarisation_matrix(direction)
            else:
                depolarisation = Calculator.initialise_sphere_depolarisation_matrix()
                direction = np.array( [] )
            direction = direction / np.linalg.norm(direction)
            directions.append(direction)
            print('direction',direction)
            depolarisations.append(depolarisation)
        # Set up the parallel processing requirements before looping over the frequencies
        call_parameters = []
        for v in np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc)):
            vau = v * wavenumber
            call_parameters.append( (v, vau, mode_list, frequencies, sigmas, oscillator_strengths,
                                     volume, epsilon_inf, drude, drude_plasma, drude_sigma) )
        number_of_processors = cpu_count()
        self.info_la.setText("Starting calculation of frequency dependent permittivity")
        p = Pool(number_of_processors, initializer=set_affinity_on_worker, maxtasksperchild=10)
        dielecv_results = p.map(Calculator.parallel_dielectric, call_parameters)
        print("Number of frequencies", len(dielecv_results))
        print("Length of call parameters", len(call_parameters))
        #nplots = len(dielecv_results)*len(scenarios)
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
        concentration = self.settings["concentration"]
        # Assemble each scenario settings
        self.xaxes = []
        self.realPermittivities = []
        self.imagPermittivities = []
        self.absorptionCoefficients = []
        self.molarAbsorptionCoefficients = []
        for i,(scenario,L) in enumerate(zip(self.notebook.scenarios,depolarisations)):
            self.info_la.setText("Starting calculation absorption for scenario {}".format(i))
            print('Scenario ',i,L)
            matrix = scenario.settings["matrix"]
            method = scenario.settings["method"].lower()
            matrix_density = scenario.settings["matrix_density"]
            matrix_permittivity = np.identity(3) * scenario.settings["matrix_permittivity"]
            mass_fraction = scenario.settings["mass_fraction"]
            volume_fraction = scenario.settings["volume_fraction"]
            particle_size_mu = scenario.settings["particle_size"]
            particle_sigma = scenario.settings["particle_sigma"]
            shape = scenario.settings["shape"].lower()
            hkl = [scenario.settings["h"], scenario.settings["k"], scenario.settings["l"]] 
            aoverb = scenario.settings["aoverb"]
            # convert the size to a dimensionless number which is 2*pi*size/lambda
            lambda_mu = 1.0E4 / (v + 1.0e-12)
            if particle_size_mu < 1.0e-12:
                particle_size_mu = 1.0e-12
            size = 2.0*PI*particle_size_mu / lambda_mu
            data = ''
            vf_type = ''
            call_parameters = []
            nplot = 0
            for v,vau,dielecv in dielecv_results:
                vau = v * wavenumber
                call_parameters.append( (v,vau,dielecv,method,volume_fraction,vf_type,particle_size_mu,particle_sigma,size,nplot,
                                         matrix_permittivity,dielecv,shape,data,L,concentration,previous_solution_shared) )
                nplot += 1
            print("Length of call parameters", len(call_parameters))
            results = p.map(Calculator.solve_effective_medium_equations, call_parameters)
            print("Method", method, "Shape",shape)
            xaxis = []
            realPermittivity = []
            imagPermittivity = []
            absorptionCoefficient = []
            molarAbsorptionCoefficient = []
            for v,nplot,method,vf_type,size_mu,size_sigma,shape,data,trace, absorption_coefficient,molar_absorption_coefficient in results:
                 xaxis.append(v)
                 realPermittivity.append(np.real(trace))
                 imagPermittivity.append(np.imag(trace))
                 absorptionCoefficient.append(absorption_coefficient)
                 molarAbsorptionCoefficient.append(molar_absorption_coefficient)
            self.xaxes.append(xaxis)
            self.realPermittivities.append(realPermittivity)
            self.imagPermittivities.append(imagPermittivity)
            self.absorptionCoefficients.append(absorptionCoefficient)
            self.molarAbsorptionCoefficients.append(molarAbsorptionCoefficient)
        p.close()
        p.join()
        self.info_la.setText("Finished calculation".format(i))

    def plot(self,xs,ys,ylabel):
        # import matplotlib.pyplot as pl
        # mp.use('Qt5Agg')
        self.figure.clf()
        plot_list = []
        plot_labels = []
        self.subplot = self.figure.add_subplot(111)
        for x,y in zip(xs,ys):
            line, = self.subplot.plot(x,y,lw=2 )
        self.subplot.set_xlabel('Frequency (cm-1)')
        self.subplot.set_ylabel(ylabel)
        self.subplot.set_title(self.settings["title"])
        self.canvas.draw_idle()
 
