import math
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets         import  QPushButton, QWidget
from PyQt5.QtWidgets         import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets         import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets         import  QSpinBox, QTabWidget, QDoubleSpinBox
from PyQt5.QtWidgets         import  QSizePolicy
from PyQt5.QtCore            import  Qt
from Python.Constants        import  wavenumber, amu, PI, avogadro_si, angstrom
from Python.Constants        import  covalent_radii, elemental_colours
# Import plotting requirements
from Python.Utilities        import Debug
from Python.Plotter          import print_strings, print_reals
from Python.GUI.OpenGLWidget import OpenGLWidget

class ViewerTab(QWidget):
    def __init__(self, parent, debug=False ):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'ViewerTab')
        self.debug = debug
        self.dirty = True
        self.subplot = None
        self.setWindowTitle('Viewer')
        self.selected_mode = 0
        self.settings = {}
        self.settings['Atom scaling'] = 0.5
        self.settings['Maximum displacement'] = 1.0
        self.settings['Bond color'] = [ 0.3, 0.3, 0.3 ]
        self.settings['Bond radius'] = 0.1
        self.settings['Cell colour'] = [ 1.0, 0.0, 0.0 ]
        self.settings['Cell radius'] = 0.1
        self.settings['Number of phase steps'] = 41
        self.light_switches = [False]*8
        self.light_switches[0] = True
        self.light_switches[1] = True
        self.plot_types = ['Animation','Arrows']
        self.plot_type_index = 1
        self.number_of_molecules = 0
        self.cell_edges = None
        self.cell_corners = None
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        self.frequencies_cm1  = self.notebook.settingsTab.frequencies_cm1
        # Create last tab - ViewerTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The selected mode
        #
        self.selected_mode_sb = QSpinBox(self)
        if self.frequencies_cm1 is not None:
            self.selected_mode_sb.setRange(0,len(self.frequencies_cm1)-1)
        else:
            self.selected_mode_sb.setRange(0,2000)
        self.selected_mode_sb.setValue(self.selected_mode)
        self.selected_mode_sb.setToolTip('Set the mode to be visualised')
        self.selected_mode_sb.valueChanged.connect(self.on_selected_mode_changed)
        label = QLabel('Select phonon mode', self)
        label.setToolTip('Set the mode to be visualised')
        form.addRow(label, self.selected_mode_sb)
        #
        # The atom scaling
        #
        self.atom_scaling_sb = QDoubleSpinBox(self) 
        self.atom_scaling_sb.setRange(0.01,5.0)
        self.atom_scaling_sb.setSingleStep(0.1)
        self.atom_scaling_sb.setDecimals(2)
        self.atom_scaling_sb.setValue(self.settings['Atom scaling'])
        self.atom_scaling_sb.setToolTip('Scale the covalent radii to determine the size of the atomic spheres')
        self.atom_scaling_sb.valueChanged.connect(self.on_atom_scaling_changed)
        #
        # The bond width
        #
        self.bond_radius_sb = QDoubleSpinBox(self) 
        self.bond_radius_sb.setRange(0.01,1.0)
        self.bond_radius_sb.setSingleStep(0.1)
        self.bond_radius_sb.setDecimals(2)
        self.bond_radius_sb.setValue(self.settings['Bond radius'])
        self.bond_radius_sb.setToolTip('Determines the size of the bonds drawn')
        self.bond_radius_sb.valueChanged.connect(self.on_bond_radius_changed)
        #
        # The cell width
        #
        self.cell_radius_sb = QDoubleSpinBox(self) 
        self.cell_radius_sb.setRange(0.01,1.0)
        self.cell_radius_sb.setSingleStep(0.1)
        self.cell_radius_sb.setDecimals(2)
        self.cell_radius_sb.setValue(self.settings['Cell radius'])
        self.cell_radius_sb.setToolTip('Determines the size of the cell outline')
        self.cell_radius_sb.valueChanged.connect(self.on_cell_radius_changed)
        #
        # The maximum displacement
        #
        self.maximum_displacement_sb = QDoubleSpinBox(self) 
        self.maximum_displacement_sb.setRange(0.1,10.0)
        self.maximum_displacement_sb.setSingleStep(0.1)
        self.maximum_displacement_sb.setDecimals(1)
        self.maximum_displacement_sb.setValue(self.settings['Maximum displacement'])
        self.maximum_displacement_sb.setToolTip('Set the size of the maximum displacement')
        self.maximum_displacement_sb.valueChanged.connect(self.on_maximum_displacement_changed)
        #
        # Add a comb box to select which type of plot
        #
        self.plottype_cb = QComboBox(self) 
        self.plottype_cb.setToolTip('The plot can either be an animation or the modes can be shown by arrow')
        for choice in self.plot_types:
            self.plottype_cb.addItem(choice)
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        self.plottype_cb.currentIndexChanged.connect(self.on_plottype_cb_changed)
        label = QLabel('Choose the plot type', self)
        label.setToolTip('The plot can either be an animation or the modes can be shown by arrow')
        form.addRow(label, self.plottype_cb)
        # 
        # Light switches
        #
        self.light_switches_cb = QComboBox(self)
        self.light_switches_cb.setToolTip('Toogle the light switches on or off')
        for index,light in enumerate(self.light_switches):
            if light:
                string = 'switch light {} off'.format(index)
            else:
                string = 'switch light {} on'.format(index)
            self.light_switches_cb.addItem(string)
        #self.light_switches_cb.currentIndexChanged.connect(self.on_light_switches_cb_changed)
        self.light_switches_cb.activated.connect(self.on_light_switches_cb_activated)
        #
        # Add a tab widget for the settings
        #
        self.settingsTab = QTabWidget(self)
        #self.settingsTab.currentChanged.connect(self.on_settingsTab_currentChanged)
        self.settingsTab.addTab(self.maximum_displacement_sb, 'Maximum Displacement')
        self.settingsTab.addTab(self.cell_radius_sb, 'Cell Width')
        self.settingsTab.addTab(self.bond_radius_sb, 'Bond Width')
        self.settingsTab.addTab(self.atom_scaling_sb, 'Atom Scaling')
        self.settingsTab.addTab(self.light_switches_cb, 'Lighting')
        label = QLabel('Settings', self)
        form.addRow(label,self.settingsTab)
        #
        # Add the opengl figure to the bottom 
        #
        self.opengl_widget = OpenGLWidget(self, debug=debug)
        self.opengl_widget.setFocusPolicy(Qt.StrongFocus)
        self.opengl_widget.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding))
        form.addRow(self.opengl_widget)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)

    def on_light_switches_cb_activated(self, index):
        self.light_switches[index] = not self.light_switches[index]
        if self.light_switches[index]:
            string = 'switch light {} off'.format(index)
        else:
            string = 'switch light {} on'.format(index)
        self.light_switches_cb.setItemText(index,string)
        self.opengl_widget.defineLights()
        self.calculate()
        self.plot()

    def on_maximum_displacement_changed(self,value):
        debugger.print('on maximum_displacement changed ', value)
        self.settings['Maximum displacement'] = value
        self.calculate()
        self.plot()
        
    def on_atom_scaling_changed(self,value):
        debugger.print('on atom_scaling changed ', value)
        self.settings['Atom scaling'] = value
        self.calculate()
        self.plot()
        
    def on_cell_radius_changed(self,value):
        debugger.print('on cell_radius changed ', value)
        self.settings['Cell radius'] = value
        self.calculate()
        self.plot()
        
    def on_bond_radius_changed(self,value):
        self.settings['Bond radius'] = value
        self.calculate()
        self.plot()
        
    def on_selected_mode_changed(self):
        self.selected_mode = self.selected_mode_sb.value()
        debugger.print('on selected_mode change ', self.selected_mode)
        self.opengl_widget.deleteArrows()
        maxR = np.max( np.abs(np.array(self.UVW[self.selected_mode])))
        arrow_scaling = self.settings['Maximum displacement'] / maxR
        for uvw in self.UVW[self.selected_mode]:
            self.opengl_widget.addArrows( uvw, arrow_scaling )
        self.calculate()
        self.plot()

    def on_plottype_cb_changed(self, index):
        self.plot_type_index = index
        debugger.print('Plot type index changed to ', self.plot_type_index)
        self.plot()

    def calculate(self):
        debugger.print('calculate')
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['Program']
        filename = settings['Output file name']
        debugger.print('calculate program file name',program, filename)
        self.reader = self.notebook.reader
        if self.reader is None:
            return
        if program is '':
            return
        if filename is '':
            return
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        debugger.print(' settingsTab settings')
        debugger.print(settings)
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        self.selected_mode_sb.setRange(0,len(self.frequencies_cm1)-1)
        debugger.print(self.frequencies_cm1)
        frequencies = np.array(self.frequencies_cm1) * wavenumber
        intensities = self.notebook.settingsTab.intensities
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = self.reader.volume*angstrom*angstrom*angstrom
        # Get the cell with whole molecules from the analysis tab
        self.unit_cell = self.notebook.analysisTab.cell_of_molecules
        if self.debug:
            self.unit_cell.printInfo()
        self.normal_modes = self.notebook.analysisTab.new_normal_modes
        self.bonds = self.unit_cell.bonds
        self.snapshot_number = 0
        self.nbonds = len(self.bonds)
        #
        self.XYZ = self.unit_cell.xyz_coordinates
        self.natoms = len(self.XYZ)
        self.number_of_modes = len(self.normal_modes)
        # get the cell edges for the bounding box
        self.cell_corners,self.cell_edges = self.unit_cell.getBoundingBox()
        self.element_names = self.unit_cell.element_names
        self.colours = []
        self.radii = []
        # set the colours (the colour is a triplet of floats from 0.0 to 1.0)
        for el in self.element_names:
           colour_list = [ colour/255.0 for colour in elemental_colours[el] ] 
           colour_tuple = tuple(colour_list)
           self.colours.append( colour_tuple )
           self.radii.append(self.settings['Atom scaling']*covalent_radii[el])
        # reorder the displacement info in the normal modes into U,V and W lists 
        self.UVW = []
        for mode,displacements in enumerate(self.normal_modes):
            uvw = []
            for i in range(0,len(displacements),3):
                uvw.append( displacements[i:i+3] )
            self.UVW.append(uvw)
        # CalculatePhasePositions stores all the sphere and bond information
        self.calculatePhasePositions()
        # Add the arrows
        maxR = np.max( np.abs(np.array(self.UVW[self.selected_mode])))
        arrow_scaling = self.settings['Maximum displacement'] / maxR
        self.opengl_widget.deleteArrows()
        for uvw in self.UVW[self.selected_mode]:
            self.opengl_widget.addArrows( uvw, arrow_scaling )
        self.opengl_widget.setRotationCentre(self.unit_cell.calculateCentreOfMass() )
        self.opengl_widget.setImageSize()
        self.dirty = False
        return

    def calculatePhasePositions(self):
        # we need the number of phase steps to be odd
        if self.settings['Number of phase steps']%2 == 0:
            self.settings['Number of phase steps'] += 1
        UVW = np.array( self.UVW[self.selected_mode] )
        debugger.print('calculate phase positions UVW', UVW)
        maxR = np.amax(np.abs(UVW))
        debugger.print('calculate phase positions maxr', maxR)
        self.scale_vibrations = self.settings['Maximum displacement'] / maxR
        debugger.print('calculate phase positions scale_vibrations', self.scale_vibrations)
        self.newXYZ       = np.zeros( (self.settings['Number of phase steps'], self.natoms, 3) )
        sign = +1.0
        phases = []
        small = 1.0E-8
        n2 = int(self.settings['Number of phase steps']/2)
        delta = 1.0 / float(n2)
        phases = np.arange(-1.0, 1.0+delta-1.0E-10, delta)
        debugger.print('Phases', phases)
        for phase_index,phase in enumerate(phases):
            self.newXYZ[phase_index] = self.XYZ+phase*self.scale_vibrations*UVW
        # end for phase_index
        #
        # Store the results in the opengl widget
        #
        self.opengl_widget.deleteSpheres()
        self.opengl_widget.deleteCylinders()
        self.opengl_widget.createArrays(len(phases))
        for phase_index, phase in enumerate(phases):
            for col, rad, xyz in zip(self.colours, self.radii, self.newXYZ[phase_index]):
                self.opengl_widget.addSphere(col, rad, xyz, phase=phase_index )
            for p in self.cell_corners:
                self.opengl_widget.addSphere(self.settings['Cell colour'], self.settings['Cell radius'], p, phase=phase_index )
            for bond in self.unit_cell.bonds:
                i,j = bond
                self.opengl_widget.addCylinder(self.settings['Bond color'], self.settings['Bond radius'], self.newXYZ[phase_index,i], self.newXYZ[phase_index,j], phase=phase_index)
            for p1,p2 in self.cell_edges:
                self.opengl_widget.addCylinder(self.settings['Cell colour'], self.settings['Cell radius'], p1, p2, phase=phase_index)
        return

    def plot(self):
        if self.reader is None:
            return
        if self.plot_type_index == 0:
            self.plot_animation()
        else:
            self.plot_arrows()

    def plot_animation(self):
        self.opengl_widget.showArrows(False)
        self.opengl_widget.update()
        self.opengl_widget.startAnimation()
        return

    def plot_arrows(self):
        self.opengl_widget.showArrows(True)
        self.opengl_widget.stopAnimation()
        self.opengl_widget.update()
        return

    def refresh(self,force=False):
        if not self.dirty and not force:
            debugger.print('refresh aborted',self.dirty,force)
            return
        debugger.print('refresh widget',force)
        self.selected_mode_sb.setValue(self.selected_mode)
        self.atom_scaling_sb.setValue(self.settings['Atom scaling'])
        self.bond_radius_sb.setValue(self.settings['Bond radius'])
        self.cell_radius_sb.setValue(self.settings['Cell radius'])
        self.maximum_displacement_sb.setValue(self.settings['Maximum displacement'])
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        for index,light in enumerate(self.light_switches):
            if light:
                string = 'switch light {} off'.format(index)
            else:
                string = 'switch light {} on'.format(index)
            self.light_switches_cb.setItemText(index,string)
        self.calculate()
        self.plot()
        return

