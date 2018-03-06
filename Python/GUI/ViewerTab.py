import math
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets         import  QPushButton, QWidget
from PyQt5.QtWidgets         import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets         import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets         import  QSpinBox
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
        self.subplot = None
        self.setWindowTitle('Viewer')
        self.selected_mode = 4
        self.atom_scaling = 0.3
        self.maxdisplacement = 1.0
        self.plot_types = ['Animation','Arrows']
        self.plot_type_index = 0
        self.number_of_molecules = 0
        self.bond_colour = [ 0.3, 0.3, 0.3 ]
        self.bond_radius = 0.1
        self.cell_edge_colour = [ 1.0, 0.0, 0.0 ]
        self.cell_edge_radius = 0.1
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
        self.atom_scaling_le = QLineEdit(self)
        self.atom_scaling_le.setText('{}'.format(self.atom_scaling))
        self.atom_scaling_le.setToolTip('Scale the covalent radii to determine the size of the atomic spheres')
        self.atom_scaling_le.textChanged.connect(self.on_atom_scaling_changed)
        label = QLabel('Atom scaling', self)
        label.setToolTip('Scale the covalent radii to determine the size of the atomic spheres')
        form.addRow(label, self.atom_scaling_le)
        #
        # The bond width
        #
        self.bond_radius_le = QLineEdit(self)
        self.bond_radius_le.setText('{}'.format(self.bond_radius))
        self.bond_radius_le.setToolTip('Determines the size of the bonds drawn')
        self.bond_radius_le.textChanged.connect(self.on_bond_radius_changed)
        label = QLabel('Bond radius', self)
        label.setToolTip('Determines the size of the bonds drawn')
        form.addRow(label, self.bond_radius_le)
        #
        #
        # The cell width
        #
        self.cell_radius_le = QLineEdit(self)
        self.cell_radius_le.setText('{}'.format(self.cell_edge_radius))
        self.cell_radius_le.setToolTip('Determines the size of the cell outline')
        self.cell_radius_le.textChanged.connect(self.on_cell_radius_changed)
        label = QLabel('Cell radius', self)
        label.setToolTip('Determines the size of the cell outline')
        form.addRow(label, self.cell_radius_le)
        #
        # The maximum displacement
        #
        self.maximum_displacement_le = QLineEdit(self)
        self.maximum_displacement_le.setText('{}'.format(self.maxdisplacement))
        self.maximum_displacement_le.setToolTip('Set the size of the maximum displacement')
        self.maximum_displacement_le.textChanged.connect(self.on_maximum_displacement_changed)
        label = QLabel('Maximum displacement', self)
        label.setToolTip('Set the size of the maximum displacement')
        form.addRow(label, self.maximum_displacement_le)
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
        # Add the opengl figure to the bottom 
        #
        self.opengl_widget = OpenGLWidget(self)
        self.opengl_widget.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding))
        form.addRow(self.opengl_widget)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)

    def on_maximum_displacement_changed(self,text):
        debugger.print('on maximum_displacement changed ', text)
        try:
          self.maxdisplacement = float(text)
        except:
          pass
        self.maximum_displacement_le.blockSignals(True)
        if self.maxdisplacement < 0.0:
            self.maxdisplacement = 0.0
            self.maximum_displacement_le.setText('{}'.format(self.maxdisplacement))
        if self.maxdisplacement > 10.0:
            self.maxdisplacement = 10.0
            self.maximum_displacement_le.setText('{}'.format(self.maxdisplacement))
        self.maximum_displacement_le.blockSignals(False)
        self.plot()
        
    def on_atom_scaling_changed(self,text):
        debugger.print('on atom_scaling_le changed ', text)
        try:
            self.atom_scaling = float(text)
        except:
          pass
        self.atom_scaling_le.blockSignals(True)
        if self.atom_scaling < 0.0:
            self.atom_scaling = 0.0
            self.atom_scaling_le.setText('{}'.format(self.atom_scaling))
        if self.atom_scaling > 10.0:
            self.atom_scaling = 10.0
            self.atom_scaling_le.setText('{}'.format(self.atom_scaling))
        self.atom_scaling_le.blockSignals(False)
        self.refresh()
        
    def on_cell_radius_changed(self,text):
        debugger.print('on cell_radius changed ', text)
        try:
            self.cell_edge_radius = float(text)
        except:
          pass
        self.cell_radius_le.blockSignals(True)
        if self.cell_edge_radius < 0.0:
            self.cell_edge_radius = 0.0
            self.cell_radius_le.setText('{}'.format(self.cell_edge_radius))
        if self.cell_edge_radius > 2.0:
            self.cell_edge_radius = 2.0
            self.cell_radius_le.setText('{}'.format(self.cell_edge_radius))
        self.cell_radius_le.blockSignals(False)
        self.refresh()
        
    def on_bond_radius_changed(self,text):
        debugger.print('on bond_radius_le changed ', text)
        try:
            self.bond_radius = float(text)
        except:
          pass
        self.bond_radius_le.blockSignals(True)
        if self.bond_radius < 0.0:
            self.bond_radius = 0.0
            self.bond_radius_le.setText('{}'.format(self.bond_radius))
        if self.bond_radius > 2.0:
            self.bond_radius = 2.0
            self.bond_radius_le.setText('{}'.format(self.bond_radius))
        self.bond_radius_le.blockSignals(False)
        self.refresh()
        
    def on_tolerance_changed(self,text):
        debugger.print('on file_tolerance_le changed ', text)
        try:
            self.bond_tolerance = float(text)
        except:
          pass
        self.tolerance_le.blockSignals(True)
        if self.bond_tolerance < 0.0:
            self.bond_tolerance = 0.0
            self.tolerance_le.setText('{}'.format(self.bond_tolerance))
        if self.bond_tolerance > 2.0:
            self.bond_tolerance = 2.0
            self.tolerance_le.setText('{}'.format(self.bond_tolerance))
        self.tolerance_le.blockSignals(False)
        self.refresh()

    def on_selected_mode_changed(self):
        self.selected_mode = self.selected_mode_sb.value()
        debugger.print('on selected_mode change ', self.selected_mode)
        self.opengl_widget.deleteDisplacements()
        for x,y,z in zip( self.U[self.selected_mode], self.V[self.selected_mode], self.W[self.selected_mode] ):
            self.opengl_widget.addDisplacement( [x,y,z] )
        self.plot()

    def refresh(self):
        debugger.print('refreshing widget')
        self.calculate()
        self.plot()
        return

    def on_plottype_cb_changed(self, index):
        self.plot_type_index = index
        debugger.print('Plot type index changed to ', self.plot_type_index)
        self.plot()

    def calculate(self):
        debugger.print('calculate')
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['program']
        filename = settings['filename']
        self.reader = self.notebook.mainTab.reader
        if self.reader is None:
            return
        if program is '':
            return
        if filename is '':
            return
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        modes_selected = self.notebook.settingsTab.modes_selected
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        mass_weighted_normal_modes = self.notebook.settingsTab.mass_weighted_normal_modes
        frequencies = np.array(self.frequencies_cm1) * wavenumber
        intensities = self.notebook.settingsTab.intensities
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = self.reader.volume*angstrom*angstrom*angstrom
        # Find the last unit cell read by the reader and its masses
        cell = self.reader.unit_cells[-1]
        atom_masses = self.reader.masses
        cell.set_atomic_masses(atom_masses)
        scale = 1.1
        tolerance = 0.2
        newcell,nmols,old_order = cell.calculate_molecular_contents(scale, tolerance, covalent_radii)
        # get the normal modes from the mass weighted ones
        self.normal_modes = Calculator.normal_modes(atom_masses, mass_weighted_normal_modes)
        # Reorder the atoms so that the mass weighted normal modes order agrees with the ordering in the new cell
        nmodes,nions,temp = np.shape(self.normal_modes)
        new_normal_modes = np.zeros( (nmodes,3*nions) )
        new_mass_weighted_normal_modes = np.zeros( (nmodes,3*nions) )
        masses = newcell.atomic_masses
        for imode,mode in enumerate(mass_weighted_normal_modes):
            for index,old_index in enumerate(old_order):
                i = index*3
                j = old_index*3
                new_mass_weighted_normal_modes[imode,i+0] = mode[old_index][0] 
                new_mass_weighted_normal_modes[imode,i+1] = mode[old_index][1] 
                new_mass_weighted_normal_modes[imode,i+2] = mode[old_index][2] 
                new_normal_modes[imode,i+0] = new_mass_weighted_normal_modes[imode,i+0] / math.sqrt(masses[index])
                new_normal_modes[imode,i+1] = new_mass_weighted_normal_modes[imode,i+1] / math.sqrt(masses[index])
                new_normal_modes[imode,i+2] = new_mass_weighted_normal_modes[imode,i+2] / math.sqrt(masses[index])

        self.normal_modes = new_normal_modes
        self.snapshot_number = 0
        self.unit_cell = newcell
        self.bonds = self.unit_cell.bonds
        self.nbonds = len(self.bonds)
        #
        self.X = self.unit_cell.xyz_coordinates[:,0]
        self.Y = self.unit_cell.xyz_coordinates[:,1]
        self.Z = self.unit_cell.xyz_coordinates[:,2]
        self.natoms = len(self.X)
        self.number_of_modes = len(self.normal_modes)
        self.number_of_modesm1 = len(self.normal_modes)-1
        # reorder the displacement info in the normal modes into U,V and W lists
        self.U, self.V, self.W = [], [], []
        for mode,displacements in enumerate(self.normal_modes):
            U = np.zeros_like( self.X )
            V = np.zeros_like( self.Y )
            W = np.zeros_like( self.Z )
            index = -1
            for i in range(0,len(displacements),3):
                index += 1
                U[index] = displacements[i + 0]
                V[index] = displacements[i + 1]
                W[index] = displacements[i + 2]
            self.U.append(U)
            self.V.append(V)
            self.W.append(W)
        # get the cell edges for the bounding box
        self.cell_edges = self.unit_cell.getBoundingBox()
        self.element_names = self.unit_cell.element_names
        self.colours = []
        self.radii = []
        # set the colours (the colour is a triplet of floats from 0.0 to 1.0)
        for el in self.element_names:
           colour_list = [ colour/255.0 for colour in elemental_colours[el] ] 
           colour_tuple = tuple(colour_list)
           self.colours.append( colour_tuple )
           self.radii.append(self.atom_scaling*covalent_radii[el])
        self.mode = 0
        self.number_of_phase_steps = 21
        self.phase_direction = 1
        self.phase_index = 0
        self.calculatePhasePositions()
        # Give the opengl_widget the information it needs
        self.opengl_widget.deleteSpheres()
        for col, rad, xyz in zip(self.colours, self.radii, self.unit_cell.xyz_coordinates):
            self.opengl_widget.addSphere(col, rad, xyz )
        self.opengl_widget.deleteCylinders()
        for bond in self.unit_cell.bonds:
            i,j = bond
            start = [ self.X[i], self.Y[i], self.Z[i] ]
            end   = [ self.X[j], self.Y[j], self.Z[j] ]
            self.opengl_widget.addCylinder(self.bond_colour, self.bond_radius, start, end)
        self.opengl_widget.deleteDisplacements()
        for u,v,w in zip( self.U[self.selected_mode], self.V[self.selected_mode], self.W[self.selected_mode] ):
            self.opengl_widget.addDisplacement( [u, v, w] )
        self.opengl_widget.setRotationCentre(self.unit_cell.calculateCentreOfMass() )
        for p1,p2 in self.cell_edges:
            self.opengl_widget.addCylinder(self.cell_edge_colour, self.cell_edge_radius, p1, p2)
        self.opengl_widget.setImageSize()
        return

    def calculatePhasePositions(self):
        # we need the number of phase steps to be odd
        if self.number_of_phase_steps%2 == 0:
            self.numbers_of_phase_steps += 1
        U = self.U[self.mode]
        V = self.V[self.mode]
        W = self.W[self.mode]
        return

    def plot(self):
         if self.reader is None:
             return
         if self.plot_type_index == 0:
             self.plot_animation()
         else:
             self.plot_arrows()

    def plot_animation(self):
        self.opengl_widget.update()

    def plot_arrows(self):
        self.opengl_widget.update()

