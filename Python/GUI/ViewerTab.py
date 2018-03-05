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
        self.title = 'Viewer'
        self.atom_scaling = 0.3
        self.maxdisplacement = 1.0
        self.cell_width = 0.1
        self.plot_types = ['Animation','Arrows']
        self.plot_type_index = 0
        self.number_of_molecules = 0
        self.bond_colour = [ 0.2, 0.2, 0.2 ]
        self.bond_radius = 0.1
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
        label = QLabel('Bond width', self)
        label.setToolTip('Determines the size of the bonds drawn')
        form.addRow(label, self.bond_radius_le)
        #
        #
        # The cell width
        #
        self.cell_width_le = QLineEdit(self)
        self.cell_width_le.setText('{}'.format(self.cell_width))
        self.cell_width_le.setToolTip('Determines the size of the cell outline')
        self.cell_width_le.textChanged.connect(self.on_cell_width_changed)
        label = QLabel('Cell width', self)
        label.setToolTip('Determines the size of the cell outline')
        form.addRow(label, self.cell_width_le)
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
        # Set the plot title         
        #
        self.title_le = QLineEdit(self) 
        self.title_le.setToolTip('Set the plot title')
        self.title_le.setText(self.title)
        self.title_le.textChanged.connect(self.on_title_changed)
        label = QLabel('Plot title', self)
        label.setToolTip('Set the plot title')
        form.addRow(label, self.title_le)
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
        
    def on_cell_width_changed(self,text):
        debugger.print('on cell_width_le changed ', text)
        try:
            self.cell_width = float(text)
        except:
          pass
        self.cell_width_le.blockSignals(True)
        if self.cell_width < 0.0:
            self.cell_width = 0.0
            self.cell_width_le.setText('{}'.format(self.cell_width))
        if self.cell_width > 2.0:
            self.cell_width = 2.0
            self.cell_width_le.setText('{}'.format(self.cell_width))
        self.cell_width_le.blockSignals(False)
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

    def on_file_store_le_changed(self,text):
        self.spreadsheet = text
        debugger.print('on file_store_le change ', self.spreadsheet)

    def on_title_changed(self,text):
        self.title = text
        if self.subplot is not None:
            self.subplot.set_title(self.title)
            self.canvas.draw_idle()
        debugger.print('on title change ', self.title)

    def on_selected_mode_changed(self):
        self.selected_mode = self.selected_mode_sb.value()
        debugger.print('on selected_mode change ', self.selected_mode)
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
        print('rot centre newcell', newcell.calculateCentreOfMass() )
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
        for col, rad, x, y, z in zip(self.colours, self.radii, self.X, self.Y, self.Z):
            self.opengl_widget.addSphere(col, rad, x, y, z)
        for bond in self.unit_cell.bonds:
            i,j = bond
            start = [ self.X[i], self.Y[i], self.Z[i] ]
            end   = [ self.X[j], self.Y[j], self.Z[j] ]
            self.opengl_widget.addCylinder(self.bond_colour, self.bond_radius, start, end)
        print('rot centre', self.unit_cell.calculateCentreOfMass() )
        self.opengl_widget.setRotationCentre(self.unit_cell.calculateCentreOfMass() )
        return

    def calculatePhasePositions(self):
        # we need the number of phase steps to be odd
        if self.number_of_phase_steps%2 == 0:
            self.numbers_of_phase_steps += 1
        U = self.U[self.mode]
        V = self.V[self.mode]
        W = self.W[self.mode]
        maxR = np.amax(U) 
        maxR = max(maxR, np.amax(V) )
        maxR = max(maxR, np.amax(W))
        self.scale_vibrations = self.maxdisplacement / maxR
        #print('MaxR is',maxR)
        self.newX    = np.zeros( (self.number_of_phase_steps,self.natoms) )
        self.newY    = np.zeros( (self.number_of_phase_steps,self.natoms) )
        self.newZ    = np.zeros( (self.number_of_phase_steps,self.natoms) )
        self.startX  = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.startY  = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.startZ  = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.dirX    = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.dirY    = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.dirZ    = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.bondx   = np.zeros( (self.number_of_phase_steps,self.natoms,2) )
        self.bondy   = np.zeros( (self.number_of_phase_steps,self.natoms,2) )
        self.bondz   = np.zeros( (self.number_of_phase_steps,self.natoms,2) )
        # work out the scaling factor to get the require maximumdisplacement
        disX = self.X - U
        sign = +1.0
        phases = []
        small = 1.0E-8
        n2 = int(self.number_of_phase_steps/2)
        delta = 1.0 / n2
        phase = -delta
        for phase_index in range(self.number_of_phase_steps):
            phase += delta*sign
            self.newX[phase_index] = self.X+U*phase*self.scale_vibrations
            self.newY[phase_index] = self.Y+V*phase*self.scale_vibrations
            self.newZ[phase_index] = self.Z+W*phase*self.scale_vibrations
            if phase >= 1.0-small:
                sign = -1.0
                phase += delta*sign
            elif phase <= -1.0+small:
                sign = +1.0
                phase += delta*sign
            phases.append(phase)
        # end for phase_index
        for phase_index,phase in enumerate(phases):
            for ibond,bond in enumerate(self.bonds):
                i,j = bond
                self.bondx[phase_index,ibond,0] = self.newX[phase_index,i]
                self.bondy[phase_index,ibond,0] = self.newY[phase_index,i]
                self.bondz[phase_index,ibond,0] = self.newZ[phase_index,i]
                self.bondx[phase_index,ibond,1] = self.newX[phase_index,j]
                self.bondy[phase_index,ibond,1] = self.newY[phase_index,j]
                self.bondz[phase_index,ibond,1] = self.newZ[phase_index,j]
                self.startX[phase_index,ibond] = self.newX[phase_index,i]
                self.startY[phase_index,ibond] = self.newY[phase_index,i]
                self.startZ[phase_index,ibond] = self.newZ[phase_index,i]
                self.dirX[phase_index,ibond] = self.newX[phase_index,j] - self.newX[phase_index,i]
                self.dirY[phase_index,ibond] = self.newY[phase_index,j] - self.newY[phase_index,i]
                self.dirZ[phase_index,ibond] = self.newZ[phase_index,j] - self.newZ[phase_index,i]
            #end for ibond,bond
        # end of phase_index,phase
        

    def animate(self):
        pass

    def refreshBonds(self, phase_index):
        pass
        return

    def refreshAtoms(self, phase_index):
        pass
        return

    def refreshDisplacements(self):
        pass
        return

    def draw(self):
        self.drawUnitCell()
        self.drawBonds()
        self.drawAtoms()
        self.drawDisplacements()

    def drawDisplacements(self):
        pass
        return

    def drawUnitCell(self):
        bondx = np.zeros( 2 )
        bondy = np.zeros( 2 )
        bondz = np.zeros( 2 )
        ones  = np.array( (1.0,1.0) )
        for edge in self.cell_edges:
            p1,p2 =edge
            bondx[0] = p1[0]
            bondy[0] = p1[1]
            bondz[0] = p1[2]
            bondx[1] = p2[0]
            bondy[1] = p2[1]
            bondz[1] = p2[2]
            #mlab.plot3d(bondx,bondy,bondz,ones, tube_radius=0.15, tube_sides=10, color=(1,0,0) )
        return

    def drawQuiverBonds(self):
        self.glyph_quiver_bonds = mlab.quiver3d(self.startX[0],self.startY[0],self.startZ[0],
                                                self.dirX[0],  self.dirY[0],  self.dirZ[0],
                                                color=(0,0,0), mode='cylinder')

    def drawBonds(self):
        bondx = np.zeros( 2 )
        bondy = np.zeros( 2 )
        bondz = np.zeros( 2 )
        ones  = np.array( (1.0,1.0) )
        self.glyph_bonds = []
        for bond in self.unit_cell.bonds:
            i,j = bond
            bondx[0] = self.X[i]
            bondy[0] = self.Y[i]
            bondz[0] = self.Z[i]
            bondx[1] = self.X[j]
            bondy[1] = self.Y[j]
            bondz[1] = self.Z[j]
            #glyph = mlab.plot3d(bondx,bondy,bondz,ones,tube_radius=0.1, tube_sides=10, color=(0,0,0) )
            #self.glyph_bonds.append(glyph)
        return

    def drawAtoms(self):
        #for x,y,z,colour,size in zip(self.X, self.Y, self.Z, self.colours, self.radii):
           #glyph = mlab.points3d([x],[y],[z],[1],color=colour,resolution=20, scale_factor=size)
           #self.glyph_atoms.append(glyph)
        return


    def plot(self):
         if self.reader is None:
             return
         if self.plot_type_index == 0:
             self.plot_animation()
         else:
             self.plot_arrows()

    def plot_animation(self):
        pass

    def plot_internal_external(self):
        pass

