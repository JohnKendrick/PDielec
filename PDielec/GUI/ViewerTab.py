import os
import math
import numpy as np
import PDielec.Calculator as Calculator
from collections import deque
from PyQt5.QtWidgets         import  QPushButton, QWidget
from PyQt5.QtWidgets         import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets         import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets         import  QSpinBox, QTabWidget, QDoubleSpinBox
from PyQt5.QtWidgets         import  QSizePolicy, QColorDialog, QMessageBox, QApplication
from PyQt5.QtCore            import  Qt
from PDielec.Constants        import  wavenumber, amu, PI, avogadro_si, angstrom
from PDielec.Constants        import  elemental_colours
# Import plotting requirements
from PDielec.Utilities        import Debug
from PDielec.GUI.OpenGLWidget import OpenGLWidget

class ViewerTab(QWidget):
    def __init__(self, parent, debug=False ):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'ViewerTab')
        self.debug = debug
        self.dirty = True
        self.setWindowTitle('Viewer')
        self.selected_mode = 0
        self.settings = {}
        self.UVW = deque()
        self.settings['Atom scaling'] = 0.5
        self.settings['Maximum displacement']  = 1.0
        self.settings['Bond colour']           = [  80,  80,  80, 255 ]
        self.settings['Bond radius']           = 0.1
        self.settings['Cell colour']           = [ 255,   0,   0, 255 ]
        self.settings['Cell radius']           = 0.1
        self.settings['Background colour']     = [ 120, 120, 120, 255 ]
        self.settings['Arrow colour']          = [   0, 255,   0, 255 ]
        self.settings['Arrow radius']          = 0.07
        self.settings['Number of phase steps'] = 41
        self.light_switches = [False]*8
        self.light_switches[0] = True
        self.light_switches[1] = True
        self.plot_types = ['Animation','Arrows']
        self.plot_type_index = 1
        self.number_of_molecules = 0
        self.unit_cell = None
        self.cell_edges = None
        self.cell_corners = None
        # element_colours is a dictionary
        self.element_colours = elemental_colours
        self.element_names = []
        self.species       = []
        self.element_coloured_buttons = []
        self.bond_cell_background_arrow_buttons = []
        self.bond_cell_background_arrow_names = ['Background','Cell','Bonds','Arrows']
        self.image_filename = ''
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
        # Add frequency of mode
        self.frequency_le = QLineEdit(self)
        self.frequency_le.setEnabled(False)
        self.frequency_le.setText('{}'.format(0.0))
        label = QLabel('Frequency (cm-1)', self)
        form.addRow(label, self.frequency_le)
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
        self.cell_radius_sb.setToolTip('Determines the radius of the cell outline')
        self.cell_radius_sb.valueChanged.connect(self.on_cell_radius_changed)
        #
        # The arrow radius
        #
        self.arrow_radius_sb = QDoubleSpinBox(self)
        self.arrow_radius_sb.setRange(0.01,1.0)
        self.arrow_radius_sb.setSingleStep(0.1)
        self.arrow_radius_sb.setDecimals(2)
        self.arrow_radius_sb.setValue(self.settings['Arrow radius'])
        self.arrow_radius_sb.setToolTip('Determines the radius of the arrow outline')
        self.arrow_radius_sb.valueChanged.connect(self.on_arrow_radius_changed)
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
        self.light_switches_cb.activated.connect(self.on_light_switches_cb_activated)
        #
        # Colours list of buttons with element colours
        #
        self.coloured_elements_widget = QWidget(self)
        self.element_coloured_buttons = []
        self.element_coloured_hbox = QHBoxLayout()
        for el in self.species:
            r,g,b,a = self.element_colours[el]
            button = QPushButton(el)
            button.setStyleSheet('background-color:rgba( {}, {}, {}, {});'.format(r,g,b,a))
            button.clicked.connect(self.on_coloured_element_clicked)
            self.element_coloured_buttons.append(button)
            self.element_coloured_hbox.addWidget(button)
        self.coloured_elements_widget.setLayout(self.element_coloured_hbox)
        #
        # Colours list of buttons with colours
        #
        self.coloured_buttons_widget = QWidget(self)
        self.bond_cell_background_arrow_buttons = []
        hbox = QHBoxLayout()
        button_colours = [  (r,g,b,a) for (r,g,b,a) in [self.settings['Background colour'], self.settings['Cell colour'], self.settings['Bond colour'],self.settings['Arrow colour']] ]
        self.bond_cell_background_arrow_buttons = []
        for col, name in zip(button_colours, self.bond_cell_background_arrow_names):
            button = QPushButton(name)
            button.setStyleSheet('background-color:rgba( {}, {}, {}, {});'.format(*col))
            button.clicked.connect(self.on_coloured_button_clicked)
            self.bond_cell_background_arrow_buttons.append(button)
            hbox.addWidget(button)
        self.coloured_buttons_widget.setLayout(hbox)
        #
        # Add a tab widget for the settings
        #
        self.settingsTab = QTabWidget(self)
        #self.settingsTab.currentChanged.connect(self.on_settingsTab_currentChanged)
        self.settingsTab.addTab(self.coloured_elements_widget, 'Elements')
        self.settingsTab.addTab(self.coloured_buttons_widget, 'Colours')
        self.settingsTab.addTab(self.atom_scaling_sb, 'Atom Scaling')
        self.settingsTab.addTab(self.light_switches_cb, 'Lighting')
        self.settingsTab.addTab(self.maximum_displacement_sb, 'Maximum Displacement')
        self.settingsTab.addTab(self.bond_radius_sb, 'Bond Radius')
        self.settingsTab.addTab(self.cell_radius_sb, 'Cell Radius')
        self.settingsTab.addTab(self.arrow_radius_sb, 'Arrow Radius')
        label = QLabel('Settings', self)
        form.addRow(label,self.settingsTab)
        #
        # Add buttons to save a png or movie
        #
        hbox = QHBoxLayout()
        self.filename_le = QLineEdit(self)
        self.filename_le.setToolTip('Give a file name in which to save the image.\nIf a return is pressed the filename is scanned to see whatyoe of file has been requested.')
        self.filename_le.setText(self.image_filename)
        self.filename_le.textChanged.connect(self.on_filename_le_changed)
        self.filename_le.returnPressed.connect(self.on_filename_le_return)
        self.png_filename_button = QPushButton('Save as png')
        self.png_filename_button.clicked.connect(self.on_filename_button_clicked)
        self.png_filename_button.setToolTip('Save the image as a png file with arrows')
        self.mp4_filename_button = QPushButton('Save as mp4')
        self.mp4_filename_button.clicked.connect(self.on_filename_button_clicked)
        self.mp4_filename_button.setToolTip('Save the animation as an mp4/avi file')
        self.gif_filename_button = QPushButton('Save as gif')
        self.gif_filename_button.clicked.connect(self.on_filename_button_clicked)
        self.gif_filename_button.setToolTip('Save the animation as a gif file')
        label = QLabel('Image file name', self)
        label.setToolTip('Give a file name in which to save the image.\nIf a return is pressed the filename is scanned to see whatyoe of file has been requested.')
        hbox.addWidget(self.filename_le)
        hbox.addWidget(self.png_filename_button)
        hbox.addWidget(self.mp4_filename_button)
        hbox.addWidget(self.gif_filename_button)
        form.addRow(label, hbox)
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

    def on_filename_le_return(self):
        debugger.print('on filename le return pressed')
        self.on_filename_button_clicked(True)

    def on_filename_le_changed(self,text):
        debugger.print('on filename le changed', text)
        self.image_filename = text

    def on_filename_button_clicked(self,boolean):
        debugger.print('on filename button clicked')
        #button = self.sender()
        #text = button.text()
        #
        # Check to see if the filename has been set
        #
        filename = self.image_filename
        root,extension = os.path.splitext(filename)
        valid_extension = ( extension == '.mp4' or extension == '.avi' or extension == '.png' or extension == '.gif' )
        if filename == '' or not valid_extension:
            debugger.print('Aborting on filename button clicked', filename)
            debugger.print('Aborting on filename button clicked', valid_extension)
            debugger.print('Aborting on filename button clicked', extension)
            QMessageBox.about(self,'Image file name', 'The file name for the image is not valid '+filename)
            return
        #
        # The filename will be opened in the directory of the output file which was read in mainTab
        #
        filename = os.path.join(self.notebook.mainTab.directory, filename)
        #
        # We need to remember the visualisation settings
        #
        old_plot_type = self.plot_type_index
        QApplication.setOverrideCursor(Qt.WaitCursor)
        if extension == '.png':
            self.opengl_widget.snapshot(filename)
        elif extension == '.avi':
            self.opengl_widget.save_movie(filename)
        elif extension == '.mp4':
            self.opengl_widget.save_movie(filename)
        elif extension == '.gif':
            self.opengl_widget.save_movie(filename)
        self.plot_type_index = old_plot_type
        QApplication.restoreOverrideCursor()
        return

    def on_coloured_element_clicked(self,boolean):
        debugger.print('on coloured elements clicked')
        button = self.sender()
        text = button.text()
        colour = QColorDialog.getColor(options=QColorDialog.DontUseNativeDialog)
        rgba = [ colour.red(), colour.green(), colour.blue(), colour.alpha() ]
        self.element_colours[text] = rgba
        self.dirty = True
        self.refresh()

    def on_coloured_button_clicked(self,boolean):
        debugger.print('on coloured button clicked')
        button = self.sender()
        text = button.text()
        colour = QColorDialog.getColor(options=QColorDialog.DontUseNativeDialog)
        rgba = [ colour.red(), colour.green(), colour.blue(), colour.alpha() ]
        if text == 'Background':
            self.settings['Background colour'] = rgba
        elif text == 'Cell':
            self.settings['Cell colour'] = rgba
        elif text == 'Bonds':
            self.settings['Bond colour'] = rgba
        elif text == 'Arrows':
            self.settings['Arrow colour'] = rgba
        self.dirty = True
        self.refresh()


    def on_light_switches_cb_activated(self, index):
        debugger.print('on_light_switches_cb_activated')
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

    def on_arrow_radius_changed(self,value):
        debugger.print('on arrow_radius changed ', value)
        self.settings['Arrow radius'] = value
        self.calculate()
        self.plot()

    def on_cell_radius_changed(self,value):
        debugger.print('on cell_radius changed ', value)
        self.settings['Cell radius'] = value
        self.calculate()
        self.plot()

    def on_bond_radius_changed(self,value):
        debugger.print('on_bond_radius_changed')
        self.settings['Bond radius'] = value
        self.calculate()
        self.plot()

    def on_selected_mode_changed(self):
        debugger.print('on_selected_mode_changed')
        self.selected_mode = self.selected_mode_sb.value()
        debugger.print('on selected_mode change ', self.selected_mode)
        self.frequency_le.setText('{:.5f}'.format(self.notebook.settingsTab.frequencies_cm1[self.selected_mode]))
        self.opengl_widget.deleteArrows()
        maxR = np.max( np.abs(np.array(self.UVW[self.selected_mode])))
        arrow_scaling = self.settings['Maximum displacement'] / maxR
        for uvw in self.UVW[self.selected_mode]:
            self.opengl_widget.addArrows( self.settings['Arrow colour'],self.settings['Arrow radius'], uvw, arrow_scaling )
        self.calculate()
        self.plot()

    def on_plottype_cb_changed(self, index):
        debugger.print('on_plottype_cb_changed')
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
        if program == '':
            return
        if filename == '':
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        self.selected_mode_sb.setRange(0,len(self.frequencies_cm1)-1)
        frequencies = np.array(self.frequencies_cm1) * wavenumber
        intensities = self.notebook.settingsTab.intensities
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = self.reader.volume*angstrom*angstrom*angstrom
        # Get the cell with whole molecules from the analysis tab
        self.unit_cell = self.notebook.analysisTab.cell_of_molecules
        # if self.debug:
        #     self.unit_cell.printInfo()
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
        self.species = self.reader.getSpecies()
        covalent_radii = self.notebook.analysisTab.element_radii
        self.radii = [self.settings['Atom scaling']*covalent_radii[el] for el in self.element_names ]
        self.colours = [ self.element_colours[el] for el in self.element_names ]
        # reorder the displacement info in the normal modes into U,V and W lists
        # Using deque here rather than a simple list as the memory allocation doesn't have to be contiguous
        self.UVW.clear()
        # self.UVW = []
        for mode,displacements in enumerate(self.normal_modes):
            uvw = deque()
            #uvw = []
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
            self.opengl_widget.addArrows( self.settings['Arrow colour'],self.settings['Arrow radius'], uvw, arrow_scaling )
        self.opengl_widget.setRotationCentre(self.unit_cell.calculateCentreOfMass() )
        self.opengl_widget.setImageSize()
        self.dirty = False
        QApplication.restoreOverrideCursor()
        return

    def setColour(self, element, colour):
        debugger.print('setcolour')
        if element == 'Background' or element == 'background':
            self.settings['Background colour'] = colour
        elif element == 'Cell' or element == 'cell':
            self.settings['Cell colour'] = colour
        elif element == 'Arrow' or element == 'arrow':
            self.settings['Arrow colour'] = colour
        else:
            self.element_colours[element] = colour
        self.plot()

    def calculatePhasePositions(self):
        debugger.print('calculatePhasePositions')
        # we need the number of phase steps to be odd
        if self.settings['Number of phase steps']%2 == 0:
            self.settings['Number of phase steps'] += 1
        UVW = np.array( self.UVW[self.selected_mode] )
        maxR = np.amax(np.abs(UVW))
        self.scale_vibrations = self.settings['Maximum displacement'] / maxR
        self.newXYZ       = np.zeros( (self.settings['Number of phase steps'], self.natoms, 3) )
        sign = +1.0
        small = 1.0E-8
        n2 = int(self.settings['Number of phase steps']/2)
        delta = 1.0 / float(n2)
        phases = np.arange(-1.0, 1.0+delta-1.0E-10, delta)
        for phase_index,phase in enumerate(phases):
            self.newXYZ[phase_index] = self.XYZ+phase*self.scale_vibrations*UVW
        # end for phase_index
        #
        # Store the results in the opengl widget
        #
        self.opengl_widget.deleteSpheres()
        self.opengl_widget.deleteCylinders()
        self.opengl_widget.createArrays(len(phases))
        debugger.print('calculatePhasePositions - adding spheres and cylinders')
        for phase_index, phase in enumerate(phases):
            for col, rad, xyz in zip(self.colours, self.radii, self.newXYZ[phase_index]):
                self.opengl_widget.addSphere(col, rad, xyz, phase=phase_index )
            for p in self.cell_corners:
                self.opengl_widget.addSphere(self.settings['Cell colour'], self.settings['Cell radius'], p, phase=phase_index )
            for bond in self.unit_cell.bonds:
                i,j = bond
                self.opengl_widget.addCylinder(self.settings['Bond colour'], self.settings['Bond radius'], self.newXYZ[phase_index,i], self.newXYZ[phase_index,j], phase=phase_index)
            for p1,p2 in self.cell_edges:
                self.opengl_widget.addCylinder(self.settings['Cell colour'], self.settings['Cell radius'], p1, p2, phase=phase_index)
        debugger.print('calculatePhasePositions - exiting')
        return

    def plot(self):
        debugger.print('plot')
        if self.reader is None:
            return
        if self.plot_type_index == 0:
            self.plot_animation()
        else:
            self.plot_arrows()

    def plot_animation(self):
        debugger.print('plot_animation')
        self.opengl_widget.showArrows(False)
        self.opengl_widget.update()
        self.opengl_widget.startAnimation()
        return

    def plot_arrows(self):
        debugger.print('plot_arrows')
        self.opengl_widget.showArrows(True)
        self.opengl_widget.stopAnimation()
        self.opengl_widget.update()
        return

    def refresh(self,force=False):
        debugger.print('refresh')
        if not self.dirty and not force and not self.notebook.visualerCalculationRequired:
            debugger.print('refresh aborted',self.dirty,force)
            return
        debugger.print('refresh widget',force)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        self.selected_mode_sb.setValue(self.selected_mode)
        try:
            self.frequency_le.setText('{:.5f}'.format(self.notebook.settingsTab.frequencies_cm1[self.selected_mode]))
        except:
            self.frequency_le.setText('{:.5f}'.format(0.0))
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
        #
        # Colours list of buttons with element colours
        #
        self.reader = self.notebook.mainTab.reader
        if self.reader is None:
            debugger.print('reader is none ')
            self.element_names = []
            self.species = []
        else:
            self.unit_cell = self.notebook.analysisTab.cell_of_molecules
            self.element_names = self.unit_cell.element_names
            self.species = self.reader.getSpecies()
            debugger.print('refresh - species ',self.species)
        count = self.element_coloured_hbox.count()
        if count == 0:
            debugger.print('initialising element colours widget')
            self.element_coloured_buttons = []
            for el in self.species:
                r,g,b,a = self.element_colours[el]
                button = QPushButton(el)
                button.setStyleSheet('background-color:rgba( {}, {}, {}, {});'.format(r,g,b,a))
                button.clicked.connect(self.on_coloured_element_clicked)
                self.element_coloured_buttons.append(button)
                self.element_coloured_hbox.addWidget(button)
        else:
            debugger.print('update element colours widget')
            for el,button in zip(self.species,self.element_coloured_buttons):
                r,g,b,a = self.element_colours[el]
                button.setStyleSheet('background-color:rgba( {}, {}, {}, {});'.format(r,g,b,a))
        #
        # Colours list of buttons with colours
        #
        button_colours = [  (r,g,b,a) for (r,g,b,a) in [self.settings['Background colour'], self.settings['Cell colour'], self.settings['Bond colour'],self.settings['Arrow colour']] ]
        for col, button in zip(button_colours, self.bond_cell_background_arrow_buttons):
            button.setStyleSheet('background-color:rgba( {}, {}, {}, {});'.format(*col))
        #
        # Filename
        #
        self.filename_le.setText(self.image_filename)
        #
        # Unlock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.calculate()
        self.plot()
        QApplication.restoreOverrideCursor()
        return
