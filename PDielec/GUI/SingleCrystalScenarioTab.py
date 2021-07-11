# -*- coding: utf8 -*-
import os.path
import os
import numpy as np
import copy
import PDielec.Calculator as Calculator
import PDielec.GTMcore as GTM
from PyQt5.QtWidgets            import QPushButton, QWidget
from PyQt5.QtWidgets            import QComboBox, QLabel, QLineEdit, QListWidget
from PyQt5.QtWidgets            import QApplication
from PyQt5.QtWidgets            import QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets            import QSpinBox,QDoubleSpinBox
from PyQt5.QtWidgets            import QSizePolicy
from PyQt5.QtCore               import QCoreApplication, Qt
from PDielec.DielectricFunction import DielectricFunction
from PDielec.Constants          import wavenumber, PI, avogadro_si, angstrom, speed_light_si
from PDielec.Utilities          import Debug
from PDielec.GUI.ScenarioTab    import ScenarioTab

def set_affinity_on_worker():
    '''When a new worker process is created, the affinity is set to all CPUs'''
    #JK print('I'm the process %d, setting affinity to all CPUs.' % os.getpid())
    #JK Commented out for the time being
    #JK os.system('taskset -p 0xff %d > /dev/null' % os.getpid())


class SingleCrystalScenarioTab(ScenarioTab):
    def __init__(self, parent, debug=False ):
        ScenarioTab.__init__(self,parent)
        global debugger
        debugger = Debug(debug,'SingleCrystalScenarioTab:')
        debugger.print('In the initialiser')
        self.dirty = True
        self.scenarioType = 'Single crystal'
        self.settings['Scenario type'] = 'Single crystal'
        self.settings['Unique direction - h'] = 0
        self.settings['Unique direction - k'] = 0
        self.settings['Unique direction - l'] = 1
        self.settings['Azimuthal angle'] = 0.0
        self.settings['Angle of incidence'] = 0.0
        self.settings['Superstrate dielectric'] = 1.0
        self.settings['Substrate dielectric'] = 1.0
        self.settings['Superstrate depth'] = 999.0
        self.settings['Substrate depth'] = 999.0
        self.settings['Film thickness'] = 100.0
        self.settings['Mode'] = 'Thick slab'
        self.settings['Frequency units'] = 'wavenumber'
        self.p_reflectance = []
        self.s_reflectance = []
        self.p_transmittance = []
        self.s_transmittance = []
        self.p_absorbtance = []
        self.s_absorbtance = []
        self.epsilon = []
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
        self.mode_cb.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two media; the incident and the crystal),\n Coherent thin film assumes there are three media; the incident, the crystal and the substrate')
        self.mode_cb.addItems( ['Thick slab','Coherent thin film'] )
        self.settings['Mode'] = 'Thick slab'
        self.mode_cb.activated.connect(self.on_mode_cb_activated)
        label = QLabel('Single crystal mode', self)
        label.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two media; the incident and the crystal),\n Coherent thin film assumes there are three media; the incident, the crystal and the substrate')
        form.addRow(label, self.mode_cb)
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
        label = QLabel('Azimuthal angle',self)
        label.setToolTip('Define the azimuthal angle (rotation of the crystal about the lab Z-axis)')
        form.addRow(label, self.azimuthal_angle_sb)
        self.angle_of_incidence_sb = QDoubleSpinBox(self)
        self.angle_of_incidence_sb.setToolTip('Define the angle of incidence, (normal incidence is 0 degrees)')
        self.angle_of_incidence_sb.setRange(0,90)
        self.angle_of_incidence_sb.setSingleStep(5)
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        self.angle_of_incidence_sb.valueChanged.connect(self.on_angle_of_incidence_sb_changed)
        label = QLabel('Angle of incidence',self)
        label.setToolTip('Define the angle of incidence, (normal incidence is 0 degrees).')
        form.addRow(label, self.angle_of_incidence_sb)
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
        self.film_thickness_sb.setToolTip('Define the thin film thickness in nanometres')
        self.film_thickness_sb.setRange(0,100000)
        self.film_thickness_sb.setSingleStep(1)
        self.film_thickness_sb.setValue(self.settings['Film thickness'])
        self.film_thickness_sb.valueChanged.connect(self.on_film_thickness_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.film_thickness_sb)
        film_thickness_label = QLabel('Film thickness (nm)',self)
        film_thickness_label.setToolTip('Define the depth of the thin film in nanometres.')
        form.addRow(film_thickness_label, hbox)
        #
        # Add a legend option
        #
        self.legend_le = QLineEdit(self)
        self.legend_le.setToolTip('The legend will be used to describe the results in the plot')
        self.legend_le.setText('Scenario legend')
        self.legend_le.textChanged.connect(self.on_legend_le_changed)
        label = QLabel('Scenario legend',self)
        label.setToolTip('The legend will be used to describe the results in the plot')
        form.addRow(label, self.legend_le)
        #
        #
        # Final buttons for changing/deleting and switching the scenarios
        #
        hbox = self.add_scenario_buttons()
        form.addRow(hbox)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()

    def on_film_thickness_sb_changed(self,value):
        debugger.print('on_film_thickness_sb', value)
        self.dirty = True
        self.settings['Film thickness'] = value

    def on_superstrate_dielectric_sb_changed(self,value):
        debugger.print('on_superstrate_dielectric_sb', value)
        self.dirty = True
        self.settings['Superstrate dielectric'] = value

    def on_substrate_dielectric_sb_changed(self,value):
        debugger.print('on_substrate_dielectric_sb', value)
        self.dirty = True
        self.settings['Substrate dielectric'] = value

    def on_azimuthal_angle_sb_changed(self,value):
        debugger.print('on_azimuthal_angl_sb_changed', value)
        self.dirty = True
        self.settings['Azimuthal angle'] = value
        self.calculate_euler_angles()

    def on_angle_of_incidence_sb_changed(self,value):
        debugger.print('on_angle_of_incidence_sb_changed', value)
        self.dirty = True
        self.settings['Angle of incidence'] = value

    def on_h_sb_changed(self,value):
        debugger.print('on_h_sb_changed', value)
        self.dirty = True
        self.settings['Unique direction - h'] = value
        self.calculate_euler_angles()

    def on_k_sb_changed(self,value):
        debugger.print('on_k_sb_changed', value)
        self.dirty = True
        self.settings['Unique direction - k'] = value
        self.calculate_euler_angles()

    def on_l_sb_changed(self,value):
        debugger.print('on_l_sb_changed', value)
        self.dirty = True
        self.settings['Unique direction - l'] = value
        self.calculate_euler_angles()

    def refresh(self,force=False):
        if not self.dirty and not force :
            debugger.print('refreshing widget aborted', self.dirty,force)
            return
        debugger.print('refreshing widget', force)
        # Change any greyed out items
        self.greyed_out()
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Now refresh values
        self.reader = self.notebook.mainTab.reader
        if self.reader is not None:
            self.cell = self.reader.unit_cells[-1]
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
            self.legend_le.setText(self.settings['Legend'])
            # Refresh the widgets that depend on the reader
            self.reader = self.notebook.reader
            self.calculate_euler_angles()
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        return

    def on_mode_cb_activated(self, index):
        debugger.print('on_mode_cb_activated')
        if index == 0:
            self.settings['Mode'] = 'Thick slab'
        elif index == 1:
            self.settings['Mode'] = 'Coherent thin film'
        else:
            self.settings['Mode'] = 'Incoherent thin film'
        self.dirty = True
        self.refresh()
        debugger.print('Mode changed to ', self.settings['Mode'])

    def calculate(self,vs_cm1):
        debugger.print('Calculate')
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
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        frequencies = np.array(frequencies_cm1) * wavenumber
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
        crystal          = GTM.Layer(thickness=crystalDepth*1e-9,    epsilon=crystalPermittivityFunction)
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
        # Get a pool of processors
        pool = Calculator.get_pool(self.notebook.ncpus, self.notebook.threading)
        # Prepare parameters for a parallel call to the layered absorption / reflection
        call_parameters = []
        results = []
        # Assemble all the parameters we need for parallel execution
        # About to call
        for v in vs_cm1:
            #jk system = copy.deepcopy(system)
            call_parameters.append( (v, angleOfIncidence, system) )
            # results.append( Calculator.solve_single_crystal_equations( (v, angleOfIncidence, system) ) )
        for result in pool.map(Calculator.solve_single_crystal_equations, call_parameters, chunksize=40):
            self.notebook.progressbars_update()
            results.append(result)
        QCoreApplication.processEvents()
        # Initialise plotting variables
        self.vs_cm1= []
        self.p_reflectance = []
        self.s_reflectance = []
        self.p_transmittance = []
        self.s_transmittance = []
        self.p_absorbtance = []
        self.s_absorbtance = []
        self.epsilon = []
        for v,r,R,t,T,epsilon in results:
            self.vs_cm1.append(v)
            self.p_reflectance.append(R[0]+R[2])
            self.s_reflectance.append(R[1]+R[3])
            self.p_transmittance.append(T[0])
            self.s_transmittance.append(T[1])
            self.p_absorbtance.append(R[0]+R[2]+T[0])
            self.s_absorbtance.append(R[1]+R[3]+T[1])
            self.epsilon.append(epsilon)
        # Close parallel processing down
        pool.close()
        pool.join()
        self.dirty = False
        QCoreApplication.processEvents()

    def write_spreadsheet(self):
        debugger.print('write_spreadsheet')
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
        sp.writeNextRow([ 'Surface definition (h)',       self.settings['Unique direction - h'] ],col=1)
        sp.writeNextRow([ 'Surface definition (k)',       self.settings['Unique direction - k'] ],col=1)
        sp.writeNextRow([ 'Surface definition (l)',       self.settings['Unique direction - l'] ],col=1)
        sp.writeNextRow([ 'Azimuthal angle',              self.settings['Azimuthal angle'] ],col=1)
        sp.writeNextRow([ 'Angle of incidence',           self.settings['Angle of incidence'] ],col=1)
        sp.writeNextRow([ 'Superstrate dielectric',       self.settings['Superstrate dielectric'] ],col=1)
        sp.writeNextRow([ 'Substrate dielectric',         self.settings['Substrate dielectric'] ],col=1)
        sp.writeNextRow([ 'Film thickness(nm)',           self.settings['Film thickness'] ],col=1)
        headings = ['R_p', 'R_s', 'T_p', 'T_s']
        self.write_results(sp, 'Crystal R&T',     self.vs_cm1, [self.p_reflectance, self.s_reflectance, self.p_transmittance, self.s_transmittance], headings)

    def write_results(self, sp, name, vs, yss, headings):
        """ 
        sp        is the spreadsheet object
        name      is the worksheet name used for writing
        vs        an np.array of the frequencies
        yss       a list of np.arrays of the reflections and transmittance ] 
        headings  the heading names for the yss
        """
        debugger.print('write_results')
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

    def get_result(self, vs_cm1, plot_type):
        """Return a particular result"""
        self.get_results(vs_cm1)
        if plot_type   == 'Crystal Reflectance (P polarisation)':
            return self.p_reflectance
        elif plot_type == 'Crystal Reflectance (S polarisation)':
            return self.s_reflectance
        elif plot_type == 'Crystal Transmittance (P polarisation)':
            return self.p_transmittance
        elif plot_type == 'Crystal Transmittance (S polarisation)':
            return self.s_transmittance
        elif plot_type == 'Crystal Absorbtance (P polarisation)':
            return self.p_absorbtance
        elif plot_type == 'Crystal Absorbtance (S polarisation)':
            return self.s_absorbtance
        else:
            print('Error in returning result from CrystalScenarioTab: ',plot_type)
            return None


    def get_results(self, vs_cm1):
        """Return the results of the effective medium theory calculation"""
        if len(vs_cm1)>0 and (self.dirty or len(self.vs_cm1) != len(vs_cm1) or self.vs_cm1[0] != vs_cm1[0] or self.vs_cm1[1] != vs_cm1[1]) :
            self.calculate(vs_cm1)
        else:
            self.notebook.progressbars_update(increment=len(vs_cm1))
        return

    def greyed_out(self):
        """Have a look through the settings and see if we need to grey anything out"""
        # If the single crystal mode is Thick Slab, there is no need for film thickness or substrate permittivity
        if self.settings['Mode'] == 'Thick slab':
            self.film_thickness_sb.setEnabled(False)
            self.substrate_dielectric_sb.setEnabled(False)
        else:
            self.film_thickness_sb.setEnabled(True)
            self.substrate_dielectric_sb.setEnabled(True)

