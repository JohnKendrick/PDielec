import sys
import os.path
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox
from PyQt5.QtCore     import  Qt
from Python.Constants import  support_matrix_db, wavenumber, amu, PI
from Python.Constants import  avogadro_si
from Python.Utilities import  Debug

class ScenarioTab(QWidget):
    def __init__(self, parent, debug=False):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'ScenarioTab:')
        self.settings = {}
        matrix = 'ptfe'
        self.settings['matrix'] = matrix
        self.settings['matrix_density'] = support_matrix_db[matrix][0]
        self.settings['matrix_permittivity'] = support_matrix_db[matrix][1]
        self.settings['mass_fraction'] = 0.1
        self.settings['volume_fraction'] = 0.1
        self.settings['particle_size'] = 0.0001
        self.settings['particle_sigma'] = 0.0
        self.settings['aoverb'] = 1.0
        self.settings['h'] = 0
        self.settings['k'] = 0
        self.settings['l'] = 1
        # get the reader from the main tab
        self.notebook = parent
        self.reader = self.notebook.mainTab.reader
        self.settings['method'] = 'Maxwell-Garnett'
        self.methods = ['Maxwell-Garnett', 'Bruggeman', 'Averaged Permittivity', 'Mie']
        self.settings['shape'] = 'Sphere'
        self.shapes = ['Sphere', 'Needle', 'Plate', 'Ellipsoid']
        self.scenarioIndex = None
        # Create a scenario tab 
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # Support matrix
        #
        self.matrix_cb = QComboBox(self)
        self.matrix_cb.setToolTip('Define the permittivity and density of the support matrix')
        for matrix in support_matrix_db:
            debugger.print('matrix:',matrix, support_matrix_db[matrix])
            self.matrix_cb.addItem(matrix)
        index = self.matrix_cb.findText(self.settings['matrix'], Qt.MatchFixedString)
        if index >=0:
            self.matrix_cb.setCurrentIndex(index)
        else:
            print('support matrix index was not 0',matrix)
        self.matrix_cb.currentIndexChanged.connect(self.on_matrix_cb_changed)
        label = QLabel('Support matrix',self)
        label.setToolTip('Define the permittivity and density of the support matrix')
        form.addRow(label, self.matrix_cb)
        #
        # Support matrix permittivity
        #
        self.density_le = QLineEdit(self) 
        self.density_le.setToolTip('Define the support matrix density. \nThis makes changes to the support density and permittivity')
        self.density_le.setText('{0:.6f}'.format(self.settings['matrix_density']))
        self.density_le.textChanged.connect(self.on_density_le_changed)
        label = QLabel('Support density', self)
        label.setToolTip('Define the support matrix density. \nThis makes changes to the support density and permittivity')
        form.addRow(label, self.density_le)
        #
        # Support matrix permittivity
        #
        self.permittivity_le = QLineEdit(self) 
        self.permittivity_le.setToolTip('Define the support matrix permittivity')
        self.permittivity_le.setText('{0:.6f}'.format(self.settings['matrix_permittivity']))
        self.permittivity_le.textChanged.connect(self.on_permittivity_le_changed)
        label = QLabel('Support permittivity', self)
        label.setToolTip('Define the support matrix permittivity')
        form.addRow(label, self.permittivity_le)
        #
        # Mass fraction of dielectric medium
        #
        self.mf_le = QLineEdit(self) 
        self.mf_le.setToolTip('The mole fraction of the dielectric medium. \nNote that volume and mole fraction are linked')
        self.mf_le.setText('{0:.6f}'.format(self.settings['mass_fraction']))
        self.mf_le.textChanged.connect(self.on_mf_le_changed)
        label = QLabel('Mole fraction of dielectric', self)
        label.setToolTip('The mole fraction of the dielectric medium. \nNote that volume and mole fraction are linked')
        form.addRow(label, self.mf_le)
        #
        # Volume fraction of dielectric medium
        #
        self.vf_le = QLineEdit(self) 
        self.vf_le.setToolTip('The volume fraction of the dielectric medium. \nNote that volume and mole fraction are linked')
        self.vf_le.textChanged.connect(self.on_vf_le_changed)
        self.vf_le.setText('{0:.6f}'.format(self.settings['volume_fraction']))
        label = QLabel('Volume fraction of dielectric', self)
        label.setToolTip('The volume fraction of the dielectric medium. \nNote that volume and mole fraction are linked')
        form.addRow(label, self.vf_le)
        #
        # Calculation method
        #
        self.methods_cb = QComboBox(self)
        self.methods_cb.setToolTip('Choose the calculation method for the effective medium theory')
        for method in self.methods:
            self.methods_cb.addItem(method)
        index = self.methods_cb.findText(self.settings['method'], Qt.MatchFixedString)
        if index >=0:
            self.methods_cb.setCurrentIndex(index)
        else:
            print('Method index was not 0',self.settings['method'])
        self.methods_cb.currentIndexChanged.connect(self.on_methods_cb_changed)
        label = QLabel('Method',self)
        label.setToolTip('Choose the calculation method for the effective medium theory')
        form.addRow(label, self.methods_cb)
        #
        # Particle size option
        #
        self.size_le = QLineEdit(self) 
        self.size_le.setToolTip('Define the particle radius of the sphere in microns. \nOnly applicable for the Mie method')
        self.size_le.setText('{0:.6f}'.format(self.settings['particle_size']))
        self.size_le.textChanged.connect(self.on_size_le_changed)
        label = QLabel('Particle size (radius in microns)',self)
        label.setToolTip('Define the particle radius of the sphere in microns. \nOnly applicable for the Mie method')
        form.addRow(label, self.size_le)
        #
        # Particle sigma option
        #
        self.sigma_le = QLineEdit(self) 
        self.sigma_le.setToolTip('Define the particle size distribution as a lognormal distribution with the give sigma. \nOnly applicable for the Mie method')
        self.sigma_le.setText('{0:.6f}'.format(self.settings['particle_sigma']))
        self.sigma_le.textChanged.connect(self.on_sigma_le_changed)
        label = QLabel('Particle sigma (radius in microns)',self)
        label.setToolTip('Define the particle size distribition as a lognormal with the given sigma. \nOnly applicable for the Mie method')
        form.addRow(label, self.sigma_le)
        #
        # Crystallite shape
        #
        self.shape_cb = QComboBox(self)
        self.shape_cb.setToolTip('Choose a particle shape. \nFor the Mie method only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [hkl].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction')
        for shape in self.shapes:
            self.shape_cb.addItem(shape)
        index = self.shape_cb.findText(self.settings['shape'], Qt.MatchFixedString)
        if index >=0:
            self.shape_cb.setCurrentIndex(index)
        else:
            print('Method index was not 0',self.settings['shape'])
        self.shape_cb.currentIndexChanged.connect(self.on_shape_cb_changed)
        label = QLabel('Shape',self)
        label.setToolTip('Choose a particle shape. \nFor the Mie method only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [hkl].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction')
        form.addRow(label, self.shape_cb)
        #
        # Particle shape information
        # unique direction (hkl) or [hkl]
        self.h_sb = QSpinBox(self)
        self.h_sb.setToolTip('Define the h dimension of the unique direction')
        self.h_sb.setRange(0,20)
        self.h_sb.setValue(self.settings['h'])
        self.h_sb.valueChanged.connect(self.on_h_sb_changed)
        self.k_sb = QSpinBox(self)
        self.k_sb.setToolTip('Define the k dimension of the unique direction')
        self.k_sb.setRange(0,20)
        self.k_sb.setValue(self.settings['k'])
        self.k_sb.valueChanged.connect(self.on_k_sb_changed)
        self.l_sb = QSpinBox(self)
        self.l_sb.setToolTip('Define the l dimension of the unique direction')
        self.l_sb.setRange(0,20)
        self.l_sb.setValue(self.settings['l'])
        self.l_sb.valueChanged.connect(self.on_l_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.h_sb)
        hbox.addWidget(self.k_sb)
        hbox.addWidget(self.l_sb)
        self.hkl_label = QLabel('Unique direction [hkl]',self)
        self.hkl_label.setToolTip('Define the unique direction by [hkl] or (hkl). \n[hkl] is used by needles and ellipsoids.  It defines the unique direction in crystallographic units. \n(hkl) is used by plates it defines a surface and the unique direction is perpendicular to it.')

        form.addRow(self.hkl_label, hbox)
        #
        # a over b ratio for ellipse
        #
        self.aoverb_le = QLineEdit(self) 
        self.aoverb_le.setToolTip('Define the ellipsoid a/b ratio or eccentricity.  \nOnly applicable for the ellipsoid shapes \na/b < 1: oblate ellipsoid \na/b > 1: prolate ellipsoid')
        self.aoverb_le.setText('{0:.6f}'.format(self.settings['aoverb']))
        self.aoverb_le.textChanged.connect(self.on_aoverb_le_changed)
        label = QLabel('Ellipsoid a/b eccentricty',self)
        label.setToolTip('Define the ellipsoid a/b ratio or eccentricity.  \nOnly applicable for the ellipsoid shapes \na/b < 1: oblate ellipsoid \na/b > 1: prolate ellipsoid')
        form.addRow(label, self.aoverb_le)
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
        # Final buttons
        #
        hbox = QHBoxLayout()
        self.pushButton1 = QPushButton('Add another scenario')
        self.pushButton1.setToolTip('Use another scenario to calculate the effect of changing the material on the absorption and permittivity')
        self.pushButton1.clicked.connect(self.pushButton1Clicked)
        hbox.addWidget(self.pushButton1)
        self.pushButton3 = QPushButton('Delete this scenario')
        self.pushButton3.setToolTip('Delete the current scenario')
        self.pushButton3.clicked.connect(self.pushButton3Clicked)
        hbox.addWidget(self.pushButton3)
        form.addRow(hbox)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        # sort out greying of boxes
        self.change_greyed_out()

    def pushButton1Clicked(self):
        # Add another scenario
        debugger.print('Button 1 pressed')
        self.notebook.addScenario(self.scenarioIndex)

    def pushButton3Clicked(self):
        # Delete a scenario
        debugger.print('Button 3 pressed')
        self.notebook.deleteScenario(self.scenarioIndex)

    def crystal_density(self):
        if not self.reader:
            return 1.0
        volume = self.reader.volume
        mass = 0.0
        for m in self.reader.masses:
            mass += m
        return m / (avogadro_si * volume * 1.0e-24)
        

    def on_h_sb_changed(self,value):
        debugger.print('on_h_sb_changed', value)
        self.settings['h'] = value

    def on_k_sb_changed(self,value):
        debugger.print('on_k_sb_changed', value)
        self.settings['k'] = value

    def on_l_sb_changed(self,value):
        debugger.print('on_l_sb_changed', value)
        self.settings['l'] = value

    def on_shape_cb_changed(self,index):
        debugger.print('on shape cb changed', index)
        self.settings['shape'] = self.shapes[index]
        self.change_greyed_out()

    def on_methods_cb_changed(self,index):
        debugger.print('on methods cb changed', index)
        self.settings['method'] = self.methods[index]
        self.change_greyed_out()

    def on_mf_le_changed(self,text):
        debugger.print('on mass fraction line edit changed', text)
        mf = float(text)
        if mf > 0.0:
            self.settings['mass_fraction'] = mf
            self.update_vf_le()

    def update_vf_le(self):
        mf1 = self.settings['mass_fraction']
        mf2 = 1.0 - mf1
        rho1 = self.crystal_density()
        rho2 = self.settings['matrix_density']
        vf1 = 1.0 / ( 1.0 + mf2/mf1 * (rho1/rho2) )
        self.settings['volume_fraction'] = vf1
        self.vf_le.blockSignals(True)
        self.vf_le.setText('{0:.6}'.format(vf1))
        self.vf_le.blockSignals(False)
        debugger.print('Update_vf_le')
        debugger.print('rho 1', rho1)
        debugger.print('rho 2', rho2)
        debugger.print('vf 1 ', vf1)
        
    def on_aoverb_le_changed(self,text):
        debugger.print('on_aoverb_le_changed',text)
        self.settings['aoverb'] = float(text)

    def on_legend_le_changed(self,text):
        debugger.print('on legend change', text)
        self.settings['legend'] = text
        self.legend_le.setText(text)

    def on_sigma_le_changed(self,text):
        debugger.print('on sigma line edit changed', text)
        self.settings['particle_sigma'] = float(text)

    def on_size_le_changed(self,text):
        debugger.print('on size line edit changed', text)
        self.settings['particle_size'] = float(text)

    def on_vf_le_changed(self,text):
        debugger.print('on volume fraction line edit changed', text)
        vf = float(text)
        if vf > 1.0e-12:
            self.settings['volume_fraction'] = vf
            self.update_mf_le()

    def update_mf_le(self):
        vf1 = self.settings['volume_fraction']
        vf2 = 1.0 - vf1
        rho1 = self.crystal_density()
        rho2 = self.settings['matrix_density']
        mf1 = 1.0 / ( 1.0 + (vf2/vf1) * (rho2/rho1) )
        self.settings['mass_fraction'] = mf1
        self.mf_le.blockSignals(True)
        self.mf_le.setText('{0:.6}'.format(mf1))
        self.vf_le.blockSignals(False)
        debugger.print('Update_mf_le')
        debugger.print('rho 1', rho1)
        debugger.print('rho 2', rho2)
        debugger.print('mf 1 ', mf1)

    def on_matrix_cb_changed(self,index):
        debugger.print('on matrix combobox changed', index)
        debugger.print('on matrix combobox changed', self.matrix_cb.currentText())
        matrix = self.matrix_cb.currentText()
        self.matrix_cb.blockSignals(True)
        self.density_le.blockSignals(True)
        self.permittivity_le.blockSignals(True)
        self.settings['matrix'] = matrix
        self.settings['matrix_density'] = support_matrix_db[matrix][0]
        self.settings['matrix_permittivity'] = support_matrix_db[matrix][1]
        self.density_le.setText('{0:.6f}'.format(self.settings['matrix_density']))
        self.permittivity_le.setText('{0:.6f}'.format(self.settings['matrix_permittivity']))
        # volume fraction takes precedence
        self.update_mf_le()
        self.update_vf_le()
        self.matrix_cb.blockSignals(False)
        self.density_le.blockSignals(False)
        self.permittivity_le.blockSignals(False)

    def on_density_le_changed(self,text):
        self.settings['density'] = float(text)
        # mass fraction taked precedence
        self.update_mf_le()
        self.update_vf_le()
        debugger.print('on density line edit changed', text)

    def on_permittivity_le_changed(self,text):
        self.settings['permittivity'] = float(text)
        debugger.print('on density line edit changed', text)

    def set_reader(self,reader):
        self.reader = reader

    def change_greyed_out(self):
        # Have a look through the settings and see if we need to grey anything out
        if self.settings['method'] == 'Mie':
            self.size_le.setEnabled(True)
            self.sigma_le.setEnabled(True)
        else:
            self.size_le.setEnabled(False)
            self.sigma_le.setEnabled(False)
        # deal with shapes
        if self.settings['shape'] == 'Ellipsoid':
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText('Unique direction [hkl]')
            self.aoverb_le.setEnabled(True)
        elif self.settings['shape'] == 'Plate':
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText('Unique direction (hkl)')
            self.aoverb_le.setEnabled(False)
        elif self.settings['shape'] == 'Needle':
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText('Unique direction [hkl]')
            self.aoverb_le.setEnabled(False)
        elif self.settings['shape'] == 'Sphere':
            self.h_sb.setEnabled(False)
            self.k_sb.setEnabled(False)
            self.l_sb.setEnabled(False)
            self.aoverb_le.setEnabled(False)
        else:
            print('ScenarioTab: Shape not recognised', self.settings['shape'])
        
    def setScenarioIndex(self,index):
        self.scenarioIndex = index
        text = self.legend_le.text()
        if text == 'Scenario legend':
            self.legend_le.setText('Scenario '+str(index + 1))
        return


    def refresh(self):
        debugger.print('refresh')
        # First see if we can get the reader from the mainTab
        self.reader = self.notebook.mainTab.reader
        # block signals to all the widgets
        self.matrix_cb.blockSignals(True)
        self.density_le.blockSignals(True)
        self.permittivity_le.blockSignals(True)
        #self.mf_le.blockSignals(True)
        #self.vf_le.blockSignals(True)
        self.methods_cb.blockSignals(True)
        self.size_le.blockSignals(True)
        self.shape_cb.blockSignals(True)
        self.h_sb.blockSignals(True)
        self.k_sb.blockSignals(True)
        self.l_sb.blockSignals(True)
        self.aoverb_le.blockSignals(True)
        # use the settings values to initialise the widgets
        index = self.matrix_cb.findText(self.settings['matrix'], Qt.MatchFixedString)
        self.matrix_cb.setCurrentIndex(index)
        self.density_le.setText('{0:.6f}'.format(self.settings['matrix_density']))
        self.permittivity_le.setText('{0:.6f}'.format(self.settings['matrix_permittivity']))
        # volume fraction takes precedence
        self.update_mf_le()
        self.update_vf_le()
        #self.vf_le.setText('{0:.6f}'.format(self.settings['volume_fraction']))
        #self.mf_le.setText('{0:.6f}'.format(self.settings['mass_fraction']))
        index = self.methods_cb.findText(self.settings['method'], Qt.MatchFixedString)
        self.methods_cb.setCurrentIndex(index)
        self.size_le.setText('{0:.6f}'.format(self.settings['particle_size']))
        index = self.shape_cb.findText(self.settings['shape'], Qt.MatchFixedString)
        self.shape_cb.setCurrentIndex(index)
        self.h_sb.setValue(self.settings['h'])
        self.k_sb.setValue(self.settings['k'])
        self.l_sb.setValue(self.settings['l'])
        self.aoverb_le.setText('{0:.6f}'.format(self.settings['aoverb']))
        self.change_greyed_out()
        # ublock signals to all the widgets
        self.matrix_cb.blockSignals(False)
        self.density_le.blockSignals(False)
        self.permittivity_le.blockSignals(False)
        self.mf_le.blockSignals(False)
        self.vf_le.blockSignals(False)
        self.methods_cb.blockSignals(False)
        self.size_le.blockSignals(False)
        self.shape_cb.blockSignals(False)
        self.h_sb.blockSignals(False)
        self.k_sb.blockSignals(False)
        self.l_sb.blockSignals(False)
        self.aoverb_le.blockSignals(False)
        return
