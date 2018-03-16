# -*- coding: utf8 -*-
import sys
import os.path
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit, QDoubleSpinBox
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
        self.dirty = True
        self.settings = {}
        self.notebook = parent
        self.notebook.newPlottingCalculationRequired = True
        matrix = 'ptfe'
        self.settings['Matrix'] = matrix
        self.settings['Matrix density'] = support_matrix_db[matrix][0]
        self.settings['Matrix permittivity'] = support_matrix_db[matrix][1]
        self.settings['Mass fraction'] = 0.1
        self.settings['Volume fraction'] = 0.1
        self.settings['Particle size(mu)'] = 0.0001
        self.settings['Particle size distribution sigma(mu)'] = 0.0
        self.settings['Ellipsoid a/b'] = 1.0
        self.settings['Unique direction - h'] = 0
        self.settings['Unique direction - k'] = 0
        self.settings['Unique direction - l'] = 1
        self.settings['Mass or volume fraction'] = 'volume'
        # get the reader from the main tab
        self.notebook = parent
        self.reader = self.notebook.mainTab.reader
        self.settings['Effective medium method'] = 'Maxwell-Garnett'
        self.methods = ['Maxwell-Garnett', 'Bruggeman', 'Averaged Permittivity', 'Mie']
        self.settings['Particle shape'] = 'Sphere'
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
        self.matrix_cb.addItems(support_matrix_db)
        index = self.matrix_cb.findText(self.settings['Matrix'], Qt.MatchFixedString)
        if index >=0:
            self.matrix_cb.setCurrentIndex(index)
        else:
            print('support matrix index was not 0',matrix)
        self.matrix_cb.activated.connect(self.on_matrix_cb_activated)
        label = QLabel('Support matrix',self)
        label.setToolTip('Define the permittivity and density of the support matrix')
        form.addRow(label, self.matrix_cb)
        #
        # Support matrix permittivity
        #
        self.density_sb = QDoubleSpinBox(self) 
        self.density_sb.setRange(0.001, 100.0)
        self.density_sb.setSingleStep(0.01)
        self.density_sb.setDecimals(3)
        self.density_sb.setToolTip('Define the support matrix density. \nThis makes changes to the support density and permittivity')
        self.density_sb.setValue(self.settings['Matrix density'])
        self.density_sb.valueChanged.connect(self.on_density_sb_changed)
        label = QLabel('Support density', self)
        label.setToolTip('Define the support matrix density. \nThis makes changes to the support density and permittivity')
        form.addRow(label, self.density_sb)
        #
        # Support matrix permittivity
        #
        self.permittivity_sb = QDoubleSpinBox(self) 
        self.permittivity_sb.setRange(0.001, 100.0)
        self.permittivity_sb.setSingleStep(0.01)
        self.permittivity_sb.setDecimals(3)
        self.permittivity_sb.setToolTip('Define the support matrix permittivity')
        self.permittivity_sb.setValue(self.settings['Matrix permittivity'])
        self.permittivity_sb.valueChanged.connect(self.on_permittivity_sb_changed)
        label = QLabel('Support permittivity', self)
        label.setToolTip('Define the support matrix permittivity')
        form.addRow(label, self.permittivity_sb)
        #
        # Mass fraction of dielectric medium
        #
        self.mf_sb = QDoubleSpinBox(self)
        self.mf_sb.setRange(0.000001, 100.0)
        self.mf_sb.setSingleStep(0.1)
        self.mf_sb.setDecimals(6)
        self.mf_sb.setToolTip('The percentage mass fraction of the dielectric medium. \nNote that volume and mass fraction are linked')
        self.mf_sb.setValue(100.0*self.settings['Mass fraction'])
        self.mf_sb.valueChanged.connect(self.on_mf_sb_changed)
        label = QLabel('% mass fraction of dielectric', self)
        label.setToolTip('The percentage mass fraction of the dielectric medium. \nNote that volume and mass fraction are linked')
        form.addRow(label, self.mf_sb)
        #
        # Volume fraction of dielectric medium
        #
        self.vf_sb = QDoubleSpinBox(self)
        self.vf_sb.setRange(0.000001, 100.0)
        self.vf_sb.setSingleStep(0.1)
        self.vf_sb.setDecimals(6)
        self.vf_sb.setToolTip('The percentage volume fraction of the dielectric medium. \nNote that volume and mass fraction are linked')
        self.vf_sb.valueChanged.connect(self.on_vf_sb_changed)
        self.vf_sb.setValue(100.0*self.settings['Volume fraction'])
        label = QLabel('% volume fraction of dielectric', self)
        label.setToolTip('The percentage volume fraction of the dielectric medium. \nNote that volume and mass fraction are linked')
        form.addRow(label, self.vf_sb)
        #
        # Calculation method
        #
        self.methods_cb = QComboBox(self)
        self.methods_cb.setToolTip('Choose the calculation method for the effective medium theory')
        self.methods_cb.addItems(self.methods)
        index = self.methods_cb.findText(self.settings['Effective medium method'], Qt.MatchFixedString)
        if index >=0:
            self.methods_cb.setCurrentIndex(index)
        else:
            print('Method index was not 0',self.settings['Effective medium method'])
        self.methods_cb.activated.connect(self.on_methods_cb_activated)
        label = QLabel('Method',self)
        label.setToolTip('Choose the calculation method for the effective medium theory')
        form.addRow(label, self.methods_cb)
        #
        # Particle size option
        #
        self.size_sb = QDoubleSpinBox(self)
        self.size_sb.setRange(0.000001, 1000.0)
        self.size_sb.setSingleStep(0.1)
        self.size_sb.setDecimals(6)
        self.size_sb.setToolTip('Define the particle radius of the sphere in μm. \nOnly applicable for the Mie method')
        self.size_sb.setValue(self.settings['Particle size(mu)'])
        self.size_sb.valueChanged.connect(self.on_size_sb_changed)
        label = QLabel('Particle radius (μm)',self)
        label.setToolTip('Define the particle radius of the sphere in μm. \nOnly applicable for the Mie method')
        form.addRow(label, self.size_sb)
        #
        # Particle sigma option
        #
        self.sigma_sb = QDoubleSpinBox(self)
        self.sigma_sb.setRange(0.0, 1000.0)
        self.sigma_sb.setSingleStep(0.1)
        self.sigma_sb.setDecimals(6)
        self.sigma_sb.setToolTip('Define the particle size distribution as a lognormal distribution with the given sigma. \nOnly applicable for the Mie method')
        self.sigma_sb.setValue(self.settings['Particle size distribution sigma(mu)'])
        self.sigma_sb.valueChanged.connect(self.on_sigma_sb_changed)
        label = QLabel('Particle sigma (μm)',self)
        label.setToolTip('Define the particle size distribition as a lognormal with the given sigma. \nOnly applicable for the Mie method')
        form.addRow(label, self.sigma_sb)
        #
        # Crystallite shape
        #
        self.shape_cb = QComboBox(self)
        self.shape_cb.setToolTip('Choose a particle shape. \nFor the Mie method only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [abc].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction')
        self.shape_cb.addItems(self.shapes)
        index = self.shape_cb.findText(self.settings['Particle shape'], Qt.MatchFixedString)
        if index >=0:
            self.shape_cb.setCurrentIndex(index)
        else:
            print('Method index was not 0',self.settings['Particle shape'])
        self.shape_cb.activated.connect(self.on_shape_cb_activated)
        label = QLabel('Particle shape',self)
        label.setToolTip('Choose a particle shape. \nFor the Mie method only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [abc].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction')
        form.addRow(label, self.shape_cb)
        #
        # Particle shape information
        # unique direction (hkl) or [abc]
        self.h_sb = QSpinBox(self)
        self.h_sb.setToolTip('Define the h dimension of the unique direction')
        self.h_sb.setRange(-20,20)
        self.h_sb.setValue(self.settings['Unique direction - h'])
        self.h_sb.valueChanged.connect(self.on_h_sb_changed)
        self.k_sb = QSpinBox(self)
        self.k_sb.setToolTip('Define the k dimension of the unique direction')
        self.k_sb.setRange(-20,20)
        self.k_sb.setValue(self.settings['Unique direction - k'])
        self.k_sb.valueChanged.connect(self.on_k_sb_changed)
        self.l_sb = QSpinBox(self)
        self.l_sb.setToolTip('Define the l dimension of the unique direction')
        self.l_sb.setRange(-20,20)
        self.l_sb.setValue(self.settings['Unique direction - l'])
        self.l_sb.valueChanged.connect(self.on_l_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.h_sb)
        hbox.addWidget(self.k_sb)
        hbox.addWidget(self.l_sb)
        self.hkl_label = QLabel('Unique direction [abc]',self)
        self.hkl_label.setToolTip('Define the unique direction by [abc] or (hkl). \n[abc] is used by needles and ellipsoids.  It defines the unique direction in crystallographic units. \n(hkl) is used by plates it defines a surface and the unique direction is perpendicular to it.')

        form.addRow(self.hkl_label, hbox)
        #
        # a over b ratio for ellipse
        #
        self.aoverb_sb = QDoubleSpinBox(self)
        self.aoverb_sb.setRange(0.0, 1000.0)
        self.aoverb_sb.setSingleStep(0.1)
        self.aoverb_sb.setDecimals(6)
        self.aoverb_sb.setToolTip('Define the ellipsoid a/b ratio or eccentricity.  \nOnly applicable for the ellipsoid shapes \na/b < 1: oblate ellipsoid \na/b > 1: prolate ellipsoid')
        self.aoverb_sb.setValue(self.settings['Ellipsoid a/b'])
        self.aoverb_sb.valueChanged.connect(self.on_aoverb_sb_changed)
        label = QLabel('Ellipsoid a/b eccentricty',self)
        label.setToolTip('Define the ellipsoid a/b ratio or eccentricity.  \nOnly applicable for the ellipsoid shapes \na/b < 1: oblate ellipsoid \na/b > 1: prolate ellipsoid')
        form.addRow(label, self.aoverb_sb)
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
        self.notebook.addScenario(copyFromIndex=self.scenarioIndex)

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
        density = mass / (avogadro_si * volume * 1.0e-24)
        return density
        

    def on_h_sb_changed(self,value):
        debugger.print('on_h_sb_changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Unique direction - h'] = value

    def on_k_sb_changed(self,value):
        debugger.print('on_k_sb_changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Unique direction - k'] = value

    def on_l_sb_changed(self,value):
        debugger.print('on_l_sb_changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Unique direction - l'] = value

    def on_shape_cb_activated(self,index):
        debugger.print('on shape cb activated', index)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Particle shape'] = self.shapes[index]
        self.change_greyed_out()

    def on_methods_cb_activated(self,index):
        debugger.print('on methods cb activated', index)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Effective medium method'] = self.methods[index]
        self.change_greyed_out()

    def on_mf_sb_changed(self,value):
        debugger.print('on mass fraction line edit changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Mass fraction'] =  value/100.0
        if self.settings['Mass or volume fraction'] == 'mass':
            self.update_vf_sb()
            self.update_mf_sb()
        else:
            self.update_mf_sb()
            self.update_vf_sb()

    def update_vf_sb(self):
        mf1 = self.settings['Mass fraction']
        mf2 = 1.0 - mf1
        rho1 = self.crystal_density()
        rho2 = self.settings['Matrix density']
        vf1 = 1.0 / ( 1.0 + mf2/mf1 * (rho1/rho2) )
        self.settings['Volume fraction'] = vf1
        self.vf_sb.blockSignals(True)
        self.vf_sb.setValue(100.0*vf1)
        self.vf_sb.blockSignals(False)
        debugger.print('Update_vf_sb')
        debugger.print('rho 1', rho1)
        debugger.print('rho 2', rho2)
        debugger.print('vf 1 ', vf1)
        
    def on_aoverb_sb_changed(self,value):
        debugger.print('on_aoverb_le_changed',value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Ellipsoid a/b'] = value

    def on_legend_le_changed(self,text):
        debugger.print('on legend change', text)
        self.dirty = True
        self.settings['Legend'] = text
        self.legend_le.setText(text)

    def on_sigma_sb_changed(self,value):
        debugger.print('on sigma line edit changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Particle size distribution sigma(mu)'] = value

    def on_size_sb_changed(self,value):
        debugger.print('on size line edit changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Particle size(mu)'] = value

    def on_vf_sb_changed(self,value):
        debugger.print('on volume fraction line edit changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.settings['Volume fraction'] = value/100.0
        if self.settings['Mass or volume fraction'] == 'volume':
            self.update_mf_sb()
            self.update_vf_sb()
        else:
            self.update_vf_sb()
            self.update_mf_sb()

    def update_mf_sb(self):
        vf1 = self.settings['Volume fraction']
        vf2 = 1.0 - vf1
        rho1 = self.crystal_density()
        rho2 = self.settings['Matrix density']
        mf1 = 1.0 / ( 1.0 + (vf2/vf1) * (rho2/rho1) )
        self.settings['Mass fraction'] = mf1
        self.mf_sb.blockSignals(True)
        self.mf_sb.setValue(100.0*mf1)
        self.mf_sb.blockSignals(False)
        debugger.print('Update_mf_sb')
        debugger.print('rho 1', rho1)
        debugger.print('rho 2', rho2)
        debugger.print('mf 1 ', mf1)

    def on_matrix_cb_activated(self,index):
        debugger.print('on matrix combobox activated', index)
        debugger.print('on matrix combobox activated', self.matrix_cb.currentText())
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        matrix = self.matrix_cb.currentText()
        self.matrix_cb.blockSignals(True)
        self.density_sb.blockSignals(True)
        self.permittivity_sb.blockSignals(True)
        self.settings['Matrix'] = matrix
        self.settings['Matrix density'] = support_matrix_db[matrix][0]
        self.settings['Matrix permittivity'] = support_matrix_db[matrix][1]
        self.density_sb.setValue(self.settings['Matrix density'])
        self.permittivity_sb.setValue(self.settings['Matrix permittivity'])
        # volume fraction takes precedence
        self.update_mf_sb()
        self.update_vf_sb()
        self.matrix_cb.blockSignals(False)
        self.density_sb.blockSignals(False)
        self.permittivity_sb.blockSignals(False)

    def on_density_sb_changed(self,value):
        self.settings['Matrix density'] = value
        # volume fraction taked precedence
        self.update_mf_sb()
        self.update_vf_sb()
        debugger.print('on density line edit changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True

    def on_permittivity_sb_changed(self,value):
        self.settings['Matrix permittivity'] = value
        debugger.print('on permittivity line edit changed', value)
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True

    def set_reader(self,reader):
        self.dirty = True
        self.notebook.newPlottingCalculationRequired = True
        self.reader = reader

    def change_greyed_out(self):
        # Have a look through the settings and see if we need to grey anything out
        if self.settings['Effective medium method'] == 'Mie':
            self.size_sb.setEnabled(True)
            self.sigma_sb.setEnabled(True)
            for i,shape in enumerate(self.shapes):
                self.shape_cb.model().item(i).setEnabled(False)
            self.settings['Particle shape'] = 'Sphere'
            self.shape_cb.setEnabled(True)
            index = self.shape_cb.findText(self.settings['Particle shape'], Qt.MatchFixedString)
            if index >=0:
                self.shape_cb.model().item(index).setEnabled(True)
                self.shape_cb.setCurrentIndex(index)
            else:
                print('Method index was not 0',self.settings['Particle shape'])
        elif self.settings['Effective medium method'] == 'Averaged Permittivity':
            self.size_sb.setEnabled(False)
            self.sigma_sb.setEnabled(False)
            self.settings['Particle shape'] = 'Sphere'
            index = self.shape_cb.findText(self.settings['Particle shape'], Qt.MatchFixedString)
            if index >=0:
                self.shape_cb.model().item(index).setEnabled(True)
                self.shape_cb.setCurrentIndex(index)
            self.shape_cb.setEnabled(False)
            for i,shape in enumerate(self.shapes):
                self.shape_cb.model().item(i).setEnabled(False)
        else:
            self.size_sb.setEnabled(False)
            self.sigma_sb.setEnabled(False)
            self.shape_cb.setEnabled(True)
            for i,shape in enumerate(self.shapes):
                self.shape_cb.model().item(i).setEnabled(True)
        # deal with shapes
        if self.settings['Particle shape'] == 'Ellipsoid':
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText('Unique direction [abc]')
            self.aoverb_sb.setEnabled(True)
        elif self.settings['Particle shape'] == 'Plate':
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText('Unique direction (hkl)')
            self.aoverb_sb.setEnabled(False)
        elif self.settings['Particle shape'] == 'Needle':
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText('Unique direction [abc]')
            self.aoverb_sb.setEnabled(False)
        elif self.settings['Particle shape'] == 'Sphere':
            self.h_sb.setEnabled(False)
            self.k_sb.setEnabled(False)
            self.l_sb.setEnabled(False)
            self.aoverb_sb.setEnabled(False)
        else:
            print('ScenarioTab: Shape not recognised', self.settings['Particle shape'])
        
    def setScenarioIndex(self,index):
        self.scenarioIndex = index
        text = self.legend_le.text()
        if text == 'Scenario legend':
            self.legend_le.setText('Scenario '+str(index + 1))
        return


    def print_settings(self):
        debugger.print('SETTINGS')
        for key in self.settings:
            debugger.print(key, self.settings[key]) 
        
    def refresh(self,force=False):
        if not self.dirty and not force:
            debugger.print('refresh aborted', self.dirty,force)
            return
        debugger.print('refresh', force)
        # Tell the main notebook that we need to recalculate any plot
        self.notebook.newPlottingCalculationRequired = True
        # First see if we can get the reader from the mainTab
        self.reader = self.notebook.mainTab.reader
        # block signals to all the widgets
        self.matrix_cb.blockSignals(True)
        self.density_sb.blockSignals(True)
        self.permittivity_sb.blockSignals(True)
        self.mf_sb.blockSignals(True)
        self.vf_sb.blockSignals(True)
        self.methods_cb.blockSignals(True)
        self.size_sb.blockSignals(True)
        self.sigma_sb.blockSignals(True)
        self.shape_cb.blockSignals(True)
        self.h_sb.blockSignals(True)
        self.k_sb.blockSignals(True)
        self.l_sb.blockSignals(True)
        self.aoverb_sb.blockSignals(True)
        # use the settings values to initialise the widgets
        index = self.matrix_cb.findText(self.settings['Matrix'], Qt.MatchFixedString)
        self.matrix_cb.setCurrentIndex(index)
        self.density_sb.setValue(self.settings['Matrix density'])
        self.permittivity_sb.setValue(self.settings['Matrix permittivity'])
        if self.settings['Mass or volume fraction'] == 'volume':
            # volume fraction takes precedence
            self.update_mf_sb()
            self.update_vf_sb()
        else:
            # mass fraction takes precedence
            self.update_vf_sb()
            self.update_mf_sb()
        #
        index = self.methods_cb.findText(self.settings['Effective medium method'], Qt.MatchFixedString)
        self.methods_cb.setCurrentIndex(index)
        self.size_sb.setValue(self.settings['Particle size(mu)'])
        self.sigma_sb.setValue(self.settings['Particle size distribution sigma(mu)'])
        index = self.shape_cb.findText(self.settings['Particle shape'], Qt.MatchFixedString)
        self.shape_cb.setCurrentIndex(index)
        self.h_sb.setValue(self.settings['Unique direction - h'])
        self.k_sb.setValue(self.settings['Unique direction - k'])
        self.l_sb.setValue(self.settings['Unique direction - l'])
        self.aoverb_sb.setValue(self.settings['Ellipsoid a/b'])
        self.change_greyed_out()
        # ublock signals to all the widgets
        self.matrix_cb.blockSignals(False)
        self.density_sb.blockSignals(False)
        self.permittivity_sb.blockSignals(False)
        self.mf_sb.blockSignals(False)
        self.vf_sb.blockSignals(False)
        self.methods_cb.blockSignals(False)
        self.size_sb.blockSignals(False)
        self.sigma_sb.blockSignals(False)
        self.shape_cb.blockSignals(False)
        self.h_sb.blockSignals(False)
        self.k_sb.blockSignals(False)
        self.l_sb.blockSignals(False)
        self.aoverb_sb.blockSignals(False)
        self.dirty = False
        return
