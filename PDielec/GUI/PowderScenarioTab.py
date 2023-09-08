# -*- coding: utf8 -*-
from PyQt5.QtWidgets         import QPushButton, QWidget, QFileDialog
from PyQt5.QtWidgets         import QComboBox, QLabel, QLineEdit, QDoubleSpinBox
from PyQt5.QtWidgets         import QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets         import QSpinBox
from PyQt5.QtCore            import Qt, QCoreApplication

from PDielec                 import __file__  as PDielec_init_filename
from PDielec.Constants       import avogadro_si
from PDielec.Utilities       import Debug
from PDielec.GUI.ScenarioTab import ScenarioTab
from multiprocessing         import Array
from PDielec.Constants       import wavenumber
from functools               import partial
import ctypes
import PDielec.Calculator as Calculator
import PDielec.DielectricFunction as DielectricFunction
from PDielec.Materials       import MaterialsDataBase
from PDielec.Materials       import Material
import numpy as np
import os

class PowderScenarioTab(ScenarioTab):
    def __init__(self, parent, debug=False):
        ScenarioTab.__init__(self,parent)
        # super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'PowderScenarioTab:')
        debugger.print('Start:: initialiser')
        self.scenarioType = 'Powder'
        self.settings['Scenario type'] = 'Powder'
        self.noCalculationsRequired = 1
        matrix = 'ptfe'
        self.settings['Matrix'] = None
        self.settings['Matrix density'] = None
        self.settings['Matrix permittivity'] = None
        self.settings['Bubble radius'] = 30.0
        self.settings['Bubble volume fraction'] = 0.0
        self.settings['Mass fraction'] = 0.1
        self.settings['Volume fraction'] = 0.1
        self.settings['Particle size(mu)'] = 0.0001
        self.settings['Particle size distribution sigma(mu)'] = 0.0
        self.settings['Ellipsoid a/b'] = 1.0
        self.settings['Unique direction - h'] = 0
        self.settings['Unique direction - k'] = 0
        self.settings['Unique direction - l'] = 1
        self.settings['Mass or volume fraction'] = 'volume'
        self.settings['ATR material refractive index'] = 4.0
        self.settings['ATR theta'] = 45.0
        self.settings['ATR S polarisation fraction'] = 0.5
        self.settings['Effective medium method'] = 'Maxwell-Garnett'
        self.settings['Particle shape'] = 'Sphere'
        PDielec_Directory = os.path.dirname(PDielec_init_filename)
        self.settings['Materials database'] = PDielec_Directory + '/MaterialsDataBase.xlsx'
        self.methods = ['Maxwell-Garnett', 'Bruggeman', 'Averaged Permittivity', 'Mie']
        self.shapes = ['Sphere', 'Needle', 'Plate', 'Ellipsoid']
        self.direction = np.array([0,0,0])
        self.depolarisation = np.array([0,0,0])
        self.scenarioIndex = None
        self.refreshRequired = True
        self.calculationRequired = False
        self.reader = None
        self.realPermittivity = []
        self.imagPermittivity = []
        self.absorptionCoefficient = []
        self.molarAbsorptionCoefficient = []
        self.sp_atr = []
        # Open the database
        self.DataBase = MaterialsDataBase(self.settings['Materials database'],debug=debug)
        # Create a scenario tab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # Option to open a database of permittivities for the support
        #
        self.openDB_button = QPushButton("Open materials' database")
        self.openDB_button.clicked.connect(self.openDB_button_clicked)
        self.openDB_button.setToolTip("Open a new materials' database (.xlsx file)")
        label = QLabel("Open materials' database")
        label.setToolTip("Open a new materials' database (.xlsx file)")
        form.addRow(label, self.openDB_button)
        #
        # Add a name for the database
        #
        self.database_le = QLineEdit(self)
        self.database_le.setToolTip("Provides information about the name of the materials' database")
        self.database_le.setText(self.settings['Materials database'])
        self.database_le.setReadOnly(True)
        label = QLabel("Current materials' database")
        label.setToolTip("Provides information about the name of the materials' database")
        form.addRow(label, self.database_le)
        #
        # Support matrix, read information from the database
        #
        hbox = QHBoxLayout()
        self.matrix_cb = QComboBox(self)
        self.matrix_cb.setToolTip('Define the permittivity and density of the support matrix')
        self.materialNames = self.DataBase.getSheetNames()
        self.materialNames.append('Define material manually')
        self.matrix_cb.addItems(self.materialNames)
        if not self.settings['Matrix'] in self.materialNames:
            self.settings['Matrix'] = self.materialNames[0]
        index = self.matrix_cb.findText(self.settings['Matrix'], Qt.MatchFixedString)
        if index >=0:
            self.matrix_cb.setCurrentIndex(index)
        else:
            print('support matrix index was not 0',matrix)
        self.matrix_cb.activated.connect(self.on_matrix_cb_activated)
        if 'Define material manually' in self.settings['Matrix']:
            permittivityObject = DielectricFunction.ConstantScalar(self.settings['Matrix permittivity'])
            self.matrixMaterial = Material('manual',permittivityObject=permittivityObject,density=self.settings['Matrix density'])
        else:
            self.matrixMaterial = self.DataBase.getMaterial(self.settings['Matrix'])
        label = QLabel('Support matrix',self)
        label.setToolTip('Define the permittivity and density of the support matrix')
        self.matrix_info_le = QLineEdit(self)
        self.matrix_info_le.setToolTip("Provides details about database entry")
        text = self.matrixMaterial.getInformation()
        self.matrix_info_le.setText(text)
        self.matrix_info_le.setReadOnly(True)
        hbox.addWidget(self.matrix_cb)
        hbox.addWidget(self.matrix_info_le)
        form.addRow(label, hbox)
        # Set the Matrix density and permittivity at 0cm-1
        self.settings['Matrix density'] = self.matrixMaterial.getDensity()
        self.matrixPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
        self.settings['Matrix permittivity'] = self.matrixPermittivityFunction(0.0)
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
        hbox = QHBoxLayout()
        self.permittivity_r_sb = QDoubleSpinBox(self)
        self.permittivity_r_sb.setRange(0.0, 1000.0)
        self.permittivity_r_sb.setSingleStep(0.01)
        self.permittivity_r_sb.setDecimals(3)
        self.permittivity_r_sb.setToolTip('Define the real component of the support matrix permittivity')
        self.permittivity_r_sb.setValue(np.real(self.settings['Matrix permittivity']))
        self.permittivity_r_sb.valueChanged.connect(self.on_permittivity_r_sb_changed)
        hbox.addWidget(self.permittivity_r_sb)
        self.permittivity_i_sb = QDoubleSpinBox(self)
        self.permittivity_i_sb.setRange(0.0, 1000.0)
        self.permittivity_i_sb.setSingleStep(0.01)
        self.permittivity_i_sb.setDecimals(3)
        self.permittivity_i_sb.setToolTip('Define imaginary component of the the support matrix permittivity')
        self.permittivity_i_sb.setValue(np.imag(self.settings['Matrix permittivity']))
        self.permittivity_i_sb.valueChanged.connect(self.on_permittivity_i_sb_changed)
        hbox.addWidget(self.permittivity_i_sb)
        label = QLabel('Support permittivity', self)
        label.setToolTip('Define the complex support matrix permittivity')
        form.addRow(label, hbox)
        #
        # Bubble volume fraction
        #
        self.bubble_vf_sb = QDoubleSpinBox(self)
        self.bubble_vf_sb.setRange(0.0, 100.0*(1.0-self.settings['Volume fraction']))
        self.bubble_vf_sb.setSingleStep(1.0)
        self.bubble_vf_sb.setDecimals(1)
        self.bubble_vf_sb.setToolTip('Define the % volume fraction of air bubble inclusions in the matrix')
        self.bubble_vf_sb.setValue(100*self.settings['Bubble volume fraction'])
        self.bubble_vf_sb.valueChanged.connect(self.on_bubble_vf_sb_changed)
        label = QLabel('% Air void volume fraction', self)
        label.setToolTip('Define the % volume fraction of air bubble inclusions in the matrix')
        form.addRow(label, self.bubble_vf_sb)
        #
        # Bubble radius in microns
        #
        self.bubble_radius_sb = QDoubleSpinBox(self)
        self.bubble_radius_sb.setRange(0.001, 1000.0)
        self.bubble_radius_sb.setSingleStep(1.0)
        self.bubble_radius_sb.setDecimals(3)
        self.bubble_radius_sb.setToolTip('Define the air bubble radius')
        self.bubble_radius_sb.setValue(self.settings['Bubble radius'])
        self.bubble_radius_sb.valueChanged.connect(self.on_bubble_radius_sb_changed)
        label = QLabel('Air void radius (μm)', self)
        label.setToolTip('Define the air void radius')
        form.addRow(label, self.bubble_radius_sb)
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
        label = QLabel('% Mass fraction of dielectric', self)
        label.setToolTip('The percentage mass fraction of the dielectric medium. \nNote that volume and mass fraction are linked')
        form.addRow(label, self.mf_sb)
        #
        # Volume fraction of dielectric medium
        #
        self.vf_sb = QDoubleSpinBox(self)
        self.vf_sb.setRange(0.000001, 100.0*(1.0-self.settings['Bubble volume fraction']))
        self.vf_sb.setSingleStep(0.1)
        self.vf_sb.setDecimals(6)
        self.vf_sb.setToolTip('The percentage volume fraction of the dielectric medium. \nNote that volume and mass fraction are linked')
        self.vf_sb.valueChanged.connect(self.on_vf_sb_changed)
        self.vf_sb.setValue(100.0*self.settings['Volume fraction'])
        label = QLabel('% Volume fraction of dielectric', self)
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
        self.size_sb.setToolTip('Define the particle radius of the sphere in μm.')
        self.size_sb.setValue(self.settings['Particle size(mu)'])
        self.size_sb.valueChanged.connect(self.on_size_sb_changed)
        label = QLabel('Particle radius (μm)',self)
        label.setToolTip('Define the particle radius of the sphere in μm.')
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
        self.shape_cb.setToolTip('Choose a particle shape. \nFor the Mie methods only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [abc].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction')
        self.shape_cb.addItems(self.shapes)
        index = self.shape_cb.findText(self.settings['Particle shape'], Qt.MatchFixedString)
        if index >=0:
            self.shape_cb.setCurrentIndex(index)
        else:
            print('Method index was not 0',self.settings['Particle shape'])
        self.shape_cb.activated.connect(self.on_shape_cb_activated)
        label = QLabel('Particle shape',self)
        label.setToolTip('Choose a particle shape. \nFor the Mie methods only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [abc].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction')
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
        # Add ATR options
        # Refractive Index
        self.atr_index_sb = QDoubleSpinBox(self)
        self.atr_index_sb.setRange(0.001, 100.0)
        self.atr_index_sb.setSingleStep(0.01)
        self.atr_index_sb.setDecimals(3)
        self.atr_index_sb.setToolTip('Define the ATR material refractive index')
        self.atr_index_sb.setValue(self.settings['ATR material refractive index'])
        self.atr_index_sb.valueChanged.connect(self.on_atr_index_sb_changed)
        label = QLabel('ATR material refractive index', self)
        label.setToolTip('Define the ATR material refractive index')
        form.addRow(label, self.atr_index_sb)
        # Incident angle in degreees
        self.atr_incident_ang_sb = QDoubleSpinBox(self)
        self.atr_incident_ang_sb.setRange(0.0, 180.0)
        self.atr_incident_ang_sb.setSingleStep(0.1)
        self.atr_incident_ang_sb.setDecimals(1)
        self.atr_incident_ang_sb.setToolTip('Define the ATR incident angle')
        self.atr_incident_ang_sb.setValue(self.settings['ATR theta'])
        self.atr_incident_ang_sb.valueChanged.connect(self.on_atr_incident_ang_sb_changed)
        label = QLabel('ATR incident angle', self)
        label.setToolTip('Define the ATR incident angle')
        form.addRow(label, self.atr_incident_ang_sb)
        # S polarisation fraction
        self.atr_spolfrac_sb = QDoubleSpinBox(self)
        self.atr_spolfrac_sb.setRange(0.0, 1.0)
        self.atr_spolfrac_sb.setSingleStep(0.01)
        self.atr_spolfrac_sb.setDecimals(3)
        self.atr_spolfrac_sb.setToolTip('Define the ATR S polarisation fraction, the rest is P polarisation')
        self.atr_spolfrac_sb.setValue(self.settings['ATR S polarisation fraction'])
        self.atr_spolfrac_sb.valueChanged.connect(self.on_atr_spolfrac_sb_changed)
        label = QLabel('ATR S polarisation fraction', self)
        label.setToolTip('Define the S polarisation fraction, the rest is P polarisation')
        form.addRow(label, self.atr_spolfrac_sb)
        #
        # Add a legend option
        #
        self.legend_le = QLineEdit(self)
        self.legend_le.setToolTip('The legend will be used to describe the results in the plot')
        self.legend_le.setText(self.settings['Legend'])
        self.legend_le.textChanged.connect(self.on_legend_le_changed)
        label = QLabel('Powder scenario legend',self)
        label.setToolTip('The legend will be used to describe the results in the plotting tab')
        form.addRow(label, self.legend_le)

        #
        # Final buttons
        #
        hbox = self.add_scenario_buttons()
        form.addRow(hbox)
        #
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        # sort out greying of boxes
        self.change_greyed_out()
        debugger.print('Finished:: initialiser')
        return

    def crystal_density(self):
        '''Find the crystal density from the current reader and return the density'''
        if not self.reader:
            debugger.print('Finished:: crystal_density - no reader')
            return 1.0
        volume = self.reader.volume
        mass = 0.0
        for m in self.reader.masses:
            mass += m
        density = mass / (avogadro_si * volume * 1.0e-24)
        debugger.print('Finished:: crystal_density',density)
        return density

    def openDB_button_clicked(self):
        '''Open a new materials' database'''
        debugger.print('Start:: openDB_button_clicked')
        selfilter = 'Spreadsheet (*.xlsx)'
        filename,myfilter = QFileDialog.getOpenFileName(self,'Open spreadsheet','','Spreadsheet (*.xls);;Spreadsheet (*.xlsx);;All files(*)',selfilter)
        # Process the filename
        if filename != '':
            self.settings['Materials database'] = filename
        self.database_le.setText(self.settings['Materials database'])
        self.DataBase = MaterialsDataBase(self.settings['Materials database'],debug=debugger.state())
        self.materialNames = self.DataBase.getSheetNames()
        self.materialNames.append('Define material manually')
        if not self.settings['Matrix'] in self.materialNames:
            self.settings['Matrix'] = self.materialNames[0]
        index = self.matrix_cb.findText(self.settings['Matrix'], Qt.MatchFixedString)
        self.matrixMaterial = self.DataBase.getMaterial(self.settings['Matrix'])
        materialPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
        self.settings["Matrix permittivity"] = materialPermittivityFunction(0.0)
        self.settings["Matrix density"] = self.matrixMaterial.getDensity()
        self.refreshRequired = True
        self.refresh()
        return

    def on_h_sb_changed(self,value):
        '''Handle a change to the h parameter of the (hkl) surface'''
        debugger.print(self.settings['Legend'],'on_h_sb_changed', value)
        self.refreshRequired = True
        self.settings['Unique direction - h'] = value
        return

    def on_k_sb_changed(self,value):
        '''Handle a change to the k parameter of the (hkl) surface'''
        debugger.print(self.settings['Legend'],'on_k_sb_changed', value)
        self.refreshRequired = True
        self.settings['Unique direction - k'] = value
        return

    def on_l_sb_changed(self,value):
        '''Handle a change to the l parameter of the (hkl) surface'''
        debugger.print(self.settings['Legend'],'on_l_sb_changed', value)
        self.refreshRequired = True
        self.settings['Unique direction - l'] = value
        return

    def on_shape_cb_activated(self,index):
        '''Handle changes to the shape of the powder particles'''
        debugger.print(self.settings['Legend'],'on shape cb activated', index)
        self.refreshRequired = True
        self.settings['Particle shape'] = self.shapes[index]
        if self.settings['Particle shape'] == 'Sphere':
            self.settings['Unique direction - h'] = 0
            self.settings['Unique direction - k'] = 0
            self.settings['Unique direction - l'] = 0
        self.change_greyed_out()
        return

    def on_methods_cb_activated(self,index):
        '''Handle changes in the calculation method for the effective medium theory'''
        debugger.print(self.settings['Legend'],'on methods cb activated', index)
        self.refreshRequired = True
        self.settings['Effective medium method'] = self.methods[index]
        if self.settings['Effective medium method'] == 'Mie':
            self.settings['Particle shape'] = 'Sphere'
        elif self.settings['Effective medium method'] == 'Anisotropic-Mie':
            self.settings['Particle shape'] = 'Sphere'
        elif self.settings['Effective medium method'] == 'Maxwell-Garnett':
            self.settings['Particle size distribution sigma(mu)'] = 0.0
        elif self.settings['Effective medium method'] == 'Bruggeman':
            self.settings['Particle size distribution sigma(mu)'] = 0.0
        elif self.settings['Effective medium method'] == 'Averaged Permittivity':
            self.settings['Particle size(mu)'] = 0.0001
            self.settings['Particle size distribution sigma(mu)'] = 0.0
        self.change_greyed_out()
        return

    def on_mf_sb_changed(self,value):
        '''Handle a mass fraction change and update the volume fraction'''
        debugger.print(self.settings['Legend'],'on mass fraction line edit changed', value)
        self.refreshRequired = True
        self.settings['Mass or volume fraction'] = 'mass'
        self.settings['Mass fraction'] =  value/100.0
        self.update_vf_sb()
        return

    def update_vf_sb(self):
        '''Update the volume fraction according to the mass fraction'''
        mf1 = self.settings['Mass fraction']
        mf2 = 1.0 - mf1
        rho1 = self.crystal_density()
        rho2 = self.settings['Matrix density']
        #
        # Avoid overflow through division by 0
        #
        mf2  = max(mf2,1.0E-18)
        rho1 = max(rho1,1.0E-18)
        vf1 = ( 1.0 - self.settings['Bubble volume fraction'] ) * (mf1/mf2)*(rho2/rho1) / ( 1 + (mf1/mf2)*(rho2/rho1))
        self.settings['Volume fraction'] = vf1
        blocking_state = self.vf_sb.signalsBlocked()
        self.vf_sb.blockSignals(True)
        self.vf_sb.setValue(100.0*vf1)
        self.bubble_vf_sb.setRange(0.0, 100.0*(1.0-self.settings['Volume fraction']))
        self.vf_sb.setRange(0.0, 100.0*(1.0-self.settings['Bubble volume fraction']))
        self.vf_sb.blockSignals(blocking_state)
        debugger.print(self.settings['Legend'],'Update_vf_sb')
        debugger.print(self.settings['Legend'],'rho 1', rho1)
        debugger.print(self.settings['Legend'],'rho 2', rho2)
        debugger.print(self.settings['Legend'],'vf 1 ', vf1)
        return

    def on_aoverb_sb_changed(self,value):
        '''Handle a change to the a/b ratio for an ellipsoid'''
        debugger.print(self.settings['Legend'],'on_aoverb_le_changed',value)
        self.refreshRequired = True
        self.settings['Ellipsoid a/b'] = value
        return

    def on_legend_le_changed(self,text):
        '''Handle a legend change'''
        debugger.print(self.settings['Legend'],'on legend change', text)
        self.refreshRequired = True
        self.settings['Legend'] = text
        return

    def on_sigma_sb_changed(self,value):
        '''Handle a particle size distribution change'''
        debugger.print(self.settings['Legend'],'on sigma line edit changed', value)
        self.refreshRequired = True
        self.settings['Particle size distribution sigma(mu)'] = value
        return

    def on_size_sb_changed(self,value):
        '''Handle a particle size change'''
        debugger.print(self.settings['Legend'],'on size line edit changed', value)
        self.refreshRequired = True
        self.settings['Particle size(mu)'] = value
        return

    def on_vf_sb_changed(self,value):
        '''Handle a volume fraction change, alter the mass fraction accordingly'''
        debugger.print(self.settings['Legend'],'Start:: on_vf_sb_changed', value)
        self.refreshRequired = True
        self.settings['Mass or volume fraction'] = 'volume'
        self.settings['Volume fraction'] = value/100.0
        self.update_mf_sb()
        debugger.print(self.settings['Legend'],'Finished:: on_vf_sb_changed', value)
        return

    def update_mf_sb(self):
        '''Update the mass fraction and according to the volume fraction'''
        debugger.print(self.settings['Legend'],'Start:: update_mf_sb')
        vf1 = self.settings['Volume fraction']
        vf2 = 1.0 - vf1 - self.settings['Bubble volume fraction']
        rho1 = self.crystal_density()
        rho2 = self.settings['Matrix density']
        mf1 = rho1*vf1 / ( rho1*vf1 + rho2*vf2 )
        self.settings['Mass fraction'] = mf1
        blocking_state = self.mf_sb.signalsBlocked()
        self.mf_sb.blockSignals(True)
        self.mf_sb.setValue(100.0*mf1)
        self.mf_sb.blockSignals(blocking_state)
        debugger.print(self.settings['Legend'],'Update_mf_sb')
        debugger.print(self.settings['Legend'],'rho 1', rho1)
        debugger.print(self.settings['Legend'],'rho 2', rho2)
        debugger.print(self.settings['Legend'],'mf 1 ', mf1)
        debugger.print(self.settings['Legend'],'Finished:: update_mf_sb')
        return

    def on_matrix_cb_activated(self,index):
        '''Handle a change to the support matrix supplied by the materials' database'''
        debugger.print(self.settings['Legend'],'on matrix combobox activated', index)
        debugger.print(self.settings['Legend'],'on matrix combobox activated', self.matrix_cb.currentText())
        # We will need to recalculate everything for a new support matrix
        self.refreshRequired = True
        # matrix is the name of the sheet in the database
        matrix = self.matrix_cb.currentText()
        # Make some of the widgets quiet as we update them
        m_blocking = self.matrix_cb.signalsBlocked()
        d_blocking = self.density_sb.signalsBlocked()
        r_blocking = self.permittivity_r_sb.signalsBlocked()
        i_blocking = self.permittivity_i_sb.signalsBlocked()
        self.matrix_cb.blockSignals(True)
        self.density_sb.blockSignals(True)
        self.permittivity_r_sb.blockSignals(True)
        self.permittivity_i_sb.blockSignals(True)
        # Store the new matrix material name
        self.settings['Matrix'] = matrix
        if 'Define material manually' in matrix:
            # The manual option has been chosen, so create a new material with the right permittivity and density
            permittivityObject = DielectricFunction.ConstantScalar(self.settings['Matrix permittivity'])
            self.matrixMaterial = Material('manual',permittivityObject=permittivityObject,density=self.settings['Matrix density'])
            self.matrixPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
        else:
            # Read the material information for permittivity and density from the data base
            self.matrixMaterial = self.DataBase.getMaterial(self.settings['Matrix'])
            self.settings['Matrix density'] = self.matrixMaterial.getDensity()
            self.matrixPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
            # The permittivity may be frequency dependent, show the value at 0 cm-1
            self.settings['Matrix permittivity'] = self.matrixPermittivityFunction(0.0)
            self.density_sb.setValue(self.settings['Matrix density'])
        # Update the matrix material information
        text = self.matrixMaterial.getInformation()
        self.matrix_info_le.setText(text)
        # Update the values of the real and imaginary permittivity
        self.permittivity_r_sb.setValue(np.real(self.settings['Matrix permittivity']))
        self.permittivity_i_sb.setValue(np.imag(self.settings['Matrix permittivity']))
        # This is a new material so reset the volume / mass fractions
        # volume fraction takes precedence
        if self.settings['Mass or volume fraction'] == 'volume':
            self.settings['Volume fraction'] = 0.1
            self.update_mf_sb()
            self.update_vf_sb()
        else:
            self.settings['Mass fraction'] = 0.1
            self.update_vf_sb()
            self.update_mf_sb()
        # Restore the signal settings on the widgets
        self.matrix_cb.blockSignals(m_blocking)
        self.density_sb.blockSignals(d_blocking)
        self.permittivity_r_sb.blockSignals(r_blocking)
        self.permittivity_i_sb.blockSignals(i_blocking)
        return

    def on_density_sb_changed(self,value):
        '''Handle a change to the matrix density'''
        self.settings['Matrix density'] = value
        # update the matrix density
        self.matrixMaterial.setDensity(value)
        # Force the matrix to be defined manually
        self.settings['Matrix'] = self.materialNames[-1]
        index = self.matrix_cb.findText(self.settings['Matrix'], Qt.MatchFixedString)
        if index >=0:
            self.matrix_cb.setCurrentIndex(index)
        else:
            print('support matrix index was not 0',matrix)
        # volume fraction taked precedence
        if self.settings['Mass or volume fraction'] == 'volume':
            self.update_mf_sb()
            self.update_vf_sb()
        else:
            self.update_vf_sb()
            self.update_mf_sb()
        debugger.print(self.settings['Legend'],'on density line edit changed', value)
        self.refreshRequired = True
        return

    def on_bubble_vf_sb_changed(self,value):
        '''Handle a change to the bubble volume fraction'''
        self.settings['Bubble volume fraction'] = value/100.0
        if self.settings['Mass or volume fraction'] == 'volume':
            self.update_mf_sb()
        else:
            self.update_vf_sb()
        debugger.print(self.settings['Legend'],'on bubble volume fraction changed', value)
        self.refreshRequired = True
        return

    def on_bubble_radius_sb_changed(self,value):
        self.settings['Bubble radius'] = value
        debugger.print(self.settings['Legend'],'on bubble raduys line edit changed', value)
        self.refreshRequired = True
        return

    def on_permittivity_i_sb_changed(self,value):
        real = np.real(self.settings['Matrix permittivity'])
        self.settings['Matrix permittivity'] = complex(real,value)
        newPermittivityObject = DielectricFunction.ConstantScalar(self.settings['Matrix permittivity'])
        self.matrixMaterial.setPermittivityObject(newPermittivityObject)
        debugger.print(self.settings['Legend'],'on imaginary permittivity line edit changed', value)
        self.refreshRequired = True
        return

    def on_permittivity_r_sb_changed(self,value):
        imaginary = np.imag(self.settings['Matrix permittivity'])
        self.settings['Matrix permittivity'] = complex(value,imaginary)
        newPermittivityObject = DielectricFunction.ConstantScalar(self.settings['Matrix permittivity'])
        self.matrixMaterial.setPermittivityObject(newPermittivityObject)
        debugger.print(self.settings['Legend'],'on permittivity line edit changed', value)
        self.refreshRequired = True
        return

    def on_atr_index_sb_changed(self,value):
        self.settings['ATR material refractive index'] = value
        debugger.print(self.settings['Legend'],'on atr index line edit changed', value)
        self.refreshRequired = True
        return

    def on_atr_incident_ang_sb_changed(self,value):
        self.settings['ATR theta'] = value
        debugger.print(self.settings['Legend'],'on atr incident angle line edit changed', value)
        self.refreshRequired = True
        return

    def on_atr_spolfrac_sb_changed(self,value):
        self.settings['ATR S polarisation fraction'] = value
        debugger.print(self.settings['Legend'],'on atr spolfraction line edit changed', value)
        self.refreshRequired = True
        return

    def change_greyed_out(self):
        # Have a look through the settings and see if we need to grey anything out
        debugger.print(self.settings['Legend'],'Start:: change_greyed_out')
        method = self.settings['Effective medium method']
        if method == 'Mie' or method == 'Anisotropic-Mie':
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
        elif method == 'Averaged Permittivity':
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
        elif method == 'Maxwell-Garnett' or method == 'Bruggeman':
            self.size_sb.setEnabled(True)
            self.sigma_sb.setEnabled(False)
            self.shape_cb.setEnabled(True)
            for i,shape in enumerate(self.shapes):
                self.shape_cb.model().item(i).setEnabled(True)
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
        debugger.print(self.settings['Legend'],'Finished:: change_greyed_out')
        return

    def calculate(self, vs_cm1):
        """Calculate the powder absorption for the range of frequencies in vs_cm1"""
        # Only allow a calculation if the plottingTab is defined
        debugger.print(self.settings['Legend'],'Start:: calculate')
        if not self.calculationRequired:
            debugger.print(self.settings['Legend'],'Finished:: calculate - immediate return because calculationRequired false')
            return None
        if self.notebook.plottingTab is None:
            debugger.print(self.settings['Legend'],'Finished:: calculate - immediate return because plottingTab unavailable')
            return None
        if self.reader is None:
            debugger.print(self.settings['Legend'],'Finished:: calculate - immediate return because reader unavailable')
            return None
        debugger.print(self.settings['Legend'],'calculate - number of frequencies',len(vs_cm1))
        cell = self.reader.unit_cells[-1]
        shape = self.settings['Particle shape']
        hkl = [self.settings['Unique direction - h'], self.settings['Unique direction - k'], self.settings['Unique direction - l']]
        if shape == 'Ellipsoid':
            self.direction = cell.convert_abc_to_xyz(hkl)
            self.depolarisation = Calculator.initialise_ellipsoid_depolarisation_matrix(self.direction,self.aoverb)
        elif shape == 'Plate':
            self.direction = cell.convert_hkl_to_xyz(hkl)
            self.depolarisation = Calculator.initialise_plate_depolarisation_matrix(self.direction)
        elif shape == 'Needle':
            self.direction = cell.convert_abc_to_xyz(hkl)
            self.depolarisation = Calculator.initialise_needle_depolarisation_matrix(self.direction)
        else:
            self.depolarisation = Calculator.initialise_sphere_depolarisation_matrix()
            self.direction = np.array( [] )
        self.direction = self.direction / np.linalg.norm(self.direction)
        # Get the crystal permittivity function from the settings tab
        crystalPermittivity = self.notebook.settingsTab.getCrystalPermittivity(vs_cm1)
        # Allocate space for the shared memory, we need twice as much as we have a complex data type
        shared_array_base = Array(ctypes.c_double, 18)
        previous_solution_shared = np.ctypeslib.as_array(shared_array_base.get_obj())
        # Convert the space allocated to complex
        previous_solution_shared.dtype = np.complex128
        # Reshape the array and fill everything with zero's
        previous_solution_shared = previous_solution_shared.reshape(3,3)
        previous_solution_shared.fill(0.0+0.0j)
        # Prepare parallel call parameters for the loop over frequencies, methods, volume fractions
        # The concentration is defined in the plottingTab, which may not exist yet
        concentration = self.notebook.plottingTab.settings['cell concentration']
        # Set the material parameters
        method = self.settings['Effective medium method'].lower()
        volume_fraction = self.settings['Volume fraction']
        particle_size_mu = self.settings['Particle size(mu)']
        particle_sigma_mu = self.settings['Particle size distribution sigma(mu)']
        shape = self.settings['Particle shape'].lower()
        atr_refractive_index = self.settings['ATR material refractive index']
        atr_theta = self.settings['ATR theta']
        atr_spolfraction = self.settings['ATR S polarisation fraction']
        bubble_vf = self.settings['Bubble volume fraction']
        bubble_radius = self.settings['Bubble radius']
        # Use the pool of processors already available
        # define a partial function to use with the pool
        partial_function = partial(Calculator.solve_effective_medium_equations, method,volume_fraction,particle_size_mu,particle_sigma_mu,self.matrixPermittivityFunction,shape,self.depolarisation,concentration,atr_refractive_index,atr_theta,atr_spolfraction,bubble_vf,bubble_radius,previous_solution_shared)
        if self.notebook.pool is None:
            self.notebook.startPool()
        debugger.print('About to use the pool to calculate effective medium equations')
        results = []
        indices = range(len(vs_cm1))
        for result in self.notebook.pool.map(partial_function, zip(vs_cm1,crystalPermittivity), chunksize=20):
            results.append(result)
            self.notebook.progressbars_update()
        QCoreApplication.processEvents()
        self.realPermittivity = []
        self.imagPermittivity = []
        self.absorptionCoefficient = []
        self.molarAbsorptionCoefficient = []
        self.sp_atr = []
        self.vs_cm1 = []
        debugger.print('Extracting results')
        for v,method,size_mu,size_sigma,shape,data,trace, absorption_coefficient,molar_absorption_coefficient,spatr in results:
             self.realPermittivity.append(np.real(trace))
             self.imagPermittivity.append(np.imag(trace))
             self.absorptionCoefficient.append(absorption_coefficient)
             self.molarAbsorptionCoefficient.append(molar_absorption_coefficient)
             self.sp_atr.append(spatr)
             self.vs_cm1.append(v)
        self.calculationRequired = False
        QCoreApplication.processEvents()
        debugger.print(self.settings['Legend'],'Finished:: calculate')
        return

    def get_result(self, vs_cm1, plot_type):
        """Return a particular result"""
        debugger.print(self.settings['Legend'],'Start:: get_result')
        self.get_results(vs_cm1)
        debugger.print(self.settings['Legend'],'Finished:: get_result')
        if plot_type   == 'Powder Molar Absorption':
            return self.molarAbsorptionCoefficient
        elif plot_type == 'Powder Absorption':
            return self.absorptionCoefficient
        elif plot_type == 'Powder Real Permittivity':
            return self.realPermittivity
        elif plot_type == 'Powder Imaginary Permittivity':
            return self.imagPermittivity
        elif plot_type == 'Powder ATR':
            return self.sp_atr
        else:
            # print('Error in returning result from PowderScenarioTab: ',plot_type)
            return None


    def get_results(self, vs_cm1):
        """Return the results of the effective medium theory calculation"""
        debugger.print(self.settings['Legend'],'Start:: get_results', len(vs_cm1))
        if len(vs_cm1) > 0 and ( self.refreshRequired or len(self.vs_cm1) != len(vs_cm1) or self.vs_cm1[0] != vs_cm1[0] or self.vs_cm1[1] != vs_cm1[1] ) :
            debugger.print(self.settings['Legend'],'get_results recalculating')
            self.refresh()
            self.calculate(vs_cm1)
        else:
            debugger.print(self.settings['Legend'],'get_results no need for recalculation')
            self.notebook.progressbars_update(increment=len(vs_cm1))
        debugger.print(self.settings['Legend'],'Finished:: get_results', len(vs_cm1))
        return


    def refresh(self,force=False):
        '''Refresh the GUI interface with up to date values
           if force is True then a refresh is forced whatever the state of the scenario
        '''
        debugger.print(self.settings['Legend'],'Start:: refresh, force =', force)
        if not self.refreshRequired and not force:
            debugger.print(self.settings['Legend'],'Finished:: refresh aborted', self.refreshRequired,force)
            return
        # Force a recalculation
        self.calculationRequired = True
        # First see if we can get the reader from the mainTab
        self.reader = self.notebook.mainTab.reader
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # use the settings values to initialise the widgets
        # Update the database name
        self.database_le.setText(self.settings['Materials database'])
        self.DataBase = MaterialsDataBase(self.settings['Materials database'],debug=debugger.state())
        # Update the possible matrix material names from the database
        self.materialNames = self.DataBase.getSheetNames()
        self.materialNames.append('Define material manually')
        # Set the matrix support material
        if 'Define material manually' in self.settings['Matrix']:
            permittivityObject = DielectricFunction.ConstantScalar(self.settings['Matrix permittivity'])
            self.matrixMaterial = Material('manual',permittivityObject=permittivityObject,density=self.settings['Matrix density'])
            self.matrixPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
        else:
           # Check to see if the permittivity or the density have been changed from the database values
           # if they have then force a 'Manual' setting.  This should allow backward compatibility with
            if self.settings['Matrix'] in self.materialNames:
                self.matrixMaterial = self.DataBase.getMaterial(self.settings['Matrix'])
                self.matrixPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
            else:
                print('Error: matrix ',self.settings['Matrix'],' not available in database')
                print('       available materials are:', self.materialNames)
                exit()
            material = self.DataBase.getMaterial(self.settings['Matrix'])
            materialPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
            materialPermittivity = materialPermittivityFunction(0.0)
            materialDensity = self.matrixMaterial.getDensity()
            if (np.abs(materialPermittivity - self.settings['Matrix permittivity'])) or (np.abs(materialDensity-self.settings['Matrix density']) > 1.0E-8):
                # It looks like the permittivity or density has been set manually
                # Creat a constant permittivity material and set the Matrix to manual
                self.settings['Matrix'] = 'Define material manually'
                permittivityObject = DielectricFunction.ConstantScalar(self.settings['Matrix permittivity'])
                self.matrixMaterial = Material('manual',permittivityObject=permittivityObject,density=self.settings['Matrix density'])
                self.matrixPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
            else:
                # Use the database material
                self.matrixMaterial = material
                self.matrixPermittivityFunction = self.matrixMaterial.getPermittivityFunction()
        # Reset the matrix combo box with new names
        self.matrix_cb.clear()
        self.matrix_cb.addItems(self.materialNames)
        index = self.matrix_cb.findText(self.settings['Matrix'], Qt.MatchFixedString)
        self.matrix_cb.setCurrentIndex(index)
        # Update the matrix material information
        text = self.matrixMaterial.getInformation()
        self.matrix_info_le.setText(text)
        # Set the matrix density widget
        self.density_sb.setValue(self.settings['Matrix density'])
        # Set the matrix permittivity widget
        self.permittivity_r_sb.setValue(np.real(self.settings['Matrix permittivity']))
        self.permittivity_i_sb.setValue(np.imag(self.settings['Matrix permittivity']))
        self.bubble_vf_sb.setValue(100*self.settings['Bubble volume fraction'])
        self.bubble_radius_sb.setValue(self.settings['Bubble radius'])
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
        self.legend_le.setText(self.settings['Legend'])
        self.aoverb = self.settings['Ellipsoid a/b']
        self.atr_index_sb.setValue(self.settings['ATR material refractive index'])
        self.atr_incident_ang_sb.setValue(self.settings['ATR theta'])
        self.atr_spolfrac_sb.setValue(self.settings['ATR S polarisation fraction'])
        self.change_greyed_out()
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.refreshRequired = False
        debugger.print(self.settings['Legend'],'Finished:: refresh, force =', force)
        return
