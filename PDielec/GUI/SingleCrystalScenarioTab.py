# -*- coding: utf8 -*-
import os.path
import os
import numpy as np
import copy
import PDielec.Calculator as Calculator
from PyQt5.QtWidgets            import QPushButton, QWidget
from PyQt5.QtWidgets            import QSpacerItem
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
from functools                  import partial
from scipy                      import signal

class SingleCrystalScenarioTab(ScenarioTab):
    def __init__(self, parent, debug=False ):
        ScenarioTab.__init__(self,parent)
        global debugger
        debugger = Debug(debug,'SingleCrystalScenarioTab:')
        debugger.print('Start:: initialiser')
        self.refreshRequired = True
        self.calculationRequired = False
        self.scenarioType = 'Single crystal'
        self.noCalculationsRequired = 1
        self.settings['Scenario type'] = 'Single crystal'
        self.settings['Unique direction - h'] = 0
        self.settings['Unique direction - k'] = 0
        self.settings['Unique direction - l'] = 1
        self.settings['Azimuthal angle'] = 0.0
        self.settings['Angle of incidence'] = 0.0
        self.settings['Superstrate dielectric'] = 1.0
        self.settings['Substrate dielectric'] = 1.0
        self.settings['Superstrate depth'] = 99999.0
        self.settings['Substrate depth'] = 99999.0
        self.settings['Film thickness'] = 10
        self.settings['Thickness unit'] = 'um'
        self.settings['Mode'] = 'Thick slab'
        self.settings['Frequency units'] = 'wavenumber'
        self.settings['Partially incoherent samples'] = 20
        self.settings['Percentage partial incoherence'] = 100
        self.settings['Filter kernel size'] = 1
        self.settings['Filter polynomial size'] = 3
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
        self.mode_cb.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two semi-infinite media; the incident and the crystal),\n Partially incoherent, incoherent and coherent thin films assume there are three media; the incident, the crystal and the substrate. Partially incoherent film mode takes a lot longer than any of the other modes.')
        self.mode_cb.addItems( ['Thick slab','Coherent thin film','Incoherent thin film','Partially incoherent thin film'] )
        self.settings['Mode'] = 'Thick slab'
        self.mode_cb.activated.connect(self.on_mode_cb_activated)
        label = QLabel('Single crystal mode', self)
        label.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two semi-infinite media; the incident and the crystal),\n Partially incoherent, incoherent and coherent thin films assume there are three media; the incident, the crystal and the substrate.  Partially incoherent film mode takes a lot longer to run than any of the other modes.')
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
        self.azimuthal_angle_sb.setToolTip('Define the slab azimuthal angle (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal in the laboratory frame can be seen in the laboratory frame information below')
        self.azimuthal_angle_sb.setRange(-180,360)
        self.azimuthal_angle_sb.setSingleStep(10)
        self.azimuthal_angle_sb.setValue(self.settings['Azimuthal angle'])
        self.azimuthal_angle_sb.valueChanged.connect(self.on_azimuthal_angle_sb_changed)
        label = QLabel('Azimuthal angle',self)
        label.setToolTip('Define the slab azimuthal angle (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal in the laboratory frame can be seen in the laboratory frame information below')
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
        # Set the units used for thickness
        #
        self.thickness_unit_cb = QComboBox(self)
        self.thickness_unit_cb.setToolTip('Set the units to be used for thickness; either nm, um, mm or cm')
        self.thickness_unit_cb.addItems( ['nm','um','mm','cm'] )
        self.settings['Thickness unit'] = 'nm'
        self.thickness_unit_cb.activated.connect(self.on_thickness_units_cb_activated)
        label = QLabel('Set thickness units', self)
        label.setToolTip('Set the units to be used for thickness: either nm, um, mm or cm')
        form.addRow(label, self.thickness_unit_cb)
        #
        # Define the Film thickness
        #
        self.film_thickness_sb = QDoubleSpinBox(self)
        self.film_thickness_sb.setToolTip('Define the thin film thickness in the defined thickness units')
        self.film_thickness_sb.setRange(0,100000)
        self.film_thickness_sb.setSingleStep(1)
        self.film_thickness_sb.setValue(self.settings['Film thickness'])
        self.film_thickness_sb.valueChanged.connect(self.on_film_thickness_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.film_thickness_sb)
        film_thickness_label = QLabel('Film thickness',self)
        film_thickness_label.setToolTip('Define the depth of the thin crystal in the defined thickness units.')
        form.addRow(film_thickness_label, hbox)
        #
        # Define the number of amount of partial incoherence to be included
        #
        self.percentage_partial_incoherence_sb = QSpinBox(self)
        self.percentage_partial_incoherence_sb.setToolTip('Define the maximum percentage changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees')
        self.percentage_partial_incoherence_sb.setRange(0,100)
        self.percentage_partial_incoherence_sb.setSingleStep(1)
        self.percentage_partial_incoherence_sb.setValue(self.settings['Percentage partial incoherence'])
        self.percentage_partial_incoherence_sb.valueChanged.connect(self.on_percentage_partial_incoherence_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.percentage_partial_incoherence_sb)
        percentage_partial_incoherence_label = QLabel('Percentage of partial incoherence',self)
        percentage_partial_incoherence_label.setToolTip('Define the maximum percentage changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees')
        form.addRow(percentage_partial_incoherence_label, hbox)
        #
        # Define the number of samples to be used in calculating the incoherent case
        #
        self.partially_incoherent_samples_sb = QSpinBox(self)
        self.partially_incoherent_samples_sb.setToolTip('Define the number of samples to be used in the calculation of an incoherent spectrum.\nA large number of samples will take a long time but will give smoother results.')
        self.partially_incoherent_samples_sb.setRange(0,10000)
        self.partially_incoherent_samples_sb.setSingleStep(1)
        self.partially_incoherent_samples_sb.setValue(self.settings['Partially incoherent samples'])
        self.partially_incoherent_samples_sb.valueChanged.connect(self.on_partially_incoherent_samples_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.partially_incoherent_samples_sb)
        partially_incoherent_samples_label = QLabel('Number of partially incoherent samples',self)
        partially_incoherent_samples_label.setToolTip('Define the number of samples to be used in the calculation of an incoherent spectrum.\nA large number of samples will take a long time but will give smoother results.')
        form.addRow(partially_incoherent_samples_label, hbox)
        #
        # Define the kernel size for the smoothing of incoherent spectra
        #
        self.partially_incoherent_kernel_sb = QSpinBox(self)
        self.partially_incoherent_kernel_sb.setToolTip('Define the kernel size for the smoothing of incoherent spectra (must be an odd number)\nIf the kernel size is less than 3, no smoothing is done.\nThe larger the number, the smoother the spectrum but beware of too much smoothing.')
        self.partially_incoherent_kernel_sb.setRange(1,1001)
        self.partially_incoherent_kernel_sb.setSingleStep(2)
        self.partially_incoherent_kernel_sb.setValue(self.settings['Filter kernel size'])
        self.partially_incoherent_kernel_sb.valueChanged.connect(self.on_partially_incoherent_kernel_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.partially_incoherent_kernel_sb)
        partially_incoherent_kernel_label = QLabel('Incoherent smoothing filter, kernel size',self)
        partially_incoherent_kernel_label.setToolTip('Define the kernel size for the smoothing of incoherent spectra (must be an odd number)\nIf the kernel size is less than 3, no smoothing is done.\nThe larger the number, the smoother the spectrum but beware of too much smoothing. ')
        form.addRow(partially_incoherent_kernel_label, hbox)
        #
        # Define the polynomial degree for the smoothing function
        #
        self.partially_incoherent_polynomial_sb = QSpinBox(self)
        self.partially_incoherent_polynomial_sb.setToolTip('Define the maximum degree of polynomial to be used in the smoothing filter')
        self.partially_incoherent_polynomial_sb.setRange(2,10)
        self.partially_incoherent_polynomial_sb.setSingleStep(1)
        self.partially_incoherent_polynomial_sb.setValue(self.settings['Filter polynomial size'])
        self.partially_incoherent_polynomial_sb.valueChanged.connect(self.on_partially_incoherent_polynomial_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.partially_incoherent_polynomial_sb)
        partially_incoherent_polynomial_label = QLabel('Incoherent smoothing filter, polynomial',self)
        partially_incoherent_polynomial_label.setToolTip('Define the maximum degree of polynomial to be used in the smoothing filter.')
        form.addRow(partially_incoherent_polynomial_label, hbox)
        #
        # Add a legend option
        #
        self.legend_le = QLineEdit(self)
        self.legend_le.setToolTip('The legend will be used to describe the results in the plot')
        self.legend_le.setText(self.settings['Legend'])
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
        #
        # Spacer item to soak up space
        #
        verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        form.addItem(verticalSpacer)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        debugger.print('Finished:: initialiser')

    def on_partially_incoherent_kernel_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_partially_incoherent_kernel_sb_changed', value)
        self.refreshRequired = True
        self.settings['Filter kernel size'] = value
        return

    def on_partially_incoherent_polynomial_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_partially_incoherent_polynomial_sb_changed', value)
        self.refreshRequired = True
        self.settings['Filter polynomial size'] = value
        return

    def on_partially_incoherent_samples_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_partially_incoherent_samples_sb_changed', value)
        self.refreshRequired = True
        self.settings['Partially incoherent samples'] = value
        self.noCalculationsRequired = value
        return

    def on_percentage_partial_incoherence_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_percentage_partial_incoherence_sb_changed', value)
        self.refreshRequired = True
        self.settings['Percentage partial incoherence'] = value
        return

    def on_film_thickness_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_film_thickness_sb_changed', value)
        self.refreshRequired = True
        self.settings['Film thickness'] = value
        return

    def on_superstrate_dielectric_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_superstrate_dielectric_sb_changed', value)
        self.refreshRequired = True
        self.settings['Superstrate dielectric'] = value
        return

    def on_substrate_dielectric_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_substrate_dielectric_sb_changed', value)
        self.refreshRequired = True
        self.settings['Substrate dielectric'] = value
        return

    def on_azimuthal_angle_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_azimuthal_angl_sb_changed', value)
        self.refreshRequired = True
        self.settings['Azimuthal angle'] = value
        self.calculate_euler_angles()
        return

    def on_angle_of_incidence_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_angle_of_incidence_sb_changed', value)
        self.refreshRequired = True
        self.settings['Angle of incidence'] = value
        return

    def on_h_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_h_sb_changed', value)
        self.refreshRequired = True
        self.settings['Unique direction - h'] = value
        self.calculate_euler_angles()
        return

    def on_k_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_k_sb_changed', value)
        self.refreshRequired = True
        self.settings['Unique direction - k'] = value
        self.calculate_euler_angles()
        return

    def on_l_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_l_sb_changed', value)
        self.refreshRequired = True
        self.settings['Unique direction - l'] = value
        self.calculate_euler_angles()
        return

    def refresh(self,force=False):
        debugger.print(self.settings['Legend'],'Start:: refresh, force =', force)
        if not self.refreshRequired and not force :
            debugger.print(self.settings['Legend'],'Finished:: refreshing widget aborted', self.refreshRequired,force)
            return
        # Force recalculation
        self.calculationRequired = True
        # Change any greyed out items
        self.greyed_out()
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Now refresh values
        # First check the mode and if it is Incoherent increase the number of spectra to be calculated
        index = self.mode_cb.findText(self.settings['Mode'], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        if self.settings['Mode'] == 'Partially incoherent thin film':
            self.noCalculationsRequired = self.settings['Partially incoherent samples']
        else:
            self.settings['Percentage partial incoherence'] = 0
            self.noCalculationsRequired = 1
        # Update the Legend
        self.legend_le.setText(self.settings['Legend'])
        # Update hkl
        self.h_sb.setValue(self.settings['Unique direction - h'])
        self.k_sb.setValue(self.settings['Unique direction - k'])
        self.l_sb.setValue(self.settings['Unique direction - l'])
        # Update angles
        self.azimuthal_angle_sb.setValue(self.settings['Azimuthal angle'])
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        # Update super and substrate permittivities
        self.superstrate_dielectric_sb.setValue(self.settings['Superstrate dielectric'])
        self.substrate_dielectric_sb.setValue(self.settings['Substrate dielectric'])
        # Set the film thickness and the units of thickness
        self.film_thickness_sb.setValue(self.settings['Film thickness'])
        index = self.thickness_unit_cb.findText(self.settings['Thickness unit'], Qt.MatchFixedString)
        self.thickness_unit_cb.setCurrentIndex(index)
        # For partial incoherent case, set percentage variation of angles and thickness and the number of samples
        self.percentage_partial_incoherence_sb.setValue(self.settings['Percentage partial incoherence'])
        self.partially_incoherent_samples_sb.setValue(self.settings['Partially incoherent samples'])
        # For partial incoherent case, set the smoothing parameters
        self.partially_incoherent_kernel_sb.setValue(self.settings['Filter kernel size'])
        self.partially_incoherent_polynomial_sb.setValue(self.settings['Filter polynomial size'])
        # Check to see if there is a reader, if there is set up the cell
        self.reader = self.notebook.mainTab.reader
        if self.reader is not None:
            self.cell = self.reader.unit_cells[-1]
            # Refresh the widgets that depend on the reader
            self.reader = self.notebook.reader
            self.calculate_euler_angles()
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        self.refreshRequired = False
        debugger.print(self.settings['Legend'],'Finished:: refresh, force =', force)
        return

    def on_thickness_units_cb_activated(self, index):
        debugger.print(self.settings['Legend'],'Start:: on_thickness_units_cb_activated')
        if index == 0:
            self.settings['Thickness unit'] = 'nm'
        elif index == 1:
            self.settings['Thickness unit'] = 'um'
        elif index == 2:
            self.settings['Thickness unit'] = 'mm'
        elif index == 3:
            self.settings['Thickness unit'] = 'cm'
        else:
            pass
        self.refreshRequired = True
        self.refresh()
        self.refreshRequired = True
        debugger.print(self.settings['Legend'],'Thickness unit changed to ', self.settings['Thickness unit'])
        debugger.print(self.settings['Legend'],'Finished:: on_thickness_unit_cb_activated')
        return

    def on_mode_cb_activated(self, index):
        debugger.print(self.settings['Legend'],'Start:: on_mode_cb_activated')
        if index == 0:
            self.settings['Mode'] = 'Thick slab'
            self.noCalculationsRequired = 1
        elif index == 1:
            self.settings['Mode'] = 'Coherent thin film'
            self.noCalculationsRequired = 1
        elif index == 2:
            self.settings['Mode'] = 'Incoherent thin film'
            self.noCalculationsRequired = 1
        elif index == 3:
            self.settings['Mode'] = 'Partially incoherent thin film'
            self.noCalculationsRequired = self.settings['Partially incoherent samples']
            if self.settings['Percentage partial incoherence'] == 0:
                self.settings['Percentage partial incoherence'] = 20
        self.refreshRequired = True
        self.refresh()
        self.refreshRequired = True
        debugger.print(self.settings['Legend'],'Mode changed to ', self.settings['Mode'])
        debugger.print(self.settings['Legend'],'Finished:: on_mode_cb_activated')
        return

    def partially_incoherent_calculator( self,
                            superstrateDielectricFunction,
                            substrateDielectricFunction,
                            crystalPermittivityFunction,
                            superstrateDepth,
                            substrateDepth,
                            crystalDepth,
                            mode,
                            theta,
                            phi,
                            psi,
                            angleOfIncidence):
        """ Calculates the incoherent component of light reflectance and transmission
            by sampling the path length in the incident medium """
        debugger.print(self.settings['Legend'],'Start:: partially_incoherent_calculator')
        #
        # Zero the arrays we will need
        #
        size = len(self.vs_cm1)
        av_p_reflectance = np.zeros(size) 
        av_s_reflectance = np.zeros(size) 
        av_p_transmittance = np.zeros(size) 
        av_s_transmittance = np.zeros(size) 
        av_s_absorbtance = np.zeros(size) 
        av_p_absorbtance = np.zeros(size) 
        av_epsilon = np.zeros((size,3,3),dtype=np.complex128)
        #
        # Loop over the number of samples requred
        #
        d = crystalDepth
        t = theta
        p  = phi
        q  = psi
        a  = angleOfIncidence
        for s in range(self.settings["Partially incoherent samples"]):
            crystalDepth = d + d*( -1 + 2*np.random.rand())*self.settings['Percentage partial incoherence']/100.0
            theta = t + np.pi/2.0*( -1 +2*np.random.rand())* self.settings['Percentage partial incoherence'] /100.0
            phi   = p + np.pi/2.0*( -1 +2*np.random.rand())* self.settings['Percentage partial incoherence'] /100.0
            psi = q + np.pi/2.0*( -1 +2*np.random.rand())* self.settings['Percentage partial incoherence'] /100.0
            angleOfIncidence = a + np.pi/2.0*( -1 +2*np.random.rand())* self.settings['Percentage partial incoherence'] /100.0
            ( p_reflectance, 
            s_reflectance, 
            p_transmittance, 
            s_transmittance, 
            s_absorbtance, 
            p_absorbtance, 
            epsilon) = self.coherent_calculator(
                                       superstrateDielectricFunction,
                                       substrateDielectricFunction,
                                       crystalPermittivityFunction,
                                       superstrateDepth,
                                       substrateDepth,
                                       crystalDepth,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence)
            av_p_reflectance   += np.array(p_reflectance) / self.settings['Partially incoherent samples']
            av_s_reflectance   += np.array(s_reflectance) / self.settings['Partially incoherent samples']
            av_p_transmittance += np.array(p_transmittance) / self.settings['Partially incoherent samples']
            av_s_transmittance += np.array(s_transmittance) / self.settings['Partially incoherent samples']
            av_s_absorbtance   += np.array(p_absorbtance) / self.settings['Partially incoherent samples']
            av_p_absorbtance   += np.array(s_absorbtance) / self.settings['Partially incoherent samples']
            av_epsilon         += np.array(epsilon) / self.settings['Partially incoherent samples']
        # Only apply the smoothing filter if the kernel is larger than 2
        k = self.settings['Filter kernel size']
        if k > 2:
            p = self.settings['Filter polynomial size']
            av_p_reflectance   = signal.savgol_filter(av_p_reflectance, k, p, mode="nearest")
            av_p_transmittance = signal.savgol_filter(av_p_transmittance, k, p, mode="nearest")
            av_p_absorbtance   = signal.savgol_filter(av_p_absorbtance, k, p, mode="nearest")
            av_s_reflectance   = signal.savgol_filter(av_s_reflectance, k, p, mode="nearest")
            av_s_transmittance = signal.savgol_filter(av_s_transmittance, k, p, mode="nearest")
            av_s_absorbtance   = signal.savgol_filter(av_s_absorbtance, k, p, mode="nearest")
        return (  av_p_reflectance.tolist(), av_s_reflectance.tolist(), av_p_transmittance.tolist(), av_s_transmittance.tolist(), av_p_absorbtance.tolist(), av_s_absorbtance.tolist(), av_epsilon.tolist() )


    def coherent_calculator( self,
                            superstrateDielectricFunction,
                            substrateDielectricFunction,
                            crystalPermittivityFunction,
                            superstrateDepth,
                            substrateDepth,
                            crystalDepth,
                            mode,
                            theta,
                            phi,
                            psi,
                            angleOfIncidence):
        """ Calculates the coherent component of light reflectance and transmission """
        debugger.print(self.settings['Legend'],'Entering the coherent_calculator function')
        #
        # Initialise the partial function to pass through to the pool
        #
        crystalIncoherence = self.settings["Percentage partial incoherence"] * np.pi / 100.0 
        partial_function = partial(Calculator.solve_single_crystal_equations,
                                       superstrateDielectricFunction,
                                       substrateDielectricFunction,
                                       crystalPermittivityFunction,
                                       superstrateDepth,
                                       substrateDepth,
                                       crystalDepth,
                                       crystalIncoherence,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence)
        results = []
        # About to call
        debugger.print(self.settings['Legend'],'About to calculate single crystal scenario using pool')
        #begin serial version 
        #for v in self.vs_cm1:
        #   results.append(Calculator.solve_single_crystal_equations( superstrateDielectricFunction, substrateDielectricFunction, crystalPermittivityFunction, superstrateDepth, substrateDepth, crystalDepth, crystalIncoherence, mode, theta, phi, psi, angleOfIncidence,v))
           #print('results',results[0:4])
        #end serial version
        if self.notebook.pool is None:
            self.notebook.startPool()
        for result in self.notebook.pool.imap(partial_function, self.vs_cm1, chunksize=20):
            self.notebook.progressbars_update()
            results.append(result)
        QCoreApplication.processEvents()
        # Initialise result variables
        vs_cm1= []
        p_reflectance = []
        s_reflectance = []
        p_transmittance = []
        s_transmittance = []
        p_absorbtance = []
        s_absorbtance = []
        epsilon = []
        debugger.print(self.settings['Legend'],'About to extract results for single crystal scenario')
        for v,r,R,t,T,eps in results:
            p_reflectance.append(R[0]+R[2])
            s_reflectance.append(R[1]+R[3])
            if self.settings['Mode'] == 'Thick slab':
                p_absorbtance.append(1.0 - R[0]-R[2])
                s_absorbtance.append(1.0 - R[1]-R[3])
                p_transmittance.append(np.zeros_like(T[0]))
                s_transmittance.append(np.zeros_like(T[1]))
            elif self.settings['Mode'] == 'Incoherent thin film':
                p_transmittance.append(T[0])
                s_transmittance.append(T[1])
                p_absorbtance.append(1.0 - p_reflectance[-1] - p_transmittance[-1])
                s_absorbtance.append(1.0 - s_reflectance[-1] - s_transmittance[-1])
            else:
                p_transmittance.append(T[0])
                s_transmittance.append(T[1])
                p_absorbtance.append(1.0 - R[0]-R[2]-T[0])
                s_absorbtance.append(1.0 - R[1]-R[3]-T[1])
            epsilon.append(eps)
        debugger.print(self.settings['Legend'],'Finished the coherent_calculator function')
        return ( p_reflectance, s_reflectance, p_transmittance, s_transmittance, p_absorbtance, s_absorbtance, epsilon )

    def calculate(self,vs_cm1):
        debugger.print(self.settings['Legend'],'Start:: calculate - number of frequencies',len(vs_cm1))
        if not self.calculationRequired:
            debugger(self.settings['Legend'],'Finished:: calculate aborted because calculationRequired false')
            return
        QCoreApplication.processEvents()
        self.vs_cm1 = vs_cm1
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['Program']
        filename = self.notebook.mainTab.getFullFileName()
        if self.reader is None:
            debugger.print(self.settings['Legend'],'Finished:: Calculate aborting - no reader')
            return
        if program == '':
            debugger.print(self.settings['Legend'],'Finished:: Calculate aborting - no program')
            return
        if filename == '':
            debugger.print(self.settings['Legend'],'Finished:: Calculate aborting - no file')
            return
        # Make sure the filter kernel size is odd, if not make it so and update the GUI
        if self.settings['Filter kernel size'] % 2 == 0:
            self.settings['Filter kernel size'] += 1
            self.partially_incoherent_kernel_sb.setValue(self.settings['Filter kernel size'])
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        frequencies = np.array(frequencies_cm1) * wavenumber
        # The dielectric variables are functions of frequency
        superstrateDielectric = self.settings['Superstrate dielectric']
        substrateDielectric   = self.settings['Substrate dielectric']
        superstrateDielectricFunction = DielectricFunction(epsType='constant',units='hz',parameters=superstrateDielectric).function()
        substrateDielectricFunction   = DielectricFunction(epsType='constant',units='hz',parameters=substrateDielectric).function()
        # The crystal dielectric has already been defined in the SettingsTab
        # Force a call to get the permittivity, if it needs recalculating it will be
        temp = self.notebook.settingsTab.get_crystal_permittivity(vs_cm1)
        # Make sure the system knows that frequency will be supplied using Hz
        self.notebook.settingsTab.CrystalPermittivity.setUnits('hz')
        # Actually use the permittivity function in what follows
        crystalPermittivityFunction     = self.notebook.settingsTab.CrystalPermittivity.function()
        # Create 3 layers, thickness is converted to metres
        if self.settings['Thickness unit'] == 'nm':
            tom = 1.0E-9
        elif self.settings['Thickness unit'] == 'um':
            tom = 1.0E-6
        elif self.settings['Thickness unit'] == 'mm':
            tom = 1.0E-3
        elif self.settings['Thickness unit'] == 'cm':
            tom = 1.0E-2
        superstrateDepth = tom * self.settings['Superstrate depth']
        superstrateDepth = tom * self.settings['Superstrate depth']
        substrateDepth   = tom * self.settings['Substrate depth']
        crystalDepth     = tom * self.settings['Film thickness']
        # Determine the euler angles
        theta,phi,psi = self.calculate_euler_angles()
        # Set the angle of incidence in radians
        angle = self.settings['Angle of incidence']
        angleOfIncidence      = np.pi / 180.0 * angle
        mode = self.settings['Mode']
        if mode == 'Thick slab' or mode == 'Coherent thin film' or mode == 'Incoherent thin film':
            calculator = self.coherent_calculator
        else:
            calculator = self.partially_incoherent_calculator
        # Call the relevant calculator
        ( self.p_reflectance, 
        self.s_reflectance, 
        self.p_transmittance, 
        self.s_transmittance, 
        self.s_absorbtance, 
        self.p_absorbtance, 
        self.epsilon) = calculator( superstrateDielectricFunction,
                                    substrateDielectricFunction,
                                    crystalPermittivityFunction,
                                    superstrateDepth,
                                    substrateDepth,
                                    crystalDepth,
                                    mode,
                                    theta,
                                    phi,
                                    psi,
                                    angleOfIncidence)
        QCoreApplication.processEvents()
        debugger.print(self.settings['Legend'],'Finished:: calculate - number of frequencies',len(vs_cm1))
        return

    def calculate_euler_angles(self):
        '''Calculate the Euler angles for the crystal to lab transformation'''
        debugger.print(self.settings['Legend'],'Start:: calculate_euler_angles')
        # Get plane specification
        hkl = [ self.settings['Unique direction - h'] , self.settings['Unique direction - k'], self.settings['Unique direction - l'] ]
        sum2 = hkl[0]*hkl[0] + hkl[1]*hkl[1] + hkl[2]*hkl[2]
        if sum2 < 1:
            debugger.print(self.settings['Legend'],'Finished:: calculate_euler_angles')
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
        self.labframe_w.addItem('crystal a-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(a[0],a[1],a[2]) )
        self.labframe_w.addItem('crystal b-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(b[0],b[1],b[2]) )
        self.labframe_w.addItem('crystal c-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(c[0],c[1],c[2]) )
        self.labframe_a = a
        self.labframe_b = b
        self.labframe_c = c
        #a = a / np.linalg.norm(a)
        #b = b / np.linalg.norm(b)
        #c = c / np.linalg.norm(c)
        # print('Projection of a,b,c onto the lab Y-axis (s-pol)', a[1],b[1],c[1])
        # print('Projection of a,b,c onto the lab X-axis (p-pol)', a[0],b[0],c[0])
        debugger.print(self.settings['Legend'],'Finished:: calculate_euler_angles')
        return (theta, phi, psi)

    def get_result(self, vs_cm1, plot_type):
        """Return a particular result"""
        debugger.print(self.settings['Legend'],'Start:: get_result',len(vs_cm1),plot_type)
        self.get_results(vs_cm1)
        debugger.print(self.settings['Legend'],'Finished:: get_result',len(vs_cm1),plot_type)
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
            # print('Error in returning result from CrystalScenarioTab: ',plot_type)
            return None


    def get_results(self, vs_cm1):
        """Return the results of the effective medium theory calculation"""
        debugger.print(self.settings['Legend'],'Start:: get_results',len(vs_cm1),self.refreshRequired)
        if len(vs_cm1)>0 and (self.refreshRequired or len(self.vs_cm1) != len(vs_cm1) or self.vs_cm1[0] != vs_cm1[0] or self.vs_cm1[1] != vs_cm1[1]) :
            debugger.print(self.settings['Legend'],'get_results recalculating')
            self.refresh()
            self.calculate(vs_cm1)
        else:
            debugger.print(self.settings['Legend'],'get_results no need for recalculation')
            self.notebook.progressbars_update(increment=len(vs_cm1))
        debugger.print(self.settings['Legend'],'Finished:: get_results',len(vs_cm1),self.refreshRequired)
        return

    def greyed_out(self):
        """Have a look through the settings and see if we need to grey anything out"""
        # If the single crystal mode is Thick Slab, there is no need for film thickness or substrate permittivity
        debugger.print(self.settings['Legend'],'Start:: greyed_out')
        if self.settings['Mode'] == 'Thick slab':
            self.film_thickness_sb.setEnabled(False)
            self.thickness_unit_cb.setEnabled(False)
            self.substrate_dielectric_sb.setEnabled(False)
            self.partially_incoherent_samples_sb.setEnabled(False)
            self.percentage_partial_incoherence_sb.setEnabled(False)
            self.partially_incoherent_polynomial_sb.setEnabled(False)
            self.partially_incoherent_kernel_sb.setEnabled(False)
        elif self.settings['Mode'] == 'Coherent thin film':
            self.film_thickness_sb.setEnabled(True)
            self.thickness_unit_cb.setEnabled(True)
            self.substrate_dielectric_sb.setEnabled(True)
            self.partially_incoherent_samples_sb.setEnabled(False)
            self.percentage_partial_incoherence_sb.setEnabled(False)
            self.partially_incoherent_polynomial_sb.setEnabled(False)
            self.partially_incoherent_kernel_sb.setEnabled(False)
        elif self.settings['Mode'] == 'Partially incoherent thin film':
            self.film_thickness_sb.setEnabled(True)
            self.thickness_unit_cb.setEnabled(True)
            self.substrate_dielectric_sb.setEnabled(True)
            self.partially_incoherent_samples_sb.setEnabled(True)
            self.percentage_partial_incoherence_sb.setEnabled(True)
            self.partially_incoherent_polynomial_sb.setEnabled(True)
            self.partially_incoherent_kernel_sb.setEnabled(True)
        elif self.settings['Mode'] == 'Incoherent thin film':
            self.film_thickness_sb.setEnabled(True)
            self.thickness_unit_cb.setEnabled(True)
            self.substrate_dielectric_sb.setEnabled(True)
            self.partially_incoherent_samples_sb.setEnabled(False)
            self.percentage_partial_incoherence_sb.setEnabled(False)
            self.partially_incoherent_polynomial_sb.setEnabled(False)
            self.partially_incoherent_kernel_sb.setEnabled(False)
        debugger.print(self.settings['Legend'],'Finished:: greyed_out')

