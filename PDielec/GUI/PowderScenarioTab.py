#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""PowderScenarioTab module."""
import ctypes
import logging
import math
import sys
from functools import partial
from multiprocessing import Array

import numpy as np
from qtpy.QtCore import QCoreApplication, Qt
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from PDielec import Calculator, DielectricFunction, Materials
from PDielec.Constants import amu, angs2bohr, wavenumber
from PDielec.GUI.ScenarioTab import ScenarioTab
from PDielec.Materials import MaterialsDataBase
from PDielec.RamanPolarCalculator import apply_particle_eo_correction
from PDielec.RamanSpectrum import stokes_prefactor, valid_stokes_mode, validate_laser_frequency

logger = logging.getLogger(__name__)
class PowderScenarioTab(ScenarioTab):
    """A class for managing the Powder Infrared Scenario Tab.

    It inherits from :class:`~PDielec.GUI.ScenarioTab`, thus utilizing its layout and properties, 
    with additional features and settings pertinent to powder infrared scenarios.

    Parameters
    ----------
    parent : QWidget
        The parent widget to attach to this scenario tab.
    debug : bool, optional
        If set to True, debug messages will be printed. Defaults to False.

    Attributes
    ----------
    spectroscopy : str
        Defines the type of spectrocopy as 'Powder Infrared'.
    settings : dict
        Dictionary holding various settings related to the scenario.
    methods : list
        List of available effective medium methods.
    shapes : list
        List of available particle shapes.
    material_names : list
        List of available material names.
    materialDefined : bool
        Indicates whether the material is defined manually or not.
    direction : ndarray
        Array defining the unique direction for particle shape.
    depolarisation : ndarray (3x3)
        Array defining the depolarisation factors.
    scenarioIndex : any
        Index associated with the current scenario.
    refresh_required : bool
        Indicates if UI refresh is required.
    calculation_required : bool
        Indicates if recalculation is required due to parameter changes.
    reader : An dft output reader
        Used to read external data.
    realPermittivity : list
        List holding real parts of permittivity values.
    imagPermittivity : list
        List holding imaginary parts of permittivity values.
    absorptionCoefficient : list
        List to hold absorption coefficient values.
    molarAbsorptionCoefficient : list
        List to hold molar absorption coefficient values.
    sp_atr : list
        List to hold ATR specific calculation results.

    Methods
    -------
    initialise_infrared_gui()
        Set the gui appropriate for an infrared calculation
    initialise_atr_gui()
        Set the gui appropriate for an atr calculation
    initialise_raman_gui()
        Set the settings appropriate for an raman calculation
    initialise_infrared_settings()
        Set the settings appropriate for an infrared calculation
    initialise_atr_settings()
        Set the settings appropriate for an atr calculation
    initialise_raman_settings()
        Set the settings appropriate for an raman calculation
    crystal_density()
        Compute and return the density of the crystal based on its volume and mass.
    open_db_button_clicked()
        Handler for the event where the database open button is clicked.
    on_h_sb_changed(value)
        Handler for changes in the 'h' dimension parameter of the crystal's unique direction.
    on_k_sb_changed(value)
        Handler for changes in the 'k' dimension parameter of the crystal's unique direction.
    on_l_sb_changed(value)
        Handler for changes in the 'l' dimension parameter of the crystal's unique direction.
    on_shape_cb_activated(index)
        Handler for changes in the selected shape of powder particles.
    on_methods_cb_activated(index)
        Handler to deal with changes in the selected calculation method for effective medium theory.
    on_mf_sb_changed(value)
        Handler for changes to the mass fraction value, updating the volume fraction value accordingly.
    update_vf_sb()
        Update the volume fraction spinbox according to the mass fraction and recalculates as necessary.
    on_aoverb_sb_changed(value)
        Handler for changes to the ellipsoid a/b ratio or eccentricity.
    on_legend_le_changed(text)
        Handler for changes to the legend label of the scenario.
    on_sigma_sb_changed(value)
        Handler for changes to the particle size distribution sigma.
    on_size_sb_changed(value)
        Handler for changes to the particle size.
    on_vf_sb_changed(value)
        Handle changes to the volume fraction spinbox, updating mass fraction accordingly.
    update_mf_sb()
        Update the mass fraction spinbox according to the volume fraction and recalculates as necessary.
    on_matrix_cb_activated(index)
        Handler for changes to the selected support matrix in the ComboBox.
    on_density_sb_changed(value)
        Handle changes to the density spinbox of the support matrix.
    on_bubble_vf_sb_changed(value)
        Handle changes to the bubble volume fraction spinbox.
    on_bubble_radius_sb_changed(value)
        Handle changes to the bubble radius spinbox.
    on_permittivity_i_sb_changed(value)
        Handle changes to the imaginary part of the support matrix permittivity.
    on_permittivity_r_sb_changed(value)
        Handle changes to the real part of the support matrix permittivity.
    on_atr_index_sb_changed(value)
        Handle changes to the ATR material refractive index spinbox.
    on_atr_incident_ang_sb_changed(value)
        Handle changes to the ATR incident angle spinbox.
    on_atr_spolfrac_sb_changed(value)
        Handle changes to the ATR s-polarisation fraction spinbox.
    on_laser_frequency_sb_changed(value)
        Handle changes to the laser frequency spin box
    on_polarisation_cb_activated(index)
        Handle changes to polarisation setting for the experiment
    on_temperature_sb_changed(value)
        Handle changes to temperature spin box
    change_greyed_out()
        Enable or disable elements of the GUI based on the current scenario settings.
    calculate(vs_cm1)
        Calls the _calculate_?? routines based on the spectroscopy 
    _calculate_atr(vs_cm1)
        Calculate the ATR spectrum
    _calculate_infrared(vs_cm1)
        Calculate the infrared powder absorption for the given frequency range.
    _calculate_raman(vs_cm1)
        Calculate the powder Raman scattering for the given frequency range.
    calculate_depolarisation_tensor()
        Calculate the depolarisation tensor
    get_result(vs_cm1, plot_type)
        Return the calculation result for the specified plot type.
    get_results(vs_cm1)
        Prepare the necessary results for displaying or processing.
    refresh(force=False)
        Refresh the GUI interface with up to date values.
    refresh_raman()
        Refresh the raman gui components
    refresh_infrared()
        Refresh the infrared gui components
    refresh_atr()
        Refresh the atr gui components

    """

    def __init__(self, parent, spectroscopy="Powder Infrared", debug=False):
        """Initialize the ScenarioTab subclass for Powder Infrared Scenario with UI and connectivity.

        Parameters
        ----------
        parent : QWidget
            The parent widget to attach to this scenario tab.
        spectroscopy : str
            Can be one of "Powder Infrared", "Powder ATR", "Powder Raman"
        debug : bool, optional
            If set to True, debug messages will be printed. Defaults to False.

        Attributes
        ----------
        spectrocopy : str
            Defines the type of spectrocopy as 'Powder Infrared', "Powder ATR", "Powder Raman".
        settings : dict
            Dictionary holding various settings related to the scenario.
        methods : list
            List of available effective medium methods.
        shapes : list
            List of available particle shapes.
        material_names : list
            List of available material names.
        materialDefined : bool
            Indicates whether the material is defined manually or not.
        direction : ndarray
            Array defining the unique direction for particle shape.
        depolarisation : ndarray (3x3)
            Array defining the depolarisation factors.
        scenarioIndex : any
            Index associated with the current scenario.
        refresh_required : bool
            Indicates if UI refresh is required.
        calculation_required : bool
            Indicates if recalculation is required due to parameter changes.
        reader : An dft output reader
            Used to read external data.
        realPermittivity : list
            List holding real parts of permittivity values.
        imagPermittivity : list
            List holding imaginary parts of permittivity values.
        absorptionCoefficient : list
            List to hold absorption coefficient values.
        molarAbsorptionCoefficient : list
            List to hold molar absorption coefficient values.
        sp_atr : list
            List to hold ATR specific calculation results.
        raman_spectrum : list
            List to hold calculate Raman spectrum

        """        
        ScenarioTab.__init__(self,parent)
        logger.debug("Start:: initialiser")
        if "Powder" not in spectroscopy:
            logger.error(f"PowderSenarioTab failed incompatible spectroscopy {spectroscopy}")
        self.spectroscopy = spectroscopy
        #
        # No testing on spectroscopy type, all settings are set independent of the spectroscopy
        #
        self.initialise_infrared_settings()
        self.initialise_atr_settings()
        self.initialise_raman_settings()
        #
        # Create a scenario tab
        #
        vbox = QVBoxLayout()
        form = QFormLayout()
        if "Powder Infrared" in self.spectroscopy:
            self.initialise_infrared_gui(vbox,form)
        elif "Powder ATR" in self.spectroscopy:
            self.initialise_infrared_gui(vbox,form)
            self.initialise_atr_gui(vbox,form)
        elif "Powder Raman" in self.spectroscopy:
            self.initialise_infrared_gui(vbox,form, include_infrared_only=False)
            self.initialise_raman_gui(vbox,form)
        #
        # Add a legend option
        #
        self.legend_le = QLineEdit(self)
        self.legend_le.setToolTip("The legend will be used to describe the results in the plot")
        self.legend_le.setText(self.settings["Legend"])
        self.legend_le.textChanged.connect(self.on_legend_le_changed)
        label = QLabel("Powder scenario legend",self)
        label.setToolTip("The legend will be used to describe the results in the plotting tab")
        form.addRow(label, self.legend_le)
        #
        # Final buttons
        #
        hbox = self.add_scenario_buttons()
        form.addRow(hbox)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        # sort out greying of boxes
        self.change_greyed_out()
        logger.debug("Finished:: initialiser")
        return

    def initialise_atr_gui(self, vbox, form):
        """Initialise the GUI for atr calculations.

        Parameters
        ----------
        vbox : QVBoxLayout
            vbox is defined by the PowderScenarioTab initialised
        form : QFormLayout
            form is defined by the PowderScenarioTab initialised

        Returns
        -------
        vbox, form

        """
        #
        # Refractive Index
        #
        self.atr_index_sb = QDoubleSpinBox(self)
        self.atr_index_sb.setRange(0.001, 100.0)
        self.atr_index_sb.setSingleStep(0.01)
        self.atr_index_sb.setDecimals(3)
        self.atr_index_sb.setToolTip("Define the ATR material refractive index")
        self.atr_index_sb.setValue(self.settings["ATR material refractive index"])
        self.atr_index_sb.valueChanged.connect(self.on_atr_index_sb_changed)
        label = QLabel("ATR material refractive index", self)
        label.setToolTip("Define the ATR material refractive index")
        form.addRow(label, self.atr_index_sb)
        #
        # Incident angle in degreees
        #
        self.atr_incident_ang_sb = QDoubleSpinBox(self)
        self.atr_incident_ang_sb.setRange(0.0, 180.0)
        self.atr_incident_ang_sb.setSingleStep(0.1)
        self.atr_incident_ang_sb.setDecimals(1)
        self.atr_incident_ang_sb.setToolTip("Define the ATR incident angle")
        self.atr_incident_ang_sb.setValue(self.settings["ATR theta"])
        self.atr_incident_ang_sb.valueChanged.connect(self.on_atr_incident_ang_sb_changed)
        label = QLabel("ATR incident angle", self)
        label.setToolTip("Define the ATR incident angle")
        form.addRow(label, self.atr_incident_ang_sb)
        #
        # S polarisation fraction
        #
        self.atr_spolfrac_sb = QDoubleSpinBox(self)
        self.atr_spolfrac_sb.setRange(0.0, 1.0)
        self.atr_spolfrac_sb.setSingleStep(0.01)
        self.atr_spolfrac_sb.setDecimals(3)
        self.atr_spolfrac_sb.setToolTip("Define the ATR S polarisation fraction, the rest is P polarisation")
        self.atr_spolfrac_sb.setValue(self.settings["ATR S polarisation fraction"])
        self.atr_spolfrac_sb.valueChanged.connect(self.on_atr_spolfrac_sb_changed)
        label = QLabel("ATR S polarisation fraction", self)
        label.setToolTip("Define the S polarisation fraction, the rest is P polarisation")
        form.addRow(label, self.atr_spolfrac_sb)
        return vbox, form

    def initialise_infrared_gui(self, vbox, form, include_infrared_only=True):
        """Initialise the GUI for atr calculations.

        Parameters
        ----------
        vbox : QVBoxLayout
            vbox is defined by the PowderScenarioTab initialised
        form : QFormLayout
            form is defined by the PowderScenarioTab initialised
        include_infrared_only : bool
            If True, include controls used only by the infrared/ATR effective-medium calculation.

        Returns
        -------
        vbox, form

        """
        #
        # Option to open a database of permittivities for the support
        # label and button are defined the parent class
        #
        form.addRow(self.openDB_label, self.openDB_button)
        #
        # Add a name for the database
        # label and button are defined the parent class
        form.addRow(self.database_le_label, self.database_le)
        #
        # Support matrix, read information from the database
        #
        hbox = QHBoxLayout()
        self.matrix_cb = QComboBox(self)
        self.matrix_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.matrix_cb.setToolTip("Define the permittivity and density of the support matrix")
        self.material_names = self.DataBase.get_sheet_names()
        self.matrix_cb.addItems(self.material_names)
        if self.settings["Matrix"] not in self.material_names:
            self.settings["Matrix"] = self.material_names[0]
        index = self.matrix_cb.findText(self.settings["Matrix"], Qt.MatchFixedString)
        if index >=0:
            self.matrix_cb.setCurrentIndex(index)
        self.matrix_cb.activated.connect(self.on_matrix_cb_activated)
        if "Material defined manually" in self.settings["Matrix"]:
            self.materialDefinedManually = True
            self.matrixMaterial = Materials.Constant("manual",permittivity=self.settings["Matrix permittivity"],density=self.settings["Matrix density"])
        else:
            self.materialDefinedManually = False
            self.matrixMaterial = self.DataBase.get_material(self.settings["Matrix"])
        label = QLabel("Support matrix",self)
        label.setToolTip("Define the permittivity and density of the support matrix")
        self.matrix_info_le = QLineEdit(self)
        self.matrix_info_le.setToolTip("Provides details about database entry")
        text = self.matrixMaterial.get_information()
        self.matrix_info_le.setText(text)
        self.matrix_info_le.setReadOnly(True)
        hbox.addWidget(self.matrix_cb)
        hbox.addWidget(self.matrix_info_le)
        form.addRow(label, hbox)
        #
        # Set the Matrix density and displayed permittivity
        #
        self.settings["Matrix density"] = self.matrixMaterial.get_density()
        self.matrixPermittivityFunction = self.matrixMaterial.get_permittivity_function()
        self.settings["Matrix permittivity"] = self._matrix_display_permittivity()
        #
        # Support matrix permittivity
        #
        self.density_sb = QDoubleSpinBox(self)
        self.density_sb.setRange(0.001, 100.0)
        self.density_sb.setSingleStep(0.01)
        self.density_sb.setDecimals(3)
        self.density_sb.setToolTip("Define the support matrix density. \nThis makes changes to the support density and permittivity")
        self.density_sb.setValue(self.settings["Matrix density"])
        self.density_sb.valueChanged.connect(self.on_density_sb_changed)
        label = QLabel("Support density", self)
        label.setToolTip("Define the support matrix density. \nThis makes changes to the support density and permittivity")
        form.addRow(label, self.density_sb)
        #
        # Support matrix permittivity
        #
        hbox = QHBoxLayout()
        self.permittivity_r_sb = QDoubleSpinBox(self)
        self.permittivity_r_sb.setRange(0.0, 1000.0)
        self.permittivity_r_sb.setSingleStep(0.01)
        self.permittivity_r_sb.setDecimals(3)
        self.permittivity_r_sb.setToolTip("Define the real component of the support matrix permittivity")
        self.permittivity_r_sb.setValue(np.real(self.settings["Matrix permittivity"]))
        self.permittivity_r_sb.valueChanged.connect(self.on_permittivity_r_sb_changed)
        hbox.addWidget(self.permittivity_r_sb)
        self.permittivity_i_sb = QDoubleSpinBox(self)
        self.permittivity_i_sb.setRange(0.0, 1000.0)
        self.permittivity_i_sb.setSingleStep(0.01)
        self.permittivity_i_sb.setDecimals(3)
        self.permittivity_i_sb.setToolTip("Define imaginary component of the the support matrix permittivity")
        self.permittivity_i_sb.setValue(np.imag(self.settings["Matrix permittivity"]))
        self.permittivity_i_sb.valueChanged.connect(self.on_permittivity_i_sb_changed)
        hbox.addWidget(self.permittivity_i_sb)
        label = QLabel("Support permittivity", self)
        label.setToolTip("Define the complex support matrix permittivity")
        form.addRow(label, hbox)
        #
        # Bubble volume fraction
        #
        if include_infrared_only:
            self.bubble_vf_sb = QDoubleSpinBox(self)
            self.bubble_vf_sb.setRange(0.0, 100.0*(1.0-self.settings["Volume fraction"]))
            self.bubble_vf_sb.setSingleStep(1.0)
            self.bubble_vf_sb.setDecimals(1)
            self.bubble_vf_sb.setToolTip("Define the % volume fraction of air bubble inclusions in the matrix")
            self.bubble_vf_sb.setValue(100*self.settings["Bubble volume fraction"])
            self.bubble_vf_sb.valueChanged.connect(self.on_bubble_vf_sb_changed)
            label = QLabel("% Air void volume fraction", self)
            label.setToolTip("Define the % volume fraction of air bubble inclusions in the matrix")
            form.addRow(label, self.bubble_vf_sb)
            #
            # Bubble radius in microns
            #
            self.bubble_radius_sb = QDoubleSpinBox(self)
            self.bubble_radius_sb.setRange(0.001, 1000.0)
            self.bubble_radius_sb.setSingleStep(1.0)
            self.bubble_radius_sb.setDecimals(3)
            self.bubble_radius_sb.setToolTip("Define the air bubble radius")
            self.bubble_radius_sb.setValue(self.settings["Bubble radius"])
            self.bubble_radius_sb.valueChanged.connect(self.on_bubble_radius_sb_changed)
            label = QLabel("Air void radius (μm)", self)
            label.setToolTip("Define the air void radius")
            form.addRow(label, self.bubble_radius_sb)
        #
        # Mass fraction of dielectric medium
        #
        self.mf_sb = QDoubleSpinBox(self)
        self.mf_sb.setRange(0.000001, 100.0)
        self.mf_sb.setSingleStep(0.1)
        self.mf_sb.setDecimals(6)
        self.mf_sb.setToolTip("The percentage mass fraction of the dielectric medium. \nNote that volume and mass fraction are linked")
        self.mf_sb.setValue(100.0*self.settings["Mass fraction"])
        self.mf_sb.valueChanged.connect(self.on_mf_sb_changed)
        label = QLabel("% Mass fraction of dielectric", self)
        label.setToolTip("The percentage mass fraction of the dielectric medium. \nNote that volume and mass fraction are linked")
        form.addRow(label, self.mf_sb)
        #
        # Volume fraction of dielectric medium
        #
        self.vf_sb = QDoubleSpinBox(self)
        self.vf_sb.setRange(0.000001, 100.0*(1.0-self.settings["Bubble volume fraction"]))
        self.vf_sb.setSingleStep(0.1)
        self.vf_sb.setDecimals(6)
        self.vf_sb.setToolTip("The percentage volume fraction of the dielectric medium. \nNote that volume and mass fraction are linked")
        self.vf_sb.valueChanged.connect(self.on_vf_sb_changed)
        self.vf_sb.setValue(100.0*self.settings["Volume fraction"])
        label = QLabel("% Volume fraction of dielectric", self)
        label.setToolTip("The percentage volume fraction of the dielectric medium. \nNote that volume and mass fraction are linked")
        form.addRow(label, self.vf_sb)
        #
        # Calculation method
        #
        if include_infrared_only:
            self.methods_cb = QComboBox(self)
            self.methods_cb.setToolTip("Choose the calculation method for the effective medium theory")
            self.methods_cb.addItems(self.methods)
            index = self.methods_cb.findText(self.settings["Effective medium method"], Qt.MatchFixedString)
            if index >=0:
                self.methods_cb.setCurrentIndex(index)
            else:
                logger.error(f"Method index was not 0 {self.settings['Effective medium method']}")
            self.methods_cb.activated.connect(self.on_methods_cb_activated)
            label = QLabel("Method",self)
            label.setToolTip("Choose the calculation method for the effective medium theory")
            form.addRow(label, self.methods_cb)
            #
            # Particle size option
            #
            self.size_sb = QDoubleSpinBox(self)
            self.size_sb.setRange(0.000001, 1000.0)
            self.size_sb.setSingleStep(0.1)
            self.size_sb.setDecimals(6)
            self.size_sb.setToolTip("Define the particle radius of the sphere in μm.")
            self.size_sb.setValue(self.settings["Particle size(mu)"])
            self.size_sb.valueChanged.connect(self.on_size_sb_changed)
            label = QLabel("Particle radius (μm)",self)
            label.setToolTip("Define the particle radius of the sphere in μm.")
            form.addRow(label, self.size_sb)
            #
            # Particle sigma option
            #
            self.sigma_sb = QDoubleSpinBox(self)
            self.sigma_sb.setRange(0.0, 1000.0)
            self.sigma_sb.setSingleStep(0.1)
            self.sigma_sb.setDecimals(6)
            self.sigma_sb.setToolTip("Define the particle size distribution as a lognormal distribution with the given sigma. \nOnly applicable for the Mie method")
            self.sigma_sb.setValue(self.settings["Particle size distribution sigma(mu)"])
            self.sigma_sb.valueChanged.connect(self.on_sigma_sb_changed)
            label = QLabel("Particle sigma (μm)",self)
            label.setToolTip("Define the particle size distribition as a lognormal with the given sigma. \nOnly applicable for the Mie method")
            form.addRow(label, self.sigma_sb)
        #
        # Crystallite shape
        #
        self.shape_cb = QComboBox(self)
        self.shape_cb.setToolTip("Choose a particle shape. \nFor the Mie methods only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [abc].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction")
        self.shape_cb.addItems(self.shapes)
        index = self.shape_cb.findText(self.settings["Particle shape"], Qt.MatchFixedString)
        if index >=0:
            self.shape_cb.setCurrentIndex(index)
        else:
            logger.error(f"Method index was not 0 {self.settings['Particle shape']}")
        self.shape_cb.activated.connect(self.on_shape_cb_activated)
        label = QLabel("Particle shape",self)
        label.setToolTip("Choose a particle shape. \nFor the Mie methods only sphere is allowed.  \nFor shapes other than sphere there is a unique direction. \nFor ellipsoidal and needle like this is a direction [abc].  \nFor a plate the perpendicular to a crystal face (hkl) is used to define the unique direction")
        form.addRow(label, self.shape_cb)
        #
        # Particle shape information
        # unique direction (hkl) or [abc]
        #
        self.h_sb = QSpinBox(self)
        self.h_sb.setToolTip("Define the h dimension of the unique direction")
        self.h_sb.setRange(-20,20)
        self.h_sb.setValue(self.settings["Unique direction - h"])
        self.h_sb.valueChanged.connect(self.on_h_sb_changed)
        self.k_sb = QSpinBox(self)
        self.k_sb.setToolTip("Define the k dimension of the unique direction")
        self.k_sb.setRange(-20,20)
        self.k_sb.setValue(self.settings["Unique direction - k"])
        self.k_sb.valueChanged.connect(self.on_k_sb_changed)
        self.l_sb = QSpinBox(self)
        self.l_sb.setToolTip("Define the l dimension of the unique direction")
        self.l_sb.setRange(-20,20)
        self.l_sb.setValue(self.settings["Unique direction - l"])
        self.l_sb.valueChanged.connect(self.on_l_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.h_sb)
        hbox.addWidget(self.k_sb)
        hbox.addWidget(self.l_sb)
        self.hkl_label = QLabel("Unique direction [abc]",self)
        self.hkl_label.setToolTip("Define the unique direction by [abc] or (hkl). \n[abc] is used by needles and ellipsoids.  It defines the unique direction in crystallographic units. \n(hkl) is used by plates it defines a surface and the unique direction is perpendicular to it.")

        form.addRow(self.hkl_label, hbox)
        #
        # a over b ratio for ellipse
        #
        self.aoverb_sb = QDoubleSpinBox(self)
        self.aoverb_sb.setRange(0.0, 1000.0)
        self.aoverb_sb.setSingleStep(0.1)
        self.aoverb_sb.setDecimals(6)
        self.aoverb_sb.setToolTip("Define the ellipsoid a/b ratio or eccentricity.  \nOnly applicable for the ellipsoid shapes \na/b < 1: oblate ellipsoid \na/b > 1: prolate ellipsoid")
        self.aoverb_sb.setValue(self.settings["Ellipsoid a/b"])
        self.aoverb_sb.valueChanged.connect(self.on_aoverb_sb_changed)
        label = QLabel("Ellipsoid a/b eccentricty",self)
        label.setToolTip("Define the ellipsoid a/b ratio or eccentricity.  \nOnly applicable for the ellipsoid shapes \na/b < 1: oblate ellipsoid \na/b > 1: prolate ellipsoid")
        form.addRow(label, self.aoverb_sb)
        return vbox, form


    def initialise_atr_settings(self):
        """Set the settings attribute for atr calculations.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.settings["ATR material refractive index"] = 4.0
        self.settings["ATR theta"] = 45.0
        self.settings["ATR S polarisation fraction"] = 0.5
        self.sp_atr = []
        return

    def initialise_infrared_settings(self):
        """Set the settings attribute for infrared calculations.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.no_calculations_required = 1
        self.settings["Matrix"] = "ptfe"
        self.settings["Matrix density"] = None
        self.settings["Matrix permittivity"] = None
        self.settings["Bubble radius"] = 30.0
        self.settings["Bubble volume fraction"] = 0.0
        self.settings["Mass fraction"] = 0.1
        self.settings["Volume fraction"] = 0.1
        self.settings["Particle size(mu)"] = 0.0001
        self.settings["Particle size distribution sigma(mu)"] = 0.0
        self.settings["Ellipsoid a/b"] = 1.0
        self.settings["Unique direction - h"] = 0
        self.settings["Unique direction - k"] = 0
        self.settings["Unique direction - l"] = 1
        self.settings["Mass or volume fraction"] = "volume"
        self.settings["Effective medium method"] = "Maxwell-Garnett"
        self.settings["Particle shape"] = "Sphere"
        self.methods = ["Maxwell-Garnett", "Bruggeman", "Averaged Permittivity", "Mie"]
        self.shapes = ["Sphere", "Needle", "Plate", "Ellipsoid"]
        self.material_names = []
        self.materialDefinedManually = False
        self.direction = np.array([0,0,0])
        self.depolarisation = np.array([0,0,0])
        self.scenarioIndex = None
        self.refresh_required = True
        self.calculation_required = False
        self.reader = None
        self.realPermittivity = []
        self.imagPermittivity = []
        self.absorptionCoefficient = []
        self.molarAbsorptionCoefficient = []
        return

    def crystal_density(self):
        """Find the crystal density from the current reader and return the density.

        Returns
        -------
        float
            The crystal density.

        """
        if not self.reader:
            logger.debug("Finished:: crystal_density - no reader")
            return 1.0
        return self.reader.get_crystal_density()

    def _matrix_display_permittivity(self):
        """Return the support-matrix permittivity shown in the Powder GUI.

        Powder Raman uses the matrix optical permittivity for particle local-field
        corrections.  Powder IR/ATR uses the matrix dielectric function at each
        IR frequency, so the editable support-matrix value should show the
        zero-frequency dielectric value instead.
        """
        if self.spectroscopy == "Powder Raman":
            return self.matrixMaterial.get_optical_permittivity()
        return self.matrixPermittivityFunction(0.0)

    def open_db_button_clicked(self):
        """Open a new materials' database.

        Returns
        -------
        None

        """
        logger.debug("Start:: open_db_button_clicked")
        self.open_data_base()
        if self.settings["Matrix"] not in self.material_names:
            self.settings["Matrix"] = self.material_names[0]
        self.matrixMaterial = self.DataBase.get_material(self.settings["Matrix"])
        # Check to see that the matrix return a scalar permittivity
        if self.matrixMaterial.is_tensor():
            logger.error("Error: matrix must have a scalar permittivity using ptfe")
            self.settings["Matrix"] = "ptfe"
            self.matrixMaterial = self.DataBase.get_material(self.settings["Matrix"])
        materialPermittivityFunction = self.matrixMaterial.get_permittivity_function()
        self.settings["Matrix permittivity"] = materialPermittivityFunction(0.0)
        self.settings["Matrix density"] = self.matrixMaterial.get_density()
        self.refresh_required = True
        self.refresh()
        return

    def on_h_sb_changed(self,value):
        """Handle a change to the h parameter of the (hkl) surface.

        Parameters
        ----------
        value : int
            The h, k or l value

        """
        logger.debug(f"{self.settings['Legend']} on_h_sb_changed {value}")
        self.refresh_required = True
        self.settings["Unique direction - h"] = value
        return

    def on_k_sb_changed(self,value):
        """Handle a change to the k parameter of the (hkl) surface.

        Parameters
        ----------
        value : int
            The h, k or l value

        """
        logger.debug(f"{self.settings['Legend']} on_k_sb_changed {value}")
        self.refresh_required = True
        self.settings["Unique direction - k"] = value
        return

    def on_l_sb_changed(self,value):
        """Handle a change to the l parameter of the (hkl) surface.

        Parameters
        ----------
        value : int
            The h, k or l value

        """
        logger.debug(f"{self.settings['Legend']} on_l_sb_changed {value}")
        self.refresh_required = True
        self.settings["Unique direction - l"] = value
        return

    def on_shape_cb_activated(self,index):
        """Handle changes to the shape of the particles.

        Parameters
        ----------
        index
            The index in the list of shapes

        """
        logger.debug(f"{self.settings['Legend']} on shape cb activated {index}")
        self.refresh_required = True
        self.settings["Particle shape"] = self.shapes[index]
        if self.settings["Particle shape"] == "Sphere":
            self.settings["Unique direction - h"] = 0
            self.settings["Unique direction - k"] = 0
            self.settings["Unique direction - l"] = 0
        self.change_greyed_out()
        return

    def on_methods_cb_activated(self,index):
        """Handle changes in the calculation method for the effective medium theory.

        Parameters
        ----------
        index
            The index in the list of methods

        """
        logger.debug(f"{self.settings['Legend']} on methods cb activated {index}")
        self.refresh_required = True
        self.settings["Effective medium method"] = self.methods[index]
        if self.settings["Effective medium method"] == "Mie" or self.settings["Effective medium method"] == "Anisotropic-Mie":
            self.settings["Particle shape"] = "Sphere"
        elif self.settings["Effective medium method"] == "Maxwell-Garnett" or self.settings["Effective medium method"] == "Bruggeman":
            self.settings["Particle size distribution sigma(mu)"] = 0.0
        elif self.settings["Effective medium method"] == "Averaged Permittivity":
            self.settings["Particle size(mu)"] = 0.0001
            self.settings["Particle size distribution sigma(mu)"] = 0.0
        self.change_greyed_out()
        return

    def on_mf_sb_changed(self,value):
        """Handle a mass fraction change and update the volume fraction.

        The mass fraction is prioritized and stored.  The volume fraction is calculated.

        Parameters
        ----------
        value : float
            The mass fraction

        """
        logger.debug(f"{self.settings['Legend']} on mass fraction line edit changed {value}")
        self.refresh_required = True
        self.settings["Mass or volume fraction"] = "mass"
        self.settings["Mass fraction"] =  value/100.0
        self.update_vf_sb()
        return

    def update_vf_sb(self):
        """Update the volume fraction according to the mass fraction.

        Parameters
        ----------
        None

        """
        mf1 = self.settings["Mass fraction"]
        mf2 = 1.0 - mf1
        rho1 = self.crystal_density()
        rho2 = self.settings["Matrix density"]
        #
        # Avoid overflow through division by 0
        # When matrix density is zero (e.g. 'none' matrix), volume fraction is undefined;
        # leave the existing volume fraction unchanged.
        #
        mf2  = max(mf2,1.0E-18)
        rho1 = max(rho1,1.0E-18)
        if rho2 < 1.0E-18:
            vf1 = self.settings["Volume fraction"]
        else:
            vf1 = ( 1.0 - self.settings["Bubble volume fraction"] ) * (mf1/mf2)*(rho2/rho1) / ( 1 + (mf1/mf2)*(rho2/rho1))
        self.settings["Volume fraction"] = vf1
        blocking_state = self.vf_sb.signalsBlocked()
        self.vf_sb.blockSignals(True)
        self.vf_sb.setValue(100.0*vf1)
        if hasattr(self, "bubble_vf_sb"):
            self.bubble_vf_sb.setRange(0.0, 100.0*(1.0-self.settings["Volume fraction"]))
        self.vf_sb.setRange(0.0, 100.0*(1.0-self.settings["Bubble volume fraction"]))
        self.vf_sb.blockSignals(blocking_state)
        logger.debug(f"{self.settings['Legend']} Update_vf_sb")
        logger.debug(f"{self.settings['Legend']} rho 1 {rho1}")
        logger.debug(f"{self.settings['Legend']} rho 2 {rho2}")
        logger.debug(f"{self.settings['Legend']} vf 1 {vf1}")
        return

    def on_aoverb_sb_changed(self,value):
        """Handle a change to the a/b ratio for an ellipsoid.

        Parameters
        ----------
        value : float
            The a/b ratio of the elllipsoid

        """
        logger.debug(f"{self.settings['Legend']} on_aoverb_le_changed {value}")
        self.refresh_required = True
        self.settings["Ellipsoid a/b"] = value
        return

    def on_legend_le_changed(self,text):
        """Handle a legend change.

        Parameters
        ----------
        text : str
            The legend

        """
        logger.debug(f"{self.settings['Legend']} on legend change {text}")
        self.refresh_required = True
        self.settings["Legend"] = text
        return

    def on_sigma_sb_changed(self,value):
        """Handle a particle size distribution change.

        Parameters
        ----------
        value : float
            The particle size distribution width in microns

        """
        logger.debug(f"{self.settings['Legend']} on sigma line edit changed {value}")
        self.refresh_required = True
        self.settings["Particle size distribution sigma(mu)"] = value
        return

    def on_size_sb_changed(self,value):
        """Handle a particle size change.

        Parameters
        ----------
        value : float
            The particle size in microns

        """
        logger.debug(f"{self.settings['Legend']} on size line edit changed {value}")
        self.refresh_required = True
        self.settings["Particle size(mu)"] = value
        return

    def on_vf_sb_changed(self,value):
        """Handle a volume fraction change, alter the mass fraction accordingly.

        Parameters
        ----------
        value : float
            The volume fraction

        """
        logger.debug(f"{self.settings['Legend']} Start:: on_vf_sb_changed {value}")
        self.refresh_required = True
        self.settings["Mass or volume fraction"] = "volume"
        self.settings["Volume fraction"] = value/100.0
        self.update_mf_sb()
        logger.debug(f"{self.settings['Legend']} Finished:: on_vf_sb_changed {value}")
        return

    def update_mf_sb(self):
        """Update the mass fraction and according to the volume fraction.

        Parameters
        ----------
        None

        """
        logger.debug(f"{self.settings['Legend']} Start:: update_mf_sb")
        vf1 = self.settings["Volume fraction"]
        vf2 = 1.0 - vf1 - self.settings["Bubble volume fraction"]
        rho1 = self.crystal_density()
        rho2 = self.settings["Matrix density"]
        denominator = rho1*vf1 + rho2*vf2
        mf1 = 1.0 if denominator < 1e-30 else rho1 * vf1 / denominator
        self.settings["Mass fraction"] = mf1
        blocking_state = self.mf_sb.signalsBlocked()
        self.mf_sb.blockSignals(True)
        self.mf_sb.setValue(100.0*mf1)
        self.mf_sb.blockSignals(blocking_state)
        logger.debug(f"{self.settings['Legend']} Update_mf_sb")
        logger.debug(f"{self.settings['Legend']} rho 1 {rho1}")
        logger.debug(f"{self.settings['Legend']} rho 2 {rho2}")
        logger.debug(f"{self.settings['Legend']} mf 1 {mf1}")
        logger.debug(f"{self.settings['Legend']} Finished:: update_mf_sb")
        return

    def on_matrix_cb_activated(self,index):
        """Handle a change to the support matrix supplied by the materials' database.

        A refresh is flagged and the matrix name set from the combobox.
        QT signals are blocked while data is being updated.
        A material object is defined on the basis of the matrix name and the database
        A refresh is performed

        Parameters
        ----------
        index : int
            The index is not used but the associated text of the combobox is.

        Returns
        -------
        None

        """
        logger.debug(f"{self.settings['Legend']} on matrix combobox activated {index}")
        logger.debug(f"{self.settings['Legend']} on matrix combobox activated {self.matrix_cb.currentText()}")
        # We will need to recalculate everything for a new support matrix
        self.refresh_required = True
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
        if "Material defined manually" in matrix:
            # The manual option has been chosen, so create a new material with the right permittivity and density
            self.materialDefinedManually = True
            self.matrixMaterial = Materials.Constant("manual",permittivity=self.settings["Matrix permittivity"],density=self.settings["Matrix density"])
            # Store the new matrix material name
            self.settings["Matrix"] = matrix
        else:
            self.materialDefinedManually = False
            if "Material defined manually" in self.material_names:
                # We don't need the manual entry any more
                self.matrix_cb.clear()
                self.material_names = self.material_names[:-1]
                self.matrix_cb.addItems(self.material_names)
            # Read the material information for permittivity and density from the data base
            matrixMaterial = self.DataBase.get_material(matrix)
            if matrixMaterial.is_scalar():
                # Only change the matrix material if it is a scalar material
                self.matrixMaterial = matrixMaterial
                # Store the new matrix material name
                self.settings["Matrix"] = matrix
            else:
                logger.error("Error: matrix material must have a scalar permittivity")
        # Use the spectroscopy-appropriate displayed permittivity of the matrix material
        self.matrixPermittivityFunction = self.matrixMaterial.get_permittivity_function()
        self.settings["Matrix permittivity"] = self._matrix_display_permittivity()
        self.settings["Matrix density"] = self.matrixMaterial.get_density()
        self.density_sb.setValue(self.settings["Matrix density"])
        # Update the matrix material information
        text = self.matrixMaterial.get_information()
        self.matrix_info_le.setText(text)
        # Update the values of the real and imaginary permittivity
        self.permittivity_r_sb.setValue(np.real(self.settings["Matrix permittivity"]))
        self.permittivity_i_sb.setValue(np.imag(self.settings["Matrix permittivity"]))
        # This is a new material so reset the volume / mass fractions
        # volume fraction takes precedence
        if self.settings["Mass or volume fraction"] == "volume":
            self.settings["Volume fraction"] = 0.1
            self.update_mf_sb()
            self.update_vf_sb()
        else:
            self.settings["Mass fraction"] = 0.1
            self.update_vf_sb()
            self.update_mf_sb()
        # Restore the signal settings on the widgets
        self.matrix_cb.blockSignals(m_blocking)
        self.density_sb.blockSignals(d_blocking)
        self.permittivity_r_sb.blockSignals(r_blocking)
        self.permittivity_i_sb.blockSignals(i_blocking)
        self.refresh()
        self.refresh_required = True
        return

    def on_density_sb_changed(self,value):
        """Handle a change to the matrix density.

        Changes in density cause changes in the volume and mass fractions.  
        These are calculated according to the mass or volume fraction having priority and the spinboxed updated.
        A refresh is forced

        Parameters
        ----------
        value : float
            The density

        """
        self.settings["Matrix density"] = value
        # update the matrix density
        self.matrixMaterial.set_density(value)
        # Force the matrix to be defined manually
        self.settings["Matrix"] = "Material defined manually"
        self.materialDefinedManually = True
        # volume fraction taked precedence
        if self.settings["Mass or volume fraction"] == "volume":
            self.update_mf_sb()
            self.update_vf_sb()
        else:
            self.update_vf_sb()
            self.update_mf_sb()
        logger.debug(f"{self.settings['Legend']} on density line edit changed {value}")
        self.refresh_required = True
        self.refresh()
        self.refresh_required = True
        return

    def on_bubble_vf_sb_changed(self,value):
        """Handle a change to the bubble volume fraction.

        A change in the bubble volume fraction causes changes in the volume and mass fractions
        These are calculated according to the mass or volume fraction having priority and the spinboxed updated.
        A refresh is forced

        Parameters
        ----------
        value : float
            The bubble volume fraction

        """
        self.settings["Bubble volume fraction"] = value/100.0
        if self.settings["Mass or volume fraction"] == "volume":
            self.update_mf_sb()
        else:
            self.update_vf_sb()
        logger.debug(f"{self.settings['Legend']} on bubble volume fraction changed {value}")
        self.refresh_required = True
        return

    def on_bubble_radius_sb_changed(self,value):
        """Update the bubble radius setting and mark refresh as required.

        This method updates the 'Bubble radius' in the settings dictionary and marks the instance for refresh.

        Parameters
        ----------
        value : float 
            The new value to update the bubble radius setting with.

        Returns
        -------
        int

        """        
        self.settings["Bubble radius"] = value
        logger.debug(f"{self.settings['Legend']} on bubble raduys line edit changed {value}")
        self.refresh_required = True
        return

    def on_permittivity_i_sb_changed(self,value):
        """Update the imaginary part of the permittivity in the matrix material settings and trigger a refresh.

        Parameters
        ----------
        value : float
            The new imaginary part of the permittivity to be set.

        Returns
        -------
        int

        Notes
        -----
        - Updates the imaginary part of the `Matrix permittivity` in the settings dictionary with the new value, keeping
          the real part unchanged. - Replaces the permittivity object of the matrix material with a new constant scalar
          permittivity object created using the updated `Matrix permittivity`. - Marks the matrix material as manually
          defined and triggers a required refresh to update any dependent calculations or displays.

        """        
        self.refresh_required = True
        real = np.real(self.settings["Matrix permittivity"])
        self.settings["Matrix permittivity"] = complex(real,value)
        new_permittivity_object = DielectricFunction.ConstantScalar(self.settings["Matrix permittivity"])
        self.matrixMaterial.set_permittivity_object(new_permittivity_object)
        self.settings["Matrix"] = "Material defined manually"
        self.materialDefinedManually = True
        logger.debug(f"{self.settings['Legend']} on imaginary permittivity line edit changed {value}")
        self.refresh()
        self.refresh_required = True
        return

    def on_permittivity_r_sb_changed(self,value):
        """Update the real part of the permittivity in the matrix material settings and refresh the material definition.

        Parameters
        ----------
        value : float
            The new real value of the matrix permittivity to be set.

        Returns
        -------
        int

        """        
        self.refresh_required = True
        imaginary = np.imag(self.settings["Matrix permittivity"])
        self.settings["Matrix permittivity"] = complex(value,imaginary)
        new_permittivity_object = DielectricFunction.ConstantScalar(self.settings["Matrix permittivity"])
        self.matrixMaterial.set_permittivity_object(new_permittivity_object)
        self.settings["Matrix"] = "Material defined manually"
        self.materialDefinedManually = True
        logger.debug(f"{self.settings['Legend']} on permittivity line edit changed {value}")
        self.refresh()
        self.refresh_required = True
        return

    def on_atr_index_sb_changed(self,value):
        """Handle the change in settings for ATR material refractive index.

        Parameters
        ----------
        value : float
            The new value for the ATR material refractive index.

        Returns
        -------
        int

        Notes
        -----
        - This function updates the 'ATR material refractive index' in the settings dictionary.
        - Marks the instance as requiring a refresh, possibly to update some UI elements or calculations.

        """        
        self.settings["ATR material refractive index"] = value
        logger.debug(f"{self.settings['Legend']} on atr index line edit changed {value}")
        self.refresh_required = True
        return

    def on_atr_incident_ang_sb_changed(self,value):
        """Handle the update to the ATR incident angle situation.

        This method updates the ATR incident angle setting based on user inputs or changes. It also triggers an update
        to ensure the new settings are reflected across the application.

        Parameters
        ----------
        value : float or int
            The new value for the ATR incident angle.

        Returns
        -------
        int

        """        
        self.settings["ATR theta"] = value
        logger.debug(f"{self.settings['Legend']} on atr incident angle line edit changed {value}")
        self.refresh_required = True
        return

    def on_atr_spolfrac_sb_changed(self,value):
        """Update the ATR S polarisation fraction setting and request a refresh as required.

        Parameters
        ----------
        value : float
            The new value for the ATR S polarisation fraction. 

        Returns
        -------
        int

        Notes
        -----
        - This function updates the settings dictionary on the `self` object with the new value for the key 'ATR S
          polarisation fraction'. - The refresh flag `self.refresh_required` is set to `True` to indicate that some
          action is needed to reflect the change in the application.

        """        
        self.settings["ATR S polarisation fraction"] = value
        logger.debug(f"{self.settings['Legend']} on atr spolfraction line edit changed {value}")
        self.refresh_required = True
        return

    def change_greyed_out(self):
        """Modify UI elements based on the selected effective medium method.

        This function updates the enabled status of UI elements such as size, sigma, and shape selection based on the
        currently selected effective medium method in the settings. It also updates the particle shape setting if
        necessary and adjusts the UI to reflect any changes.

        Parameters
        ----------
        None

        Returns
        -------
        int

        Notes
        -----
        - The method checks the 'Effective medium method' and 'Particle shape' in the settings. - For 'Mie' and
          'Anisotropic-Mie' methods, it enables size and sigma spinboxes and sets the particle shape to 'Sphere',
          disabling all other shapes. - For 'Averaged Permittivity', it disables size and sigma spinboxes and sets the
          particle shape to 'Sphere'. - For 'Maxwell-Garnett' and 'Bruggeman', it enables the size spinbox and disables
          the sigma spinbox, allowing all shapes to be selected. - If the particle shape is set to 'Ellipsoid', 'Plate',
          or 'Needle', it enables direction (h, k, l) spinboxes and adjusts the display label accordingly. For 'Sphere',
          it disables these spinboxes.

        """        
        logger.debug(f"{self.settings['Legend']} Start:: change_greyed_out")
        if self.settings["Matrix"] == "none":
            # No matrix selected: disable all EMT and matrix-related controls
            if hasattr(self, "size_sb"):
                self.size_sb.setEnabled(False)
            if hasattr(self, "sigma_sb"):
                self.sigma_sb.setEnabled(False)
            self.shape_cb.setEnabled(False)
            if hasattr(self, "methods_cb"):
                self.methods_cb.setEnabled(False)
            self.density_sb.setEnabled(False)
            self.permittivity_r_sb.setEnabled(False)
            self.permittivity_i_sb.setEnabled(False)
            self.vf_sb.setEnabled(False)
            self.mf_sb.setEnabled(False)
            if hasattr(self, "bubble_vf_sb"):
                self.bubble_vf_sb.setEnabled(False)
            if hasattr(self, "bubble_radius_sb"):
                self.bubble_radius_sb.setEnabled(False)
            self.h_sb.setEnabled(False)
            self.k_sb.setEnabled(False)
            self.l_sb.setEnabled(False)
            self.aoverb_sb.setEnabled(False)
            logger.debug(f"{self.settings['Legend']} Finished:: change_greyed_out")
            return
        if hasattr(self, "methods_cb"):
            method = self.settings["Effective medium method"]
            if method in ( "Mie",  "Anisotropic-Mie" ):
                self.size_sb.setEnabled(True)
                self.sigma_sb.setEnabled(True)
                for i in range(len(self.shapes)):
                    self.shape_cb.model().item(i).setEnabled(False)
                self.settings["Particle shape"] = "Sphere"
                self.shape_cb.setEnabled(True)
                index = self.shape_cb.findText(self.settings["Particle shape"], Qt.MatchFixedString)
                if index >=0:
                    self.shape_cb.model().item(index).setEnabled(True)
                    self.shape_cb.setCurrentIndex(index)
                else:
                    logger.error(f"Method index was not 0 {self.settings['Particle shape']}")
            elif method == "Averaged Permittivity":
                self.size_sb.setEnabled(False)
                self.sigma_sb.setEnabled(False)
                self.settings["Particle shape"] = "Sphere"
                index = self.shape_cb.findText(self.settings["Particle shape"], Qt.MatchFixedString)
                if index >=0:
                    self.shape_cb.model().item(index).setEnabled(True)
                    self.shape_cb.setCurrentIndex(index)
                self.shape_cb.setEnabled(False)
                for i in range(len(self.shapes)):
                    self.shape_cb.model().item(i).setEnabled(False)
            elif method in ( "Maxwell-Garnett", "Bruggeman" ):
                self.size_sb.setEnabled(True)
                self.sigma_sb.setEnabled(False)
                self.shape_cb.setEnabled(True)
                for i in range(len(self.shapes)):
                    self.shape_cb.model().item(i).setEnabled(True)
            else:
                self.size_sb.setEnabled(False)
                self.sigma_sb.setEnabled(False)
                self.shape_cb.setEnabled(True)
                for i in range(len(self.shapes)):
                    self.shape_cb.model().item(i).setEnabled(True)
        else:
            self.shape_cb.setEnabled(True)
            for i in range(len(self.shapes)):
                self.shape_cb.model().item(i).setEnabled(True)
        # deal with shapes
        if self.settings["Particle shape"] == "Ellipsoid":
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText("Unique direction [abc]")
            self.aoverb_sb.setEnabled(True)
        elif self.settings["Particle shape"] == "Plate":
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText("Unique direction (hkl)")
            self.aoverb_sb.setEnabled(False)
        elif self.settings["Particle shape"] == "Needle":
            self.h_sb.setEnabled(True)
            self.k_sb.setEnabled(True)
            self.l_sb.setEnabled(True)
            self.hkl_label.setText("Unique direction [abc]")
            self.aoverb_sb.setEnabled(False)
        elif self.settings["Particle shape"] == "Sphere":
            self.h_sb.setEnabled(False)
            self.k_sb.setEnabled(False)
            self.l_sb.setEnabled(False)
            self.aoverb_sb.setEnabled(False)
        else:
            logger.warning(f"ScenarioTab: Shape not recognised {self.settings['Particle shape']}")
        # Grey out ATR widgets when not in ATR mode
        is_atr = self.spectroscopy == "Powder ATR"
        if is_atr:
            self.atr_index_sb.setEnabled(is_atr)
            self.atr_incident_ang_sb.setEnabled(is_atr)
            self.atr_spolfrac_sb.setEnabled(is_atr)
        logger.debug(f"{self.settings['Legend']} Finished:: change_greyed_out")
        return

    def calculate(self, vs_cm1):
        """Calculate the powder spectrum for the range of frequencies in vs_cm1.

        Dispatches to the appropriate calculation method based on the spectroscopy
        type set in the settings tab.

        Parameters
        ----------
        vs_cm1 : array_like
            Array of frequencies for which to calculate the powder spectrum.

        Returns
        -------
        None

        """
        spectroscopy_type = self.notebook.settingsTab.settings.get("Spectroscopy type", "Powder Infrared")
        if spectroscopy_type == "Powder Infrared":
            self._calculate_infrared(vs_cm1)
        elif spectroscopy_type == "Powder ATR":
            self._calculate_atr(vs_cm1)
        elif spectroscopy_type == "Powder Raman":
            self._calculate_raman(vs_cm1)
        else:
            logger.error(f"{self.settings['Legend']} calculate: unknown spectroscopy type: {spectroscopy_type}")

    def _calculate_atr(self, vs_cm1):
        """Calculate the ATR powder spectrum for the range of frequencies in vs_cm1.

        Calls _calculate_infrared to obtain the effective medium permittivity, then
        computes the ATR reflectance from the stored complex permittivity values.

        Parameters
        ----------
        vs_cm1 : array_like
            Array of frequencies for which to calculate the ATR spectrum.

        Returns
        -------
        None

        """
        self._calculate_infrared(vs_cm1)
        atr_refractive_index = self.settings["ATR material refractive index"]
        atr_theta = self.settings["ATR theta"]
        atr_spolfraction = self.settings["ATR S polarisation fraction"]
        self.sp_atr = []
        for real_perm, imag_perm in zip(self.realPermittivity, self.imagPermittivity):
            trace = complex(real_perm, imag_perm)
            ri = Calculator.calculate_refractive_index_scalar(trace)
            spatr = Calculator.reflectance_atr(ri, atr_refractive_index, atr_theta, atr_spolfraction)
            self.sp_atr.append(spatr)

    def _calculate_raman(self, vs_cm1):
        """Calculate the powder Raman spectrum for the range of frequencies in vs_cm1.

        Implements the macroscopic powder approach in ``Raman-Theory.tex`` for
        small ellipsoidal particles embedded in a non-absorbing matrix.

        Two code paths are used depending on particle shape:

        * **Sphere** — the depolarisation tensor L = I/3 is rotationally invariant, so
          R_eff(Ω) = R @ R_eff_crystal @ R.T is a pure rank-2 tensor rotation
          and the analytical rotational invariants (``eq-invariants1`` and
          ``eq-Intensities``) are exact.
        * **Non-sphere** — D^{particle} is diagonalised once in the crystal frame
          (ΔD is orientation-invariant); for each of ``n_samples`` SO(3) orientations
          the optical internal-field tensor is recomputed and used to evaluate the
          effective bulk Raman tensor R_eff(Ω), whose contribution is
          accumulated into the spectrum.

        Parameters
        ----------
        vs_cm1 : array_like
            Array of frequencies for which to calculate the Raman spectrum.

        Returns
        -------
        None

        Notes
        -----
        Results are stored in ``self.raman_spectrum`` as a list of intensities in
        arbitrary units, one per frequency in ``vs_cm1``.

        """
        logger.debug(f"{self.settings['Legend']} Start:: _calculate_raman")
        if not self.calculation_required:
            logger.debug(f"{self.settings['Legend']} Finished:: _calculate_raman - not required")
            return
        if self.notebook.plottingTab is None:
            logger.debug(f"{self.settings['Legend']} Finished:: _calculate_raman - plottingTab unavailable")
            return
        if self.reader is None:
            logger.debug(f"{self.settings['Legend']} Finished:: _calculate_raman - reader unavailable")
            return

        raman_tensors = self.reader.get_raman_tensors()
        if raman_tensors is None or len(raman_tensors) == 0:
            logger.warning(f"{self.settings['Legend']} _calculate_raman: no Raman tensors available")
            self.raman_spectrum = list(np.zeros(len(vs_cm1)))
            self.vs_cm1 = list(np.array(vs_cm1))
            self.calculation_required = False
            self.notebook.progressbars_update(increment=len(vs_cm1))
            return

        # Mode information from settings tab
        settings_tab = self.notebook.settingsTab
        frequencies_cm1 = settings_tab.frequencies_cm1
        sigmas_cm1 = settings_tab.sigmas_cm1
        modes_selected = settings_tab.modes_selected

        # Crystal optical permittivity ε_i^∞ (3×3)
        epsilon_inf_i = np.array(self.notebook.settingsTab.settings["Optical permittivity"], dtype=complex)
        I3 = np.eye(3, dtype=complex)

        is_none_matrix = (self.settings["Matrix"] == "none")
        if is_none_matrix:
            # No matrix: skip all EMT/field corrections.
            # N = I (no internal field), sphere path, no particle frequency correction.
            N = I3.copy()
            L = I3 / 3.0
            epsilon_e = None
            has_correction_data = False
            is_sphere = True
        else:
            # Matrix optical permittivity ε_e^∞ (scalar)
            epsilon_e = float(np.real(self.matrixMaterial.get_optical_permittivity()))
            # Depolarisation tensor L from particle shape (same logic as _calculate_infrared)
            L = self.calculate_depolarisation_tensor()
            # Internal field tensor N from eq-internal_external_fields3.
            N = Calculator.compute_internal_field_tensor(L, epsilon_inf_i, epsilon_e)
            # Particle phonon frequencies from eq-particle_dynamical and eq-particle_eigenvalues:
            # When Born charges and the hessian are available, diagonalise D^particle
            # to obtain shifted frequencies and transformed Raman tensors.
            has_hessian = hasattr(self.reader, "hessian") and self.reader.hessian is not None
            has_born = len(self.reader.born_charges) > 0
            has_normal_modes = np.any(self.reader.mass_weighted_normal_modes)
            has_correction_data = has_hessian and has_born and has_normal_modes
            is_sphere = np.allclose(np.real(L), (1.0 / 3.0) * np.eye(3), atol=1e-12)

        # Volume fraction of crystal in the sample (scales total scattering intensity).
        volume_fraction = self.settings["Volume fraction"]

        # Raman experiment parameters
        laser_nm = self.settings["Raman laser frequency"]
        validate_laser_frequency(laser_nm)
        nu_L = 1.0e7 / laser_nm          # laser frequency in cm^-1
        polarisation = self.settings["Raman laser polarisation"]
        temperature = self.settings["Raman temperature"]
        n_samples = self.settings["Raman orientation samples"]
        chi2 = getattr(self.reader, "nonlinear_optical_susceptibility", None)
        eo_requested = self.settings.get("Raman electro-optic term", True)
        include_eo = (
            eo_requested
            and not is_none_matrix
            and has_correction_data
            and chi2 is not None
        )
        if eo_requested and not include_eo:
            logger.debug(
                "%s _calculate_raman: particle EO disabled "
                "(matrix=%s, correction_data=%s, chi2=%s)",
                self.settings["Legend"],
                not is_none_matrix,
                has_correction_data,
                chi2 is not None,
            )

        vs_cm1 = np.array(vs_cm1, dtype=float)

        # Non-sphere: numerical SO(3) averaging.  The finite-particle EO
        # correction is constructed once in the particle frame and does not
        # by itself require orientation sampling.
        if not is_sphere and has_correction_data:
            logger.debug(f"{self.settings['Legend']} _calculate_raman: orientation numerical average ({n_samples} samples)")
            n_freqs = len(vs_cm1)
            _accumulated = [0.0]
            _updated = [0]
            _increment_per_sample = n_freqs / n_samples

            def _progress_callback():
                _accumulated[0] += _increment_per_sample
                new_int = int(_accumulated[0])
                if new_int > _updated[0]:
                    self.notebook.progressbars_update(increment=new_int - _updated[0])
                    _updated[0] = new_int
                    QCoreApplication.processEvents()

            spectrum = self._compute_orientation_sampled_spectrum(
                L, epsilon_e, epsilon_inf_i, I3,
                raman_tensors, frequencies_cm1, sigmas_cm1, modes_selected,
                polarisation, nu_L, temperature, n_samples, vs_cm1,
                no_matrix=is_none_matrix,
                include_eo=include_eo,
                chi2_repsilon=chi2,
                progress_callback=_progress_callback)
            # Flush any remaining fractional increments
            remaining = n_freqs - _updated[0]
            if remaining > 0:
                self.notebook.progressbars_update(increment=remaining)
            self.raman_spectrum = (spectrum * volume_fraction).tolist()
            self.vs_cm1 = list(vs_cm1)
            self.calculation_required = False
            QCoreApplication.processEvents()
            logger.debug(f"{self.settings['Legend']} Finished:: _calculate_raman")
            return

        # Sphere path (or non-sphere fallback when correction data unavailable):
        # compute particle modes once in crystal frame, then use analytical invariants.
        if has_correction_data:
            logger.debug(f"{self.settings['Legend']} _calculate_raman: applying particle frequency correction")
            loop_freqs, loop_raman, loop_sigmas, loop_selected = self._compute_particle_modes(
                N, L, epsilon_e, epsilon_inf_i, I3,
                raman_tensors, frequencies_cm1, sigmas_cm1, modes_selected,
                include_eo=include_eo, chi2_repsilon=chi2)
        else:
            logger.debug(f"{self.settings['Legend']} _calculate_raman: no hessian/born charges, using bulk TO frequencies")
            loop_freqs = np.array(frequencies_cm1)
            loop_raman = [np.array(rt, dtype=complex) for rt in raman_tensors]
            loop_sigmas = list(sigmas_cm1)
            loop_selected = list(modes_selected)

        spectrum = np.zeros(len(vs_cm1))

        # Fixed shape/crystal axes in an isotropic host rotate as one rank-two
        # effective tensor: eq-invariants1 and eq-Intensities also cover ellipsoids.
        if not is_sphere:
            logger.warning("%s _calculate_raman: polar response data unavailable; using TO tensors "
                           "with ellipsoid optical fields and analytic orientation invariants. "
                           "Particle frequency and EO corrections are omitted.", self.settings["Legend"])
        for freq, sigma, selected, R_eps in zip(loop_freqs, loop_sigmas, loop_selected, loop_raman):
            if not selected or not valid_stokes_mode(freq, nu_L, 1.0):
                continue

            # Effective particle Raman tensor.
            if is_none_matrix:
                R_eff = np.array(R_eps, dtype=complex)
            else:
                R_eff = Calculator.compute_effective_raman_tensor(R_eps, N, N)

            # Powder-averaged scattering intensity for the chosen polarisation.
            if polarisation in ("VV", "VH", "HV"):
                vv, vh = Calculator.compute_powder_raman_intensities(R_eff)
                intensity_factor = vv if polarisation == "VV" else vh
            else:  # Unpolarised: 45α² + 7γ² + 5κ² (≠ VV+VH for antisymmetric tensors)
                alpha = np.trace(R_eff) / 3.0
                gamma_t = 0.5 * (R_eff + R_eff.T) - alpha * I3
                kappa_t = 0.5 * (R_eff - R_eff.T)
                alpha2 = float(np.real(alpha * np.conj(alpha)))
                gamma2 = 3.0 / 2.0 * float(np.real(np.sum(gamma_t * np.conj(gamma_t))))
                kappa2 = 3.0 / 2.0 * float(np.real(np.sum(kappa_t * np.conj(kappa_t))))
                intensity_factor = 45.0 * alpha2 + 7.0 * gamma2 + 5.0 * kappa2

            S_m = stokes_prefactor(freq, nu_L, temperature) * intensity_factor

            # Add Lorentzian contribution to the spectrum from eq-raman-intensity.
            spectrum += S_m * sigma / ((vs_cm1 - freq) ** 2 + sigma ** 2)

        self.notebook.progressbars_update(increment=len(vs_cm1))
        self.raman_spectrum = (spectrum * volume_fraction).tolist()
        self.vs_cm1 = list(vs_cm1)
        self.calculation_required = False
        QCoreApplication.processEvents()
        logger.debug(f"{self.settings['Legend']} Finished:: _calculate_raman")

    def _compute_particle_modes(self, N_bg, L, epsilon_e, epsilon_inf_i, I3,
                                raman_tensors, frequencies_cm1, sigmas_cm1, modes_selected,
                                include_eo=False, chi2_repsilon=None):
        """Compute particle phonon frequencies and Raman tensors in the particle normal-mode basis.

        Implements ``eq-particle_dynamical`` and ``eq-particle_eigenvalues``.
        The particle dynamical matrix

        .. math::

            D^{\\mathrm{particle}} = D^{TO} + \\frac{4\\pi}{\\varepsilon_e V}
                \\left(Z^{mw}\\right)^T N_{bg} L Z^{mw}

        where :math:`Z^{mw}_{\\alpha,\\kappa\\beta} = Z_{\\alpha,\\kappa\\beta}/\\sqrt{M_\\kappa}`
        (``eq-polarisation-born-mass-weighted``) carries **no volume factor**, and
        :math:`V` is the unit-cell volume (explicit in the :math:`1/V` prefactor).
        The matrix is diagonalised to give phonon frequencies and normal modes
        appropriate for a small particle with depolarisation tensor L embedded in a medium
        with permittivity ε_e.  The bulk Raman tensors are transformed into the resulting
        particle normal-mode basis.

        Parameters
        ----------
        N_bg : ndarray, shape (3, 3)
            Background internal field tensor, analogous to ``eq-internal_external_fields3``.
        L : ndarray, shape (3, 3)
            Depolarisation tensor for the particle shape.
        epsilon_e : float
            Optical permittivity of the external medium (scalar).
        epsilon_inf_i : ndarray, shape (3, 3)
            Optical permittivity tensor of the inclusion.
        I3 : ndarray, shape (3, 3)
            3×3 identity matrix.
        raman_tensors : list of ndarray
            Bulk TO Raman tensors (one per TO mode, each 3×3).
        frequencies_cm1 : list of float
            Bulk TO frequencies in cm^{-1}.
        sigmas_cm1 : list of float
            Lorentzian half-widths in cm^{-1} (one per TO mode).
        modes_selected : list of bool
            Mode selection flags (one per TO mode).
        include_eo : bool, optional
            Apply the finite-particle electro-optic correction in the particle-mode basis.
        chi2_repsilon : ndarray, shape (3, 3, 3), optional
            Electro-optic susceptibility in the internal ``R_epsilon`` convention.

        Returns
        -------
        particle_freqs : ndarray, shape (3N,)
            Particle phonon frequencies in cm^{-1}.
        particle_raman : list of ndarray
            Raman tensors in the particle normal-mode basis (one 3×3 array per mode).
        particle_sigmas : list of float
            Lorentzian half-widths inherited from the dominant TO component.
        particle_selected : list of bool
            Mode selection flags inherited from the dominant TO component.

        Notes
        -----
        The correction uses the Gaussian-unit convention ``1/ε₀ → 4π`` that matches
        the non-analytic correction in :func:`~PDielec.Calculator.longitudinal_modes`.
        The Born charge matrix Z' = Z*/√M carries **no volume factor**; the volume
        enters only as ``1/V`` in the prefactor ``4π/(ε_e V)``.  This is consistent
        with the NAC prefactor ``4π/V`` in ``longitudinal_modes``.  The Born charges
        are in units of the elementary charge (dimensionless), masses in electron-mass
        units, and the volume in Bohr³.

        """
        nAtoms = self.reader.nions
        n_modes = 3 * nAtoms

        # Unit-cell volume in Bohr³ (reader stores Å³)
        volume_au = self.reader.volume * angs2bohr ** 3

        # Atomic masses in atomic units (electron masses)
        masses_au = np.array(self.reader.masses) * amu

        # Born effective charges Z*[κ, α, β]: (nAtoms, 3, 3)
        born_charges = np.array(self.reader.born_charges)

        # Z' matrix (3 × 3N): Z'[α, κβ] = Z*[κ, α, β] / √M_κ_au  (no volume factor)
        # This matches the Born-charge convention in Calculator.longitudinal_modes.
        # The volume enters only in the prefactor below.
        Z_mat = np.zeros((3, n_modes))
        for kappa in range(nAtoms):
            inv_sqrtM = 1.0 / math.sqrt(masses_au[kappa])
            for beta in range(3):
                Z_mat[:, kappa * 3 + beta] = born_charges[kappa, :, beta] * inv_sqrtM

        # Correction to dynamical matrix from eq-particle_dynamical:
        # ΔD = (4π / (ε_e V)) Z'^T (N_bg L) Z'
        # Derived from F^mw = (e²/(ε₀ε_e V)) Z'^T N_bg L Z' x in atomic units
        # (1/ε₀ → 4π, same Gaussian convention as Calculator.longitudinal_modes).
        # N_bg and L are real for non-absorbing media; take real part to be safe.
        NbgL = np.real(N_bg) @ np.real(L)
        delta_D = (4.0 * np.pi / (epsilon_e * volume_au)) * (Z_mat.T @ NbgL @ Z_mat)

        # Bulk TO dynamical matrix (mass-weighted hessian in atomic units)
        D_TO = np.array(self.reader.hessian, dtype=float)

        # Particle dynamical matrix and eigendecomposition.
        eig_val, eig_vec = np.linalg.eigh(D_TO + delta_D)

        # Particle frequencies in cm^{-1}; preserve sign for dynamically unstable modes
        particle_freqs = np.array([
            (math.sqrt(abs(ev)) / wavenumber) * (1.0 if ev >= 0.0 else -1.0)
            for ev in eig_val
        ])

        # Build U_TO: rows = TO eigenvectors flattened to (n_to_modes × 3N)
        n_to_modes = len(self.reader.mass_weighted_normal_modes)
        U_TO = np.zeros((n_to_modes, n_modes))
        for imode, mode in enumerate(self.reader.mass_weighted_normal_modes):
            col = 0
            for atom in mode:
                U_TO[imode, col:col + 3] = atom
                col += 3

        # Overlap matrix C[n, m] = <u_n^TO | u_m^particle>
        # eig_vec columns are the particle eigenvectors (shape 3N × 3N)
        C = U_TO @ eig_vec  # shape (n_to_modes, 3N)

        # All mode-indexed lists must have the same length; warn if not.
        n_rt = len(raman_tensors)
        n_sigma = len(sigmas_cm1)
        n_sel = len(modes_selected)
        if not (n_to_modes == n_rt == n_sigma == n_sel == n_modes):
            logger.warning(
                f"_compute_particle_modes: inconsistent list lengths — "
                f"n_to_modes={n_to_modes}, n_raman_tensors={n_rt}, "
                f"n_sigmas={n_sigma}, n_selected={n_sel}, n_modes={n_modes}"
            )

        particle_raman = []
        particle_sigmas = []
        particle_selected = []

        # loop over the particle modes
        for p_idx in range(n_modes):

            # Raman tensor in particle mode basis: R^(m,particle) = Σ_n C[n,m] R_eps^(n)
            R_p = np.zeros((3, 3), dtype=complex)
            for n_to in range(n_to_modes):
                R_p += C[n_to, p_idx] * np.array(raman_tensors[n_to], dtype=complex)

            # Store the particle raman tensor
            particle_raman.append(R_p)

            # Assign sigma and selected flag from the dominant TO mode
            # In the case of degeneracy this approach might fail
            # But any other approach seems fraught with dangers too.
            dominant_to = int(np.argmax(np.abs(C[:, p_idx])))
            particle_sigmas.append(sigmas_cm1[dominant_to])
            particle_selected.append(modes_selected[dominant_to])

        if include_eo and chi2_repsilon is not None:
            K_particle = NbgL / epsilon_e
            particle_raman = apply_particle_eo_correction(
                particle_raman, chi2_repsilon, K_particle, Z_mat, eig_vec)

        return particle_freqs, particle_raman, particle_sigmas, particle_selected

    @staticmethod
    def _get_sobol_rotations(n_samples):
        """Return ``n_samples`` rotations in SO(3) drawn uniformly via a Sobol sequence.

        Delegates to :func:`PDielec.Calculator.sobol_rotations` with the fixed
        seed used for reproducible powder averaging.

        Parameters
        ----------
        n_samples : int
            Number of orientations required.

        Returns
        -------
        list of ndarray, each shape (3, 3)
            Orthogonal rotation matrices (det = +1).

        """
        return Calculator.sobol_rotations(n_samples, seed=42)

    def _compute_orientation_sampled_spectrum(
            self, L, epsilon_e, epsilon_inf_i, I3,
            raman_tensors, frequencies_cm1, sigmas_cm1, modes_selected,
            polarisation, nu_L, temperature, n_samples, vs_cm1,
            no_matrix=False, include_eo=False, chi2_repsilon=None,
            progress_callback=None):
        """Compute the powder Raman spectrum for non-spherical particles by numerical SO(3) averaging.

        The particle dynamical matrix correction ΔD is orientation-invariant (proof: rotating
        L, ε_inf, and Z to the lab frame and back always recovers ΔD_crystal), so D^{particle}
        is diagonalised **once** before the orientation loop.  The finite-particle EO
        correction is also constructed once in the particle frame using the same
        electrostatic kernel as ΔD.  For each sampled orientation Ω ∈ SO(3), the
        optical internal-field tensor is recomputed and supplies exactly one incident
        and one scattered local-field factor.  Scattering strengths from
        ``eq-ramanefficiency_depolarised`` are accumulated as Lorentzian
        contributions from ``eq-raman-intensity`` and normalised by ``n_samples``.

        Parameters
        ----------
        L : ndarray, shape (3, 3)
            Crystal-frame depolarisation tensor.
        epsilon_e : float
            Optical permittivity of the external medium (scalar).
        epsilon_inf_i : ndarray, shape (3, 3)
            Optical permittivity tensor of the inclusion (crystal frame).
        I3 : ndarray, shape (3, 3)
            3×3 identity matrix.
        raman_tensors : list of ndarray
            Bulk TO Raman tensors (crystal frame), one 3×3 array per mode.
        frequencies_cm1 : list of float
            Bulk TO frequencies in cm⁻¹.
        sigmas_cm1 : list of float
            Lorentzian half-widths in cm⁻¹ (one per TO mode).
        modes_selected : list of bool
            Mode selection flags (one per TO mode).
        polarisation : str
            One of ``"VV"``, ``"VH"``, ``"HV"``, ``"Unpolarised"``.
        nu_L : float
            Laser frequency in cm⁻¹.
        temperature : float
            Sample temperature in K.
        n_samples : int
            Number of SO(3) orientations to sample.
        vs_cm1 : ndarray, shape (n_freqs,)
            Frequency axis for the spectrum in cm⁻¹.
        no_matrix : bool, optional
            If True, skip particle frequency and local-field corrections.
        include_eo : bool, optional
            If True, include the finite-particle electro-optic χ² correction.
        chi2_repsilon : ndarray or None, optional
            Reader χ² tensor in the internal R_epsilon convention.
        progress_callback : callable or None, optional
            If provided, called once after each orientation sample to allow the
            caller to update a progress bar.  The callable takes no arguments.

        Returns
        -------
        spectrum : ndarray, shape (n_freqs,)
            Accumulated and normalised Raman spectrum.

        Notes
        -----
        Polarisation vectors follow the backscattering geometry with the laser along Z:
        VV uses ``e_L = e_S = [0, 1, 0]``; VH uses ``e_S = [1, 0, 0]``; Unpolarised
        sums both.  The rotation ``R`` maps crystal-frame coordinates to the lab frame.

        The Born-charge α-index transforms as a vector under rotation, i.e.
        ``Z_lab[α, κβ] = Σ_γ R[α, γ] Z_crystal[γ, κβ]``, so
        ``Z_mat_lab = R @ Z_mat_crystal``.

        The Lorentzian width and mode-selection flags are inherited from the dominant
        TO component of each particle mode (same heuristic as ``_compute_particle_modes``).

        """
        # Build all orientation-independent particle-mode data through the same
        # path used by the analytical sphere calculation.
        if no_matrix:
            part_freqs = np.asarray(frequencies_cm1, dtype=float)
            R_eps_cryst_list = [np.asarray(tensor, dtype=complex) for tensor in raman_tensors]
            part_sigmas = list(sigmas_cm1)
            part_selected = list(modes_selected)
        else:
            N_phonon_bg_crystal = Calculator.compute_internal_field_tensor(
                np.real(L), np.real(epsilon_inf_i), epsilon_e)
            part_freqs, R_eps_cryst_list, part_sigmas, part_selected = self._compute_particle_modes(
                N_phonon_bg_crystal, L, epsilon_e, epsilon_inf_i, I3,
                raman_tensors, frequencies_cm1, sigmas_cm1, modes_selected,
                include_eo=include_eo, chi2_repsilon=chi2_repsilon)
        n_modes = len(part_freqs)

        # Per-mode orientation-independent scalars: (freq, sigma, spectral_weight)
        # None marks modes that should be skipped entirely.
        mode_data = []
        for p_idx in range(n_modes):
            freq = part_freqs[p_idx]
            if not valid_stokes_mode(freq, nu_L, 1.0):
                mode_data.append(None)
                continue
            if not part_selected[p_idx]:
                mode_data.append(None)
                continue
            sigma = part_sigmas[p_idx]
            mode_data.append((freq, sigma, stokes_prefactor(freq, nu_L, temperature)))

        # ── Orientation loop ──────────────────────────────────────────────────────────
        # Lab-frame polarisation vectors (backscattering, laser along Z)
        e_L  = np.array([0.0, 1.0, 0.0])
        e_VV = np.array([0.0, 1.0, 0.0])
        e_VH = np.array([1.0, 0.0, 0.0])
        spectrum = np.zeros(len(vs_cm1))
        
        # Loop over a random set of rotations
        for R in self._get_sobol_rotations(n_samples):

            # Rotate depolarisation tensor and ε_inf_i to the lab frame
            L_lab   = R @ np.real(L) @ R.T
            eps_lab = R @ np.real(epsilon_inf_i) @ R.T

            # Optical internal-field tensor.  The present approximation uses
            # the same optical permittivity at the laser and Stokes
            # frequencies, while the helper preserves their separate roles.
            N_optical_lab = I3 if no_matrix else Calculator.compute_internal_field_tensor(
                L_lab, eps_lab, epsilon_e)

            # Loop over particle modes (skip inactive ones)
            for p_idx in range(n_modes):
                if mode_data[p_idx] is None:
                    continue
                freq, sigma, spectral_weight = mode_data[p_idx]

                # Rotate crystal-frame Raman tensor to the lab frame
                R_eps_lab = R @ R_eps_cryst_list[p_idx] @ R.T

                # Effective particle Raman tensor in the lab frame.
                if no_matrix:
                    R_eff_lab = R_eps_lab
                else:
                    R_eff_lab = Calculator.compute_effective_raman_tensor(
                        R_eps_lab, N_optical_lab, N_optical_lab)

                # Polarisation-specific intensity
                if polarisation == "VV":
                    intensity_factor = 45.0 * abs(e_VV @ R_eff_lab @ e_L) ** 2
                elif polarisation in ("VH", "HV"):
                    intensity_factor = 45.0 * abs(e_VH @ R_eff_lab @ e_L) ** 2
                else:  # Unpolarised
                    intensity_factor = 45.0 * (
                        abs(e_VV @ R_eff_lab @ e_L) ** 2
                        + abs(e_VH @ R_eff_lab @ e_L) ** 2
                    )

                # Scattering strength accumulated as a Lorentzian line contribution.
                S_m = spectral_weight * intensity_factor
                spectrum += S_m * sigma / ((vs_cm1 - freq) ** 2 + sigma ** 2)

            if progress_callback is not None:
                progress_callback()

        # Normalise by number of orientations
        spectrum /= n_samples
        return spectrum

    def _calculate_infrared(self, vs_cm1):
        """Calculate the powder infrared absorption for the range of frequencies in vs_cm1.

        Parameters
        ----------
        vs_cm1 : array_like
            Array of frequencies for which to calculate the powder absorption.

        Returns
        -------
        array_like
            The calculated powder absorption for the given range of frequencies.

        """
        # Only allow a calculation if the plottingTab is defined
        logger.debug(f"{self.settings['Legend']} Start:: calculate")
        if not self.calculation_required:
            logger.debug(f"{self.settings['Legend']} Finished:: calculate - immediate return because calculation_required false")
            return
        if self.notebook.plottingTab is None:
            logger.debug(f"{self.settings['Legend']} Finished:: calculate - immediate return because plottingTab unavailable")
            return
        if self.reader is None:
            logger.debug(f"{self.settings['Legend']} Finished:: calculate - immediate return because reader unavailable")
            return
        logger.debug(f"{self.settings['Legend']} calculate - number of frequencies {len(vs_cm1)}")
        # Bypass EMT when 'none' matrix is selected: use raw DFT permittivity
        if self.settings["Matrix"] == "none":
            self._calculate_infrared_none(vs_cm1)
            return

        # Calculate the depolarisation tensor and the unique direction
        # (sets self.depolarisation and self.direction as side effects)
        self.calculate_depolarisation_tensor()

        # Get the crystal permittivity function from the settings tab
        crystalPermittivity = self.notebook.settingsTab.get_crystal_permittivity(vs_cm1)

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
        concentration = self.notebook.plottingTab.settings["cell concentration"]

        # Set the material parameters
        method = self.settings["Effective medium method"].lower()
        volume_fraction = self.settings["Volume fraction"]
        particle_size_mu = self.settings["Particle size(mu)"]
        particle_sigma_mu = self.settings["Particle size distribution sigma(mu)"]
        shape = self.settings["Particle shape"].lower()
        bubble_vf = self.settings["Bubble volume fraction"]
        bubble_radius = self.settings["Bubble radius"]
        # Use the pool of processors already available
        # define a partial function to use with the pool
        partial_function = partial(Calculator.solve_effective_medium_equations, method,volume_fraction,particle_size_mu,particle_sigma_mu,self.matrixPermittivityFunction,shape,self.depolarisation,concentration,bubble_vf,bubble_radius,previous_solution_shared)
        if self.notebook.pool is None:
            self.notebook.start_pool()
        logger.debug("About to use the pool to calculate effective medium equations")
        results = []
        for result in self.notebook.pool.imap(partial_function, zip(vs_cm1,crystalPermittivity), chunksize=20):
            results.append(result)
            self.notebook.progressbars_update()
        QCoreApplication.processEvents()
        self.realPermittivity = []
        self.imagPermittivity = []
        self.absorptionCoefficient = []
        self.molarAbsorptionCoefficient = []
        self.vs_cm1 = []
        logger.debug("Extracting results")
        for v,_method,_size_mu,_size_sigma,_shape,_data,trace,absorption_coefficient,molar_absorption_coefficient in results:
             self.realPermittivity.append(np.real(trace))
             self.imagPermittivity.append(np.imag(trace))
             self.absorptionCoefficient.append(absorption_coefficient)
             self.molarAbsorptionCoefficient.append(molar_absorption_coefficient)
             self.vs_cm1.append(v)
        self.calculation_required = False
        QCoreApplication.processEvents()
        logger.debug(f"{self.settings['Legend']} Finished:: calculate")
        return

    def _calculate_infrared_none(self, vs_cm1):
        """Calculate the powder infrared spectrum directly from DFT data, bypassing EMT.

        Used when the matrix material is 'none'.  The isotropic average of the
        DFT crystal permittivity tensor (trace/3) is used directly as the
        effective permittivity; no effective medium correction is applied.
        The absorption coefficient is derived from the imaginary refractive
        index in the same way as for the normal EMT path.

        Parameters
        ----------
        vs_cm1 : array_like
            Array of frequencies in cm^-1.

        Returns
        -------
        None

        """
        logger.debug(f"{self.settings['Legend']} Start:: _calculate_infrared_none")
        # Ensure self.depolarisation is a valid 3×3 matrix for spreadsheet writing.
        # For the "none" path there is no particle shape, so use sphere (L = I/3).
        self.depolarisation = Calculator.initialise_sphere_depolarisation_matrix()
        self.direction = np.array([])
        crystalPermittivity = self.notebook.settingsTab.get_crystal_permittivity(vs_cm1)
        concentration = self.notebook.plottingTab.settings["cell concentration"]
        volume_fraction = self.settings["Volume fraction"]
        self.realPermittivity = []
        self.imagPermittivity = []
        self.absorptionCoefficient = []
        self.molarAbsorptionCoefficient = []
        self.vs_cm1 = []
        for v_cm1, crystal_perm in zip(vs_cm1, crystalPermittivity):
            trace = (crystal_perm[0, 0] + crystal_perm[1, 1] + crystal_perm[2, 2]) / 3.0
            self.realPermittivity.append(float(np.real(trace)))
            self.imagPermittivity.append(float(np.imag(trace)))
            # Use Im(ε) directly rather than Im(√ε): Im(ε) from the DFT Lorentzian model
            # is a pure sum of Lorentzians centred at the TO frequencies.  Computing
            # Im(√ε) instead introduces the full Kramers-Kronig reststrahlen lineshape
            # (asymmetric peaks, frequency shifts) which is not the intended baseline.
            # Scale by volume_fraction so that the raw absorption_coefficient represents
            # the mixture (consistent with the EMT path), and divide molar by vf too
            # so the intrinsic molar absorption coefficient is vf-independent.
            absorption_coefficient = v_cm1 * 4 * np.pi * np.imag(trace) * math.log10(math.e) * volume_fraction
            molar_absorption_coefficient = absorption_coefficient / concentration / volume_fraction if concentration > 1e-30 else 0.0
            self.absorptionCoefficient.append(absorption_coefficient)
            self.molarAbsorptionCoefficient.append(molar_absorption_coefficient)
            self.vs_cm1.append(v_cm1)
            self.notebook.progressbars_update()
        self.calculation_required = False
        QCoreApplication.processEvents()
        logger.debug(f"{self.settings['Legend']} Finished:: _calculate_infrared_none")
        return

    def get_result(self, vs_cm1, plot_type):
        """Return a result for plotting.

        Parameters
        ----------
        vs_cm1 : list
            A list of frequencies in cm-1.
        plot_type : str
            A string indicating the type of result to return. One of: 'Powder Molar Absorption', 'Powder Absorption',
            'Powder Real Permittivity', 'Powder Imaginary Permittivity', 'Powder ATR'.


        Returns
        -------
        list
            A list of the required results.

        """
        logger.debug(f"{self.settings['Legend']} Start:: get_result")
        self.get_results(vs_cm1)
        logger.debug(f"{self.settings['Legend']} Finished:: get_result")
        return { "Powder Molar Absorption"      : self.molarAbsorptionCoefficient,
                 "Powder Absorption"            : self.absorptionCoefficient,
                 "Powder Real Permittivity"     : self.realPermittivity,
                 "Powder Imaginary Permittivity": self.imagPermittivity,
                 "Powder ATR"                   : self.sp_atr,
                 "Powder Raman"                 : self.raman_spectrum,
                }.get(plot_type)


    def get_results(self, vs_cm1):
        """Return the results of the effective medium theory calculation.

        This routine jut causes a refresh and calculation if it is needed.
        See get_result() for accessing the results themselves

        Parameters
        ----------
        vs_cm1 : list
            A list of frequencies in cm-1

        Returns
        -------
        None

        """
        logger.debug(f"{self.settings['Legend']} Start:: get_results {len(vs_cm1)}")
        if len(vs_cm1) > 0 and ( self.refresh_required or len(self.vs_cm1) != len(vs_cm1) or self.vs_cm1[0] != vs_cm1[0] or self.vs_cm1[1] != vs_cm1[1] ) :
            logger.debug(f"{self.settings['Legend']} get_results recalculating")
            self.refresh()
            self.calculate(vs_cm1)
        else:
            logger.debug(f"{self.settings['Legend']} get_results no need for recalculation")
            self.notebook.progressbars_update(increment=len(vs_cm1))
        logger.debug(f"{self.settings['Legend']} Finished:: get_results {len(vs_cm1)}")
        return


    def refresh(self,force=False):
        """Refresh the GUI interface with up to date values.

        Parameters
        ----------
        force : bool
            If True, then a refresh is forced regardless of the state of the scenario.

        Returns
        -------
        None

        """
        logger.debug(f"{self.settings['Legend']} Start:: refresh, force = {force}")
        if not self.refresh_required and not force:
            logger.debug(f"{self.settings['Legend']} Finished:: refresh aborted {self.refresh_required} {force}")
            return
        # Force a recalculation
        self.calculation_required = True
        # Sync spectroscopy from global spectroscopy type setting
        self.spectroscopy = self.notebook.settingsTab.settings.get("Spectroscopy type", "Powder Infrared")
        # First see if we can get the reader from the mainTab
        self.reader = self.notebook.mainTab.reader
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        self.refresh_infrared()
        if self.spectroscopy == "Powder ATR":
           self.refresh_atr()
        if self.spectroscopy == "Powder Raman":
           self.refresh_raman()
        self.change_greyed_out()
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.refresh_required = False
        logger.debug(f"{self.settings['Legend']} Finished:: refresh, force = {force}")
        return

    def refresh_atr(self):
        """Refresh the atr settings in the GUI.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.atr_index_sb.setValue(self.settings["ATR material refractive index"])
        self.atr_incident_ang_sb.setValue(self.settings["ATR theta"])
        self.atr_spolfrac_sb.setValue(self.settings["ATR S polarisation fraction"])
        return

    def on_laser_frequency_sb_changed(self, value):
        """Handle a change to the Raman laser frequency.

        Parameters
        ----------
        value : float
            The laser wavelength in nm.

        """
        logger.debug(f"{self.settings['Legend']} on_laser_frequency_sb_changed {value}")
        self.refresh_required = True
        self.settings["Raman laser frequency"] = value

    def on_polarisation_cb_activated(self, index):
        """Handle a change to the Raman polarisation configuration.

        Parameters
        ----------
        index : int
            Index into the polarisation list ["VV", "VH", "HV", "Unpolarised"].

        """
        polarisations = ["VV", "VH", "HV", "Unpolarised"]
        logger.debug(f"{self.settings['Legend']} on_polarisation_cb_activated {polarisations[index]}")
        self.refresh_required = True
        self.settings["Raman laser polarisation"] = polarisations[index]

    def on_temperature_sb_changed(self, value):
        """Handle a change to the sample temperature.

        Parameters
        ----------
        value : float
            The temperature in K.

        """
        logger.debug(f"{self.settings['Legend']} on_temperature_sb_changed {value}")
        self.refresh_required = True
        self.settings["Raman temperature"] = value

    def on_eo_term_cb_toggled(self, checked):
        """Handle a toggle of the Raman electro-optic tensor contribution."""
        logger.debug(f"{self.settings['Legend']} on_eo_term_cb_toggled {checked}")
        self.refresh_required = True
        self.calculation_required = True
        self.settings["Raman electro-optic term"] = checked

    def on_orientation_samples_cb_activated(self, index):
        """Handle a change to the number of orientation samples.

        Parameters
        ----------
        index : int
            Index into the orientation sample count list.

        """
        counts = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
        value = counts[index]
        logger.debug(f"{self.settings['Legend']} on_orientation_samples_cb_activated {value}")
        self.refresh_required = True
        self.settings["Raman orientation samples"] = value

    def refresh_raman(self):
        """Refresh the raman settings in the GUI.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.laser_frequency_sb.setValue(self.settings["Raman laser frequency"])
        polarisations = ["VV", "VH", "HV", "Unpolarised"]
        pol = self.settings["Raman laser polarisation"]
        if pol in polarisations:
            self.polarisation_cb.setCurrentIndex(polarisations.index(pol))
        self.temperature_sb.setValue(self.settings["Raman temperature"])
        self.eo_term_cb.setChecked(self.settings.get("Raman electro-optic term", True))
        counts = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
        current = self.settings["Raman orientation samples"]
        idx = counts.index(current) if current in counts else counts.index(256)
        self.orientation_samples_cb.setCurrentIndex(idx)
        return

    def refresh_infrared(self):
        """Refresh the infrared settings in the GUI.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        #
        # Update the database 
        #
        if self.settings["Materials database"] != self.DataBase.get_file_name():
            self.DataBase = MaterialsDataBase(self.settings["Materials database"])
            self.settings["Materials database"] = self.DataBase.get_file_name()
            self.database_le.setText(self.settings["Materials database"])
            # Update the possible  material names from the database
            self.material_names = self.DataBase.get_sheet_names()
        self.materialDefinedManually = False
        if self.settings["Matrix"] == "Material defined manually":
            self.materialDefinedManually = True
        if self.materialDefinedManually:
            self.material_names.append("Material defined manually")
            self.matrixMaterial = Materials.Constant("manual",permittivity=self.settings["Matrix permittivity"],density=self.settings["Matrix density"])
            self.matrixPermittivityFunction = self.matrixMaterial.get_permittivity_function()
        elif self.settings["Matrix"] in self.material_names:
            self.matrixMaterial = self.DataBase.get_material(self.settings["Matrix"])
            self.matrixPermittivityFunction = self.matrixMaterial.get_permittivity_function()
            self.settings["Matrix permittivity"] = self._matrix_display_permittivity()
            self.settings["Matrix density"] = self.matrixMaterial.get_density()
        else:
            logger.error(f"Error: matrix {self.settings['Matrix']} not available in database; available materials are: {self.material_names}")
            sys.exit()
        #
        # Reset the matrix combo box with new names
        #
        self.matrix_cb.clear()
        self.matrix_cb.addItems(self.material_names)
        index = self.matrix_cb.findText(self.settings["Matrix"], Qt.MatchFixedString)
        self.matrix_cb.setCurrentIndex(index)
        #
        # Update the matrix material information
        #
        text = self.matrixMaterial.get_information()
        self.matrix_info_le.setText(text)
        #
        # Set the matrix density widget
        #
        self.density_sb.setValue(self.settings["Matrix density"])
        #
        # Set the matrix permittivity widget
        #
        self.permittivity_r_sb.setValue(np.real(self.settings["Matrix permittivity"]))
        self.permittivity_i_sb.setValue(np.imag(self.settings["Matrix permittivity"]))
        #
        # Set the bubble data
        #
        if hasattr(self, "bubble_vf_sb"):
            self.bubble_vf_sb.setValue(100*self.settings["Bubble volume fraction"])
        if hasattr(self, "bubble_radius_sb"):
            self.bubble_radius_sb.setValue(self.settings["Bubble radius"])
        #
        # Set the volume/mass fraction
        #
        if self.settings["Mass or volume fraction"] == "volume":
            # volume fraction takes precedence
            self.update_mf_sb()
            self.update_vf_sb()
        else:
            # mass fraction takes precedence
            self.update_vf_sb()
            self.update_mf_sb()
        #
        # Set the effective medium theory method
        #
        if hasattr(self, "methods_cb"):
            index = self.methods_cb.findText(self.settings["Effective medium method"], Qt.MatchFixedString)
            self.methods_cb.setCurrentIndex(index)
        #
        # Set the particle size and distribution
        #
        if hasattr(self, "size_sb"):
            self.size_sb.setValue(self.settings["Particle size(mu)"])
        if hasattr(self, "sigma_sb"):
            self.sigma_sb.setValue(self.settings["Particle size distribution sigma(mu)"])
        #
        # Set the particle shape
        #
        index = self.shape_cb.findText(self.settings["Particle shape"], Qt.MatchFixedString)
        self.shape_cb.setCurrentIndex(index)
        self.h_sb.setValue(self.settings["Unique direction - h"])
        self.k_sb.setValue(self.settings["Unique direction - k"])
        self.l_sb.setValue(self.settings["Unique direction - l"])
        self.aoverb_sb.setValue(self.settings["Ellipsoid a/b"])
        self.legend_le.setText(self.settings["Legend"])
        self.aoverb = self.settings["Ellipsoid a/b"]
        return

    def initialise_raman_gui(self, vbox, form):
        """Initialise the GUI for atr calculations.

        Parameters
        ----------
        vbox : QVBoxLayout
            vbox is defined by the PowderScenarioTab initialised
        form : QFormLayout
            form is defined by the PowderScenarioTab initialised

        Returns
        -------
        vbox, form

        """
        #
        # Laser light frequency
        #
        self.laser_frequency_sb = QDoubleSpinBox(self)
        self.laser_frequency_sb.setRange(0.001, 2000.0)
        self.laser_frequency_sb.setSingleStep(0.01)
        self.laser_frequency_sb.setDecimals(3)
        self.laser_frequency_sb.setToolTip("Define the laser frequency in nm")
        self.laser_frequency_sb.setValue(self.settings["Raman laser frequency"])
        self.laser_frequency_sb.valueChanged.connect(self.on_laser_frequency_sb_changed)
        label = QLabel("Laser frequency (nm)", self)
        label.setToolTip("Define the laser frequency in nm")
        form.addRow(label, self.laser_frequency_sb)
        #
        # Raman polarisation configuration
        #
        self.polarisation_cb = QComboBox(self)
        self.polarisation_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.polarisation_cb.setToolTip("Define the Raman polarisation geometry: VV (parallel), VH/HV (cross), Unpolarised")
        self.polarisation_cb.addItems( ["VV", "VH", "HV", "Unpolarised" ] )
        index = self.polarisation_cb.findText(self.settings["Raman laser polarisation"], Qt.MatchFixedString)
        if index >=0:
            self.polarisation_cb.setCurrentIndex(index)
        self.polarisation_cb.activated.connect(self.on_polarisation_cb_activated)
        label = QLabel("Raman laser polarisation", self)
        label.setToolTip("Define the experimental polarisation configuration VV/VH/HV or Unpolarised")
        form.addRow(label, self.polarisation_cb)
        #
        # Sample temperature for Bose-Einstein factor
        #
        self.temperature_sb = QDoubleSpinBox(self)
        self.temperature_sb.setRange(0.01, 10000.0)
        self.temperature_sb.setSingleStep(10.0)
        self.temperature_sb.setDecimals(1)
        self.temperature_sb.setToolTip("Sample temperature in K (used for the Bose-Einstein occupation factor)")
        self.temperature_sb.setValue(self.settings["Raman temperature"])
        self.temperature_sb.valueChanged.connect(self.on_temperature_sb_changed)
        label = QLabel("Temperature (K)", self)
        label.setToolTip("Sample temperature in K (used for the Bose-Einstein occupation factor)")
        form.addRow(label, self.temperature_sb)
        #
        # Electro-optic Raman tensor correction
        #
        self.eo_term_cb = QCheckBox(self)
        self.eo_term_cb.setChecked(self.settings["Raman electro-optic term"])
        self.eo_term_cb.toggled.connect(self.on_eo_term_cb_toggled)
        self.eo_term_cb.setToolTip(
            "Include the finite-particle electro-optic χ⁽²⁾ contribution using the same "
            "depolarisation-field kernel as the particle phonon frequencies."
        )
        label = QLabel("Include electro-optic term", self)
        label.setToolTip(self.eo_term_cb.toolTip())
        form.addRow(label, self.eo_term_cb)
        #
        # Orientation samples for numerical powder averaging (non-spherical particles)
        #
        counts = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
        self.orientation_samples_cb = QComboBox(self)
        for c in counts:
            self.orientation_samples_cb.addItem(str(c))
        current = self.settings["Raman orientation samples"]
        idx = counts.index(current) if current in counts else counts.index(256)
        self.orientation_samples_cb.setCurrentIndex(idx)
        self.orientation_samples_cb.setToolTip(
            "Number of SO(3) orientations used for numerical powder averaging "
            "(applies to non-spherical particles only; powers of 2 give optimal Sobol coverage)"
        )
        self.orientation_samples_cb.activated.connect(self.on_orientation_samples_cb_activated)
        label = QLabel("Orientation samples", self)
        label.setToolTip(
            "Number of SO(3) orientations for numerical powder averaging "
            "(non-spherical particles only; 256 is usually sufficient)"
        )
        form.addRow(label, self.orientation_samples_cb)
        return vbox, form

    def initialise_raman_settings(self):
        """Set the settings attribute for raman calculations.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.settings["Raman laser frequency"] = 785
        self.settings["Raman laser polarisation"] = "VH"
        self.settings["Raman temperature"] = 298.0
        self.settings["Raman orientation samples"] = 512
        self.settings["Raman electro-optic term"] = False
        self.raman_spectrum = []
        return

    def calculate_depolarisation_tensor(self):
        """Calculate the depolarisation tensor.

        Parameters
        ----------
        None

        Returns
        -------
        L : np.array 3x3

        """
        cell = self.reader.get_unit_cell()
        shape = self.settings["Particle shape"]
        hkl = [self.settings["Unique direction - h"],
               self.settings["Unique direction - k"],
               self.settings["Unique direction - l"]]
        if shape == "Ellipsoid":
            self.direction = cell.convert_abc_to_xyz(hkl)
            self.direction = self.direction / np.linalg.norm(self.direction)
            self.depolarisation = Calculator.initialise_ellipsoid_depolarisation_matrix(
                self.direction, self.settings["Ellipsoid a/b"])
        elif shape == "Plate":
            self.direction = cell.convert_hkl_to_xyz(hkl)
            self.direction = self.direction / np.linalg.norm(self.direction)
            self.depolarisation = Calculator.initialise_plate_depolarisation_matrix(self.direction)
        elif shape == "Needle":
            self.direction = cell.convert_abc_to_xyz(hkl)
            self.direction = self.direction / np.linalg.norm(self.direction)
            self.depolarisation = Calculator.initialise_needle_depolarisation_matrix(self.direction)
        else:
            self.depolarisation = Calculator.initialise_sphere_depolarisation_matrix()
            self.direction = np.array([])
        return self.depolarisation
