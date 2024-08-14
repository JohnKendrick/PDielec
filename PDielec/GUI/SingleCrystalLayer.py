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
"""SingleCrystalLayer module."""
# -*- coding: utf8 -*-
import sys

import numpy as np
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from PDielec.Utilities import Debug


class SingleCrystalLayer:
    """A class representing a single crystal layer.

    Parameters
    ----------
    material : Material
        An instance of Material object (see :class:`~PDielec.Materials.Material`), which contains:
          - name (str)
          - density (float)
          - permittivity (float)
          - unit cell (UnitCell)
    hkl : list of int
        A list of 3 integers defining the plane of interest.
    azimuthal : float
        The azimuthal angle of rotation of the crystal about the z-axis.
    thickness : float
        The thickness of the layer in the specified thickness units.
    thicknessUnits : {'ang','nm', 'um', 'mm', 'cm'}
        The units of thickness measurement.
    incoherentOption : {'Coherent', 'Incoherent (intensity)', 'Incoherent (phase cancelling)', 'Incoherent (phase averaging)', 'Incoherent non-reflective'}
        The option for handling incoherent scattering. 
    dielectricFlag : bool
        True if the layer material is the dielectric being studied (ie it was read in as a DFT calculation), False otherwise.

    Attributes
    ----------
    material : Material
        The material of the layer. (See :class:`~PDielec.Materials.Material`)
    hkl : list of int
        The Miller indices of the plane of interest in the crystal.
    incoherentOption : str
        The incoherent option for the layer.
    azimuthal : float
        The azimuthal angle of rotation of the crystal about the z-axis, in degrees.
    thickness : float
        The thickness of the layer in the specified units.
    thicknessUnit : str
        The units used for the thickness measurement.
    dielectricFlag : bool
        A flag indicating if the layer material is the dielectric material.
    euler : ndarray
        The Euler rotation matrix for the crystal to laboratory frame transformation, stored as a 3x3 numpy array with `dtype=np.longdouble`.
    euler_inverse : ndarray
        The inverse of the Euler rotation matrix, stored as a 3x3 numpy array with `dtype=np.longdouble`.
    labframe_w : 
        Placeholder attribute for GUI interactions, not initialized within the class definition.
    labframe : ndarray
        The lattice vectors of the crystal in the laboratory frame, stored as a numpy array.
    phaseShift : float
        The phase shift for the layer.

    Methods
    -------
    print()
        Prints the main attributes of the layer.
    getPhaseShift()
        Returns the phase shift for the layer.
    setPhaseShift(phaseShift)
        Sets the phase shift for the layer.
    setAzimuthal(angle)
        Sets the azimuthal angle of the layer and recalculates the Euler matrix.
    isCoherent()
        Returns True if the layer is coherent, False otherwise.
    getIncoherentOption()
        Returns the incoherent interaction option for the layer.
    setIncoherentOption(option)
        Sets the incoherent interaction option for the layer.
    getAzimuthal()
        Returns the azimuthal angle.
    setThickness(thickness)
        Sets the thickness of the layer.
    getThickness()
        Returns the thickness of the layer.
    getThicknessInMetres()
        Converts and returns the thickness of the layer in meters.
    setThicknessUnit(thicknessUnit)
        Sets the unit of thickness measurement.
    getThicknessUnit()
        Returns the unit of thickness measurement.
    setHKL(hkl)
        Sets the Miller indices for the plane of interest in the crystal.
    getHKL()
        Returns the Miller indices of the plane of interest.
    getName()
        Returns the name of the material.
    getMaterial()
        Returns the material object associated with the layer.
    getPermittivityFunction()
        Returns the permittivity tensor for the crystal at the surface and azimuthal angle specified.
    isTensor()
        Returns True if the material of the layer is a tensor material, False otherwise.
    isScalar()
        Returns True if the material of the layer is a scalar material, False otherwise.
    isDielectric()
        Returns True if the material of the layer is the dielectric being studied.
    calculate_euler_matrix()
        Calculates the Euler angles for the crystal to lab frame transformation.
    azimuthalRotationMatrix(angle)
        Calculates a rotation matrix for the azimuthal angle.
    changeLabFrameInfo()
        Changes the Lab Frame Information, used with GUI interactions.
    getLabFrame()
        Returns the information of the crystal axes in the laboratory frame.
    calculateLabFrameEpsilon(v)
        Computes the permittivity tensor for the material at a given frequency and transforms it to the laboratory frame.
    invert(m)
        Calculates the inverse of a given 3x3 complex matrix.

    Notes
    -----
    The class heavily relies on numpy for various calculations including matrix manipulations and rotation calculations.

    """

    def __init__(self,material,hkl=None,azimuthal=0.0,thickness=0.0,thicknessUnit="nm",incoherentOption="Coherent",dielectricFlag=False):
        """Initialise a single crystal layer.

        Parameters
        ----------
        material : Material
            An instance of :class:`~PDielec.Materials.Material`, which contains:
              - name (str)
              - density (float)
              - permittivity (float)
              - unit cell (UnitCell)
        hkl : list of ints
            A list of 3 integers defining the plane of interest. Default is [0,0,1]
            If the material is isotropic then [0,0,0] is used
        azimuthal : float
            The azimuthal angle of rotation of the crystal about the z-axis.
        thickness : float
            The thickness of the layer in the specified thickness units.
        thicknessUnit : {'ang','nm', 'um', 'mm', 'cm'}
            The units of thickness measurement.
        incoherentOption : {'Coherent', 'Incoherent (intensity)', 'Incoherent (phase cancelling)', 'Incoherent (phase averaging)', 'Incoherent non-reflective'}
            The option for handling incoherent scattering. 
        dielectricFlag : bool
            True if the layer material is the dielectric being studied (ie it was read in as a DFT calculation), False otherwise.

        """
        if hkl is None:
            hkl = [0,0,1]
        self.material = material
        if material.isScalar():
            hkl = [0,0,0]
        self.hkl = hkl
        self.incoherentOption = incoherentOption
        self.azimuthal = azimuthal
        self.thickness = thickness
        self.thicknessUnit = thicknessUnit
        self.dielectricFlag = dielectricFlag
        self.euler = np.zeros((3,3),dtype=np.longdouble)
        self.euler_inverse = np.zeros((3,3),dtype=np.longdouble)
        self.labframe_w = None
        self.labframe = None
        self.phaseShift = 0.0
        self.calculate_euler_matrix()

    def print(self):
        """Print the main attributes of the layer.

        Parameters
        ----------
        None

        """
        print("------------------ ")
        print("Material         : ", self.material)
        print("HKL              :      ", self.hkl)
        print("Azimuthal        :", self.azimuthal)
        print("Thickness        :", self.thickness)
        print("Thickness unit   :", self.thicknessUnit)
        print("Dielectric flag  :", self.dielectricFlag)
        print("Incoherent option:", self.incoherentOption)
        print("Phase shift      :", self.phaseShift)

    def getPhaseShift(self):
        """Get the phase shift for the layer.

        Parameters
        ----------
        None

        Returns
        -------
        float
            The phase shift used in the averaged phase shift method

        """
        return self.phaseShift

    def setPhaseShift(self, phaseShift):
        """Set the phase shift for the layer.

        Parameters
        ----------
        phaseShift : float
           The phase shift

        Returns
        -------
        None

        """
        self.phaseShift = phaseShift
        return

    def setAzimuthal(self, angle):
        """Set the azimuthal angle of the layer and calculate the Euler matrix.

        Parameters
        ----------
        angle : float
           The azimuthal angle in degrees

        Returns
        -------
        None

        """
        self.azimuthal = angle
        self.calculate_euler_matrix()
        return 

    def isCoherent(self):
        """Return True if this is a coherent layer, False otherwise.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if this is a coherent layer, False otherwise.

        """
        result = False
        if self.incoherentOption == "Coherent":
            result = True
        return result

    def getIncoherentOption(self):
        """Return the value of the incoherent option attribute.

        This will return 'Coherent' for a coherent layer

        Parameters
        ----------
        None

        Returns
        -------
        type
            The value of the `incoherentOption` attribute of the object.

        """        
        return self.incoherentOption

    def setIncoherentOption(self,option):
        """Set the incoherent option.

        Parameters
        ----------
        option : string
            Sets the incoherent option.  Can be one of;
            - "Coherent"
            - "Incoherent (intensity)"
            - "Incoherent (phase averaging)"
            - "Incoherent (phase cancelling)"
            - "Incoherent (non-reflective)"

        Returns
        -------
        None

        """
        self.incoherentOption = option
        return 

    def getAzimuthal(self):
        """Get the azimuthal angle.

        Parameters
        ----------
        None

        Returns
        -------
        float
            The azimuthal angle.

        """
        return self.azimuthal

    def setThickness(self, thickness):
        """Set the thickness.

        Parameters
        ----------
        thickness : float
            The value of the layer thickness in the current thickness units

        Returns
        -------
        None

        """
        self.thickness = thickness
        return 

    def getThickness(self):
        """Get the thickness.

        Parameters
        ----------
        None

        Returns
        -------
        float

        """
        return self.thickness

    def getThicknessInMetres(self):
        """Get the thickness in metres.

        Parameters
        ----------
        None

        Returns
        -------
        float

        """
        thickness_conversion_factors = {"ang":1.0E-10, "nm":1.0E-9, "um":1.0E-6, "mm":1.0E-3, "cm":1.0E-2}
        tom = thickness_conversion_factors[self.thicknessUnit]
        return tom*self.thickness

    def setThicknessUnit(self, thicknessUnit):
        """Set the thickness unit.

        Parameters
        ----------
        thicknessUnit : str
            The unit of thickness. Can be one of "ang", "nm", "um", "mm", or "cm"

        Returns
        -------
        None

        """
        self.thicknessUnit = thicknessUnit
        return 

    def getThicknessUnit(self):
        """Get the thickness unit.

        Parameters
        ----------
        None

        Returns
        -------
        str

        """
        return self.thicknessUnit

    def setHKL(self, hkl):
        """Set the hkl and recalculate the Euler matrix.

        Parameters
        ----------
        hkl : list of 3 ints
           The hkl parameters of the surface

        Returns
        -------
        None

        """
        self.hkl = hkl
        self.calculate_euler_matrix()
        return 

    def getHKL(self):
        """Get the hkl.

        Parameters
        ----------
        None

        Returns
        -------
        list of 3 ints

        """
        return self.hkl

    def getName(self):
        """Return the material name.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The name of the material.

        """
        return self.material.getName()

    def getMaterial(self):
        """Return the material.

        Parameters
        ----------
        None

        Returns
        -------
        an instance of a Materials object

        """
        return self.material

    def getPermittivityFunction(self):
        """Return the permittivity function associated with this material.

        Parameters
        ----------
        None

        Returns
        -------
        A permittivity function

        """
        return self.calculateLabFrameEpsilon

    def isTensor(self):
        """Return true if the material of the layer is a tensor material.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the material of the layer is a tensor material, False otherwise.

        """
        return self.getMaterial().isTensor()

    def isScalar(self):
        """Return true if the material of the layer is a scalar material.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the material of the layer is a scalar material, False otherwise.

        """
        return self.getMaterial().isScalar()

    def isDielectric(self):
        """Return true if the material of the layer is the dielectric material.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the material of the layer is the dielectric material, otherwise False.

        """
        return self.dielectricFlag

    def calculate_euler_matrix(self):
        """Calculate the Euler angles for the crystal to lab transformation.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # Get plane specification
        hkl = self.hkl
        sum2 = hkl[0]*hkl[0] + hkl[1]*hkl[1] + hkl[2]*hkl[2]
        if sum2 < 1:
            return 
        x = 0
        y = 1
        z = 2
        # convert normal to plane to a direction in xyz coordinates
        planez = self.material.cell.convert_hkl_to_xyz(hkl)
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
        self.euler = rotation.T
        # Rotate by azimuthal angle
        self.euler = np.matmul(self.azimuthalRotationMatrix(self.azimuthal),self.euler)
        # Calculate and keep euler inverse
        self.euler_inverse = self.invert(self.euler)
        # Calculate lab frame
        self.labframe = np.matmul(self.euler,self.material.cell.lattice.T).T
        normal_to_plane_lab = np.matmul(self.euler,plane[z])
        if normal_to_plane_lab[2] < 0.9999 and normal_to_plane_lab[2] > -0.9999:
            print("Error in Euler rotations - surface normal is not along Z-axis", normal_to_plane_lab)
            sys.exit()
        return 

    def azimuthalRotationMatrix(self,angle):
        """Calculate a rotation matrix for the azimuthal angle (in degrees).

        Parameters
        ----------
        angle : float
            The azimuthal angle in degrees

        Returns
        -------
        ndarray 3x3
            The rotation matrix for the given azimuthal angle.

        """
        angle = np.radians(angle)
        matrix = np.eye(3)
        # in numpy the first index is the column, the second is the row
        matrix[0,0] = np.cos(angle)
        matrix[1,1] = np.cos(angle)
        matrix[0,1] = -np.sin(angle)
        matrix[1,0] = +np.sin(angle)
        return matrix

    def changeLabFrameInfo(self):
        """Change the Lab Frame Info panel.

        Information about the relationship between the crystal and labframe coordinate systems is presented.
        """
        if self.labframe_w is None:
            return
        a,b,c = self.getLabFrame()
        self.labframe_w.clear()
        self.labframe_w.addItem(f"crystal a-axis in lab frame: {a[0]: 3.5f}, {a[1]: 3.5f}, {a[2]: 3.5f}" )
        self.labframe_w.addItem(f"crystal b-axis in lab frame: {b[0]: 3.5f}, {b[1]: 3.5f}, {b[2]: 3.5f}" )
        self.labframe_w.addItem(f"crystal c-axis in lab frame: {c[0]: 3.5f}, {c[1]: 3.5f}, {c[2]: 3.5f}" )
        return

    def getLabFrame(self):
        """Return the a, b, and c cell axes in the laboratory frame coordinates.

        Parameters
        ----------
        None

        Returns
        -------
        (a, b, c) : 3 lists of floats each with 3 elements
            - a, The a cell axis in the laboratory frame.
            - b, The b cell axis in the laboratory frame.
            - c, The c cell axis in the laboratory frame.

        """
        return self.labframe

    def calculateLabFrameEpsilon(self, v):
        """Get the permittivity for this material at frequency v and transform it to the laboratory frame.

        Parameters
        ----------
        v : float
            Frequency at which the permittivity is evaluated.

        Returns
        -------
        float
            The permittivity of the material in the laboratory frame.

        """
        epsilon_xstal = self.material.getPermittivityFunction()(v)
        if self.material.isScalar():
            self.labFrameEpsilon = epsilon_xstal * np.eye(3)
        else:
            self.labFrameEpsilon = np.matmul(self.euler, np.matmul(epsilon_xstal,self.euler_inverse))
        return self.labFrameEpsilon

    def invert(self,m):
        """Calculate the inverse of m.

        Parameters
        ----------
        m : complex 3x3 array
            The input matrix for which to calculate the inverse.

        Returns
        -------
        numpy.ndarray
            The inverse of `m` as a 3x3 complex numpy array.

        """
        m1, m2, m3, m4, m5, m6, m7, m8, m9 = m.flatten()
        determinant = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9  
        return np.array([[m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5],
                         [m6*m7-m4*m9, m1*m9-m3*m7, m3*m4-m1*m6],
                         [m4*m8-m5*m7, m2*m7-m1*m8, m1*m5-m2*m4]])/determinant

class ShowLayerWindow(QDialog):
    """A GUI window for displaying and editing the properties of a layer in a project.

    Parameters
    ----------
    layer : object
        The layer object to be displayed or acted upon in the window. (see :class:`SingleCrystalLayer`)
    message : str, optional
        Custom message to display in the window. The default is an empty string.
    parent : QWidget, optional
        The parent widget of this window. The default is None.
    debug : bool, optional
        Flag to enable or disable debug mode. The default is False.

    Attributes
    ----------
    layer : object
        The :class:`SingleCrystalLayer` instance that is displayed or being edited in this window.
    message : str
        The message displayed at the top of the window.
    layout : QVBoxLayout
        The main layout of the window.
    buttonBox : QDialogButtonBox
        The dialog box containing Ok and Cancel buttons.

    Methods
    -------
    getLayer()
        Return the :class:`SingleCrystalLayer` instance being edited or displayed.
    drawLayerWidget()
        Create a QWidget displaying all information about the layer.
    drawLayerWidgetLine1()
        Add widgets to display the first line of layer properties.
    drawLayerWidgetLine2()
        Add widgets to display the second line of layer properties related to the crystallographic direction (hkl).
    drawLayerWidgetLine3()
        Add widgets for displaying the third line of layer properties related to lab frame information.
    changeLabFrameInfo()
        Update the displayed lab frame information based on the layer's current state.
    on_film_thickness_sb_changed(value)
        Respond to changes in the film thickness spin box.
    on_azimuthal_angle_sb_changed(value)
        Respond to changes in the azimuthal angle spin box.
    on_hkl_sb_changed(value, hkorl)
        Respond to changes in the h, k, or l spin boxes.
    on_thickness_units_cb_activated(index)
        Respond to thickness unit selection changes.

    """

    def __init__(self, layer, message = "", parent=None, debug=False ):
        """Initialize a ShowLayerWindow instance.

        Parameters
        ----------
        layer : object
            The layer object to be displayed or acted upon in the window.
        message : str, optional
            Custom message to display in the window. The default is an empty string.
        parent : QWidget, optional
            The parent widget of this window. The default is None.
        debug : bool, optional
            Flag to enable or disable debug mode. The default is False.

        Notes
        -----
        This constructor initializes the ShowLayerWindow by setting up the UI elements
        including buttons and the layout. It also configures the debug mode according
        to the provided argument.

        """        
        super().__init__(parent)
        global debugger
        debugger = Debug(debug,"ShowLayerWindow")
        debugger.print("Start:: initialiser")
        # Set up the buttons of the button box
        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        # Initialize the layer
        self.layer = layer
        self.message = message
        # The dialog will have a vertical layout
        self.layout = QVBoxLayout(self)
        # Create the layer widget to display
        layerWidget = self.drawLayerWidget()
        # Add the widget showing the layer information
        self.layout.addWidget(layerWidget)
        # Add the button box
        self.layout.addWidget(self.buttonBox)
        debugger.print("Finished:: initialiser")

    def getLayer(self):
        """Return the edited layer.

        Parameters
        ----------
        None

        Returns
        -------
        self.layer : An instance of a single crystal layer

        """
        return self.layer

    def drawLayerWidget(self):
        """Create a layer widget showing all the information about the layer.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        debugger.print("drawLayerWidget")
        widget = QWidget()
        form = QFormLayout()
        label = QLabel("Layer type:")
        form.addRow(label,QLabel(self.message))
        label = QLabel("Material:")
        material = self.layer.getMaterial()
        materialName = material.getName()
        form.addRow(label,QLabel(materialName))
        label,layout = self.drawLayerWidgetLine1()
        form.addRow(label,layout)
        if self.layer.getMaterial().isTensor():
            label,layout = self.drawLayerWidgetLine2()
            form.addRow(label,layout)
            label,layout = self.drawLayerWidgetLine3()
            form.addRow(label,layout)
        widget.setLayout(form)
        return widget

    def drawLayerWidgetLine1(self):
        """Add the first line of the layer description.

        Shows the material name and the film thickness and units. 

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        hbox = QHBoxLayout()
        debugger.print("drawLayerWidgetLine1")
        # Define material thickness
        materialThickness = self.layer.getThickness()
        # Handle thickness
        film_thickness_sb = QDoubleSpinBox(self)
        film_thickness_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        film_thickness_sb.setToolTip("Define the thin film thickness in the defined thickness units")
        film_thickness_sb.setRange(0,100000)
        film_thickness_sb.setSingleStep(0.01)
        film_thickness_sb.setValue(materialThickness)
        film_thickness_sb.valueChanged.connect(self.on_film_thickness_sb_changed)
        # Handle thickness units
        thicknessUnit = self.layer.getThicknessUnit()
        thickness_units_cb = QComboBox(self)
        thickness_units_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thickness_units_cb.setToolTip("Set the units to be used for thickness; either nm, um, mm or cm")
        thickness_units_cb.addItems( ["nm","um","mm","cm"] )
        index = thickness_units_cb.findText(thicknessUnit, Qt.MatchFixedString)
        thickness_units_cb.setCurrentIndex(index)
        thickness_units_cb.activated.connect(self.on_thickness_units_cb_activated)
        thicknessLabel = QLabel("Thickness:")
        thicknessLabel.setToolTip("Define the depth of the thin crystal in the defined thickness units.")
        # Create the line of widgets
        hbox.addWidget(film_thickness_sb)
        hbox.addWidget(thickness_units_cb)
        return thicknessLabel,hbox

    def drawLayerWidgetLine2(self):
        """Add the second line of the layer description.

        Shows the hkl specification and the azimuthal angle

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        debugger.print("drawLayerWidgetLine2")
        hbox = QHBoxLayout()
        # define hkl
        h_sb = QSpinBox(self)
        h_sb.setToolTip("Define the h dimension of the unique direction")
        h_sb.setRange(-20,20)
        h_sb.setSingleStep(1)
        h_sb.setValue(self.layer.getHKL()[0])
        h_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,0))
        k_sb = QSpinBox(self)
        k_sb.setToolTip("Define the k dimension of the unique direction")
        k_sb.setRange(-20,20)
        k_sb.setSingleStep(1)
        k_sb.setValue(self.layer.getHKL()[1])
        k_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,1))
        l_sb = QSpinBox(self)
        l_sb.setToolTip("Define the l dimension of the unique direction")
        l_sb.setRange(-20,20)
        l_sb.setSingleStep(1)
        l_sb.setValue(self.layer.getHKL()[2])
        l_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,2))
        hklLabel = QLabel("hkl:")
        hklLabel.setToolTip("Define the crystal surface (hkl). Defines the unique direction in crystallographic units.")
        # define azimuthal angle
        azimuthal = self.layer.getAzimuthal()
        azimuthal_angle_sb = QDoubleSpinBox(self)
        azimuthal_angle_sb.setToolTip("Define the slab azimuthal angle (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal in the laboratory frame can be seen in the laboratory frame information below")
        azimuthal_angle_sb.setRange(-180,360)
        azimuthal_angle_sb.setSingleStep(10)
        azimuthal_angle_sb.setValue(azimuthal)
        azimuthal_angle_sb.valueChanged.connect(self.on_azimuthal_angle_sb_changed)
        azimuthalLabel = QLabel("Azimuthal:")
        # Create the line of widgets
        hbox.addWidget(h_sb)
        hbox.addWidget(k_sb)
        hbox.addWidget(l_sb)
        hbox.addWidget(azimuthalLabel)
        hbox.addWidget(azimuthal_angle_sb)
        return hklLabel,hbox

    def drawLayerWidgetLine3(self):
        """Add the third line of the layer description.

        Shows the relationship between crystal and labframe coordinate systems

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        debugger.print("drawLayerWidgetLine3")
        hbox = QHBoxLayout()
        label = QLabel("Lab frame\ninformation", self)
        label.setToolTip("The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y")
        self.labframe_w = QListWidget(self)
        fm = self.labframe_w.fontMetrics()
        h = fm.ascent() + fm.descent()
        self.labframe_w.setMaximumHeight(6*h)
        self.labframe_w.setToolTip("The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y")
        self.changeLabFrameInfo()
        hbox.addWidget(self.labframe_w)
        return label,hbox

    def changeLabFrameInfo(self):
        """Update the lab frame information in the labframe widget.

        Shows a, b and c crystal axis in the labframe.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        a,b,c = self.layer.getLabFrame()
        self.labframe_w.clear()
        self.labframe_w.addItem(f"crystal a-axis in lab frame: {a[0]: 3.5f}, {a[1]: 3.5f}, {a[2]: 3.5f}" )
        self.labframe_w.addItem(f"crystal b-axis in lab frame: {b[0]: 3.5f}, {b[1]: 3.5f}, {b[2]: 3.5f}" )
        self.labframe_w.addItem(f"crystal c-axis in lab frame: {c[0]: 3.5f}, {c[1]: 3.5f}, {c[2]: 3.5f}" )
        return

    def on_film_thickness_sb_changed(self,value):
        """Handle film thickness spin box change.

        Parameters
        ----------
        value : float
            The new thickness

        Returns
        -------
        None

        """
        debugger.print("on_film_thickness_sb_changed", value)
        self.layer.setThickness(value)
        return

    def on_azimuthal_angle_sb_changed(self,value):
        """Handle azimuthal spin box change.

        Parameters
        ----------
        value : float
            The new azimuthal angle in degrees

        Returns
        -------
        None

        """
        debugger.print("on_azimuthal_angl_sb_changed", value)
        self.layer.setAzimuthal(value)
        self.layer.changeLabFrameInfo()
        self.changeLabFrameInfo()
        return

    def on_hkl_sb_changed(self,value,hkorl):
        """Handle h, k, or l change.

        Parameters
        ----------
        value : int
           The value of h, k or l
        hkorl : int
            0, 1, or 2 for h, k or l respectively.

        Returns
        -------
        None

        """
        debugger.print("on_h_sb_changed", value)
        hkl = self.layer.getHKL()
        hkl[hkorl] = value
        self.layer.setHKL(hkl)
        self.layer.changeLabFrameInfo()
        self.changeLabFrameInfo()
        return

    def on_thickness_units_cb_activated(self, index):
        """Handle activation of thickness units combobox.

        Parameters
        ----------
        index : int
            The index of the selected item in the combobox.
            One of 'nm', 'um', 'mm', 'cm'

        Returns
        -------
        None

        """        
        debugger.print("Start:: on_thickness_units_cb_activated",index)
        units = ["nm","um","mm","cm"]
        unit = units[index]
        self.layer.setThicknessUnit(unit)
        return
