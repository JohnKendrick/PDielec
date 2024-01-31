# -*- coding: utf8 -*-
import numpy as np
from PyQt5.QtWidgets            import QPushButton, QWidget, QFrame, QDialog
from PyQt5.QtWidgets            import QComboBox, QLabel, QLineEdit, QListWidget
from PyQt5.QtWidgets            import QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets            import QSpinBox,QDoubleSpinBox, QCheckBox
from PyQt5.QtWidgets            import QSizePolicy, QDialogButtonBox
from PyQt5.QtCore               import QCoreApplication, Qt
from PDielec.Utilities          import Debug
import PDielec.Calculator as Calculator

class Layer():
    def __init__(self,material,hkl=None,azimuthal=0.0,thickness=0.0,thicknessUnit='nm',incoherentOption='Coherent',dielectricFlag=False):
        '''A Layer class to handle layers
           material          is an instance of Material object, the material object has the following;
                             - name
                             - density
                             - permittivity
                             - unit cell
           hkl               is a list of 3 integers defining the plane of interest
           azimuthal         is the azimuthal angle of rotation of the crystal about z
           thickness         is the thickness of the layer in the specfied thickness units
           thicknessUnits    can be either nm, um, mm or cm 
           incoherentOption  can be 'Coherent', 'Incoherent intensity', 'Incoherent phase cancelling' or 'Thick slab'
           dielectricFlag    this is true if the layer material is the dielectric being studied'''
        self.material = material
        self.hkl = hkl
        if material.isTensor() and hkl is None:
            hkl = [0,0,1]
        else:
            hkl = [0,0,0]
        self.incoherentOption = incoherentOption
        self.azimuthal = azimuthal
        self.thickness = thickness
        self.thicknessUnit = thicknessUnit
        self.dielectricFlag = dielectricFlag
        self.euler = np.zeros((3,3),dtype=np.longdouble)
        self.euler_inverse = np.zeros((3,3),dtype=np.longdouble)
        self.labframe_w = None
        self.labframe = None
        self.calculate_euler_matrix()

    def print(self):
        '''Print the main attributes of the layer'''
        print('------------------ ')
        print('Material         : ', self.material)
        print('HKL              :      ', self.hkl)
        print('Azimuthal        :', self.azimuthal)
        print('Thickness        :', self.thickness)
        print('Thickness unit   :', self.thicknessUnit)
        print('Dielectric flag  :', self.dielectricFlag)
        print('Incoherent option:', self.incoherentOption)

    def setAzimuthal(self, angle):
        '''Set the azimuthal angle of the layer and calculate the Euler matrix'''
        self.azimuthal = angle
        self.calculate_euler_matrix()
        return 

    def isCoherent(self):
        '''Returns True if this is a coherent layer, False otherwise'''
        result = False
        if incoherentOption == 'Coherent':
            result = True
        return result

    def getIncoherentOption(self):
        return self.incoherentOption

    def setIncoherentOption(self,option):
        '''Set coherent to False'''
        self.incoherentOption = option
        return 

    def getAzimuthal(self):
        '''Get the azimuthal angle'''
        return self.azimuthal

    def setThickness(self, thickness):
        '''Set the thickness'''
        self.thickness = thickness
        return 

    def getThickness(self):
        '''Get the thickness'''
        return self.thickness

    def getThicknessInMetres(self):
        '''Get the thickness'''
        thicknessUnits = {'nm':1.0E-9, 'um':1.0E-6, 'mm':1.0E-3, 'cm':1.0E-2}
        tom = thicknessUnits[self.thicknessUnit]
        return tom*self.thickness

    def setThicknessUnit(self, thicknessUnit):
        '''Set the thicknessUnit'''
        self.thicknessUnit = thicknessUnit
        return 

    def getThicknessUnit(self):
        '''Get the thicknessUnit'''
        return self.thicknessUnit

    def setHKL(self, hkl):
        '''Set the hkl'''
        self.hkl = hkl
        self.calculate_euler_matrix()
        return 

    def getHKL(self):
        '''Get the hkl'''
        return self.hkl

    def getName(self):
        '''Return the material name'''
        return self.material.getName()

    def getMaterial(self):
        '''Return the material'''
        return self.material

    def getPermittivityFunction(self):
        '''Return the permittivity tensor for the surface and azimuthal angle specified'''
        return self.calculateLabFrameEpsilon

    def isTensor(self):
        '''Return true if the material of the layer is a tensor material'''
        return self.getMaterial().isTensor()

    def isScalar(self):
        '''Return true if the material of the layer is a scalar material'''
        return self.getMaterial().isScalar()

    def isDielectric(self):
        '''Return true if the material of the layer is the dielectric material'''
        return self.dielectricFlag

    def calculate_euler_matrix(self):
        '''Calculate the Euler angles for the crystal to lab transformation'''
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
            print('Error in Euler rotations - surface normal is not along Z-axis', normal_to_plane_lab)
            exit()
        return 

    def azimuthalRotationMatrix(self,angle):
        '''Calculate a rotation matrix for the azimuthal angle (in degrees)'''
        angle = np.radians(angle)
        matrix = np.eye(3)
        # in numpy the first index is the column, the second is the row
        matrix[0,0] = np.cos(angle)
        matrix[1,1] = np.cos(angle)
        matrix[0,1] = -np.sin(angle)
        matrix[1,0] = +np.sin(angle)
        return matrix

    def changeLabFrameInfo(self):
        '''Change the Lab Frame Info panel'''
        if self.labframe_w is None:
            return
        a,b,c = self.getLabFrame()
        self.labframe_w.clear()
        self.labframe_w.addItem('crystal a-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(a[0],a[1],a[2]) )
        self.labframe_w.addItem('crystal b-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(b[0],b[1],b[2]) )
        self.labframe_w.addItem('crystal c-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(c[0],c[1],c[2]) )
        return

    def getLabFrame(self):
        '''Return the a,b and c cell axes in the laboratory frame coordinates'''
        return self.labframe

    def calculateLabFrameEpsilon(self, v):
        '''Get the permittivity for this material at frequency v and transform it to the laboratory frame'''
        epsilon_xstal = self.material.getPermittivityFunction()(v)
        if self.material.isScalar():
            self.labFrameEpsilon = epsilon_xstal * np.eye(3)
        else:
            self.labFrameEpsilon = np.matmul(self.euler, np.matmul(epsilon_xstal,self.euler_inverse))
        return self.labFrameEpsilon

    def invert(self,m):
        """Calculate the inverse of m

        Parameters
        ----------
        m : complex 3x3 array

        Returns
        -------
        the inverse as a 3x3 complex np array
        """
        m1, m2, m3, m4, m5, m6, m7, m8, m9 = m.flatten()
        determinant = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9  
        return np.array([[m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5],
                         [m6*m7-m4*m9, m1*m9-m3*m7, m3*m4-m1*m6],
                         [m4*m8-m5*m7, m2*m7-m1*m8, m1*m5-m2*m4]])/determinant

class ShowLayerWindow(QDialog):
    def __init__(self, layer, parent=None, debug=False ):
        super(ShowLayerWindow,self).__init__(parent)
        global debugger
        debugger = Debug(debug,'ShowLayerWindow')
        debugger.print('Start:: initialiser')
        # Set up the buttons of the button box
        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        # Initialize the layer
        self.layer = layer
        # The dialog will have a vertical layout
        self.layout = QVBoxLayout(self)
        # Create the layer widget to display
        layerWidget = self.drawLayerWidget()
        # Add the widget showing the layer information
        self.layout.addWidget(layerWidget)
        # Add the button box
        self.layout.addWidget(self.buttonBox)
        debugger.print('Finished:: initialiser')

    def getLayer(self):
        '''Return the edited layer'''
        return self.layer

    def drawLayerWidget(self):
        '''Creat a layer widget showing all the information about the layer'''
        debugger.print('drawLayerWidget')
        widget = QWidget()
        form = QFormLayout()
        layout = self.drawLayerWidgetLine1()
        form.addRow(layout)
        if self.layer.getMaterial().isTensor():
            layout = self.drawLayerWidgetLine2()
            form.addRow(layout)
            layout = self.drawLayerWidgetLine3()
        form.addRow(layout)
        widget.setLayout(form)
        return widget

    def drawLayerWidgetLine1(self):
        '''Add the first line of the layer description'''
        hbox = QHBoxLayout()
        debugger.print('drawLayerWidgetLine1')
        # Define material name
        material = self.layer.getMaterial()
        materialName = material.getName()
        materialLabel = QLabel(materialName,self)
        materialThickness = self.layer.getThickness()
        # Handle thickness
        film_thickness_sb = QDoubleSpinBox(self)
        film_thickness_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        film_thickness_sb.setToolTip('Define the thin film thickness in the defined thickness units')
        film_thickness_sb.setRange(0,100000)
        film_thickness_sb.setSingleStep(0.01)
        film_thickness_sb.setValue(materialThickness)
        film_thickness_sb.valueChanged.connect(self.on_film_thickness_sb_changed)
        # Handle thickness units
        thicknessUnit = self.layer.getThicknessUnit()
        thickness_units_cb = QComboBox(self)
        thickness_units_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thickness_units_cb.setToolTip('Set the units to be used for thickness; either nm, um, mm or cm')
        thickness_units_cb.addItems( ['nm','um','mm','cm'] )
        index = thickness_units_cb.findText(thicknessUnit, Qt.MatchFixedString)
        thickness_units_cb.setCurrentIndex(index)
        thickness_units_cb.activated.connect(self.on_thickness_units_cb_activated)
        thicknessLabel = QLabel('Thickness',self)
        thicknessLabel.setToolTip('Define the depth of the thin crystal in the defined thickness units.')
        # Create the line of widgets
        hbox.addWidget(materialLabel)
        hbox.addWidget(thicknessLabel)
        hbox.addWidget(film_thickness_sb)
        hbox.addWidget(thickness_units_cb)
        return hbox

    def drawLayerWidgetLine2(self):
        '''Add the second line of the layer description'''
        debugger.print('drawLayerWidgetLine2')
        hbox = QHBoxLayout()
        # define hkl
        h_sb = QSpinBox(self)
        h_sb.setToolTip('Define the h dimension of the unique direction')
        h_sb.setRange(-20,20)
        h_sb.setSingleStep(1)
        h_sb.setValue(self.layer.getHKL()[0])
        h_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,0))
        k_sb = QSpinBox(self)
        k_sb.setToolTip('Define the k dimension of the unique direction')
        k_sb.setRange(-20,20)
        k_sb.setSingleStep(1)
        k_sb.setValue(self.layer.getHKL()[1])
        k_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,1))
        l_sb = QSpinBox(self)
        l_sb.setToolTip('Define the l dimension of the unique direction')
        l_sb.setRange(-20,20)
        l_sb.setSingleStep(1)
        l_sb.setValue(self.layer.getHKL()[2])
        l_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,2))
        hklLabel = QLabel('hkl:')
        hklLabel.setToolTip('Define the crystal surface (hkl). Defines the unique direction in crystallographic units.')
        # define azimuthal angle
        azimuthal = self.layer.getAzimuthal()
        azimuthal_angle_sb = QDoubleSpinBox(self)
        azimuthal_angle_sb.setToolTip('Define the slab azimuthal angle (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal in the laboratory frame can be seen in the laboratory frame information below')
        azimuthal_angle_sb.setRange(-180,360)
        azimuthal_angle_sb.setSingleStep(10)
        azimuthal_angle_sb.setValue(azimuthal)
        azimuthal_angle_sb.valueChanged.connect(self.on_azimuthal_angle_sb_changed)
        azimuthalLabel = QLabel('Azimuthal:')
        # Create a checkbox for incoherence
        incoherence_cb = QCheckBox(self)
        incoherence_cb.setToolTip('Treat layer as incoherent using phase cancellation')
        incoherence_cb.setText('Incoherent?')
        incoherence_cb.setLayoutDirection(Qt.RightToLeft)
        if self.layer.isCoherent():
            incoherence_cb.setCheckState(Qt.Unchecked)
        else:
            incoherence_cb.setCheckState(Qt.Checked)
        incoherence_cb.stateChanged.connect(self.on_incoherence_cb_changed)
        # Create the line of widgets
        hbox.addWidget(hklLabel)
        hbox.addWidget(h_sb)
        hbox.addWidget(k_sb)
        hbox.addWidget(l_sb)
        hbox.addWidget(azimuthalLabel)
        hbox.addWidget(azimuthal_angle_sb)
        hbox.addWidget(incoherence_cb)
        return hbox

    def drawLayerWidgetLine3(self):
        '''Add the third line of the layer description'''
        debugger.print('drawLayerWidgetLine3')
        hbox = QHBoxLayout()
        label = QLabel('Lab frame\ninformation', self)
        label.setToolTip('The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y')
        self.labframe_w = QListWidget(self)
        fm = self.labframe_w.fontMetrics()
        h = fm.ascent() + fm.descent()
        self.labframe_w.setMaximumHeight(6*h)
        self.labframe_w.setToolTip('The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y')
        self.changeLabFrameInfo()
        hbox.addWidget(label)
        hbox.addWidget(self.labframe_w)
        return hbox

    def changeLabFrameInfo(self):
        '''Report the lab frame information in the labframe widget'''
        a,b,c = self.layer.getLabFrame()
        self.labframe_w.clear()
        self.labframe_w.addItem('crystal a-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(a[0],a[1],a[2]) )
        self.labframe_w.addItem('crystal b-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(b[0],b[1],b[2]) )
        self.labframe_w.addItem('crystal c-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(c[0],c[1],c[2]) )
        return

    def on_incoherence_cb_changed(self,value):
        '''The incoherence checkbox has changed'''
        debugger.print('on_incoherence_cb_changed', value)
        state = self.layer.isCoherent()
        if state:
            self.layer.setIncoherent()
        else:
            self.layer.setCoherent()
        return

    def on_film_thickness_sb_changed(self,value):
        '''Handle film thickness spin box change'''
        debugger.print('on_film_thickness_sb_changed', value)
        self.layer.setThickness(value)
        return

    def on_azimuthal_angle_sb_changed(self,value):
        '''Handle azimuthal spin box change'''
        debugger.print('on_azimuthal_angl_sb_changed', value)
        self.layer.setAzimuthal(value)
        self.layer.changeLabFrameInfo()
        self.changeLabFrameInfo()
        return

    def on_hkl_sb_changed(self,value,hkorl):
        '''Handle h, k or l change'''
        debugger.print('on_h_sb_changed', value)
        hkl = self.layer.getHKL()
        hkl[hkorl] = value
        self.layer.setHKL(hkl)
        self.layer.changeLabFrameInfo()
        self.changeLabFrameInfo()
        return

    def on_thickness_units_cb_activated(self, index):
        debugger.print('Start:: on_thickness_units_cb_activated',index)
        units = ['nm','um','mm','cm']
        unit = units[index]
        self.layer.setThicknessUnit(unit)
        return
