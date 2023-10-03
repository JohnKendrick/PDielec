# -*- coding: utf8 -*-
import numpy as np
from PyQt5.QtWidgets            import QPushButton, QWidget, QFrame, QDialog
from PyQt5.QtWidgets            import QComboBox, QLabel, QLineEdit, QListWidget
from PyQt5.QtWidgets            import QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets            import QSpinBox,QDoubleSpinBox
from PyQt5.QtWidgets            import QSizePolicy
from PyQt5.QtCore               import QCoreApplication, Qt
from PDielec.Utilities          import Debug
import PDielec.Calculator as Calculator

class Layer():
    def __init__(self,material,hkl=None,azimuthal=0.0,thickness=0.0,thicknessUnit='nm',dielectricFlag=False):
        '''A Layer class to handle layers
           material    is an instance of Material object, the material object has the following;
                       - name
                       - density
                       - permittivity
                       - unit cell
           hkl         is a list of 3 integers defining the plane of interest
           azimuthal   is the azimuthal angle of rotation of the crystal about z
           thickness   is the thickness of the layer in the specfied thickness units
           thicknessUnits can be either nm, um, mm or cm 
           dielectricFlag this is true if the layer material is the dielectric being studied'''
        self.material = material
        self.hkl = hkl
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
        print('---------------- ')
        print('Material       : ', self.material)
        print('HKL            :      ', self.hkl)
        print('Azimuthal      :', self.azimuthal)
        print('Thickness      :', self.thickness)
        print('Thickness unit :', self.thicknessUnit)
        print('Dielectric flag:', self.dielectricFlag)

    def setAzimuthal(self, angle):
        '''Set the azimuthal angle'''
        self.azimuthal = angle
        self.calculate_euler_matrix()
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
        # debugger.print(self.settings['Legend'],'Start:: calculate_euler_matrix')
        # Get plane specification
        hkl = self.hkl
        sum2 = hkl[0]*hkl[0] + hkl[1]*hkl[1] + hkl[2]*hkl[2]
        if sum2 < 1:
            # debugger.print(self.settings['Legend'],'Finished:: calculate_euler_matrix')
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
        self.labframe = np.matmul(self.euler,self.material.cell.lattice).T
        normal_to_plane_lab = np.matmul(self.euler,plane[z])
        if normal_to_plane_lab[2] < 0.9999 and normal_to_plane_lab[2] > -0.9999:
            print('Error in Euler rotations - surface normal is not along Z-axis', normal_to_plane_lab)
            exit()
        # debugger.print(self.settings['Legend'],'Finished:: calculate_euler_matrix')
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

class EditLayerWindow(QDialog):
    def __init__(self, layers,database,parent=None,debug=False ):
        super(EditLayerWindow,self).__init__(parent)
        global debugger
        debugger = Debug(debug,'EditLayerWindow')
        debugger.print('Start:: initialiser')
        self.DataBase = database
        self.layers = layers
        self.newMaterialName = None
        self.materialNames = None
        mainLayout = QVBoxLayout()
        self.setLayout(mainLayout)
        self.redraw(clearLayout = False)
        debugger.print('Finished:: initialiser')

    def clearLayout(self,layout):
        '''Clear the existing layout of all widgets and layouts'''
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clearLayout(item.layout())

    def redraw(self,clearLayout=False):
        '''Redraw the edit layers window'''
        if clearLayout:
            self.clearLayout(self.layout())
            QCoreApplication.processEvents()
        form = QFormLayout()
        # Add each layer to the window
        for i,layer in enumerate(self.layers):
            layout = self.drawLayerWidget(i,layer)
            form.addRow(layout)
            line = QFrame()
            line.setFrameShape(QFrame.HLine)
            form.addRow(line)
        layout = self.drawAddLayerWidget()
        form.addRow(layout)
        # Add the form to the existing layout
        self.layout().addLayout(form)
        self.update()
        for i in range(0,10):
            QCoreApplication.processEvents()
        self.resize(self.minimumSizeHint())
        return

    def drawAddLayerWidget(self):
        '''Add widget with the option to add more layers'''
        hbox = QHBoxLayout()
        material_cb = QComboBox(self)
        material_cb.setToolTip('Choose the material to be added as a new layer.')
        self.materialNames = self.DataBase.getSheetNames()
        if self.newMaterialName not in self.materialNames:
            self.newMaterialName = self.materialNames[0]
        material_cb.clear()
        material_cb.addItems(self.materialNames)
        index = material_cb.findText(self.newMaterialName, Qt.MatchFixedString)
        material_cb.setCurrentIndex(index)
        material_cb.activated.connect(self.on_material_cb_activated)
        add_button = QPushButton('Add layer of.....')
        add_button.clicked.connect(lambda x: self.addAnotherLayer(x))
        close_button = QPushButton('Close window')
        close_button.clicked.connect(lambda x: self.hideWindow(x))
        hbox.addWidget(add_button)
        hbox.addWidget(material_cb)
        hbox.addWidget(close_button)
        return hbox

    def on_material_cb_activated(self,index):
        self.newMaterialName = self.materialNames[index]
        return
    
    def hideWindow(self,x):
        '''Hide the window'''
        self.hide()

    def addAnotherLayer(self,x):
        '''Add another layer to the window'''
        print('Ready to add another layer')
        hkl = [0,0,0]
        newMaterial = self.DataBase.getMaterial(self.newMaterialName)
        new_layer = Layer(newMaterial,hkl=hkl,azimuthal=0.0,thickness=1.0,thicknessUnit='um')
        self.layers.append(new_layer)
        self.redraw(clearLayout = True)
        return

    def drawLayerWidget(self,sequenceNumber,layer):
        '''Draw a layer widget'''
        form = QFormLayout()
        layout = self.drawLayerWidgetLine1(sequenceNumber,layer)
        form.addRow(layout)
        if layer.getMaterial().isTensor():
            layout = self.drawLayerWidgetLine2(sequenceNumber,layer)
            form.addRow(layout)
            layout = self.drawLayerWidgetLine3(sequenceNumber,layer)
        form.addRow(layout)
        return form

    def drawLayerWidgetLine1(self,sequenceNumber,layer):
        hbox = QHBoxLayout()
        # Add sequence number
        sequenceLabel = QLabel('No. '+str(sequenceNumber),self)
        # Define material name
        material = layer.getMaterial()
        materialName = material.getName()
        materialLabel = QLabel(materialName,self)
        materialThickness = layer.getThickness()
        # Handle thickness units
        thicknessUnit = layer.getThicknessUnit()
        film_thickness_sb = QDoubleSpinBox(self)
        film_thickness_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        film_thickness_sb.setToolTip('Define the thin film thickness in the defined thickness units')
        film_thickness_sb.setRange(0,100000)
        film_thickness_sb.setSingleStep(1)
        film_thickness_sb.setValue(materialThickness)
        film_thickness_sb.valueChanged.connect(lambda x: self.on_film_thickness_sb_changed(x,layer))
        thickness_unit_cb = QComboBox(self)
        thickness_unit_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thickness_unit_cb.setToolTip('Set the units to be used for thickness; either nm, um, mm or cm')
        thickness_unit_cb.addItems( ['nm','um','mm','cm'] )
        index = thickness_unit_cb.findText(thicknessUnit, Qt.MatchFixedString)
        thickness_unit_cb.setCurrentIndex(index)
        thickness_unit_cb.activated.connect(lambda x: self.on_thickness_units_cb_activated(x, layer))
        thicknessLabel = QLabel('Film thickness',self)
        thicknessLabel.setToolTip('Define the depth of the thin crystal in the defined thickness units.')
        # Add move up and down buttons
        move_up_button = QPushButton('Up')
        move_up_button.clicked.connect(lambda x: self.moveLayerUp(x, sequenceNumber))
        move_down_button = QPushButton('Down')
        move_down_button.clicked.connect(lambda x: self.moveLayerDown(x, sequenceNumber))
        # Add delete button
        delete_button = QPushButton('Del')
        delete_button.setToolTip('Delete this layer')
        delete_button.clicked.connect(lambda x: self.deleteLayer(x, sequenceNumber,layer))
        # Create the line of widgets
        hbox.addWidget(sequenceLabel)
        hbox.addWidget(materialLabel)
        hbox.addWidget(thicknessLabel)
        hbox.addWidget(film_thickness_sb)
        hbox.addWidget(thickness_unit_cb)
        hbox.addWidget(move_up_button)
        hbox.addWidget(move_down_button)
        hbox.addWidget(delete_button)
        return hbox

    def drawLayerWidgetLine2(self,sequenceNumber,layer):
        hbox = QHBoxLayout()
        # define hkl
        h_sb = QSpinBox(self)
        h_sb.setToolTip('Define the h dimension of the unique direction')
        h_sb.setRange(-20,20)
        h_sb.setSingleStep(1)
        h_sb.setValue(layer.getHKL()[0])
        h_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,0,layer))
        k_sb = QSpinBox(self)
        k_sb.setToolTip('Define the k dimension of the unique direction')
        k_sb.setRange(-20,20)
        k_sb.setSingleStep(1)
        k_sb.setValue(layer.getHKL()[1])
        k_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,1,layer))
        l_sb = QSpinBox(self)
        l_sb.setToolTip('Define the l dimension of the unique direction')
        l_sb.setRange(-20,20)
        l_sb.setSingleStep(1)
        l_sb.setValue(layer.getHKL()[2])
        l_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,2,layer))
        hklLabel = QLabel('hkl:')
        hklLabel.setToolTip('Define the crystal surface (hkl). Defines the unique direction in crystallographic units.')
        # define azimuthal angle
        azimuthal = layer.getAzimuthal()
        azimuthal_angle_sb = QDoubleSpinBox(self)
        azimuthal_angle_sb.setToolTip('Define the slab azimuthal angle (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal in the laboratory frame can be seen in the laboratory frame information below')
        azimuthal_angle_sb.setRange(-180,360)
        azimuthal_angle_sb.setSingleStep(10)
        azimuthal_angle_sb.setValue(azimuthal)
        azimuthal_angle_sb.valueChanged.connect(lambda x: self.on_azimuthal_angle_sb_changed(x,layer))
        azimuthalLabel = QLabel('Azimuthal:')
        # Create the line of widgets
        hbox.addWidget(hklLabel)
        hbox.addWidget(h_sb)
        hbox.addWidget(k_sb)
        hbox.addWidget(l_sb)
        hbox.addWidget(azimuthalLabel)
        hbox.addWidget(azimuthal_angle_sb)
        return hbox

    def drawLayerWidgetLine3(self,sequenceNumber,layer):
        '''Add the lab frame box'''
        hbox = QHBoxLayout()
        label = QLabel('Lab frame information', self)
        label.setToolTip('The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y')
        layer.labframe_w = QListWidget(self)
        fm = layer.labframe_w.fontMetrics()
        h = fm.ascent() + fm.descent()
        layer.labframe_w.setMaximumHeight(6*h)
        layer.labframe_w.setToolTip('The normal to the surface defines the Z-axis in the  lab frame\nThe incident and reflected light lie in the XZ plane\nThe p-polarization is direction lies in the XZ plane, s-polarisation is parallel to Y')
        a,b,c = layer.getLabFrame()
        layer.labframe_w.clear()
        layer.labframe_w.addItem('crystal a-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(a[0],a[1],a[2]) )
        layer.labframe_w.addItem('crystal b-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(b[0],b[1],b[2]) )
        layer.labframe_w.addItem('crystal c-axis in lab frame: {: 3.5f}, {: 3.5f}, {: 3.5f}'.format(c[0],c[1],c[2]) )
        hbox.addWidget(label)
        hbox.addWidget(layer.labframe_w)
        return hbox

    def getLayers(self):
        '''Return the list of layers'''
        return self.layers

    def deleteLayer(self,i,sequenceNumber,layer):
        '''Handle a delete layer button press'''
        # Only delete a layer if it doesn't have the dielectricFlag set
        if not layer.isDielectric():
            del self.layers[sequenceNumber]
        self.redraw(clearLayout = True)

    def moveLayerUp(self,x,old):
        '''Move a layer up (sequence number gets smaller by 1)'''
        print('MoveLayerup:',x,old)
        if old < 1:
            return
        new = old - 1
        item = self.layers[old]
        self.layers.pop(old)
        self.layers.insert(new, item)
        self.redraw(clearLayout = True)
        return

    def moveLayerDown(self,x,old):
        '''Move a layer down (sequence number gets larger by 1)'''
        print('MoveLayerup:',x,old)
        last = len(self.layers) - 1
        if old >= last:
            return
        new = old + 1
        item = self.layers[old]
        self.layers.pop(old)
        self.layers.insert(new, item)
        self.redraw(clearLayout = True)
        return

    def on_film_thickness_sb_changed(self,value,layer):
        '''Handle film thickness spin box change'''
        debugger.print('on_film_thickness_sb_changed', value, layer.getName())
        layer.setThickness(value)
        return

    def on_azimuthal_angle_sb_changed(self,value,layer):
        '''Handle azimuthal spin box change'''
        debugger.print('on_azimuthal_angl_sb_changed', value, layer.getName())
        layer.setAzimuthal(value)
        layer.changeLabFrameInfo()
        return

    def on_hkl_sb_changed(self,value,hkorl,layer):
        '''Handle h, k or l change'''
        debugger.print('on_h_sb_changed', value)
        hkl = layer.getHKL()
        hkl[hkorl] = value
        layer.setHKL(hkl)
        layer.changeLabFrameInfo()
        return

    def on_thickness_units_cb_activated(self, index, layer):
        debugger.print('Start:: on_thickness_units_cb_activated',index,layer.getName())
        units = ['nm','um','mm','cm']
        unit = units[index]
        layer.setThicknessUnit(unit)
        return
