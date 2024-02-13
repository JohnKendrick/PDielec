# -*- coding: utf8 -*-
import os.path
import os
import numpy as np
import copy
import PDielec.Calculator         as Calculator
import PDielec.DielectricFunction as DielectricFunction
import PDielec.Materials          as Materials
import PDielec.GTMcore            as GTM
from PyQt5.QtWidgets             import QToolBar, QHeaderView, QAction
from PyQt5.QtWidgets             import QPushButton, QWidget, QFrame, QToolButton
from PyQt5.QtWidgets             import QSpacerItem, QCheckBox
from PyQt5.QtWidgets             import QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets             import QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets             import QSpinBox,QDoubleSpinBox
from PyQt5.QtWidgets             import QSizePolicy, QApplication, QStyle
from PyQt5.QtWidgets             import QTableWidget, QTableWidgetItem
from PyQt5.QtCore                import QCoreApplication, Qt, QSize
from PyQt5.QtGui                 import QIcon
from PDielec.Constants           import wavenumber, PI, avogadro_si, angstrom, speed_light_si
from PDielec.Utilities           import Debug
from PDielec.GUI.ScenarioTab     import ScenarioTab
from PDielec.Materials           import MaterialsDataBase
from PDielec.GUI.EditLayerWindow import ShowLayerWindow
from PDielec.GUI.EditLayerWindow import Layer
from functools                   import partial
from scipy                       import signal

thickness_conversion_factors = {'ang':1.0E-10, 'nm':1.0E-9, 'um':1.0E-6, 'mm':1.0E-3, 'cm':1.0E-2}
thickness_units = list(thickness_conversion_factors.keys())
incoherentOptions = ['Coherent','Incoherent (intensity)','Incoherent (phase cancelling)','Incoherent (non-reflective)'] 
gtmMethods ={'Coherent':GTM.CoherentLayer,
             'Incoherent (intensity)':GTM.IncoherentIntensityLayer,
             'Incoherent (phase cancelling)':GTM.IncoherentPhaseLayer,
             'Incoherent (non-reflective)':GTM.IncoherentThickLayer} 

def trace3x3(m):
    if isinstance(m,np.ndarray):
        return (m[0,0] + m[1,1] + m[2,2])/3.0
    else:
        return m

def solve_single_crystal_equations( 
        layers                        ,
        mode                          ,
        theta                         ,
        phi                           ,
        psi                           ,
        angleOfIncidence              ,
        sliceThickness                ,
        exponent_threshold            ,
        v                             ,
        ):
    """
        This is a parallel call to the single crystal equation solver, system is a GTM system   
        layers                        a list of material layers (their permittivity functions)
        mode                          'Transfer matrix' or 'Scattering matrix'
        theta                         the theta angle of the sla,
        phi                           the angle of the slab
        psi                           the psi angle of the slab
        angleOfIncidence              the angle of incidence
        sliceThickness                a thickness in m, used to subdivide thicker films
                                      if zero then the full fat film is used
        exponent_threshold            largest exponent allowed in calculation of the propagation matrix
        v                             the frequency of the light in cm-1

    """
    # Create superstrate from the first layer
    superstrate = layers[0]
    superstrateDepth = superstrate.getThicknessInMetres()
    superstrateDielectricFunction = superstrate.getPermittivityFunction()
    superstrate      = GTM.CoherentLayer(thickness=superstrateDepth,epsilon=superstrateDielectricFunction,exponent_threshold=exponent_threshold)
    # Create substrate from the last layer
    substrate = layers[-1]
    substrateDepth = substrate.getThicknessInMetres()
    substrateDielectricFunction = substrate.getPermittivityFunction()
    substrate        = GTM.CoherentLayer(thickness=substrateDepth,  epsilon=substrateDielectricFunction,exponent_threshold=exponent_threshold)
    selectedLayers = layers[1:-1]
    gtmLayers = []
    # Create layers from all the layers between first and last
    for layer in selectedLayers:
        permittivityFunction = layer.getPermittivityFunction()
        depth = layer.getThicknessInMetres()
        incoherentOption = layer.getIncoherentOption()
        if sliceThickness != 0 and depth > sliceThickness:
            # Slice if the slicing is chosen and if the layer thickness is larger than the threshold
            no_of_layers = int(depth / sliceThickness) + 1
            sliced_depth = depth / no_of_layers
            for i in range(no_of_layers):
                gtmLayers.append(gtmMethods[incoherentOption](thickness=sliced_depth, epsilon=permittivityFunction, exponent_threshold=exponent_threshold))
        else:
            # No slicing is necessary
            gtmLayers.append(gtmMethods[incoherentOption](thickness=depth, epsilon=permittivityFunction, exponent_threshold=exponent_threshold))
    # Creat the system with the layers 
    system = GTM.System(substrate=substrate, superstrate=superstrate, layers=gtmLayers)
    # Rotate the dielectric constants to the laboratory frame
    # This is a global rotation of all the layers.
    system.substrate.set_euler(theta, phi, psi)
    system.superstrate.set_euler(theta, phi, psi)
    for layer in system.layers:
        layer.set_euler(theta, phi, psi)
    # 
    # convert cm-1 to frequency
    #
    freq = v * speed_light_si * 1e2
    system.initialize_sys(freq)
    zeta_sys = np.sin(angleOfIncidence)*np.sqrt(system.superstrate.epsilon[0,0])
    Sys_Gamma = system.calculate_GammaStar(freq, zeta_sys)
    r, R, t, T = system.calculate_r_t(zeta_sys)
    if len(system.layers) > 0:
        epsilon = system.layers[0].epsilon
    else:
        epsilon = system.substrate.epsilon
    errors,largest_exponent = system.overflowErrors()
    return v,r,R,t,T,epsilon,errors,largest_exponent

class SingleCrystalScenarioTab(ScenarioTab):
    def __init__(self, parent, debug=False ):
        ScenarioTab.__init__(self,parent)
        global debugger
        debugger = Debug(debug,'SingleCrystalScenarioTab:')
        debugger.print('Start:: initialiser')
        self.refreshRequired = True
        self.calculationRequired = True
        self.scenarioType = 'Single crystal'
        self.noCalculationsRequired = 1
        self.settings['Scenario type'] = 'Single crystal'
        self.settings['Global azimuthal angle'] = 0.0
        self.settings['Angle of incidence'] = 0.0
        self.settings['Mode'] = 'Transfer matrix'
        self.settings['Frequency units'] = 'wavenumber'
        self.settings['Partially incoherent samples'] = 20
        self.settings['Percentage partial incoherence'] = 0
        self.settings['Filter kernel size'] = 1
        self.settings['Filter polynomial size'] = 3
        # Define a default superstrate/dielectric/substrate system
        self.settings['Layer material names']     = ['air',     'Dielectric layer','air'     ]
        self.settings['Layer hkls']               = [ [0,0,0],  [0,0,1],           [0,0,0]   ]
        self.settings['Layer azimuthals']         = [ 0,         0,                 0        ]
        self.settings['Layer thicknesses']        = [ 1,         1,                 1        ]
        self.settings['Layer thickness units']    = ['um',      'um',              'um'      ]
        self.settings['Layer incoherent options'] = ['Coherent','Coherent',        'Coherent']
        self.settings['Layer dielectric flags']   = [ False,     True,              False    ]
        # The maximum allowed thickness of a layer in metres
        # used to subdivide thicker films into many thinner films
        # if zero no subdivision is performed
        self.settings['Slice thickness'] = 0
        self.settings['Slice thickness unit'] = 'um'
        self.materialNames = []
        self.p_reflectance = []
        self.s_reflectance = []
        self.p_transmittance = []
        self.s_transmittance = []
        self.p_absorbtance = []
        self.s_absorbtance = []
        self.epsilon = []
        self.layers = []
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Get the last unit cell in the reader
        self.cell = None
        if self.reader is not None:
            self.cell = self.reader.unit_cells[-1]
        # Set the exponent threshold to be used by GTM
        self.exponent_threshold = 11000    
        # Open the database and get the material names
        self.DataBase = MaterialsDataBase(self.settings['Materials database'],debug=debugger.state())
        self.settings['Materials database'] = self.DataBase.getFileName()
        self.materialNames = self.setMaterialNames()
        # Create the layers - superstrate / dielectric / substrate from the defaults layer settings
        if self.reader is not None:
            self.settings2Layers()
        # Create last tab - SingleCrystalTab
        vbox = QVBoxLayout()
        self.form = QFormLayout()
        #
        # Option to open a database of permittivities for the support
        # label and button are defined the parent class
        self.form.addRow(self.openDB_label, self.openDB_button)
        #
        # Add a name for the database
        # label and button are defined the parent class
        self.form.addRow(self.database_le_label, self.database_le)
        #
        # Chose mode of operation
        #
        self.mode_cb = QComboBox(self)
        self.mode_cb.setToolTip('Set the method for calculating light transmission and reflectanceb;\n Transfer matrix.  This method is fast but can be numerically unstable.\n Scattering matrix. This method is slow but is numerically stable')
        self.mode_cb.addItems( ['Transfer matrix','Scattering matrix'] )
        index = self.mode_cb.findText(self.settings['Mode'], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        self.mode_cb.activated.connect(self.on_mode_cb_activated)
        label = QLabel('Single crystal methodology', self)
        label.setToolTip('Set the method for calculating light transmission and reflectanceb;\n Transfer matrix.  This method is fast but can be numerically unstable.\n Scattering matrix. This method is slow but is numerically stable')
        self.form.addRow(label, self.mode_cb)
        #
        # Define the global azimuthal angle widget
        #
        label,layout = self.globalAzimuthalWidget()
        self.form.addRow(label, layout)
        #
        # Define the angle of incidence widget
        #
        label, layout = self.angleOfIncidenceWidget()
        self.form.addRow(label,layout)
        #
        # Layer information widget
        #
        label = QLabel('Layer information')
        line  = QFrame()
        line.setFrameShape(QFrame.HLine)
        hbox = QHBoxLayout()
        hbox.addWidget(line)
        hbox.setAlignment(Qt.AlignVCenter)
        self.form.addRow(label,hbox)
        self.form.addRow(self.drawLayerTable())
        label = QLabel('    ')
        line  = QFrame()
        line.setFrameShape(QFrame.HLine)
        hbox = QHBoxLayout()
        hbox.addWidget(line)
        hbox.setAlignment(Qt.AlignVCenter)
        self.form.addRow(label,hbox)
        #
        # Widgets for setting the slice thickness
        #
        label,layout = self.sliceThicknessWidget()
        self.form.addRow(label, layout)
        #
        # Partial incoherence widget
        #
        label,layout = self.partialIncoherenceWidget()
        self.form.addRow(label, layout)
        #
        # Smoothing widget
        #
        label,layout = self.smoothingWidget()
        self.form.addRow(label, layout)
        #
        # Add a legend option
        #
        self.legend_le = QLineEdit(self)
        self.legend_le.setToolTip('The legend will be used to describe the results in the plot')
        self.legend_le.setText(self.settings['Legend'])
        self.legend_le.textChanged.connect(self.on_legend_le_changed)
        label = QLabel('Scenario legend')
        label.setToolTip('The legend will be used to describe the results in the plot')
        self.form.addRow(label, self.legend_le)
        #
        #
        # Final buttons for changing/deleting and switching the scenarios
        #
        hbox = self.add_scenario_buttons()
        self.form.addRow(hbox)
        #
        # Spacer item to soak up space
        #
        verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.form.addItem(verticalSpacer)
        vbox.addLayout(self.form)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        debugger.print('Finished:: initialiser')

    def redrawLayerTable(self):
        '''Redraw the layer table widget'''
        self.layerTable_tw.setRowCount(1)
        rowCount = 0
        for sequenceNumber,layer in enumerate(self.layers):
            rowCount += 1
            firstLayer = False
            lastLayer = False
            if layer == self.layers[0]:
                firstLayer = True
            elif layer == self.layers[-1]:
                lastLayer = True
            self.layerTable_tw.setRowCount(rowCount)
            self.redrawLayerTableRow(sequenceNumber,layer,rowCount,firstLayer,lastLayer)
        # Add a 'create new layer' button
        rowCount += 1
        newLayer_cb = self.newLayerWidget()
        newLayer_cb.setStyleSheet('Text-align:left')
        self.layerTable_tw.setRowCount(rowCount)
        self.layerTable_tw.setCellWidget(rowCount-1,0,newLayer_cb)

    def redrawLayerTableRow(self,sequenceNumber,layer,rowCount,firstLayer,lastLayer):
        '''Draw a row of the layer table
           Need a seperate routine for this as there are problems with the
           lambda code only keeping the last in a list
        '''
        # Create a layer button
        material = layer.getMaterial()
        materialName = material.getName()
        layer_button = QPushButton(materialName)
        layer_button.setToolTip('Show all the material properties in a new window')
        layer_button.setStyleSheet('Text-align:left')
        layer_button.clicked.connect(lambda x: self.on_layer_button_clicked(x,layer,sequenceNumber))
        self.layerTable_tw.setCellWidget(sequenceNumber,0,layer_button)
        # Handle thickness 
        materialThickness = layer.getThickness()
        thicknessUnit = layer.getThicknessUnit()
        film_thickness_sb = QDoubleSpinBox(self)
        film_thickness_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        film_thickness_sb.setToolTip('Define the thin film thickness in the defined thickness units')
        film_thickness_sb.setRange(0,100000)
        film_thickness_sb.setDecimals(3)
        film_thickness_sb.setSingleStep(0.001)
        film_thickness_sb.setValue(materialThickness)
        film_thickness_sb.valueChanged.connect(lambda x: self.on_film_thickness_sb_changed(x,layer))
        self.layerTable_tw.setCellWidget(sequenceNumber,1,film_thickness_sb)
        # thickness unit
        thickness_unit_cb = QComboBox(self)
        thickness_unit_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thickness_unit_cb.setToolTip('Set the units to be used for thickness; either angs nm, um, mm or cm')
        thickness_unit_cb.addItems( thickness_units )
        index = thickness_unit_cb.findText(thicknessUnit, Qt.MatchFixedString)
        thickness_unit_cb.setCurrentIndex(index)
        thickness_unit_cb.activated.connect(lambda x: self.on_thickness_units_cb_activated(x, layer))
        self.layerTable_tw.setCellWidget(sequenceNumber,2,thickness_unit_cb)
        # define hkl
        h_sb = QSpinBox(self)
        h_sb.setToolTip('Define the h dimension of the unique direction')
        h_sb.setRange(-20,20)
        h_sb.setSingleStep(1)
        h_sb.setValue(layer.getHKL()[0])
        h_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,0,layer))
        h_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        k_sb = QSpinBox(self)
        k_sb.setToolTip('Define the k dimension of the unique direction')
        k_sb.setRange(-20,20)
        k_sb.setSingleStep(1)
        k_sb.setValue(layer.getHKL()[1])
        k_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,1,layer))
        k_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        l_sb = QSpinBox(self)
        l_sb.setToolTip('Define the l dimension of the unique direction')
        l_sb.setRange(-20,20)
        l_sb.setSingleStep(1)
        l_sb.setValue(layer.getHKL()[2])
        l_sb.valueChanged.connect(lambda x: self.on_hkl_sb_changed(x,2,layer))
        l_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.layerTable_tw.setCellWidget(sequenceNumber,3,h_sb)
        self.layerTable_tw.setCellWidget(sequenceNumber,4,k_sb)
        self.layerTable_tw.setCellWidget(sequenceNumber,5,l_sb)
        # define azimuthal angle
        azimuthal = layer.getAzimuthal()
        azimuthal_angle_sb = QDoubleSpinBox(self)
        azimuthal_angle_sb.setToolTip('Define the slab azimuthal angle (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal in the laboratory frame can be seen in the laboratory frame information below')
        azimuthal_angle_sb.setRange(-180,360)
        azimuthal_angle_sb.setSingleStep(10)
        azimuthal_angle_sb.setValue(azimuthal)
        azimuthal_angle_sb.valueChanged.connect(lambda x: self.on_azimuthal_angle_sb_changed(x,layer))
        azimuthal_angle_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.layerTable_tw.setCellWidget(sequenceNumber,6,azimuthal_angle_sb)
        # Create a checkbox for coherence/incoherence
        option_cb = QComboBox(self)
        option_cb.setToolTip('Change optional settings for the layer')
        option_cb.addItems( incoherentOptions )
        index = option_cb.findText(layer.getIncoherentOption(), Qt.MatchFixedString)
        option_cb.setCurrentIndex(index)
        option_cb.activated.connect(lambda x: self.on_option_cb_activated(x,layer))
        option_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.layerTable_tw.setCellWidget(sequenceNumber,7,option_cb)
        # Create a toolbar for up down delete
        toolbar = self.createToolBar(layer,sequenceNumber,len(self.layers))
        self.layerTable_tw.setCellWidget(sequenceNumber,8,toolbar)
        # Add a Print option if debug is on
        if debugger.state():
            printButton = QPushButton('Print')
            printButton.setToolTip('Print the permittivity')
            printButton.clicked.connect(lambda x: self.on_print_button_clicked(x,layer))
            self.layerTable_tw.setCellWidget(sequenceNumber,9,printButton)
        if layer.isScalar():
            h_sb.setEnabled(False)
            k_sb.setEnabled(False)
            l_sb.setEnabled(False)
            azimuthal_angle_sb.setEnabled(False)
        if firstLayer or lastLayer:
            film_thickness_sb.setEnabled(False)
            thickness_unit_cb.setEnabled(False)
        return

    def drawLayerTable(self):
        '''Draw a table with all the layers in it '''
        self.layerTable_tw = QTableWidget()
        self.layerTable_tw.setToolTip('Define the layers in the system')
        self.layerTable_tw.itemChanged.connect(self.on_layerTable_itemChanged)
        self.layerTable_tw.setStyleSheet("QTableWidget::item {padding-left: 0px; border; 0px}")
        self.layerTable_tw.verticalHeader().setVisible(False)
        self.layerTable_tw.setShowGrid(False)
        headers = ['Material', 'Thickness', 'Units', 'H', 'K', 'L', 'Azimuthal', 'Options', 'Move']
        if debugger.state():
            headers.append('Print')
        self.layerTable_tw.setRowCount(1)
        self.layerTable_tw.setColumnCount(len(headers))
        self.layerTable_tw.setHorizontalHeaderLabels(headers)
        header = self.layerTable_tw.horizontalHeader()
        # Material 
        header.setSectionResizeMode(0,QHeaderView.Stretch)
        # Thickness and units 
        header.setSectionResizeMode(1,QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2,QHeaderView.ResizeToContents)
        # HKL
        header.setSectionResizeMode(3,QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4,QHeaderView.ResizeToContents)
        header.setSectionResizeMode(5,QHeaderView.ResizeToContents)
        # Azimuthal
        header.setSectionResizeMode(6,QHeaderView.ResizeToContents)
        # Options
        header.setSectionResizeMode(7,QHeaderView.ResizeToContents)
        # Move
        header.setSectionResizeMode(8,QHeaderView.ResizeToContents)
        if 'Print' in headers:
            # Print
            header.setSectionResizeMode(9,QHeaderView.ResizeToContents)
        if self.reader is not None:
            self.redrawLayerTable()
        return self.layerTable_tw

    def deleteLayer(self,x,layer,layerIndex):
        '''Handle a delete layer button press'''
        new = layerIndex + 1
        if layerIndex == 0 and layer[new].isTensor():
            #  Only allow scalar materials as the superstrate
            print('New superstrate material must be a scalar dielectric')
            return
        # Delete the layer
        del self.layers[layerIndex]
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def moveLayerUp(self,x,layer,layerIndex):
        '''Move a layer up (sequence number gets smaller by 1)'''
        if layerIndex < 1:
            return
        if layerIndex == 1 and layer.isTensor():
            #  Only allow scalar materials as the superstrate
            print('New superstrate material must be a scalar dielectric')
            return
        new = layerIndex - 1
        item = self.layers[layerIndex]
        self.layers.pop(layerIndex)
        self.layers.insert(new, item)
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def moveLayerDown(self,x,layer,layerIndex):
        '''Move a layer down (sequence number gets larger by 1)'''
        last = len(self.layers) - 1
        if layerIndex >= last:
            return
        new = layerIndex + 1
        if layerIndex == 0 and self.layers[new].isTensor():
            #  Only allow scalar materials as the superstrate
            print('New superstrate material must be a scalar dielectric')
            return
        item = self.layers[layerIndex]
        self.layers.pop(layerIndex)
        self.layers.insert(new, item)
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def createToolBarMoveUpButton(self,layer,layerIndex,nLayers):
        '''Create the move up button as part of the layer toolbar'''
        moveUpButton = QPushButton()
        moveUpButton.setIcon(QApplication.style().standardIcon(QStyle.SP_ArrowUp))
        moveUpButton.clicked.connect(lambda x: self.moveLayerUp(x,layer,layerIndex))
        moveUpButton.setFixedSize(20,20)
        moveUpButton.setIconSize(QSize(20,20))
        moveUpButton.setStyleSheet('border: none;')
        moveUpButton.setToolTip('Move this layer up the list of layers')
        return moveUpButton

    def createToolBarMoveDownButton(self,layer,layerIndex,nLayers):
        '''Create the move down button as part of the layer toolbar'''
        moveDownButton = QPushButton()
        moveDownButton.setIcon(QApplication.style().standardIcon(QStyle.SP_ArrowDown))
        moveDownButton.clicked.connect(lambda x: self.moveLayerDown(x,layer,layerIndex))
        moveDownButton.setFixedSize(20,20)
        moveDownButton.setIconSize(QSize(20,20))
        moveDownButton.setStyleSheet('border: none;')
        moveDownButton.setToolTip('Move this layer down the list of layers')
        return moveDownButton

    def createToolBarDeleteButton(self,layer,layerIndex,nLayers):
        '''Create the delete button as part of the layer toolbar'''
        deleteButton = QPushButton()
        deleteButton.setIcon(QApplication.style().standardIcon(QStyle.SP_DialogCloseButton))
        deleteButton.clicked.connect(lambda x: self.deleteLayer(x,layer,layerIndex))
        deleteButton.setFixedSize(20,20)
        deleteButton.setIconSize(QSize(20,20))
        deleteButton.setStyleSheet('border: none;')
        deleteButton.setToolTip('Delete this layer')
        return deleteButton

    def createToolBar(self,layer,layerIndex,nLayers):
        '''Create the tool bar used for the material layer
           layer is the layer concerned
           layerIndex is its index in the list
           nLayers is the number of layers in the list'''
        if nLayers <= 1:
            return
        frame = QFrame()
        frame_layout = QHBoxLayout()
        frame.setLayout(frame_layout)
        # Create the buttons in different routines because of the lambda function usage
        moveUpButton   =  self.createToolBarMoveUpButton(layer,layerIndex,nLayers)
        moveDownButton =  self.createToolBarMoveDownButton(layer,layerIndex,nLayers)
        deleteButton   =  self.createToolBarDeleteButton(layer,layerIndex,nLayers)
        nextIndex = layerIndex+1
        previousIndex = layerIndex-1
        # disable any buttons that are irrelevant to the layer
        if layerIndex == 0:
            moveUpButton.setEnabled(False)
        if layerIndex == nLayers-1:
            moveDownButton.setEnabled(False)
        if layerIndex == 0 and self.layers[nextIndex].isTensor():
            moveDownButton.setEnabled(False)
            deleteButton.setEnabled(False)
        if layerIndex == 1 and layer.isTensor():
            moveUpButton.setEnabled(False)
        # Add the buttons to the frame and return the frame
        frame_layout.addWidget(moveUpButton)
        frame_layout.addWidget(moveDownButton)
        frame_layout.addWidget(deleteButton)
        return frame
        
    def newLayerWidget(self):
        newLayer_cb = QComboBox()
        newLayer_cb.setToolTip('Create a new layer')
        materialNames = ['New layer...']
        materialNames += self.materialNames
        newLayer_cb.addItems(materialNames)
        newLayer_cb.setCurrentIndex(0)
        newLayer_cb.activated.connect(self.on_newLayer_cb_activated)
        return newLayer_cb

    def on_layerTable_itemChanged(self,item):
        '''Handle a change to the layer table'''
        print('on_layerTable_itemChanged: ',item)
        return

    def on_newLayer_cb_activated(self,index):
        '''Handle a new layer button click'''
        if index == 0:
            return
        # Subtract 1 from the index because the widget thinks the list includes 'New layer...' at the start
        newMaterialName = self.materialNames[index-1]
        if 'manual' in newMaterialName:
            return
        newMaterial = self.getMaterialFromDataBase(newMaterialName)
        hkl = [0,0,0]
        if newMaterial.isTensor():
            hkl = [0,0,1]
        new_layer = Layer(newMaterial,hkl=hkl,azimuthal=0.0,thickness=1.0,thicknessUnit='um')
        self.layers.append(new_layer)
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired = True
        return

    def on_print_button_clicked(self,x,layer):
        '''Print the permittivity for the layer'''
        material = layer.getMaterial()
        permittivityObject = material.getPermittivityObject()
        name = material.getName()
        name = name.replace(' ','_')
        name += '_permittivity.csv'
        print('Printing permittivity information to',name)
        permittivityObject.print(0.0,2000.0,1.0,file=name)
        return

    def on_option_cb_activated(self,index,layer):
        '''The incoherence option combo box has been activated'''
        debugger.print('on_incoherence_cb_activated', index,layer.getName())
        option = incoherentOptions[index]
        layer.setIncoherentOption(option)
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_film_thickness_sb_changed(self,value,layer):
        '''Handle film thickness spin box change'''
        debugger.print('on_film_thickness_sb_changed', value, layer.getName())
        layer.setThickness(value)
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_thickness_units_cb_activated(self, index, layer):
        debugger.print('Start:: on_thickness_units_cb_activated',index,layer.getName())
        unit = thickness_units[index]
        layer.setThicknessUnit(unit)
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_azimuthal_angle_sb_changed(self,value,layer):
        '''Handle azimuthal spin box change'''
        debugger.print('on_azimuthal_angl_sb_changed', value, layer.getName())
        layer.setAzimuthal(value)
        layer.changeLabFrameInfo()
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def on_hkl_sb_changed(self,value,hkorl,layer):
        '''Handle h, k or l change'''
        debugger.print('on_hkl_sb_changed', value)
        hkl = layer.getHKL()
        hkl[hkorl] = value
        layer.setHKL(hkl)
        layer.changeLabFrameInfo()
        self.generateLayerSettings()
        self.refreshRequired = True
        return

    def generateLayerSettings(self):
        '''Generate a list of settings for each layer'''
        debugger.print(self.settings['Legend'],'generateLayerSettings')
        self.settings['Layer material names']  = []
        self.settings['Layer hkls']            = []
        self.settings['Layer azimuthals']      = []
        self.settings['Layer thicknesses']     = []
        self.settings['Layer thickness units'] = []
        self.settings['Layer dielectric flags'] = []
        self.settings['Layer incoherent options'] = []
        for layer in self.layers:
            self.settings['Layer material names'].append(layer.getMaterial().getName())
            self.settings['Layer hkls'].append(layer.getHKL())
            self.settings['Layer azimuthals'].append(layer.getAzimuthal())
            self.settings['Layer thicknesses'].append(layer.getThickness())
            self.settings['Layer thickness units'].append(layer.getThicknessUnit())
            self.settings['Layer dielectric flags'].append(layer.isDielectric())
            self.settings['Layer incoherent options'].append(layer.getIncoherentOption())
        return

    def setMaterialNames(self):
        # Get the list of material names from the database
        materialNames = self.DataBase.getSheetNames()
        materialNames.append('Dielectric layer')
        return materialNames

    def printLayerSettings(self,message):
        print(message)
        print(self.settings['Layer material names'])
        print(self.settings['Layer hkls'])
        print(self.settings['Layer azimuthals'])
        print(self.settings['Layer thicknesses'])
        print(self.settings['Layer thickness units'])
        print(self.settings['Layer dielectric flags'])
        print(self.settings['Layer incoherent options'])
        print()

    def settings2Layers(self):
        '''Read the layer settings and generate a list of layers'''
        debugger.print(self.settings['Legend'],'settings2Layers')
        self.layers = []
        self.materialNames = self.setMaterialNames()
        # Process the settings information and append each layer to the list
        for index, (name, hkl, azimuthal, thickness, thicknessUnit, dielectricFlag, incoherentOption) in enumerate(zip(
                          self.settings['Layer material names'],
                          self.settings['Layer hkls'],
                          self.settings['Layer azimuthals'],
                          self.settings['Layer thicknesses'],
                          self.settings['Layer thickness units'],
                          self.settings['Layer dielectric flags'],
                          self.settings['Layer incoherent options'])):
            if name not in self.materialNames:
                print('Error material ', name, ' not available ', self.materialNames)
                name = 'air'
            material = self.getMaterialFromDataBase(name)
            self.layers.append(Layer(material,hkl=hkl,azimuthal=azimuthal,
                                     thickness=thickness,thicknessUnit=thicknessUnit,
                                     incoherentOption=incoherentOption,dielectricFlag=dielectricFlag))
        return

    def getMaterialFromDataBase(self,name,permittivity=None):
        '''Get the material from the database
           name         is the name of the spreadsheet
                        it can also be 'Dielectric layer' or 'Material defined manually'
           permittivity is the permittivity of a 'Defined manually material' material
        '''
        if name == 'Dielectric layer':
            # Create the dielectric material
            crystalPermittivityObject = self.notebook.settingsTab.getCrystalPermittivityObject()
            material = Materials.External(name,permittivityObject=crystalPermittivityObject,cell=self.cell)
        elif name == 'Material defined manually':
            material = Materials.Constant('Material defined manually',permittivity=permittivity)
        else:
            # Get the material from the data base
            # set the units for frequency to Hz for all materials
            material = self.DataBase.getMaterial(name)
        return material

    def getDielectricLayerIndex(self):
        '''Return the index of the dielectric layer in the list of layers'''
        for index,layer in enumerate(self.layers):
            if layer.isDielectric():
                self.dielectricLayer = layer
                return index
        return None

    def on_layer_button_clicked(self,x,layer,layerIndex):
        '''Handle a click on the show layer widget'''
        # Create the dialog box with all the information on the layer, work on a copy of the layer
        showLayerWindow = ShowLayerWindow(copy.copy(layer),debug=debugger.state())
        if showLayerWindow.exec():
            # The 'Ok' button was pressed
            # get the new Layer and replace the old one
            self.layers[layerIndex] = showLayerWindow.getLayer()
            self.generateLayerSettings()
            self.refresh(force=True)
            self.refreshRequired = True
        return

    def globalAzimuthalWidget(self):
        '''Create an global azimuthal angle widget'''
        self.global_azimuthal_angle_sb = QDoubleSpinBox(self)
        self.global_azimuthal_angle_sb.setToolTip('Define the global azimuthal angle for all layers (rotation of the crystal about the lab Z-axis).\nThe orientation of the crystal layers in the laboratory frame is shown in the edit layers window')
        self.global_azimuthal_angle_sb.setRange(-180,360)
        self.global_azimuthal_angle_sb.setSingleStep(10)
        self.global_azimuthal_angle_sb.setValue(self.settings['Global azimuthal angle'])
        self.global_azimuthal_angle_sb.valueChanged.connect(self.on_global_azimuthal_angle_sb_changed)
        label = QLabel('Global azimuthal angle')
        label.setToolTip('Define the global azimuthal angle for all layers (rotation of the layers about the lab Z-axis).\nThe orientation of the crystal layers in the laboratory frame is shown in the edit layers window')
        return label,self.global_azimuthal_angle_sb

    def angleOfIncidenceWidget(self):
        '''Create an angle of incidence widget'''
        self.angle_of_incidence_sb = QDoubleSpinBox(self)
        self.angle_of_incidence_sb.setToolTip('Define the angle of incidence, (normal incidence is 0 degrees)')
        self.angle_of_incidence_sb.setRange(0,90)
        self.angle_of_incidence_sb.setSingleStep(5)
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        self.angle_of_incidence_sb.valueChanged.connect(self.on_angle_of_incidence_sb_changed)
        label = QLabel('Angle of incidence')
        label.setToolTip('Define the angle of incidence, (normal incidence is 0 degrees).')
        return label,self.angle_of_incidence_sb

    def partialIncoherenceWidget(self):
        '''Create a partial incoherence widget'''
        hbox = QHBoxLayout()
        self.percentage_partial_incoherence_sb = QSpinBox(self)
        self.percentage_partial_incoherence_sb.setToolTip('Define the maximum percentage changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees\nIf the value of the percentage incoherence is zero, no incoherence is calculated')
        self.percentage_partial_incoherence_sb.setRange(0,100)
        self.percentage_partial_incoherence_sb.setSingleStep(1)
        self.percentage_partial_incoherence_sb.setValue(self.settings['Percentage partial incoherence'])
        self.percentage_partial_incoherence_sb.valueChanged.connect(self.on_percentage_partial_incoherence_sb_changed)
        self.partially_incoherent_samples_sb = QSpinBox(self)
        self.partially_incoherent_samples_sb.setToolTip('Define the number of samples to be used in the calculation of an incoherent spectrum.\nA large number of samples will take a long time but will give smoother results.')
        self.partially_incoherent_samples_sb.setRange(0,10000)
        self.partially_incoherent_samples_sb.setSingleStep(1)
        self.partially_incoherent_samples_sb.setValue(self.settings['Partially incoherent samples'])
        self.partially_incoherent_samples_sb.valueChanged.connect(self.on_partially_incoherent_samples_sb_changed)
        hbox.addWidget(self.percentage_partial_incoherence_sb)
        hbox.addWidget(self.partially_incoherent_samples_sb)
        label = QLabel('Partial incoherence (percentage & samples)')
        label.setToolTip('Define the percentage partial incoherence and the number of samples to be used in the calculation of an incoherent spectrum.\nThe percentage reflects changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence.\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees.\nA large number of samples will take a long time but will give smoother results.\nIf the value of the percentage incoherence is zero, no incoherence is calculated')
        return label,hbox

    def smoothingWidget(self):
        '''Create a smoothing widget'''
        hbox = QHBoxLayout()
        self.partially_incoherent_kernel_sb = QSpinBox(self)
        self.partially_incoherent_kernel_sb.setRange(1,1001)
        self.partially_incoherent_kernel_sb.setSingleStep(2)
        self.partially_incoherent_kernel_sb.setValue(self.settings['Filter kernel size'])
        self.partially_incoherent_kernel_sb.valueChanged.connect(self.on_partially_incoherent_kernel_sb_changed)
        self.partially_incoherent_kernel_sb.setToolTip('Define the kernel size for the smoothing of incoherent spectra (must be an odd number)\nIf the kernel size is less than 3, no smoothing is done.\nThe larger the number, the smoother the spectrum but beware of too much smoothing.')
        self.partially_incoherent_polynomial_sb = QSpinBox(self)
        self.partially_incoherent_polynomial_sb.setToolTip('Define the maximum degree of polynomial to be used in the smoothing filter')
        self.partially_incoherent_polynomial_sb.setRange(2,10)
        self.partially_incoherent_polynomial_sb.setSingleStep(1)
        self.partially_incoherent_polynomial_sb.setValue(self.settings['Filter polynomial size'])
        self.partially_incoherent_polynomial_sb.valueChanged.connect(self.on_partially_incoherent_polynomial_sb_changed)
        hbox.addWidget(self.partially_incoherent_kernel_sb)
        hbox.addWidget(self.partially_incoherent_polynomial_sb)
        label = QLabel('Smoothing filter (kernel & polynomial size)')
        label.setToolTip('Define the kernel size for the smoothing of incoherent spectra (must be an odd number)\nIf the kernel size is less than 3, no smoothing is done.\nThe larger the number, the smoother the spectrum but beware of too much smoothing.\nAlso defines the polynomial size for the fitting of the points in the kernel')
        return label,hbox

    def sliceThicknessWidget(self):
        self.slice_thickness_sb = QSpinBox(self)
        self.slice_thickness_sb.setToolTip('Define a slice thickness to subdivide thick films\nA value of zero means no slicing is performed')
        self.slice_thickness_sb.setRange(0,10000)
        self.slice_thickness_sb.setSingleStep(1)
        self.slice_thickness_sb.setValue(self.settings['Slice thickness'])
        self.slice_thickness_sb.valueChanged.connect(self.on_slice_thickness_sb_changed)
        self.slice_thickness_sb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.slice_thickness_unit_cb = QComboBox(self)
        self.slice_thickness_unit_cb.setToolTip('Set the units to be used for thickness; either nm, um, mm or cm')
        self.slice_thickness_unit_cb.addItems( thickness_units )
        self.slice_thickness_unit_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thicknessUnit = self.settings['Slice thickness unit']
        index = self.slice_thickness_unit_cb.findText(thicknessUnit, Qt.MatchFixedString)
        self.slice_thickness_unit_cb.setCurrentIndex(index)
        self.slice_thickness_unit_cb.activated.connect(self.on_slice_thickness_unit_cb_activated)
        hbox = QHBoxLayout()
        hbox.addWidget(self.slice_thickness_sb)
        hbox.addWidget(self.slice_thickness_unit_cb)
        label = QLabel('Slice thickness')
        label.setToolTip('Thick films can be sliced into thinner films to overcome numerical problems\nSlicing is not performed if the slice thickness is 0')
        return label,hbox

    def on_slice_thickness_sb_changed(self,value):
        # We will need to recalculate everything for the new slice
        self.refreshRequired = True
        self.settings['Slice thickness'] = value
        return

    def on_slice_thickness_unit_cb_activated(self,index):
        # We will need to recalculate everything for the new slice
        self.refreshRequired = True
        self.settings['Slice thickness unit'] = thickness_units[index]
        return

    def openDB_button_clicked(self):
        '''Open a new materials' database'''
        debugger.print('Start:: openDB_button_clicked')
        self.openDataBase()
        self.refresh(force=True)
        self.refreshRequired = True
        return

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

    def on_global_azimuthal_angle_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_global_azimuthal_angl_sb_changed', value)
        self.refreshRequired = True
        self.settings['Global azimuthal angle'] = value
        return

    def on_angle_of_incidence_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_angle_of_incidence_sb_changed', value)
        self.refreshRequired = True
        self.settings['Angle of incidence'] = value
        return

    def refresh(self,force=False):
        debugger.print(self.settings['Legend'],'Start:: refresh, force =', force)
        if not self.refreshRequired and not force :
            debugger.print(self.settings['Legend'],'Finished:: refreshing widget aborted', self.refreshRequired,force)
            return
        # Check to see if there is a new reader, if there is set up the cell
        self.reader = self.notebook.reader
        if self.reader is None:
            return
        #
        # Determine the exponent_threshold to be used by GTM
        #
        if 'Linux' in self.notebook.mainTab.settings['Compatibility mode']:
            self.exponent_threshold = 11000
        elif 'Windows' in self.notebook.mainTab.settings['Compatibility mode']:
            self.exponent_threshold = 700
        else:
            self.exponent_threshold = 11000
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Initialise the cell
        self.cell = self.reader.unit_cells[-1]
        self.DataBase = MaterialsDataBase(self.settings['Materials database'],debug=debugger.state())
        self.settings['Materials database'] = self.DataBase.getFileName()
        self.database_le.setText(self.settings['Materials database'])
        # Update the possible  material names from the database
        self.materialNames = self.setMaterialNames()
        # Generate the layers from the settings
        self.settings2Layers()
        self.generateLayerSettings()
        # Force recalculation
        self.calculationRequired = True
        # Change any greyed out items
        self.greyed_out()
        #
        # Now refresh values that need updating
        #
        # Check the mode 
        index = self.mode_cb.findText(self.settings['Mode'], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        # Check partial incoherence settings
        if self.settings['Percentage partial incoherence'] > 0:
            self.noCalculationsRequired = self.settings['Partially incoherent samples']
        else:
            self.noCalculationsRequired = 1
        # Update the slice information
        self.slice_thickness_sb.setValue(self.settings['Slice thickness'])
        thicknessUnit = self.settings['Slice thickness unit']
        index = self.slice_thickness_unit_cb.findText(thicknessUnit, Qt.MatchFixedString)
        self.slice_thickness_unit_cb.setCurrentIndex(index)
        # Update the Legend widget
        self.legend_le.setText(self.settings['Legend'])
        # Update angle widgets
        self.global_azimuthal_angle_sb.setValue(self.settings['Global azimuthal angle'])
        self.angle_of_incidence_sb.setValue(self.settings['Angle of incidence'])
        # For partial incoherent case, set percentage variation of angles and thickness and the number of samples
        self.percentage_partial_incoherence_sb.setValue(self.settings['Percentage partial incoherence'])
        self.partially_incoherent_samples_sb.setValue(self.settings['Partially incoherent samples'])
        # For partial incoherent case, set the smoothing parameters
        self.partially_incoherent_kernel_sb.setValue(self.settings['Filter kernel size'])
        self.partially_incoherent_polynomial_sb.setValue(self.settings['Filter polynomial size'])
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        # Redraw the layer information widget
        self.redrawLayerTable()
        for i in range(20):
            QCoreApplication.processEvents()
        self.refreshRequired = False
        debugger.print(self.settings['Legend'],'Finished:: refresh, force =', force)
        return

    def on_mode_cb_activated(self, index):
        debugger.print(self.settings['Legend'],'Start:: on_mode_cb_activated')
        if index == 0:
            self.settings['Mode'] = 'Transfer matrix'
        elif index == 1:
            self.settings['Mode'] = 'Scattering matrix'
        if self.settings['Percentage partial incoherence'] > 0:
            self.noCalculationsRequired = self.settings['Partially incoherent samples']
        else:
            self.noCalculationsRequired = 1 
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired = True
        debugger.print(self.settings['Legend'],'Mode changed to ', self.settings['Mode'])
        debugger.print(self.settings['Legend'],'Finished:: on_mode_cb_activated')
        return

    def partially_incoherent_calculator( self,
                            layers,
                            mode,
                            theta,
                            phi,
                            psi,
                            angleOfIncidence,
                            sliceThickness, exponent_threshold):
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
        av_epsilon = np.zeros((size,3,3),dtype=np.cdouble)
        # Work out which of the layers is the crystal dielectric
        crystalLayer = None
        for layer in layers:
            if layer.isDielectric():
                crystalLayer = layer
        #
        # Loop over the number of samples requred
        #
        crystalDepth = crystalLayer.getThickness()
        keepCrystalDepth = crystalDepth
        d = crystalDepth
        t = theta
        p  = phi
        q  = psi
        a  = angleOfIncidence
        fractionalIncoherence = self.settings['Percentage partial incoherence']/100.0
        for s in range(self.settings["Partially incoherent samples"]):
            crystalDepth = d + d*( -1 + 2*np.random.rand())*fractionalIncoherence
            crystalLayer.setThickness(crystalDepth)
            theta = t + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            phi   = p + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            psi = q + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            angleOfIncidence = a + np.pi/2.0*( -1 +2*np.random.rand())*fractionalIncoherence
            ( p_reflectance, 
            s_reflectance, 
            p_transmittance, 
            s_transmittance, 
            s_absorbtance, 
            p_absorbtance, 
            epsilon) = self.coherent_calculator(
                                       layers,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence,
                                       sliceThickness, exponent_threshold)
            av_p_reflectance   += np.array(p_reflectance) / self.settings['Partially incoherent samples']
            av_s_reflectance   += np.array(s_reflectance) / self.settings['Partially incoherent samples']
            av_p_transmittance += np.array(p_transmittance) / self.settings['Partially incoherent samples']
            av_s_transmittance += np.array(s_transmittance) / self.settings['Partially incoherent samples']
            av_s_absorbtance   += np.array(p_absorbtance) / self.settings['Partially incoherent samples']
            av_p_absorbtance   += np.array(s_absorbtance) / self.settings['Partially incoherent samples']
            av_epsilon         += np.array(epsilon) / self.settings['Partially incoherent samples']
        crystalLayer.setThickness(keepCrystalDepth)
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
                            layers,
                            mode,
                            theta,
                            phi,
                            psi,
                            angleOfIncidence,
                            sliceThickness, exponent_threshold):
        """ Calculates the coherent component of light reflectance and transmission """
        debugger.print(self.settings['Legend'],'Entering the coherent_calculator function')
        #
        # Initialise the partial function to pass through to the pool
        #
        #print('coherent_calculator')
        #print('layers', layers)
        #print('mode', mode)
        #print('theta', theta)
        #print('phi', phi)
        #print('psi', psi)
        #print('angleOfIncidence', angleOfIncidence)
        #print('sliceThickness, exponent_threshold', sliceThickness, exponent_threshold)
        #
        partial_function = partial(solve_single_crystal_equations,
                                       layers,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence,
                                       sliceThickness, exponent_threshold)
        results = []
        # About to call
        debugger.print(self.settings['Legend'],'About to calculate single crystal scenario using pool')
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
        for v,r,R,t,T,eps,errors,largest_exponent in results:
            if errors > 0:
                print('Warning exponential overflow occured at frequency',v,errors,largest_exponent)
            p_reflectance.append(R[0]+R[2])
            s_reflectance.append(R[1]+R[3])
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
            debugger.print(self.settings['Legend'],'Finished:: calculate aborted because calculationRequired false')
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
        # Set the sliceThickness in metres
        tom = thickness_conversion_factors[self.settings['Slice thickness unit']]
        sliceThickness   = tom * self.settings['Slice thickness']
        # The euler angles are set to zero, apart from the global azimuthal angle
        # the rotation of each layer is now handled by the layer class.
        theta = 0.0
        phi = 0.0
        # global azimuthal angle in radians
        psi = np.radians(self.settings['Global azimuthal angle'])
        # Set the angle of incidence in radians
        angleOfIncidence = np.radians(self.settings['Angle of incidence'])
        # Tell each layer to calculate the euler matrix for rotation to the lab frame
        for layer in self.layers:
            if layer.isTensor():
                layer.calculate_euler_matrix()
        # Choose the calculator function depending on the mode of operation
        mode = self.settings['Mode']
        if self.settings['Percentage partial incoherence'] > 0:
            calculator = self.partially_incoherent_calculator
        else:
            calculator = self.coherent_calculator
        exponent_threshold = self.exponent_threshold
        # Call the relevant calculator
        ( self.p_reflectance, 
        self.s_reflectance, 
        self.p_transmittance, 
        self.s_transmittance, 
        self.s_absorbtance, 
        self.p_absorbtance, 
        self.epsilon) = calculator( self.layers,
                                    mode,
                                    theta,
                                    phi,
                                    psi,
                                    angleOfIncidence,
                                    sliceThickness, exponent_threshold)
        debugger.print(self.settings['Legend'],'Finished:: calculate - number of frequencies',len(vs_cm1))
        return

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
        # At the moment it appears there is nothing to do.
        debugger.print(self.settings['Legend'],'Finished:: greyed_out')
        return

