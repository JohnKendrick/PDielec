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

def solve_single_crystal_equations( 
        superstrateDielectricFunction ,
        substrateDielectricFunction   ,
        superstrateDepth              ,
        substrateDepth                ,
        layers                        ,
        crystalIncoherence            ,
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
        superstrateDielectricFunction is a dielectric function providing the superstrate dielectric
        substrateDielectricFunction   is a dielectric function providing the substrate dielectric
        superstrateDepth              the depth of the superstrate (m)
        substrateDepth                the depth of the substrate (m)
        layers                        a list of material layers (their permittivity functions)
        crystalIncoherence            a parameter between 0 and pi giving the amount of incoherence
        mode                          'thick slab' indicates the substrate will be the crystal so semi-infinite
                                      'coherent thin film' a standard GTM call with a single layer
                                      'incoherent thin film' multiple GTM calls with random phase
        theta                         the theta angle of the sla,
        phi                           the angle of the slab
        psi                           the psi angle of the slab
        angleOfIncidence              the angle of incidence
        sliceThickness                a thickness in m, used to subdivide thicker films
                                      if zero then the full fat film is used
        exponent_threshold            largest exponent allowed in calculation of the propagation matrix
        v                             the frequency of the light in cm-1

    """
    # Create layers, thickness is in metres
    superstrate      = GTM.CoherentLayer(thickness=superstrateDepth,epsilon1=superstrateDielectricFunction,exponent_threshold=exponent_threshold)
    substrate        = GTM.CoherentLayer(thickness=substrateDepth,  epsilon1=substrateDielectricFunction,exponent_threshold=exponent_threshold)
    selectedLayers = layers
    if mode == 'Thick slab':
        # For a thick slab the last layer is used as the thick layer
        # so redefined the substrate and remove the last layer from the list of layers
        selectedLayers = layers[:-1]
        substrateDepth = layers[-1].getThicknessInMetres()
        substrateDielectricFunction = layers[-1].getPermittivityFunction()
        substrate = GTM.CoherentLayer(thickness=substrateDepth,  epsilon=substrateDielectricFunction, exponent_threshold=exponent_threshold)
    gtmLayers = []
    for layer in selectedLayers:
        permittivityFunction = layer.getPermittivityFunction()
        depth = layer.getThicknessInMetres()
        coherent = layer.isCoherent()
        if sliceThickness != 0 and depth > sliceThickness:
            no_of_layers = int(depth / sliceThickness) + 1
            newdepth = depth / no_of_layers
            for i in range(no_of_layers):
                if coherent:
                    gtmLayers.append(GTM.CoherentLayer(thickness=newdepth, epsilon=permittivityFunction, exponent_threshold=exponent_threshold))
                else:
                    gtmLayers.append(GTM.IncoherentLayer(thickness=newdepth, epsilon=permittivityFunction, exponent_threshold=exponent_threshold))
        else:
            if coherent:
                gtmLayers.append(GTM.CoherentLayer(thickness=depth, epsilon=permittivityFunction, exponent_threshold=exponent_threshold))
            else:
                gtmLayers.append(GTM.IncoherentLayer(thickness=depth, epsilon=permittivityFunction, exponent_threshold=exponent_threshold))
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
    if mode == 'Incoherent thin film':
        Sys_Gamma = system.calculate_incoherent_GammaStar(freq, zeta_sys)
    else:
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
        self.settings['Superstrate'] = 'air'
        self.settings['Substrate'] = 'air'
        self.settings['Superstrate permittivity'] = 1.0
        self.settings['Substrate permittivity'] = 1.0
        self.settings['Superstrate depth'] = 99999.0
        self.settings['Substrate depth'] = 99999.0
        self.settings['Superstrate & substrate thickness unit'] = 'mm'
        self.settings['Mode'] = 'Thick slab'
        self.settings['Frequency units'] = 'wavenumber'
        self.settings['Partially incoherent samples'] = 20
        self.settings['Percentage partial incoherence'] = 100
        self.settings['Filter kernel size'] = 1
        self.settings['Filter polynomial size'] = 3
        self.settings['Layer material names'] = []
        self.settings['Layer hkls'] = []
        self.settings['Layer azimuthals'] = []
        self.settings['Layer thicknesses'] = []
        self.settings['Layer thickness units'] = []
        self.settings['Layer dielectric flags'] = []
        self.settings['Layer coherent flags'] = []
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
        self.dielectricLayer = None
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Get the last unit cell in the reader
        if self.reader is not None:
            self.cell = self.reader.unit_cells[-1]
        # Set the exponent threshold to be used by GTM
        self.exponent_threshold = 11000    
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
        self.mode_cb.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two semi-infinite media; the incident and the crystal),\n Partially incoherent, incoherent and coherent thin films assume there are three media; the incident, the crystal and the substrate. Partially incoherent film mode takes a lot longer than any of the other modes.')
        self.mode_cb.addItems( ['Thick slab','Coherent thin film','Incoherent thin film','Partially incoherent thin film'] )
        self.settings['Mode'] = 'Thick slab'
        self.mode_cb.activated.connect(self.on_mode_cb_activated)
        label = QLabel('Single crystal mode', self)
        label.setToolTip('Set the mode of operation for this tab;\n Thick slab means that only reflections are significant (the film thickness has no effect and there are only two semi-infinite media; the incident and the crystal),\n Partially incoherent, incoherent and coherent thin films assume there are three media; the incident, the crystal and the substrate.  Partially incoherent film mode takes a lot longer to run than any of the other modes.')
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
        # Define the superstrate material widget
        #
        label,layout = self.superstrateWidget()
        self.form.addRow(label, layout)
        #
        # Define the superstrate real and imaginary components
        #
        label,layout = self.superstratePermittivityWidget()
        self.form.addRow(label, layout)
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
        self.layerInformationFormCount = self.form.rowCount()-1
        label = QLabel('    ')
        line  = QFrame()
        line.setFrameShape(QFrame.HLine)
        hbox = QHBoxLayout()
        hbox.addWidget(line)
        hbox.setAlignment(Qt.AlignVCenter)
        self.form.addRow(label,hbox)
        #
        # Define the substrate material widget
        #
        label,layout = self.substrateWidget()
        self.form.addRow(label, layout)
        #
        # Define the substrate real and imaginary components
        #
        label,layout = self.substratePermittivityWidget()
        self.form.addRow(label, layout)
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
        #self.layerTable_tw.clearContents()
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
        film_thickness_sb.setSingleStep(0.01)
        film_thickness_sb.setValue(materialThickness)
        film_thickness_sb.valueChanged.connect(lambda x: self.on_film_thickness_sb_changed(x,layer))
        self.layerTable_tw.setCellWidget(sequenceNumber,1,film_thickness_sb)
        # thickness unit
        thickness_unit_cb = QComboBox(self)
        thickness_unit_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        thickness_unit_cb.setToolTip('Set the units to be used for thickness; either nm, um, mm or cm')
        thickness_unit_cb.addItems( ['nm','um','mm','cm'] )
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
        # Create a checkbox for incoherence
        option_cb = QCheckBox(self)
        option_cb.setToolTip('Change optional settings for the layer')
        option_cb.setText('Options')
        option_cb.setLayoutDirection(Qt.RightToLeft)
        if layer.isCoherent():
            option_cb.setCheckState(Qt.Unchecked)
        else:
            option_cb.setCheckState(Qt.Checked)
        option_cb.stateChanged.connect(lambda x: self.on_option_cb_changed(x,layer))
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
        return

    def drawLayerTable(self):
        '''Draw a table with all the layers in it '''
        self.layerTable_tw = QTableWidget()
        self.layerTable_tw.setToolTip('Define the layers in the system')
        self.layerTable_tw.itemChanged.connect(self.on_layerTable_itemChanged)
        self.layerTable_tw.setStyleSheet("QTableWidget::item {padding-left: 0px; border; 0px}")
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
        self.redrawLayerTable()
        return self.layerTable_tw

    def deleteLayer(self,x,layer,layerIndex):
        '''Handle a delete layer button press'''
        # Only delete a layer if it doesn't have the dielectricFlag set
        if not layer.isDielectric():
            del self.layers[layerIndex]
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def moveLayerUp(self,x,layer,layerIndex):
        '''Move a layer up (sequence number gets smaller by 1)'''
        print('JK30 move layer up',x,layer,layerIndex)
        if layerIndex < 1:
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
        item = self.layers[layerIndex]
        self.layers.pop(layerIndex)
        self.layers.insert(new, item)
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired=True
        return

    def createToolBarMoveUpButton(self,layer,layerIndex,nLayers):
        '''Create the move up button as part of the layer toolbar'''
        print('JK20','Move up',layer,layerIndex,nLayers)
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
        frame = QFrame()
        frame_layout = QHBoxLayout()
        frame.setLayout(frame_layout)
        # Create the buttons in different routines because of the lambda function usage
        moveUpButton   =  self.createToolBarMoveUpButton(layer,layerIndex,nLayers)
        moveDownButton =  self.createToolBarMoveDownButton(layer,layerIndex,nLayers)
        deleteButton   =  self.createToolBarDeleteButton(layer,layerIndex,nLayers)
        # disable any buttons that are irrelevant to the layer
        if layerIndex == 0:
            moveUpButton.setEnabled(False)
        if layerIndex == nLayers-1:
            moveDownButton.setEnabled(False)
        if nLayers == 1:
            moveUpButton.setEnabled(False)
            moveDownButton.setEnabled(False)
            deleteButton.setEnabled(False)
        if layer.isDielectric():
            deleteButton.setEnabled(False)
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
        newMaterialName = self.materialNames[index-1]
        newMaterial = self.DataBase.getMaterial(newMaterialName)
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

    def on_option_cb_changed(self,value,layer):
        '''The incoherence checkbox has changed'''
        debugger.print('on_incoherence_cb_changed', value,layer.getName())
        state = layer.isCoherent()
        if state:
            layer.setIncoherent()
        else:
            layer.setCoherent()
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
        units = ['nm','um','mm','cm']
        unit = units[index]
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
        self.settings['Layer coherent flags'] = []
        for layer in self.layers:
            self.settings['Layer material names'].append(layer.getMaterial().getName())
            self.settings['Layer hkls'].append(layer.getHKL())
            self.settings['Layer azimuthals'].append(layer.getAzimuthal())
            self.settings['Layer thicknesses'].append(layer.getThickness())
            self.settings['Layer thickness units'].append(layer.getThicknessUnit())
            self.settings['Layer dielectric flags'].append(layer.isDielectric())
            self.settings['Layer coherent flags'].append(layer.isCoherent())
        return

    def addDielectricLayer(self,name,hkl,azimuthal,thickness,thicknessUnit,coherentFlag,forTheFirstTime=False):
        '''Add a dielectric layer to the list of layers'''
        debugger.print(self.settings['Legend'],'addDielectricLayer',name,hkl,azimuthal,thicknessUnit,forTheFirstTime)
        # Update the notebook permittivity if it needs to be updated
        temp = self.notebook.settingsTab.getCrystalPermittivity(self.vs_cm1)
        # Make sure the system knows that frequency will be supplied using Hz
        crystalPermittivityObject = self.notebook.settingsTab.getCrystalPermittivityObject()
        crystalPermittivityObject.setUnits('hz')
        # Create the dielectric material
        dielectricMaterial = Materials.External(name,permittivityObject=crystalPermittivityObject,cell=self.cell)
        if forTheFirstTime:
            # There are no dielectric layers so we have to add one
            # If older versions (before version 8.0.0) are used then
            # Different settings are used to define the layer.  If these settings are encountered then
            # use them once and then delete them
            debugger.print(self.settings['Legend'],'processing for the first time')
            if 'Film thickness' in self.settings:
                thickness = self.settings['Film thickness']
                del(self.settings['Film thickness'])
            if 'Thickness unit' in self.settings:
                thicknessUnit = self.settings['Thickness unit']
                del(self.settings['Thickness unit'])
            if 'Azimuthal angle' in self.settings:
                azimuthal = self.settings['Azimuthal angle']
                del(self.settings['Azimuthal angle'])
            if 'Unique direction - h' in self.settings:
                hkl = [ self.settings['Unique direction - h'], self.settings['Unique direction - k'], self.settings['Unique direction - l']]
                del(self.settings['Unique direction - h'])
                del(self.settings['Unique direction - k'])
                del(self.settings['Unique direction - l'])
        # Append the dielectric material to the list of layers
        self.layers.append(Layer(dielectricMaterial,hkl=hkl,azimuthal=azimuthal,thickness=thickness,thicknessUnit=thicknessUnit,coherentFlag=coherentFlag,dielectricFlag=True))
        self.generateLayerSettings()
        return

    def settings2Layers(self):
        '''Read the layer settings and generate a list of layers'''
        debugger.print(self.settings['Legend'],'settings2Layers')
        self.layers = []
        # Get the list of material names from the database
        self.materialNames = self.DataBase.getSheetNames()
        # If there is no layer information in settings then just add the dielectric
        if len(self.settings['Layer material names']) == 0:
            # Defaults come from version 7.3 of the code
            name = 'Dielectric film'
            hkl = [0, 0, 1]
            azimuthal = 0.0
            thickness = 10.0
            thicknessUnit = 'um'
            coherentFlag = True
            self.addDielectricLayer(name,hkl,azimuthal,thickness,thicknessUnit,coherentFlag,forTheFirstTime=True)
            return
        # Process the settings information and append each layer to the list
        for name, hkl, azimuthal, thickness, thicknessUnit, dielectricFlag, coherentFlag in zip(
                          self.settings['Layer material names'],
                          self.settings['Layer hkls'],
                          self.settings['Layer azimuthals'],
                          self.settings['Layer thicknesses'],
                          self.settings['Layer thickness units'],
                          self.settings['Layer dielectric flags'],
                          self.settings['Layer coherent flags']):
            # If the dielectric flag is set then initialise the dielectric and create a new layer
            # Otherwise use the materials database
            if dielectricFlag:
                self.addDielectricLayer(name,hkl,azimuthal,thickness,thicknessUnit,coherentFlag)
            else:
                if name not in self.materialNames:
                    print('Error material ', name, ' not available ', self.materialNames)
                    name = 'air'
                material = self.DataBase.getMaterial(name)
                self.layers.append(Layer(material,hkl=hkl,azimuthal=azimuthal,thickness=thickness,thicknessUnit=thicknessUnit,coherentFlag=coherentFlag,dielectricFlag=dielectricFlag))
        return

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
        print('on_layer_button_clicked',layer,layerIndex)
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
        self.percentage_partial_incoherence_sb.setToolTip('Define the maximum percentage changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees')
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
        label.setToolTip('Define the percentage partial incoherence and the number of samples to be used in the calculation of an incoherent spectrum.\nThe percentage reflects changes in the slab geometric parameters (thickness, alpha, beta and gamma euler angles) and the angle of incidence.\nFor thickness it is a percentage of the required thickness.  For angles it is a percentage of 90 degrees.\nA large number of samples will take a long time but will give smoother results.')
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

    def superstrateWidget(self):
        '''Create the superstrate widget '''
        hbox = QHBoxLayout()
        self.superstrate_cb = QComboBox(self)
        self.superstrate_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.superstrate_cb.setToolTip('Define the incident medium')
        self.materialNames = self.DataBase.getSheetNames()
        if 'manually' in self.settings['Superstrate']:
            self.materialNames.append('Material defined manually')
        self.superstrate_cb.clear()
        self.superstrate_cb.addItems(self.materialNames)
        if not self.settings['Superstrate'] in self.materialNames:
            if 'air' in self.materialNames:
                self.settings['Superstrate'] = 'air'
            else:
                self.settings['Superstrate'] = self.materialNames[0]
        index = self.superstrate_cb.findText(self.settings['Superstrate'], Qt.MatchFixedString)
        self.superstrate_cb.setCurrentIndex(index)
        self.superstrate_cb.activated.connect(self.on_superstrate_cb_activated)
        if 'Material defined manually' in self.settings['Superstrate']:
            self.superstrateMaterial = Materials.Constant('manual',permittivity=self.settings['Superstrate permittivity'])
        else:
            self.superstrateMaterial = self.DataBase.getMaterial(self.settings['Superstrate'])
        label = QLabel('Superstrate material')
        label.setToolTip('Define the incident medium permittivity')
        self.superstrate_info_le = QLineEdit(self)
        self.superstrate_info_le.setToolTip("Provides details about database entry")
        text = self.superstrateMaterial.getInformation()
        self.superstrate_info_le.setText(text)
        self.superstrate_info_le.setReadOnly(True)
        hbox.addWidget(self.superstrate_cb)
        hbox.addWidget(self.superstrate_info_le)
        return label,hbox

    def superstratePermittivityWidget(self):
        ''' Create a superstrate permittivity widget'''
        self.superstrate_permittivity_r_sb = QDoubleSpinBox(self)
        self.superstrate_permittivity_r_sb.setToolTip('Define the incident medium permittivity (real component)')
        self.superstrate_permittivity_r_sb.setRange(0,1000)
        self.superstrate_permittivity_r_sb.setSingleStep(0.1)
        self.superstrate_permittivity_r_sb.setValue(self.settings['Superstrate permittivity'].real)
        self.superstrate_permittivity_r_sb.valueChanged.connect(self.on_superstrate_permittivity_r_sb_changed)
        self.superstrate_permittivity_i_sb = QDoubleSpinBox(self)
        self.superstrate_permittivity_i_sb.setToolTip('Define the incident medium permittivity (imaginary component)')
        self.superstrate_permittivity_i_sb.setRange(0,1000)
        self.superstrate_permittivity_i_sb.setSingleStep(0.1)
        self.superstrate_permittivity_i_sb.setValue(self.settings['Superstrate permittivity'].imag)
        self.superstrate_permittivity_i_sb.valueChanged.connect(self.on_superstrate_permittivity_i_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.superstrate_permittivity_r_sb)
        hbox.addWidget(self.superstrate_permittivity_i_sb)
        label = QLabel('Superstrate permittivity')
        label.setToolTip('Define the complex permittivity of the superstrate')
        return label,hbox

    def substrateWidget(self):
        ''' Create a substrate widget'''
        hbox = QHBoxLayout()
        self.substrate_cb = QComboBox(self)
        self.substrate_cb.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        self.substrate_cb.setToolTip('Define the substrate material')
        self.materialNames = self.DataBase.getSheetNames()
        if 'manually' in self.settings['Substrate']:
            self.materialNames.append('Material defined manually')
        self.substrate_cb.clear()
        self.substrate_cb.addItems(self.materialNames)
        if not self.settings['Substrate'] in self.materialNames:
            if 'air' in self.materialNames:
                self.settings['Substrate'] = 'air'
            else:
                self.settings['Substrate'] = self.materialNames[0]
        index = self.substrate_cb.findText(self.settings['Substrate'], Qt.MatchFixedString)
        self.substrate_cb.setCurrentIndex(index)
        self.substrate_cb.activated.connect(self.on_substrate_cb_activated)
        if 'Material defined manually' in self.settings['Substrate']:
            self.substrateMaterial = Materials.Constant('manual',permittivity=self.settings['Substrate permittivity'])
        else:
            self.substrateMaterial = self.DataBase.getMaterial(self.settings['Substrate'])
        label = QLabel('Substrate material')
        label.setToolTip('Define the substrate permittivity')
        self.substrate_info_le = QLineEdit(self)
        self.substrate_info_le.setToolTip("Provides details about database entry")
        text = self.substrateMaterial.getInformation()
        self.substrate_info_le.setText(text)
        self.substrate_info_le.setReadOnly(True)
        hbox.addWidget(self.substrate_cb)
        hbox.addWidget(self.substrate_info_le)
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
        self.slice_thickness_unit_cb.addItems( ['nm','um','mm','cm'] )
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
        units = ['nm', 'um', 'mm', 'cm']
        self.settings['Slice thickness unit'] = units[index]
        return

    def substratePermittivityWidget(self):
        ''' Create a substrate permittivity widget'''
        self.substrate_permittivity_r_sb = QDoubleSpinBox(self)
        self.substrate_permittivity_r_sb.setToolTip('Define the substrate permittivity (real component)')
        self.substrate_permittivity_r_sb.setRange(0,1000)
        self.substrate_permittivity_r_sb.setSingleStep(0.1)
        self.substrate_permittivity_r_sb.setValue(self.settings['Substrate permittivity'].real)
        self.substrate_permittivity_r_sb.valueChanged.connect(self.on_substrate_permittivity_r_sb_changed)
        self.substrate_permittivity_i_sb = QDoubleSpinBox(self)
        self.substrate_permittivity_i_sb.setToolTip('Define the substrate permittivity (imaginary component)')
        self.substrate_permittivity_i_sb.setRange(0,1000)
        self.substrate_permittivity_i_sb.setSingleStep(0.1)
        self.substrate_permittivity_i_sb.setValue(self.settings['Substrate permittivity'])
        self.substrate_permittivity_i_sb.valueChanged.connect(self.on_substrate_permittivity_i_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.substrate_permittivity_r_sb)
        hbox.addWidget(self.substrate_permittivity_i_sb)
        label = QLabel('Substrate permittivity')
        label.setToolTip('Define the complex permittivity of the substrate')
        return label,hbox

    def openDB_button_clicked(self):
        '''Open a new materials' database'''
        debugger.print('Start:: openDB_button_clicked')
        self.openDataBase()
        if not self.settings['Substrate'] in self.materialNames:
            if 'air' in self.materialNames:
               self.settings['Substrate'] = 'air'
            else:
               self.settings['Substrate'] = self.materialNames[0]
        if not self.settings['Superstrate'] in self.materialNames:
            if 'air' in self.materialNames:
               self.settings['Superstrate'] = 'air'
            else:
               self.settings['Superstrate'] = self.materialNames[0]
        self.substrateMaterial = self.DataBase.getMaterial(self.settings['Substrate'])
        self.superstrateMaterial = self.DataBase.getMaterial(self.settings['Superstrate'])
        # Check to see that the materials return a scalar permittivity
        if self.substrateMaterial.isTensor():
            print('Error: substrate must have a scalar permittivity, using air')
            self.settings['Substrate'] = 'air'
            self.substrateMaterial = self.DataBase.getMaterial(self.settings['Substrate'])
        if self.superstrateMaterial.isTensor():
            print('Error: superstrate must have a scalar permittivity, using air')
            self.settings['Superstrate'] = 'air'
            self.superstrateMaterial = self.DataBase.getMaterial(self.settings['Superstrate'])
        self.settings['Substrate permittivity'] = self.substrateMaterial.getPermittivityFunction()(0.0)
        self.settings['Superstrate permittivity'] = self.superstrateMaterial.getPermittivityFunction()(0.0)
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

    def on_superstrate_permittivity_r_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_superstrate_permittivity_r_sb_changed', value)
        self.refreshRequired = True
        old_imag = self.settings['Superstrate permittivity'].imag
        self.settings['Superstrate permittivity'] = complex(value,old_imag)
        newPermittivityObject = DielectricFunction.ConstantScalar(self.settings['Superstrate permittivity'])
        self.superstrateMaterial.setPermittivityObject(newPermittivityObject)
        self.settings['Superstrate'] = 'Material defined manually'
        if 'Material defined manually' not in self.materialNames:
            self.materialNames.append('Material defined manually')
        self.refresh(force=True)
        self.refreshRequired = True
        return

    def on_superstrate_permittivity_i_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_superstrate_permittivity_i_sb_changed', value)
        self.refreshRequired = True
        old_real = self.settings['Superstrate permittivity'].real
        self.settings['Superstrate permittivity'] = complex(old_real,value)
        newPermittivityObject = DielectricFunction.ConstantScalar(self.settings['Superstrate permittivity'])
        self.superstrateMaterial.setPermittivityObject(newPermittivityObject)
        self.settings['Superstrate'] = 'Material defined manually'
        if 'Material defined manually' not in self.materialNames:
            self.materialNames.append('Material defined manually')
        self.refresh(force=True)
        self.refreshRequired = True
        return

    def on_substrate_permittivity_r_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_substrate_permittivity_r_sb_changed', value)
        self.refreshRequired = True
        old_imag = self.settings['Substrate permittivity'].imag
        self.settings['Substrate permittivity'] = complex(value,old_imag)
        newPermittivityObject = DielectricFunction.ConstantScalar(self.settings['Substrate permittivity'])
        self.substrateMaterial.setPermittivityObject(newPermittivityObject)
        self.settings['Substrate'] = 'Material defined manually'
        if 'Material defined manually' not in self.materialNames:
            self.materialNames.append('Material defined manually')
        self.refresh(force=True)
        self.refreshRequired = True
        return

    def on_substrate_permittivity_i_sb_changed(self,value):
        debugger.print(self.settings['Legend'],'on_substrate_permittivity_i_sb_changed', value)
        self.refreshRequired = True
        old_real = self.settings['Substrate permittivity'].real
        self.settings['Substrate permittivity'] = complex(old_real,value)
        newPermittivityObject = DielectricFunction.ConstantScalar(self.settings['Substrate permittivity'])
        self.substrateMaterial.setPermittivityObject(newPermittivityObject)
        self.settings['Substrate'] = 'Material defined manually'
        if 'Material defined manually' not in self.materialNames:
            self.materialNames.append('Material defined manually')
        self.refresh(force=True)
        self.refreshRequired = True
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
        # Update the possible super/substrate material names from the database
        self.materialNames = self.DataBase.getSheetNames()
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
        self.superstrate_cb.clear()
        self.superstrate_cb.addItems(self.materialNames)
        self.substrate_cb.clear()
        self.substrate_cb.addItems(self.materialNames)
        # If the permittivity has been edited manually in the GUI then update the combo box
        if 'manually' in self.settings['Superstrate'] or 'manually' in self.settings['Substrate']:
            self.materialNames.append('Material defined manually')
            if 'manually' in self.settings['Superstrate']:
                self.superstrate_cb.clear()
                self.superstrate_cb.addItems(self.materialNames)
            if 'manually' in self.settings['Substrate']:
                self.substrate_cb.clear()
                self.substrate_cb.addItems(self.materialNames)
        # Check the mode and if it is Incoherent increase the number of spectra to be calculated
        index = self.mode_cb.findText(self.settings['Mode'], Qt.MatchFixedString)
        self.mode_cb.setCurrentIndex(index)
        if self.settings['Mode'] == 'Partially incoherent thin film':
            self.noCalculationsRequired = self.settings['Partially incoherent samples']
        else:
            self.settings['Percentage partial incoherence'] = 0
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
        # Set the superstrate material
        self.setSuperstrate()
        # Set the substrate material
        self.setSubstrate()
        # Update super and substrate permittivities widgets
        index = self.superstrate_cb.findText(self.settings['Superstrate'], Qt.MatchFixedString)
        self.superstrate_cb.setCurrentIndex(index)
        index = self.substrate_cb.findText(self.settings['Substrate'], Qt.MatchFixedString)
        self.substrate_cb.setCurrentIndex(index)
        debugger.print('Updating permittivity widgets')
        self.superstrate_permittivity_r_sb.setValue(self.settings['Superstrate permittivity'].real)
        self.superstrate_permittivity_i_sb.setValue(self.settings['Superstrate permittivity'].imag)
        self.substrate_permittivity_r_sb.setValue(self.settings['Substrate permittivity'].real)
        self.substrate_permittivity_i_sb.setValue(self.settings['Substrate permittivity'].imag)
        # For partial incoherent case, set percentage variation of angles and thickness and the number of samples
        self.percentage_partial_incoherence_sb.setValue(self.settings['Percentage partial incoherence'])
        self.partially_incoherent_samples_sb.setValue(self.settings['Partially incoherent samples'])
        # For partial incoherent case, set the smoothing parameters
        self.partially_incoherent_kernel_sb.setValue(self.settings['Filter kernel size'])
        self.partially_incoherent_polynomial_sb.setValue(self.settings['Filter polynomial size'])
        # Redraw the layer information widget
        self.redrawLayerTable()
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        for i in range(20):
            QCoreApplication.processEvents()
        self.refreshRequired = False
        debugger.print(self.settings['Legend'],'Finished:: refresh, force =', force)
        return

    def setSuperstrate(self):
        ''' Use the superstrate material to define a superstrate and permittivity'''
        debugger.print(self.settings['Legend'],'setSuperstrate')
        if 'Material defined manually' in self.settings['Superstrate']:
            # Handle the manual case
            debugger.print(self.settings['Legend'],'set to manual')
            self.superstrateMaterial = Materials.Constant('manual',permittivity=self.settings['Superstrate permittivity'])
            self.superstratePermittivityFunction = self.superstrateMaterial.getPermittivityFunction()
        else:
            # check to see if the superstrate is available in the database
            if self.settings['Superstrate'] in self.materialNames:
                debugger.print(self.settings['Legend'],'set superstrate to',self.settings['Superstrate'])
                self.superstrateMaterial = self.DataBase.getMaterial(self.settings['Superstrate'])
                self.superstratePermittivityFunction = self.superstrateMaterial.getPermittivityFunction()
            else:
                print('Error: superstrate ',self.settings['Superstrate'],' not available in database')
                print('       available materials are:', self.materialNames)
                exit()
            materialPermittivityFunction = self.superstrateMaterial.getPermittivityFunction()
            materialPermittivity = materialPermittivityFunction(0.0)
            # Check to see if the permittivity has changed from the database values
            # if it has then force a 'Manual' setting.  This should allow backward compatibility
            if (np.abs(materialPermittivity - self.settings['Superstrate permittivity'])) > 1.0E-8:
                debugger.print(self.settings['Legend'],'setting superstrate to manual permittivity has changed')
                # It looks like the permittivity or density has been set manually
                # Creat a constant permittivity material and set the Superstrate to manual
                self.settings['Superstrate'] = 'Material defined manually'
                self.superstrateMaterial = Materials.Constant('manual',permittivity=self.settings['Superstrate permittivity'])
                self.superstratePermittivityFunction = self.superstrateMaterial.getPermittivityFunction()
            else:
                # Use the database material
                debugger.print(self.settings['Legend'],'setting superstrate to database material')
                self.superstratePermittivityFunction = self.superstrateMaterial.getPermittivityFunction()
        return

    def setSubstrate(self):
        ''' Use the substrate material to define a substrate and permittivity'''
        debugger.print(self.settings['Legend'],'setSubstrate')
        if 'Material defined manually' in self.settings['Substrate']:
            # Handle the manual case
            debugger.print(self.settings['Legend'],'set to manual')
            self.substrateMaterial = Materials.Constant('manual',permittivity=self.settings['Substrate permittivity'])
            self.substratePermittivityFunction = self.substrateMaterial.getPermittivityFunction()
        else:
            # check to see if the substrate is available in the database
            if self.settings['Substrate'] in self.materialNames:
                debugger.print(self.settings['Legend'],'set substrate to',self.settings['Substrate'])
                self.substrateMaterial = self.DataBase.getMaterial(self.settings['Substrate'])
                self.substratePermittivityFunction = self.substrateMaterial.getPermittivityFunction()
            else:
                print('Error: substrate ',self.settings['Substrate'],' not available in database')
                print('       available materials are:', self.materialNames)
                exit()
            materialPermittivityFunction = self.substrateMaterial.getPermittivityFunction()
            materialPermittivity = materialPermittivityFunction(0.0)
            # Check to see if the permittivity has changed from the database values
            # if it has then force a 'Manual' setting.  This should allow backward compatibility
            if (np.abs(materialPermittivity - self.settings['Substrate permittivity'])) > 1.0E-8:
                debugger.print(self.settings['Legend'],'setting superstrate to manual permittivity has changed')
                # It looks like the permittivity or density has been set manually
                # Creat a constant permittivity material and set the Substrate to manual
                self.settings['Substrate'] = 'Material defined manually'
                self.substrateMaterial = Materials.Constant('manual',permittivity=self.settings['Substrate permittivity'])
                self.substratePermittivityFunction = self.substrateMaterial.getPermittivityFunction()
            else:
                # Use the database material
                debugger.print(self.settings['Legend'],'setting superstrate to database material')
                self.substratePermittivityFunction = self.substrateMaterial.getPermittivityFunction()
        return

    def on_superstrate_cb_activated(self, index):
        '''Handle a change to the superstrate  material'''
        debugger.print(self.settings['Legend'],'on superstrate combobox activated', index)
        debugger.print(self.settings['Legend'],'on superstrate combobox activated', self.superstrate_cb.currentText())
        # We will need to recalculate everything for the new superstrate
        self.refreshRequired = True
        # superstrate is the name of the sheet in the database
        superstrate = self.superstrate_cb.currentText()
        # Make some of the widgets quiet as we update them
        m_blocking = self.superstrate_cb.signalsBlocked()
        r_blocking = self.superstrate_permittivity_r_sb.signalsBlocked()
        i_blocking = self.superstrate_permittivity_i_sb.signalsBlocked()
        self.superstrate_cb.blockSignals(True)
        self.superstrate_permittivity_r_sb.blockSignals(True)
        self.superstrate_permittivity_i_sb.blockSignals(True)
        # Store the new superstrate material
        if superstrate == 'Material defined manually':
            # The manual option has been chosen, so create a new material with the right permittivity
            self.superstrateMaterial = Materials.Constant('manual',permittivity=self.settings['Superstrate permittivity'])
            self.superstratePermittivityFunction = self.superstrateMaterial.getPermittivityFunction()
            self.settings['Superstrate'] = superstrate
        else:
            # Read the material information for permittivity and density from the data base
            superstrateMaterial = self.DataBase.getMaterial(superstrate)
            if superstrateMaterial.isScalar():
                # Only switch if it is a scalar material
                self.superstrateMaterial = superstrateMaterial
                self.settings['Superstrate'] = superstrate
            else:
                print('Error superstrate must have a scalar permittivity')
            self.superstratePermittivityFunction = self.superstrateMaterial.getPermittivityFunction()
            # The permittivity may be frequency dependent, show the value at 0 cm-1
            self.settings['Superstrate permittivity'] = self.superstratePermittivityFunction(0.0)
        # Update the superstrate material information
        text = self.superstrateMaterial.getInformation()
        self.superstrate_info_le.setText(text)
        # Update the values of the real and imaginary permittivity
        self.superstrate_permittivity_r_sb.setValue(self.settings['Superstrate permittivity'].real)
        self.superstrate_permittivity_i_sb.setValue(self.settings['Superstrate permittivity'].imag)
        # Restore the signal settings on the widgets
        self.superstrate_cb.blockSignals(m_blocking)
        self.superstrate_permittivity_r_sb.blockSignals(r_blocking)
        self.superstrate_permittivity_i_sb.blockSignals(i_blocking)
        self.refresh(force=True)
        self.refreshRequired = True
        return

    def on_substrate_cb_activated(self, index):
        '''Handle a change to the substrate  material'''
        debugger.print(self.settings['Legend'],'on substrate combobox activated', index)
        debugger.print(self.settings['Legend'],'on substrate combobox activated', self.substrate_cb.currentText())
        # We will need to recalculate everything for the new substrate
        self.refreshRequired = True
        # substrate is the name of the sheet in the database
        substrate = self.substrate_cb.currentText()
        # Make some of the widgets quiet as we update them
        m_blocking = self.substrate_cb.signalsBlocked()
        r_blocking = self.substrate_permittivity_r_sb.signalsBlocked()
        i_blocking = self.substrate_permittivity_i_sb.signalsBlocked()
        self.substrate_cb.blockSignals(True)
        self.substrate_permittivity_r_sb.blockSignals(True)
        self.substrate_permittivity_i_sb.blockSignals(True)
        # Store the new substrate material
        if 'Material defined manually' in substrate:
            # The manual option has been chosen, so create a new material with the right permittivity
            self.substrateMaterial = Materials.Constant('manual',permittivity=self.settings['Substrate permittivity'])
            self.substratePermittivityFunction = self.substrateMaterial.getPermittivityFunction()
            self.settings['Substrate'] = substrate
        else:
            # Read the material information for permittivity and density from the data base
            substrateMaterial = self.DataBase.getMaterial(substrate)
            if substrateMaterial.isScalar():
                # only switch if it is a scalar material
                self.substrateMaterial = substrateMaterial
                self.settings['Substrate'] = substrate
            else:
                print('Error substrate must have a scalar permittivity')
            self.substratePermittivityFunction = self.substrateMaterial.getPermittivityFunction()
            # The permittivity may be frequency dependent, show the value at 0 cm-1
            self.settings['Substrate permittivity'] = self.substratePermittivityFunction(0.0)
        # Update the substrate material information
        text = self.substrateMaterial.getInformation()
        self.substrate_info_le.setText(text)
        # Update the values of the real and imaginary permittivity
        self.substrate_permittivity_r_sb.setValue(np.real(self.settings['Substrate permittivity']))
        self.substrate_permittivity_i_sb.setValue(np.imag(self.settings['Substrate permittivity']))
        # Restore the signal settings on the widgets
        self.substrate_cb.blockSignals(m_blocking)
        self.substrate_permittivity_r_sb.blockSignals(r_blocking)
        self.substrate_permittivity_i_sb.blockSignals(i_blocking)
        self.refresh(force=True)
        self.refreshRequired = True
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
        self.generateLayerSettings()
        self.refresh(force=True)
        self.refreshRequired = True
        debugger.print(self.settings['Legend'],'Mode changed to ', self.settings['Mode'])
        debugger.print(self.settings['Legend'],'Finished:: on_mode_cb_activated')
        return

    def partially_incoherent_calculator( self,
                            superstrateDielectricFunction,
                            substrateDielectricFunction,
                            superstrateDepth,
                            substrateDepth,
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
        for s in range(self.settings["Partially incoherent samples"]):
            crystalDepth = d + d*( -1 + 2*np.random.rand())*self.settings['Percentage partial incoherence']/100.0
            crystalLayer.setThickness(crystalDepth)
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
                                       superstrateDepth,
                                       substrateDepth,
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
                            superstrateDielectricFunction,
                            substrateDielectricFunction,
                            superstrateDepth,
                            substrateDepth,
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
        crystalIncoherence = self.settings["Percentage partial incoherence"] * np.pi / 100.0 
        partial_function = partial(solve_single_crystal_equations,
                                       superstrateDielectricFunction,
                                       substrateDielectricFunction,
                                       superstrateDepth,
                                       substrateDepth,
                                       layers,
                                       crystalIncoherence,
                                       mode,
                                       theta,
                                       phi,
                                       psi,
                                       angleOfIncidence,
                                       sliceThickness, exponent_threshold)
        results = []
        # About to call
        debugger.print(self.settings['Legend'],'About to calculate single crystal scenario using pool')
        #begin serial version 
        #for v in self.vs_cm1:
        #   results.append(Calculator.solve_single_crystal_equations( superstrateDielectricFunction, substrateDielectricFunction, crystalPermittivityFunction, superstrateDepth, substrateDepth, crystalDepth, crystalIncoherence, mode, theta, phi, psi, angleOfIncidence,sliceThickness, exponent_threshold,v))
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
        for v,r,R,t,T,eps,errors,largest_exponent in results:
            if errors > 0:
                print('Warning exponential overflow occured at frequency',v,errors,largest_exponent)
            p_reflectance.append(R[0]+R[2])
            s_reflectance.append(R[1]+R[3])
            if self.settings['Mode'] == 'Thick slab':
                p_absorbtance.append(1.0 - R[0]-R[2])
                s_absorbtance.append(1.0 - R[1]-R[3])
                p_transmittance.append(np.zeros_like(T[0]))
                s_transmittance.append(np.zeros_like(T[1]))
#            elif self.settings['Mode'] == 'Incoherent thin film':
#                p_transmittance.append(T[0])
#                s_transmittance.append(T[1])
#                p_absorbtance.append(1.0 - p_reflectance[-1] - p_transmittance[-1])
#                s_absorbtance.append(1.0 - s_reflectance[-1] - s_transmittance[-1])
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
        # The dielectric variables are functions of frequency
        superstrateDielectric = self.settings['Superstrate permittivity']
        substrateDielectric   = self.settings['Substrate permittivity']
        superstrateDielectricFunction = DielectricFunction.ConstantScalar(superstrateDielectric).function()
        substrateDielectricFunction   = DielectricFunction.ConstantScalar(substrateDielectric).function()
        # make sure the units are correct for thickness
        thicknessUnits = {'nm':1.0E-9, 'um':1.0E-6, 'mm':1.0E-3, 'cm':1.0E-2}
        tom = thicknessUnits[self.settings['Superstrate & substrate thickness unit']]
        superstrateDepth = tom * self.settings['Superstrate depth']
        superstrateDepth = tom * self.settings['Superstrate depth']
        substrateDepth   = tom * self.settings['Substrate depth']
        tom = thicknessUnits[self.settings['Slice thickness unit']]
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
        if mode == 'Thick slab' or mode == 'Coherent thin film' or mode == 'Incoherent thin film':
            calculator = self.coherent_calculator
        else:
            calculator = self.partially_incoherent_calculator
        exponent_threshold = self.exponent_threshold
        # Call the relevant calculator
        ( self.p_reflectance, 
        self.s_reflectance, 
        self.p_transmittance, 
        self.s_transmittance, 
        self.s_absorbtance, 
        self.p_absorbtance, 
        self.epsilon) = calculator( superstrateDielectricFunction,
                                    substrateDielectricFunction,
                                    superstrateDepth,
                                    substrateDepth,
                                    self.layers,
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
        # If the single crystal mode is Thick Slab, there is no need for substrate permittivity
        debugger.print(self.settings['Legend'],'Start:: greyed_out')
        if self.settings['Mode'] == 'Thick slab':
            self.substrate_cb.setEnabled(False)
            self.substrate_info_le.setEnabled(False)
            self.substrate_permittivity_r_sb.setEnabled(False)
            self.substrate_permittivity_i_sb.setEnabled(False)
            self.partially_incoherent_samples_sb.setEnabled(False)
            self.percentage_partial_incoherence_sb.setEnabled(False)
            self.partially_incoherent_polynomial_sb.setEnabled(False)
            self.partially_incoherent_kernel_sb.setEnabled(False)
        elif self.settings['Mode'] == 'Coherent thin film':
            self.substrate_cb.setEnabled(True)
            self.substrate_info_le.setEnabled(True)
            self.substrate_permittivity_r_sb.setEnabled(True)
            self.substrate_permittivity_i_sb.setEnabled(True)
            self.partially_incoherent_samples_sb.setEnabled(False)
            self.percentage_partial_incoherence_sb.setEnabled(False)
            self.partially_incoherent_polynomial_sb.setEnabled(False)
            self.partially_incoherent_kernel_sb.setEnabled(False)
        elif self.settings['Mode'] == 'Partially incoherent thin film':
            self.substrate_cb.setEnabled(True)
            self.substrate_info_le.setEnabled(True)
            self.substrate_permittivity_r_sb.setEnabled(True)
            self.substrate_permittivity_i_sb.setEnabled(True)
            self.partially_incoherent_samples_sb.setEnabled(True)
            self.percentage_partial_incoherence_sb.setEnabled(True)
            self.partially_incoherent_polynomial_sb.setEnabled(True)
            self.partially_incoherent_kernel_sb.setEnabled(True)
        elif self.settings['Mode'] == 'Incoherent thin film':
            self.substrate_cb.setEnabled(True)
            self.substrate_info_le.setEnabled(True)
            self.substrate_permittivity_r_sb.setEnabled(True)
            self.substrate_permittivity_i_sb.setEnabled(True)
            self.partially_incoherent_samples_sb.setEnabled(False)
            self.percentage_partial_incoherence_sb.setEnabled(False)
            self.partially_incoherent_polynomial_sb.setEnabled(False)
            self.partially_incoherent_kernel_sb.setEnabled(False)
        debugger.print(self.settings['Legend'],'Finished:: greyed_out')

