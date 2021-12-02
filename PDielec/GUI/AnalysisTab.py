import math
import numpy as np
import PDielec.Calculator as Calculator
from PyQt5.QtWidgets  import  QWidget, QApplication
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox, QDoubleSpinBox
from PyQt5.QtWidgets  import  QSizePolicy, QTableWidgetItem
from PyQt5.QtCore     import  Qt, QCoreApplication
from PDielec.Constants import  covalent_radii
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PDielec.Utilities import Debug
from PDielec.GUI.SettingsTab import FixedQTableWidget

class AnalysisTab(QWidget):
    def __init__(self, parent, debug=False ):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'AnalysisTab')
        self.settings = {}
        self.subplot = None
        self.setWindowTitle('Analysis')
        self.settings['Minimum frequency'] = -1
        self.settings['Maximum frequency'] = 400
        self.settings['title'] = 'Analysis'
        self.settings['Covalent radius scaling'] = 1.1
        self.settings['Bonding tolerance'] = 0.1
        self.settings['Bar width'] = 0.5
        self.refreshRequired = True
        self.plot_types = ['Internal vs External','Molecular Composition']
        self.plot_type_index = 0
        self.number_of_molecules = 0
        self.frequency_units = None
        self.cell_of_molecules = None
        self.original_atomic_order = None
        self.frequencies_cm1 = []
        self.mode_energies = []
        self.element_radii = covalent_radii
        self.species = []
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Create last tab - AnalysisTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The minimum frequency
        #
        self.vmin_sb = QDoubleSpinBox(self)
        self.vmin_sb.setRange(-100,9000)
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmin_sb.setToolTip('Set the minimum frequency to be considered)')
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        label = QLabel('Minimum frequency:', self)
        label.setToolTip('Set the minimum frequency to be considered)')
        form.addRow(label, self.vmin_sb)
        #
        # The maximum frequency
        #
        self.vmax_sb = QDoubleSpinBox(self)
        self.vmax_sb.setRange(0,9000)
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.vmax_sb.setToolTip('Set the maximum frequency to be considered)')
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        label = QLabel('Maximum frequency:', self)
        label.setToolTip('Set the maximum frequency to be considered)')
        form.addRow(label, self.vmax_sb)
        #
        # The bonding tolerance and scaling of radii frequency
        #
        hbox = QHBoxLayout()
        self.scale_sp = QDoubleSpinBox(self)
        self.scale_sp.setRange(0.01,10.0)
        self.scale_sp.setSingleStep(0.01)
        self.scale_sp.setDecimals(2)
        self.scale_sp.setValue(self.settings['Covalent radius scaling'])
        self.scale_sp.setToolTip('Scale the covalent radii to determine bonding')
        self.scale_sp.valueChanged.connect(self.on_scale_changed)
        hbox.addWidget(self.scale_sp)
        self.tolerance_sp = QDoubleSpinBox(self)
        self.tolerance_sp.setRange(0.01,2.0)
        self.tolerance_sp.setSingleStep(0.01)
        self.tolerance_sp.setDecimals(2)
        self.tolerance_sp.setValue(self.settings['Bonding tolerance'])
        self.tolerance_sp.setToolTip('Tolerance for bonding is determined from scale*(radi+radj)+toler')
        self.tolerance_sp.valueChanged.connect(self.on_tolerance_changed)
        hbox.addWidget(self.tolerance_sp)
        label = QLabel('Bonding scale and tolerance', self)
        label.setToolTip('Bonding is determined from scale*(radi+radj)+toler')
        form.addRow(label, hbox)
        # Add a table of covalent radii
        self.element_radii_tw = FixedQTableWidget(parent=self)
        self.element_radii_tw.setToolTip('Individual covalent radii used to determine bonding can be set here')
        self.element_radii_tw.itemClicked.connect(self.on_element_radii_tw_itemClicked)
        self.element_radii_tw.itemChanged.connect(self.on_element_radii_tw_itemChanged)
        self.element_radii_tw.setRowCount(1)
        self.element_radii_tw.blockSignals(False)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        self.element_radii_tw.setSizePolicy(sizePolicy)
        form.addRow(QLabel('Atomic radii', self), self.element_radii_tw)
        # Add number of molecules found
        self.molecules_le = QLineEdit(self)
        self.molecules_le.setEnabled(False)
        self.molecules_le.setText('{}'.format(self.number_of_molecules))
        self.molecules_le.setToolTip('The bonding tolerances can change the number of molecules found')
        label = QLabel('Number of molecules found', self)
        label.setToolTip('The bonding tolerances can change the number of molecules found')
        form.addRow(label, self.molecules_le)
        #
        # The plotting width of bar
        #
        self.width_sp = QDoubleSpinBox(self)
        self.width_sp.setRange(0.01,2.0)
        self.width_sp.setSingleStep(0.01)
        self.width_sp.setDecimals(2)
        self.width_sp.setValue(self.settings['Bar width'])
        self.width_sp.setToolTip('Change the width of the bars - should be between 0 and 1')
        self.width_sp.valueChanged.connect(self.on_width_changed)
        label = QLabel('Bar width', self)
        label.setToolTip('Change the width of the bars - should be between 0 and 1')
        form.addRow(label, self.width_sp)
        #
        # Set the plot title
        #
        self.title_le = QLineEdit(self)
        self.title_le.setToolTip('Set the plot title')
        self.title_le.setText(self.settings['title'])
        self.title_le.textChanged.connect(self.on_title_changed)
        label = QLabel('Plot title', self)
        label.setToolTip('Set the plot title')
        form.addRow(label, self.title_le)
        #
        # Add a comb box to select which type of plot
        #
        self.plottype_cb = QComboBox(self)
        self.plottype_cb.setToolTip('The energy can be decomposed either according to internal vs external motion or into a molecular based decompostion')
        self.plottype_cb.addItems(self.plot_types)
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        self.plottype_cb.currentIndexChanged.connect(self.on_plottype_cb_changed)
        label = QLabel('Choose the plot type', self)
        label.setToolTip('The energy can be decomposed either according to internal vs external motion or into a molecular based decompostion')
        form.addRow(label, self.plottype_cb)
        #
        # Add the matplotlib figure to the bottom
        #
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding))
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        form.addRow(self.canvas)
        form.addRow(self.toolbar)
        vbox.addLayout(form)
        # finalise the layout
        self.setLayout(vbox)
        QCoreApplication.processEvents()
        #if self.notebook.spreadsheet is not None:
        #    self.writeSpreadsheet()
        #QCoreApplication.processEvents()

    def on_element_radii_tw_itemClicked(self,item):
        self.element_radii_tw.blockSignals(False)

    def on_element_radii_tw_itemChanged(self,item):
        if self.reader is None:
            return
        col = item.column()
        try:
            debugger.print('Changing the element radius',col,item.text())
            self.element_radii[self.species[col]] = float(item.text())
            self.calculate()
            self.plot()
            if self.notebook.viewerTab is not None:
                self.notebook.viewerTab.requestRefresh()
        except:
            debugger.print('Failed Changing the element radius',col,item.text())
            pass

    def set_radii_tw(self):
        if self.reader is None:
            return
        self.element_radii_tw.blockSignals(True)
        self.species = self.reader.getSpecies()
        radii = [ self.element_radii[el] for el in self.species ]
        self.element_radii_tw.setColumnCount(len(self.species))
        self.element_radii_tw.setHorizontalHeaderLabels(self.species)
        self.element_radii_tw.setVerticalHeaderLabels([''])
        for i,(radius,element) in enumerate(zip(radii,self.species)):
            qw = QTableWidgetItem()
            qw.setText('{0:.6f}'.format(radius))
            qw.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.element_radii_tw.setItem(0,i, qw )
        self.element_radii_tw.blockSignals(False)
        # end if
        return

    def setCovalentRadius(self,element,radius):
        self.element_radii[element] = radius
        self.set_radii_tw()
        self.calculate()
        self.plot()

    def writeSpreadsheet(self):
        if self.notebook.spreadsheet is None:
            return
        sp = self.notebook.spreadsheet
        sp.selectWorkSheet('Analysis')
        sp.delete()
        sp.writeNextRow(['Analysis of the vibrational modes into percentage contributions for molecules and internal/external modes'], row=0,col=1)
        headers = ['Mode','Frequency (cm-1)', 'Centre of mass %','Rotational %', 'Vibrational %']
        for mol in range(self.number_of_molecules):
            headers.append('Molecule '+str(mol)+' %')
        #
        sp.writeNextRow(headers,col=1)
        for imode,(freq,energies) in enumerate(zip(self.frequencies_cm1,self.mode_energies)):
           tote,cme,rote,vibe, molecular_energies = energies
           tote = max(tote,1.0E-8)
           output = [ imode+1, freq, 100*cme/tote, 100*rote/tote, 100*vibe/tote ]
           for e in molecular_energies:
               output.append(100*e/tote)
           sp.writeNextRow(output,col=1,check=1)


    def on_width_changed(self,value):
        debugger.print('on width changed ', value)
        self.settings['Bar width'] = value
        self.plot()

    def on_scale_changed(self,value):
        debugger.print('on scale_le changed ', value)
        self.settings['Covalent radius scaling'] = value
        self.refreshRequired = True
        self.calculate()
        self.plot()

    def on_tolerance_changed(self,value):
        debugger.print('on_tolerance_le changed ', value)
        self.settings['Bonding tolerance'] = value
        self.refreshRequired = True
        self.calculate()
        self.plot()

    def on_title_changed(self,text):
        self.settings['title'] = text
        if self.subplot is not None:
            self.subplot.set_title(self.settings['title'])
            self.canvas.draw_idle()
        debugger.print('on title change ', self.settings['title'])

    def on_vmin_changed(self):
        self.vmin_sb.blockSignals(True)
        vmin = self.vmin_sb.value()
        vmax = self.vmax_sb.value()
        if vmin < vmax:
            self.settings['Minimum frequency'] = vmin
            debugger.print('on_vmin_changed new value', self.settings['Minimum frequency'])
        else:
            self.vmin_sb.setValue(self.settings['Maximum frequency']-1)
            self.settings['Minimum frequency'] = self.settings['Maximum frequency']-1
            debugger.print('on_vmin_changed restricting value to', self.settings['Minimum frequency'])
        self.plot()
        self.vmin_sb.blockSignals(False)
        return

    def on_vmax_changed(self):
        self.vmin_sb.blockSignals(True)
        vmin = self.vmin_sb.value()
        vmax = self.vmax_sb.value()
        if vmax > vmin:
            self.settings['Maximum frequency'] = vmax
            debugger.print('on_vmax_changed new value', self.settings['Maximum frequency'])
        else:
            self.vmin_sb.setValue(self.settings['Minimum frequency']+1)
            self.settings['Maximum frequency'] = self.settings['Minimum frequency']+1
            debugger.print('on_vmax_changed restricting value to', self.settings['Maximum frequency'])
        self.plot()
        self.vmin_sb.blockSignals(False)
        return

    def requestRefresh(self):
        self.refreshRequired = True

    def refresh(self, force=False):
        if not self.refreshRequired and not force:
            debugger.print('return with no refresh', self.refreshRequired, force)
            return
        debugger.print('Refreshing widget')
        #
        # Block signals during refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.scale_sp.setValue(self.settings['Covalent radius scaling'])
        self.tolerance_sp.setValue(self.settings['Bonding tolerance'])
        self.molecules_le.setText('{}'.format(self.number_of_molecules))
        self.width_sp.setValue(self.settings['Bar width'])
        self.title_le.setText(self.settings['title'])
        self.plottype_cb.setCurrentIndex(self.plot_type_index)
        self.calculate()
        self.set_radii_tw()
        self.plot()
        #
        # Unlock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        return

    def on_plottype_cb_changed(self, index):
        self.plot_type_index = index
        debugger.print('Plot type index changed to ', self.plot_type_index)
        self.plot()

    def calculate(self):
        debugger.print('calculate')
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['Program']
        filename = self.notebook.mainTab.getFullFileName()
        self.reader = self.notebook.mainTab.reader
        if self.reader is None:
            return
        if program == '':
            return
        if filename == '':
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        mass_weighted_normal_modes = self.notebook.settingsTab.mass_weighted_normal_modes
        scale = self.settings['Covalent radius scaling']
        tolerance = self.settings['Bonding tolerance']
        # Find the last unit cell read by the reader and its masses
        cell = self.reader.unit_cells[-1]
        atom_masses = self.reader.masses
        cell.set_atomic_masses(atom_masses)
        self.cell_of_molecules,nmols,self.original_atomic_order = cell.calculate_molecular_contents(scale, tolerance, covalent_radii)
        # if the number of molecules has changed then tell the viewerTab that the cell has changed
        if self.number_of_molecules != nmols:
            if self.notebook.viewerTab is not None:
                self.notebook.viewerTab.requestRefresh()
            self.number_of_molecules = nmols
        self.molecules_le.setText('{}'.format(self.number_of_molecules))
        # get the normal modes from the mass weighted ones
        normal_modes = Calculator.normal_modes(atom_masses, mass_weighted_normal_modes)
        # Reorder the atoms so that the mass weighted normal modes order agrees with the ordering in the cell_of_molecules cell
        nmodes,nions,temp = np.shape(normal_modes)
        self.new_normal_modes = np.zeros( (nmodes,3*nions) )
        self.new_mass_weighted_normal_modes = np.zeros( (nmodes,3*nions) )
        masses = self.cell_of_molecules.atomic_masses
        for imode,mode in enumerate(mass_weighted_normal_modes):
            for index,old_index in enumerate(self.original_atomic_order):
                i = index*3
                j = old_index*3
                self.new_mass_weighted_normal_modes[imode,i+0] = mode[old_index][0]
                self.new_mass_weighted_normal_modes[imode,i+1] = mode[old_index][1]
                self.new_mass_weighted_normal_modes[imode,i+2] = mode[old_index][2]
                self.new_normal_modes[imode,i+0] = self.new_mass_weighted_normal_modes[imode,i+0] / math.sqrt(masses[index])
                self.new_normal_modes[imode,i+1] = self.new_mass_weighted_normal_modes[imode,i+1] / math.sqrt(masses[index])
                self.new_normal_modes[imode,i+2] = self.new_mass_weighted_normal_modes[imode,i+2] / math.sqrt(masses[index])
        # Calculate the distribution in energy for the normal modes
        mode_energies = Calculator.calculate_energy_distribution(self.cell_of_molecules, self.frequencies_cm1, self.new_mass_weighted_normal_modes)
        # Deal with degeneracies
        degenerate_list = [ [] for f in self.frequencies_cm1]
        for i,fi in enumerate(self.frequencies_cm1):
            for j,fj in enumerate(self.frequencies_cm1):
                if abs(fi-fj) < 1.0E-5:
                    degenerate_list[i].append(j)
        self.mode_energies = []
        for i,fi in enumerate(self.frequencies_cm1):
            tote,cme,rote,vibe,mole = mode_energies[i]
            tote = max(tote,1.0E-8)
            sums = [0.0]*5
            sume = [0.0]*len(mole)
            degeneracy = len(degenerate_list[i])
            for j in degenerate_list[i]:
                tote,cme,rote,vibe,mole = mode_energies[j]
                tote = max(tote,1.0E-8)
                sums[0] += tote / degeneracy
                sums[1] += cme / degeneracy
                sums[2] += rote / degeneracy
                sums[3] += vibe / degeneracy
                for k,e in enumerate(mole):
                    sume[k] += e / degeneracy
                sums[4] = sume
            self.mode_energies.append(sums)
        # Store the results in the spread shee
        # if self.notebook.spreadsheet is not None:
        #     self.writeSpreadsheet()
        # Flag that a recalculation is not needed
        self.refreshRequired = False
        QApplication.restoreOverrideCursor()

    def plot(self):
         if self.reader is None:
             return
         if len(self.frequencies_cm1) <= 0:
             return
         if self.plot_type_index == 0:
             self.plot_internal_external()
         else:
             self.plot_molecular()

    def plot_molecular(self):
        self.subplot = None
        self.figure.clf()
        self.subplot = self.figure.add_subplot(111)
        self.subplot.xaxis.set_major_locator(MaxNLocator(integer=True))
        xlabel = 'Mode Number'
        ylabel = 'Percentage energy'
        # Decide which modes to analyse
        mode_list = []
        mode_list_text = []
        vmin = self.settings['Minimum frequency']
        vmax = self.settings['Maximum frequency']
        tote,cme,rote,vibe,mole = self.mode_energies[0]
        tote = max(tote,1.0E-8)
        mol_energies = [ [] for _ in range(self.number_of_molecules) ]
        mol_bottoms  = [ [] for _ in range(self.number_of_molecules) ]
        for imode, frequency in enumerate(self.frequencies_cm1):
            if frequency >= vmin and frequency <= vmax:
                mode_list.append(imode+1)
                mode_list_text.append(str(imode+1))
                tote,cme,rote,vibe,mole = self.mode_energies[imode]
                tote = max(tote,1.0E-8)
                for i,mol in enumerate(mole):
                    mol_energies[i].append(100.0*mol/tote)
                    if i == 0:
                        mol_bottoms[i].append(0.0)
                    else:
                        mol_bottoms[i].append(mol_bottoms[i-1][-1]+mol_energies[i-1][-1])
        if len(mode_list) < 3:
            return
        width = self.settings['Bar width']
        plots = []
        colours = ['y','b','r','c','m','k']
        for i,(energies,bottoms) in enumerate(zip(mol_energies, mol_bottoms )):
            plots.append(self.subplot.bar(mode_list,energies,width, bottom=bottoms,color=colours[i%6]))
        legends = []
        for i in range(self.number_of_molecules):
            legends.append('Molecule '+str(i))
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        self.subplot.legend( plots, legends)
        self.subplot.set_title('Molecular Composition of Vibrational Energy')
        self.canvas.draw_idle()

    def plot_internal_external(self):
        self.subplot = None
        self.figure.clf()
        self.subplot = self.figure.add_subplot(111)
        self.subplot.xaxis.set_major_locator(MaxNLocator(integer=True))
        xlabel = 'Mode Number'
        ylabel = 'Percentage energy'
        # Decide which modes to analyse
        mode_list = []
        mode_list_text = []
        cme_energy = []
        rot_energy = []
        vib_energy = []
        vib_bottom = []
        mol_energy = []
        vmin = self.settings['Minimum frequency']
        vmax = self.settings['Maximum frequency']
        colours = ['y','b','r','c','m','k']
        for imode, frequency in enumerate(self.frequencies_cm1):
            if frequency >= vmin and frequency <= vmax:
                mode_list.append(imode+1)
                mode_list_text.append(str(imode+1))
                tote,cme,rote,vibe,molecular_energies = self.mode_energies[imode]
                molecular_energies = np.array(molecular_energies)
                tote = max(tote,1.0E-8)
                cme_energy.append(cme/tote*100.0)
                rot_energy.append(rote/tote*100.0)
                vib_energy.append(vibe/tote*100.0)
                vib_bottom.append( (cme+rote)/tote*100.0 )
                mol_energy.append(molecular_energies/tote*100)
        if len(mode_list) < 3:
            return
        width = self.settings['Bar width']
        p1 = self.subplot.bar(mode_list,cme_energy,width,color=colours[0])
        p2 = self.subplot.bar(mode_list,rot_energy,width,bottom=cme_energy,color=colours[1])
        p3 = self.subplot.bar(mode_list,vib_energy,width,bottom=vib_bottom,color=colours[2])
        plots = ( p1[0], p2[0], p3[0] )
        legends = ('translation','rotation','vibration')
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        self.subplot.legend( plots, legends)
        self.subplot.set_title('Internal-External Composition of Vibrational Energy')
        self.canvas.draw_idle()

