import math
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox
from PyQt5.QtWidgets  import  QSizePolicy
from PyQt5.QtCore     import  Qt
from Python.Constants import  wavenumber, amu, PI, avogadro_si, angstrom
from Python.Constants import  covalent_radii
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from Python.Utilities import Debug
from Python.Plotter import print_strings, print_reals

class AnalysisTab(QWidget):
    def __init__(self, parent, debug=False ):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'AnalysisTab')
        self.settings = {}
        self.subplot = None
        self.setWindowTitle('Analysis')
        self.settings['vmin'] = 0
        self.settings['vmax'] = 400
        self.settings['title'] = 'Analysis'
        self.settings['bond_scaling'] = 1.1
        self.settings['bond_tolerance'] = 0.1
        self.settings['bar_width'] = 0.5
        self.settings['plot_types'] = ['Internal vs External','Molecular Composition']
        self.plot_type_index = 0
        self.number_of_molecules = 0
        self.frequency_units = None
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
        self.vmin_sb = QSpinBox(self)
        self.vmin_sb.setRange(-100,9000)
        self.vmin_sb.setValue(self.settings['vmin'])
        self.vmin_sb.setToolTip('Set the minimum frequency to be considered)')
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        label = QLabel('Minimum frequency:', self)
        label.setToolTip('Set the minimum frequency to be considered)')
        form.addRow(label, self.vmin_sb)
        #
        # The maximum frequency
        #
        self.vmax_sb = QSpinBox(self)
        self.vmax_sb.setRange(0,9000)
        self.vmax_sb.setValue(self.settings['vmax'])
        self.vmax_sb.setToolTip('Set the maximum frequency to be considered)')
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        label = QLabel('Maximum frequency:', self)
        label.setToolTip('Set the maximum frequency to be considered)')
        form.addRow(label, self.vmax_sb)
        #
        # The bond tolerance and scaling of radii frequency
        #
        hbox = QHBoxLayout()
        self.scale_le = QLineEdit(self)
        self.scale_le.setText('{}'.format(self.settings['bond_scaling']))
        self.scale_le.setToolTip('Scale the covalent radii to determine bonding')
        self.scale_le.textChanged.connect(self.on_scale_changed)
        hbox.addWidget(self.scale_le)
        self.tolerance_le = QLineEdit(self)
        self.tolerance_le.setText('{}'.format(self.settings['bond_tolerance']))
        self.tolerance_le.setToolTip('Tolerance for bonding is determined from scale*(radi+radj)+toler')
        self.tolerance_le.textChanged.connect(self.on_tolerance_changed)
        hbox.addWidget(self.tolerance_le)
        label = QLabel('Bonding scale and tolerance', self)
        label.setToolTip('Bonding is determined from scale*(radi+radj)+toler')
        form.addRow(label, hbox)
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
        self.width_le = QLineEdit(self)
        self.width_le.setText('{}'.format(self.settings['bar_width']))
        self.width_le.setToolTip('Change the width of the bars - should be between 0 and 1')
        self.width_le.textChanged.connect(self.on_width_changed)
        label = QLabel('Bar width', self)
        label.setToolTip('Change the width of the bars - should be between 0 and 1')
        form.addRow(label, self.width_le)
        # 
        # Store results in a file?
        #
        #self.file_store_le = QLineEdit(self) 
        #self.file_store_le.setToolTip('Store the results in a .csv or .xlsx file')
        #self.file_store_le.setText(self.settings['spreadsheet'])
        #self.file_store_le.textChanged.connect(self.on_file_store_le_changed)
        #form.addRow(QLabel('Output spreadsheet', self), self.file_store_le)
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
        for choice in self.settings["plot_types"]:
            self.plottype_cb.addItem(choice)
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

    def on_width_changed(self,text):
        debugger.print('on width changed ', text)
        try:
          self.settings['bar_width'] = float(text)
        except:
          pass
        self.width_le.blockSignals(True)
        if self.settings['bar_width'] < 0.0:
            self.settings['bar_width'] = 0.0
            self.width_le.setText('{}'.format(self.settings['bar_width']))
        if self.settings['bar_width'] > 1.0:
            self.settings['bar_width'] = 1.0
            self.width_le.setText('{}'.format(self.settings['bar_width']))
        self.width_le.blockSignals(False)
        self.plot()
        
    def on_scale_changed(self,text):
        debugger.print('on scale_le changed ', text)
        try:
            self.settings['bond_scaling'] = float(text)
        except:
          pass
        self.width_le.blockSignals(True)
        if self.settings['bond_scaling'] < 0.0:
            self.settings['bond_scaling'] = 0.0
            self.scale_le.setText('{}'.format(self.settings['bond_scaling']))
        if self.settings['bond_scaling'] > 2.0:
            self.settings['bond_scaling'] = 2.0
            self.scale_le.setText('{}'.format(self.settings['bond_scaling']))
        self.width_le.blockSignals(False)
        self.refresh()
        
    def on_tolerance_changed(self,text):
        debugger.print('on file_tolerance_le changed ', text)
        try:
            self.settings['bond_tolerance'] = float(text)
        except:
          pass
        self.width_le.blockSignals(True)
        if self.settings['bond_tolerance'] < 0.0:
            self.settings['bond_tolerance'] = 0.0
            self.tolerance_le.setText('{}'.format(self.settings['bond_tolerance']))
        if self.settings['bond_tolerance'] > 2.0:
            self.settings['bond_tolerance'] = 2.0
            self.tolerance_le.setText('{}'.format(self.settings['bond_tolerance']))
        self.width_le.blockSignals(False)
        self.refresh()

    def on_file_store_le_changed(self,text):
        self.settings['spreadsheet'] = text
        debugger.print('on file_store_le change ', self.settings['spreadsheet'])

    def on_title_changed(self,text):
        self.settings['title'] = text
        if self.subplot is not None:
            self.subplot.set_title(self.settings['title'])
            self.canvas.draw_idle()
        debugger.print('on title change ', self.settings['title'])

    def on_vmin_changed(self):
        self.settings['vmin'] = self.vmin_sb.value()
        debugger.print('on vmin change ', self.settings['vmin'])
        vmin = self.settings['vmin']
        vmax = self.settings['vmax']
        if vmax > vmin:
            self.plot()

    def on_vmax_changed(self):
        self.settings['vmax'] = self.vmax_sb.value()
        debugger.print('on vmax change ', self.settings['vmax'])
        vmin = self.settings['vmin']
        vmax = self.settings['vmax']
        if vmax > vmin:
            self.plot()

    def refresh(self):
        debugger.print('refreshing widget')
        self.calculate()
        self.plot()
        return

    def on_plottype_cb_changed(self, index):
        self.plot_type_index = index
        debugger.print('Plot type index changed to ', self.plot_type_index)
        self.plot()

    def calculate(self):
        debugger.print('calculate')
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['program']
        filename = settings['filename']
        self.reader = self.notebook.mainTab.reader
        if self.reader is None:
            return
        if program is '':
            return
        if filename is '':
            return
        # Assemble the settingsTab settings
        settings = self.notebook.settingsTab.settings
        eckart = settings['eckart']
        neutral = settings['neutral']
        hessian_symm = settings['hessian_symmetrisation']
        epsilon_inf = np.array(settings['optical_permittivity'])
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        mass_weighted_normal_modes = self.notebook.settingsTab.mass_weighted_normal_modes
        frequencies = np.array(self.frequencies_cm1) * wavenumber
        intensities = self.notebook.settingsTab.intensities
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = self.reader.volume*angstrom*angstrom*angstrom
        vmin = self.settings['vmin']
        vmax = self.settings['vmax']
        scale = self.settings['bond_scaling']
        tolerance = self.settings['bond_tolerance']
        # Find the last unit cell read by the reader and its masses
        cell = self.reader.unit_cells[-1]
        atom_masses = self.reader.masses
        cell.set_atomic_masses(atom_masses)
        newcell,nmols,old_order = cell.calculate_molecular_contents(scale, tolerance, covalent_radii)
        #newcell.printInfo()
        self.number_of_molecules = nmols
        self.molecules_le.setText('{}'.format(self.number_of_molecules))
        # get the normal modes from the mass weighted ones
        normal_modes = Calculator.normal_modes(atom_masses, mass_weighted_normal_modes)
        # Reorder the atoms so that the mass weighted normal modes order agrees with the ordering in the new cell
        nmodes,nions,temp = np.shape(normal_modes)
        new_normal_modes = np.zeros( (nmodes,3*nions) )
        new_mass_weighted_normal_modes = np.zeros( (nmodes,3*nions) )
        masses = newcell.atomic_masses
        for imode,mode in enumerate(mass_weighted_normal_modes):
            for index,old_index in enumerate(old_order):
                i = index*3
                j = old_index*3
                new_mass_weighted_normal_modes[imode,i+0] = mode[old_index][0] 
                new_mass_weighted_normal_modes[imode,i+1] = mode[old_index][1] 
                new_mass_weighted_normal_modes[imode,i+2] = mode[old_index][2] 
                new_normal_modes[imode,i+0] = new_mass_weighted_normal_modes[imode,i+0] / math.sqrt(masses[index])
                new_normal_modes[imode,i+1] = new_mass_weighted_normal_modes[imode,i+1] / math.sqrt(masses[index])
                new_normal_modes[imode,i+2] = new_mass_weighted_normal_modes[imode,i+2] / math.sqrt(masses[index])
        # Calculate the distribution in energy for the normal modes
        self.mode_energies = Calculator.calculate_energy_distribution(newcell, self.frequencies_cm1, new_mass_weighted_normal_modes)
        # Output the final result
        #title = ['Freq(cm-1)','%mol-cme','%mol-rot','%vib']
        #for i in range(self.number_of_molecules):
        #    title.append('%mol-'+str(i))
        #print_strings('Percentage energies in vibrational modes',title,format="{:>10}")
        fd_csvfile = None
        #if not fd_csvfile is None:
        #    print_strings('Percentage energies in vibrational modes',title,format="{:>10}",file=fd_csvfile,separator=",")
        #for freq,energies in zip(self.frequencies_cm1,self.mode_energies):
        #       tote,cme,rote,vibe,molecular_energies = energies
        #       output = [ freq, 100*cme/tote, 100*rote/tote, 100*vibe/tote ]
        #       for e in molecular_energies:
        #           output.append(100*e/tote)
        #       print_reals('',output,format="{:10.2f}")
        #       if not fd_csvfile is None:
        #           print_reals('',output,format="{:10.2f}",file=fd_csvfile,separator=",")
        # if using a excel file write out the results
        fd_excelfile = None
        if not fd_excelfile is None:
            excel_row += 2; worksheet.write(excel_row,0,'Percentage energies in vibrational modes')
            excel_row += 1; [ worksheet.write(excel_row,col,f) for col,f in enumerate(title) ]
            for freq,energies in zip(self.frequencies_cm1,self.mode_energies):
               tote,cme,rote,vibe, molecular_energies = energies
               output = [ freq, 100*cme/tote, 100*rote/tote, 100*vibe/tote ]
               for e in molecular_energies:
                   output.append(100*e/tote)
               excel_row += 1; [ worksheet.write(excel_row,col,f) for col,f in enumerate( output ) ]
        if fd_excelfile is not None:
            workbook.close()
        if fd_csvfile is not None:
            fd_csvfile.close()
        #

    def plot(self):
         if self.plot_type_index == 0:
             self.plot_internal_external()
         else:
             self.plot_molecular()

    def plot_molecular(self):
        self.subplot = None
        self.figure.clf()
        self.subplot = self.figure.add_subplot(111)
        xlabel = 'Molecule'
        ylabel = 'Percentage energy'
        # Decide which modes to analyse
        mode_list = []
        mode_list_text = []
        vmin = self.settings['vmin']
        vmax = self.settings['vmax']
        tote,cme,rote,vibe,mole = self.mode_energies[0]
        mol_energies = None
        mol_energies = [ [] for _ in range(self.number_of_molecules) ]
        mol_bottoms  = [ [] for _ in range(self.number_of_molecules) ]
        for imode, frequency in enumerate(self.frequencies_cm1):
            if frequency >= vmin and frequency <= vmax:
                mode_list.append(imode)
                mode_list_text.append(str(imode))
                tote,cme,rote,vibe,mole = self.mode_energies[imode]
                for i,mol in enumerate(mole):
                    mol_energies[i].append(100.0*mol/tote)
                    if i == 0:
                        mol_bottoms[i].append(0.0)
                    else:
                        mol_bottoms[i].append(mol_bottoms[i-1][-1]+mol_energies[i-1][-1])
        width = self.settings['bar_width']
        plots = []
        for i,(energies,bottoms) in enumerate(zip(mol_energies, mol_bottoms )):
            plots.append(self.subplot.bar(mode_list,energies,width, bottom=bottoms))
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
        vmin = self.settings['vmin']
        vmax = self.settings['vmax']
        for imode, frequency in enumerate(self.frequencies_cm1):
            if frequency >= vmin and frequency <= vmax:
                mode_list.append(imode)
                mode_list_text.append(str(imode))
                tote,cme,rote,vibe,molecular_energies = self.mode_energies[imode]
                cme_energy.append(cme/tote*100.0)
                rot_energy.append(rote/tote*100.0)
                vib_energy.append(vibe/tote*100.0)
                vib_bottom.append( (cme+rote)/tote*100.0 )
                mol_energy.append(molecular_energies/tote*100)
        width = self.settings['bar_width']
        p1 = self.subplot.bar(mode_list,cme_energy,width)
        p2 = self.subplot.bar(mode_list,rot_energy,width,bottom=cme_energy)
        p3 = self.subplot.bar(mode_list,vib_energy,width,bottom=vib_bottom)
        plots = ( p1[0], p2[0], p3[0] )
        legends = ('translation','rotation','vibration')
        #self.subplot.xticks(mode_list,mode_list_text)
        #self.subplot.yticks(np.arrange(0,101,10))
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        #self.subplot.legend( ploc='best')
        self.subplot.legend( plots, legends)
        self.subplot.set_title('Internal-External Composition of Vibrational Energy')
        self.canvas.draw_idle()

