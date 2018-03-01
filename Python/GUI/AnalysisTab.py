import sys
import os.path
import os
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QFileDialog, QProgressBar
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox
from PyQt5.QtWidgets  import  QSizePolicy
from PyQt5.QtCore     import  Qt
from Python.Constants import  wavenumber, amu, PI, avogadro_si, angstrom
from Python.Constants import  average_masses, isotope_masses
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from Python.Utilities import Debug

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
        self.vmin_sb.setRange(0,9000)
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
        # Set the x-axis frequency units
        #
        self.funits_cb = QComboBox(self) 
        self.funits_cb.setToolTip('Set the frequency units for the x-axis')
        for choice in ['wavenumber','THz']:
            self.funits_cb.addItem(choice)
        self.frequency_units = 'wavenumber'
        self.funits_cb.currentIndexChanged.connect(self.on_funits_cb_changed)
        label = QLabel('Frequency units for the x-axis', self)
        label.setToolTip('Set the frequency units for the x-axis')
        form.addRow(label, self.funits_cb)
        # Add the matplotlib figure to the bottom 
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
        self.notebook.newCalculationRequired = True
        self.progressbar.setValue(0)
        debugger.print('on vmin change ', self.settings['vmin'])

    def on_vmax_changed(self):
        self.settings['vmax'] = self.vmax_sb.value()
        self.notebook.newCalculationRequired = True
        self.progressbar.setValue(0)
        debugger.print('on vmax change ', self.settings['vmax'])

    def refresh(self):
        debugger.print('refreshing widget')
        self.calculate()
        return

    def on_funits_cb_changed(self, index):
        if index == 0:
            self.frequency_units = 'wavenumber'
        else:
            self.frequency_units = 'THz'
        self.replot()
        debugger.print('Frequency units changed to ', self.frequency_units)

    def calculate(self):
        debugger.print('calculate')
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['program']
        filename = settings['filename']
        reader = self.notebook.mainTab.reader
        if reader is None:
            return
        if program is '':
            return
        if filename is '':
            return
        # Assemble the settingsTab settings
        self.notebook.newCalculationRequired = False
        settings = self.notebook.settingsTab.settings
        eckart = settings['eckart']
        neutral = settings['neutral']
        hessian_symm = settings['hessian_symmetrisation']
        epsilon_inf = np.array(settings['optical_permittivity'])
        sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        sigmas = np.array(sigmas_cm1) * wavenumber
        modes_selected = self.notebook.settingsTab.modes_selected
        frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        frequencies = np.array(frequencies_cm1) * wavenumber
        intensities = self.notebook.settingsTab.intensities
        oscillator_strengths = self.notebook.settingsTab.oscillator_strengths
        volume = reader.volume*angstrom*angstrom*angstrom
        vmin = self.settings['vmin']
        vmax = self.settings['vmax']
        cell = reader.unit_cells[-1]
        atom_masses = self.reader.masses
        cell.set_atomic_masses(atom_masses)
        newcell,nmols,old_order = cell.calculate_molecular_contents(scale, toler, covalent_radii)
        newcell.printInfo()
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
        mode_energies = Calculator.calculate_energy_distribution(newcell, frequencies_cm1, new_mass_weighted_normal_modes)
        # Output the final result
        title = ['Freq(cm-1)','%mol-cme','%mol-rot','%vib']
        for i in range(nmols):
            title.append('%mol-'+str(i))
        print_strings('Percentage energies in vibrational modes',title,format="{:>10}")
        fd_csvfile = None
        if not fd_csvfile is None:
            print_strings('Percentage energies in vibrational modes',title,format="{:>10}",file=fd_csvfile,separator=",")
        for freq,energies in zip(frequencies_cm1,mode_energies):
               tote,cme,rote,vibe,molecular_energies = energies
               output = [ freq, 100*cme/tote, 100*rote/tote, 100*vibe/tote ]
               for e in molecular_energies:
                   output.append(100*e/tote)
               print_reals('',output,format="{:10.2f}")
               if not fd_csvfile is None:
                   print_reals('',output,format="{:10.2f}",file=fd_csvfile,separator=",")
        # if using a excel file write out the results
        fd_excelfile = None
        if not fd_excelfile is None:
            excel_row += 2; worksheet.write(excel_row,0,'Percentage energies in vibrational modes')
            excel_row += 1; [ worksheet.write(excel_row,col,f) for col,f in enumerate(title) ]
            for freq,energies in zip(frequencies_cm1,mode_energies):
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
        self.xaxes = []

    def plot(self,xs,ys,ylabel):
        # import matplotlib.pyplot as pl
        # mp.use('Qt5Agg')
        self.subplot = None
        self.remember_xs = xs
        self.remember_ys = ys
        self.remember_ylabel = ylabel
        self.figure.clf()
        if self.frequency_units == 'wavenumber':
            xlabel = r'Frequency $\mathdefault{(cm^{-1})}}$'
            scale = 1.0
        else:
            xlabel = r'THz'
            scale = 0.02998
        self.subplot = self.figure.add_subplot(111)
        for scenario,x,y in zip(self.scenarios,xs,ys):
            x = np.array(x)
            legend = scenario.settings['legend']
            line, = self.subplot.plot(scale*x,y,lw=2, label=legend )
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(ylabel)
        self.subplot.legend(loc='best')
        self.subplot.set_title(self.settings['title'])
        self.canvas.draw_idle()

    def replot(self):
        if self.subplot is not None:
            self.plot(self.remember_xs, self.remember_ys, self.remember_ylabel)
            
 
