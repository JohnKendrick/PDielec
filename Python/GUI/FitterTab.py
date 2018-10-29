# -*- coding: utf8 -*-
import sys
import os.path
import numpy as np
import Python.Calculator as Calculator
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit, QDoubleSpinBox
from PyQt5.QtWidgets    import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets    import  QSpinBox
from PyQt5.QtWidgets    import  QFileDialog
from PyQt5.QtWidgets    import  QTableWidget
from PyQt5.QtWidgets    import  QTableWidgetItem
from PyQt5.QtWidgets    import  QSizePolicy
from PyQt5.QtWidgets    import  QCheckBox
from PyQt5.QtCore       import  Qt
from PyQt5.QtCore       import  QCoreApplication
from Python.Utilities   import  Debug
from Python.GUI.SettingsTab import  FixedQTableWidget
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from Python.Utilities import Debug
import time
from openpyxl import load_workbook
from scipy.interpolate import interp1d
from scipy.optimize import minimize

class FitterTab(QWidget):
    def __init__(self, parent, debug=False):   
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'FitterTab:')
        self.dirty = True
        self.settings = {}
        self.notebook = parent
        self.notebook.plottingCalculationRequired = True
        self.settings['Excel file name'] = ''
        self.settings['Plot title'] = 'Experimental and Calculated Spectral Comparison'
        self.settings['Plot type'] = 'Molar absorption'
        self.plot_type_definitions = ['Molar absorption', 'Absorption', 'ATR','Real','Imaginary']
        self.settings['Number of iterations'] = 20
        self.settings['Frequency scaling'] = False
        self.settings['Frequency scaling factor'] = 1.0
        self.settings['Absorption scaling'] = False
        self.settings['Absorption scaling factor'] = 1.0
        self.plot_frequency_shift = False
        self.xcorr0=0.0
        self.xcorr1=0.0
        self.lag=0.0
        self.excel_file_has_been_read = False
        self.sigmas_cm1 = []
        self.modes_fitted = []
        self.modes_selected = []
        self.excel_frequencies = []
        self.excel_absorption = []
        self.experimental_absorption = []
        self.frequency_units = 'wavenumber'
        # Create a scenario tab 
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # Ask for the excel spread sheet
        #
        self.spectrafile_le = QLineEdit(self)
        self.spectrafile_le.setToolTip('Provide the name of an .xlsx file of spectrum')
        self.spectrafile_le.setText(self.settings['Excel file name'])
        self.spectrafile_le.returnPressed.connect(self.on_spectrafile_le_return)
        self.spectrafile_le.textChanged.connect(self.on_spectrafile_le_changed)
        label = QLabel('Excel spread sheet with spectra')
        label.setToolTip('Provide the name of an .xlsx containing the experimental spectrum.  The spreadsheet should have two columns, with no headings.  The first column contains the frequencies in cm-1.  The second contains the experimental spectrum.')
        form.addRow(label, self.spectrafile_le)
        # 
        # Select the type of plot we are going to use
        #
        self.plot_type_cb = QComboBox(self)
        self.plot_type_cb.setToolTip('What type of data are we going to fit to?')
        self.plot_type_cb.addItems(self.plot_type_definitions)
        self.plot_type_cb.activated.connect(self.on_plot_type_cb_activated)
        self.plot_type_cb.setCurrentIndex(self.plot_type_definitions.index(self.settings['Plot type']))
        self.settings['Plot type'] = self.plot_type_definitions[0]
        label = QLabel('Plot and data type')
        label.setToolTip('What type of data is stored in the experiment spread-sheet')
        form.addRow(label,self.plot_type_cb)
        # get initial sigmas from the settings tab
        self.sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        self.modes_selected = self.notebook.settingsTab.modes_selected
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        self.intensities = self.notebook.settingsTab.intensities
        if len(self.modes_fitted) == 0 and len(self.modes_selected) > 0:
            self.modes_selected = [ False  for _ in self.modes_selected ]
        self.sigmas_tw = QTableWidget(self)
        #self.sigmas_tw = FixedQTableWidget(self)
        self.sigmas_tw.setToolTip('Choose the sigmas which will be used in the fitting')
        self.sigmas_tw.itemChanged.connect(self.on_sigmas_tw_itemChanged)
        self.sigmas_tw.setRowCount(len(self.sigmas_cm1))
        self.sigmas_tw.setColumnCount(3)
        self.sigmas_tw.setHorizontalHeaderLabels(['   Sigma   \n(cm-1)', ' Frequency \n(cm-1)', '  Intensity  \n(Debye2/Å2/amu)'])
        self.redraw_sigmas_tw()
        label = QLabel('Lorentzian widths:')
        label.setToolTip('The Lorentzian widths can be edited here.  If checked they will also be used in the optimisation of the cross-correlation between the experiment and calculated spectra')
        form.addRow(label,self.sigmas_tw)
        #
        # See if we want frequency scaling
        #
        self.frequency_scaling_cb = QCheckBox(self)
        self.frequency_scaling_cb.setToolTip('Frequency scaling is applied during the fitting process')
        self.frequency_scaling_cb.setText('')
        self.frequency_scaling_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings['Frequency scaling']:
            self.frequency_scaling_cb.setCheckState(Qt.Checked)
        else:
            self.frequency_scaling_cb.setCheckState(Qt.Unchecked)
        self.frequency_scaling_cb.stateChanged.connect(self.on_frequency_scaling_cb_changed)
        self.frequency_scaling_factor_sb = QDoubleSpinBox(self)
        self.frequency_scaling_factor_sb.setRange(0.000001,10000000.0)
        self.frequency_scaling_factor_sb.setSingleStep(0.1)
        self.frequency_scaling_factor_sb.setValue(self.settings['Frequency scaling factor'])
        self.frequency_scaling_factor_sb.setToolTip('Set the value for scaling the frequency axis of the calculated spectrum')
        self.frequency_scaling_factor_sb.valueChanged.connect(self.on_frequency_scaling_factor_sb_changed)
        hbox = QHBoxLayout()
        hbox.addWidget(self.frequency_scaling_cb)
        hbox.addWidget(self.frequency_scaling_factor_sb)
        label = QLabel('Set frequency scaling factor to be optimised')
        label.setToolTip('Set frequency scaling factor to be optimised.  The frequency scaling factor given will always be used, even if not included in the optimisation.  The plot does not reflect the frequency scaling.')
        form.addRow(label, hbox)
        #
        # See if we want absorption scaling
        #
#        self.absorption_scaling_cb = QCheckBox(self)
#        self.absorption_scaling_cb.setToolTip('Absorption scaling is applied during the fitting process')
#        self.absorption_scaling_cb.setText('')
#        self.absorption_scaling_cb.setLayoutDirection(Qt.RightToLeft)
#        if self.settings['Absorption scaling']:
#            self.absorption_scaling_cb.setCheckState(Qt.Checked)
#        else:
#            self.absorption_scaling_cb.setCheckState(Qt.Unchecked)
#        self.absorption_scaling_cb.stateChanged.connect(self.on_absorption_scaling_cb_changed)
#        self.absorption_scaling_factor_sb = QDoubleSpinBox(self)
#        self.absorption_scaling_factor_sb.setRange(0.000001,10000000.0)
#        self.absorption_scaling_factor_sb.setSingleStep(0.1)
#        self.absorption_scaling_factor_sb.setValue(self.settings['Absorption scaling factor'])
#        self.absorption_scaling_factor_sb.setToolTip('Set the value for scaling the absorption axis of the calculated spectrum')
#        self.absorption_scaling_factor_sb.valueChanged.connect(self.on_absorption_scaling_factor_sb_changed)
#        hbox = QHBoxLayout()
#        hbox.addWidget(self.absorption_scaling_cb)
#        hbox.addWidget(self.absorption_scaling_factor_sb)
#        label = QLabel('Set absorption scaling factor to be optimised')
#        label.setToolTip('Set absorption scaling factor to be optimised.  The absorption scaling factor given will always be used, even if not included in the optimisation.  The plot does not reflect the absorption scaling.')
#        form.addRow(label, hbox)
        # Add the number of iterations
        self.iterations_sb = QSpinBox(self)
        self.iterations_sb.setRange(1,900)
        self.iterations_sb.setValue(self.settings['Number of iterations'])
        self.iterations_sb.setToolTip('Set the number of iterations to be used to optimise the cross-correlation coefficient')
        self.iterations_sb.valueChanged.connect(self.on_iterations_sb_changed)
        label = QLabel('Number of iterations:')
        label.setToolTip('Set the number of iterations to be used to optimise the cross-correlation coefficient')
        form.addRow(label,self.iterations_sb)
        # Add output of the cross correlation coefficient
        hbox = QHBoxLayout()
        self.cross_correlation_le = QLineEdit(self)
        self.cross_correlation_le.setEnabled(False)
        self.cross_correlation_le.setText('{}'.format(0.0))
        self.lag_frequency_le = QLineEdit(self)
        self.lag_frequency_le.setEnabled(False)
        self.lag_frequency_le.setText('{}'.format(0.0))
        hbox.addWidget(self.cross_correlation_le)
        hbox.addWidget(self.lag_frequency_le)
        label = QLabel('Cross correlation value and shift')
        label.setToolTip('The highest cross-correlation value and its associated frequency shift is shown')
        form.addRow(label, hbox)
        # Add a replot and recalculate button
        hbox = QHBoxLayout()
        self.replotButton1 = QPushButton('Replot')
        self.replotButton1.setToolTip('Recalculate the absorption with the new sigma values')
        self.replotButton1.clicked.connect(self.replotButton1Clicked)
        hbox.addWidget(self.replotButton1)
        self.replotButton2 = QPushButton('Replot with frequency shift')
        self.replotButton2.setToolTip('Recalculate the absorption with the new sigma values, including a shft in the frequencies to maximise the cross-correlation')
        self.replotButton2.clicked.connect(self.replotButton2Clicked)
        hbox.addWidget(self.replotButton2)
        # Add a fitting button
        self.fittingButton = QPushButton('Perform fitting')
        self.fittingButton.setToolTip('Attempt to fit the calculated spectrum to the experimental one')
        self.fittingButton.clicked.connect(self.fittingButtonClicked)
        hbox.addWidget(self.fittingButton)
        form.addRow(hbox)
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
        self.dirty = True

    def on_iterations_sb_changed(self):
         self.settings['Number of iterations'] = self.iterations_sb.value()

    def on_absorption_scaling_cb_changed(self,value):
        #print('on_absorption_scaling_cb_changed',value)
        self.settings['Absorption scaling'] = self.absorption_scaling_cb.isChecked()
        return
 
    def on_absorption_scaling_factor_sb_changed(self,value):
        #print('on_absorption_scaling_factor_sb_changed',value)
        self.settings['Absorption scaling factor'] = float(value)
        return

    def on_frequency_scaling_cb_changed(self,value):
        #print('on_frequency_scaling_cb_changed',value)
        self.settings['Frequency scaling'] = self.frequency_scaling_cb.isChecked()
        return
 
    def on_frequency_scaling_factor_sb_changed(self,value):
        #print('on_frequency_scaling_factor_sb_changed',value)
        self.settings['Frequency scaling factor'] = float(value)
        return

    def replotButton1Clicked(self):
        self.dirty = True
        self.notebook.plottingTab.refresh(force=True)
        self.plot_frequency_shift = False
        self.refresh(force=True)
        return

    def replotButton2Clicked(self):
        self.dirty = True
        self.notebook.plottingTab.refresh(force=True)
        self.plot_frequency_shift = True
        self.refresh(force=True)
        return

    def plot(self,experiment,xs,ys,legends,label):
        # Plot the experimental values on the left y-axis
        # Plot all the others in xs, ys on the right x-axis
        self.subplot1 = None
        self.figure.clf()
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap("tab10")
        if self.frequency_units == 'wavenumber':
            xlabel = r'Frequency $\mathdefault{(cm^{-1})}}$'
            scale = 1.0
        else:
            xlabel = r'THz'
            scale = 0.02998
        self.subplot1 = self.figure.add_subplot(111)
        self.subplot2 = self.subplot1.twinx()
        cmap_index = 0
        lines = []
        if self.plot_frequency_shift:
            lag = float(self.lag_frequency_le.text())
            scale_calc = scale * self.settings['Frequency scaling factor']
        else:
            lag = 0.0
            scale_calc = scale
        for x,y,legend in zip(xs,ys,legends):
            x = np.array(x)
            line, = self.subplot1.plot(lag+scale_calc*x,y,lw=2, color=cmap(cmap_index), label=legend )
            lines.append(line)
            cmap_index += 1
        if len(experiment) > 0:
            # Use the x variables from the previous xs, ys
            line, = self.subplot2.plot(scale*x,experiment,lw=2, color=cmap(cmap_index), label='Experiment' )
            lines.append(line)
        labels = [l.get_label() for l in lines]
        self.subplot2.set_ylabel('Experiment')
        self.subplot1.set_ylabel('Calculated '+self.settings['Plot type'] )
        self.subplot1.legend(lines, labels, loc='best')
        self.subplot1.set_title(self.settings['Plot title'])
        self.canvas.draw_idle()

    def fittingButtonClicked(self):
        self.refresh()
        #print('Cross correlations', lag,xcorr0,xcorr1)
        final_point = self.optimiseFit()
        return

    def optimiseFit(self):
        # Optimise the fit of the first scenario to the experimental data
        # First determine who many variables we have
        self.fit_list = []
        for mode,fitted in enumerate(self.modes_fitted):
            if fitted:
                self.fit_list.append(mode)
        initial_point = [ self.sigmas_cm1[i] for i in self.fit_list ]
        # Append a scaling option
        if self.settings['Frequency scaling']:
            initial_point.append(self.settings['Frequency scaling factor'])
        nvariables = len(initial_point)
        if nvariables > 0:
            final_point = minimize(self.optimiseFunction, initial_point, method='nelder-mead', options={'xtol':1.0, 'disp':False, 'maxiter':self.settings['Number of iterations']} )
        else: 
             print('No sigmas have been selected for optimisation')
             final_point = []
        return final_point

    def optimiseFunction(self,variables) :
        # Determine the function to be optimised (minimised)
        #print('Current point in optimisation')
        #print(variables)
        if self.settings['Frequency scaling']:
            sigmas = variables[:-1]
            scaling_factor = variables[-1]
        else:
            sigmas = variables
            scaling_factor = self.settings['Frequency scaling factor']
        #print('sigmas',sigmas)
        #print('scaling_factor',scaling_factor)
        for index,sigma in zip(self.fit_list,sigmas):
            self.sigmas_cm1[index] = sigma
            self.redraw_sigmas_tw()
            self.notebook.settingsTab.sigmas_cm1[index] = sigma
        self.notebook.settingsTab.redraw_output_tw()
        self.notebook.plottingTab.refresh(force=True)
        self.refresh(force=True)
        #print(' ')
        # Returning the best correlation but made negative because we need to minimise
        return -1.0*self.xcorr0

    def calculateCrossCorrelation(self,scaling_factor):
        # Calculate the cross correlation coefficient between the experimental and the first scenario
        if len(self.experimental_absorption) == 0:
            return (0.0,0.0,0.0)
        col1 = np.array(self.experimental_absorption)
        col1 = ( col1 - np.mean(col1)) / ( np.std(col1) * np.sqrt(len(col1)) )
        col2 = np.array(self.calculated_absorptions[0])
        f = interp1d(scaling_factor*np.array(self.xaxis), col2, kind='cubic',fill_value='extrapolate')
        col2 = f(self.xaxis)
        col2 = ( col2 - np.mean(col2)) / ( np.std(col2) * np.sqrt(len(col2)) )
        correlation = np.correlate(col1, col2,  mode='full')
        lag = np.argmax(correlation) - (len(correlation)-1)/2
        lag = (self.xaxis[1] - self.xaxis[0]) * lag
        return (lag,np.max(correlation),correlation[int((len(correlation)-1)/2)])


    def redraw_sigmas_tw(self):
        if len(self.sigmas_cm1) <= 0:
            return
        self.sigmas_tw.blockSignals(True)
        self.sigmas_tw.setRowCount(len(self.sigmas_cm1))
        for i,(f,sigma,intensity) in enumerate(zip(self.frequencies_cm1, self.sigmas_cm1, self.intensities)):
            # Sigma and check / unchecked column
            items = []
            itemFlags = []
            item = QTableWidgetItem('{0:.2f}'.format(sigma))
            if self.modes_selected[i]:
                if self.modes_fitted[i]:
                    item.setCheckState(Qt.Checked)
                else:
                    item.setCheckState(Qt.Unchecked)
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable )
                otherFlags = item.flags() & Qt.NoItemFlags | Qt.ItemIsEnabled
            else:
                #itemFlags.append( item.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled & ~Qt.ItemIsSelectable & ~Qt.ItemIsEditable )
                itemFlags.append( item.flags() & Qt.NoItemFlags | Qt.ItemIsUserCheckable | Qt.ItemIsEnabled )
                item.setCheckState(Qt.Unchecked)
                otherFlags = item.flags() & Qt.NoItemFlags
            items.append(item)
            # Frequency column cm-1
            items.append(QTableWidgetItem('{0:.4f}'.format(f) ) )
            itemFlags.append( otherFlags )
            # Intensity column Debye2/Angs2/amu
            items.append(QTableWidgetItem('{0:.4f}'.format(intensity) ) )
            itemFlags.append( otherFlags )
            for j,(item,flag) in enumerate(zip(items,itemFlags)):
                item.setFlags(flag)
                item.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
                self.sigmas_tw.setItem(i, j, item )
        self.sigmas_tw.resizeColumnsToContents()
        # Release the block on signals for the frequency output table
        self.sigmas_tw.blockSignals(False)
        QCoreApplication.processEvents()

    def on_plot_type_cb_activated(self,index):
        self.dirty = True
        self.settings['Plot type'] = self.plot_type_definitions[index]

    def on_spectrafile_le_return(self):
        # Handle a return in the excel file name line editor
        file_name = self.spectrafile_le.text()
        if not os.path.isfile(file_name):
            qfd = QFileDialog(self)
            qfd.setDirectory(self.notebook.mainTab.directory)
            file_name, _ = qfd.getOpenFileName(self,'Open Excel file','','Excel (*.xlsx);; All Files (*)')
        if not os.path.isfile(file_name):
            return
        self.settings['Excel file name'] = file_name
        self.spectrafile_le.setText(self.settings['Excel file name'])
        debugger.print('new file name', self.settings['Excel file name'])
        self.dirty = True
        self.excel_file_has_been_read = False
        return

    def read_excel_file(self):
        # 
        file_name = self.settings['Excel file name']
        if not os.path.isfile(file_name):
            return
        if self.excel_file_has_been_read:
            return
        # 
        # Open the work book
        #
        self.wb = load_workbook(filename=self.settings['Excel file name'], read_only=True)
        self.ws = self.wb.worksheets[0]
        self.excel_frequencies = []
        self.excel_absorption = []
        for row in self.ws.rows:
            self.excel_frequencies.append(row[0].value)
            self.excel_absorption.append(row[1].value)
        self.excel_file_has_been_read = True
        return

    def on_spectrafile_le_changed(self,text):
        debugger.print('on_spectrafile_le_changed', text)
        text = self.spectrafile_le.text()
        self.dirty = True
        self.settings['Excel file name'] = text
        self.excel_file_has_been_read = False

    def on_sigmas_tw_itemChanged(self, item):
        self.sigmas_tw.blockSignals(True)
        debugger.print('on_sigmas_tw_itemChanged)', item.row(), item.column() )
        col = item.column()
        row = item.row()
        if col == 0:
            # If this is the first column alter the check status but reset the sigma value
            if item.checkState() == Qt.Checked:
                self.modes_fitted[row] = True
            else:
                 self.modes_fitted[row] = False
            new_value = float(item.text())
            self.sigmas_cm1[row] = new_value
            self.redraw_sigmas_tw()
            self.notebook.settingsTab.sigmas_cm1[row] = new_value
            self.notebook.settingsTab.redraw_output_tw()
        elif col == 1:
            self.redraw_sigmas_tw()
        else:
            self.redraw_sigmas_tw()
        self.notebook.plottingCalculationRequired = True
        self.notebook.analysisCalculationRequired = True
        self.notebook.fittingCalculationRequired = True
        QCoreApplication.processEvents()


    def print_settings(self):
        debugger.print('SETTINGS')
        for key in self.settings:
            debugger.print(key, self.settings[key]) 
        
    def refresh(self,force=False):
        if not self.dirty and not force and not self.notebook.fittingCalculationRequired:
            debugger.print('refresh aborted', self.dirty,force)
            return
        debugger.print('refresh', force)
        #
        # Block signals during refresh
        # 
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # 
        # use the settings values to initialise the widgets
        #
        self.spectrafile_le.setText(self.settings['Excel file name'])
        self.read_excel_file()
        self.plot_type_cb.setCurrentIndex(self.plot_type_definitions.index(self.settings['Plot type']))
        self.iterations_sb.setValue(self.settings['Number of iterations'])
        # 
        # If the sigmas are not set return
        #
        self.sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        if len(self.sigmas_cm1) < 1:
            return
        self.modes_selected = self.notebook.settingsTab.modes_selected
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        self.intensities = self.notebook.settingsTab.intensities
        if len(self.modes_fitted) == 0 and len(self.modes_selected) > 0:
            self.modes_fitted = [ False  for _ in self.modes_selected ]
        self.redraw_sigmas_tw()
        # Resample the spectrum
        self.calculated_frequencies = self.notebook.plottingTab.xaxes
        if len(self.excel_absorption) > 0:
            self.resample_spectrum()
        index = self.plot_type_cb.currentIndex()
        self.settings['Plot type'] = self.plot_type_definitions[index]
        if index == 0:
            self.calculated_absorptions = self.notebook.plottingTab.molarAbsorptionCoefficients
            label = 'Molar absorption coefficient'
        elif index == 1:
            self.calculated_absorptions = self.notebook.plottingTab.absorptionCoefficients
            label = 'Absorption coefficient'
        elif index == 2:
            self.calculated_absorptions = self.notebook.plottingTab.sp_atrs
            label = 'Absorption coefficient'
        elif index == 3:
            self.calculated_absorptions = self.notebook.plottingTab.realPermittivities
            label = 'Real Permittivity'
        elif index == 4:
            self.calculated_absorptions = self.notebook.plottingTab.imagPermittivities
            label = 'Imaginary Permittivity'
        legends = self.notebook.plottingTab.legends
        scaling_factor = self.settings['Frequency scaling factor']
        self.lag,self.xcorr0,self.xcorr1 = self.calculateCrossCorrelation(scaling_factor)
        self.cross_correlation_le.setText('{}'.format(self.xcorr0))
        self.lag_frequency_le.setText('{}'.format(self.lag))
        self.plot(self.experimental_absorption,self.calculated_frequencies,self.calculated_absorptions,legends,label)
        #
        # Unblock signals after refresh
        # 
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.dirty = False
        self.notebook.fittingCalculationRequired = False
        return


    def resample_spectrum(self):
    #
    #  The experimental spectrum needs to be in the same range as the calculated spectrum
    #  It also needs to be calculated at the same frequencies as the calculated spectrum
    #
        self.xaxis = self.calculated_frequencies[0]
        f = interp1d(self.excel_frequencies, self.excel_absorption, kind='cubic',fill_value='extrapolate')
        self.experimental_absorption = f(self.xaxis)
        return
    
