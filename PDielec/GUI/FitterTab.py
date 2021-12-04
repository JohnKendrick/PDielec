# -*- coding: utf8 -*-
import os.path
import numpy as np
import PDielec.Calculator as Calculator
from PyQt5.QtWidgets    import  QPushButton, QWidget, QProgressBar
from PyQt5.QtWidgets    import  QComboBox, QLabel, QLineEdit, QDoubleSpinBox
from PyQt5.QtWidgets    import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets    import  QSpinBox
from PyQt5.QtWidgets    import  QFileDialog
from PyQt5.QtWidgets    import  QTableWidgetItem
from PyQt5.QtWidgets    import  QSizePolicy
from PyQt5.QtWidgets    import  QCheckBox
from PyQt5.QtWidgets    import  QTabWidget
from PyQt5.QtCore       import  Qt
from PyQt5.QtCore       import  QCoreApplication
from PDielec.Utilities   import  Debug
from PDielec.GUI.SettingsTab import  FixedQTableWidget
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from scipy.interpolate import interp1d
from scipy.optimize import minimize

class FitterTab(QWidget):
    def __init__(self, parent, debug=False):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'FitterTab:')
        debugger.print('Start:: Initialising')
        self.refreshRequired = True
        self.calculationInProgress = False
        self.settings = {}
        self.notebook = parent
        self.settings['Excel file name'] = ''
        self.settings['Plot title'] = 'Experimental and Calculated Spectral Comparison'
        self.settings['Fitting type'] = 'Minimise x-correlation'
        self.fitting_type_definitions = ['Minimise x-correlation', 'Minimise spectral difference']
        self.settings['Number of iterations'] = 20
        self.settings['Frequency scaling factor'] = 1.0
        self.settings['Optimise frequency scaling'] = False
        self.settings['Spectrum scaling'] = False
        self.settings['Spectrum scaling factor'] = 1.0
        self.settings['Independent y-axes'] = True
        self.settings['Spectral difference threshold'] = 0.05
        self.settings['HPFilter lambda'] = 7.0
        self.settings['Baseline removal'] = False
        self.settings['Scenario index'] = len(self.notebook.scenarios) - 1
        self.scenario_legends = [ scenario.settings['Legend'] for scenario in self.notebook.scenarios ]
        self.lastButtonPressed = self.replotButton1Clicked
        self.plot_frequency_shift = False
        self.xcorr0=0.0
        self.xcorr1=0.0
        self.lag=0.0
        self.excel_file_has_been_read = False
        self.sigmas_cm1 = []
        self.modes_fitted = []
        self.modes_selected = []
        self.excel_frequencies = []
        self.excel_spectrum = []
        self.experimental_spectrum = []
        self.frequency_units = 'wavenumber'
        # Create a tab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # Ask for the excel spread sheet
        #
        self.spectrafile_le = QLineEdit(self)
        self.spectrafile_le.setToolTip('Provide the name of an .xlsx containing the experimental spectrum.  The spreadsheet should have two columns, with no headings.  The first column contains the frequencies in cm-1.  The second contains the experimental spectrum.')
        self.spectrafile_le.setText(self.settings['Excel file name'])
        self.spectrafile_le.returnPressed.connect(self.on_spectrafile_le_return)
        self.spectrafile_le.textChanged.connect(self.on_spectrafile_le_changed)
        #
        # Select the type of fitting we are going to use
        #
        self.fitting_type_cb = QComboBox(self)
        self.fitting_type_cb.setToolTip('What type of fitting?')
        self.fitting_type_cb.addItems(self.fitting_type_definitions)
        self.fitting_type_cb.activated.connect(self.on_fitting_type_cb_activated)
        self.fitting_type_cb.setCurrentIndex(self.fitting_type_definitions.index(self.settings['Fitting type']))
        #
        # Select the type of plot we are going to use
        #
        self.scenario_cb = QComboBox(self)
        self.scenario_cb.setToolTip('What type of data is stored in the experimental spreadsheet?')
        self.scenario_cb.addItems(self.scenario_legends)
        self.scenario_cb.activated.connect(self.on_scenario_cb_activated)
        self.scenario_cb.setCurrentIndex(self.settings['Scenario index'])
        #
        # See if we want frequency scaling
        #
        self.frequency_scaling_factor_sb = QDoubleSpinBox(self)
        self.frequency_scaling_factor_sb.setToolTip('Frequency scaling factor')
        self.frequency_scaling_factor_sb.setRange(0.000001,10000000.0)
        self.frequency_scaling_factor_sb.setDecimals(4)
        self.frequency_scaling_factor_sb.setSingleStep(0.1)
        self.frequency_scaling_factor_sb.setValue(self.settings['Frequency scaling factor'])
        self.frequency_scaling_factor_sb.setToolTip('Set the value for scaling the frequency axis of the calculated spectrum')
        self.frequency_scaling_factor_sb.valueChanged.connect(self.on_frequency_scaling_factor_sb_changed)
        #
        # See if we want base line removal
        #
        self.baseline_cb = QCheckBox(self)
        self.baseline_cb.setToolTip('Apply base line correction')
        self.baseline_cb.setText('')
        self.baseline_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings['Baseline removal']:
            self.baseline_cb.setCheckState(Qt.Checked)
        else:
            self.baseline_cb.setCheckState(Qt.Unchecked)
        self.baseline_cb.stateChanged.connect(self.on_baseline_cb_changed)

        # Hodrick-Prescott Filter Lambda
        self.hpfilter_lambda_sb = QDoubleSpinBox(self)
        self.hpfilter_lambda_sb.setRange(0.0,10000000.0)
        self.hpfilter_lambda_sb.setDecimals(1)
        self.hpfilter_lambda_sb.setSingleStep(1)
        self.hpfilter_lambda_sb.setValue(self.settings['HPFilter lambda'])
        self.hpfilter_lambda_sb.setToolTip('The Hodrick-Prescott filter lambda value')
        self.hpfilter_lambda_sb.valueChanged.connect(self.on_hpfilter_lambda_sb_changed)
        # Spectral difference threshold
        self.spectral_difference_threshold_sb = QDoubleSpinBox(self)
        self.spectral_difference_threshold_sb.setRange(0.000001,10.0)
        self.spectral_difference_threshold_sb.setDecimals(2)
        self.spectral_difference_threshold_sb.setSingleStep(0.1)
        self.spectral_difference_threshold_sb.setValue(self.settings['Spectral difference threshold'])
        self.spectral_difference_threshold_sb.setToolTip('Set the value the spectral difference threshold')
        self.spectral_difference_threshold_sb.valueChanged.connect(self.on_spectral_difference_threshold_sb_changed)
        # Add the number of iterations
        self.iterations_sb = QSpinBox(self)
        self.iterations_sb.setRange(1,900)
        self.iterations_sb.setValue(self.settings['Number of iterations'])
        self.iterations_sb.setToolTip('Set the number of iterations to be used to optimise the cross-correlation coefficient')
        self.iterations_sb.valueChanged.connect(self.on_iterations_sb_changed)
        # Independent y-axes
        self.independent_yaxes_cb = QCheckBox(self)
        self.independent_yaxes_cb.setToolTip('Check to use indpendent y-axes in the plot')
        self.independent_yaxes_cb.setText('')
        self.independent_yaxes_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings['Independent y-axes']:
            self.independent_yaxes_cb.setCheckState(Qt.Checked)
        else:
            self.independent_yaxes_cb.setCheckState(Qt.Unchecked)
        self.independent_yaxes_cb.stateChanged.connect(self.on_independent_yaxes_cb_changed)
        # Optimise frequency scaling?
        self.optimise_frequency_scaling_cb = QCheckBox(self)
        self.optimise_frequency_scaling_cb.setToolTip('Use frequency scaling in optimisation')
        self.optimise_frequency_scaling_cb.setText('')
        self.optimise_frequency_scaling_cb.setLayoutDirection(Qt.RightToLeft)
        if self.settings['Optimise frequency scaling']:
            self.optimise_frequency_scaling_cb.setCheckState(Qt.Checked)
        else:
            self.optimise_frequency_scaling_cb.setCheckState(Qt.Unchecked)
        self.optimise_frequency_scaling_cb.stateChanged.connect(self.on_optimise_frequency_scaling_cb_changed)
        #
        # Add a tab widget for the settings ######################################################################################
        #
        self.settingsTab = QTabWidget(self)
        self.settingsTab.addTab(self.spectrafile_le, 'Experimental spectrum')
        self.settingsTab.addTab(self.scenario_cb, 'Scenario')
        self.settingsTab.addTab(self.frequency_scaling_factor_sb, 'Frequency scaling factor')
        self.settingsTab.addTab(self.iterations_sb, 'No. of iterations')
        self.settingsTab.addTab(self.independent_yaxes_cb, 'Independent y-axes')
        self.settingsTab.addTab(self.fitting_type_cb, 'Fitting type')
        self.settingsTab.addTab(self.optimise_frequency_scaling_cb, 'Optimise scaling')
        self.settingsTab.addTab(self.spectral_difference_threshold_sb, 'Spectral difference threshold')
        self.settingsTab.addTab(self.baseline_cb, 'Baseline removal?')
        self.settingsTab.addTab(self.hpfilter_lambda_sb, 'HP Filter Lambda')
        label = QLabel('Options:', self)
        form.addRow(label,self.settingsTab)
        # END OF THE SETTINGS TAB #################################################################################################
        # Add Lorentzian widths table
        # get initial sigmas from the settings tab
        self.sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        self.modes_selected = self.notebook.settingsTab.modes_selected
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        self.intensities = self.notebook.settingsTab.intensities
        if len(self.modes_fitted) == 0 and len(self.modes_selected) > 0:
            self.modes_selected = [ False  for _ in self.modes_selected ]
        self.sigmas_tw = FixedQTableWidget(parent=self,rows=7)
        #jk sizePolicy = QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        #jk self.sigmas_tw.setSizePolicy(sizePolicy)
        self.sigmas_tw.setToolTip('Choose the sigmas which will be used in the fitting')
        self.sigmas_tw.itemChanged.connect(self.on_sigmas_tw_itemChanged)
        self.sigmas_tw.setRowCount(max(8,len(self.sigmas_cm1)))
        self.sigmas_tw.setColumnCount(3)
        self.sigmas_tw.setHorizontalHeaderLabels(['   Sigma   \n(cm-1)', ' Frequency \n(cm-1)', '  Intensity  \n(Debye2/Å2/amu)'])
        self.redraw_sigmas_tw()
        label = QLabel('Lorentzian widths:')
        label.setToolTip('The Lorentzian widths can be edited here.  If checked they will also be used in the optimisation of the cross-correlation between the experiment and calculated spectra')
        form.addRow(label,self.sigmas_tw)
        # Add a replot and recalculate button
        hbox = QHBoxLayout()
        self.replotButton1 = QPushButton('Replot')
        self.replotButton1.setToolTip('Recalculate the spectrum with the new sigma values')
        self.replotButton1.clicked.connect(self.replotButton1Clicked)
        hbox.addWidget(self.replotButton1)
        self.replotButton2 = QPushButton('Replot with frequency shift')
        self.replotButton2.setToolTip('Recalculate the spectrum with the new sigma values, including a shft in the frequencies to maximise the cross-correlation')
        self.replotButton2.clicked.connect(self.replotButton2Clicked)
        hbox.addWidget(self.replotButton2)
        # Add a fitting button
        self.fittingButton = QPushButton('Perform fitting')
        self.fittingButton.setToolTip('Attempt to fit the calculated spectrum to the experimental one')
        self.fittingButton.clicked.connect(self.fittingButtonClicked)
        hbox.addWidget(self.fittingButton)
        form.addRow(hbox)
        # Add a progress bar
        label = QLabel('Calculation progress', self)
        label.setToolTip('Show the progress of any calculations')
        self.progressbar = QProgressBar(self)
        form.addRow(label,self.progressbar)
        self.notebook.progressbars_add(self.progressbar)
        # Add output of the cross correlation coefficient
        hbox = QHBoxLayout()
        self.cross_correlation_le = QLineEdit(self)
        self.cross_correlation_le.setEnabled(False)
        self.cross_correlation_le.setText('{}'.format(0.0))
        self.lag_frequency_le = QLineEdit(self)
        self.lag_frequency_le.setEnabled(False)
        self.lag_frequency_le.setText('{}'.format(0.0))
        self.frequency_scaling_le = QLineEdit(self)
        self.frequency_scaling_le.setEnabled(False)
        self.frequency_scaling_le.setText('{}'.format(self.settings['Frequency scaling factor']))
        self.rmse_le = QLineEdit(self)
        self.rmse_le.setEnabled(False)
        self.rmse_le.setText('{}'.format(0.0))
        hbox.addWidget(self.cross_correlation_le)
        hbox.addWidget(self.lag_frequency_le)
        hbox.addWidget(self.frequency_scaling_le)
        hbox.addWidget(self.rmse_le)
        label = QLabel('X-correlation, shift/lag, frequency scale and rmse')
        label.setToolTip('The highest cross-correlation value and its associated frequency shift is shown followed by the spectral error')
        form.addRow(label, hbox)
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
        self.refreshRequired = True
        debugger.print('Finished:: Initialising')

    def on_iterations_sb_changed(self):
        debugger.print('on_iterations_sb_changed')
        self.settings['Number of iterations'] = self.iterations_sb.value()
        self.refreshRequired = True
        return

    def on_independent_yaxes_cb_changed(self,value):
        debugger.print('independent_yaxes_cb_changed',value)
        self.settings['Independent y-axes'] = self.independent_yaxes_cb.isChecked()
        self.refreshRequired = True
        return

    def on_optimise_frequency_scaling_cb_changed(self,value):
        debugger.print('optimise_frequency_scaling_cb_changed',value)
        self.settings['Optimise frequency scaling'] = self.optimise_frequency_scaling_cb.isChecked()
        self.refreshRequired = True
        return

    def on_spectrum_scaling_cb_changed(self,value):
        debugger.print('on_spectrum_scaling_cb_changed',value)
        self.settings['Spectrum scaling'] = self.spectrum_scaling_cb.isChecked()
        self.refreshRequired = True
        return

    def on_hpfilter_lambda_sb_changed(self,value):
        debugger.print('hpfilter_lambda_sb_changed',value)
        self.refreshRequired = True
        try:
            self.settings['HPFilter lambda'] = float(value)
        except:
            print('Failed to convert to float', value)
        return

    def on_spectrum_scaling_factor_sb_changed(self,value):
        debugger.print('on_spectrum_scaling_factor_cb_changed',value)
        self.refreshRequired = True
        try:
            self.settings['Spectrum scaling factor'] = float(value)
        except:
            print('Failed to convert to float', value)
        return

    def on_baseline_cb_changed(self,value):
        debugger.print('on_baseline_cb_changed',value)
        self.refreshRequired = True
        self.settings['Baseline removal'] = self.baseline_cb.isChecked()
        return

    def on_spectral_difference_threshold_sb_changed(self,value):
        debugger.print('on_spectral_difference_threshold_sb_changed',value)
        self.refreshRequired = True
        try:
            self.settings['Spectral difference threshold'] = float(value)
        except:
            print('Failed to convert to float', value)
        return

    def on_frequency_scaling_factor_sb_changed(self,value):
        debugger.print('on_frequency_scaling_factor_cb_changed',value)
        self.refreshRequired = True
        try:
            self.settings['Frequency scaling factor'] = float(value)
        except:
            print('Failed to convert to float', value)
        return

    def replotButton1Clicked(self):
        debugger.print('Start:: replotButton1Clicked')
        self.refreshRequired = True
        self.plot_frequency_shift = False
        self.lastButtonPressed = self.replotButton1Clicked
        self.refresh()
        debugger.print('Finished:: replotButton1Clicked')
        return

    def replotButton2Clicked(self):
        debugger.print('Start:: replotButton2Clicked')
        self.refreshRequired = True
        self.plot_frequency_shift = True
        self.lastButtonPressed = self.replotButton2Clicked
        self.refresh()
        debugger.print('Finished:: replotButton2Clicked')
        return

    def plot(self,experiment,xs,ys,legends,label):
        # Plot the experimental values on the left y-axis
        # Plot all the others in xs, ys on the right x-axis
        debugger.print('Start:: plot')
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
        if self.settings['Independent y-axes']:
            self.subplot2 = self.subplot1.twinx()
        else:
            self.subplot2 = self.subplot1
        cmap_index = 0
        lines = []
        if self.plot_frequency_shift:
            lag = float(self.lag_frequency_le.text())
            scale_calc = scale * self.settings['Frequency scaling factor']
        else:
            lag = 0.0
            # scale_calc = scale
            scale_calc = scale * self.settings['Frequency scaling factor']
        for x,y,legend in zip(xs,ys,legends):
            if y is not None:
               x = np.array(x)
               line, = self.subplot1.plot(lag+scale_calc*x,y,lw=2, color=cmap(cmap_index), label=legend )
               lines.append(line)
               cmap_index += 1
        if len(experiment) > 0:
            # Use the x variables from the previous xs, ys
            line, = self.subplot2.plot(scale*x,experiment,lw=2, color=cmap(cmap_index), label='Experiment' )
            lines.append(line)
        labels = [l.get_label() for l in lines]
        if self.settings['Independent y-axes']:
            self.subplot2.set_ylabel('Experiment')
            self.subplot2.set_ylim(bottom=0.0)
            self.subplot1.set_ylabel('Calculated '+self.plot_label )
        else:
            self.subplot1.set_ylabel(self.plot_label )
        self.subplot1.set_xlabel(xlabel)
        self.subplot1.set_ylim(bottom=0.0)
        self.subplot1.legend(lines, labels, loc='best')
        self.subplot1.set_title(self.settings['Plot title'])
        self.canvas.draw_idle()
        debugger.print('Finished:: plot')

    def fittingButtonClicked(self):
        debugger.print('Start:: fittingButtonClicked')
        self.refreshRequired = True
        if self.calculationInProgress:
            self.fittingButton.setText('Perform fitting')
            self.calculationInProgress = False
        else:
            self.fittingButton.setText('Interupt fitting')
            self.calculationInProgress = True
        debugger.print('replotButton2Clicked',self.refreshRequired)
        self.refresh()
        self.replot()
        final_point = self.optimiseFit()
        self.fittingButton.setText('Perform fitting')
        self.calculationInProgress = False
        debugger.print('Finished:: fittingButtonClicked')
        return

    def optimiseFit(self):
        # Optimise the fit of the first scenario to the experimental data
        # First determine who many variables we have
        debugger.print('Start:: optimiseFit')
        self.functionCalls = 0
        self.fit_list = []
        for mode,fitted in enumerate(self.modes_fitted):
            if fitted:
                self.fit_list.append(mode)
        initial_point = [ self.sigmas_cm1[i] for i in self.fit_list ]
        # Append a scaling option
        if self.settings['Optimise frequency scaling']:
            initial_point.append(self.settings['Frequency scaling factor'])
        nvariables = len(initial_point)
        if nvariables > 0:
            final_point = minimize(self.optimiseFunction, initial_point, method='nelder-mead', options={'xatol':0.01, 'disp':True, 'maxfev':nvariables+nvariables*self.settings['Number of iterations']} )
        else: 
            print('No sigmas have been selected for optimisation')
            final_point = []
        debugger.print('Finished:: optimiseFit')
        return final_point

    def optimiseFunction(self,variables) :
        # Determine the function to be optimised (minimised)
        # print('Current point in optimisation')
        # print(variables)
        if not self.calculationInProgress:
            return -9.9E99
        debugger.print('Start:: optimiseFunction',variables)
        self.functionCalls += 1
        nvariables = len(variables)
        self.fittingButton.setText('Interrupt fitting ({}/{})'.format(self.functionCalls,nvariables+1+nvariables*self.settings['Number of iterations']))
        if self.settings['Optimise frequency scaling']:
            sigmas = variables[:-1]
            scaling_factor = variables[-1]
            self.settings['Frequency scaling factor'] = scaling_factor
        else:
            sigmas = variables
            scaling_factor = self.settings['Frequency scaling factor']
        for index,sigma in zip(self.fit_list,sigmas):
            self.sigmas_cm1[index] = sigma
            self.redraw_sigmas_tw()
            self.notebook.settingsTab.sigmas_cm1[index] = sigma
            self.notebook.settingsTab.redraw_output_tw()
            self.notebook.settingsTab.requestRefresh()
        self.refresh(force=True)
        self.replot()
        # Returning the best correlation but made negative because we need to minimise
        if self.settings['Fitting type'] == 'Minimise x-correlation':
            function_value = -1.0*self.xcorr0
        elif self.settings['Fitting type'] == 'Minimise spectral difference':
            function_value = self.rmse
        debugger.print('optimiseFunction - xcorr0,rmse',self.xcorr0,self.rmse)
        debugger.print('Finished:: optimiseFunction',function_value)
        return function_value

    def calculateSpectralDifference(self,scaling_factor):
        debugger.print('calculateSpectralDifference',scaling_factor)
        # Calculate the spectral difference  between the experimental and the first scenario
        if len(self.experimental_spectrum) == 0:
            return 0.0
        debugger.print('Start:: optimiseFunction')
        # col1 contains the experimental spectrum
        col1 = np.array(self.experimental_spectrum)
        maxcol1 = np.max(col1)
        col1 = col1 / maxcol1
        col1[ col1< self.settings['Spectral difference threshold'] ] = 0.0
#        for i,val in enumerate(col1):
#            if val < self.settings['Spectral difference threshold']:
#                col1[i] = 0.0
        # col2 contains the calculated spectrum
        col2 = np.array(self.calculated_spectrum)
        # The new xaxis for the calculated spectrum is scaling_factor*xaxis
        f = interp1d(scaling_factor*np.array(self.xaxis), col2, kind='cubic',fill_value='extrapolate')
        col2 = f(self.xaxis)
        col2 = col2 / maxcol1
        col1[ col1< self.settings['Spectral difference threshold'] ] = 0.0
        diff = col1 - col2
        rmse = np.sqrt(np.dot(diff,diff)/len(col2))
        debugger.print('rmse',rmse)
        debugger.print('Finished:: optimiseFunction')
        return rmse

    def calculateCrossCorrelation(self,scaling_factor):
        debugger.print('Start:: calculateCrossCorrelation',scaling_factor)
        # Calculate the cross correlation coefficient between the experimental and the first scenario
        if len(self.experimental_spectrum) == 0:
            debugger.print('calculateCrossCorrelation experimental_spectrum is not defined')
            debugger.print('Finshed:: calculateCrossCorrelation',scaling_factor)
            return (0.0,0.0,0.0)
        # col1 contains the experimental spectrum
        col1 = np.array(self.experimental_spectrum)
        col1 = ( col1 - np.mean(col1)) / ( np.std(col1) * np.sqrt(len(col1)) )
        # col2 contains the calculated spectrum
        col2 = np.array(self.calculated_spectrum)
        # The new xaxis for the calculated spectrum is scaling_factor*xaxis
        f = interp1d(scaling_factor*np.array(self.xaxis), col2, kind='cubic',fill_value='extrapolate')
        col2 = f(self.xaxis)
        col2 = ( col2 - np.mean(col2)) / ( np.std(col2) * np.sqrt(len(col2)) )
        correlation = np.correlate(col1, col2,  mode='full')
        lag = np.argmax(correlation) - (len(correlation)-1)/2
        lag = (self.xaxis[1] - self.xaxis[0]) * lag
        debugger.print('lag , max(corr), index', lag,np.max(correlation),correlation[int((len(correlation)-1)/2)])
        debugger.print('Finshed:: calculateCrossCorrelation',scaling_factor)
        return (lag,np.max(correlation),correlation[int((len(correlation)-1)/2)])


    def redraw_sigmas_tw(self):
        debugger.print('Start:: redraw_sigmas_tw')
        if len(self.sigmas_cm1) <= 0:
            debugger.print('Finished:: redraw_sigmas_tw')
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
                item.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
                self.sigmas_tw.setItem(i, j, item )
        self.sigmas_tw.resizeColumnsToContents()
        # Release the block on signals for the frequency output table
        self.sigmas_tw.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print('Finished:: redraw_sigmas_tw')

    def on_fitting_type_cb_activated(self,index):
        # Change in fitting type
        debugger.print('on_fitting_type_cb_activated', index)
        self.refreshRequired = True
        self.settings['Fitting type'] = self.fitting_type_definitions[index]

    def on_scenario_cb_activated(self,index):
        # Change in Scenario to be used for fitting
        debugger.print('on_scenario_cb_activated', index)
        self.refreshRequired = True
        self.settings['Scenario index'] = index

    def on_spectrafile_le_return(self):
        # Handle a return in the excel file name line editor
        debugger.print('Start:: on_spectrafile_le_return')
        file_name = self.spectrafile_le.text()
        if not os.path.isfile(file_name):
            qfd = QFileDialog(self)
            qfd.setDirectory(self.notebook.mainTab.directory)
            file_name, _ = qfd.getOpenFileName(self,'Open Excel file','','Excel (*.xlsx);; All Files (*)')
        if not os.path.isfile(file_name):
            debugger.print('Finished:: on_spectrafile_le_return')
            return
        self.settings['Excel file name'] = file_name
        self.spectrafile_le.setText(self.settings['Excel file name'])
        debugger.print('new file name', self.settings['Excel file name'])
        self.refreshRequired = True
        self.excel_file_has_been_read = False
        # redo the plot if a return is pressed
        self.lastButtonPressed()
        debugger.print('Finished:: on_spectrafile_le_return')
        return

    def read_excel_file(self):
        # 
        from openpyxl import load_workbook
        debugger.print('Start:: read_excel_file',self.settings['Excel file name'])
        file_name = self.settings['Excel file name']
        if not os.path.isfile(file_name):
            debugger.print('Finished:: read_excel_file does not exist',self.settings['Excel file name'])
            return
        if self.excel_file_has_been_read:
            debugger.print('Finished:: read_excel_file has been read')
            return
        # 
        # Open the work book
        #
        self.wb = load_workbook(filename=self.settings['Excel file name'], read_only=True)
        self.ws = self.wb.worksheets[0]
        self.excel_frequencies = []
        self.excel_spectrum = []
        for row in self.ws.rows:
            if isinstance(row[0].value, (int, float, complex)):
                self.excel_frequencies.append(row[0].value)
                self.excel_spectrum.append(row[1].value)
        self.excel_file_has_been_read = True
        self.refreshRequired = True
        debugger.print('Finished:: read_excel_file')
        return

    def on_spectrafile_le_changed(self,text):
        debugger.print('on_spectrafile_le_changed', text)
        text = self.spectrafile_le.text()
        self.refreshRequired = True
        self.settings['Excel file name'] = text
        self.excel_file_has_been_read = False

    def on_sigmas_tw_itemChanged(self, item):
        debugger.print('Start:: on_sigmas_tw_itemChanged', item)
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
            try:
                new_value = float(item.text())
                self.sigmas_cm1[row] = new_value
                self.redraw_sigmas_tw()
                self.notebook.settingsTab.sigmas_cm1[row] = new_value
                self.notebook.settingsTab.requestRefresh()
                self.notebook.settingsTab.redraw_output_tw()
            except:
                 print('Failed to convert to float',item.txt())
        elif col == 1:
            self.redraw_sigmas_tw()
        else:
            self.redraw_sigmas_tw()
        self.refreshRequired = True
        QCoreApplication.processEvents()
        debugger.print('Finished:: on_sigmas_tw_itemChanged', item)
        return


    def print_settings(self):
        debugger.print('print_settings')
        for key in self.settings:
            debugger.print(key, self.settings[key])
        return

    def replot(self):
        debugger.print('Start:: replot')
        if len(self.excel_spectrum) > 0:
            self.resample_experimental_spectrum()
        self.plot(self.experimental_spectrum,self.calculated_frequencies,self.calculated_spectra,self.scenario_legends,self.plot_label)
        debugger.print('Finished:: replot')
        return

    def refresh(self,force=False):
        debugger.print('Start:: refresh', force)
        if not self.refreshRequired and not force:
            debugger.print('refresh aborted', self.refreshRequired,force)
            debugger.print('Finished:: refresh', force)
            return
        self.frequencies_cm1 = self.notebook.settingsTab.frequencies_cm1
        if np.sum(self.frequencies_cm1) < 1.0E-10:
            debugger.print('refresh aborted there are no frequencies')
            debugger.print('Finished:: refresh', force)
            return
        #
        # Flag all the scenarios as needing an update
        #
        for scenario in self.notebook.scenarios:
            scenario.requestRefresh()
        self.notebook.settingsTab.requestRefresh()
        #
        # Now refresh the plotting tab 
        #
        self.notebook.plottingTab.refresh()
        #
        # Block signals during refresh
        # 
        #for w in self.findChildren(QWidget):
        #    w.blockSignals(True)
        # 
        # use the settings values to initialise the widgets
        #
        self.spectrafile_le.setText(self.settings['Excel file name'])
        self.read_excel_file()
        self.plot_label = self.notebook.plottingTab.settings['Plot type']
        self.scenario_legends = [ scenario.settings['Legend'] for scenario in self.notebook.scenarios ]
        self.scenario_cb.clear()
        self.scenario_cb.addItems(self.scenario_legends)
        self.scenario_cb.setCurrentIndex(self.settings['Scenario index'])
        vs_cm1 = self.notebook.plottingTab.vs_cm1
        self.calculated_spectra = [ scenario.get_result(vs_cm1, self.plot_label) for scenario in self.notebook.scenarios ]
        self.calculated_spectrum = self.calculated_spectra[self.settings['Scenario index']]
        debugger.print('refresh scenario index' , self.settings['Scenario index'])
        self.iterations_sb.setValue(self.settings['Number of iterations'])
        self.frequency_scaling_factor_sb.setValue(self.settings['Frequency scaling factor'])
        if self.settings['Independent y-axes']:
            self.independent_yaxes_cb.setCheckState(Qt.Checked)
        else:
            self.independent_yaxes_cb.setCheckState(Qt.Unchecked)
        # 
        # If the sigmas are not set return
        #
        self.sigmas_cm1 = self.notebook.settingsTab.sigmas_cm1
        if len(self.sigmas_cm1) < 1:
            return
        self.modes_selected = self.notebook.settingsTab.modes_selected
        self.intensities = self.notebook.settingsTab.intensities
        if len(self.modes_fitted) == 0 and len(self.modes_selected) > 0:
            self.modes_fitted = [ False  for _ in self.modes_selected ]
        self.redraw_sigmas_tw()
        # Resample the spectrum
        self.calculated_frequencies = [ vs_cm1 for scenario in self.notebook.scenarios ]
        if len(self.excel_spectrum) > 0:
            self.resample_experimental_spectrum()
        scaling_factor = self.settings['Frequency scaling factor']
        self.lag,self.xcorr0,self.xcorr1 = self.calculateCrossCorrelation(scaling_factor)
        self.rmse = self.calculateSpectralDifference(scaling_factor)
        self.cross_correlation_le.setText('{:6.4f}'.format(self.xcorr0))
        self.lag_frequency_le.setText('{:8.2f}'.format(self.lag))
        self.frequency_scaling_le.setText('{:8.2f}'.format(self.settings['Frequency scaling factor']))
        self.rmse_le.setText('{:.2e}'.format(self.rmse))
        self.replot()
        #
        # Unblock signals after refresh
        # 
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        self.refreshRequired = False
        debugger.print('Finished:: refresh', force)
        return

    def requestRefresh(self):
        """Request a refresh of the interface"""
        debugger.print('Start:: requestRefresh')
        self.refreshRequired = True
        debugger.print('Finished:: requestRefresh')

    def resample_experimental_spectrum(self):
    #
    #  The experimental spectrum needs to be in the same range as the calculated spectrum
    #  It also needs to be calculated at the same frequencies as the calculated spectrum
    #
        debugger.print('Start:: Resample_experimental_spectrum')
        self.xaxis = np.array(self.calculated_frequencies[0])
        # If the experimental frequencies starts at a higher frequency 
        # than the calculated frequencies then add new frequencies to pad the range out
        #indices = np.where (self.xaxis < self.excel_frequencies[0])
        indices = self.xaxis < self.excel_frequencies[0]
        padded_xaxis = self.xaxis[indices]
        padded_yaxis = np.array([ 0 for index in indices if index ])
        # Add the experimental frequencies
        padded_xaxis = np.append(padded_xaxis,self.excel_frequencies)
        padded_yaxis = np.append(padded_yaxis,self.excel_spectrum)
        # If the experimental frequencies ends at a lower frequency 
        # than the calculated frequencies then add new frequencies to pad the range out
        indices = self.xaxis > self.excel_frequencies[-1]
        padded_xaxis = np.append(padded_xaxis,self.xaxis[indices])
        padded_yaxis = np.append(padded_yaxis, np.array([ 0 for index in indices if index ]) )
        # 
        # Create a function using the padded frequencies to calculate the spectrum at the calculated frequencies
        f = interp1d(padded_xaxis, padded_yaxis, kind='cubic',fill_value='extrapolate')
        # Store the experimental spectrum at the calculated frequencies
        experimental_spectrum = f(self.xaxis)
        if self.settings['Baseline removal']:
            # Apply a Hodrick-Prescott filter to remove the background
            self.experimental_spectrum = Calculator.hodrick_prescott_filter(
                                          experimental_spectrum, 0.01,
                                          self.settings['HPFilter lambda'], 10)
        else:
            self.experimental_spectrum = experimental_spectrum
        debugger.print('Finished:: Resample_experimental_spectrum')
        return

