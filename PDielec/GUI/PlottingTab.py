import os.path
import os
import numpy as np
from PyQt5.QtWidgets  import  QPushButton, QWidget
from PyQt5.QtWidgets  import  QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets  import  QProgressBar, QApplication
from PyQt5.QtWidgets  import  QVBoxLayout, QHBoxLayout, QFormLayout
from PyQt5.QtWidgets  import  QSpinBox,QDoubleSpinBox
from PyQt5.QtWidgets  import  QSizePolicy
from PyQt5.QtCore     import  QCoreApplication, Qt
from PDielec.Constants import  PI, avogadro_si, angstrom
# Import plotting requirements
import matplotlib
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PDielec.Utilities import Debug

class PlottingTab(QWidget):
    def __init__(self, parent, debug=False ):
        super(QWidget, self).__init__(parent)
        global debugger
        debugger = Debug(debug,'PlottingTab')
        debugger.print('Start:: Plotting tab initialisation')
        self.settings = {}
        self.refreshRequired = True
        self.subplot = None
        self.setWindowTitle('Plotting')
        self.settings['Minimum frequency'] = 0
        self.settings['Maximum frequency'] = 200
        self.settings['Frequency increment'] = 0.2
        self.molar_definitions = ['Unit cells','Atoms','Molecules']
        self.settings['Molar definition'] = 'Unit cells'
        self.settings['Number of atoms'] = 1
        self.settings['Plot type'] = 'Powder Molar Absorption'
        # self.settings['Plot title'] = 'Plot Title'
        self.legends = []
        self.vs_cm1 = []
        self.frequency_units = None
        self.molar_cb_current_index = 0
        # store the notebook
        self.notebook = parent
        # get the reader from the main tab
        self.reader = self.notebook.mainTab.reader
        # Create last tab - PlottingTab
        vbox = QVBoxLayout()
        form = QFormLayout()
        #
        # The minimum frequency
        #
        self.vmin_sb = QDoubleSpinBox(self)
        self.vmin_sb.setRange(0,9000)
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmin_sb.setToolTip('Set the minimum frequency to be considered)')
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        #
        # The maximum frequency
        #
        self.vmax_sb = QDoubleSpinBox(self)
        self.vmax_sb.setRange(0,9000)
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.vmax_sb.setToolTip('Set the maximum frequency to be considered)')
        self.vmax_sb.valueChanged.connect(self.on_vmax_changed)
        #
        # Choose a suitable increment
        #
        self.vinc_sb = QDoubleSpinBox(self)
        self.vinc_sb.setRange(0.0001,5.0)
        self.vinc_sb.setSingleStep(0.1)
        self.vinc_sb.setDecimals(4)
        self.vinc_sb.setToolTip('Choose an increment for the frequency when plotting')
        self.vinc_sb.setValue(self.settings['Frequency increment'])
        self.vinc_sb.valueChanged.connect(self.on_vinc_changed)
        #
        label = QLabel('Frequency min, max and increment', self)
        label.setToolTip('Choose minimum, maximum and increment for the frequency when plotting')
        #
        hbox = QHBoxLayout()
        hbox.addWidget(self.vmin_sb)
        hbox.addWidget(self.vmax_sb)
        hbox.addWidget(self.vinc_sb)
        form.addRow(label, hbox)
        #
        # Define molar quantity
        #
        self.molar_cb = QComboBox(self)
        self.molar_cb.setToolTip('Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given')
        self.molar_cb.addItems(self.molar_definitions)
        try:
            self.molar_cb_current_index = self.molar_definitions.index(self.settings["Molar definition"])
        except:
            self.molar_cb_current_index = 0
            self.settings["Molar definition"] = self.molar_definitions[self.molar_cb_current_index]
        self.molar_cb.setCurrentIndex(self.molar_cb_current_index)
        self.molar_cb.activated.connect(self.on_molar_cb_activated)
        label = QLabel('Molar definition', self)
        label.setToolTip('Define what a mole is.  \nIn the case of Molecules, the number of atoms in a molecule must be given')
        form.addRow(label, self.molar_cb)
        #
        # Number of atoms in a molecule
        #
        self.natoms_sb = QSpinBox(self)
        self.natoms_sb.setToolTip('Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed')
        self.natoms_sb.setRange(1,500)
        self.natoms_sb.setValue(self.settings['Number of atoms'])
        self.natoms_sb.valueChanged.connect(self.on_natoms_changed)
        self.natoms_sb.setEnabled(False)
        label = QLabel('Number of atoms per molecule', self)
        label.setToolTip('Set the number of atoms in a molecule. \nOnly need this if moles of molecules is needed')
        form.addRow(label, self.natoms_sb)
        #jk #
        #jk # Set the plot title
        #jk #
        #jk self.title_le = QLineEdit(self)
        #jk self.title_le.setToolTip('Set the plot title')
        #jk self.title_le.setText(self.settings['Plot title'])
        #jk self.title_le.textChanged.connect(self.on_title_changed)
        #jk label = QLabel('Plot title', self)
        #jk label.setToolTip('Set the plot title')
        #jk form.addRow(label, self.title_le)
        #
        # Set the x-axis frequency units
        #
        self.funits_cb = QComboBox(self)
        self.funits_cb.setToolTip('Set the frequency units for the x-axis')
        self.funits_cb.addItems( ['wavenumber','THz'] )
        self.frequency_units = 'wavenumber'
        self.funits_cb.activated.connect(self.on_funits_cb_activated)
        label = QLabel('Frequency units for the x-axis', self)
        label.setToolTip('Set the frequency units for the x-axis')
        form.addRow(label, self.funits_cb)
        #
        # Final button
        #
        self.plot_type_cb = QComboBox(self)
        self.plot_type_cb.setToolTip('Choose the which data to plot')
        self.plot_types = [
                            'Powder Molar Absorption',
                            'Powder Absorption',
                            'Powder Real Permittivity',
                            'Powder Imaginary Permittivity',
                            'Powder ATR',
                            'Crystal Reflectance (P polarisation)',
                            'Crystal Reflectance (S polarisation)',
                            'Crystal Transmittance (P polarisation)',
                            'Crystal Transmittance (S polarisation)',
                            'Crystal Absorbtance (P polarisation)',
                            'Crystal Absorbtance (S polarisation)',
                          ]
        self.plot_ylabels = {
                     'Powder Molar Absorption': r'Molar Absorption Coefficient $\mathdefault{(L mole^{-1} cm^{-1})}$',
                           'Powder Absorption': r'Absorption Coefficient $\mathdefault{(cm^{-1})}$',
                    'Powder Real Permittivity': r'Real Component of Permittivity',
               'Powder Imaginary Permittivity': r'Imaginary Component of Permittivity',
                                  'Powder ATR': r'ATR absorption',
        'Crystal Reflectance (P polarisation)': r'Fraction of p-polarised reflectance',
        'Crystal Reflectance (S polarisation)': r'Fraction of s-polarised reflectance',
      'Crystal Transmittance (P polarisation)': r'Fraction of p-polarised transmitted',
      'Crystal Transmittance (S polarisation)': r'Fraction of s-polarised transmitted',
        'Crystal Absorbtance (P polarisation)': r'Fraction of p-polarised absorbtance',
        'Crystal Absorbtance (S polarisation)': r'Fraction of s-polarised absorbtance',
                            }

        self.plot_type_cb.activated.connect(self.on_plot_type_cb_activated)
        self.plot_type_cb.addItems( self.plot_types )
        label = QLabel('Choose plot type', self)
        label.setToolTip('Choose the plot type')
        index = self.plot_type_cb.findText(self.settings['Plot type'], Qt.MatchFixedString)
        self.plot_type_cb.setCurrentIndex(index)
        plot_button = QPushButton('Update plot')
        plot_button.clicked.connect(self.refresh)
        plot_button.setToolTip('Update the plot')
        hbox = QHBoxLayout()
        hbox.addWidget(self.plot_type_cb)
        hbox.addWidget(plot_button)
        form.addRow(label, hbox)
        # Add a progress bar
        self.progressbar = QProgressBar(self)
        self.progressbar.setToolTip('Show the progress of any calculations')
        # Append the progress bar to the list of progress bars managed by the notebook
        self.notebook.progressbars_add(self.progressbar)
        self.notebook.progressbars_set_maximum(0)
        label = QLabel('Calculation progress', self)
        label.setToolTip('Show the progress of any calculations')
        form.addRow(label,self.progressbar)
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
        QCoreApplication.processEvents()
        # Create the plot
        debugger.print('Finished:: Plotting tab initialisation')
        return

    def requestRefresh(self):
        debugger.print('Start:: requestRefresh')
        self.refreshRequired
        debugger.print('Finished:: requestRefresh')
        return

    def requestScenarioRefresh(self):
        debugger.print('Start:: requestScenarioRefresh')
        self.notebook.settingsTab.requestRefresh()
        for scenario in self.notebook.scenarios:
            scenario.requestRefresh()
        debugger.print('Finished:: requestScenarioRefresh')
        return

    def on_vinc_changed(self,value):
        debugger.print('Start:: on_vinc_changed', value)
        self.vinc_sb.blockSignals(True)
        value = self.vinc_sb.value()
        self.settings['Frequency increment'] = value
        self.notebook.fitterTab.requestRefresh()
        self.refreshRequired = True
        self.requestScenarioRefresh()
        debugger.print('on_vinc_change ', self.settings['Frequency increment'])
        self.vinc_sb.blockSignals(False)
        debugger.print('Finished:: on_vinc_changed', value)

    def on_vmin_changed(self):
        debugger.print('Start:: on_vmin_changed')
        self.vmin_sb.blockSignals(True)
        vmin = self.vmin_sb.value()
        vmax = self.vmax_sb.value()
        self.settings['Minimum frequency'] = vmin
        debugger.print('on_vmin_changed setting vmin to', self.settings['Minimum frequency'])
        self.notebook.fitterTab.requestRefresh()
        self.refreshRequired = True
        self.requestScenarioRefresh()
        self.vmin_sb.blockSignals(False)
        debugger.print('Finished:: on_vmin_changed')

    def on_vmax_changed(self):
        debugger.print('Start:: on_vmax_changed')
        self.vmax_sb.blockSignals(True)
        vmin = self.vmin_sb.value()
        vmax = self.vmax_sb.value()
        self.settings['Maximum frequency'] = vmax
        debugger.print('on_vmax_changed setting vmax to ', self.settings['Maximum frequency'])
        self.notebook.fitterTab.requestRefresh()
        self.refreshRequired = True
        self.requestScenarioRefresh()
        self.vmax_sb.blockSignals(False)
        debugger.print('Finished:: on_vmax_changed')

    def refresh(self,force=False):
        debugger.print('Start:: refresh', force)
        if not self.refreshRequired and not force:
            debugger.print('Finished:: refreshing widget not required')
            return
        #
        # Block signals during refresh
        #
        self.greyed_out()
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Now refresh values
        if self.settings['Maximum frequency'] < self.settings['Minimum frequency']:
            self.settings['Maximum frequency'] = self.settings['Minimum frequency']+1
        if self.settings['Frequency increment'] > self.settings['Maximum frequency'] - self.settings['Minimum frequency']:
            self.settings['Frequency increment'] = (self.settings['Maximum frequency'] - self.settings['Minimum frequency'])/2
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmax_sb.setValue(self.settings['Maximum frequency'])
        self.vinc_sb.setValue(self.settings['Frequency increment'])
        index = self.plot_type_cb.findText(self.settings['Plot type'], Qt.MatchFixedString)
        self.plot_type_cb.setCurrentIndex(index)
        try:
            self.molar_cb_current_index = self.molar_definitions.index(self.settings["Molar definition"])
        except:
            self.molar_cb_current_index = 0
            self.settings["Molar definition"] = self.molar_definitions[self.molar_cb_current_index]
        self.molar_cb.setCurrentIndex(self.molar_cb_current_index)
        self.natoms_sb.setValue(self.settings['Number of atoms'])
        # Refresh the widgets that depend on the reader
        self.reader = self.notebook.reader
        if self.reader is not None:
            self.set_concentrations()
        # Reset the progress bar
        self.notebook.progressbars_set_maximum(0)
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        debugger.print('calling plot from refresh')
        self.plot()
        refreshRequired = False
        debugger.print('Finished:: refresh', force)
        return

    def on_natoms_changed(self, value):
        debugger.print('Start:: on_natoms_changed', value)
        self.settings['Number of atoms'] = value
        debugger.print('on natoms changed ', self.settings['Number of atoms'])
        self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 * self.settings['Number of atoms'] / self.reader.nions)
        debugger.print('The concentration has been set', self.settings['Molar definition'], self.settings['concentration'])
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refresh()
        debugger.print('Finished:: on_natoms_changed', value)

    def on_plot_type_cb_activated(self, index):
        debugger.print('Start:: on_plot_type_cb_activated', index)
        self.settings['Plot type'] = self.plot_type_cb.currentText()
        debugger.print('Changed plot type to ', self.settings['Plot type'])
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refresh()
        debugger.print('Finished:: on_plot_type_cb_activated', index)

    def on_funits_cb_activated(self, index):
        debugger.print('Start:: on_funits_cb_activated', index)
        if index == 0:
            self.frequency_units = 'wavenumber'
        else:
            self.frequency_units = 'THz'
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refresh()
        debugger.print('Frequency units changed to ', self.frequency_units)
        debugger.print('Finished:: on_funits_cb_activated', index)

    def on_molar_cb_activated(self, index):
        debugger.print('Start:: on_molar_cb_activated', index)
        self.molar_cb_current_index = index
        self.settings['Molar definition'] = self.molar_definitions[index]
        self.set_concentrations()
        self.refreshRequired = True
        self.notebook.fitterTab.requestRefresh()
        self.refresh()
        debugger.print('The concentration has been set', self.settings['Molar definition'], self.settings['concentration'])
        debugger.print('Finished:: on_molar_cb_activated', index)
        return

    def set_concentrations(self):
        debugger.print('Start:: set_concentration')
        if self.settings['Molar definition'] == 'Molecules':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 * self.settings['Number of atoms'] / self.reader.nions)
            self.natoms_sb.setEnabled(True)
        elif self.settings['Molar definition'] == 'Unit cells':
            self.settings['concentration']      = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
            self.settings['cell concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
            self.natoms_sb.setEnabled(False)
        elif self.settings['Molar definition'] == 'Atoms':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 / self.reader.nions)
            self.natoms_sb.setEnabled(False)
        debugger.print('Finished:: set_concentration')
        return

    def writeSpreadsheet(self):
        debugger.print('Start::writeSpreadsheet')
        if self.notebook.spreadsheet is None:
            debugger.print('Finished::writeSpreadsheet spreadsheet is None')
            return
        # make sure the plottingTab is up to date
        self.refresh()
        # Handle powder plots
        molarAbsorptionCoefficients = []
        absorptionCoefficients      = []
        realPermittivities          = []
        imagPermittivities          = []
        sp_atrs                     = []
        R_ps                        = []
        R_ss                        = []
        T_ps                        = []
        T_ss                        = []
        powder_legends              = []
        crystal_legends             = []
        # Deal with Scenarios 
        sp = self.notebook.spreadsheet
        sp.selectWorkSheet('Scenarios')
        sp.delete()
        sp.writeNextRow(['A list of the scenarios used the calculation of the effective medium'],col=1)
        for index,scenario in enumerate(self.notebook.scenarios):
            if scenario.scenarioType == 'Powder':
                direction = scenario.direction
                depolarisation = scenario.depolarisation
                sp.writeNextRow([''],col=1)
                sp.writeNextRow(['Scenario '+str(index)],col=1,check=1)
                settings = scenario.settings
                for key in sorted(settings,key=str.lower):
                    sp.writeNextRow([key, settings[key]],col=1,check=1)
                sp.writeNextRow(['Normalised unique direction']+direction.tolist(), col=1,check=1)
                sp.writeNextRow(['Depolarisation matrix'], col=1,check=1)
                sp.writeNextRow(depolarisation[0].tolist(), col=2, check=1)
                sp.writeNextRow(depolarisation[1].tolist(), col=2, check=1)
                sp.writeNextRow(depolarisation[2].tolist(), col=2, check=1)
                molarAbsorptionCoefficients.append( scenario.get_result(self.vs_cm1,self.plot_types[0] ) )
                absorptionCoefficients.append( scenario.get_result(self.vs_cm1,self.plot_types[1] ) )
                realPermittivities.append( scenario.get_result(self.vs_cm1,self.plot_types[2] ) )
                imagPermittivities.append( scenario.get_result(self.vs_cm1,self.plot_types[3] ) )
                sp_atrs.append( scenario.get_result(self.vs_cm1,self.plot_types[4] ) )
                powder_legends.append(scenario.settings['Legend'])
            else:
                sp.writeNextRow([''],col=1)
                sp.writeNextRow(['Scenario '+str(index)],col=1,check=1)
                settings = scenario.settings
                for key in sorted(settings,key=str.lower):
                    sp.writeNextRow([key, settings[key]],col=1,check=1)
                sp.writeNextRow(['Crystal axes in the laboratory frame'], col=1,check=1)
                sp.writeNextRow(scenario.labframe_a.tolist(), col=2, check=1)
                sp.writeNextRow(scenario.labframe_b.tolist(), col=2, check=1)
                sp.writeNextRow(scenario.labframe_c.tolist(), col=2, check=1)
                # Store the reflectance and transmittance
                R_ps.append( scenario.get_result(self.vs_cm1,self.plot_types[5] ) )
                R_ss.append( scenario.get_result(self.vs_cm1,self.plot_types[6] ) )
                T_ps.append( scenario.get_result(self.vs_cm1,self.plot_types[7] ) )
                T_ss.append( scenario.get_result(self.vs_cm1,self.plot_types[8] ) )
                crystal_legends.append(scenario.settings['Legend'])
        # Single crystal Permittivity
        dielecv = self.notebook.settingsTab.get_crystal_permittivity(self.vs_cm1)
        # Powder results
        # Work out what molar units we are using
        if len(molarAbsorptionCoefficients) > 0:
            if self.settings['Molar definition'] == 'Molecules':
                sheet_name = 'Powder Molar Absorption (mols)'
            elif self.settings['Molar definition'] == 'Unit cells':
                sheet_name = 'Powder Molar Absorption (cells)'
            elif self.settings['Molar definition'] == 'Atoms':
                sheet_name = 'Powder Molar Absorption (atoms)'
            # Always write out the moles of cell
            self.write_powder_results(sp, 'Powder Molar Absorption (cells)', self.vs_cm1, powder_legends, molarAbsorptionCoefficients)
            if not self.settings['Molar definition'] == 'Unit cells':
                # If some other molar definition has been used then write that out too
                molarAbsorptionCoefficients_mols = []
                molar_scaling = self.settings['cell concentration']/self.settings['concentration']
                for absorption in molarAbsorptionCoefficients:
                    molarAbsorptionCoefficients_mols.append(molar_scaling * np.array(absorption))
                self.write_powder_results(sp, sheet_name,                      self.vs_cm1, powder_legends, molarAbsorptionCoefficients_mols)
            # end if
            self.write_powder_results(sp, 'Powder Absorption',             self.vs_cm1, powder_legends, absorptionCoefficients)
            self.write_powder_results(sp, 'Powder Real Permittivity',      self.vs_cm1, powder_legends, realPermittivities)
            self.write_powder_results(sp, 'Powder Imaginary Permittivity', self.vs_cm1, powder_legends, imagPermittivities)
            self.write_powder_results(sp, 'Powder ATR Reflectance',        self.vs_cm1, powder_legends, sp_atrs)
        # Single Crystal results
        if len(R_ps) > 0:
            self.write_crystal_results(sp, 'Crystal R_p', self.vs_cm1, crystal_legends, R_ps)
            self.write_crystal_results(sp, 'Crystal R_s', self.vs_cm1, crystal_legends, R_ss)
            self.write_crystal_results(sp, 'Crystal T_p', self.vs_cm1, crystal_legends, T_ps)
            self.write_crystal_results(sp, 'Crystal T_s', self.vs_cm1, crystal_legends, T_ss)

        if len(dielecv) > 0:
            self.write_eps_results(sp, self.vs_cm1, dielecv)
        debugger.print('Finished::writeSpreadsheet')
        return

    def write_eps_results(self, sp, vs, dielecv):
        debugger.print('Start:: write_eps_results length vs',len(vs))
        sp.selectWorkSheet('Real Crystal Permittivity')
        sp.delete()
        headers = ['frequencies (cm-1)', 'xx', 'yy', 'zz', 'xy', 'xz', 'yz' ]
        sp.writeNextRow(headers,row=0, col=1)
        for v,eps in zip(vs,dielecv):
            eps_xx_r = np.real(eps[0][0])
            eps_yy_r = np.real(eps[1][1])
            eps_zz_r = np.real(eps[2][2])
            eps_xy_r = np.real(eps[0][1])
            eps_xz_r = np.real(eps[0][2])
            eps_yz_r = np.real(eps[1][2])
            output = [v, eps_xx_r, eps_yy_r, eps_zz_r, eps_xy_r, eps_xz_r, eps_yz_r ]
            sp.writeNextRow(output, col=1,check=1)
        sp.selectWorkSheet('Imag Crystal Permittivity')
        sp.delete()
        sp.writeNextRow(headers,row=0, col=1)
        for v,eps in zip(vs,dielecv):
            eps_xx_i = np.imag(eps[0][0])
            eps_yy_i = np.imag(eps[1][1])
            eps_zz_i = np.imag(eps[2][2])
            eps_xy_i = np.imag(eps[0][1])
            eps_xz_i = np.imag(eps[0][2])
            eps_yz_i = np.imag(eps[1][2])
            output = [v, eps_xx_i, eps_yy_i, eps_zz_i, eps_xy_i, eps_xz_i, eps_yz_i ]
            sp.writeNextRow(output, col=1,check=1)
        debugger.print('Finished:: write_eps_results length vs',len(vs))
        return

    def write_crystal_results(self, sp, name, vs, legends, yss):
        """ 
        sp        is the spreadsheet object
        name      is the worksheet name used for writing
        vs        an np.array of the frequencies
        yss       a list of np.arrays of the reflections and transmittance ] 
        headings  the heading names for the yss
        """
        debugger.print('Start:: write_crystal_results')
        debugger.print('write_crystal_results name',name)
        debugger.print('write_crystal_results legends',legends)
        debugger.print('write_crystal_results length vs',len(vs))
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ['frequencies (cm-1)']
        headers.extend(legends)
        sp.writeNextRow(headers,row=0, col=1)
        for iv,v in enumerate(vs):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)
        debugger.print('Finished:: write_crystal_results')
        return


    def write_powder_results(self, sp, name, vs, legends, yss):
        debugger.print('Start:: write powder results')
        debugger.print('write_powder_results name',name)
        debugger.print('write_powder_results legends',legends)
        debugger.print('write_powder_results length vs',len(vs))
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ['frequencies (cm-1)']
        #for isc,ys in enumerate(yss):
        #    headers.append('Scenario'+str(isc))
        headers.extend(legends)
        sp.writeNextRow(headers,row=0, col=1)
        for iv,v in enumerate(vs):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)
        debugger.print('Finished:: write powder results')
        return

    def plot(self):
        # import matplotlib.pyplot as pl
        # mp.use('Qt5Agg')
        debugger.print('Start:: plot')
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['Program']
        filename = self.notebook.mainTab.getFullFileName()
        reader = self.notebook.mainTab.reader
        if reader is None:
            debugger.print('Finished:: plot aborting because reader is NONE')
            return
        if program == '':
            debugger.print('Finished:: plot aborting because program is not set')
            return
        if filename == '':
            debugger.print('Finished:: plot aborting because filename is not set')
            return
        if self.notebook.settingsTab.CrystalPermittivity is None:
            debugger.print('Finished:: plot aborting because settingTab.CrystalPermittivity is not set')
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        vmin = self.settings['Minimum frequency']
        vmax = self.settings['Maximum frequency']
        vinc = self.settings['Frequency increment']
        self.vs_cm1 = np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc))
        self.subplot = None
        self.figure.clf()
        if self.frequency_units == 'wavenumber':
            xlabel = r'Frequency $\mathdefault{(cm^{-1})}}$'
            xscale = 1.0
        else:
            xlabel = r'THz'
            xscale = 0.02998
        x = np.array(self.vs_cm1)
        self.subplot = self.figure.add_subplot(111)
        n = len(self.notebook.scenarios)
        if self.notebook.settingsTab.refreshRequired:
            n += 1
        self.notebook.progressbars_set_maximum(n*len(x))
        self.legends = []
        plots = 0
        for scenario in self.notebook.scenarios:
            legend = scenario.settings['Legend']
            self.legends.append(legend)
            y = scenario.get_result(self.vs_cm1,self.settings['Plot type'])
            if y is not None and len(y) > 0:
                y = np.array(y)
                if self.settings['Plot type'] == 'Powder Molar Absorption':
                    y = y * self.settings['cell concentration']/self.settings['concentration']
                plots += 1
                line, = self.subplot.plot(xscale*x,y,lw=2, label=legend )
        if plots > 0:
            self.subplot.set_xlabel(xlabel)
            self.subplot.set_ylabel(self.plot_ylabels[self.settings['Plot type']])
            self.subplot.legend(loc='best')
            self.subplot.set_title(self.settings['Plot type'])
            self.canvas.draw_idle()
        QApplication.restoreOverrideCursor()
        debugger.print('Finished:: plot')

    def greyed_out(self):
        """Handle items that should be greyed out if they are not needed"""
        debugger.print('Start:: greyed_out')
        powder_scenarios_present = False
        crystal_scenarios_present = False
        for scenario in self.notebook.scenarios:
            if scenario.scenarioType == 'Powder':
                powder_scenarios_present = True
            else:
                crystal_scenarios_present = True
        # end of for loop
        # 
        # Disable any plot types that are not needed
        #
        self.plot_type_cb.model().item(0).setEnabled(True)
        self.plot_type_cb.model().item(1).setEnabled(True)
        self.plot_type_cb.model().item(2).setEnabled(True)
        self.plot_type_cb.model().item(3).setEnabled(True)
        self.plot_type_cb.model().item(4).setEnabled(True)
        self.plot_type_cb.model().item(5).setEnabled(True)
        self.plot_type_cb.model().item(6).setEnabled(True)
        self.plot_type_cb.model().item(7).setEnabled(True)
        self.plot_type_cb.model().item(8).setEnabled(True)
        self.plot_type_cb.model().item(9).setEnabled(True)
        self.plot_type_cb.model().item(10).setEnabled(True)
        index = self.plot_type_cb.findText(self.settings['Plot type'], Qt.MatchFixedString)
        if not powder_scenarios_present:
            self.plot_type_cb.model().item(0).setEnabled(False)
            self.plot_type_cb.model().item(1).setEnabled(False)
            self.plot_type_cb.model().item(2).setEnabled(False)
            self.plot_type_cb.model().item(3).setEnabled(False)
            self.plot_type_cb.model().item(4).setEnabled(False)
            if index < 5:
                self.plot_type_cb.setCurrentIndex(5)
                self.settings['Plot type'] = self.plot_type_cb.currentText()
        if not crystal_scenarios_present:
            self.plot_type_cb.model().item(5).setEnabled(False)
            self.plot_type_cb.model().item(6).setEnabled(False)
            self.plot_type_cb.model().item(7).setEnabled(False)
            self.plot_type_cb.model().item(8).setEnabled(False)
            self.plot_type_cb.model().item(9).setEnabled(False)
            self.plot_type_cb.model().item(10).setEnabled(False)
            if index >= 5:
                self.plot_type_cb.setCurrentIndex(0)
                self.settings['Plot type'] = self.plot_type_cb.currentText()
        debugger.print('Finished:: greyed_out')
