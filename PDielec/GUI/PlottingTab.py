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
        debugger.print('Plotting tab initialisaton')
        self.settings = {}
        self.subplot = None
        self.setWindowTitle('Plotting')
        self.settings['Minimum frequency'] = 0
        self.settings['Maximum frequency'] = 200
        self.settings['Frequency increment'] = 0.2
        self.molar_definitions = ['Unit cells','Atoms','Molecules']
        self.settings['Molar definition'] = 'Unit cells'
        self.settings['Number of atoms'] = 1
        self.settings['Plot title'] = 'Plot Title'
        self.legends = []
        self.directions = []
        self.depolarisations = []
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
        self.vmin_sb = QSpinBox(self)
        self.vmin_sb.setRange(0,9000)
        self.vmin_sb.setValue(self.settings['Minimum frequency'])
        self.vmin_sb.setToolTip('Set the minimum frequency to be considered)')
        self.vmin_sb.valueChanged.connect(self.on_vmin_changed)
        #
        # The maximum frequency
        #
        self.vmax_sb = QSpinBox(self)
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
        #
        # Set the plot title
        #
        self.title_le = QLineEdit(self)
        self.title_le.setToolTip('Set the plot title')
        self.title_le.setText(self.settings['Plot title'])
        self.title_le.textChanged.connect(self.on_title_changed)
        label = QLabel('Plot title', self)
        label.setToolTip('Set the plot title')
        form.addRow(label, self.title_le)
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
        self.settings['Plot type'] = self.plot_types[0]
        index = self.plot_type_cb.findText(self.settings['Plot type'], Qt.MatchFixedString)
        self.plot_type_cb.setCurrentIndex(index)
        plot_button = QPushButton('Update plot')
        plot_button.clicked.connect(self.plot)
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
        debugger.print('Calling plot() from initialiser')
        self.plot()

    def on_vinc_changed(self,value):
        self.settings['Frequency increment'] = value
        debugger.print('on vinc change ', self.settings['Frequency increment'])

    def on_vmin_changed(self):
        self.settings['Minimum frequency'] = self.vmin_sb.value()
        self.notebook.fitterTab.dirty = True
        debugger.print('on vmin change ', self.settings['Minimum frequency'])

    def on_vmax_changed(self):
        self.settings['Maximum frequency'] = self.vmax_sb.value()
        self.notebook.fitterTab.dirty = True
        debugger.print('on vmax change ', self.settings['Maximum frequency'])

    def refresh(self,force=False):
        debugger.print('refreshing widget', force)
        #
        # Block signals during refresh
        #
        self.greyed_out()
        for w in self.findChildren(QWidget):
            w.blockSignals(True)
        # Now refresh values
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
        self.title_le.setText(self.settings['Plot title'])
        # Refresh the widgets that depend on the reader
        self.reader = self.notebook.reader
        if self.reader is not None:
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
        # Reset the progress bar
        self.notebook.progressbars_set_maximum(0)
        self.plot()
        #
        # Unblock signals after refresh
        #
        for w in self.findChildren(QWidget):
            w.blockSignals(False)
        QCoreApplication.processEvents()
        return

    def on_natoms_changed(self, value):
        self.settings['Number of atoms'] = value
        debugger.print('on natoms changed ', self.settings['Number of atoms'])
        self.notebook.fitterTab.dirty = True

    def on_plot_type_cb_activated(self, index):
        self.settings['Plot type'] = self.plot_type_cb.currentText()
        debugger.print('Changed plot type to ', self.settings['Plot type'])
        self.notebook.fitterTab.dirty = True
        self.plot()

    def on_funits_cb_activated(self, index):
        if index == 0:
            self.frequency_units = 'wavenumber'
        else:
            self.frequency_units = 'THz'
        self.notebook.fitterTab.dirty = True
        debugger.print('Frequency units changed to ', self.frequency_units)
        self.plot()

    def on_molar_cb_activated(self, index):
        self.molar_cb_current_index = index
        self.settings['Molar definition'] = self.molar_definitions[index]
        if self.settings['Molar definition'] == 'Molecules':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 * self.settings['Number of atoms'] / self.reader.nions)
            self.natoms_sb.setEnabled(True)
        elif self.settings['Molar definition'] == 'Unit cells':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24)
            self.natoms_sb.setEnabled(False)
        elif self.settings['Molar definition'] == 'Atoms':
            self.settings['concentration'] = 1000.0 / (avogadro_si * self.reader.volume * 1.0e-24 / self.reader.nions)
            self.natoms_sb.setEnabled(False)
        self.notebook.fitterTab.dirty = True
        debugger.print('The concentration has been set', self.settings['Molar definition'], self.settings['concentration'])

    def write_spreadsheet(self):
        debugger.print('write spreadsheet')
        if self.notebook.spreadsheet is None:
            return
        # make sure the plottingTab is up to date
        self.notebook.plottingTab.refresh()
        # Handle powder plots
        molarAbsorptionCoefficients = []
        absorptionCoefficients      = []
        realPermittivities          = []
        imagPermittivities          = []
        sp_atrs                     = []
        vss                         = []
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
                vss.append(self.vs_cm1)
        # Now deal with Molar absorption, absorption, real and imaginary permittivity
        if len(molarAbsorptionCoefficients) > 0:
            self.write_results(sp, 'Molar Absorption', vss, molarAbsorptionCoefficients)
            self.write_results(sp, 'Absorption', vss, absorptionCoefficients)
            self.write_results(sp, 'Real Permittivity', vss, realPermittivities)
            self.write_results(sp, 'Imaginary Permittivity', vss, imagPermittivities)
            self.write_results(sp, 'ATR Reflectance', vss, sp_atrs)
        # Handle Single Crystal Plots (A temporary fix for one crystal tab only)
        for index,scenario in enumerate(self.notebook.scenarios):
            if scenario.scenarioType == 'Single crystal':
                sp.selectWorkSheet('Single Crystal')
                sp.delete()
                sp.writeNextRow(['Settings for the single crystal calculation of absorption and reflection'],col=1)
                sp.writeNextRow([''],col=1)
                sp.writeNextRow([ 'Single crystal mode',          scenario.settings['Mode'] ],col=1)
                sp.writeNextRow([ 'Minimum frequency',            self.notebook.plottingTab.settings['Minimum frequency'] ],col=1)
                sp.writeNextRow([ 'Maximum frequency',            self.notebook.plottingTab.settings['Maximum frequency'] ],col=1)
                sp.writeNextRow([ 'Frequency increment',          self.notebook.plottingTab.settings['Frequency increment'] ],col=1)
                sp.writeNextRow([ 'Surface definition (h)',       scenario.settings['Unique direction - h'] ],col=1)
                sp.writeNextRow([ 'Surface definition (k)',       scenario.settings['Unique direction - k'] ],col=1)
                sp.writeNextRow([ 'Surface definition (l)',       scenario.settings['Unique direction - l'] ],col=1)
                sp.writeNextRow([ 'Azimuthal angle',              scenario.settings['Azimuthal angle'] ],col=1)
                sp.writeNextRow([ 'Angle of incidence',           scenario.settings['Angle of incidence'] ],col=1)
                sp.writeNextRow([ 'Superstrate dielectric',       scenario.settings['Superstrate dielectric'] ],col=1)
                sp.writeNextRow([ 'Substrate dielectric',         scenario.settings['Substrate dielectric'] ],col=1)
                sp.writeNextRow([ 'Film thickness(nm)',           scenario.settings['Film thickness'] ],col=1)
                headings = ['R_p', 'R_s', 'T_p', 'T_s']
                self.write_crystal_results(sp, 'Crystal R&T',     self.vs_cm1, [scenario.p_reflectance, scenario.s_reflectance, scenario.p_transmittance, scenario.s_transmittance], headings )

    def write_crystal_results(self, sp, name, vs, yss, headings):
        """ 
        sp        is the spreadsheet object
        name      is the worksheet name used for writing
        vs        an np.array of the frequencies
        yss       a list of np.arrays of the reflections and transmittance ] 
        headings  the heading names for the yss
        """
        debugger.print('write_crystal_results')
        debugger.print('write_crystal_results name',name)
        debugger.print('write_crystal_results headings',headings)
        debugger.print('write_crystal_results length vs',len(vs))
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ['frequencies (cm-1)']
        headers.extend(headings)
        sp.writeNextRow(headers,row=0, col=1)
        for ys in yss:
        for iv,v in enumerate(vs):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)


    def write_results(self, sp, name, vss, yss):
        debugger.print('write results')
        sp.selectWorkSheet(name)
        sp.delete()
        headers = ['frequencies (cm-1)']
        #for isc,ys in enumerate(yss):
        #    headers.append('Scenario'+str(isc))
        headers.extend(self.legends)
        sp.writeNextRow(headers,row=0, col=1)
        for iv,v in enumerate(vss[0]):
           output = [v]
           for ys in yss:
               output.append(ys[iv])
           sp.writeNextRow(output, col=1,check=1)

    def plot(self):
        # import matplotlib.pyplot as pl
        # mp.use('Qt5Agg')
        debugger.print('plot')
        QApplication.setOverrideCursor(Qt.WaitCursor)
        # Assemble the mainTab settings
        settings = self.notebook.mainTab.settings
        program = settings['Program']
        filename = settings['Output file name']
        reader = self.notebook.mainTab.reader
        if reader is None:
            return
        if program == '':
            return
        if filename == '':
            return
        vmin = self.settings['Minimum frequency']
        vmax = self.settings['Maximum frequency']
        vinc = self.settings['Frequency increment']
        self.vs_cm1 = np.arange(float(vmin), float(vmax)+0.5*float(vinc), float(vinc))
        debugger.print('plot')
        self.subplot = None
        self.figure.clf()
        if self.frequency_units == 'wavenumber':
            xlabel = r'Frequency $\mathdefault{(cm^{-1})}}$'
            scale = 1.0
        else:
            xlabel = r'THz'
            scale = 0.02998
        x = np.array(self.vs_cm1)
        self.subplot = self.figure.add_subplot(111)
        self.notebook.progressbars_set_maximum(len(x)*len(self.notebook.scenarios))
        self.legends = []
        for scenario in self.notebook.scenarios:
            legend = scenario.settings['Legend']
            self.legends.append(legend)
            y = scenario.get_result(self.vs_cm1,self.settings['Plot type'])
            if y is not None and len(y) > 0:
                line, = self.subplot.plot(scale*x,y,lw=2, label=legend )
        self.subplot.set_xlabel(xlabel)
        self.subplot.set_ylabel(self.plot_ylabels[self.settings['Plot type']])
        self.subplot.legend(loc='best')
        self.subplot.set_title(self.settings['Plot type'])
        self.canvas.draw_idle()
        QApplication.restoreOverrideCursor()


    def on_title_changed(self,text):
        self.settings['Plot title'] = text
        if self.subplot is not None:
            self.subplot.set_title(self.settings['Plot title'])
            self.canvas.draw_idle()
        debugger.print('on title change ', self.settings['Plot title'])

    def greyed_out(self):
        """Handle items that should be greyed out if they are not needed"""
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
