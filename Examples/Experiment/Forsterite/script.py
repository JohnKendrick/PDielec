#
# Change the first scenario tab to single crystal
#
tab = self.notebook
tab.switchScenario(0,scenarioType='Single crystal')
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'experiment'
tab.settings['Output file name'] = 'forsterite.exp'
tab.settings['Excel file name'] = ''
tab.settings['Script file name'] = 'temp.py'
tab.settings['QM program'] = ''
tab.settings['Hessian symmetrisation'] = 'symm'
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  True
tab.settings['Neutral Born charges'] =  False
tab.settings['Sigma value'] =  5
tab.settings['Mass definition'] = 'average'
tab.settings['Optical permittivity edited'] =  False
tab.sigmas_cm1 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
#
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = 'a-axis (001) phi=90'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Azimuthal angle'] =  90.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  100.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
#
#
self.notebook.addScenario()
tab = self.notebook.scenarios[1]
tab.settings['Legend'] = 'a-axis (010) phi=90'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  1
tab.settings['Unique direction - l'] =  0
tab.settings['Azimuthal angle'] =  90.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  100.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
#
#
self.notebook.addScenario()
tab = self.notebook.scenarios[2]
tab.settings['Legend'] = 'b-axis (001) phi=0'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  100.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
#
#
self.notebook.addScenario()
tab = self.notebook.scenarios[3]
tab.settings['Legend'] = 'b-axis (100) phi=90'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Unique direction - h'] =  1
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Azimuthal angle'] =  90.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  100.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
#
#
self.notebook.addScenario()
tab = self.notebook.scenarios[4]
tab.settings['Legend'] = 'c-axis (100) phi=0'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Unique direction - h'] =  1
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  100.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
#
#
self.notebook.addScenario()
tab = self.notebook.scenarios[5]
tab.settings['Legend'] = 'c-axis (010) phi=0'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  1
tab.settings['Unique direction - l'] =  0
tab.settings['Azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  100.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
#
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  350
tab.settings['Maximum frequency'] =  600
tab.settings['Frequency increment'] =  0.1
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Number of atoms'] =  1
tab.settings['Plot type'] = 'Crystal Reflectance (S polarisation)'
tab.settings['Plot title'] = 'Plot Title'
tab.settings['concentration'] =  5.751258957038493
#
#
tab = self.notebook.analysisTab
tab.settings['Minimum frequency'] =  -1
tab.settings['Maximum frequency'] =  400
tab.settings['title'] = 'Analysis'
tab.settings['Covalent radius scaling'] =  1.1
tab.settings['Bonding tolerance'] =  0.1
tab.settings['Bar width'] =  0.5
#
#
tab = self.notebook.viewerTab
tab.settings['Atom scaling'] =  0.5
tab.settings['Maximum displacement'] =  1.0
tab.settings['Bond colour'] =  [80, 80, 80, 255]
tab.settings['Bond radius'] =  0.1
tab.settings['Cell colour'] =  [255, 0, 0, 255]
tab.settings['Cell radius'] =  0.1
tab.settings['Background colour'] =  [120, 120, 120, 255]
tab.settings['Arrow colour'] =  [0, 255, 0, 255]
tab.settings['Arrow radius'] =  0.07
tab.settings['Number of phase steps'] =  41
tab.settings['Super Cell'] =  [1, 1, 1]
#
#
tab = self.notebook.fitterTab
tab.settings['Excel file name'] = ''
tab.settings['Plot title'] = 'Experimental and Calculated Spectral Comparison'
tab.settings['Fitting type'] = 'Minimise x-correlation'
tab.settings['Number of iterations'] =  20
tab.settings['Frequency scaling factor'] =  1.0
tab.settings['Optimise frequency scaling'] =  False
tab.settings['Spectrum scaling'] =  False
tab.settings['Spectrum scaling factor'] =  1.0
tab.settings['Independent y-axes'] =  True
tab.settings['Spectral difference threshold'] =  0.05
tab.settings['HPFilter lambda'] =  7.0
tab.settings['Baseline removal'] =  False
tab.settings['Scenario index'] =  0
