#
# Handle the special case of the first scenario
#
self.notebook.switchScenario(0,scenarioType="Powder")
#
#
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'Phonopy'
tab.settings['Output file name'] = 'phonopy.yaml'
tab.settings['Excel file name'] = 'results.xlsx'
tab.settings['Script file name'] = 'script.py'
tab.settings['Compatibility mode'] = 'Linux'
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  True
tab.settings['Neutral Born charges'] =  False
tab.settings['Sigma value'] =  5
tab.settings['Mass definition'] = 'average'
tab.settings['Optical permittivity edited'] =  False
tab.sigmas_cm1 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
#
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = 'Scenario 1'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Materials database'] = '../../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  (2+0j)
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.16650079135886944
tab.settings['Volume fraction'] =  0.09999999999999999
tab.settings['Particle size(mu)'] =  0.0001
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  1.0
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
#
#
tab = self.notebook.analysisTab
tab.settings['Radii'] =  [1.21, 0.66]
tab.settings['Minimum frequency'] =  -1
tab.settings['Maximum frequency'] =  400
tab.settings['title'] = 'Analysis'
tab.settings['Covalent radius scaling'] =  1.1
tab.settings['Bonding tolerance'] =  0.1
tab.settings['Bar width'] =  0.5
#
#
tab = self.notebook.viewerTab
tab.settings['Selected mode'] =  4
tab.settings['Atom scaling'] =  0.5
tab.settings['Maximum displacement'] =  1.0
tab.settings['Bond colour'] =  [80, 80, 80, 255]
tab.settings['Bond radius'] =  0.1
tab.settings['Cell colour'] =  [255, 0, 0, 255]
tab.settings['Cell radius'] =  0.1
tab.settings['Text colour'] =  [255, 255, 0, 255]
tab.settings['Text size'] = '24'
tab.settings['Background colour'] =  [120, 120, 120, 255]
tab.settings['Arrow colour'] =  [0, 255, 0, 255]
tab.settings['Arrow radius'] =  0.07
tab.settings['Number of phase steps'] =  41
tab.settings['Super Cell'] =  [1, 1, 1]
tab.settings['Transform'] =  [['1', '0', '0'], ['0', '1', '0'], ['0', '0', '1']]
tab.settings['hkl'] =  (0, 0, 1)
tab.settings['uvw'] =  (1, 0, 0)
tab.settings['Element colours'] =  None
tab.settings['Element palette'] = 'Jmol'
tab.settings['Toggle states'] =  [True, True, True, True, True]
#
#
tab = self.notebook.fitterTab
tab.settings['Experimental file name'] = ''
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
#
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  1.0
tab.settings['Maximum frequency'] =  900.0
tab.settings['Frequency increment'] =  0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Number of atoms'] =  1
tab.settings['Plot type'] = 'Powder Molar Absorption'
tab.settings['Frequency unit'] = 'wavenumber'
tab.settings['concentration'] =  19.39594919292005
tab.settings['cell concentration'] =  19.39594919292005
