#
# Handle the special case of the first scenario
#
self.notebook.switchScenario(0,scenarioType="Powder")
#
#
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'vasp'
tab.settings['Output file name'] = 'OUTCAR'
tab.settings['Excel file name'] = 'application_note.xlsx'
tab.settings['Script file name'] = 'application_note.py'
tab.settings['QM program'] = 'vasp'
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  False
tab.settings['Neutral Born charges'] =  False
tab.settings['Sigma value'] =  2.0
tab.settings['Mass definition'] = 'average'
tab.settings['Optical permittivity edited'] =  False
tab.sigmas_cm1 = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
#
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = 'Averaged permittivity'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  2.0
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.13646688033570456
tab.settings['Volume fraction'] =  0.09999999999999999
tab.settings['Particle size(mu)'] =  0.0001
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  2.0
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Averaged permittivity'
tab.settings['Particle shape'] = 'Sphere'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[1]
tab.settings['Legend'] = 'Ellipsoid [001] a/b=2'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  2.0
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.13646688033570456
tab.settings['Volume fraction'] =  0.09999999999999999
tab.settings['Particle size(mu)'] =  0.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  2.0
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Ellipsoid'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[2]
tab.settings['Legend'] = 'Needle [001]'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  2.0
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.13646688033570456
tab.settings['Volume fraction'] =  0.09999999999999999
tab.settings['Particle size(mu)'] =  0.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  2.0
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Needle'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[3]
tab.settings['Legend'] = 'Plate (100)'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  2.0
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.13646688033570456
tab.settings['Volume fraction'] =  0.09999999999999999
tab.settings['Particle size(mu)'] =  0.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  2.0
tab.settings['Unique direction - h'] =  1
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Plate'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[4]
tab.settings['Legend'] = 'Sphere'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  2.0
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.13646688033570456
tab.settings['Volume fraction'] =  0.09999999999999999
tab.settings['Particle size(mu)'] =  0.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  2.0
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
#
#
tab = self.notebook.analysisTab
tab.settings['Minimum frequency'] =  -1
tab.settings['Maximum frequency'] =  800
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
#
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  900.0
tab.settings['Maximum frequency'] =  1100.0
tab.settings['Frequency increment'] =  0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Number of atoms'] =  1
tab.settings['Plot type'] = 'Powder Molar Absorption'
tab.settings['Plot title'] = 'VASP - F-Apatite'
tab.settings['concentration'] =  3.102361108197526
tab.settings['cell concentration'] =  3.102361108197526