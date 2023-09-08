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
tab.settings['Excel file name'] = 'results.xlsx'
tab.settings['Script file name'] = 'script.py'
tab.settings['QM program'] = 'vasp'
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  False
tab.settings['Neutral Born charges'] =  True
tab.settings['Sigma value'] =  5.0
tab.settings['Mass definition'] = 'program'
tab.settings['Optical permittivity edited'] =  False
tab.sigmas_cm1 = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
#
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = 'Fused Silica - Maxwell Garnett (10%)'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'Fused Silica'
tab.settings['Matrix density'] =  2.201
tab.settings['Matrix permittivity'] =  (3.498916942719922+0j)
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.11615971773044519
tab.settings['Volume fraction'] =  0.10000000000000003
tab.settings['Particle size(mu)'] =  0.0001
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  0.5
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Materials database'] = '/home/john/Software/PDielec/PDielec/MaterialsDataBase.xlsx'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[1]
tab.settings['Legend'] = 'Manual (rho=2.2,eps=3+0.1j) - Maxwell Garnett (10%)'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'Define material manually'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  (3+0.1j)
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.11620638192008519
tab.settings['Volume fraction'] =  0.09999999999999996
tab.settings['Particle size(mu)'] =  0.0001
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  0.5
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Materials database'] = '/home/john/Software/PDielec/PDielec/MaterialsDataBase.xlsx'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[2]
tab.settings['Legend'] = 'Water (25C) - Maxwell Garnett (10%)'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'Water(25C)'
tab.settings['Matrix density'] =  1.0
tab.settings['Matrix permittivity'] =  (5.15619999999997+2.6243579999999938j)
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.2243666379592112
tab.settings['Volume fraction'] =  0.1
tab.settings['Particle size(mu)'] =  10.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  0.5
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Materials database'] = '/home/john/Software/PDielec/PDielec/MaterialsDataBase.xlsx'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[3]
tab.settings['Legend'] = 'Water (19C) - Maxwell Garnett (10%)'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'Water(19C)'
tab.settings['Matrix density'] =  1.0
tab.settings['Matrix permittivity'] =  (12.550267025473484+16.940581205830732j)
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.2243666379592112
tab.settings['Volume fraction'] =  0.1
tab.settings['Particle size(mu)'] =  10.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  0.5
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Materials database'] = '/home/john/Software/PDielec/PDielec/MaterialsDataBase.xlsx'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[4]
tab.settings['Legend'] = 'PVC - Maxwell Garnett (10%)'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'PVC'
tab.settings['Matrix density'] =  1.38
tab.settings['Matrix permittivity'] =  (2.459384803400675+0.01756442892715219j)
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.17329079239199838
tab.settings['Volume fraction'] =  0.10000000000000003
tab.settings['Particle size(mu)'] =  10.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  0.5
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Materials database'] = '/home/john/Software/PDielec/PDielec/MaterialsDataBase.xlsx'
#
#
self.notebook.addScenario(scenarioType="Powder")
tab = self.notebook.scenarios[5]
tab.settings['Legend'] = 'Test - Maxwell Garnett (10%)'
tab.settings['Scenario type'] = 'Powder'
tab.settings['Matrix'] = 'Test'
tab.settings['Matrix density'] =  2.0
tab.settings['Matrix permittivity'] =  (10+0.1j)
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.12635865193552115
tab.settings['Volume fraction'] =  0.10000000000000003
tab.settings['Particle size(mu)'] =  10.0
tab.settings['Particle size distribution sigma(mu)'] =  0.0
tab.settings['Ellipsoid a/b'] =  0.5
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['ATR material refractive index'] =  4.0
tab.settings['ATR theta'] =  45.0
tab.settings['ATR S polarisation fraction'] =  0.5
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Materials database'] = '/home/john/Software/PDielec/Examples/Vasp/Na2SO42_v/TestMaterialsDataBase.xlsx'
#
#
tab = self.notebook.analysisTab
tab.settings['Radii'] =  [1.66, 1.05, 0.66]
tab.settings['Minimum frequency'] =  -1
tab.settings['Maximum frequency'] =  800
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
tab.settings['Background colour'] =  [120, 120, 120, 255]
tab.settings['Arrow colour'] =  [0, 255, 0, 255]
tab.settings['Arrow radius'] =  0.07
tab.settings['Number of phase steps'] =  41
tab.settings['Super Cell'] =  [1, 1, 1]
tab.settings['Element colours'] =  None
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
tab.settings['Maximum frequency'] =  500.0
tab.settings['Frequency increment'] =  0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Number of atoms'] =  1
tab.settings['Plot type'] = 'Powder Absorption'
tab.settings['Plot title'] = 'VASP - Na2(SO4)2'
tab.settings['concentration'] =  10.933228754034278
tab.settings['cell concentration'] =  10.933228754034278
