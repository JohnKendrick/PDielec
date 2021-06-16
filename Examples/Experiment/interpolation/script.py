#
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'experiment'
tab.settings['Output file name'] = 'interpolation.exp'
tab.settings['Excel file name'] = 'results.xlsx'
tab.settings['QM program'] = ''
tab.settings['Hessian symmetrisation'] = 'symm'
tab.refresh(force=True)
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  True
tab.settings['Neutral Born charges'] =  False
tab.settings['Sigma value'] =  5
tab.settings['Mass definition'] = 'average'
tab.settings['Optical permittivity edited'] =  False
tab.sigmas_cm1 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
tab.refresh(force=True)
#
#
tab = self.notebook.scenarios[0]
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix density'] =  2.2
tab.settings['Matrix permittivity'] =  2.0
tab.settings['Bubble radius'] =  30.0
tab.settings['Bubble volume fraction'] =  0.0
tab.settings['Mass fraction'] =  0.11120635086212491
tab.settings['Volume fraction'] =  0.1
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
tab.settings['Legend'] = 'Scenario 1'
tab.refresh(force=True)
#
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  0
tab.settings['Maximum frequency'] =  400
tab.settings['Frequency increment'] =  0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Number of atoms'] =  1
tab.settings['Plot title'] = 'Plot Title'
tab.settings['concentration'] =  13.743950329674377
tab.refresh(force=True)
#
#
tab = self.notebook.singleCrystalTab
tab.settings['Minimum frequency'] =  100
tab.settings['Maximum frequency'] =  400
tab.settings['Frequency increment'] =  0.2
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  45.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  10.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Plot title'] = 'Single crystal Thick slab'
tab.settings['Use default plot title'] =  True
tab.settings['Frequency units'] = 'wavenumber'
tab.refresh(force=True)
#
#
tab = self.notebook.analysisTab
tab.settings['Minimum frequency'] =  -1
tab.settings['Maximum frequency'] =  400
tab.settings['title'] = 'Analysis'
tab.settings['Covalent radius scaling'] =  1.1
tab.settings['Bonding tolerance'] =  0.1
tab.settings['Bar width'] =  0.5
tab.refresh(force=True)
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
tab.refresh(force=True)
#
#
tab = self.notebook.fitterTab
tab.settings['Excel file name'] = ''
tab.settings['Plot title'] = 'Experimental and Calculated Spectral Comparison'
tab.settings['Plot type'] = 'Molar absorption'
tab.settings['Fitting type'] = 'Minimise x-correlation'
tab.settings['Number of iterations'] =  20
tab.settings['Frequency scaling factor'] =  1.0
tab.settings['Optimise frequency scaling'] =  False
tab.settings['Absorption scaling'] =  False
tab.settings['Absorption scaling factor'] =  1.0
tab.settings['Independent y-axes'] =  True
tab.settings['Spectral difference threshold'] =  0.05
tab.settings['HPFilter lambda'] =  7.0
tab.settings['Baseline removal'] =  False
tab.refresh(force=True)
