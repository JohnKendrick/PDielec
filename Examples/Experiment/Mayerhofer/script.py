#
# Handle the special case of the first scenario
#
self.notebook.switchScenario(0,scenarioType="Single crystal")
#
#
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'experiment'
tab.settings['Output file name'] = 'mayerhofer.exp'
tab.settings['Excel file name'] = ''
tab.settings['Script file name'] = 'new.py'
tab.settings['QM program'] = ''
tab.settings['Compatibility mode'] = 'Linux'
tab.settings['Hessian symmetrisation'] = 'symm'
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  True
tab.settings['Neutral Born charges'] =  False
tab.settings['Sigma value'] =  5
tab.settings['Mass definition'] = 'average'
tab.settings['Optical permittivity edited'] =  False
tab.sigmas_cm1 = [5, 5, 5, 5, 5, 5]
#
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = 'a-axis along X'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1]]
tab.settings['Layer azimuthals'] =  [0, 0]
tab.settings['Layer thicknesses'] =  [1, 1.0]
tab.settings['Layer thickness units'] =  ['um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True]
tab.settings['Slice thickness'] =  0
tab.settings['Slice thickness unit'] = 'um'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Azimuthal angle'] =  0.0
tab.settings['Film thickness'] =  100.0
#
#
self.notebook.addScenario(scenarioType="Single crystal")
tab = self.notebook.scenarios[1]
tab.settings['Legend'] = 'b-axis along X'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1]]
tab.settings['Layer azimuthals'] =  [0, 90.0]
tab.settings['Layer thicknesses'] =  [1, 1.0]
tab.settings['Layer thickness units'] =  ['um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True]
tab.settings['Slice thickness'] =  0
tab.settings['Slice thickness unit'] = 'um'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Azimuthal angle'] =  90.0
tab.settings['Film thickness'] =  100.0
#
#
self.notebook.addScenario(scenarioType="Single crystal")
tab = self.notebook.scenarios[2]
tab.settings['Legend'] = 'c-axis along X'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 1, 0]]
tab.settings['Layer azimuthals'] =  [0, 90.0]
tab.settings['Layer thicknesses'] =  [1, 1]
tab.settings['Layer thickness units'] =  ['um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True]
tab.settings['Slice thickness'] =  0
tab.settings['Slice thickness unit'] = 'um'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  1
tab.settings['Unique direction - l'] =  0
tab.settings['Azimuthal angle'] =  90.0
tab.settings['Film thickness'] =  100.0
#
#
tab = self.notebook.analysisTab
tab.settings['Radii'] =  [1.11, 0.66]
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
tab.settings['Excel file name'] = ''
#
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  100
tab.settings['Maximum frequency'] =  1300
tab.settings['Frequency increment'] =  0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Number of atoms'] =  1
tab.settings['Plot type'] = 'Crystal Reflectance (P polarisation)'
tab.settings['Frequency unit'] = 'wavenumber'
tab.settings['concentration'] =  1660.538783162726
tab.settings['cell concentration'] =  1660.538783162726
