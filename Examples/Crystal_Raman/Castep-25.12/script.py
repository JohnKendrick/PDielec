#
# Handle the special case of the first scenario
#
self.notebook.switch_scenario(0,scenarioType="Crystal Raman")
#
#
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'Castep'
tab.settings['Output file name'] = 'raman.castep'
tab.settings['Excel file name'] = 'results.xlsx'
tab.settings['Script file name'] = 'script.py'
tab.settings['QM program'] = ''
tab.settings['Compatibility mode'] = 'Linux'
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  True
tab.settings['Neutral Born charges'] =  False
tab.settings['Sigma value'] =  5
tab.settings['Mass definition'] = 'average'
tab.settings['Optical permittivity edited'] =  False
tab.settings['Spectroscopy type'] = 'Crystal Raman'
tab.sigmas_cm1 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
#
# Scenario 1: Coherent depth integration - thin film 1 um, p-in unpolarised
# Tests: per-layer Coherent mode, Coherent amplitude depth, backscattering
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = 'Coherent depth - thin film 1um p-in'
tab.settings['Scenario type'] = 'Crystal Raman'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Scattering matrix'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer', 'air']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
tab.settings['Layer azimuthals'] =  [0, 0, 0]
tab.settings['Layer thicknesses'] =  [1, 1.0, 1]
tab.settings['Layer thickness units'] =  ['um', 'um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True, False]
tab.settings['Laser wavelength nm'] =  532.0
tab.settings['Incident polarisation'] = 'p'
tab.settings['Detected polarisation'] = 'unpolarised'
tab.settings['Temperature K'] =  298.0
tab.settings['GL point density'] =  20.0
tab.settings['Collection side'] = 'superstrate'
tab.settings['Collection angle'] =  -1.0
tab.settings['Depth coherence'] = 'Coherent amplitude'
tab.settings['Layer combination'] = 'Incoherent intensities'
tab.settings['Approximate ES'] =  False
tab.settings['Coalesce equivalent layers'] =  True
tab.settings['Raman electro-optic term'] =  True
tab.settings['Layer NAC mode'] = 'none'
tab.settings['Modal pair combination'] = 'Group q channels'
tab.settings['Porto notation'] = ''
#
# Scenario 2: Incoherent depth integration - thick slab 100 um, Incoherent(non-reflective) layer
# Tests: Incoherent (non-reflective) layer mode, Incoherent intensity depth
#
self.notebook.add_scenario(scenarioType="Crystal Raman")
tab = self.notebook.scenarios[1]
tab.settings['Legend'] = 'Incoherent depth - thick slab 100um'
tab.settings['Scenario type'] = 'Crystal Raman'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Scattering matrix'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer', 'air']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
tab.settings['Layer azimuthals'] =  [0, 0, 0]
tab.settings['Layer thicknesses'] =  [1, 100.0, 1]
tab.settings['Layer thickness units'] =  ['um', 'um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Incoherent (non-reflective)', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True, False]
tab.settings['Laser wavelength nm'] =  532.0
tab.settings['Incident polarisation'] = 'p'
tab.settings['Detected polarisation'] = 'unpolarised'
tab.settings['Temperature K'] =  298.0
tab.settings['GL point density'] =  20.0
tab.settings['Collection side'] = 'superstrate'
tab.settings['Collection angle'] =  -1.0
tab.settings['Depth coherence'] = 'Incoherent intensity'
tab.settings['Layer combination'] = 'Incoherent intensities'
tab.settings['Approximate ES'] =  False
tab.settings['Coalesce equivalent layers'] =  True
tab.settings['Raman electro-optic term'] =  True
tab.settings['Layer NAC mode'] = 'none'
tab.settings['Modal pair combination'] = 'Group q channels'
tab.settings['Porto notation'] = ''
#
# Scenario 3: s-polarisation in, s-polarisation out - thin film 1 um
# Tests: different polarisation channel selection
#
self.notebook.add_scenario(scenarioType="Crystal Raman")
tab = self.notebook.scenarios[2]
tab.settings['Legend'] = 's-in s-out polarisation 1um'
tab.settings['Scenario type'] = 'Crystal Raman'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Scattering matrix'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer', 'air']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
tab.settings['Layer azimuthals'] =  [0, 0, 0]
tab.settings['Layer thicknesses'] =  [1, 1.0, 1]
tab.settings['Layer thickness units'] =  ['um', 'um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True, False]
tab.settings['Laser wavelength nm'] =  532.0
tab.settings['Incident polarisation'] = 's'
tab.settings['Detected polarisation'] = 's'
tab.settings['Temperature K'] =  298.0
tab.settings['GL point density'] =  20.0
tab.settings['Collection side'] = 'superstrate'
tab.settings['Collection angle'] =  -1.0
tab.settings['Depth coherence'] = 'Coherent amplitude'
tab.settings['Layer combination'] = 'Incoherent intensities'
tab.settings['Approximate ES'] =  False
tab.settings['Coalesce equivalent layers'] =  True
tab.settings['Raman electro-optic term'] =  True
tab.settings['Layer NAC mode'] = 'none'
tab.settings['Modal pair combination'] = 'Group q channels'
tab.settings['Porto notation'] = ''
#
# Scenario 4: Forward scattering (substrate collection side) - thin film 1 um
# Tests: substrate collection (transmission geometry)
#
self.notebook.add_scenario(scenarioType="Crystal Raman")
tab = self.notebook.scenarios[3]
tab.settings['Legend'] = 'Forward scattering - substrate collection 1um'
tab.settings['Scenario type'] = 'Crystal Raman'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Scattering matrix'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer', 'air']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
tab.settings['Layer azimuthals'] =  [0, 0, 0]
tab.settings['Layer thicknesses'] =  [1, 1.0, 1]
tab.settings['Layer thickness units'] =  ['um', 'um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True, False]
tab.settings['Laser wavelength nm'] =  532.0
tab.settings['Incident polarisation'] = 'p'
tab.settings['Detected polarisation'] = 'unpolarised'
tab.settings['Temperature K'] =  298.0
tab.settings['GL point density'] =  20.0
tab.settings['Collection side'] = 'substrate'
tab.settings['Collection angle'] =  -1.0
tab.settings['Depth coherence'] = 'Coherent amplitude'
tab.settings['Layer combination'] = 'Incoherent intensities'
tab.settings['Approximate ES'] =  False
tab.settings['Coalesce equivalent layers'] =  True
tab.settings['Raman electro-optic term'] =  True
tab.settings['Layer NAC mode'] = 'none'
tab.settings['Modal pair combination'] = 'Group q channels'
tab.settings['Porto notation'] = ''
#
# Scenario 5: Oblique incidence 45 degrees - thin film 1 um
# Tests: non-zero angle of incidence
#
self.notebook.add_scenario(scenarioType="Crystal Raman")
tab = self.notebook.scenarios[4]
tab.settings['Legend'] = 'Oblique incidence 45deg p-in 1um'
tab.settings['Scenario type'] = 'Crystal Raman'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  45.0
tab.settings['Mode'] = 'Scattering matrix'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer', 'air']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
tab.settings['Layer azimuthals'] =  [0, 0, 0]
tab.settings['Layer thicknesses'] =  [1, 1.0, 1]
tab.settings['Layer thickness units'] =  ['um', 'um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True, False]
tab.settings['Laser wavelength nm'] =  532.0
tab.settings['Incident polarisation'] = 'p'
tab.settings['Detected polarisation'] = 'unpolarised'
tab.settings['Temperature K'] =  298.0
tab.settings['GL point density'] =  20.0
tab.settings['Collection side'] = 'superstrate'
tab.settings['Collection angle'] =  -1.0
tab.settings['Depth coherence'] = 'Coherent amplitude'
tab.settings['Layer combination'] = 'Incoherent intensities'
tab.settings['Approximate ES'] =  False
tab.settings['Coalesce equivalent layers'] =  True
tab.settings['Raman electro-optic term'] =  True
tab.settings['Layer NAC mode'] = 'none'
tab.settings['Modal pair combination'] = 'Group q channels'
tab.settings['Porto notation'] = ''
#
# Scenario 6: NAC geometry correction - thin film 1 um
# Tests: Layer NAC mode = 'geometry' (Snell's law phonon q)
#
self.notebook.add_scenario(scenarioType="Crystal Raman")
tab = self.notebook.scenarios[5]
tab.settings['Legend'] = 'NAC geometry correction 1um'
tab.settings['Scenario type'] = 'Crystal Raman'
tab.settings['Materials database'] = '../../../PDielec/MaterialsDataBase.xlsx'
tab.settings['Global azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  0.0
tab.settings['Mode'] = 'Scattering matrix'
tab.settings['Frequency units'] = 'wavenumber'
tab.settings['Partially incoherent samples'] =  20
tab.settings['Percentage partial incoherence'] =  0
tab.settings['Filter kernel size'] =  1
tab.settings['Filter polynomial size'] =  3
tab.settings['Layer material names'] =  ['air', 'Dielectric layer', 'air']
tab.settings['Layer hkls'] =  [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
tab.settings['Layer azimuthals'] =  [0, 0, 0]
tab.settings['Layer thicknesses'] =  [1, 1.0, 1]
tab.settings['Layer thickness units'] =  ['um', 'um', 'um']
tab.settings['Layer incoherent options'] =  ['Coherent', 'Coherent', 'Coherent']
tab.settings['Layer dielectric flags'] =  [False, True, False]
tab.settings['Laser wavelength nm'] =  532.0
tab.settings['Incident polarisation'] = 'p'
tab.settings['Detected polarisation'] = 'unpolarised'
tab.settings['Temperature K'] =  298.0
tab.settings['GL point density'] =  20.0
tab.settings['Collection side'] = 'superstrate'
tab.settings['Collection angle'] =  -1.0
tab.settings['Depth coherence'] = 'Coherent amplitude'
tab.settings['Layer combination'] = 'Incoherent intensities'
tab.settings['Approximate ES'] =  False
tab.settings['Coalesce equivalent layers'] =  True
tab.settings['Raman electro-optic term'] =  True
tab.settings['Layer NAC mode'] = 'geometry'
tab.settings['Modal pair combination'] = 'Group q channels'
tab.settings['Porto notation'] = ''
#
#
tab = self.notebook.analysisTab
tab.settings['Radii'] =  [0.66, 1.22]
tab.settings['Minimum frequency'] =  -1
tab.settings['Maximum frequency'] =  700
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
tab.settings['Maximum frequency'] =  700.0
tab.settings['Frequency increment'] =  0.5
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Number of atoms'] =  1
tab.settings['Plot type'] = 'Crystal Raman'
tab.settings['Frequency unit'] = 'wavenumber'
