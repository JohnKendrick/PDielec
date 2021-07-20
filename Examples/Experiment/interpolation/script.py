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
