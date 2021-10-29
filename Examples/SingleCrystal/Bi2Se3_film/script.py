#
# Handle the special case of the first scenario
#
self.notebook.switchScenario(0,scenarioType="Single crystal")
#
#
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'vasp'
tab.settings['Output file name'] = 'OUTCAR'
tab.settings['Excel file name'] = ''
tab.settings['QM program'] = 'vasp'
tab.settings['Hessian symmetrisation'] = 'symm'
#
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] =  True
tab.settings['Neutral Born charges'] =  False
tab.settings['Sigma value'] =  5
tab.settings['Mass definition'] = 'average'
tab.settings['Optical permittivity edited'] =  False
tab.sigmas_cm1 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
#
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = '(001) theta  0'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  60.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  100.0
tab.settings['Mode'] = 'Coherent thin film'
tab.settings['Frequency units'] = 'wavenumber'
for az  in [ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
    self.notebook.addScenario(scenarioType="Single crystal")
    tab = self.notebook.scenarios[-1]
    tab.settings['Legend'] = '(001) theta {}'.format(az)
    tab.settings['Azimuthal angle'] =  az
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  0
tab.settings['Maximum frequency'] =  250
tab.settings['Frequency increment'] =  0.2
tab.settings['Plot type'] = 'Crystal Reflectance (S polarisation)'
