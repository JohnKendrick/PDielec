#
# mainTab
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'crystal'
tab.settings['Output file name'] = 'qua_hf_2d_f_ir-int.out'
tab.settings['Hessian symmetrisation'] = 'crystal'
#
# SettingsTab
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] = False
tab.settings['Neutral Born charges'] = False
tab.settings['Sigma value'] = 5
tab.settings['Mass definition'] = 'program'
tab.settings['Optical permittivity'] = [ [2.018, 0, 0], [0, 2.018, 0], [ 0,0,2.071] ]
#
# 0th Scenario tabs
#
tab = self.notebook.scenarios[0]
tab.settings['Matrix'] = 'ptfe'
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['Volume fraction'] = 0.1
tab.settings['Ellipsoid a/b'] = 0.5
tab.settings['Unique direction - h'] = 0
tab.settings['Unique direction - k'] = 0
tab.settings['Unique direction - l'] = 1
tab.settings['Effective medium method'] = 'Averaged Permittivity'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Legend'] = 'Averaged permittivity'
# Add new scenarios
methods = [ 'Maxwell-Garnett']
shapes = ['Needle', 'Ellipsoid', 'Plate']
hkls   = [[0,0,1],  [0,0,1]    , [1,0,0,]]
for method in methods:
    for shape,hkl in zip(shapes,hkls):
        self.notebook.addScenario()
        tab = self.notebook.scenarios[-1]
        tab.settings['Particle shape'] = shape
        tab.settings['Effective medium method'] = method
        tab.settings['Unique direction - h'] = hkl[0]
        tab.settings['Unique direction - k'] = hkl[1]
        tab.settings['Unique direction - l'] = hkl[2]
        tab.settings['Legend'] = method + ' ' + shape + ' ' +str(hkl)
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] = 0
tab.settings['Maximum frequency'] = 300
tab.settings['Frequency increment'] = 0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Plot title'] = 'Crystal14 - Quartz'
#
# Analysis Tab
#
tab = self.notebook.analysisTab
tab.settings['Minimum frequency'] = -1
tab.settings['Maximum frequency'] = 300
tab.settings['title'] = 'Analysis'
tab.settings['Covalent radius scaling'] = 1.1
tab.settings['Bonding tolerance'] = 0.1
tab.settings['Bar width'] = 0.5
#
