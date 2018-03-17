#
# SettingsTab
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] = False
tab.settings['Neutral Born charges'] = False
tab.settings['Sigma value'] = 2
tab.settings['Mass definition'] = 'average'
#
# 0th Scenario tabs
#
tab = self.notebook.scenarios[0]
tab.settings['Matrix'] = 'ptfe'
tab.settings['Volume fraction'] = 0.1
tab.settings['Ellipsoid a/b'] = 0.5
tab.settings['Unique direction - h'] = 0
tab.settings['Unique direction - k'] = 0
tab.settings['Unique direction - l'] = 1
tab.settings['Effective medium method'] = 'Averaged Permittivity'
tab.settings['Particle shape'] = 'Sphere'
# Add new scenarios
method = 'Maxwell-Garnett'
shapes = ['Sphere', 'Plate', 'Plate','Ellipsoid']
hkls    = [[0,0,1],  [1,0,0], [0,0,1], [0,0,1]]
for shape,hkl in zip(shapes,hkls):
        self.notebook.addScenario()
        tab = self.notebook.scenarios[-1]
        tab.settings['Particle shape'] = shape
        tab.settings['Effective medium method'] = method
        tab.settings['Unique direction - h'] = hkl[0]
        tab.settings['Unique direction - k'] = hkl[1]
        tab.settings['Unique direction - l'] = hkl[2]
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] = 0
tab.settings['Maximum frequency'] = 4000
tab.settings['Frequency increment'] = 0.2
tab.molar_definitions = ['Unit cells','Atoms','Molecules']
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Plot title'] = 'Script Title'
#
# Analysis Tab
#
tab = self.notebook.analysisTab
tab.settings['Minimum frequency'] = -1
tab.settings['Maximum frequency'] = 400
tab.settings['title'] = 'Analysis'
tab.settings['Covalent radius scaling'] = 1.1
tab.settings['Bonding tolerance'] = 0.1
tab.settings['Bar width'] = 0.5
#
