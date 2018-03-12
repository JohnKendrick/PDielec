#
# SettingsTab
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] = True
tab.settings['Neutral Born charges'] = True
tab.settings['Sigma value'] = 5
tab.settings['Mass definition'] = 'program'
#
# 0th Scenario tabs
#
tab = self.notebook.scenarios[0]
tab.settings['Matrix'] = 'ptfe'
tab.settings['Matrix permittivity'] = 3.0
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['Volume fraction'] = 0.1
tab.settings['Ellipsoid a/b'] = 0.5
tab.settings['Unique direction - h'] = 0
tab.settings['Unique direction - k'] = 0
tab.settings['Unique direction - l'] = 1
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
# Add new scenarios
#methods = [ 'Maxwell-Garnett', 'Bruggeman','Averaged permittivity' ]
#shapes = ['Sphere', 'Needle']
#for method in methods:
#    for shape in shapes:
#        self.notebook.addScenario()
#        tab = self.notebook.scenarios[-1]
#        tab.settings['Particle shape'] = shape
#        tab.settings['Effective medium method'] = method
#self.notebook.deleteScenario(0)
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] = 300
tab.settings['Maximum frequency'] = 2800
tab.settings['Frequency increment'] = 0.2
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
