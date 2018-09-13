#
# SettingsTab
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] = False
tab.settings['Neutral Born charges'] = False
tab.settings['Sigma value'] = 5
tab.settings['Mass definition'] = 'average'
#
# 0th Scenario tabs
#
tab = self.notebook.scenarios[0]
tab.settings['Matrix'] = 'air'
tab.settings['Mass or volume fraction'] = 'volume'
tab.settings['Volume fraction'] = 0.5
tab.settings['Ellipsoid a/b'] = 0.5
tab.settings['Unique direction - h'] = 0
tab.settings['Unique direction - k'] = 0
tab.settings['Unique direction - l'] = 1
tab.settings['Effective medium method'] = 'Maxwell-Garnett'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Legend'] = 'theta=0'
tab.settings['ATR theta'] = 10.0
tab.settings['ATR S polarisation fraction'] = 0.5
tab.settings['ATR material refractive index'] = 4.0
# Add new scenarios
angles = [10, 20, 30, 40, 50, 60 , 70, 80 ]
for angle in angles:
    self.notebook.addScenario()
    tab = self.notebook.scenarios[-1]
    tab.settings['ATR theta'] = angle
    tab.settings['Legend'] = 'theta='+str(angle)
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] = 200
tab.settings['Maximum frequency'] = 500
tab.settings['Frequency increment'] = 0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Plot title'] = 'AbInit AlAs Calculation'
#
# Analysis Tab
#
tab = self.notebook.analysisTab
tab.settings['Minimum frequency'] = -1
tab.settings['Maximum frequency'] = 500
tab.settings['title'] = 'Analysis'
tab.settings['Covalent radius scaling'] = 1.1
tab.settings['Bonding tolerance'] = 0.1
tab.settings['Bar width'] = 0.5
#
