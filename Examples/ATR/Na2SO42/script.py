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
tab.settings['Legend'] = 'sfrac=0.0'
tab.settings['ATR theta'] = 45.0
tab.settings['ATR S polarisation fraction'] = 0.0
tab.settings['ATR material refractive index'] = 4.0
# Add new scenarios
sfracs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6 , 0.7, 0.8, 0.9, 1.0 ]
for sfrac in sfracs:
    self.notebook.addScenario()
    tab = self.notebook.scenarios[-1]
    tab.settings['ATR S polarisation fraction'] = sfrac
    tab.settings['Legend'] = 'sfrac='+str(sfrac)
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] = 0
tab.settings['Maximum frequency'] = 1200
tab.settings['Frequency increment'] = 0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Plot title'] = 'Vasp Na2(SO4)2 Calculation'
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
