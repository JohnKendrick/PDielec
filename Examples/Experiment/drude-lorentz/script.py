#
# Change the first scenario tab to single crystal
#
tab = self.notebook
tab.switchScenario(0,scenarioType='Single crystal')
#
tab = self.notebook.scenarios[0]
tab.settings['Unique direction - h'] =  0
tab.settings['Unique direction - k'] =  0
tab.settings['Unique direction - l'] =  1
tab.settings['Angle of incidence'] =  10.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Superstrate depth'] =  999.0
tab.settings['Substrate depth'] =  999.0
tab.settings['Film thickness'] =  10000.0
tab.settings['Mode'] = 'Thick slab'
tab.settings['Legend'] = 'Angle 10'
#
#
for degs in [ 20, 30, 40, 50, 60, 70, 80]:
    self.notebook.addScenario(scenarioType='Single crystal')
    tab = self.notebook.scenarios[-1]
    tab.settings['Angle of incidence'] =  degs
    tab.settings['Legend'] = 'Angle {}'.format(degs)
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  100
tab.settings['Maximum frequency'] =  800
tab.settings['Frequency increment'] =  0.2
tab.settings['Plot title'] = 'Test of Drude-Lorentz Dielectric for MgO'
tab.settings['Plot type'] = 'Crystal Reflectance (P polarisation)'
#
#
