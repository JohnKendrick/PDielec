#
# Change the first scenario tab to single crystal
#
tab = self.notebook
tab.switchScenario(0,scenarioType='Single crystal')
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'experiment'
tab.settings['Output file name'] = 'constant.exp'
tab.settings['QM program'] = ''
#
tab = self.notebook.scenarios[0]
tab.settings['Legend'] = 'Thick slab'
tab.settings['Scenario type'] = 'Single crystal'
tab.settings['Azimuthal angle'] =  0.0
tab.settings['Angle of incidence'] =  45.0
tab.settings['Superstrate dielectric'] =  1.0
tab.settings['Substrate dielectric'] =  1.0
tab.settings['Mode'] = 'Thick slab'
tab.refresh(force=True)
#
#
self.notebook.addScenario(scenarioType='Single crystal')
tab = self.notebook.scenarios[1]
tab.settings['Legend'] = 'Film 10um'
tab.settings['Mode'] = 'Coherent thin film'
tab.settings['Film thickness'] =  10000.0
tab.refresh(force=True)
#
#
self.notebook.addScenario(scenarioType='Single crystal')
tab = self.notebook.scenarios[2]
tab.settings['Legend'] = 'Film 20um'
tab.settings['Mode'] = 'Coherent thin film'
tab.settings['Film thickness'] =  20000.0
tab.refresh(force=True)
#
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =  100
tab.settings['Maximum frequency'] =  400
tab.settings['Frequency increment'] =  0.2
tab.settings['Plot title'] = 'Test of Constant Dielectric'
tab.settings['Plot type'] = 'Crystal Reflectance (P polarisation)'
tab.refresh(force=True)
#
#
