#
# mainTab
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'castep'
tab.settings['Output file name'] = 'phonon.castep'
#
# SettingsTab
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] = False
tab.settings['Neutral Born charges'] = False
tab.settings['Sigma value'] = 10
tab.settings['Mass definition'] = 'average'
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
tab.settings['Legend'] = 'Maxwell-Garnett sphere'
# Add new scenarios
methods = [ 'Maxwell-Garnett', 'Bruggeman']
shapes = ['Sphere', 'Needle']
vfs    = [ 0.1, 0.2 ]
hkl = [0,0,1]
for method in methods:
    for vf in vfs:
        for shape in shapes:
            self.notebook.addScenario()
            tab = self.notebook.scenarios[-1]
            tab.settings['Particle shape'] = shape
            tab.settings['Effective medium method'] = method
            tab.settings['Volume fraction'] = vf
            tab.settings['Legend'] = method+' vf='+str(vf)+' '+shape
            if shape != 'Sphere':
                tab.settings['Legend'] = method+' vf='+str(vf)+' '+shape+' '+str(hkl)
tab.settings['Effective medium method'] = 'Averaged Permittivity'
tab.settings['Particle shape'] = 'Sphere'
for vf in vfs:
    self.notebook.addScenario()
    tab = self.notebook.scenarios[-1]
    tab.settings['Volume fraction'] = vf
    tab.settings['Legend'] = 'AP vf='+str(vf)
self.notebook.deleteScenario(0)
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] = 300
tab.settings['Maximum frequency'] = 800
tab.settings['Frequency increment'] = 0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Plot title'] = 'Castep MgO'
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
