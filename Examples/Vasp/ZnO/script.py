#
# SettingsTab
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] = False
tab.settings['Neutral Born charges'] = False
tab.settings['Sigma value'] = 5.0
tab.settings['Mass definition'] = 'average'
#
# 0th Scenario tabs
#
tab = self.notebook.scenarios[0]
tab.settings['Matrix'] = 'ptfe'
tab.settings['Mass or volume fraction'] = 'mass'
tab.settings['Mass fraction'] = 0.01
tab.settings['Ellipsoid a/b'] = 2.0
tab.settings['Unique direction - h'] = 0
tab.settings['Unique direction - k'] = 0
tab.settings['Unique direction - l'] = 0
tab.settings['Effective medium method'] = 'Averaged permittivity'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Particle size(mu)'] = 0.0001
tab.settings['Legend'] = 'Averaged permittivity'
# Add new scenarios
methods = ['Maxwell-Garnett','Bruggeman']
shapes  = ['Sphere','Plate','Needle']
hkls    = [[0,0,0], [0,0,1],[0,0,1]]
mfs     = [0.01]
sizes   = [0.0]
for method in methods:
    for shape,hkl in zip(shapes,hkls):
        for mf in mfs:
            for size in sizes:
                self.notebook.addScenario()
                tab = self.notebook.scenarios[-1]
                tab.settings['Mass fraction'] = mf
                tab.settings['Particle shape'] = shape
                tab.settings['Particle size(mu)'] = size
                tab.settings['Effective medium method'] = method
                tab.settings['Unique direction - h'] = hkl[0]
                tab.settings['Unique direction - k'] = hkl[1]
                tab.settings['Unique direction - l'] = hkl[2]
                #tab.settings['Legend'] = method + ' ' + shape + ' vf='+str(vf)+' size='+str(size)
                tab.settings['Legend'] = method + ' ' + shape + str(hkl)
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =    0.0
tab.settings['Maximum frequency'] = 1000.0
tab.settings['Frequency increment'] = 0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Plot title'] = 'ZnO - VASP'
#
# Analysis Tab
#
tab = self.notebook.analysisTab
tab.settings['Minimum frequency'] = -1
tab.settings['Maximum frequency'] = 800
tab.settings['title'] = 'Analysis'
tab.settings['Covalent radius scaling'] = 1.1
tab.settings['Bonding tolerance'] = 0.1
tab.settings['Bar width'] = 0.5
#
