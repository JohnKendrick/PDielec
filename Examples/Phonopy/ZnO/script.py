#
# mainTab
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'phonopy'
tab.settings['Phonopy QM program'] = 'vasp'
tab.settings['Output file name'] = 'OUTCAR.born'
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
tab.settings['Matrix'] = 'ptfe'
tab.settings['Mass or volume fraction'] = 'mass'
tab.settings['Mass fraction'] = 0.01
tab.settings['Ellipsoid a/b'] = 0.5
tab.settings['Unique direction - h'] = 0
tab.settings['Unique direction - k'] = 0
tab.settings['Unique direction - l'] = 0
tab.settings['Effective medium method'] = 'Averaged permittivity'
tab.settings['Particle shape'] = 'Sphere'
tab.settings['Particle size(mu)'] = 0.0
tab.settings['Particle size distribution sigma(mu)'] = 0.0
tab.settings['Legend'] = 'Average permittivity'
# Add new scenarios
methods = ['Maxwell-Garnett', 'Bruggeman']
shapes  = ['Sphere','Plate','Needle']
hkls    = [[0,0,0], [0,0,1],[0,0,1]]
vfs     = [0.1]
sizes   =  [0.0001]
sigmas   = [0.0]
for method in methods:
    for shape,hkl in zip(shapes,hkls):
        for vf in vfs:
            for size,sigma in zip(sizes,sigmas):
                self.notebook.addScenario()
                tab = self.notebook.scenarios[-1]
                tab.settings['Mass or volume fraction'] = 'volume'
                tab.settings['Volume fraction'] = vf
                tab.settings['Particle shape'] = shape
                tab.settings['Particle size(mu)'] = size
                tab.settings['Effective medium method'] = method
                tab.settings['Particle size distribution sigma(mu)'] = sigma
                tab.settings['Unique direction - h'] = hkl[0]
                tab.settings['Unique direction - k'] = hkl[1]
                tab.settings['Unique direction - l'] = hkl[2]
                #tab.settings['Legend'] = method + ' ' + shape + ' vf='+str(vf)+' size='+str(size)+' sigma=',str(sigma)
                tab.settings['Legend'] = method + ' ' + shape + ' hkl='+str(hkl)
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] =   0.0
tab.settings['Maximum frequency'] = 1000.0
tab.settings['Frequency increment'] = 0.2
tab.settings['Molar definition'] = 'Unit cells'
tab.settings['Plot title'] = 'Phonopy ZnO'
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
