#
# MainTab
#
tab = self.notebook.mainTab
tab.settings['Program'] = 'abinit'
tab.settings['QM program'] = ''
tab.settings['Output file name'] = 'temp.out'
tab.settings['Excel filename'] = 'temp.xlsx'
#
# SettingsTab
#
tab = self.notebook.settingsTab
tab.settings['Eckart flag'] = True
tab.settings['Neutral Born charges'] = False
tab.settings['Hessian symmetrisation'] = 'symm'              # crystal/symm
tab.settings['Sigma value'] = 5
tab.settings['Mass definition'] = 'average'                  # average/program/isotope
tab.settings['Optical permittivity'] = [[1,0,0],[0,1,0],[0,0,1]]
#tab.setElementMass['H',2.0]
#
# Scenario tabs
#
tab = self.notebook.scenarios[0]
tab.settings['Matrix'] = 'ptfe'                              # ptfe/air/kbr/nujol/ldpe
tab.settings['Matrix density'] = 1.0
tab.settings['Matrix permittivity'] = 1.0
tab.settings['Mass or volume fraction'] = 'volume'           # mass or volume
tab.settings['Volume fraction'] = 0.1                        # by default volume fraction has priority
tab.settings['Mass fraction'] = 0.1
tab.settings['Particle size(mu)'] = 0.0001
tab.settings['Particle size distribution sigma(mu)'] = 0.0
tab.settings['Ellipsoid a/b'] = 1.0
tab.settings['Unique direction - h'] = 0
tab.settings['Unique direction - k'] = 0
tab.settings['Unique direction - l'] = 1
tab.settings['Effective medium method'] = 'Maxwell-Garnett' # Averaged permittivity/Bruggeman/Mie/Maxwell-Garnett
tab.settings['Particle shape'] = 'Sphere'                   # Sphere/Plate/Ellipsoid/Needle
# Add new scenarios
self.notebook.addScenario(1)
tab = self.notebook.scenarios[1]
tab.settings['Matrix permittivity'] = 2.0
#
# Plotting Tab
#
tab = self.notebook.plottingTab
tab.settings['Minimum frequency'] = 0
tab.settings['Maximum frequency'] = 400
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
# GUI Tab
#
tab = self.notebook.viewerTab
tab.settings['Bond colour'] = [ 0.3, 0.3, 0.3, 1.0 ]
tab.settings['Bond radius'] = 0.1
tab.settings['Number of phase steps'] = 41
tab.settings['Cell colour'] = [ 1.0, 0.0, 0.0, 1.0 ]
tab.settings['Cell radius'] = 0.1
tab.settings['Background colour'] = [ 0.0, 0.0, 0.0, 1.0 ]
tab.light_switches[0] = False
tab.setColour('Na', [255, 255, 0, 255])
#
