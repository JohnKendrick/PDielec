"""ZnO spheres in PTFE: compare the finite-particle electro-optic term off/on."""

self.notebook.switch_scenario(0, scenarioType="Powder Raman")

tab = self.notebook.mainTab
tab.settings["Program"] = "Finite_field"
tab.settings["Output file name"] = "finite_difference.json"
tab.settings["Excel file name"] = "results.xlsx"
tab.settings["Script file name"] = "script.py"

tab = self.notebook.settingsTab
tab.settings["Eckart flag"] = True
tab.settings["Neutral Born charges"] = False
tab.settings["Sigma value"] = 5.0
tab.settings["Mass definition"] = "average"
tab.settings["Optical permittivity edited"] = False
tab.settings["Spectroscopy type"] = "Powder Raman"
tab.sigmas_cm1 = [5.0] * 12

tab = self.notebook.scenarios[0]
tab.settings["Legend"] = "ZnO spheres in PTFE, EO off"
tab.settings["Materials database"] = "../../../PDielec/MaterialsDataBase.xlsx"
tab.settings["Matrix"] = "ptfe"
tab.settings["Matrix density"] = 2.2
tab.settings["Matrix permittivity"] = 2.0 + 0.0j
tab.settings["Volume fraction"] = 0.1
tab.settings["Mass or volume fraction"] = "volume"
tab.settings["Particle shape"] = "Sphere"
tab.settings["Effective medium method"] = "Maxwell-Garnett"
tab.settings["Raman laser frequency"] = 785
tab.settings["Raman laser polarisation"] = "VV"
tab.settings["Raman temperature"] = 298.0
tab.settings["Raman electro-optic term"] = False

self.notebook.add_scenario(copyFromIndex=0)
tab = self.notebook.scenarios[1]
tab.settings["Legend"] = "ZnO spheres in PTFE, EO on"
tab.settings["Raman electro-optic term"] = True

tab = self.notebook.plottingTab
tab.settings["Minimum frequency"] = 1.0
tab.settings["Maximum frequency"] = 700.0
tab.settings["Frequency increment"] = 1.0
tab.settings["Plot type"] = "Powder Raman"
tab.settings["Frequency unit"] = "wavenumber"
