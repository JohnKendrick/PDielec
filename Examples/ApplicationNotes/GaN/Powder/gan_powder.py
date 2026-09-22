"""PDGui session comparing bulk tensor averages with GaN spheres in air.

Run from this directory with pdgui -script gan_powder.py. Matrix=none is
the uncorrected TO reference; air enables particle and local-field corrections.
"""

self.notebook.switch_scenario(0, scenarioType="Powder Raman")
tab = self.notebook.mainTab
tab.settings["Program"] = "Crystal"
tab.settings["Output file name"] = "../Raman/raman.log"
tab.settings["Excel file name"] = "Data/gan_powder.xlsx"
tab.settings["Script file name"] = "gan_powder.py"

tab = self.notebook.settingsTab
tab.settings["Eckart flag"] = True
tab.settings["Neutral Born charges"] = False
tab.settings["Sigma value"] = 5
tab.settings["Mass definition"] = "average"
tab.settings["Optical permittivity edited"] = False
tab.settings["Spectroscopy type"] = "Powder Raman"
tab.settings["Raman activity units"] = "polarizability"
tab.sigmas_cm1 = [5] * 12
tab.modes_selected = [False, False, False, True, True, False,
                      True, True, True, True, True, False]

for index, (matrix, polarisation) in enumerate(
    [(matrix, pol) for matrix in ("none", "air") for pol in ("VV", "VH")]
):
    if index:
        self.notebook.add_scenario(scenarioType="Powder Raman")
    tab = self.notebook.scenarios[index]
    tab.settings["Scenario type"] = "Powder Raman"
    tab.settings["Legend"] = f"{matrix} {polarisation}"
    tab.settings["Matrix"] = matrix
    tab.settings["Mass or volume fraction"] = "volume"
    tab.settings["Volume fraction"] = 0.1
    tab.settings["Particle shape"] = "Sphere"
    tab.settings["Particle size(mu)"] = 0.0001
    tab.settings["Particle size distribution sigma(mu)"] = 0.0
    tab.settings["Bubble volume fraction"] = 0.0
    tab.settings["Raman laser wavelength"] = 488.0
    tab.settings["Raman laser polarisation"] = polarisation
    tab.settings["Raman temperature"] = 298.0
    tab.settings["Raman electro-optic term"] = False

tab = self.notebook.plottingTab
tab.settings["Minimum frequency"] = 100.0
tab.settings["Maximum frequency"] = 800.0
tab.settings["Frequency increment"] = 0.5
tab.settings["Plot type"] = "Powder Raman"
tab.settings["Frequency unit"] = "wavenumber"
tab.settings["Spectrum renormalisation"] = "none"
