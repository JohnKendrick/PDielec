"""PDGui scenarios matching the optical geometries of Irmer et al. Figs. 9--11.

This is the q->0 NAC/control half of the validation.  It deliberately does
not claim to calculate finite-q polariton dispersion; ``polariton_response.py``
performs that complementary calculation.
"""

LASER_WAVELENGTH_NM = 514.5
SCREEN_DISTANCE_MM = 80.0
Y_VALUES_MM = [2.0, 3.5, 4.5, 6.0, 8.0, 10.0, 12.5, 15.0]
Z_VALUES_MM = [3.0, 4.0, 5.0, 7.0, 10.0]


self.notebook.switch_scenario(0, scenarioType="Crystal Raman")

tab = self.notebook.mainTab
tab.settings["Program"] = "Crystal"
tab.settings["Output file name"] = "../../Raman/raman.log"
tab.settings["Excel file name"] = "Data/irmer_pdgui.xlsx"
tab.settings["Script file name"] = "irmer_pdgui.py"
tab.settings["Compatibility mode"] = "Linux"

tab = self.notebook.settingsTab
tab.settings["Eckart flag"] = True
tab.settings["Neutral Born charges"] = False
tab.settings["Sigma value"] = 5
tab.settings["Mass definition"] = "average"
tab.settings["Optical permittivity edited"] = False
tab.settings["Spectroscopy type"] = "Crystal Raman"
tab.settings["Raman activity units"] = "polarizability"
tab.sigmas_cm1 = [5] * 12
tab.modes_selected = [
    False, False, False, True, True, False,
    True, True, True, True, True, False,
]


def configure_scenario(
    tab,
    legend,
    offset_mm,
    azimuth_deg,
    incident_pol,
    detected_pol,
):
    """Configure one forward-scattering aperture position."""
    from math import atan, degrees

    tab.settings["Legend"] = legend
    tab.settings["Scenario type"] = "Crystal Raman"
    tab.settings["Layer material names"] = ["air", "Dielectric layer", "air"]
    # The paper uses an m-plane entrance face: surface normal x, c axis z.
    tab.settings["Layer hkls"] = [[0, 0, 0], [1, 0, 0], [0, 0, 0]]
    tab.settings["Layer azimuthals"] = [0.0, 0.0, 0.0]
    tab.settings["Layer thicknesses"] = [1.0, 10.0, 1.0]
    tab.settings["Layer thickness units"] = ["um", "um", "um"]
    tab.settings["Layer incoherent options"] = ["Coherent", "Coherent", "Coherent"]
    tab.settings["Layer dielectric flags"] = [False, True, False]

    # azimuth=0: scattering plane xy, s=z and p=y.
    # azimuth=90: scattering plane xz, p=z and s=y.
    tab.settings["Global azimuthal angle"] = float(azimuth_deg)
    tab.settings["Angle of incidence"] = 0.0
    tab.settings["Collection side"] = "substrate"
    tab.settings["Collection angle"] = degrees(
        atan(float(offset_mm) / 80.0)
    )
    tab.settings["Mode"] = "Scattering matrix"
    tab.settings["Frequency units"] = "wavenumber"
    tab.settings["Laser wavelength nm"] = 514.5
    tab.settings["Temperature K"] = 298.0
    tab.settings["Incident polarisation"] = incident_pol
    tab.settings["Detected polarisation"] = detected_pol

    # The bulk-like setting avoids treating the arbitrary 10 um sampling depth
    # as a phase-stable thin film.  It is the same incoherent convention used
    # in the preceding GaN modal-pairs validation.
    tab.settings["GL point density"] = 10.0
    tab.settings["Depth coherence"] = "Incoherent intensity"
    tab.settings["Layer combination"] = "Incoherent intensities"
    tab.settings["Approximate ES"] = False
    tab.settings["Coalesce equivalent layers"] = True
    tab.settings["Layer NAC mode"] = "modal_pairs"
    tab.settings["Modal pair combination"] = "Group q channels"
    tab.settings["Modal pair final-state model"] = "Bulk phase matched"
    tab.settings["Modal pair include zero q"] = None
    tab.settings["Modal pair q-angle tolerance"] = 90.0
    tab.settings["Raman electro-optic term"] = False
    tab.settings["Azimuthal sweep points"] = 1


scenarios = []

# Figure 9: x(zy)x, scattering in xy.  At azimuth zero this is s -> p.
for y_mm in Y_VALUES_MM:
    scenarios.append((f"Fig9 Y={y_mm:g}mm x(zy)x", y_mm, 0.0, "s", "p"))

# Figure 10: x(yz)x, same ray geometry with reciprocal polarization p -> s.
for y_mm in Y_VALUES_MM:
    scenarios.append((f"Fig10 Y={y_mm:g}mm x(yz)x", y_mm, 0.0, "p", "s"))

# Figure 11: scattering in xz.  Rotating the plane makes p=z and s=y.
for z_mm in Z_VALUES_MM:
    scenarios.append((f"Fig11a Z={z_mm:g}mm x(zz)x", z_mm, 90.0, "p", "p"))
for z_mm in Z_VALUES_MM:
    scenarios.append((f"Fig11b Z={z_mm:g}mm x(yy)x", z_mm, 90.0, "s", "s"))

for index, scenario in enumerate(scenarios):
    if index > 0:
        self.notebook.add_scenario(scenarioType="Crystal Raman")
    configure_scenario(self.notebook.scenarios[index], *scenario)

tab = self.notebook.analysisTab
tab.settings["Minimum frequency"] = -1
tab.settings["Maximum frequency"] = 850
tab.settings["title"] = "GaN Irmer 2013 geometry verification"
tab.settings["Covalent radius scaling"] = 1.1
tab.settings["Bonding tolerance"] = 0.1
tab.settings["Bar width"] = 0.5

tab = self.notebook.viewerTab
tab.settings["Selected mode"] = 8
tab.settings["hkl"] = (1, 0, 0)
tab.settings["uvw"] = (0, 0, 1)

tab = self.notebook.plottingTab
tab.settings["Minimum frequency"] = 150.0
tab.settings["Maximum frequency"] = 780.0
tab.settings["Frequency increment"] = 0.5
tab.settings["Plot type"] = "Crystal Raman"
tab.settings["Frequency unit"] = "wavenumber"
tab.settings["Spectrum renormalisation"] = "none"
