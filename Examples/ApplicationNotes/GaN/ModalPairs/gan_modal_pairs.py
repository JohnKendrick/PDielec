"""PDGui automation for the GaN near-forward modal-pairs comparison.

The geometry follows Mina, Alhaddad, and Pagès, Appl. Phys. Lett. 128,
092101 (2026): propagation through an m-plane edge, c in the sample plane and
perpendicular to the incidence/detection plane, and fixed forward collection.
PDGui s and p then represent extraordinary (e, electric field parallel to c)
and ordinary (o) optical channels, respectively.
"""

angles_deg = [0.0, 2.0]

nac_modes = [
    ("geometry", "geometry"),
    ("dominant", "dominant_mode"),
    ("modal_pairs", "modal_pairs"),
]

polarisation_channels = [
    ("e-e", "s", "s"),
    ("o-o", "p", "p"),
    ("e-o", "s", "p"),
    ("o-e", "p", "s"),
]

self.notebook.switch_scenario(0, scenarioType="Crystal Raman")

tab = self.notebook.mainTab
tab.settings["Program"] = "Crystal"
tab.settings["Output file name"] = "../Raman/raman.log"
tab.settings["Excel file name"] = "../ModalPairs/Data/gan_modal_pairs.xlsx"
tab.settings["Script file name"] = "gan_modal_pairs.py"
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
# Acoustic modes 1-3 and silent B1 modes 6 and 12 are excluded.
tab.modes_selected = [
    False, False, False, True, True, False,
    True, True, True, True, True, False,
]


def configure_scenario(
    tab, legend, angle_deg, nac_mode, incident_pol, detected_pol
):
    """Configure one m-plane, near-forward GaN Raman scenario."""
    tab.settings["Legend"] = legend
    tab.settings["Scenario type"] = "Crystal Raman"
    # Air / 10 um GaN / air models transmission through an edge face.  The
    # (100) m-plane normal is the propagation direction.  At azimuth 0 deg,
    # c is perpendicular to the incidence/detection plane, hence s=e and p=o.
    tab.settings["Layer material names"] = [
        "air", "Dielectric layer", "air"
    ]
    tab.settings["Layer hkls"] = [[0, 0, 0], [1, 0, 0], [0, 0, 0]]
    tab.settings["Layer azimuthals"] = [0.0, 0.0, 0.0]
    tab.settings["Layer thicknesses"] = [1.0, 10.0, 1.0]
    tab.settings["Layer thickness units"] = ["um", "um", "um"]
    tab.settings["Layer incoherent options"] = [
        "Coherent", "Coherent", "Coherent"
    ]
    tab.settings["Layer dielectric flags"] = [False, True, False]
    tab.settings["Global azimuthal angle"] = 0.0

    # The experiment varies incidence at the rear while detection remains
    # normal to the front face.  PDGui angles are external angles; the paper's
    # quoted theta values are internal scattering angles.
    tab.settings["Angle of incidence"] = float(angle_deg)
    tab.settings["Collection side"] = "substrate"
    tab.settings["Collection angle"] = 0.0
    tab.settings["Mode"] = "Scattering matrix"
    tab.settings["Frequency units"] = "wavenumber"
    tab.settings["Laser wavelength nm"] = 488.0
    tab.settings["Temperature K"] = 298.0
    tab.settings["Incident polarisation"] = incident_pol
    tab.settings["Detected polarisation"] = detected_pol

    # A bulk-like, phase-matched comparison suppresses arbitrary interference
    # across the 10 um sampling depth while retaining local optical fields.
    tab.settings["GL point density"] = 10.0
    tab.settings["Depth coherence"] = "Incoherent intensity"
    tab.settings["Layer combination"] = "Incoherent intensities"
    tab.settings["Approximate ES"] = False
    tab.settings["Coalesce equivalent layers"] = True

    tab.settings["Layer NAC mode"] = nac_mode
    tab.settings["Modal pair combination"] = "Group q channels"
    tab.settings["Modal pair final-state model"] = "Bulk phase matched"
    tab.settings["Modal pair include zero q"] = None
    tab.settings["Modal pair q-angle tolerance"] = 90.0
    tab.settings["Raman electro-optic term"] = False
    tab.settings["Azimuthal sweep points"] = 1


scenarios = []
for angle_deg in angles_deg:
    for nac_label, nac_mode in nac_modes:
        for channel, incident_pol, detected_pol in polarisation_channels:
            legend = f"theta_ext={angle_deg:g} {nac_label} {channel}"
            scenarios.append(
                (legend, angle_deg, nac_mode, incident_pol, detected_pol)
            )

for index, scenario in enumerate(scenarios):
    if index > 0:
        self.notebook.add_scenario(scenarioType="Crystal Raman")
    configure_scenario(self.notebook.scenarios[index], *scenario)

tab = self.notebook.analysisTab
tab.settings["Minimum frequency"] = -1
tab.settings["Maximum frequency"] = 850
tab.settings["title"] = "GaN modal-pairs validation"
tab.settings["Covalent radius scaling"] = 1.1
tab.settings["Bonding tolerance"] = 0.1
tab.settings["Bar width"] = 0.5

tab = self.notebook.viewerTab
tab.settings["Selected mode"] = 8
tab.settings["hkl"] = (1, 0, 0)
tab.settings["uvw"] = (0, 0, 1)

tab = self.notebook.plottingTab
tab.settings["Minimum frequency"] = 400.0
tab.settings["Maximum frequency"] = 800.0
tab.settings["Frequency increment"] = 0.5
tab.settings["Plot type"] = "Crystal Raman"
tab.settings["Frequency unit"] = "wavenumber"
tab.settings["Spectrum renormalisation"] = "none"
