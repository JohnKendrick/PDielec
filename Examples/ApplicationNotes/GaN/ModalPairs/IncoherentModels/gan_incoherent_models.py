"""Compare the two incoherent-depth modal-pairs final-state models."""

angles_deg = [0.0, 2.0]

model_configs = [
    ("bulk", "Bulk phase matched", 0.0),
    ("bulk", "Bulk phase matched", 90.0),
    ("local", "Local incoherent", 0.0),
    ("local", "Local incoherent", 90.0),
]

# The GaN c-axis is perpendicular to the incidence/detection plane, so s=e
# and p=o for propagation through the (100) m-plane edge.
polarisation_channels = [
    ("e-e", "s", "s"),
    ("o-o", "p", "p"),
    ("e-o", "s", "p"),
    ("o-e", "p", "s"),
]

self.notebook.switch_scenario(0, scenarioType="Crystal Raman")

tab = self.notebook.mainTab
tab.settings["Program"] = "Crystal"
tab.settings["Output file name"] = "../../Raman/raman.log"
tab.settings["Excel file name"] = (
    "Data/gan_incoherent_models.xlsx"
)
tab.settings["Script file name"] = "gan_incoherent_models.py"
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
    angle_deg,
    final_state_model,
    tolerance_deg,
    incident_pol,
    detected_pol,
):
    """Configure one incoherent-depth final-state comparison."""
    tab.settings["Legend"] = legend
    tab.settings["Scenario type"] = "Crystal Raman"
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

    tab.settings["Angle of incidence"] = float(angle_deg)
    tab.settings["Collection side"] = "substrate"
    tab.settings["Collection angle"] = 0.0
    tab.settings["Mode"] = "Scattering matrix"
    tab.settings["Frequency units"] = "wavenumber"
    tab.settings["Laser wavelength nm"] = 488.0
    tab.settings["Temperature K"] = 298.0
    tab.settings["Incident polarisation"] = incident_pol
    tab.settings["Detected polarisation"] = detected_pol

    tab.settings["GL point density"] = 10.0
    tab.settings["Depth coherence"] = "Incoherent intensity"
    tab.settings["Layer combination"] = "Incoherent intensities"
    tab.settings["Approximate ES"] = False
    tab.settings["Coalesce equivalent layers"] = True

    tab.settings["Layer NAC mode"] = "modal_pairs"
    tab.settings["Modal pair combination"] = "Group q channels"
    tab.settings["Modal pair final-state model"] = final_state_model
    tab.settings["Modal pair include zero q"] = None
    tab.settings["Modal pair q-angle tolerance"] = float(tolerance_deg)
    tab.settings["Raman electro-optic term"] = False
    tab.settings["Azimuthal sweep points"] = 1


scenarios = []
for angle_deg in angles_deg:
    for model_label, final_state_model, tolerance_deg in model_configs:
        for channel, incident_pol, detected_pol in polarisation_channels:
            legend = (
                f"theta_ext={angle_deg:g} model={model_label} "
                f"tol={tolerance_deg:g} channel={channel}"
            )
            scenarios.append(
                (
                    legend,
                    angle_deg,
                    final_state_model,
                    tolerance_deg,
                    incident_pol,
                    detected_pol,
                )
            )

for index, scenario in enumerate(scenarios):
    if index > 0:
        self.notebook.add_scenario(scenarioType="Crystal Raman")
    configure_scenario(self.notebook.scenarios[index], *scenario)

tab = self.notebook.analysisTab
tab.settings["Minimum frequency"] = -1
tab.settings["Maximum frequency"] = 850
tab.settings["title"] = "GaN incoherent modal-pairs models"
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
