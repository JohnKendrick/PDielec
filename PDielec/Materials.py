import numpy as np
from PDielec.DielectricFunction import DielectricFunction

class Material():
    def __init__(self, name, density=None, epsType=None, epsinf=0.0, volume=0, permittivity=None):
        # Manage a database of permittivities and densities
        self.permittivity = permittivity
        self.density = density
        self.name = name
        self.epsinf = epsinf
        self.permittivity_function = DielectricFUnction(epsType=espType,volume=volume,parameters=permittivity)
        self.permittivity_function.setEpsilonInfinity(epsinf)

    def setPermittivity(self):
        return permittvity

    def setDensity(self):
        return density

    def getPermittivityFunction(self):
        return self.permittivityFunction

    def getDensity(self):
        return density

materials_database = {}
materials_database['air'] = Material('air', density=1.225e-3, epsType=constant, permittivity=1.0)
materials_database['vacuum'] = Material('vacuum', density=0.0, epsType=constant, permittivity=1.0)
materials_database['ptfe']  = Material('ptfe',density=2.2,epsType=constant, permittivity=2.0)
materials_database['kbr']  = Material('kbr',density=2.75,epsType=constant, permittivity=2.25)
materials_database['nujol']  = Material('nujol',density=0.838,epsType=constant, permittivity=2.155)
materials_database['ldpe']  = Material('ldpe',density=0.925,epsType=constant, permittivity=2.25)
materials_database['mdpe']  = Material('mdpe',density=0.933,epsType=constant, permittivity=2.25)
materials_database['silicon']  = Material()
materials_database['aluminium']  = Material()
materials_database['SnO2']  = Material()
materials_database['water']  = Material()
materials_database['water']  = Material()
