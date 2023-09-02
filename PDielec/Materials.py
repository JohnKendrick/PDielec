import numpy as np
from PDielec.DielectricFunction import DielectricFunction
from PDielec.Utilities          import Debug
import openpyxl as xl


class MaterialsDataBase():
    def __init__(self,filename, debug=False):
        """Initialise a database of material properties using the excel spreadsheet filename
           The class has the following variables
           filename:        The filename of the spreadsheet
           sheets:          a list of the sheet names in the spreadsheet
        """
        global debugger
        debugger = Debug(debug,'MaterialsDataBase')
        debugger.print('Start:: initialise')
        if os.path.isfile(filename):
            self.filename = filename
            workbook = xl.load_workbook(self.filename)
            self.sheets = self.workbook.sheetnames
            debugger.print('Sheet names:: ',self.sheets)
            # Close the work book while it is not in use
            workbook.close()
        else:
             self.filename = None
             self.sheetnames = None
             print('  Error: MaterialsDataBase filename not valid',filename)
    
        debugger.print('Finished:: initialise')
        return

    def getSheetNames(self):
        """Return a list of the sheetnames in the database"""
        return self.sheetnames

    def getMaterial(self,sheet):
        """Return a material object based on the data in sheet (an excel sheet)"""
            material = Material('undefined')
            if not sheet in self.sheets:
                print('Error in getMaterial sheet ',sheet,' not in self.sheets',self.sheets,file=sys.stderr)
                return material
            workbook = xl.load_workbook(self.filename)
            worksheet = workbook[sheet]
            entry = worksheet['H1'].value
            density = worksheet['H2'].value
            if density == '':
                print('Error in getMaterial density was not defined',file=sys.stderr)
                return material
            if 'Constant' in entry and 'refractive' in entry:
                # Constant refractive index
                n = worksheet['C2'].value
                k = worksheet['D2'].value
                nk = complex(n, k)
                permittivity = Calculator.calculate_permittivity(nk)
                material = Material(sheet,epsType='constant',permittivityValues=permittivity,density=density)
            elif 'Constant' in entry and ('permitt' in entry or 'dielec' in entry):
                # Constant permittivity
                eps_r = worksheet['C2'].value
                eps_i = worksheet['D2'].value
                permittivity = complex(eps_r, eps_i)
                material = Material(sheet,epsType='constant',permittivityValues=permittivity,density=density)
            elif 'Tabulated' in entry and ('permitt' in entry or 'dielec' in entry):
                # Tabulated refractive index
                parameters = []
                for a, c, d in zip(worksheet['A'] ,worksheet['C'] , worksheet['D']):
                    n = c.value
                    k = d.value
                    nk = complex(n, k)
                    permittivity = Calculator.calculate_permittivity(nk)
                    permittivities.append(permittivity)
                material = Material(sheet,epsType='interpolate',permittivityValues=permittivity,density=density)

            workbook.close()
            return material
        
        

class Material():
    def __init__(self, name, density=None, epsinf=0.0, volume=0):
        """A Material has the following properties
           density:            The density of the material
           permitivity:        A function to evaluate the permittivity at a given value
           name:               The name of the material

           The material object is created from;
           name:               The name of the material
           density:            The material density
           epsType:            The type of permittivity function that is used to define permittivity
           epsInf:             Epsilon infinity if needed
           volume:             The volume of the unit cell if needed

        """
        if density is None:
            print('Error in Material',name,'density not specified')
            return
        if permittivityValues is None:
            print('Error in Material',name,'permittivityValues not specified')
            return
        self.density = density
        self.name = name
        self.volume = volume
        self.type = 'Base Class'
        self.permittivity.setEpsilonInfinity(epsinf)

    def getPermittivity(self):
        """Return the permittivity function"""
        return permittvity

    def getDensity(self):
        """Return the permittivity function"""
        return density

class Constant(Material):
    def __init__(self, name, permittivity=None, density=None, epsinf=0.0)
        """Create an instance of a material with a constant permittivity
           permittivity is the value of the permittivy and can be complex
           the Dielectric Function always returns a tensor
           The required parameters are;
           name:             The name of the material
           permittivity      The permittivity  
           density           in g/ml
        """
        self.type = 'Constant permittivity'
        self.permittivity = DielectricFunction.Constant(permittivity)
        Material.__init__(name, density=density, epsinf=epsinf)


class Tabulated(Material):
    def __init__(self, name, vs_cm1=None, permittivities=None, density=None, epsinf=0.0)
        """Create an instance of a material with a constant permittivity
           permittivity is the value of the permittivy and can be complex
           the Dielectric Function always returns a tensor
           The required parameters are;
           name:             The name of the material
           vs_cm1            The list of tabulated frequencies in cm-1
           permittivities    The permittivites, either a single (n) vector or a (3,n) vector
           density           in g/ml
        """
        vs = np.array(vs_cm1)
        eps = np.array(permittivities)
        if len(np.shape(eps)) ==2:
            if np.shape(eps)[0] == 3:
                self.permittivity = DielectricFunction.Tabulate3(vs,eps)
            elif np.shape(eps)[0] == 6:
                self.permittivity = DielectricFunction.Tabulate6(vs,eps)
            else:
                print('Error in Tabulated, shape of parameters is wrong')
         else:
             self.permittivity = DielectricFunction.Tabulate1(vs,eps)
        self.type = 'Tabulated permittivity'
        Material.__init__(name, density=density, epsinf=epsinf)

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
