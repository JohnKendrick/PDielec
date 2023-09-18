import numpy as np
import PDielec.DielectricFunction as DielectricFunction
import PDielec.Calculator as Calculator
from PDielec.Utilities          import Debug
import openpyxl as xl
import os
import sys


class MaterialsDataBase():
    def __init__(self,filename, debug=False):
        '''Initialise a database of material properties using the excel spreadsheet filename
           The class has the following variables
           filename:        The filename of the spreadsheet
           sheetNames:      a list of the sheet names in the spreadsheet
        '''
        global debugger
        debugger = Debug(debug,'MaterialsDataBase')
        debugger.print('Start:: initialise')
        if os.path.isfile(filename):
            self.filename = filename
            workbook = xl.load_workbook(self.filename,data_only=True)
            self.sheetNames = workbook.sheetnames
            debugger.print('Sheet names:: ',self.sheetNames)
            # Close the work book while it is not in use
            workbook.close()
        else:
             self.filename = None
             self.sheetNames = None
             print('  Error: MaterialsDataBase filename not valid',filename)
    
        debugger.print('Finished:: initialise')
        return

    def valid(self):
        '''Test to see if the spreadsheet is a valid materials database'''
        result = False
        if 'Information' in self.sheetNames[0]:
            result = True
        return result

    def getSheetNames(self):
        '''Return a list of the sheetnames in the database'''
        debugger.print('getSheetNames:: ',self.sheetNames)
        return self.sheetNames[1:]

    def getMaterial(self,sheet):
        '''Return a material object based on the data in sheet (an excel sheet)'''
        debugger.print('getMaterial:: ',sheet)
        if not sheet in self.sheetNames:
            print('Error in getMaterial sheet ',sheet,' not in self.sheetNames',self.sheetNames,file=sys.stderr)
            return material
        workbook = xl.load_workbook(self.filename,data_only=True)
        worksheet = workbook[sheet]
        entry = worksheet['H1'].value
        density = float(worksheet['H2'].value)
        if density is None :
            print('Error in getMaterial density was not defined',file=sys.stderr)
            debugger.print('type(density): ',type(density))
            return material
        if 'Constant' in entry and 'refractive' in entry:
            # Constant refractive index
            n = float(worksheet['C2'].value)
            k = float(worksheet['D2'].value)
            nk = complex(n, k)
            permittivity = Calculator.calculate_permittivity(nk)
            debugger.print('Constant refractive:: ',nk,permittivity,density)
            material = Constant(sheet,permittivity=permittivity,density=density)
        elif 'Constant' in entry and ('permitt' in entry or 'dielec' in entry):
            # Constant permittivity
            eps_r = float(worksheet['C2'].value)
            eps_i = float(worksheet['D2'].value)
            permittivity = complex(eps_r, eps_i)
            debugger.print('Constant permittivity:: ',permittivity,density)
            material = Constant(sheet,permittivity=permittivity,density=density)
        elif 'Tabulated' in entry and 'refractive' in entry:
            # Tabulated refractive index
            permittivities = []
            vs_cm1 = []
            for a, c, d in zip(worksheet['A'][1:] ,worksheet['C'][1:] , worksheet['D'][1:]):
                try:
                   v = float(a.value)
                   n = float(c.value)
                   k = float(d.value)
                except:
                    print('Error in Tabulated: ',a.value,c.value,d.value)
                nk = complex(n, k)
                permittivity = Calculator.calculate_permittivity(nk)
                permittivities.append(permittivity)
                vs_cm1.append(v)
            material = Tabulated(sheet,vs_cm1,permittivities=permittivities,density=density)
        elif 'Tabulated' in entry and ('permitt' in entry or 'dielec' in entry):
            # Tabulated permittivity
            permittivities = []
            vs_cm1 = []
            for a, c, d in zip(worksheet['A'][1:] ,worksheet['C'][1:] , worksheet['D'][1:]):
                v = float(a.value)
                eps_r = float(c.value)
                eps_i = float(d.value)
                permittivity = complex(eps_r, eps_i)
                permittivities.append(permittivity)
                vs_cm1.append(v)
            material = Tabulated(sheet,vs_cm1,permittivities=permittivities,density=density)
        elif 'Lorentz' in entry and 'Drude' in entry:
            # Lorentz-Drude mode for permittivity
            epsilon_infinity = np.zeros( (3,3) )
            directions = [[], [], []]
            epsinfs = [[], [], []]
            omegas = [[], [], []]
            strengths = [[], [], []]
            gammas = [[], [], []]
            for a, b, c, d, e in zip(worksheet['A'][1:] ,worksheet['B'][1:] , worksheet['C'][1:], worksheet['D'][1:], worksheet['E'][1:]) :
                print('Error in Lorentz-Drude: ',a.value,b.value,c.value,d.value,e.value)
#                try:
                if a.value is not None:
                    direction = a.value
                index = ['xx','yy','zz'].index(direction)
                if b.value is not None:
                    epsilon_infinity[[index],[index]] = float(b.value)
                if c.value is not None:
                    omegas[index].append(float(c.value))
                if d.value is not None:
                    strengths[index].append(float(d.value))
                if e.value is not None:
                    gammas[index].append(float(e.value))
#                except:
#                    print('Error in Lorentz-Drude: ',a.value,b.value,c.value,d.value,e.value)
#                    return
            material = DrudeLorentz(sheet,epsilon_infinity,omegas,strengths,gammas,density=density)
        # Close the work book
        workbook.close()
        return material
        
class Material():
    def __init__(self, name, density=None, permittivityObject=None, volume=0):
        '''A Material has the following properties
           density:            The density of the material
           permitivityObject:  An instance of a DielectricFunction
                               The permittivity object is created by the children if Material
           name:               The name of the material

           The material object is created from;
           name:               The name of the material
           density:            The material density
           volume:             The volume of the unit cell if needed

        '''
        self.density = density
        self.name = name
        self.volume = volume
        self.type = 'Base Class'
        self.permittivityObject = permittivityObject

    def getInformation(self):
        '''Return information about the material'''
        result = self.type
        if 'Tabulate' in self.type:
            low = self.permittivityObject.getLowestFrequency()
            high = self.permittivityObject.getHighestFrequency()
            result += ' freq range {:.0f}-{:.0f}'.format(low,high)
        return result

    def setPermittivityObject(self,permittivityObject):
        '''Return the permittivity function'''
        self.permittivityObject = permittivityObject
        return

    def getPermittivityFunction(self):
        '''Return the permittivity function'''
        return self.permittivityObject.function()

    def setDensity(self, value):
        '''Return the density'''
        self.density = value
        return

    def getDensity(self):
        '''Return the density'''
        return self.density

    def print(self):
        '''Return the permittivity function'''
        print('Material ',self.name)
        print('   type: ',self.type)
        print('   volume: ',self.volume)
        print('   density: ',self.density)
        print('   permittivity: ',self.permittivityObject)
        return self.density

class Constant(Material):
    def __init__(self, name, permittivity=None, density=None):
        '''Create an instance of a material with a constant scalar permittivity
           permittivity is the value of the permittivity and can be complex
           The required parameters are;
           name:             The name of the material
           permittivity      The permittivity  
           density           in g/ml
        '''
        super().__init__(name, density=density, permittivityObject=DielectricFunction.ConstantScalar(permittivity))
        self.type = 'Constant permittivity'


class DrudeLorentz(Material):
    def __init__(self, name,epsinf,omegas,strengths,gammas,density=None):
        '''Create an instance of a material with a Lorentz Drude model permittivity
           permittivity is the value of the permittivy and can be complex
           The returned permittivityObject can generate either a scalar or a tensor.
           For defining a support matrix material a scalar is used
           The required parameters are;
           name:             The name of the material
           epsinf            Epsilon infinity either a 3x3 list or a 3x3 array
           omegas            The TO frequencies as a list 
           strengths         The absorption strengths as a list 
           gammas            The absorption widths as a list
           density           in g/ml
        '''
        epsilon_infinity = np.array(epsinf)
        permittivityObject = DielectricFunction.DrudeLorentz( omegas, strengths, gammas,units='cm-1')
        permittivityObject.setEpsilonInfinity(epsilon_infinity)
        super().__init__(name, density=density, permittivityObject=permittivityObject)
        self.type = 'Drude-Lorentz'


class Tabulated(Material):
    def __init__(self, name, vs_cm1=None, permittivities=None, density=None):
        '''Create an instance of a material with a constant permittivity
           permittivity is the value of the permittivy and can be complex
           The returned permittivityObject can generate either a scalar or a tensor.
           For defining a support matrix material a scalar is used
           The required parameters are;
           name:             The name of the material
           vs_cm1            The list of tabulated frequencies in cm-1
           permittivities    The permittivites, either a single (n) vector or a (3,n) vector
           density           in g/ml
        '''
        vs = np.array(vs_cm1)
        eps = np.array(permittivities)
        if len(np.shape(eps)) ==2:
            if np.shape(eps)[0] == 3:
                permittivityObject = DielectricFunction.Tabulate3(vs,eps)
            elif np.shape(eps)[0] == 6:
                permittivityObject = DielectricFunction.Tabulate6(vs,eps)
            else:
                print('Error in Tabulated, shape of parameters is wrong')
        else:
            permittivityObject = DielectricFunction.TabulateScalar(vs,eps)
        super().__init__(name, density=density, permittivityObject=permittivityObject)
        self.type = 'Tabulated permittivity'

