import numpy as np
import PDielec.DielectricFunction as DielectricFunction
import PDielec.Calculator as Calculator
from PDielec.Utilities          import Debug
from PDielec.UnitCell           import UnitCell
from PDielec                    import __file__  as PDielec_init_filename
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
        if len(filename)> 5 and (filename.endswith('xlsx') or filename.endswith('XLSX')) and os.path.isfile(filename):
            self.filename = os.path.relpath(filename)
            workbook = xl.load_workbook(self.filename,data_only=True)
            self.sheetNames = workbook.sheetnames
            debugger.print('Sheet names:: ',self.sheetNames)
            # Close the work book while it is not in use
            workbook.close()
        else:
            # Try opening the default database
            PDielec_Directory = os.path.dirname(PDielec_init_filename)
            filename  = os.path.join(PDielec_Directory, 'MaterialsDataBase.xlsx')
            filename  = os.path.relpath(filename)
            if os.path.isfile(filename):
                self.filename = filename
                workbook = xl.load_workbook(self.filename,data_only=True)
                self.sheetNames = workbook.sheetnames
                debugger.print('Sheet names from default database ',self.sheetNames)
                # Close the work book while it is not in use
                workbook.close()
            else:
                self.filename = None
                self.sheetNames = None
                print('  Error: MaterialsDataBase filename not valid',filename)
    
        debugger.print('Finished:: initialise')
        return

    def getFileName(self):
        '''Return the filename'''
        return self.filename

    def valid(self):
        '''Test to see if the spreadsheet is a valid materials database'''
        result = False
        if 'Information' in self.sheetNames[0]:
            result = True
        return result

    def getSheetNames(self):
        '''Return a list of the sheetnames in the database'''
        # First take a copy of the sheetnames ignoring the first (Information)
        fullList = []
        if self.sheetNames is not None:
            fullList = self.sheetNames[1:].copy()
        # Append any in-built materials
        # This list is taken from the original powder code before version 8.0
        if 'air' not in fullList:
            fullList.append('air')
        if 'vacuum' not in fullList:
            fullList.append('vacuum')
        if 'ptfe' not in fullList:
            fullList.append('ptfe')
        if 'ldpe' not in fullList:
            fullList.append('ldpe')
        if 'mdpe' not in fullList:
            fullList.append('mdpe')
        if 'kbr' not in fullList:
            fullList.append('kbr')
        if 'nujol' not in fullList:
            fullList.append('nujol')
        debugger.print('getSheetNames:: ',fullList)
        return sorted(fullList, key=lambda s: s.casefold())

    def getMaterial(self,sheet):
        '''Return a material object based on the data in sheet (an excel sheet)'''
        debugger.print('getMaterial:: ',sheet)
        # Define a set of back-up materials that the program can use even if the sheet name is not in the spreadsheet
        if self.sheetNames is None or sheet not in self.sheetNames:
            if sheet == 'air':
                material = Constant('air',permittivity=1.0,density=0.001225)
            elif sheet == 'vacuum':
                material = Constant('vacuum',permittivity=1.0,density=0.0)
            elif sheet == 'ptfe':
                material = Constant('ptfe',permittivity=1.0,density=2.2)
            elif sheet == 'ldpe':
                material = Constant('ldpe',permittivity=2.25,density=0.925)
            elif sheet == 'mdpe':
                material = Constant('mdpe',permittivity=2.25,density=0.933)
            elif sheet == 'kbr':
                material = Constant('kbr',permittivity=2.25,density=2.75)
            elif sheet == 'nujol':
                material = Constant('nujol',permittivity=2.155,density=0.838)
            else:
                print('Error in getMaterial sheet ',sheet,' not in self.sheetNames',self.sheetNames,file=sys.stderr)
                material = Constant('vacuum',permittivity=1.0,density=0.0)
            return material
        # Carry on with the spreadsheet
        workbook = xl.load_workbook(self.filename,data_only=True)
        worksheet = workbook[sheet]
        unitCell = None
        avector = bvector = cvector = None
        a = b = c = alpha = beta = gamma = None
        for i in range(20):
            cell1 = 'G'+str(i+1)
            cell2 = 'H'+str(i+1)
            token = worksheet[cell1].value
            if token is not None:
                token = token.lower()
                if 'entry' in token:
                    entry = worksheet[cell2].value.lower()
                elif 'density' in token:
                    density = float(worksheet[cell2].value)
                elif 'a_vector' in token:
                    avector = [ float(cell.value) for cell in [ worksheet['I'+str(i+1)], worksheet['J'+str(i+1)], worksheet['K'+str(i+1)] ] ]
                elif 'b_vector' in token:
                    bvector = [ float(cell.value) for cell in [ worksheet['I'+str(i+1)], worksheet['J'+str(i+1)], worksheet['K'+str(i+1)] ] ]
                elif 'c_vector' in token:
                    cvector = [ float(cell.value) for cell in [ worksheet['I'+str(i+1)], worksheet['J'+str(i+1)], worksheet['K'+str(i+1)] ] ]
                elif 'a:' == token:
                    a = float(worksheet[cell2].value)
                elif 'b:' == token:
                    b = float(worksheet[cell2].value)
                elif 'c:' == token:
                    c = float(worksheet[cell2].value)
                elif 'alpha' in token:
                    alpha = float(worksheet[cell2].value)
                elif 'beta' in token:
                    beta = float(worksheet[cell2].value)
                elif 'gamma' in token:
                    gamma = float(worksheet[cell2].value)
        #
        if avector is not None and bvector is not None and cvector is not None:
            unitCell = UnitCell(a=avector,b=bvector,c=cvector)
        elif a is not None and b is not None and c is not None and alpha is not None and beta is not None and gamma is not None:
            unitCell = UnitCell(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)
        # Process the entry type
        if 'constant' in entry and 'refractive' in entry:
            material = self.readConstantRefractiveIndex(sheet,worksheet,density)
        elif 'constant' in entry and ('permitt' in entry or 'dielec' in entry):
            material = self.readConstantPermittivity(sheet,worksheet,density)
        elif 'tabulated' in entry and 'refractive' in entry:
            material = self.readTabulatedRefractiveIndex(sheet,worksheet,density)
        elif 'tabulated' in entry and ('permitt' in entry or 'dielec' in entry):
            material = self.readTabulatedPermittivity(sheet,worksheet,density)
        elif 'lorentz' in entry and 'drude' in entry:
            material = self.readLorentzDrude(sheet,worksheet,density,unitCell)
        elif 'fpsq' in entry:
            material = self.readFPSQ(sheet,worksheet,density,unitCell)
        elif 'sellmeier' in entry:
            material = self.readSellmeier(sheet,worksheet,density,unitCell)
        # Close the work book
        workbook.close()
        return material

    def readConstantRefractiveIndex(self,sheet,worksheet,density):
        '''Read constant refractive index from the spreadsheet
           sheet is the worksheet name
           worksheet is the worksheet
           density is the density of the material'''
        # Constant refractive index
        n = float(worksheet['C2'].value)
        k = float(worksheet['D2'].value)
        nk = complex(n, k)
        permittivity = Calculator.calculate_permittivity(nk)
        debugger.print('Constant refractive:: ',nk,permittivity,density)
        material = Constant(sheet,permittivity=permittivity,density=density)
        return material

    def readConstantPermittivity(self,sheet,worksheet,density):
        '''Read constant permittivity data from the spreadsheet
           sheet is the worksheet name
           worksheet is the worksheet
           density is the density of the material'''
        # Constant permittivity
        eps_r = float(worksheet['C2'].value)
        eps_i = float(worksheet['D2'].value)
        permittivity = complex(eps_r, eps_i)
        debugger.print('Constant permittivity:: ',permittivity,density)
        material = Constant(sheet,permittivity=permittivity,density=density)
        return material

    def readTabulatedRefractiveIndex(self,sheet,worksheet,density):
        '''Read tabulated refractive index data from the spreadsheet
           sheet is the worksheet name
           worksheet is the worksheet
           density is the density of the material'''
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
        return material

    def readTabulatedPermittivity(self,sheet,worksheet,density):
        '''Read tabulated permittivity data from the spreadsheet
           sheet is the worksheet name
           worksheet is the worksheet
           density is the density of the material'''
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
        return material

    def readLorentzDrude(self,sheet,worksheet,density,unitCell):
        '''Read Drude-Lorentz data from the spreadsheet
           sheet is the worksheet name
           worksheet is the worksheet
           density is the density of the material 
           unitCell unit cell'''
        # Lorentz-Drude model for permittivity
        epsilon_infinity = np.zeros( (3,3) )
        directions = [[], [], []]
        epsinfs = [[], [], []]
        omegas = [[], [], []]
        strengths = [[], [], []]
        gammas = [[], [], []]
        for a, b, c, d, e in zip(worksheet['A'][1:] ,worksheet['B'][1:] , worksheet['C'][1:], worksheet['D'][1:], worksheet['E'][1:]) :
            try:
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
            except:
                print('Error in Lorentz-Drude: ',a.value,b.value,c.value,d.value,e.value)
                return
        material = DrudeLorentz(sheet,epsilon_infinity,omegas,strengths,gammas,density=density,cell=unitCell)
        return material

    def readFPSQ(self,sheet,worksheet,density,unitCell):
        '''Read FPSQ data from the spreadsheet
           sheet is the worksheet name
           worksheet is the worksheet
           density is the density of the material 
           unitCell unit cell'''
        # FPSQ model for permittivity
        epsilon_infinity = np.zeros( (3,3) )
        directions = [[], [], []]
        omega_tos = [[], [], []]
        gamma_tos = [[], [], []]
        omega_los = [[], [], []]
        gamma_los = [[], [], []]
        for a, b, c, d, e, f in zip(worksheet['A'][1:] ,worksheet['B'][1:] , worksheet['C'][1:], worksheet['D'][1:], worksheet['E'][1:], worksheet['F'][1:]) :
            try:
                if a.value is not None:
                    direction = a.value
                index = ['xx','yy','zz'].index(direction)
                if b.value is not None:
                    epsilon_infinity[[index],[index]] = float(b.value)
                if c.value is not None:
                    omega_tos[index].append(float(c.value))
                if d.value is not None:
                    gamma_tos[index].append(float(d.value))
                if e.value is not None:
                    omega_los[index].append(float(e.value))
                if f.value is not None:
                    gamma_los[index].append(float(f.value))
            except:
                print('Error in FPSQ: ',a.value,b.value,c.value,d.value,e.value,f.value)
                return
        material = FPSQ(sheet,epsilon_infinity,omega_tos,gamma_tos,omega_los,gamma_los,density=density,cell=unitCell)
        return material

    def readSellmeier(self,sheet,worksheet,density,unitCell):
        '''Read Sellmeier data from the spreadsheet
           sheet is the worksheet name
           worksheet is the worksheet
           density is the density of the material 
           unitCell unit cell'''
        # Sellmeier model for refractive index
        directions = []
        Bs = []
        Cs = []
        for b, c in zip(worksheet['A'][1:] , worksheet['B'][1:] ) :
            try:
                if b.value is not None:
                    Bs.append(float(b.value))
                if c.value is not None:
                    Cs.append(float(c.value))
            except:
                print('Error in Sellmeier: ',b.value,c.value)
                return
        material = Sellmeier(sheet,Bs,Cs,density=density,cell=unitCell)
        return material
        
class Material():
    def __init__(self, name, density=None, permittivityObject=None, cell=None):
        '''A Material has the following properties
           density:            The density of the material
           permitivityObject:  An instance of a DielectricFunction
                               The permittivity object is created by the children if Material
           name:               The name of the material

           The material object is created from;
           name:               The name of the material
           density:            The material density
           cell:               The unit cell

        '''
        self.density            = density
        self.cell               = cell
        self.name               = name
        self.type               = 'Base Class'
        self.permittivityObject = permittivityObject
        if self.density is None and self.cell is not None:
            self.density = self.cell.calculate_density()

    def getName(self):
        return self.name

    def getInformation(self):
        '''Return information about the material'''
        result = self.type
        if 'Tabulate' in self.type:
            low = self.permittivityObject.getLowestFrequency()
            high = self.permittivityObject.getHighestFrequency()
            result += ' freq range {:.0f}-{:.0f}'.format(low,high) # + ' value at 0 {}'.format(self.permittivityObject.function()(0))
        # result += ' value at 0 {}'.format(self.permittivityObject.function()(0))
        return result

    def print(self):
        '''Print information about the material'''
        print('Material name:',self.name)
        print('Material density:',self.density)
        print('Material type:',self.type)
        print('Material is scalar?:',self.isScalar())
        print('Material is tensor?:',self.isTensor())
        print('Material permittivity:',self.getInformation())
        if self.cell is not None:
            print('Material unit cell')
            self.cell.print()
        return

    def isScalar(self):
        '''Returns true of the material returns a scalar permittivity'''
        return self.permittivityObject.isScalar()

    def isTensor(self):
        '''Returns true of the material returns a tensor permittivity'''
        return self.permittivityObject.isTensor()

    def setPermittivityObject(self,permittivityObject):
        '''Return the permittivity function'''
        self.permittivityObject = permittivityObject
        return

    def getPermittivityObject(self):
        '''Return the permittivity function'''
        return self.permittivityObject

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

    def setCell(self, cell):
        '''Return the cell'''
        self.cell = cell
        if self.density is None and self.cell is not None:
            self.density = self.cell.calculate_density()
        return

    def getCell(self):
        '''Return the cell'''
        return self.cell

class Constant(Material):
    def __init__(self, name, permittivity=None, density=None, cell=None):
        '''Create an instance of a material with a constant scalar permittivity
           permittivity is the value of the permittivity and can be complex
           The required parameters are;
           name:             The name of the material
           permittivity      The permittivity  
           density           in g/ml
           cell              the unit cell
        '''
        super().__init__(name, density=density, permittivityObject=DielectricFunction.ConstantScalar(permittivity), cell=cell)
        self.type = 'Constant permittivity'


class External(Material):
    def __init__(self, name, permittivityObject=None, density=None, cell=None):
        '''Create an instance of a material which has the permittivity object specified externally
           permittivity is the value of the permittivity and can be complex
           The required parameters are;
           name:              The name of the material
           permittivityObject The permittivity  
           density            in g/ml
           cell               the unit cell
        '''
        super().__init__(name, density=density, permittivityObject=permittivityObject, cell=cell)
        self.type = 'External permittivity'


class DrudeLorentz(Material):
    def __init__(self, name,epsinf,omegas,strengths,gammas,density=None,cell=None):
        '''Create an instance of a material with a Lorentz Drude model permittivity
           permittivity is the value of the permittivy and can be complex
           The required parameters are;
           name:             The name of the material
           epsinf            Epsilon infinity either a 3x3 list or a 3x3 array
           omegas            The TO frequencies as a list 
           strengths         The absorption strengths as a list 
           gammas            The absorption widths as a list
           density           in g/ml
           cell               the unit cell
        '''
        epsilon_infinity = np.array(epsinf)
        permittivityObject = DielectricFunction.DrudeLorentz( omegas, strengths, gammas)
        permittivityObject.setEpsilonInfinity(epsilon_infinity)
        super().__init__(name, density=density, permittivityObject=permittivityObject,cell=cell)
        self.type = 'Drude-Lorentz'

class FPSQ(Material):
    def __init__(self, name,epsinf,omega_tos,gamma_tos,omega_los,gamma_los,density=None,cell=None):
        '''Create an instance of a material with an FPSQ model permittivity
           permittivity is the value of the permittivy and can be complex
           The required parameters are;
           name:             The name of the material
           eps0              Epsilon infinity either a 3x3 list or a 3x3 array
           omega_tos         The TO frequencies as a list 
           gamma_tos         The TO absorption widths as a list
           omega_los         The LO frequencies as a list 
           gamma_los         The LO absorption widths as a list
           density           in g/ml
           cell              the unit cell
        '''
        epsilon_infinity = np.array(epsinf)
        permittivityObject = DielectricFunction.FPSQ( omega_tos, gamma_tos, omega_los, gamma_los)
        permittivityObject.setEpsilonInfinity(epsilon_infinity)
        super().__init__(name, density=density, permittivityObject=permittivityObject,cell=cell)
        self.type = 'FPSQ'

class Sellmeier(Material):
    def __init__(self, name,Bs,Cs,density=None,cell=None):
        '''Create an instance of a material with an Sellmeier model permittivity
           permittivity is the value of the permittivy and should be real for Sellmeier
           The required parameters are;
           name:             The name of the material
           Bs                The B parameters of the Sellmeier equation
           Cs                The C parameters of the Sellmeier equation
           density           in g/ml
           cell               the unit cell
        '''
        permittivityObject = DielectricFunction.Sellmeier( Bs, Cs)
        super().__init__(name, density=density, permittivityObject=permittivityObject,cell=cell)
        self.type = 'Sellmeier'


class Tabulated(Material):
    def __init__(self, name, vs_cm1=None, permittivities=None, density=None, cell=None):
        '''Create an instance of a material with a constant permittivity
           permittivity is the value of the permittivy and can be complex
           The returned permittivityObject can generate either a scalar or a tensor.
           For defining a support matrix material a scalar is used
           The required parameters are;
           name:             The name of the material
           vs_cm1            The list of tabulated frequencies in cm-1
           permittivities    The permittivites, either a single (n) vector or a (3,n) vector
           density           in g/ml
           cell               the unit cell
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
        super().__init__(name, density=density, permittivityObject=permittivityObject,cell=cell)
        self.type = 'Tabulated permittivity'

