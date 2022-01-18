#!/usr/bin/env python
"""Replace the old Makefile system for testing and installing """

import os
import sys
from PDielec.preader import main as main_preader
from PDielec.checkcsv import main as main_checkcsv
import contextlib
from termcolor import colored

test_preader = [
	'Castep/preader',
	'Vasp/preader',
	'Phonopy/preader',
	'Gulp/preader',
	'Crystal/preader',
	'AbInit/preader',
	'QE/preader',
    ]

test_p2cif = [
	'P2Cif/',
    ]

test_pdgui = [
	'Castep/AsparticAcid',
	'Castep/Isoleucine',
	'Castep/MgO',
	'Castep/Bubbles',
	'Castep/Na2SO42',
	'Castep/Castep17',
	'Vasp/F-Apatite',
	'Vasp/ZnO',
	'Vasp/Na2SO42',
	'Gulp/calcite',
	'Gulp/Na2SO42',
	'Crystal/Na2SO42',
	'Crystal/Leucine',
	'Crystal/Quartz',
	'Crystal/ZnO/Default',
	'Crystal/ZnO/CPHF',
	'Crystal/ZnO/NoEckart',
	'AbInit/AlAs',
	'AbInit/BaTiO3',
	'AbInit/Na2SO42',
	'QE/ZnO',
	'QE/Na2SO42',
	'QE/Cocaine',
	'Phonopy/ZnO',
	'Phonopy/Na2SO42',
	'SizeEffects/BaTiO3',
	'SizeEffects/MgO',
	'SizeEffects/ZnO',
	'Mie/MgO',
	'Mie/MgO_lognormal',
	'ATR/AlAs',
	'ATR/Na2SO42',
	'ATR/Na2SO42_fit',
	'Experiment/constant',
	'Experiment/fpsq',
	'Experiment/drude-lorentz',
	'Experiment/interpolation',
	'Experiment/Mayerhofer',
	'SingleCrystal/Bi2Se3',
	'SingleCrystal/Bi2Se3_film',
    ]

test_vibanalysis = [
	'VibAnalysis/AsparticAcid',
	'VibAnalysis/BaTiO3',
	'VibAnalysis/Cocaine',
	'VibAnalysis/F-Apatite',
	'VibAnalysis/Isoleucine',
	'VibAnalysis/Leucine',
    ]

benchmarks = [
	'Castep/AsparticAcid',
	'Castep/Isoleucine',
	'Castep/MgO',
	'Castep/Bubbles',
	'Castep/Na2SO42',
	'Vasp/F-Apatite',
	'Vasp/ZnO',
	'Vasp/Na2SO42',
	'Gulp/calcite',
	'Gulp/Na2SO42',
	'Crystal/Na2SO42',
	'Crystal/Leucine',
	'Crystal/Quartz',
	'Crystal/ZnO/Default',
	'Crystal/ZnO/CPHF',
	'Crystal/ZnO/NoEckart',
	'AbInit/AlAs',
	'AbInit/BaTiO3',
	'AbInit/Na2SO42',
	'QE/ZnO',
	'QE/Na2SO42',
	'QE/Cocaine',
	'Phonopy/ZnO',
	'Phonopy/Na2SO42',
	'SizeEffects/BaTiO3',
	'SizeEffects/MgO',
	'SizeEffects/ZnO',
	'Mie/MgO',
	'Mie/MgO_lognormal',
	'ATR/AlAs',
	'ATR/Na2SO42',
    ]

def usage():
    print('pdmake:', file=sys.stderr)
    print('  test' , file=sys.stderr)
    print('  tests' , file=sys.stderr)
    print('           run all the tests               ' , file=sys.stderr)
    print('  test-preader' , file=sys.stderr)
    print('           run all the preader tests' , file=sys.stderr)
    print('  test-pdgui' , file=sys.stderr)
    print('           run all the pdgui tests' , file=sys.stderr)
    print('  test-p2cif' , file=sys.stderr)
    print('           run all the p2cif tests' , file=sys.stderr)
    print('  test-vibanalysis' , file=sys.stderr)
    print('           run all the vibanalysis tests' , file=sys.stderr)
    print('  benchmarks      ' , file=sys.stderr)
    print('           run the benchmarks   ' , file=sys.stderr)
    print('  regenerate      ' , file=sys.stderr)
    print('           regenerate all the test reference data' , file=sys.stderr)
    print('  root directory  ' , file=sys.stderr)
    print('           sets the root directory for the installation' , file=sys.stderr)
    print('  pypi      ' , file=sys.stderr)
    print('           Only use on linux installations...          ' , file=sys.stderr)
    print('           prepare for pypi uploading' , file=sys.stderr)
    print('  pyinstaller' , file=sys.stderr)
    print('           Only use on linux installations...          ' , file=sys.stderr)
    print('           prepare for pypinstaller' , file=sys.stderr)
    print('  clean' , file=sys.stderr)
    print('           Only use on linux installations...          ' , file=sys.stderr)
    print('           clean up' , file=sys.stderr)
    print('  install' , file=sys.stderr)
    print('           Only use on linux installations...          ' , file=sys.stderr)
    print('           install the executables to the scripts directory' , file=sys.stderr)
    print('  scripts directory' , file=sys.stderr)
    print('           Only use on linux installations...          ' , file=sys.stderr)
    print('           Set the scripts directory to be used for installation' , file=sys.stderr)
    print('           the default is ~/bin                       ' , file=sys.stderr)
    exit()

def install(scripts):
    print('Performing installation of scripts into ',scripts)
    os.popen('echo cp -P preader {}'.format(scripts))
    os.popen('echo cp -P p1reader {}'.format(scripts))
    os.popen('echo cp -P pdgui {}'.format(scripts))
    os.popen('echo cp -P p2cif {}'.format(scripts))
    os.popen('echo cp -P pdcompare {}'.format(scripts))
    os.popen('echo cp -P graphdatagenerator {}'.format(scripts))
    os.popen('echo cp -P vibanalysis {}'.format(scripts))
    os.popen('echo mkdir -p {}/PDielec'.format(scripts))
    os.popen('echo mkdir -p {}/PDielec/GUI'.format(scripts))
    os.popen('echo cp -r PDielec/*.py {}/PDielec'.format(scripts))
    os.popen('echo cp -r PDielec/GUI/*.py {}/PDielec/GUI/'.format(scripts))
    os.popen('echo cp -r PDielec/GUI/*.png {}/PDielec/GUI/'.format(scripts))

@contextlib.contextmanager
def redirect(file):
    # capture all outputs to a log file while still printing it
    class Logger:
        def __init__(self, file):
            fd = open(file,'w')
            self.terminal = sys.stdout
            self.log = fd

        def write(self, message):
            #self.terminal.write(message)
            self.log.write(message)

        def __getattr__(self, attr):
            return getattr(self.terminal, attr)

    logger = Logger(file)

    sys.stdout = logger
    sys.stderr = open('/dev/null','w')
    try:
        yield logger.log
    finally:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__



def runPreaderTest(title,instructions, regenerate):
    '''Run a preader test
    title is the title in the pdmake file
    instructions are the command line parameters
    regenerate is set to true if the reference file is overwritten'''
    sys.argv = []
    sys.argv.append('preader')
    sys.argv.extend(instructions)
    outputfile = 'command.csv'
    if regenerate:
        outputfile = 'command.ref.csv'
    with redirect(outputfile):
        result = main_preader()
    # If not doing a regeneration perform a check
    if not regenerate:
        sys.argv = ['checkcsv', 'command.ref.csv','command.csv']
        with redirect('/dev/null'):
            result = main_checkcsv()
        nerrors,keep_line_number,keep_word1,keep_word2,max_percentage_error = result
        if nerrors > 0:
            print(title+colored(' ERRORS:','red')+"LARGEST ON LINE {} -- max %error={}".format(nerrors, keep_line_number, max_percentage_error))
        else:
            print(title+colored(' OK:','blue')+" -- max %error={}" .format(max_percentage_error))
        # end if
    # end
    return

def runTests(testlist, testType, regenerate):
    '''Run the tests given in the directories stored in testlist
       testType can be pdgui, preader, p2cif or benchmark
       regenerate causes the test data to be regenerated'''
    homedir = os.getcwd()
    from PDielec.pdgui import main as main_pdgui
    from PDielec.p2cif import main as main_p2cif
    from PDielec.VibAnalysis import main as main_vibanalysis
    arguments_pdgui = [ 'pdgui','-nosplash', '-excel results.xlsx' '-script script.py', '-exit']
    arguments_preader = [ 'preader','-nosplash', '-excel results.xlsx' '-script script.py', '-exit']
    for directory in testlist:
        print('Testing {} in directory {}'.format(testType,directory))
        if not os.path.isdir(directory):
            print('Error: command needs to be executed in the PDielec home directory')
            print('       current directory is {}'.format(homedir))
            print('       required directory is {}'.format(directory))
            exit()
        mychdir(directory)
        if testType == 'pdgui':
            main_pdgui( arguments )
        elif testType == 'preader':
            main_preader( arguments)
        elif testType == 'p2cif':
            main_p2cif( arguments)
        elif testType == 'vibanalysis':
            main_vibanalysis( arguments)
        elif testType == 'benchmarks':
            main_pdgui( arguments )
        else:
            print('Error in runTests, unkown testType')
            exit()
        mychdir(homedir)
    # end for
    print('Tests complete')

def readPdMakefile(filename):
    # Return a dictionary of instructions
    # The first word is the dictionary key, the rest is its value
    instructions = {}
    with open(filename,'r') as fd:
        line = fd.readline()
        title = line[:-1].ljust(30)
        line = fd.readline()
        while line:
            if len(line)>0 and not line.startswith('#'):
                splits = line.split()
                key = splits[0]
                value = splits[1:]
                instructions[key] = value
            # end if
            line = fd.readline()
        # end while
    #end open
    return title,instructions

def mychdir(directory):
    '''Local version of os.chdir with some testing'''
    if directory != '':
        if os.path.isdir(directory):
            result = os.chdir(directory)
            newdir = os.getcwd()
        else:
            print('Error directory does not exist',directory)
        # end if
    # end if
    return

def runPdMakefile(pdmakefile,regenerate):
    '''Run specific pdMakefile'''
    homedir = os.getcwd()
    directory,filename = os.path.split(pdmakefile)
    mychdir(directory)
    title,instructions = readPdMakefile(filename)
    for key in instructions:
        parameters = instructions[key]
        if key == 'preader':
            runPreaderTest(title,parameters,regenerate)
    mychdir(homedir)

def runClean():
    print('Cleaning old results from the examples directory')
    os.popen('echo find . -name results.xlsx -exec rm -f {} \;')
    os.popen('echo find . -name results.csv -exec rm -f {} \;')
    os.popen('echo find . -name all.cif -exec rm -f {} \;')
    os.popen('echo find . -name \*.nma -exec rm -f {} \;')
    print('Cleaning complete')

def runPyInstaller():
    print('Creating pyinstaller files in dis/pdgui         ')
    os.popen('pyinstaller pdgui.spec -y')
    os.popen('cp -r dist/pdgui/PyQt5/Qt/plugins/platforms dist/pdgui')

def runPyPi():
    print('Creating PyPi distribution files')
    os.popen('rm -rf build dist PDielec.egg-info; python setup.py sdist bdist_wheel')

def main():
    # Start processing the directories
    originalDirectory = os.getcwd()
    if len(sys.argv) <= 1 :
        usage()
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    regenerate = False
    scriptsDirectory='~/bin'
    command = sys.argv[0]
    rootDirectory,command = os.path.split(command)
    actions = []
    pdmakefiles = []
    # Loop over the tokens on the command line
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token == 'test' or token == 'tests':
            actions.append('test all')
            itoken += 1
        elif token == 'test-preader':
            actions.append('test preader')
            itoken += 1
        elif token == 'test-p2cif':
            actions.append('test preader')
            itoken += 1
        elif token == 'test-pdgui':
            actions.append('test pdgui')
            itoken += 1
        elif token == 'test-vibanalysis':
            actions.append('test vibanalysis')
            itoken += 1
        elif token == 'benchmarks':
            actions.append('run benchmarks')
            itoken += 1
        elif token == 'root':
            itoken += 1
            rootDirectory = tokens[itoken]
            itoken += 1
        elif token == 'regenerate':
            regenerate = True
            itoken += 1
        elif token == 'pypi':
            actions.append('run pypi')
            runPyPi()
            itoken += 1
        elif token == 'runPyInstaller':
            actions.append('run pypinstaller')
            runPyInstaller()
            itoken += 1
        elif token == 'runClean':
            actions.append('run clean')
            itoken += 1
        elif token == 'install':
            actions.append('install')
            install(scriptsDirectory)
            itoken += 1
        elif token == 'scripts':
            itoken += 1
            scriptsDirectory = tokens[itoken]
            itoken += 1
        elif token.endswith('.pdmake'):
            pdmakefiles.append(token)
            itoken += 1
        else:
            usage()
            exit();
        # endif
    # end while
    #
    # Change to the rootDirectory
    if not os.path.isdir(rootDirectory):
        print('Error: The root directory is not valid',rootDirectory)
        exit()
    mychdir(rootDirectory)
    # Change directory to the Examples directory
    if not os.path.isdir('Examples'):
        print('Error: command needs to be executed in the PDielec home directory')
        exit()
    #
    # Handle any pdmake files there might be
    #
    for pdmakefile in pdmakefiles:
        runPdMakefile(os.path.join(originalDirectory,pdmakefile),regenerate)
    #
    # Now loop over the actions
    #
    mychdir('Examples')
    for action in actions:
        if action == 'test all':
            runTests(test_preader    ,'preader'    ,regenerate)
            runTests(test_p2cif      ,'p2cif'      ,regenerate)
            runTests(test_pdgui      ,'pdgui'      ,regenerate)
            runTests(test_vibanalysis,'vibanalysis',regenerate)
        elif action == 'test preader':
            runTests(test_preader    ,'preader'    ,regenerate)
        elif action == 'test p2cif':
            runTests(test_p2cif      ,'p2cif'      ,regenerate)
        elif action == 'test pdgui':
            runTests(test_pdgui      ,'pdgui'      ,regenerate)
        elif action == 'test vibanalysis':
            runTests(test_vibanalysis,'vibanalysis',regenerate)
        elif action == 'run benchmarks':
            runTests(benchmarks,'benchmarks',regenerate)
        elif action == 'run pypi':
            runPyPi()
        elif action == 'run pyinstaller':
            runPyInstaller()
        elif action == 'run clean':
            runClean()
        elif action == 'install':
            install(scriptsDirectory)
        else:
            usage()
            exit();
        #end if
    #end for

# end of def main

if __name__ == '__main__':
    main(sys)
