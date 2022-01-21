#!/usr/bin/env python
"""Replace the old Makefile system for testing and installing """

import os
import sys
from PDielec.preader     import main as main_preader
from PDielec.pdgui       import main as main_pdgui
from PDielec.checkcsv    import main as main_checkcsv
from PDielec.checkexcel  import main as main_checkexcel
from PDielec.p2cif       import main as main_p2cif
from PDielec.vibanalysis import main as main_vibanalysis
import numpy as np
import contextlib
from shutil import copyfile
import subprocess
import time

settings = {}
settings['padding'] = 50
settings['title'] = 'title'

debug = False
viewing = False
useLocal = True

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
    print('  --regenerate      ' , file=sys.stderr)
    print('           regenerate all the test reference data' , file=sys.stderr)
    print('  --root directory  ' , file=sys.stderr)
    print('           sets the root directory for the installation' , file=sys.stderr)
    print('  --directory', file=sys.stderr)
    print('           Use the directory when listing the tests' , file=sys.stderr)
    print('  --padding 60', file=sys.stderr)
    print('           change the width of the output title field' , file=sys.stderr)
    print('  --view  ' , file=sys.stderr)
    print('           run pdgui with windowing on' , file=sys.stderr)
    print('  --usesystem ' , file=sys.stderr)
    print('           Force the use of system executables' , file=sys.stderr)
    print('  --debug ' , file=sys.stderr)
    print('           run with debug options switched on' , file=sys.stderr)
    if os.name != 'nt':
        print('  pypi      ' , file=sys.stderr)
        print('           Only use on linux installations...          ' , file=sys.stderr)
        print('           prepare for pypi uploading' , file=sys.stderr)
        print('  clean' , file=sys.stderr)
        print('           Only use on linux installations...          ' , file=sys.stderr)
        print('           clean up' , file=sys.stderr)
        print('  install' , file=sys.stderr)
        print('           Only use on linux installations...          ' , file=sys.stderr)
        print('           install the executables to the scripts directory' , file=sys.stderr)
        print('  --scripts directory' , file=sys.stderr)
        print('           Only use on linux installations...          ' , file=sys.stderr)
        print('           Set the scripts directory to be used for installation' , file=sys.stderr)
        print('           the default is ~/bin                       ' , file=sys.stderr)
    exit()

def checkLocalExecutables():
    '''Check to see if the rootDirectory has the executables we need'''
    os.chdir(rootDirectory)
    if not os.path.isfile('preader'):
        return False
    if not os.path.isfile('p2cif'):
        return False
    if not os.path.isfile('pdgui'):
        return False
    if not os.path.isfile('vibanalysis'):
        return False
    return True

def install(scripts):
    print('Performing installation of scripts into ',scripts)
    global rootDirectory
    mychdir(rootDirectory)
    print('Installing preader',end='')
    subprocess.run('cp -P preader {}'.format(scripts),shell=True)
    print(',p1reader',end='')
    subprocess.run('cp -P p1reader {}'.format(scripts),shell=True)
    print(',pdgui',end='')
    subprocess.run('cp -P pdgui {}'.format(scripts),shell=True)
    print(',p2cif',end='')
    subprocess.run('cp -P p2cif {}'.format(scripts),shell=True)
    print(',pdcompare',end='')
    subprocess.run('cp -P pdcompare {}'.format(scripts),shell=True)
    print(',graphdatagenerator',end='')
    subprocess.run('cp -P graphdatagenerator {}'.format(scripts),shell=True)
    print(',vibanalysis',end='')
    subprocess.run('cp -P vibanalysis {}'.format(scripts),shell=True)
    print(',pdmake',end='')
    subprocess.run('cp -P pdmake {}'.format(scripts),shell=True)
    print(',python scripts')
    subprocess.run('mkdir -p {}/PDielec'.format(scripts),shell=True)
    subprocess.run('mkdir -p {}/PDielec/GUI'.format(scripts),shell=True)
    subprocess.run('cp -r PDielec/*.py {}/PDielec'.format(scripts),shell=True)
    subprocess.run('cp -r PDielec/GUI/*.py {}/PDielec/GUI/'.format(scripts),shell=True)
    subprocess.run('cp -r PDielec/GUI/*.png {}/PDielec/GUI/'.format(scripts),shell=True)
    print('Finished installation of scripts into ',scripts)

@contextlib.contextmanager
def redirect(file):
    # capture all outputs to a log file while still printing it
    class Logger:
        def __init__(self, file):
            fd = open(file,'w')
            self.terminal = sys.stdout
            self.log = fd

        def write(self, message):
            global debug
            if debug:
                self.terminal.write(message)
            self.log.write(message)

        def __getattr__(self, attr):
            return getattr(self.terminal, attr)

    logger = Logger(file)

    sys.stdout = logger
    if not debug:
        sys.stderr = open(os.devnull,'w')
    try:
        yield logger.log
    finally:
        sys.stdout = sys.__stdout__
        if not debug:
            sys.stderr = sys.__stderr__

def readNmaFile(file):
    nmodes = 0
    frequencies = []
    rsquared = []
    explained = []
    with open(file,'r') as fd:
        line = fd.readline()
        while line:
            splits = line.split()
            if len(splits) > 0:
                if splits[0] == 'Mode':
                    if splits[1] != 'too':
                        nmodes += 1
                        frequencies.append(float(splits[2]))
                elif splits[0] == 'R**2':
                    rsquared.append(float(splits[2]))
                elif splits[0] == 'Explained':
                    explained.append(float(splits[3]))
                # end if
            # end if
            line = fd.readline()
        #end while
    #end with open
    f = np.array(frequencies[3:])
    r = np.array(rsquared[3:])
    e = np.array(explained[3:])
    return np.sum(f),np.sum(r),np.sum(e)

def compareNmaFiles(file1,file2):
    f1,r1,e1 = readNmaFile(file1)
    f2,r2,e2 = readNmaFile(file2)
    nerrors = 0
    if abs(f1 - f2) > 1.0e-6:
        nerrors += 1
    if abs(r1 - r2) > 1.0e-6:
        nerrors += 1
    if abs(e1 - e2) > 1.0e-6:
        nerrors += 1
    return nerrors

def compareFiles(file1,file2):
    global debug
    fd1 = open(file1,'r')
    fd2 = open(file2,'r')
    lines1 = fd1.readlines()
    lines2 = fd2.readlines()
    nerrors = 0
    for line1,line2 in  zip(lines1,lines2):
        for word1,word2 in zip(line1.split(),line2.split()):
            word1 = word1.replace('(','')
            word1 = word1.replace(')','')
            word1 = word1.replace('%','')
            word2 = word2.replace('(','')
            word2 = word2.replace(')','')
            word2 = word2.replace('%','')
            try:
                float1 = float(word1)
                float2 = float(word2)
                if abs(float1)+abs(float2) > 1.e-12:
                    if 2*abs(float1-float2)/(abs(float1)+abs(float2)) > 1.0e-4:
                        if debug:
                            print('Float difference',float1,float2)
                            print('Line 1',line1)
                            print('Line 2',line2)
                        nerrors += 1
            except:
                if word1 != word2:
                    if debug:
                        print('Word difference',word1,word2)
                        print('Line 1',line1)
                        print('Line 2',line2)
                    nerrors += 1
    return nerrors


def runP2CifTest(title, instructions, regenerate):
    '''Run a test on p2cif'''
    global rootDirectory
    global useLocal
    if useLocal:
        sys.argv = ['python']
        sys.argv.append(os.path.join(rootDirectory,'p2cif'))
    else:
        sys.argv = ['p2cif']
    sys.argv.extend(instructions)
    outputfile = 'all.cif'
    if regenerate:
        outputfile = 'all.ref.cif'
    with open(outputfile,'w') as stdout:
        if debug:
            result = subprocess.run(sys.argv,stdout=stdout)
        else:
            result = subprocess.run(sys.argv,stdout=stdout,stderr=subprocess.DEVNULL)
        # end if debug
    # end with open
    # If not doing a regeneration perform a check
    if not regenerate:
        nerrors = compareFiles('all.cif', 'all.ref.cif')
        if nerrors > 0:
            print(title+' {} ERRORS:'.format(nerrors))
        else:
            print(title+' OK:')
        # end if
    else:
        print(title+' Regenerated:')
    # end if
    return

def runVibAnalysis(title, instructions, regenerate):
    '''Run a vibanalysis test
    title is the title in the pdmake file
    instructions are the command line parameters
    regenerate is set to true if the reference file is overwritten'''
    global rootDirectory
    global useLocal
    if useLocal:
        sys.argv = ['python']
        sys.argv.append(os.path.join(rootDirectory,'vibanalysis'))
    else:
        sys.argv = ['vibanalysis']
    sys.argv.extend(instructions)
    filename = None
    for option in instructions:
        if not option.startswith('-'):
            filename = option
        # end if
    # end for
    if filename is None:
        print('Error in runVibAnalysis: there is no filename specified')
    header = os.path.splitext(filename)[0]
    nmafile = header+'.nma'
    reffile = header+'.nma.ref'
    if debug:
        result = subprocess.run(sys.argv)
    else:
        result = subprocess.run(sys.argv,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    # If not doing a regeneration perform a check
    if not regenerate:
        nerrors = compareNmaFiles(nmafile, reffile)
        if nerrors > 0:
            print(title+' {} ERRORS:'.format(nerrors))
        else:
            print(title+' OK:')
        # end if
    else:
        copyfile(nmafile,reffile)
        print(title+' Regenerated:')
    # end if
    return

def runPreaderTest(title, instructions, regenerate):
    '''Run a preader test
    title is the title in the pdmake file
    instructions are the command line parameters
    regenerate is set to true if the reference file is overwritten'''
    global debug
    global rootDirectory
    global useLocal
    global debug
    if debug:
        print('runPreaderTest:',title,instructions)
    if useLocal:
        sys.argv = ['python']
        sys.argv.append(os.path.join(rootDirectory,'preader'))
    else:
        sys.argv = ['preader']
    sys.argv.extend(instructions)
    outputfile = 'command.csv'
    if regenerate:
        outputfile = 'command.ref.csv'
    with open(outputfile,'w') as fd:
        if debug:
            result = subprocess.run(sys.argv,stdout=fd)
        else:
            result = subprocess.run(sys.argv,stdout=fd,stderr=subprocess.DEVNULL)
    # end with
    # If not doing a regeneration perform a check
    if not regenerate:
        sys.argv = ['checkcsv', 'command.ref.csv','command.csv']
        with redirect(os.devnull):
            result = main_checkcsv()
        nerrors,keep_line_number,keep_word1,keep_word2,max_percentage_error = result
        if nerrors > 0:
            print(title+' ERRORS:'+"LARGEST ON LINE {} - max %error={}".format(nerrors, keep_line_number, max_percentage_error))
        else:
            print(title+' OK:'+" - max %error={}" .format(max_percentage_error))
        # end if
    else:
        print(title+' Regenerated:')
    # end if
    return

def runPDGuiTest(title, instructions, regenerate, benchmarks=False):
    '''Run a pdgui test
    title is the title in the pdmake file
    instructions are the command line parameters
    regenerate is set to true if the reference file is overwritten'''
    global viewing
    global start_time
    global debug
    global rootDirectory
    global useLocal
    if useLocal:
        sys.argv = ['python']
        sys.argv.append(os.path.join(rootDirectory,'pdgui'))
    else:
        sys.argv = ['pdgui']
    sys.argv.extend(instructions)
    if debug:
        result = subprocess.run(sys.argv)
    else:
        result = subprocess.run(sys.argv,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    if not viewing:
        # If not doing a regeneration perform a check
        if not regenerate and not benchmarks:
            sys.argv = ['checkexcel', 'results.ref.xlsx','results.xlsx']
            with redirect(os.devnull):
                result = main_checkexcel()
            nerrors,row,col,sheet,file1,value1,file2,value2,max_percentage_error = result
            if nerrors > 0:
                print(title+'{} ERRORS:'.format(nerrors)+"{}@{},{} - max %error={}".format(sheet,row,col, max_percentage_error))
            else:
                print(title+' OK:'+" - max %error={}" .format(max_percentage_error))
            # end if
        elif benchmarks:
            end_time = time.time()
            elapsed_time = end_time - start_time
            start_time = end_time
            print(title+' OK:'+" - elapsed time {:.3f}s" .format(elapsed_time))
        elif not benchmarks:
            # If we asked for a benchmarking then don't do a regenerate
            copyfile('results.xlsx','results.ref.xlsx')
            print(title+' Regenerated:')
        # end if not regenerate
    #end if not viewing
    return

def changePadding(all):
    maxlen = 0
    for d in all:
        maxlen = max(maxlen,len(d))
    if maxlen < settings['padding']:
        maxlen = settings['padding']
    return maxlen

def runTests(testlist, testType, regenerate):
    '''Run the tests given in the directories stored in testlist
       testType can be pdgui, preader, p2cif or benchmark
       regenerate causes the test data to be regenerated'''
    global start_time
    global debug
    mychdir(os.path.join(rootDirectory,'Examples'))
    homedir = os.getcwd()
    convertTestType={}
    convertTestType['pdgui'] = 'PDGui tests'
    convertTestType['vibanalysis'] = 'VibAnalysis tests'
    convertTestType['preader'] = 'PReader tests'
    convertTestType['p2cif'] = 'P2Cif tests'
    convertTestType['benchmarks'] = 'BenchMarks'
    print('')
    print('--------------------------------------------------')
    print(convertTestType[testType],'starting' )
    start_time = time.time()
    test_start_time = start_time
    if debug:
        print('RunTests: homedir=',homedir)
    pdmakefile = 'command.pdmake'
    for directory in testlist:
        if debug:
            print('RunTest: directory=',directory)
        if not os.path.isdir(directory):
            print('Error: command needs to be executed in the PDielec home directory')
            print('       current directory is {}'.format(homedir))
            print('       required directory is {}'.format(directory))
            exit()
        if testType == 'benchmarks':
            runPdMakefile(directory,pdmakefile,regenerate,benchmarks=True)
        else:
            runPdMakefile(directory,pdmakefile,regenerate)
    # end for
    elapsed_time = time.time() - test_start_time
    print('--------------------------------------------------')
    print(convertTestType[testType],'completed in {:.3f}s'.format(elapsed_time))

def readPdMakefile(directory,filename):
    # Return a dictionary of instructions
    # The first word is the dictionary key, the rest is its value
    global settings
    global debug
    instructions = {}
    if debug:
        print('ReadPdMakefile:',directory,filename)
    with open(filename,'r') as fd:
        line = fd.readline()[:-1]
        if settings['title'] == 'title':
            title = line.ljust(settings['padding'])
        else:
            if len(directory)>1:
                title = directory.ljust(settings['padding'])
            else:
                title = os.getcwd().ljust(settings['padding'])
        line = fd.readline()[:-1]
        while line:
            if len(line)>0 and not line.startswith('#'):
                full_line = line.replace('\\','')
                while line.endswith('\\'):
                    line = fd.readline()[:-1]
                    full_line += line.replace('\\','')
                splits = full_line.split()
                key = splits[0]
                value = splits[1:]
                instructions[key] = value
            # end if
            line = fd.readline()[:-1]
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
            print('mychdir: Error directory does not exist',directory)
            print('mychdir: Current directory is          ',os.getcwd())
        # end if
    # end if
    return

def runPdMakefile(directory,pdmakefile,regenerate,benchmarks=False):
    '''Run specific pdMakefile'''
    global debug
    homedir = os.getcwd()
    if debug:
        print('RunPdMakefile: directory =',directory)
        print('RunPdMakefile: pdmakefile =',pdmakefile)
        print('RunPdMakefile: cwd =',homedir)
    #directory_,filename = os.path.split(pdmakefile)
    #if debug:
    #    print('RunPdMakefile: directory_ =',directory_)
    #    print('RunPdMakefile: filename =',filename)
    #mychdir(directory_)
    mychdir(directory)
    title,instructions = readPdMakefile(directory,pdmakefile)
    if debug:
        print('RunPdMakefile: title =',title)
        print('RunPdMakefile: instructions =',instructions)
    for key in instructions:
        parameters = instructions[key]
        if key == 'preader':
            runPreaderTest(title,parameters,regenerate)
        elif key == 'vibanalysis':
            runVibAnalysis(title,parameters,regenerate)
        elif key == 'p2cif':
            runP2CifTest(title,parameters,regenerate)
        elif key == 'pdgui':
            if viewing:
                parameters.extend(['-script', 'script.py'])
            elif benchmarks:
                parameters.extend(['-nosplash', '-exit', '-script', 'script.py'])
            else:
                parameters.extend(['-nosplash', 'results.xlsx', '-exit', '-script', 'script.py'])
            if debug:
                parameters.append('-d')
            runPDGuiTest(title,parameters,regenerate,benchmarks=benchmarks)
    mychdir(homedir)

def runClean():
    if os.name == 'nt':
        print('Unable to clean installation in Windows')
        return
    print('Cleaning old results from the Examples directory')
    os.chdir(rootDirectory)
    subprocess.run('find . -name results.xlsx -exec rm -f {} \;',shell=True)
    subprocess.run('find . -name results.csv -exec rm -f {} \;',shell=True)
    subprocess.run('find . -name all.cif -exec rm -f {} \;',shell=True)
    subprocess.run('find . -name \*.nma -exec rm -f {} \;',shell=True)
    print('Cleaning complete')

def runPyInstaller():
    if os.name == 'nt':
        print('Unable to create pyInstaller installation in Windows')
        return
    print('Creating pyinstaller files in dis/pdgui         ')
    os.chdir(rootDirectory)
    subprocess.run('pyinstaller pdgui.spec -y',shell=True)
    subprocess.run('cp -r dist/pdgui/PyQt5/Qt/plugins/platforms dist/pdgui',shell=True)

def runPyPi():
    if os.name == 'nt':
        print('Unable to create PyPi installation in Windows')
        return
    print('Creating PyPi distribution files')
    os.chdir(rootDirectory)
    subprocess.run('rm -rf build dist PDielec.egg-info; python setup.py sdist bdist_wheel',shell=True)

def testForRootDirectory(path):
    test = os.path.join(path,'Examples')
    if not os.path.isdir(test):
        return False
#    test = os.path.join(path,'PDielec')
#    if not os.path.isdir(test):
#        return False
#    test = os.path.join(path,'Sphinx')
#    if not os.path.isdir(test):
#        return False
    return True

def findRootDirectory(start):
    '''Find the root directory starting from start'''
    lastpath = start
    path = os.path.join(start,'Examples')
    while path != lastpath:
        lastpath = path
        path,tail=os.path.split(path)
        if testForRootDirectory(path):
            return path
    print('Unable to locate root directory')
    exit()

def main():
    global debug
    global viewing
    global settings
    global rootDirectory
    global useLocal
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
    #rootDirectory,command = os.path.split(command)
    #if rootDirectory == '':
    #    rootDirectory = originalDirectory
    rootDirectory = originalDirectory
    rootDirectory = findRootDirectory(rootDirectory)
    actions = []
    pdmakefiles = []
    # Loop over the tokens on the command line
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token == 'test' or token == 'tests':
            actions.append('test all')
        elif token == 'test-preader':
            actions.append('test preader')
        elif token == 'test-p2cif':
            actions.append('test p2cif')
        elif token == 'test-pdgui':
            actions.append('test pdgui')
        elif token == 'test-vibanalysis':
            actions.append('test vibanalysis')
        elif token == 'benchmarks':
            actions.append('run benchmarks')
        elif token == '-root' or token == '--root':
            itoken += 1
            rootDirectory = tokens[itoken]
        elif token == '-debug' or token == '--debug' or token == '-d':
            debug = True
        elif token == '-usesystem':
            useLocal = False
        elif token == '-view' or token == '--view' or token == '-v':
            viewing = True
        elif token == '-regenerate' or token == '--regenerate':
            regenerate = True
        elif token == '-padding' or token == '--padding':
            itoken += 1
            settings['padding'] = int(tokens[itoken])
        elif token == '-directory' or token == '--directory':
            itoken += 1
            settings['title'] = 'directory'
        elif token == 'pypi':
            actions.append('run pypi')
#        elif token == 'pyinstaller':
#            actions.append('run pyinstaller')
        elif token == 'clean':
            actions.append('run clean')
        elif token == 'install':
            actions.append('install')
        elif token == 'scripts' or token == '-scripts':
            itoken += 1
            scriptsDirectory = tokens[itoken]
        elif token.endswith('.pdmake'):
            pdmakefiles.append(token)
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
    #
    # Check to see if the required executables are present
    # If the are then use them
    # 
    if useLocal:
        if checkLocalExecutables():
            useLocal = True
        else:
            useLocal = False
    # Change directory to the Examples directory
    if not os.path.isdir('Examples'):
        print('Error: command needs to be executed in the PDielec home directory')
        exit()
    #
    # Handle any pdmake files there might be
    #
    for pdmakefile in pdmakefiles:
        directory = originalDirectory
        runPdMakefile(directory,pdmakefile,regenerate)
    #
    # Now loop over the actions
    #
    mychdir('Examples')
    if settings['title'] != 'title':
        test_all = []
        test_all.extend(test_p2cif)
        test_all.extend(test_preader)
        test_all.extend(test_pdgui)
        test_all.extend(test_vibanalysis)
        settings['padding'] = changePadding(test_all)
    #
    for action in actions:
        if action == 'test all':
            runTests(test_p2cif      ,'p2cif'      ,regenerate)
            runTests(test_preader    ,'preader'    ,regenerate)
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
