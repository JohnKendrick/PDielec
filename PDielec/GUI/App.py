import sys
import os.path
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtCore    import QCoreApplication
from PDielec.GUI.NoteBook import NoteBook
import PDielec.Utilities  as Utilities
from PDielec.Utilities import Debug
import PDielec.__init__
version = PDielec.__init__.__version__

class App(QMainWindow):

    def __init__(self, args, progressbar):
        super().__init__()
        program = ''
        filename = ''
        spreadsheet_name = ''
        self.program_exit = False
        self.debug = False
        self.scripting = False
        nosplash = False
        self.scriptname = ''
        # Manage options
        tokens = args[1:]
        ntokens = len(tokens)
        itoken = 0
        ncpus = 0
        threading = False
        default_scenario = 'powder'
        # Look at the environment to see if the number of cpus is specified
        token = os.getenv('PDIELEC_NUM_PROCESSORS')
        if token is not None:
            ncpus = int(token)
        else:
            ncpus = 0
        # Look at the environment to see if the number of threading is specified
        token = os.getenv('PDIELEC_THREADING')
        if token is not None:
            if token == 'FALSE' or token == 'false':
                threading = False
            else:
                threading = True
        else:
            threading = False
        parameters = []
        # Process any instructions on the input line
        while itoken < ntokens:
            token = tokens[itoken]
            if token == '-d' or token == '-debug' or token == '--debug':
                self.debug = True
            elif token == '-nosplash' or token == '--nosplash':
                nosplash = True
            elif token == '-exit' or token == '--exit' or token == '-quit' or token == '--quit':
                self.program_exit = True
            elif token == '-script' or token == '--script':
                itoken += 1
                self.scripting = True
                self.scriptname = tokens[itoken]
            elif token == '-spreadsheet' or token == '--spreadsheet' or token == '-xls' or token == '--xls':
                itoken += 1
                spreadsheet_name = tokens[itoken]
            elif token == '-program' or token == '--program':
                itoken += 1
                program = tokens[itoken]
            elif token == '-threading' or token == '--threading' or token == '-threads' or token == '--threads':
                threading = True
            elif token == '-cpus' or token == '--cpus':
                itoken += 1
                ncpus = int(tokens[itoken])
            elif token == '-scenario' or token == '--scenario':
                itoken += 1
                default_scenario = tokens[itoken]
                if default_scenario != 'powder' and default_scenario != 'crystal':
                    print('Error in default scenario: must be \'powder\' or \'crystal\'')
                    self.print_usage()
                    exit()
            elif token[0:0] == '-' or token == '-h' or token == '-help' or token == '--help':
                self.print_usage()
                exit()
            else:
                parameters.append(token)
            itoken += 1
        #
        # Lets see how many command line parameters there were
        #
        if len(parameters) == 0:
            pass
        elif len(parameters) == 1:
            filename = parameters[0]
            program = Utilities.find_program_from_name(filename)
        elif len(parameters) == 2:
            program = parameters[0]
            filename = parameters[1]
        elif len(parameters) == 3:
            program = parameters[0]
            filename = parameters[1]
            spreadsheet_name = parameters[2]
        else:
            self.print_usage()
            exit()
        #
        # Continue
        #
        self.version = version
        self.title = 'PDGui - Using PDielec library {} '.format(self.version)
        self.left = 10
        self.top = 30
        self.width = 500
        self.height = 800
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        QCoreApplication.processEvents()
        global debugger
        debugger = Debug(self.debug,'App:')
        debugger.print('Start:: Initialising')
        debugger.print('About to open the notebook')
        debugger.print('Program is', program)
        debugger.print('Filename is', filename)
        debugger.print('Spreadsheet is', spreadsheet_name)
        debugger.print('Script is', self.scriptname)
        debugger.print('The default scenario is', default_scenario)
        debugger.print('No. of cpus is', ncpus)
        debugger.print('Threading is', threading)
        self.notebook = NoteBook(self, program, filename, spreadsheet_name, scripting=self.scripting, progressbar=progressbar, debug=self.debug, ncpus=ncpus, threading=threading, default_scenario=default_scenario)
        debugger.print('About to call setCentralWidget')
        self.setCentralWidget(self.notebook)
        debugger.print('Finished call setCentralWidget')
        if self.scripting:
            debugger.print('Processing script',self.scriptname)
            self.readScript(self.scriptname,spreadsheet_name=spreadsheet_name)
        if self.program_exit:
            if spreadsheet_name != '':
                debugger.print('Writing spreadsheeet on exit',spreadsheet_name)
                self.notebook.writeSpreadsheet()
                self.notebook.spreadsheet.close()
            debugger.print('Exiting with sys.exit call')
            self.notebook.pool.close()
            self.notebook.pool.join()
            sys.exit()
        #self.show()
        debugger.print('Finished:: Initialising')

    def print_usage(self):
        print('pdgui - graphical user interface to the PDielec package')
        print('pdgui [program] filename [spreadsheet] [options]')
        print('     program      The name of the program which created the outputfile')
        print('                  Should be one of; vasp, phonopy, gulp, castep, abinit or qe')
        print('                  If the program is not given a best guess is made from the output filename')
        print('    filename      The name of the output file')
        print(' spreadsheet file The optional name of a spreadsheet (file must end with .xlsx')
        print('                  If this option is used program, filename must also be specified')
        print('   -scenario type Change the default scenario to \"type\"; either \"powder\" to \"crystal\"')
        print('-spreadsheet file An alternative way of specifying the spread sheet')
        print('    -program      An alternative way of specifying the program')
        print('     -script file The initial commands are read from a script file')
        print('   -nosplash      No splash screen is presented (useful for batch running)')
        print('  -threading      Uses threads rather than multiprocessing')
        print('       -cpus 0    Specify the number of processors or tasks (0 means all available')
        print('    -version      Print the version of the code')
        print('       -exit      Exit the program after executing any script')
        print('       -help      Prints out this help information')
        print('      -debug      Switches on debugging information')
        return

    def setMyWindowTitle(self,title):
        self.title = 'PDGui {}  - '.format(self.version) + title
        self.setWindowTitle(self.title)
        return

    def readScript(self,scriptname,spreadsheet_name=''):
        debugger.print('Start:: readScript')
        self.notebook.scripting = True
        directory = os.path.dirname(scriptname)
        # chdir to the directory that the script is in
        if directory != '':
            os.chdir(directory)
        scriptname = os.path.basename(scriptname)
        # If a script is used there are no prompts for overwriting files etc.
        self.notebook.overwriting = True
        with open(scriptname,'r') as fd:
            exec(fd.read())
        debugger.print('readScript finished reading script')
        self.notebook.scripting = False
        debugger.print('readScript notebook scripting set to False')
        if not self.program_exit:
            self.notebook.overwriting = False
            debugger.print('readScript notebook overwriting set to False')
        # The command line excel file overrides that in the script
        if spreadsheet_name != '':
            debugger.print('readScript overwriting spread sheet name:',spreadsheet_name)
            self.notebook.mainTab.settings['Excel file name'] = spreadsheet_name
        debugger.print('readScript notebook refresh')
        self.notebook.refresh()
        QCoreApplication.processEvents()
        debugger.print('Finished:: readScript')

    def closeEvent(self, event):
        # Make sure any spread sheet is closed
        debugger.print('Close event has been captured')
        self.notebook.pool.close()
        self.notebook.pool.join()
        super(App, self).closeEvent(event)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App(sys.argv)
    ex.show()
    self.notebook.pool.close()
    self.notebook.pool.join()
    sys.exit(app.exec_())
