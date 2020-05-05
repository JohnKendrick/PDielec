import sys
import os.path
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtCore    import QCoreApplication
from PDielec.GUI.NoteBook import NoteBook

class App(QMainWindow):

    def __init__(self, args, progressbar):
        super().__init__()
        program = ''
        filename = ''
        spreadsheet = ''
        exit = False
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
        while itoken < ntokens:
            token = tokens[itoken]
            if token == '-d' or token == '-debug' or token == '--debug':
                self.debug = True
            elif token == '-nosplash' or token == '--nosplash':
                nosplash = True
            elif token == '-exit' or token == '--exit' or token == '-quit' or token == '--quit':
                exit = True
            elif token == '-script' or token == '--script':
                itoken += 1
                self.scripting = True
                self.scriptname = tokens[itoken]
            elif token == '-spreadsheet' or token == '--spreadsheet':
                itoken += 1
                spreadsheet = tokens[itoken]
            elif token == '-threading' or token == '--threading' or token == '-threads' or token == '--threads':
                threading = True
            elif token == '-cpus' or token == '--cpus':
                itoken += 1
                ncpus = int(tokens[itoken])
            elif token[0:0] == '-' or token == '-h' or token == '-help' or token == '--help':
                print('pdgui - graphical user interface to the PDielec package')
                print('pdgui [-help] [-debug] [program] [filename] [spreadsheet] [-script scriptname] [-nosplash] [-threading] [-cpus 0]')
                exit()
            elif program == '':
                program = token
            elif filename == '':
                filename = token
            elif spreadsheet == '':
                spreadsheet = token
            itoken += 1
        #
        self.title = 'PDielec GUI '
        self.left = 10
        self.top = 30
        self.width = 500
        self.height = 800
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        QCoreApplication.processEvents()
        if self.debug:
            print('About to open the notebook')
            print('Program is', program)
            print('Filename is', filename)
            print('Spreadsheet is', spreadsheet)
            print('Script is', self.scriptname)
            print('Ncpus is', ncpus)
            print('Threading is', threading)
        self.notebook = NoteBook(self, program, filename, spreadsheet, scripting=self.scripting, progressbar=progressbar, debug=self.debug, ncpus=ncpus, threading=threading)
        self.setCentralWidget(self.notebook)
        if self.scripting:
            if self.debug:
                print('Processing script',self.scriptname)
            self.readScript(self.scriptname)
        if exit:
            if self.notebook.spreadsheet is not None:
                if self.debug:
                    print('Closing spreadsheeet on exit',self.scriptname)
                self.notebook.write_spreadsheet()
                self.notebook.spreadsheet.close()
            sys.exit()
        #self.show()

    def readScript(self,scriptname):
        self.notebook.scripting = True
        with open(scriptname,'r') as fd:
            exec(fd.read())
        self.notebook.scripting = False
        self.notebook.refresh()
        QCoreApplication.processEvents()

    def closeEvent(self, event):
        # Make sure any spread sheet is closed
        self.notebook.write_spreadsheet()
        if self.debug:
            print('Close event has been captured')
        if self.notebook.spreadsheet is not None:
            self.notebook.spreadsheet.close()
        else:
            if self.debug:
                print('Spreadsheet was not set')
        super(App, self).closeEvent(event)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App(sys.argv)
    ex.show()
    sys.exit(app.exec_())
