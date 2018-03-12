import sys
import os.path
from PyQt5.QtWidgets import QMainWindow, QApplication
from Python.GUI.NoteBook import NoteBook
 
class App(QMainWindow):
 
    def __init__(self, args):
        super().__init__()
        program = ''
        filename = ''
        spreadsheet = ''
        exit = False
        debug = False
        script = False
        scriptname = ''
        # Manage options
        tokens = args[1:]
        ntokens = len(tokens)
        itoken = 0
        while itoken < ntokens:
            token = tokens[itoken]
            if token == '-d' or token == '-debug' or token == '--debug':
                debug = True
            elif token == '-exit' or token == '--exit' or token == '-quit' or token == '--quit':
                exit = True
            elif token == '-script' or token == '--script':
                itoken += 1
                script = True
                scriptname = tokens[itoken]
            elif token == '-spreadsheet' or token == '--spreadsheet':
                itoken += 1
                spreadsheet = tokens[itoken]
            elif token == '-h' or token == '-help' or token == '--help':
                print('help is true')
                print('pdgui - graphical user interface to the PDielec package')
                print('pdgui [-help] [-debug] [program] [filename] [spreadsheet] [-script scriptname]')
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
 
        self.notebook = NoteBook(self, program, filename, spreadsheet, debug=debug)
        self.setCentralWidget(self.notebook)
        self.show()
        if script:
            self.readScript(scriptname)
        if exit:
            sys.exit()

    def readScript(self,scriptname):
        with open(scriptname,'r') as fd:
            exec(fd.read())
        self.notebook.refresh()

    def closeEvent(self, event):
        # Make sure any spread sheet is closed
        if self.notebook.spreadsheet is not None:
            self.notebook.spreadsheet.close()
        super(App, self).closeEvent(event)
    
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App(sys.argv)
    sys.exit(app.exec_())
