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
        exit = ''
        debug = False
        # Manage options
        for arg in args[1:]:
            if arg == '-d' or arg == '-debug' or arg == '--debug':
                debug = True
            elif arg == '-h' or arg == '-help' or arg == '--help':
                print('help is true')
                print('pdgui - graphical user interface to the PDielec package')
                print('pdgui [-help] [-debug] [program] [filename] [spreadsheet]')
                exit()
            elif program == '':
                program = arg
            elif filename == '':
                filename = arg
            elif spreadsheet == '':
                spreadsheet = arg
            elif exit == '':
                exit = arg
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
        if exit == 'exit':
            sys.exit()

    def closeEvent(self, event):
        # Make sure any spread sheet is closed
        print('captured close event')
        if self.notebook.spreadsheet is not None:
            print('closing spreadsheet')
            self.notebook.spreadsheet.close()
        super(App, self).closeEvent(event)
    
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App(sys.argv)
    sys.exit(app.exec_())
