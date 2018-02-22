import sys
import os.path
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5.QtWidgets import  QPushButton, QWidget, QAction, QTabWidget
from PyQt5.QtWidgets import  QListWidget, QComboBox, QLabel, QLineEdit
from PyQt5.QtWidgets import  QFileDialog, QCheckBox
from PyQt5.QtWidgets import  QVBoxLayout, QHBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtCore import Qt
from NoteBook import NoteBook
 
class App(QMainWindow):
 
    def __init__(self, args):
        super().__init__()
        program = ""
        filename = ""
        if len(args) >= 3:
            program = args[1]
            filename = args[2]
        self.title = 'PDielec:'
        self.left = 10
        self.top = 10
        self.width = 500
        self.height = 700
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
 
        self.notebook = NoteBook(self, program, filename)
        self.setCentralWidget(self.notebook)
 
        self.show()
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App(sys.argv)
    sys.exit(app.exec_())
