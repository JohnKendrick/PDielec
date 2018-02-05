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
 
    def __init__(self):
        super().__init__()
        self.title = 'PDielec:'
        self.left = 0
        self.top = 0
        self.width = 300
        self.height = 200
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
 
        self.notebook = NoteBook(self)
        self.setCentralWidget(self.notebook)
 
        self.show()
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
