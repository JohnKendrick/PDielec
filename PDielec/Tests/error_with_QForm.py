import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget
from PyQt5.QtWidgets import QFormLayout, QLabel, QPushButton
from PyQt5.QtGui import QPalette, QColor

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("My App")
        layout = QFormLayout()
        layout.addRow(QLabel('red'),Color('red'))
        layout.addRow(QLabel('green'),Color('green'))
        insideLayout = QFormLayout()
        insideLayout.addRow(QLabel('label 1'),QPushButton('Label 1'))
        insideLayout.addRow(QLabel('label 2'),QPushButton('Label 2'))
        layout.addRow(QLabel('inside'),insideLayout)
        row = layout.rowCount()-1
        layout.addRow(QLabel('blue'),Color('blue'))
        layout.removeRow(row)

        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

class Color(QWidget):
    def __init__(self, color):
        super(Color, self).__init__()
        self.setAutoFillBackground(True)

        palette = self.palette()
        palette.setColor(QPalette.Window, QColor(color))
        self.setPalette(palette)


app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()
