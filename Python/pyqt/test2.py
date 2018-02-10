from PyQt5 import QtWidgets, QtCore

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent = None):
        super(MainWindow, self).__init__()
        self.__setup__()

    def __setup__(self):
        self.resize(400, 400)
        tabWidget = TabWidget(self)
        self.setCentralWidget(tabWidget)

        options = Options(self)
        optionsDock = QtWidgets.QDockWidget()
        optionsDock.setWidget(options)
        optionsDock.setWindowTitle("Options")
        self.addDockWidget(QtCore.Qt.TopDockWidgetArea, optionsDock)
        tabWidget.tabWidget.currentChanged.connect(lambda index: options.spinBox_columns.
                                                   setValue(tabWidget.tabWidget.widget(index).columnCount()))

        options.spinBox_columns.valueChanged.connect(lambda value: tabWidget.tabWidget.currentWidget().
                                                     setColumnCount(value))

class Options(QtWidgets.QWidget):
    def __init__(self, parent):
        super(Options, self).__init__(parent)
        self.__setup__()

    def __setup__(self):
        self.spinBox_columns = QtWidgets.QSpinBox()
        self.spinBox_columns.setValue(1)
        self.spinBox_columns.setMinimum(1)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.spinBox_columns)
        self.setLayout(layout)

class TabWidget(QtWidgets.QWidget):
    def __init__(self, parent):
        super(TabWidget, self).__init__(parent)
        self.__setup__()

    def __setup__(self):
        self.tabWidget = QtWidgets.QTabWidget()

        for i in range(3):
            widget = QtWidgets.QTableWidget()
            widget.setColumnCount(1)
            widget.setRowCount(3)
            self.tabWidget.addTab(widget, "Column " + str(i))

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.tabWidget)
        self.setLayout(layout)

def main():
    import sys
    app = QtWidgets.QApplication(sys.argv) 
    window = MainWindow() 
    window.show() 
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
