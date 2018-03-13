#!/usr/bin/env python
#
# Copyright 2015 John Kendrick
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License
# along with this program, if not see https://opensource.org/licenses/MIT
#
"""PDgui driver program to calculate dielectric response at infrared and THz frequencies"""
from __future__ import print_function
import os
import sys
import time
from Python.GUI.App  import App
from PyQt5.QtGui     import QPixmap
from PyQt5.QtWidgets import QApplication, QSplashScreen, QProgressBar
from multiprocessing import freeze_support

def main(sys):
    app = QApplication(sys.argv)
    dirname = os.path.dirname(os.path.realpath(sys.argv[0]))
    splashfile = os.path.join(dirname, 'Python/GUI/splash.png')
    pixmap = QPixmap(splashfile)
    splash = QSplashScreen(pixmap)
    progressbar = QProgressBar(splash)
    splash.show()
    ex = App(sys.argv, progressbar)
    ex.show()
    splash.finish(ex)
    sys.exit(app.exec_())

if __name__ == '__main__':
    freeze_support()
    main(sys)
