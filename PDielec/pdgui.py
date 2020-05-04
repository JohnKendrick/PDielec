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
from PDielec.GUI.App  import App
from PyQt5.QtGui     import QPixmap
from PyQt5.QtWidgets import QApplication, QSplashScreen, QProgressBar
from multiprocessing import freeze_support

def main():
    app = QApplication(sys.argv)
    show_splash = True
    for token in sys.argv:
        if token == '-nosplash' or token == '--nosplash':
            show_splash = False
        elif token == '-h' or token == '-help' or token == '--help':
            print('pdgui - graphical user interface to the PDielec package')
            print('pdgui [-help] [-debug] [program] [filename] [spreadsheet] [-script scriptname] [-nosplash] [-threads] [-cpus ncpus] [-exit]')
            exit()

    if show_splash:
        dirname = os.path.dirname(os.path.realpath(sys.argv[0]))
        splashfile = os.path.join(dirname, 'PDielec/GUI/splash.png')
        pixmap = QPixmap(splashfile)
        splash = QSplashScreen(pixmap)
        progressbar = QProgressBar(splash)
        splash.show()
    else:
        progressbar = QProgressBar()
        progressbar = None
    ex = App(sys.argv, progressbar)
    ex.show()
    if show_splash:
        splash.finish(ex)
    sys.exit(app.exec_())

if __name__ == '__main__':
    freeze_support()
    main(sys)
