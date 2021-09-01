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
from PDielec.GUI.App  import App
from PyQt5.QtGui      import QPixmap
from PyQt5.QtWidgets  import QApplication, QSplashScreen, QProgressBar
from multiprocessing  import freeze_support
import PDielec.__init__
version = PDielec.__init__.__version__


def main():
    show_splash = True
    for token in sys.argv:
        if token == '-nosplash' or token == '--nosplash':
            show_splash = False
        elif token == '-v' or token == '-version' or token == '--version':
            print('Version ',version)
            exit()

    app = QApplication(sys.argv)

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
    main()
