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
    force_exit = False
    debug = False
    for token in sys.argv:
        if token == '-nosplash' or token == '--nosplash':
            show_splash = False
        elif token == '-v' or token == '-version' or token == '--version':
            print('Version ',version)
            exit()
        elif token == '-exit' or token == '--exit':
            force_exit = True
        elif token == '-d' or token == '-debug' or token == '--debug':
            debug = True


    app = QApplication(sys.argv)

    if show_splash:
        dirname = os.path.dirname(os.path.realpath(sys.argv[0]))
        splashfile = os.path.join(dirname, 'PDielec/GUI/splash.png')
        if debug:
            print('Looking for ',splashfile)
        if not os.path.isfile(splashfile):
            splashfile = os.path.join(dirname, '../../lib/site-packages/PDielec/GUI/splash.png')
            if debug:
                print('Looking for ',splashfile)
        if not os.path.isfile(splashfile):
            for f in sys.path:
                splashfile = os.path.join(f,'PDielec/GUI/splash.png')
                if debug:
                    print('Looking for ',splashfile)
                if os.path.isfile(splashfile):
                    break
        if debug:
            print('pdgui: splashfile ',dirname, splashfile)
        pixmap = QPixmap(splashfile)
        splash = QSplashScreen(pixmap)
        progressbar = QProgressBar(splash)
        splash.show()
    else:
        progressbar = QProgressBar()
        progressbar = None
    ex = App(sys.argv, progressbar)
    if not force_exit:
        ex.show()
        if show_splash:
            splash.finish(ex)
        app.exec_()

if __name__ == '__main__':
    freeze_support()
    main()
