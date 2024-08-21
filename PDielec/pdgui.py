#!/usr/bin/env python
#
# Copyright 2024 John Kendrick & Andrew Burnett
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
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""PDgui driver program to calculate dielectric response at infrared and THz frequencies.

Process command line arguments and start a GUI application with optional splash screen.

No explicit input parameters are taken by the function. Instead, it parses command line arguments (sys.argv) to configure application behavior. The command line arguments can control displaying a splash screen, application debugging, forcibly exiting, or showing the application version.

Notes
-----
    - The -nosplash flag disables the splash screen.
    - The -v, -version, or --version flags print the program version and exit.
    - The -exit or --exit flag forces the program to exit immediately after processing arguments.
    - The -d, -debug, or --debug flags enable debug messages.
    - The function looks for a splash image file in several locations, starting with the executable's directory, and displays it if found unless the splash screen is disabled.
    - If the splash screen is shown, it also includes a progress bar.

"""
import os
import sys
from multiprocessing import freeze_support

from qtpy.QtGui import QPixmap
from qtpy.QtWidgets import QApplication, QProgressBar, QSplashScreen

import PDielec.__init__
from PDielec.GUI.App import App

version = PDielec.__init__.__version__


def main():
    """Process command line arguments and start a GUI application with optional splash screen.

    No explicit input parameters are taken by the function. Instead, it parses command line arguments (sys.argv) to configure application behavior. The command line arguments can control displaying a splash screen, application debugging, forcibly exiting, or showing the application version.

    Notes
    -----
    - The -nosplash flag disables the splash screen.
    - The -v, -version, or --version flags print the program version and exit.
    - The -exit or --exit flag forces the program to exit immediately after processing arguments.
    - The -d, -debug, or --debug flags enable debug messages.
    - The function looks for a splash image file in several locations, starting with the executable's directory, and displays it if found unless the splash screen is disabled.
    - If the splash screen is shown, it also includes a progress bar.

    Side Effects
    ------------
    - May print the application version and exit if the corresponding flags are detected.
    - May terminate the program execution before initializing the GUI if force_exit is true.
    - Prints debug messages if the debug flag is enabled.
    - Initializes and possibly displays a GUI application, potentially with a splash screen and progress bar.

    Examples
    --------
    Command line usage examples:
    - To run the application with a splash screen and debug messages:

    >>>    pdgui dft_output_file -d

    - To run the application without a splash screen:

    >>>    pdgui dft_output_file -nosplash

    - To print the version of the application and exit:

    >>>    pdgui dft_output_file --version

    - To run a script:

    >>>    pdgui -script script.py

    """    
    show_splash = True
    force_exit = False
    debug = False
    apply_default_style = True
    for token in sys.argv:
        if token in ( "-nosplash", "--nosplash" ):
            show_splash = False
        elif token in ( "-v", "-version", "--version" ):
            print("Version ",version)
            sys.exit()
        elif token in ( "-exit", "--exit" ):
            force_exit = True
        elif token in ( "-d", "-debug", "--debug" ):
            debug = True
        elif token in ( "-style" ):
            apply_default_style = False

    if apply_default_style:
        sys.argv.append("-style")
        sys.argv.append("fusion")

    app = QApplication(sys.argv)

    if show_splash:
        dirname = os.path.dirname(os.path.realpath(sys.argv[0]))
        splashfile = os.path.join(dirname, "PDielec/GUI/splash.png")
        if debug:
            print("Looking for ",splashfile)
        if not os.path.isfile(splashfile):
            splashfile = os.path.join(dirname, "../../lib/site-packages/PDielec/GUI/splash.png")
            if debug:
                print("Looking for ",splashfile)
        if not os.path.isfile(splashfile):
            for f in sys.path:
                splashfile = os.path.join(f,"PDielec/GUI/splash.png")
                if debug:
                    print("Looking for ",splashfile)
                if os.path.isfile(splashfile):
                    break
        if debug:
            print("pdgui: splashfile ",dirname, splashfile)
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
    return 

if __name__ == "__main__":
    freeze_support()
    main()
