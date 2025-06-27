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
"""App Module."""

import os.path
import sys

from qtpy.QtCore import QCoreApplication
from qtpy.QtWidgets import QApplication, QMainWindow

import PDielec.__init__
from PDielec import Utilities
from PDielec.Calculator import set_no_of_threads
from PDielec.GUI.NoteBook import NoteBook
from PDielec.Utilities import Debug

version = PDielec.__init__.__version__


class App(QMainWindow):
    """A class representing the main application window.

    This class initializes the main application window with various configurations based on command line arguments and environment variables. It includes functionalities to read scripts, handle command line inputs, set up multiprocessing or threading as needed, and manage application events.

    Parameters
    ----------
    args : list
        List of command line arguments passed to the application.
    progressbar : QProgressBar
        A progress bar object to display progress of operations in the application.

    Attributes
    ----------
    program_exit : bool
        Flag indicating whether the application should exit after executing any script.
    debug : bool
        Debugging mode flag.
    scripting : bool
        Flag indicating whether the application is running in scripting mode.
    scriptname : str
        The name of the script file to be executed if scripting mode is enabled.
    version : str
        The version of the PDielec package.
    title : str
        The window title.
    left : int
        The x-coordinate of the window's position.
    top : int
        The y-coordinate of the window's position.
    width : int
        The width of the window.
    height : int
        The height of the window.
    notebook : NoteBook
        The notebook widget that acts as the central widget of the application.

    Methods
    -------
    print_usage()
        Prints usage information for the command line interface.
    setMyWindowTitle(title)
        Sets the window title to a formatted string including the version and the given title.
    readScript(scriptname, spreadsheet_name='')
        Executes the commands from a script file and optionally sets the spreadsheet file name.
    closeEvent(event)
        Handles the close event, ensuring that multiprocessing pools are properly closed.

    """

    def __init__(self, args, progressbar):
        """Initialize the GUI with user-specified options from command-line arguments.

        Parameters
        ----------
        args : list
            Command-line arguments passed to the program.
        progressbar : ProgressBar or similar
            Progress bar object to display progress visually in the GUI.

        Raises
        ------
        SystemExit
            If command-line arguments specify an exit condition or if invalid options are
            provided that prevent the program from running correctly.

        Notes
        -----
        This constructor processes command-line arguments to configure various aspects
        of the GUI, such as whether debugging is enabled, the number of CPUs to use,
        whether to enable threading, and specifies files and scenarios to be used. It
        initiates the GUI with these settings and handles any requests for help or
        program exit specified via command-line arguments.

        Examples
        --------
        Command-line usage could involve arguments like:

        - To specify the program should not show splash: `--nosplash`
        - To exit after processing: `--exit`
        - To specify a script to run: `--script my_script.py`
        - To set the number of CPUs used: `--cpus 4`
        - For help: `--help`

        """
        super().__init__()
        program = ""
        filename = ""
        spreadsheet_name = ""
        self.program_exit = False
        self.debug = False
        self.scripting = False
        self.scriptname = ""
        # Manage options
        tokens = args[1:]
        ntokens = len(tokens)
        itoken = 0
        ncpus = 0
        threading = False
        default_scenario = "powder"
        # Look at the environment to see if the number of cpus is specified
        token = os.getenv("PDIELEC_NUM_PROCESSORS")
        ncpus = 0 if token is None else int(token)
        # Look at the environment to see if the number of threads is specified
        token = os.getenv("PDIELEC_NUM_THREADS")
        nthreads = 1 if token is None else int(token)
        # Look at the environment to see if the number of threading is specified
        token = os.getenv("PDIELEC_THREADING")
        if token is None:
            threading = False
        elif token in ("True", "true", "TRUE"):
            threading = True
        elif token in ("False", "false", "FALSE"):
            threading = False
        else:
            threading = False
        parameters = []
        program_has_been_specified = False
        # Process any instructions on the input line
        while itoken < ntokens:
            token = tokens[itoken]
            if token in ("-d", "-debug", "--debug"):
                self.debug = True
            elif token in ("-nosplash", "--nosplash"):
                pass
            elif token in ("-exit", "--exit", "-quit", "--quit"):
                self.program_exit = True
            elif token in ("-script", "--script"):
                itoken += 1
                self.scripting = True
                self.scriptname = tokens[itoken]
            elif token in ("-spreadsheet", "--spreadsheet", "-xls", "--xls"):
                itoken += 1
                spreadsheet_name = tokens[itoken]
            elif token in ("-program", "--program"):
                itoken += 1
                program = tokens[itoken]
                program_has_been_specified = True
            elif token in ("-threading", "--threading", "-threads", "--threads"):
                threading = True
            elif token in ("-threads", "--threads"):
                itoken += 1
                nthreads = int(tokens[itoken])
            elif token in ("-cpus", "--cpus"):
                itoken += 1
                ncpus = int(tokens[itoken])
            elif token in ("-scenario", "--scenario"):
                itoken += 1
                default_scenario = tokens[itoken]
                if default_scenario not in ("powder", "crystal"):
                    print("Error in default scenario: must be 'powder' or 'crystal'")
                    self.print_usage()
                    sys.exit()
            elif token.startswith("-"):
                self.print_usage()
                sys.exit()
            else:
                parameters.append(token)
            itoken += 1
        #
        # Lets see how many command line parameters there were
        #
        if len(parameters) == 0:
            pass
        elif len(parameters) == 1:
            filename = parameters[0]
            if not program_has_been_specified:
                program = Utilities.find_program_from_name(filename)
        elif len(parameters) == 2:
            if program_has_been_specified:
                print("Warning: program has been specified twice")
            program = parameters[0]
            filename = parameters[1]
        elif len(parameters) == 3:
            if program_has_been_specified:
                print("Warning: program has been specified twice")
            program = parameters[0]
            filename = parameters[1]
            spreadsheet_name = parameters[2]
        else:
            self.print_usage()
            sys.exit()
        #
        # Continue
        #
        self.version = version
        self.title = f"PDGui - Using PDielec library {self.version} "
        self.left = 10
        self.top = 30
        self.width = 900
        self.height = 800
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        QCoreApplication.processEvents()
        global debugger
        debugger = Debug(self.debug, "App:")
        debugger.print("Start:: Initialising")
        debugger.print("About to open the notebook")
        debugger.print("Program is", program)
        debugger.print("Filename is", filename)
        debugger.print("Spreadsheet is", spreadsheet_name)
        debugger.print("Script is", self.scriptname)
        debugger.print("The default scenario is", default_scenario)
        debugger.print("No. of cpus is", ncpus)
        debugger.print("Threading is", threading)
        # Set the number of threads before NUMPY is loaded
        if threading:
            # Threading is used instead of cpu multiprocessing
            # If threading is true the number of threads is the same as the number of processors
            set_no_of_threads(ncpus)
        else:
            # Threading is False, so spawning separate processes, but each one will still have nthreads threads
            set_no_of_threads(nthreads)
        self.notebook = NoteBook(
            self,
            program,
            filename,
            spreadsheet_name,
            scripting=self.scripting,
            progressbar=progressbar,
            debug=self.debug,
            ncpus=ncpus,
            threading=threading,
            default_scenario=default_scenario,
        )
        debugger.print("About to call setCentralWidget")
        self.setCentralWidget(self.notebook)
        debugger.print("Finished call setCentralWidget")
        if self.scripting:
            debugger.print("Processing script", self.scriptname)
            self.readScript(self.scriptname, spreadsheet_name=spreadsheet_name)
        if self.program_exit:
            if spreadsheet_name != "":
                debugger.print("Writing spreadsheeet on exit", spreadsheet_name)
                self.notebook.writeSpreadsheet()
                self.notebook.spreadsheet.close()
            debugger.print("Exiting with sys.exit call")
            self.notebook.pool.close()
            self.notebook.pool.join()
            self.notebook.pool = None
        debugger.print("Finished:: Initialising")
        return

    def print_usage(self):
        """Print usage information for the PDielec package's graphical user interface.

        This function prints the general usage information, available commands, and options for the graphical user interface to the PDielec package.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        The usage information includes:

        - `program`: The name of the program which created the output file. Supported programs are 'vasp', 'phonopy', 'gulp', 'castep', 'abinit', and 'qe'. The program is guessed from the filename if not specified.
        - `filename`: The name of the output file.
        - `spreadsheet file`: An optional name of a spreadsheet file (must end with .xlsx). If provided, both program and filename must be specified.
        - `-scenario type`: Changes the default scenario to "type", which can be either "powder" or "crystal".
        - `-spreadsheet file`: An alternative way to specify the spreadsheet file.
        - `-program`: An alternative way to specify the program.
        - `-script file`: Specifies that initial commands are read from a script file.
        - `-nosplash`: No splash screen is presented, which is useful for batch running.
        - `-threading`: Use threads instead of multiprocessing.
        - `-threads 1`: Specify the number of threads to use for each cpu
        - `-cpus 0`: Specify the number of processors or tasks; 0 uses all available.
        - `-version`: Prints the version of the code.
        - `-exit`: Exit the program after executing any script.
        - `-help`: Prints out help information.
        - `-debug`: Switches on debugging information.

        """
        print("pdgui - graphical user interface to the PDielec package")
        print("pdgui [program] filename [spreadsheet] [options]")
        print("     program      The name of the program which created the outputfile")
        print("                  Should be one of; vasp, phonopy, gulp, castep, abinit or qe")
        print("                  If the program is not given a best guess is made from the output filename")
        print("    filename      The name of the output file")
        print(" spreadsheet file The optional name of a spreadsheet (file must end with .xlsx")
        print("                  If this option is used program, filename must also be specified")
        print('   -scenario type Change the default scenario to "type"; either "powder" to "crystal"')
        print("-spreadsheet file An alternative way of specifying the spread sheet")
        print("    -program      Specify the program used to generate the output")
        print("     -script file The initial commands are read from a script file")
        print("   -nosplash      No splash screen is presented (useful for batch running)")
        print("  -threading      Uses threads rather than multiprocessing")
        print("  -threads   1    Uses 1 thread in the BLAS ")
        print("       -cpus 0    Specify the number of processors or tasks (0 means all available")
        print("    -version      Print the version of the code")
        print("       -exit      Exit the program after executing any script")
        print("       -help      Prints out this help information")
        print("      -debug      Switches on debugging information")
        return

    def setMyWindowTitle(self, title):
        """Set the window title with the provided title appended to the PDGui version.

        Parameters
        ----------
        title : str
            The window title

        Returns
        -------
        None
            The title to be appended after the PDGui version.

        Notes
        -----
        This function modifies the window title attribute of the instance and then updates the actual window title to reflect this change. The version of the PDGui is prefixed to the given title.

        """
        self.title = f"PDGui {self.version}  - " + title
        self.setWindowTitle(self.title)
        return

    def readScript(self, scriptname, spreadsheet_name=""):
        """Read and execute a script, optionally changing the working directory to the script's location and optionally setting a spreadsheet name.

        Parameters
        ----------
        scriptname : str
            The file path of the script to be executed.
        spreadsheet_name : str, optional
            The name of the spreadsheet, default is blank

        Returns
        -------
        None
            The name of the spreadsheet to set in the notebook settings, by default an empty string which implies no spreadsheet name will be set.

        Notes
        -----
        This function changes the current working directory to the directory of the script if its directory part is non-empty. It executes the script in the current Python environment using `exec()`. It also sets various notebook flags such as `scripting` and `overwriting`.

        After executing the script, it potentially updates the spreadsheet name in the notebook's mainTab settings if a non-empty `spreadsheet_name` is provided. It refreshes the notebook and processes pending Qt events with `QCoreApplication.processEvents()`.

        """
        debugger.print("Start:: readScript")
        self.notebook.scripting = True
        directory = os.path.dirname(scriptname)
        # chdir to the directory that the script is in
        if directory != "":
            os.chdir(directory)
        scriptname = os.path.basename(scriptname)
        # If a script is used there are no prompts for overwriting files etc.
        self.notebook.overwriting = True
        with open(scriptname) as fd:
            # lines = fd.readlines()
            script = fd.read()
            exec(script)
            # line_no = 0
            # for line in lines:
            #    line_no += 1
            #    debugger.print('line: ',line_no,line)
            #    exec(line)
        debugger.print("readScript finished reading script")
        self.notebook.scripting = False
        debugger.print("readScript notebook scripting set to False")
        if not self.program_exit:
            self.notebook.overwriting = False
            debugger.print("readScript notebook overwriting set to False")
        # The command line excel file overrides that in the script
        if spreadsheet_name != "":
            debugger.print("readScript overwriting spread sheet name:", spreadsheet_name)
            self.notebook.mainTab.settings["Excel file name"] = spreadsheet_name
        debugger.print("readScript notebook refresh")
        self.notebook.refresh(force=True)
        QCoreApplication.processEvents()
        debugger.print("Finished:: readScript")

    def closeEvent(self, event):
        # Make sure any spread sheet is closed
        """Handle the close event of the application.

        This method is called automatically when the close event is fired, typically when
        the user tries to close the application window. It ensures that the multiprocessing
        pool is properly terminated before the application shuts down.

        Parameters
        ----------
        event : QCloseEvent
            A close event

        Returns
        -------
        None
            The close event.

        Notes
        -----
        This function must be a method of a class that inherits from a PyQt or PySide
        widget which has a closeEvent method to override, such as QMainWindow, QDialog, etc.

        The `self.notebook.pool` is an instance of a multiprocessing pool, which needs
        to be closed and joined properly to ensure all processes are terminated cleanly before the application exits.

        The `super(App, self).closeEvent(event)` call makes sure that any close event operations defined in the
        base class (from which the current class is derived) are also executed.

        When the application is about to close, this method ensures clean termination
        of multiprocessing resources and performs any additional base class close event handling.

        """
        debugger.print("Close event has been captured")
        self.notebook.pool.close()
        self.notebook.pool.join()
        super().closeEvent(event)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = App(sys.argv)
    ex.show()
    sys.exit(app.exec_())
