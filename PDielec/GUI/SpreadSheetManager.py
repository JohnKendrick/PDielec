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
"""SpreadSheetMamager module."""
import numpy as np
import xlsxwriter as xlsx


class SpreadSheetManager:
    """A manager for handling operations on an Excel workbook via xlsxwriter.

    This class provides methods for managing a workbook with various tabs
    (named worksheets), writing data to these worksheets, and handling complex
    data types specifically tailored for scientific or engineering applications.

    Parameters
    ----------
    filename : str
        The name of the file to which the workbook will be saved.

    Attributes
    ----------
    workbook : Workbook
        The workbook object from xlsxwriter this manager operates on.
    closed : bool
        Boolean flag indicating whether the workbook has been closed.
    tab_names : list of str
        A list containing the names of the tabs (worksheets) to be managed.
    worksheets : dict
        A dictionary mapping tab names to worksheet objects.
    positions : dict
        A dictionary mapping tab names to their current write position (row, column).
    max_col : dict
        A dictionary mapping tab names to the maximum column index used.
    max_row : dict
        A dictionary mapping tab names to the maximum row index used.
    opened : dict
        A dictionary indicating whether a tab (worksheet) has been opened for writing.
    name : str
        The name of the currently selected work tab.

    Methods
    -------
    openWorkSheet(tab)
        Opens a new worksheet for the given tab if not already opened.
    selectWorkSheet(name)
        Selects the specified worksheet for further operations.
    writeNextRow(items, row=None, col=None, check='')
        Writes a series of items to the next row of the currently selected worksheet.
    write(row, col, item)
        Writes an item to the specified location in the currently selected worksheet.
    delete()
        Deletes all content from the currently selected worksheet.
    close()
        Closes the workbook, finalizing it for output.

    """

    def __init__(self, filename):
        # Manage a spread sheet for PDielec / PDGui
        """Initialize the instance with the given filename.

        This constructor initializes an Excel workbook with predefined tab names and settings. It prepares the workbook for data entry and analysis across various predefined categories related to powder and crystal property analysis.

        Parameters
        ----------
        filename : str
            The name of the file to which the workbook will be saved.

        Attributes
        ----------
        workbook : xlsx.Workbook
            An instance of an Excel workbook.
        closed : bool
            A flag indicating whether the workbook is closed.
        tab_names : list
            A list of strings indicating the names of the workbook tabs.
        worksheets : dict
            A dictionary to keep track of worksheet objects by tab name.
        positions : dict
            A dictionary storing the positions of the last written row and column for each tab.
        max_col : dict
            A dictionary tracking the maximum column number used in each tab.
        max_row : dict
            A dictionary tracking the maximum row number used in each tab.
        opened : dict
            A dictionary indicating whether a tab is open for editing.
        name : str
            The name of the currently active worksheet.

        Notes
        -----
        This constructor also automatically opens the worksheet named 'Main' to ready it for data entry or manipulation.

        Examples
        --------
        >>> my_excel_file = ExcelFileHandler("example.xlsx")
        >>> my_excel_file.add_data_to_sheet("Main", data)

        """        
        self.workbook = xlsx.Workbook(filename)
        self.closed = False
        self.tab_names = ["Main", "Settings", "Scenarios",
                "Powder Molar Absorption (cells)",
                "Powder Molar Absorption (atoms)",
                "Powder Molar Absorption (mols)",
                "Powder Absorption",
                "Powder Real Permittivity",
                "Powder Imaginary Permittivity",
                "Powder ATR Reflectance",
                "Analysis",
                "Crystal R_p",
                "Crystal R_s",
                "Crystal T_p",
                "Crystal T_s",
                "Crystal A_p",
                "Crystal A_s",
                "Real Crystal Permittivity",
                "Imag Crystal Permittivity" ]
        self.worksheets = {}
        # Positions points to where we write to next
        self.positions  = {}
        self.max_col    = {}
        self.max_row    = {}
        self.opened     = {}
        for tab in self.tab_names:
            self.opened[tab] = False
        self.name = "Main"
        self.openWorkSheet(self.name)

    def openWorkSheet(self,tab):
        """Open a new worksheet in a workbook.

        Parameters
        ----------
        tab : str
            The name of the worksheet to be opened.

        Returns
        -------
        None

        Notes
        -----
        If the worksheet identified by `tab` is already opened, this method does nothing. 
        Otherwise, it creates a new worksheet with the name `tab`, initializes its position, 
        maximum column, maximum row, and marks it as opened within the workbook data structures 
        handled by the instance. 

        """        
        if self.opened[tab]:
            return
        self.worksheets[tab] = self.workbook.add_worksheet(tab)
        self.positions[tab] = (0,0)
        self.max_col[tab] = 0
        self.max_row[tab] = 0
        self.opened[tab] = True

    def selectWorkSheet(self,name):
        """Select or opens a worksheet by name.

        Parameters
        ----------
        name : str
            The name of the worksheet to select or open.

        Returns
        -------
        None

        Notes
        -----
        This method selects a worksheet if it is already opened. If the worksheet is not
        currently open, it attempts to open the worksheet by calling `openWorkSheet` with
        the worksheet name.

        """        
        self.name = name
        if not self.opened[name]:
            self.openWorkSheet(self.name)

    def writeNextRow(self,items, row=None, col=None, check=""):
        """Write a sequence of items as a row in a spread sheet, starting from a specified row and column into a grid structure.

        Parameters
        ----------
        items : list
            A list of items to be written. Items can be strings, numbers, or nested lists where each sub-item is written in subsequent columns.
        row : int, optional
            The row index to start writing from. If not specified, uses the current row position of the object.
        col : int, optional
            The column index to start writing from. If not specified, uses the current column position of the object.
        check : str, optional
            A check string written to the beginning of the row specified. Defaults to an empty string if not specified.

        Returns
        -------
        None

        Notes
        -----
        - If `col` is 0, a warning message is printed.
        - The `row` and `col` are updated as items are written, and `row` is incremented after writing all items.
        - For items that are lists, each element is written in subsequent columns. This is applied recursively for nested lists.
        - `self.positions` is a dictionary holding the current positions (row, column) for different names, and `self.name` accesses the current object's name.

        """        
        oldRow,oldCol = self.positions[self.name]
        if col is None:
            col = oldCol
        if row is None:
            row = oldRow
        if col == 0:
            print("We have a problem, col is 0")
        self.write(row,0,check)
        for item in items:
            if isinstance(item,list):
                for i in item:
                    if isinstance(i, list):
                        for j in i:
                            self.write(row,col,j)
                            col += 1
                    else:
                        self.write(row,col,i)
                        col += 1
            else:
                self.write(row, col, item)
                col += 1
        row += 1

    def write(self,row,col,item):
        """Write an item to a specific position in the active worksheet, handling complex numbers.

        Parameters
        ----------
        row : int
            The row index where the item will be written.
        col : int
            The column index where the item will be written.
        item : int, float, complex, str
            The item to be written to the worksheet. If the item is a complex number,
            its real and imaginary parts are written to adjacent columns.

        Returns
        -------
        None

        Notes
        -----
        - This method updates the worksheet's max column and row based on the provided indices.
        - When writing a complex number, `col` is incremented before writing the imaginary part,
          and the updated column is used for adjusting the max column.
        - Updates the `positions` attribute with the next available position in the worksheet.

        """        
        if isinstance(item,complex) or np.iscomplexobj(item):
            self.worksheets[self.name].write(row,col,item.real)
            col += 1
            self.worksheets[self.name].write(row,col,item.imag)
            self.max_col[self.name] = max(self.max_col[self.name], col)
            self.max_row[self.name] = max(self.max_row[self.name], row)
            self.positions[self.name] = (row+1,col+1)
        else:
            self.worksheets[self.name].write(row,col,item)
            self.max_col[self.name] = max(self.max_col[self.name], col)
            self.max_row[self.name] = max(self.max_row[self.name], row)
            self.positions[self.name] = (row+1,col+1)

    def delete(self):
        """Delete all the data in the current worksheet.

        Clears all entries in the worksheet represented by `self` by setting each cell's value to an empty string.
        After clearing the data, it resets the current position, maximum column, and maximum row counters for the worksheet.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        for row in range(0,self.max_row[self.name]):
            for col in range(0, self.max_col[self.name]):
                self.worksheets[self.name].write(row,col,"")
        self.positions[self.name] = (-1,0)
        self.max_col[self.name] = 0
        self.max_row[self.name] = 0

    def close(self):
        """Close the workbook if it hasn't been closed already.

        Closes the workbook associated with an instance and sets its status as closed to prevent multiple close operations.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        if not self.closed:
            self.workbook.close()
        self.closed = True
