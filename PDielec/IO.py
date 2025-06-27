#!/usr/bin/python
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
"""Simple IO interface.

This aims to reproduce the python readline() io method, with limited functionality but much faster.

"""
class pdielec_io:
    """A class for reading lines from a file and providing an interface similar to file objects.

    The class opens a file upon initialization, reads its lines into memory,
    and allows sequential reading of these lines through the readline method.
    It also provides a close method to reset the reading process.
    This speeds up the reading of a file using the usual python readline() method.

    Parameters
    ----------
    filename : str
        The path to the file to be read.
    mode : str
        The mode in which the file should be opened (e.g., 'r' for read).

    Attributes
    ----------
    filename : str
        The path to the file that was opened.
    lines : list of str
        The lines read from the file.
    index : int
        The current position (line) in the file being read, initialized at 0.
    end : int
        The total number of lines in the file.

    Methods
    -------
    readline()
        Reads the next line from the file. If the end of the file has been reached,
        returns an empty string.
    close()
        Resets the reading process by resetting the current line index and
        deleting all lines read from the file from memory.

    """

    def __init__(self, filename, mode):
        """Initialize the object with a file's content.

        Parameters
        ----------
        filename : str
            The path to the file to be read.
        mode : str
            The mode in which the file should be opened (e.g., 'r' for read).

        """        
        self.filename = filename
        with open(filename,mode) as fd:
            self.lines = fd.readlines()
        self.index = 0
        self.end = len(self.lines)

    def readline(self):
        """Read and return the next line from the stored lines.

        This method reads a line from the current position within an internal list of lines
        and then advances the position so that subsequent calls will return subsequent lines.
        If the end of the list is reached, it returns an empty string.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The next line in the internal list of lines, or an empty string if the end of the list is reached.

        """        
        if self.index >= self.end:
            line = ""
        else:
            line = self.lines[self.index]
            self.index += 1
        return line

    def close(self):
        """Close the current object session.

        Resets the index to 0 and deletes the lines attribute from the object.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method is designed to perform cleanup tasks, such as resetting object state and freeing resources by deleting attributes. Calling this method on an object may render it unusable for certain operations that depend on the deleted attributes.

        """        
        self.index = 0
        del(self.lines)
