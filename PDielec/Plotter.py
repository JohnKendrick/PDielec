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
"""Do some pretty printing for the dielectric functions."""
import sys


def print_ints(title, ints, no_per_line=8, format="{:9d}", file=sys.stdout, separator=" "):
    """Print ints data.

    Parameters
    ----------
    title : str
       A title for the print
    ints : a list of ints
       A list of integers to be printed
    no_per_line : int, optional
       The number of  integers per line of output
    format : str, optional 
        a format descriptor (defaults to 9d)
    file : int, optional 
        a file descriptor (defaults to stdout)
    separator : str
        the separator to use between ints (defaults to space)

    Returns
    -------
    None

    """
    #
    # Print out a list of ints prettily
    #
    len_ints = len(ints)
    if title != "":
        print(" ", file=file)
        print(title, file=file)
    nlines = int((len_ints - 1) / no_per_line) + 1
    start = 0
    for _i in range(nlines):
        end = start + no_per_line
        end = min(end, len_ints)
        print(" " + separator.join(format.format(r) for r in ints[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return

def print_strings(title, strings, no_per_line=8, format="{:9s}", file=sys.stdout, separator=" "):
    """Print strings data.

    Parameters
    ----------
    title : str
       A title for the print
    strings : a list of strings
       A list of strings to be printed
    no_per_line : int, optional
       The number of  integers per line of output
    format : str, optional 
        a format descriptor (defaults to 9s)
    file : int, optional 
        a file descriptor (defaults to stdout)
    separator : str
        the separator to use between ints (defaults to space)

    Returns
    -------
    None

    """
    #
    # Print out a list of strings prettily
    #
    len_strings = len(strings)
    if title != "":
        print(" ", file=file)
        print(title, file=file)
    nlines = int((len_strings - 1) / no_per_line) + 1
    start = 0
    for _i in range(nlines):
        end = start + no_per_line
        end = min(end, len_strings)
        print(" " + separator.join(format.format(r) for r in strings[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return


def print_reals(title, reals, no_per_line=8, format="{:9.2f}", file=sys.stdout, separator=" "):
    """Print reals data.

    Parameters
    ----------
    title : str
       A title for the print
    reals : a list of reals
       A list of reals to be printed
    no_per_line : int, optional
       The number of  integers per line of output (default is 8)
    format : str, optional 
        a format descriptor (defaults to 9.2f)
    file : int, optional 
        a file descriptor (defaults to stdout)
    separator : str
        the separator to use between ints (defaults to space)

    Returns
    -------
    None

    """
    #
    # Print out a list of reals prettily
    #
    len_reals = len(reals)
    if title != "":
        print(" ", file=file)
        print(title, file=file)
    nlines = int((len_reals - 1) / no_per_line) + 1
    start = 0
    for _i in range(nlines):
        end = start + no_per_line
        end = min(end, len_reals)
        print(" " + separator.join(format.format(r) for r in reals[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return


def print3x3(title, array, format="{:14.6f}", file=sys.stdout, separator=" "):
    """Print a 3x3 matrix.

    Parameters
    ----------
    title : str
        A title for the print
    array : tensor 3x3
        The array to be printed
    format : str, optional 
        a format descriptor (defaults to 9.2f)
    file : int, optional 
        a file descriptor (defaults to stdout)
    separator : str
        the separator to use between ints (defaults to space)

    Returns
    -------
    None

    """
    #
    # Print out a 3x3 tensor matrix
    #
    print(" ", file=file)
    if title != "":
        print(title, file=file)
    for i in range(3):
        print("      "+" ".join(format.format(p) for p in array[i]), file=file)
    # end for i
    return
