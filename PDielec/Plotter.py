#!/usr/bin/python
"""Do some plotting for the dielectric functions"""
from __future__ import print_function, division
import sys
import numpy as np


def print_ints(title, ints, no_per_line=8, format="{:9d}", file=sys.stdout, separator=" "):
    """Print ints data """
    #
    # Print out a list of ints prettily
    #
    len_ints = len(ints)
    if title != "":
        print(" ", file=file)
        print(title, file=file)
    nlines = int((len_ints - 1) / no_per_line) + 1
    start = 0
    for i in range(nlines):
        end = start + no_per_line
        if end > len_ints:
            end = len_ints
        print(" " + separator.join(format.format(r) for r in ints[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return

def print_strings(title, strings, no_per_line=8, format="{:9s}", file=sys.stdout, separator=" "):
    """Print strings data """
    #
    # Print out a list of strings prettily
    #
    len_strings = len(strings)
    if title != "":
        print(" ", file=file)
        print(title, file=file)
    nlines = int((len_strings - 1) / no_per_line) + 1
    start = 0
    for i in range(nlines):
        end = start + no_per_line
        if end > len_strings:
            end = len_strings
        print(" " + separator.join(format.format(r) for r in strings[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return


def print_reals(title, reals, no_per_line=8, format="{:9.2f}", file=sys.stdout, separator=" "):
    """Print reals data """
    #
    # Print out a list of reals prettily
    #
    len_reals = len(reals)
    if title != "":
        print(" ", file=file)
        print(title, file=file)
    nlines = int((len_reals - 1) / no_per_line) + 1
    start = 0
    for i in range(nlines):
        end = start + no_per_line
        if end > len_reals:
            end = len_reals
        print(" " + separator.join(format.format(r) for r in reals[start:end]), file=file)
        start = start + no_per_line
    # end for i
    return


def print3x3(title, array, format="{:14.6f}", file=sys.stdout, separator=" "):
    """Print 3x3 matrix"""
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
