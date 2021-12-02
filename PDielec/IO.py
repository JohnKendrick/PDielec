#!/usr/bin/python
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
"""Somple IO interface """
from __future__ import print_function
import sys
#
class pdielec_io:
    def __init__(self, filename, mode):
        self.filename = filename
        with open(filename,mode) as fd:
            self.lines = fd.readlines()
        self.index = 0
        self.end = len(self.lines)

    def readline(self):
        if self.index >= self.end:
            line = ''
        else:
            line = self.lines[self.index]
            self.index += 1
        return line

    def close(self):
        self.index = 0
        del(self.lines)
