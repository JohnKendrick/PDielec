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
from Python.GUI.App import App
from PyQt5.QtWidgets import QApplication

app = QApplication(sys.argv)
ex = App(sys.argv)
sys.exit(app.exec_())
