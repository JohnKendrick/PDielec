import setuptools
import os

# Read the version number from the pyproject.toml file
version = 'Unkown'
with open('pyproject.toml') as fh:
    for line in fh:
        line = line.lower()
        line = line.replace("="," ")
        line = line.replace('"'," ")
        line = line.replace("'"," ")
        line_list = line.split()
        if len(line_list) > 1:
            if line_list[0] == "version":
              version = line_list[1]

# Read in the license file
with open("LICENSE.txt", "r") as fh:
    licence = fh.readlines()

# Create an __init_.py file for each module
init_lines = ["# Created by setup.py, whilst creating a new PyPi release\n","__version__ = \"{}\"\n".format(version)]
with open("PDielec/__init__.py", "w") as fh:
    fh.writelines(['"""'])
    fh.writelines(['PDielec: library for processing DFT calculations for infrared and terahertz spectroscopy.\n'])
    fh.writelines(['\n'])
    fh.writelines(licence)
    fh.writelines(['"""\n'])
    fh.writelines(init_lines)
with open("PDielec/GUI/__init__.py", "w") as fh:
    fh.writelines(['"""'])
    fh.writelines(['PDGui: GUI for processing DFT calculations for infrared and terahertz spectroscopy.\n'])
    fh.writelines(['\n'])
    fh.writelines(licence)
    fh.writelines(['"""\n'])
    fh.writelines(init_lines)

setuptools.setup()
