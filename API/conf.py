# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PDielec API'
copyright = '2024, John Kendrick & Andrew Burnett'
author = 'John Kendrick & Andrew Burnett'
release = '8.1.0'

import os
import sys
sys.path.insert(0, os.path.abspath('../PDielec'))  # Source code dir relative to this file

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'autoapi.extension',
]

rst_prolog = """
.. role:: summarylabel
"""

html_css_files = [
    "css/custom.css",
]

def contains(seq, item):
    return item in seq

def prepare_jinja_env(jinja_env) -> None:
    jinja_env.tests["contains"] = contains
autoapi_prepare_jinja_env = prepare_jinja_env
autoapi_member_order = "groupwise"
autoapi_dirs = ['../PDielec']
autoapi_type = "python"
autoapi_template_dir="_template/autoapi"
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autodoc_typehints = "signature"
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'Python'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'cloud'
html_theme = 'furo'
html_static_path = ['_static']
