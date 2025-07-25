import os
import sys

from funki import __version__
from funki.common import _colors

sys.path.insert(0, os.path.abspath('../../src/'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FUNKI'
author = 'Nicolàs Palacio-Escat'
copyright = '2024, %s' % author
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'myst_parser']

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'classic'
html_title = 'FUNKI documentation'
html_static_path = ['_static']
html_favicon = '../../assets/funki_favicon.ico'
html_theme_options = {
    'relbarbgcolor': _colors['blue'],
    'headbgcolor': _colors['aqua'],
}
