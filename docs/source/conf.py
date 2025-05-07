# Configuration file for the Sphinx documentation builder.

import os
import sys
import datetime

# -- Path setup --------------------------------------------------------------

# Add project root to sys.path to enable autodoc
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------

project = 'pysplot'
author = 'pysplot contributors'
copyright = f'{datetime.datetime.now().year}, {author}'
release = '0.1.0'  # Update as needed

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'myst_parser',  # Markdown support if needed
]

templates_path = ['_templates']
exclude_patterns = []

# Support both .rst and .md files
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']