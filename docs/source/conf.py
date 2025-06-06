import os
import sys

project = 'moldenViz'
copyright = '2025, Felipe Faria'
author = 'Felipe Faria'

# Import the project version
sys.path.insert(0, os.path.abspath('../../src'))
try:
    from moldenViz.__about__ import __version__

    release = __version__
    version = '.'.join(release.split('.')[:2])  # e.g., "0.1" from "0.1.4"
except ImportError:
    release = '0.0.0'  # Fallback or error
    version = '0.0'
    print('Warning: Could not import __version__ from moldenViz.__about__.')
    print("Make sure '../src' is in sys.path and moldenViz/__about__.py exists.")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',  # Core library for html generation from docstrings
    'sphinx.ext.autosummary',  # Create neat summary tables
    'sphinx.ext.napoleon',  # Support for NumPy and Google style docstrings
    'sphinx.ext.viewcode',  # Add links to highlighted source code
    'sphinx.ext.intersphinx',  # Link to other projects' documentation
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_show_sourcelink = False

html_theme_options = {
    'use_edit_page_button': False,
    'show_toc_level': 2,
    'navbar_align': 'left',  # Left, right, content
    'github_url': 'https://github.com/Faria22/moldenviz',  # Add your repo URL
    'navbar_end': ['navbar-icon-links.html', 'search-field.html'],
}

autodoc_member_order = 'bysource'
autosummary_generate = True  # Uncomment if you want to auto-generate stubs

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}
