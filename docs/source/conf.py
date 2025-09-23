import logging  # noqa: D100
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

project = 'moldenViz'
_copyright = '2025, Felipe Faria'
author = 'Felipe Faria'

# Import the project version
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
try:
    # Try to read version directly from file to avoid importing the full module
    version_file = ROOT / 'src' / 'moldenViz' / '__about__.py'
    version_globals = {}
    exec(version_file.read_text(), version_globals)
    __version__ = version_globals['__version__']

    release = __version__
    version = '.'.join(release.split('.')[:2])  # e.g., "0.1" from "0.1.4"
except (ImportError, FileNotFoundError, KeyError) as e:
    release = '0.0.0'  # Fallback or error
    version = '0.0'
    logger.warning(
        f'Warning: Could not read version from __about__.py: {e}\nUsing fallback version.',
    )

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',  # Core library for html generation from docstrings
    'sphinx.ext.autodoc.typehints',  # Render annotated types in descriptions
    'sphinx.ext.autosummary',  # Create neat summary tables
    'sphinx.ext.napoleon',  # Support for NumPy and Google style docstrings
    'sphinx.ext.viewcode',  # Add links to highlighted source code
    'sphinx.ext.intersphinx',  # Link to other projects' documentation
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


napoleon_use_ivar = True
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

html_sidebars = {'**': []}

html_css_files = ['custom.css']

autodoc_member_order = 'bysource'
autosummary_generate = False

autodoc_typehints = 'description'

# Mock imports for modules that require GUI or other system dependencies
autodoc_mock_imports = [
    'tkinter',
    'pyvista',
    'pyvistaqt',
    'PySide6',
]

# Make autodoc more permissive about import failures
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
    'ignore-module-all': True,
}

# Don't be strict about import errors
autodoc_inherit_docstrings = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}
