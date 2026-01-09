"""Sphinx configuration for GECKO documentation."""

project = "GECKO"
copyright = "2025, ITPcp"
author = "ITPcp"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
]

# MyST configuration
myst_enable_extensions = [
    "dollarmath",
    "colon_fence",
    "deflist",
]

# Bibliography
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

# Napoleon settings (NumPy-style docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "trimesh": ("https://trimesh.org/", None),
}

# HTML output
html_theme = "furo"
html_title = "GECKO"
html_static_path = ["_static"]

# Source settings
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"
