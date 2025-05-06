# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import re
from datetime import datetime

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "MadDM"
copyright = f"{datetime.now().year}, MadDM developers"
author = "Chiara Arina & Andrew Cheek & Mattia Di Mauro & Jan Heisig & Gian Marco Lucchetti & Fabio Maltoni & Olivier Mattelaer & Daniele Massaro"

package_root_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
# load version number from setup.cfg
with open(os.path.join(package_root_path, "version"), "r") as f:
    fstring = str(f.read())
    version = re.search(r"MadDM v(?P<version>\d+\.\d+\.\d+)", fstring).group(1)

release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.githubpages",
    "sphinxcontrib.mermaid",
    "sphinx_gallery.gen_gallery"
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Sphinx Gallery
sphinx_gallery_conf = {
    "examples_dirs": os.path.join(package_root_path, "examples"),
    "gallery_dirs": "tutorial_and_examples",
    "download_all_examples": False,
    "remove_config_comments": True,
}

pygments_style = "colorful"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_logo = "_static/logo.png"
html_favicon = "_static/favicon.ico"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "logo": {"text": "MadDM"},
    "github_url": "https://github.com/maddmhep/maddm",
    "secondary_sidebar_items": ["page-toc", "sg_download_links", "sg_launcher_links"],
    "show_toc_level": 3,
}
