# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import shutil
import glob

sys.path.insert(0, os.path.abspath('../../calphy/'))

#copy necessary figures
#cdir = os.getcwd()
#figs = glob.glob("../../examples/*/*.png")

#for fig in figs:
#    shutil.copy(fig, "../_static/")

def skip(app, what, name, obj, would_skip, options):
    if name in ( '__init__',):
        return False
    return would_skip
def setup(app):
    app.connect('autodoc-skip-member', skip)

#copy examples
#---------------------------------
#copy ipynb here
if os.path.exists("examples"):
    shutil.rmtree("examples")
shutil.copytree("../../examples", "examples")


# -- Project information -----------------------------------------------------

project = 'calphy'
copyright = '2021, Sarath Menon, Yury Lysogorskiy, Ralf Drautz'
author = 'Sarath Menon'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'm2r2',
    'sphinx_markdown_tables',
    'nbsphinx',
]

#html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'
#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_logo = "../_static/calphy_logo.png"
html_theme_options = {
    'logo_only' : True,
    'canonical_url' : 'https://calphy.readthedocs.io/',
}

html_extra_path = ['../_static' ]

source_suffix = ['.rst', '.md']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


#html_static_path = ['../_static']
#def setup(app):
#    #app.add_stylesheet("theme_extra.css")
#    app.add_css_file("theme_extra.css")
