# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os

project = 'treeprofiler'
copyright = '2024, Ziqi Deng & Jaime Huerta-Cepas'
author = 'Ziqi Deng & Jaime Huerta-Cepas'
release = '[1.2.2-beta]'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Enable markdown support

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['static']

# Update the output directory to use 'static' instead of '_static'
def setup(app):
    app.connect('build-finished', rename_static)

def rename_static(app, exception):
    if app.builder.name == 'html':
        static_dir = os.path.join(app.outdir, '_static')
        new_static_dir = os.path.join(app.outdir, 'static')
        if os.path.exists(static_dir):
            os.rename(static_dir, new_static_dir)
            # Run the post-processing script to update paths in HTML files
            os.system('python post_process_html.py')