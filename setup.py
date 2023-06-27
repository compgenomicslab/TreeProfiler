#!/usr/bin/env python

from setuptools import setup, find_packages

VERSION = '1.0.0-beta'
install_requires = ['ete4', 'selenium', 'biopython','scipy']
setup(
    name='TreeProfiler',
    version=VERSION,
    # metadata for upload to PyPI
    description='TreeProfiler is command-line tool for profiling metadata table into phylogenetic tree with descriptive analysis and output visualization',
    author='Ziqi Deng, Jaime Huerta-Cepas',
    author_email='dengziqi1234@gmail.com, jhcepas@gmail.com',
    maintainer = 'Ziqi Deng',
    maintainer_email = 'dengziqi1234@gmail.com',
    url="https://github.com/compgenomicslab/MetaTreeDrawer",
    packages=find_packages(),
    install_requires=install_requires,
    keywords = "tree annotation, tree visualization, phylogeny, phylogenetics, phylogenomics",
    
)