#!/usr/bin/env python

from setuptools import setup, find_packages
import os

def get_files(directory):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            file_list.append(file_path)
    return file_list

VERSION = '1.2.4.1'
install_requires = ['ete4', 'selenium', 'biopython','scipy']
example_files = get_files('examples/')
test_files = get_files('tests/')
long_description = open("README.md").read()
setup(
    name='TreeProfiler',
    version=VERSION,
    # metadata for upload to PyPI
    description='TreeProfiler is command-line tool for profiling metadata table into phylogenetic tree with descriptive analysis and output visualization',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Ziqi Deng, Jaime Huerta-Cepas',
    author_email='dengziqi1234@gmail.com, jhcepas@gmail.com',
    maintainer = 'Ziqi Deng',
    maintainer_email = 'dengziqi1234@gmail.com',
    url="https://github.com/compgenomicslab/MetaTreeDrawer",
    
    #package_dir = {'treeprofiler' : '' },
    packages=find_packages(),
    package_data = { 
        'treeprofiler' : [
            'treeprofiler/*',
            'layouts/pfam2color.json',
            'layouts/smart2color.json',
            ],
        },
    #scripts=['treeprofiler.py'],
    entry_points = {
        'console_scripts': ['treeprofiler=treeprofiler.treeprofiler:main']
    },
    data_files=[
        ('examples', example_files),
        ('tests', test_files)
    ],
    install_requires=install_requires,
    keywords = "tree annotation, tree visualization, phylogeny, phylogenetics, phylogenomics",
)

