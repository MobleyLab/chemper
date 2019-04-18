"""
This is a semi-standard setup.py file
It began as a copy of the skeleton setup file by MolSSI in their python template repository at:
https://github.com/MolSSI/python_template
Then was extended using the Open Force Field Initiative Toolkits setup
file as inspiration.
https://github.com/openforcefield/openforcefield
Note - Like all openforcefield related projects,
this file is slowly being transitioned
to more closely resemble the MolSSI cookiecutter format:
    https://github.com/MolSSI/cookiecutter-cms
"""

import os
from os.path import relpath, join
from setuptools import setup
import versioneer

def read(fname):
    """
    Taken from openforcefield/setup.py
    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def find_package_data(data_root, package_root):
    """
    Taken from openforcefield/setup.py

    It allows non-python files to be included in the installation
    """
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files


setup(
    name = 'chemper',        # Make sure to change to match your library name
    author = 'Caitlin C. Bannan',          # add your name to author category
    author_email = 'bannanc@uci.edu',    # add your e-mail
    description = ("ChemPer"),
    license = 'MIT',      # should match license in your repo
    keywords = "chemical perception, SMARTS, SMIRKS, SMIRNOFF, Open Force Field, forcefield",
    url = "http://github.com/mobleylab/chemper",
    packages = [
        'chemper',
        'chemper/data',
        'chemper/graphs',
        'chemper/mol_toolkits',
        'chemper/tests',
        ],
    long_description = read('README.md'),
    package_data = {'chemper': find_package_data('chemper/data', 'chemper')},
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass()
)
