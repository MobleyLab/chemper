"""
This is a semi-standard setup.py file
It began as a copy of the skeleton setup file by MolSSI in their python template repository at:
https://github.com/MolSSI/python_template
Then was extended using the Open Force Field Initiative Toolkits setup
file as inspiration.
https://github.com/openforcefield/openforcefield
"""

import setuptools
import os
from os.path import relpath, join

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


if __name__ == "__main__":
    setuptools.setup(
        name='chemper',        # Make sure to change to match your library name
        version="0.1.0",    # you should keep track of versions
        description='A python package for automatically sampling chemical perception',     # add a description
        long_description=read('README.md'),
        author='Caitlin C. Bannan',          # add your name to author category
        author_email='bannanc@uci.edu',    # add your e-mail
        url="https://github.com/MobleyLab/chemper",             # add github URL
        license='MIT',      # should match license in your repo
        packages=setuptools.find_packages()+['tests', 'chemper/data'],
        install_requires=[
            'numpy>=1.7',
            # do you need any other libraries?
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',  # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },

        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],

        classifiers=[
            'Development Status :: 4 - Alpha',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
        package_data={'chemper': find_package_data('chemper/data', 'chemper')}
    )
