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

import setuptools
import os
from os.path import relpath, join
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


if __name__ == "__main__":

    short_description = __doc__.split("\n")
    try:
        with open("README.md", "r") as handle:
            long_description = handle.read()
    except:
        long_description = "\n".join(short_description[2:])

    setuptools.setup(
        name='chemper',        # Make sure to change to match your library name
        version="0.1.0",    # you should keep track of versions
        description='A python package for automatically sampling chemical perception',     # add a description
        long_description=read('README.md'),
        author='Caitlin C. Bannan',          # add your name to author category
        author_email='bannanc@uci.edu',    # add your e-mail
        description=short_description[0],
        long_description=long_description,
        long_description_content_type="text/markdown",
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        license='MIT',      # should match license in your repo

        packages=['chemper', 'chemper.tests'],
        # look for other package data, not just the python modules
        # this should install all data in the chemper/data/ folders
        package_data={'chemper': ["data/*"] },
        # This was previously
        # package_data={'chemper': find_package_data('chemper/data', 'chemper')}
        # TODO: permanently delete line above and relevant function if it works
        # these were optional in cookiecutter, I added for completeness
        url="https://github.com/MobleyLab/chemper",             # add github URL
        python_requires=">=3.5",
        zip_safe=True,
    )
