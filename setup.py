import setuptools


if __name__ == "__main__":
    setuptools.setup(
        name='chemical_perception',        # Make sure to change to match your library name
        version="0.0.0",    # you should keep track of versions
        description='A python package for automatically sampling chemical perception',     # add a description
        author='Caitlin C. Bannan',          # add your name to author category
        author_email='bannanc@uci.edu',    # add your e-mail
        url="https://github.com/MobleyLab/chemical_perception",             # add github URL
        license='MIT',      # should match license in your repo
        packages=setuptools.find_packages(),
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
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
    )
