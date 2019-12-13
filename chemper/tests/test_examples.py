"""
test_examples.py

This script scrapes README files and the example folder for
python examples and then runs those to make sure they are functional.

These follow the tests of the same name in the openforcefield toolkit,
which are on GitHub here:
https://github.com/openforcefield/openforcefield/blob/master/openforcefield/utils/utils.py
"""


import os
import glob
import re
import subprocess
import textwrap
import pytest

#======================================================================
# TEST UTILITIES
#======================================================================


ROOT_DIR_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')


def run_script_file(file_path):
    """Run through the shell a python script."""
    cmd = ['python', file_path]
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        raise Exception('Example {file_path} failed'.format(file_path=file_path))


def run_script_str(script_str):
    """
    Execute a Python string through the shell in a temporary directory.

    Parameters
    ----------
    script_str : str
                python script to be tested
    """
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_file_path = os.path.join(tmp_dir, 'temp.py')
        # Create temporary python script.
        with open(temp_file_path, 'w') as f:
            f.write(script_str)
        # Run the Python script.
        cmd = ['python', temp_file_path]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            script_str = textwrap.indent(script_str, '    ')
            raise Exception('The following script failed:\n'+script_str)


def find_python_examples(dir_name='examples', sub_dirs=0):
    """Find all examples in the specified directory.
    That is look for *.py files to be tested.

    Parameters
    ----------
    dir_name : str
              directory path from root directory with python files
    sub_dirs : idx
        How many sub-directories should we look for python files in?

    Returns
    -------
    example_file_paths : List[str]
        List of python scripts to execute.
    """
    examples_dir_path = os.path.join(ROOT_DIR_PATH, dir_name)

    example_file_paths = []
    # check provided directory
    test_pys = glob.glob(os.path.join(examples_dir_path, '*.py'))

    # look in subdirectories
    current_dirs = examples_dir_path
    for idx in range(sub_dirs):
        current_dirs = os.path.join(current_dirs, '*')
        test_pys += glob.glob(os.path.join(current_dirs, '*.py'))

    for example_file_path in test_pys:
        example_file_path = os.path.relpath(example_file_path)
        example_file_paths.append(example_file_path)
    return example_file_paths


def scrape_md_file(md_path):
    """
    Yield the Python scripts and URLs in the md_file in path.

    Parameters
    ----------
    md_path : str
             path to md file to scrape

    Returns
    -------
    python_examples : List[str]
        The list of Python scripts included in the provided file.
    urls :
    """
    # check there is a README in that folder
    if not os.path.isfile(md_path):
        return [], []

    with open(md_path, 'r') as f:
        readme_content = f.read()

    pythons = re.findall('```python(.*?)```', readme_content, flags=re.DOTALL)
    urls = re.findall('http[s]?://(?:[0-9a-zA-Z]|[-/.%:_])+', readme_content)
    return pythons, urls

#======================================================================
# ACTUAL TESTS
#======================================================================


# find ALL readme examples:
py_examples = list()
url_examples = list()
for dir, sub_dirs, files in os.walk(ROOT_DIR_PATH):
    for f in files:
        if '.md' in f:
            pys, us = scrape_md_file(os.path.join(dir, f))
            py_examples += pys
            url_examples += us


@pytest.mark.parametrize('py_str', py_examples)
def test_py_examples(py_str):
    """
    Test python example
    """
    run_script_str(py_str)


# TODO: figure out how to make these tests pass
# @pytest.mark.parametrize('link', url_examples)
# def test_urls(link):
#     """
#     Test the URLs from files run successfully
#     """
#     from urllib.request import Request, urlopen
#     # Some websites do not accept requests that don't specify the
#     # client and the type of accepted documents so we add fake info
#     # to avoid the response being an error.
#     headers = {'User-Agent':'Mozilla',
#                'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',}
#     request = Request(link, headers=headers)
#     urlopen(request)
#


py_scripts = find_python_examples(sub_dirs=2)
py_scripts += find_python_examples('docs', 2)


@pytest.mark.parametrize('py_script', py_scripts)
def test_script_examples(py_script):
    """
    Test that *.py files from examples run without errors
    """
    f = open(py_script, 'r')
    text = f.read()
    f.close()

    run_script_str(text)
    #run_script_file(py_script)
