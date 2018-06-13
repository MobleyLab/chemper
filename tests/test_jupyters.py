"""
Tests Jupyter notebooks provided in examples
by running the code in each notebook.

The exe_scriptified_ipynb was adjusted from psi4numpy test utils:
https://github.com/psi4/psi4numpy/blob/049f0a535283eabafbbd2901bd89c8941d68d260/tests/utils.py
"""

import glob
import shutil
import tempfile
import os
import re
import pytest
from chemper.mol_toolkits import mol_toolkit

# from https://stackoverflow.com/a/31499114
def sed_inplace(filename, pattern, repl):
    """Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
    `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
    Examples
    --------
    sed_inplace('/etc/apt/sources.list', r'^\# deb', 'deb')
    """
    # For efficiency, precompile the passed regular expression.
    pattern_compiled = re.compile(pattern)

    # For portability, NamedTemporaryFile() defaults to mode "w+b" (i.e., binary
    # writing with updating). This is usually a good thing. In this case,
    # however, binary writing imposes non-trivial encoding constraints trivially
    # resolved by switching to text writing. Let's do that.
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
        with open(filename) as src_file:
            for line in src_file:
                tmp_file.write(pattern_compiled.sub(repl, line))

    # Overwrite the original file with the munged temporary file in a
    # manner preserving file attributes (e.g., permissions).
    shutil.copystat(filename, tmp_file.name)
    shutil.move(tmp_file.name, filename)


def exe_scriptified_ipynb(workspace, tdir, ipynb):
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = base_dir + '/' + tdir + '/' + ipynb + '.ipynb'
    path = workspace.workspace
    workspace.run('jupyter nbconvert --to script ' + script + ' --output-dir=' + path)
    script_py = path + '/' + ipynb + '.py'
    sed_inplace(script_py,
                r"""get_ipython\(\).magic\(u?'matplotlib inline'\)""",
                """# <<<  Jupyter magic  >>>  get_ipython().magic('matplotlib inline')\nimport matplotlib as mpl; mpl.use('Agg')""")
    sed_inplace(script_py,
                r"""get_ipython\(\).magic\(u?'matplotlib notebook'\)""",
                """# <<<  Jupyter magic  >>>  get_ipython().magic('matplotlib notebook')\nimport matplotlib as mpl; mpl.use('Agg')""")
    sed_inplace(script_py,
                r"""get_ipython\(\).magic\(u?['"]timeit """,
                """# <<<  Jupyter magic  >>>""")
        #workspace.run('cp %s %s' % (os.environ['OE_LICENSE'], path))
        #from pytest_shutil import env
        #if 'OE_LICENSE' not in os.environ:
        #    os.environ['OE_LICENSE'] = '/home/oe_liceense.txt'
        #else:
        #    print(os.environ['OE_LICENSE'])
        #env.set_env('OE_LICENSE', os.environ['OE_LICENSE'])
    workspace.run('python ' + script_py)


notebooks = ['single_mol_smirks', 'smirks_from_molecules']
@pytest.mark.parametrize('notebook_name', notebooks)
def test_example_notebooks(workspace, notebook_name):
    print(os.environ['HOME'])
    if 'openeye' in mol_toolkit.__name__:
        if os.environ['HOME'] == '/home/travis':
            oe_f = '/home/travis/oe_license.txt'
            if os.path.isfile(oe_f):
                from shutil import copyfile
                copyfile(oe_f, path+'/oe_license.txt')
                print("found ", oe_f)

    exe_scriptified_ipynb(workspace, 'examples', notebook_name)
