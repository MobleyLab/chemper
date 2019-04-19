"""
test_complex_typing.py

This script tests if chemper can generate SMIRKS
matching the complexity found in the smirnoff99Frosst force field.
This is done by typing a set of molecules in chemper/data/molecules
with the SMIRKS stored at chemper/data/smarts_files.
Then the clusters based on those SMIRKS are used to generate a
hierarchical list of SMIRKS patterns using the SMIRKSifier.
"""

from copy import deepcopy
import pytest
from itertools import product
from chemper.smirksify import SMIRKSifier
from chemper.chemper_utils import get_full_path, create_tuples_for_clusters, get_typed_molecules
from chemper.mol_toolkits.mol_toolkit import mols_from_mol2
from chemper.chemper_utils import get_data_path


def parse_smarts_file(smarts_file_name):
    """
    Parameters
    ----------
    smarts_file_name: str
        smarts file located in chemper/data/smarts_files/

    Returns
    -------
    smirks_list: list of tuples (label, smirks)
        This is the ordered list of SMIRKS from the smarts file
        if a label is provided in the file then it is assigned, otherwise
        the indices from the file is used
    """
    import os
    smarts_folder = get_data_path(os.path.join('smarts_files', smarts_file_name))
    fn = get_full_path(file_path)
    f = open(fn)
    lines = f.readlines()
    f.close()

    lines = [l.split() for l in lines]
    if len(lines[0]) > 1:
        return [(label, smirks) for smirks, label in lines]

    return [(str(i), l[0]) for i, l in enumerate(lines)]



mol_files = ['AlkEthOH_filtered_tripos.mol2']
fragments = ['nonbond', 'bond', 'angle', 'proper_torsion']
pairs = list(product(mol_files, fragments))

@pytest.mark.parametrize('mol_file, frag', pairs)
def test_complex_clusters(mol_file, frag):
    # load molecules and set aromaticity
    mols_list = mols_from_mol2(mol_file)
    for m in mols_list:
        m.set_aromaticity_mdl()

    # load smirks list
    smirks_list = parse_smarts_file("%s_smirks.smarts" % frag)

    # type molecules with these smirks
    clusters = create_tuples_for_clusters(smirks_list, mols_list)

    # smirksify these clusters
    create_frag = SMIRKSifier(mols_list, clusters, verbose=False, max_layers=2)

    # try reducing these SMIRKS
    create_frag.reduce(10)
