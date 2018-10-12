"""
This script is used to test the SMIRKSifier Class and the class methods it contains
"""

from chemper.mol_toolkits import mol_toolkit
from chemper.optimize_smirks.smirksify import SMIRKSifier
import pytest

smiles_list = ['C', 'N', 'C=C', 'C#C' 'c1ccccc1']
@pytest.mark.parametrize('smiles', smiles_list)
def test_max_reduction(smiles):
    """
    starting from a single atom with no layers,
    you should get [*:1] in a minimum of 7 steps
    """
    mol = mol_toolkit.MolFromSmiles(smiles)
    cluster_lists = [('1', [[(0,)]])]
    # create reducer
    red = SMIRKSifier([mol], cluster_lists, layers=0)
    smirks_list = red.reduce(10)
    final_smirks = smirks_list[0][1]
    assert final_smirks == '[*:1]'

def test_more_complex_reducer():
    """
    Check that all SMIRKSifier class functions at least work
    """
    smiles = ['CC', 'C=C', 'C#C']
    mols = [mol_toolkit.MolFromSmiles(s) for s in smiles]
    c1 = [[(0, 1)]]*len(smiles)
    c2 = [[(0, 2)]]*len(smiles)
    cluster_lists = [('1', c1), ('2', c2)]
    # create reducer
    red = SMIRKSifier(mols, cluster_lists, verbose=False)
    # make sure printing runs:
    red.print_smirks()
    # run for a long time (assumed to hit all possible methods)
    smirks_list = red.reduce(2000)

def test_explicitly_check_methods():
    """
    Due to the random nature of this method, we should
    explicitly check each method
    """
    # practice reducer
    mol = mol_toolkit.MolFromSmiles('C')
    cluster_lists = [('1', [[(0,)]])]
    # create reducer
    red = SMIRKSifier([mol], cluster_lists, layers=0)
    red.print_smirks()

    # check generic SMIRKS output
    out_smirks, changed = red.remove_decorator("[*:1]~[*:2]")
    assert out_smirks == "[*:1]~[*:2]"
    assert not changed

    # check explicit output for removing "OR"
    new, changed = red.remove_or([])
    assert not changed

    new, changed = red.remove_or([('#6', ['X4'])])
    assert changed
    assert new == [('#6', [])]

    # check explicit output for removing "AND"
    new, changed = red.remove_and([])
    assert not changed


expected_change = [("[#6:1]~[*:2]", "[*:1]~[*:2]" ),
                   ("[#6X4:1]~[*:2]", "[#6:1]~[*:2]"),
                   ("[*;A:1]~[*:2]", "[*:1]~[*:2]"),
                   ("[*:1]#[*:2]", "[*:1]~[*:2]")
                   ]
@pytest.mark.parametrize('in_smirks, out_smirks', expected_change)
def check_expected_removal(in_smirks, out_smirks):
    mol = mol_toolkit.MolFromSmiles('C')
    cluster_lists = [('1', [[(0,)]])]
    # create reducer
    red = SMIRKSifier([mol], cluster_lists, layers=0)

    # check only OR base to remove
    red_smirks, changed = red.remove_decorator(in_smirks)
    while not changed:
        red_smirks, changed = red.remove_decorator(in_smirks)
    assert red_smirks == out_smirks
