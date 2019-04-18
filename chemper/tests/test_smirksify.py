"""
This script is used to test the SMIRKSifier Class and the class methods it contains
"""

from chemper.mol_toolkits import mol_toolkit
from chemper.smirksify import SMIRKSifier, print_smirks, ClusteringError, Reducer
import pytest
import copy

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
    red = SMIRKSifier([mol], cluster_lists, max_layers=0)
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
    print_smirks(red.current_smirks)
    # run for a long time (assumed to hit all possible methods)
    smirks_list = red.reduce(2000)

def test_reducer_methods():
    """
    Due to the random nature of this method, we should
    explicitly check each method
    """
    # practice reducer
    mol = mol_toolkit.MolFromSmiles('C')
    red = Reducer([('a','[*:1]~[*:2]')], [mol])

    # check generic SMIRKS output
    out_smirks, changed = red.remove_decorator("[*:1]~[*:2]")
    assert out_smirks == "[*:1]~[*:2]"
    assert not changed

    # check explicit output for removing "OR"
    new, changed = red.remove_or([])
    assert not changed

    # check explicit output for removing "AND"
    new, changed = red.remove_and([])
    assert not changed

    # check is_bond option for remove_or
    new, changed = red.remove_or([('-', [])], True)
    assert changed
    assert new == []

    # check top method
    new, changed = red.remove_or([('#6', ['X4'])])
    assert changed

    input_ors = [('#6', ['X4']), ('#7', ['X3']) ]
    fun_output = [
        (red.remove_all_dec_type(copy.deepcopy(input_ors)), [('#6', []), ('#7', []) ]),
        (red.remove_all_bases(copy.deepcopy(input_ors)), [('*', ['X4']), ('*', ['X3']) ]),
        (red.remove_ref(copy.deepcopy(input_ors), 0), [('#7', ['X3'])]),
        (red.remove_ref_sub_decs(copy.deepcopy(input_ors), 0), [('#6', []), ('#7', ['X3'])]),
        (red.remove_one_sub_dec(copy.deepcopy(input_ors), 0), [('#6', []), ('#7', ['X3'])])
    ]

    for output, expected in fun_output:
        assert output == expected


expected_change = ["[#6:1]~[*:2]",
                   "[*X4:1]~[*:2]",
                   "[*;A:1]~[*:2]",
                   "[*:1]#[*:2]",
                   "[*:1]~[*:2]~[*]",
                   ]
@pytest.mark.parametrize('in_smirks', expected_change)
def test_expected_removal(in_smirks):
    # create reducer
    mol = mol_toolkit.MolFromSmiles('C')
    smirks_list = [('a', '[#6AH4X4x0!r+0:1]-;!@[#1AH0X1x0!r+0:2]')]
    red = Reducer(smirks_list, [mol])

    # check only OR base to remove
    red_smirks, changed = red.remove_decorator(in_smirks)
    while not changed:
        red_smirks, changed = red.remove_decorator(in_smirks)
    assert red_smirks == "[*:1]~[*:2]"

def test_failed_layers():
    mol1 = mol_toolkit.MolFromSmiles('CCC')
    mol2 = mol_toolkit.MolFromSmiles('CCCC')
    # cluster 1 has bond 0-1 in mol1
    # cluster 2 has bond 0-1 and 2-3 in mol2
    # the only way to distinguish these is
    # with more than 1 layer so max_layers = 0 will fail
    clusters = [('1', [[(0,1)], []]),
                ('2', [[], [(0,1), (3,2)]])]
    with pytest.raises(ClusteringError):
        red = SMIRKSifier([mol1, mol2], clusters, max_layers=0)
        print(red.current_smirks)
        print_smirks(red.current_smirks)

@pytest.mark.parametrize('layers',[1,2,3,5,10,100])
def test_layer_choice(layers):
    mol1 = mol_toolkit.MolFromSmiles('CCC')
    mol2 = mol_toolkit.MolFromSmiles('CCCC')
    # cluster 1 has bond 0-1 in mol1
    # cluster 2 has bond 0-1 and 2-3 in mol2
    # this should always need 1 layer to work
    clusters = [('1', [[(0,1)], []]),
                ('2', [[], [(0,1), (2,3)]])]
    red = SMIRKSifier([mol1, mol2], clusters, max_layers=layers)
    assert red.layers == 1
