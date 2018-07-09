"""
This is a rather straight forward set of tests intended to make sure classes in
cluster_graph.py and fragment_graph.py run without failures.
"""

from chemper.graphs.cluster_graph import ClusterGraph
from chemper.graphs.fragment_graph import ChemPerGraphFromMol, ChemPerGraph
from chemper.mol_toolkits import mol_toolkit
import pytest


@pytest.mark.parametrize('graph_method', [ClusterGraph, ChemPerGraph])
def test_empty_graph(graph_method):
    """
    Test basic function of empty graphs
    """
    c = graph_method()
    assert c.as_smirks() is None

    new_atom = c.add_atom(None)
    assert c.as_smirks() == '[*]'

    newer_atom = c.add_atom(None, None, new_atom)
    assert c.as_smirks() == '[*]~[*]'

    c = graph_method()
    c.add_atom(None, new_smirks_index=1)
    assert c.as_smirks() == '[*:1]'


# Make sure tests run without failures for a variety of simple molecules
smiles_set = ['C', 'c1ccccc1', 'C1CC1F', 'OC(=O)CC', 'c1ccccc1', '[O-1]c1ccccc1']
layers_options = [1,3,'all']

frag_combos = [(s, l) for s in smiles_set for l in layers_options]
@pytest.mark.parametrize('smile,layers',frag_combos)
def test_no_fail_fragment(smile, layers):
    mol = mol_toolkit.MolFromSmiles(smile)
    smirks_dict = {1: 0, 2: 1}
    c = ChemPerGraphFromMol(mol, smirks_dict, layers)
    assert c.add_atom(None) is None


cluster_combos = [([smiles_set[i], smiles_set[i+1]], l)
                  for i in range(len(smiles_set)-1) for l in layers_options]
@pytest.mark.parametrize('smiles_list,layers',cluster_combos)
def test_no_fail_cluster(smiles_list, layers):
    smirks_dict_list = [[{1: 0, 2: 1}, {1:1, 2:2}]] * len(smiles_list)
    mols_list = [mol_toolkit.MolFromSmiles(s) for s in smiles_list]
    c = ClusterGraph(mols_list, smirks_dict_list, layers=layers)
    assert c.add_atom(None) is None

def test_mols_mismatch():
    mols_list = [mol_toolkit.MolFromSmiles('CC')]
    smirks_dict_list = [[{1:0, 2:1}], [{1:2, 2:2}]]
    with pytest.raises(Exception):
        ClusterGraph(mols_list, smirks_dict_list)
