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
    c.add_atom(None, new_label=1)
    assert c.as_smirks() == '[*:1]'


# Make sure tests run without failures for a variety of simple molecules
# including a variety of possible atom indices
smiles_set = ['C', 'c1ccccc1', 'C1CC1F', 'OC(=O)CC', 'c1ccccc1', '[O-1]c1ccccc1']
layers_options = [1,3,'all']

frag_combos = [(s, l) for s in smiles_set for l in layers_options]
@pytest.mark.parametrize('smile,layers',frag_combos)
def test_no_fail_fragment(smile, layers):
    mol = mol_toolkit.Mol.from_smiles(smile)
    smirks_atoms = (0, 1)
    c = ChemPerGraphFromMol(mol, smirks_atoms, layers)
    assert c.add_atom(None) is None
    smirks_atoms = (0,)
    c = ChemPerGraphFromMol(mol, smirks_atoms, layers)


cluster_combos = [([smiles_set[i], smiles_set[i+1]], l)
                  for i in range(len(smiles_set)-1) for l in layers_options]
@pytest.mark.parametrize('smiles_list,layers',cluster_combos)
def test_no_fail_cluster(smiles_list, layers):
    smirks_atom_lists1 = [ [(0,1), (1,2)] ] * len(smiles_list)
    smirks_atom_lists2 = [ [(0,), (1,), (2,) ] ] * len(smiles_list)
    mols_list = [mol_toolkit.Mol.from_smiles(s) for s in smiles_list]
    c1 = ClusterGraph(mols_list, smirks_atom_lists1, layers=layers)
    c2 = ClusterGraph(mols_list, smirks_atom_lists2, layers=layers)
    assert c1.add_atom(None) is None
    assert c2.add_atom(None) is None

def test_mols_mismatch():
    """
    tests that an exception is raised when the number of molecules
    and the number of smirks dictionaries is not equal
    """
    mols_list = [mol_toolkit.Mol.from_smiles('CC')]
    smirks_atom_lists = [[ (0, 1) ], [(1, 2)]]
    with pytest.raises(Exception):
        ClusterGraph(mols_list, smirks_atom_lists)
