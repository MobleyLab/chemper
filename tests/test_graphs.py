"""
This test is created to check how a single molecule SMIRKS graph works
"""

from chemper.graphs.fragment_graph import ChemPerGraphFromMol
from chemper.graphs.cluster_graph import ClusterGraph
from chemper.mol_toolkits import mol_toolkit
import pytest


def make_frag_graph(smiles, layers):
    mol = mol_toolkit.MolFromSmiles(smiles)
    smirks_dict = {1:0, 2:1}
    return ChemPerGraphFromMol(mol, smirks_dict, layers)

def make_cluster_graph(smiles_list):
    smirks_dict_list = [{1:0, 2:1}]*len(smiles_list)
    mols_list = [mol_toolkit.MolFromSmiles(smiles) for smiles in smiles_list]
    return ClusterGraph(mols_list, smirks_dict_list)


graph_data = [
    (make_frag_graph('C', 0), '[#6AH4X4x0r0+0:1]-!@[#1AH0X1x0r0+0:2]'), # no layers
    (make_frag_graph('C#C', 1), # one layer
     '[#6AH1X2x0r0+0:1](-!@[#1AH0X1x0r0+0])#!@[#6AH1X2x0r0+0:2]-!@[#1AH0X1x0r0+0]'),
    (make_frag_graph('CO', 'all'), # infinite layers
     '[#6AH3X4x0r0+0:1](-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])'\
     '-!@[#8AH1X2x0r0+0:2]-!@[#1AH0X1x0r0+0]'
     ),
    (make_cluster_graph(['CC']), "[#6AH3X4x0r0+0:1]-;!@[#6AH3X4x0r0+0:2]"),
    (make_cluster_graph(['CC', 'C=C']),
     "[#6AH2X3x0r0+0,#6AH3X4x0r0+0:1]-,=;!@[#6AH2X3x0r0+0,#6AH3X4x0r0+0:2]"),
    (make_cluster_graph(['CC', 'C=C', 'C1CC1']),
     "[#6AH2X3x0r0+0,#6AH2X4x2r3+0,#6AH3X4x0r0+0:1]-,=[#6AH2X3x0r0+0,#6AH2X4x2r3+0,#6AH3X4x0r0+0:2]")
]


@pytest.mark.parametrize('graph,expected', graph_data)
def test_smirks_frag_graph(graph, expected):
    smirks = graph.as_smirks()
    print(smirks)
    assert smirks == expected

@pytest.mark.parametrize('graph', [g for g,e in graph_data])
def test_other_cluster_graph(graph):
    bonds = graph.get_bonds()
    atom = graph.get_atoms()[0]
    neighbor = graph.get_neighbors(atom)[0]
    bond = graph.get_connecting_bond(atom, neighbor)
    assert bond is not None
