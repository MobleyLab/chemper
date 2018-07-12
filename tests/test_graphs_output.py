"""
These tests are intended to test the SMIRKS output for ClusterGraph and ChemPerGraph objects.
Here all automatically generated SMIRKS patterns are compared to expected SMIRKS patterns.
These SMIRKS patterns are VERY complicated and NOT human readable.
"""

from chemper.graphs.cluster_graph import ClusterGraph
from chemper.graphs.fragment_graph import ChemPerGraphFromMol, ChemPerGraph
from chemper.mol_toolkits import mol_toolkit
import pytest


def make_frag_graph(smiles, layers):
    """
    Generates a chemper Mol from the provided smiles and then
    uses that Mol to build a ChemPerGraph where atom 0 is assigned SMIRKS index 1
    and atom 1 is assigned SMIRKS index 2.

    The variable layers is used to set the number of atoms away from the indexed atoms to include.
    For example if layers is 0 then only the SMIRKS indexed atoms are included in the graph;
    and if layers is 1 then atoms 1 bond away from the indexed atoms are included, and so forth.
    Layers can also be "all" which will lead to all atoms in the molecule being added to the graph.
    """
    mol = mol_toolkit.MolFromSmiles(smiles)
    smirks_dict = {1:0, 2:1}
    return ChemPerGraphFromMol(mol, smirks_dict, layers)

def make_cluster_graph(smiles_list, layers=0):
    """
    Generates a chemper Mol for each of the smiles in smiles_list and then
    uses those Mols to build a ClusterGraph where the same smirks indices are used for all Mols.
    Specifically, atom 0 is assigned SMIRKS index 1 and atom 1 is assigned SMIRKS index 2.

    The variable layers is used to set the number of atoms away from the indexed atoms to include.
    For example if layers is 0 then only the SMIRKS indexed atoms are included in the graph;
    and if layers is 1 then atoms 1 bond away from the indexed atoms are included, and so forth.
    Layers can also be "all" which will lead to all atoms in the molecule being added to the graph.
    """
    smirks_dict_list = [[{1:0, 2:1}]]*len(smiles_list)
    mols_list = [mol_toolkit.MolFromSmiles(smiles) for smiles in smiles_list]
    return ClusterGraph(mols_list, smirks_dict_list, layers=layers)


# Check for expected output
graph_data = [
    (make_cluster_graph(['CC']), "[#6AH3X4x0r0+0:1]-;!@[#6AH3X4x0r0+0:2]"),
    (make_cluster_graph(['CC', 'C=C']),
     "[#6AH2X3x0r0+0,#6AH3X4x0r0+0:1]-,=;!@[#6AH2X3x0r0+0,#6AH3X4x0r0+0:2]"
     ),
    (make_cluster_graph(['CC', 'C=C', 'C1CC1']),
     "[#6AH2X3x0r0+0,#6AH2X4x2r3+0,#6AH3X4x0r0+0:1]-,=[#6AH2X3x0r0+0,#6AH2X4x2r3+0,#6AH3X4x0r0+0:2]"
     ),
    (make_cluster_graph(['CC'], layers=1),
     "[#6AH3X4x0r0+0:1](-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])" \
     "-;!@[#6AH3X4x0r0+0:2](-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])-;!@[#1AH0X1x0r0+0]"
     ),
    (make_cluster_graph(['CC'], layers='all'),
     "[#6AH3X4x0r0+0:1](-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])" \
     "-;!@[#6AH3X4x0r0+0:2](-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])-;!@[#1AH0X1x0r0+0]"
     ),
    (make_cluster_graph(['C#C'], 1), # one layer
     '[#6AH1X2x0r0+0:1](-;!@[#1AH0X1x0r0+0])#;!@[#6AH1X2x0r0+0:2]-;!@[#1AH0X1x0r0+0]'
     ),
    (make_cluster_graph(['CO'], 'all'), # infinite layers
     '[#6AH3X4x0r0+0:1](-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])' \
     '-;!@[#8AH1X2x0r0+0:2]-;!@[#1AH0X1x0r0+0]'
     ),
    (make_cluster_graph(['C#CC'], 3), # one layer
     '[#6AH1X2x0r0+0:1](-;!@[#1AH0X1x0r0+0])#;!@[#6AH0X2x0r0+0:2]-;!@[#6AH3X4x0r0+0]' \
     "(-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])-;!@[#1AH0X1x0r0+0]"
     ),
    (make_cluster_graph(['CC', 'CCC'], layers=1),
     "[#6AH3X4x0r0+0:1](-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])" \
     "-;!@[#6AH2X4x0r0+0,#6AH3X4x0r0+0:2](-;!@[#1AH0X1x0r0+0,#6AH3X4x0r0+0])" \
     "(-;!@[#1AH0X1x0r0+0])-;!@[#1AH0X1x0r0+0]"
     ),
    (make_cluster_graph(['C1CCCC1', 'C1=CNC=C1', 'CO'], layers=2),
     "[#6AH2X4x2r5+0,#6AH3X4x0r0+0,#6aH1X3x2r5+0:1](-,:[#1AH0X1x0r0+0,#6AH2X4x2r5+0,#6aH1X3x2r5+0]" \
     "(-;!@[#1AH0X1x0r0+0])(-;!@[#1AH0X1x0r0+0])-,:;@[#6AH2X4x2r5+0,#6aH1X3x2r5+0])(-;!@[#1AH0X1x0r0+0])" \
     "(-;!@[#1AH0X1x0r0+0])-,:[#6AH2X4x2r5+0,#6aH1X3x2r5+0,#8AH1X2x0r0+0:2](-;!@[#1AH0X1x0r0+0])" \
     "(-;!@[#1AH0X1x0r0+0])-,:;@[#6AH2X4x2r5+0,#7aH1X3x2r5+0](-;!@[#1AH0X1x0r0+0])-;!@[#1AH0X1x0r0+0]"
     ),
    # Make single molecule ChemPerGraphs
    (make_frag_graph('C', 0), '[#6AH4X4x0r0+0:1]-!@[#1AH0X1x0r0+0:2]'), # no layers
    (make_frag_graph('C#C', 1), # one layer
     '[#6AH1X2x0r0+0:1](-!@[#1AH0X1x0r0+0])#!@[#6AH1X2x0r0+0:2]-!@[#1AH0X1x0r0+0]'
     ),
    (make_frag_graph('C#CC', 3), # three layers
     '[#6AH1X2x0r0+0:1](-!@[#1AH0X1x0r0+0])#!@[#6AH0X2x0r0+0:2]-!@[#6AH3X4x0r0+0]' \
     "(-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])-!@[#1AH0X1x0r0+0]"
     ),
    (make_frag_graph('CO', 'all'), # infinite layers
     '[#6AH3X4x0r0+0:1](-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])' \
     '-!@[#8AH1X2x0r0+0:2]-!@[#1AH0X1x0r0+0]'
     ),
]


@pytest.mark.parametrize('graph,expected', graph_data)
def test_smirks_frag_graph(graph, expected):
    """
    Checking the smirks pattern is the easiest way to check that the
    graphs were built correctly and that new changes have not affected the output
    """
    smirks = graph.as_smirks()
    print(smirks)
    assert smirks == expected

@pytest.mark.parametrize('graph', [g for g,e in graph_data])
def test_other_cluster_graph(graph):
    """
    Check some other graph functionality
    """
    bonds = graph.get_bonds()
    atom = graph.get_atoms()[0]
    neighbor = graph.get_neighbors(atom)[0]
    bond = graph.get_connecting_bond(atom, neighbor)
    assert bond is not None

