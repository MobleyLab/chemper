"""
These tests are intended to test the SMIRKS output for ClusterGraph and ChemPerGraph objects.
Here all automatically generated SMIRKS patterns are compared to expected SMIRKS patterns.
These SMIRKS patterns are VERY complicated and NOT human readable.
"""

from chemper.graphs.cluster_graph import ClusterGraph
from chemper.graphs.fragment_graph import ChemPerGraphFromMol, ChemPerGraph
from chemper.mol_toolkits import mol_toolkit
from chemper.chemper_utils import is_valid_smirks
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
    mol = mol_toolkit.Mol.from_smiles(smiles)
    smirks_atoms = (0,1)
    return ChemPerGraphFromMol(mol, smirks_atoms, layers)

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
    smirks_atom_lists = [ [(0,1)] ] * len(smiles_list)
    mols_list = [mol_toolkit.Mol.from_smiles(smiles) for smiles in smiles_list]
    return ClusterGraph(mols_list, smirks_atom_lists, layers=layers)


# Check for expected output
# These have the form (graph, long expected, short expected
graph_data = [
    # 0
    (make_cluster_graph(['CC']), "[#6H3X4x0!r+0A:1]-;!@[#6H3X4x0!r+0A:2]",
     "[#6H3X4x0!r+0A:1]-;!@[#6H3X4x0!r+0A:2]"),
    # 1
    (make_cluster_graph(['CC', 'C=C']),
     "[#6H2X3x0!r+0A,#6H3X4x0!r+0A:1]-,=;!@[#6H2X3x0!r+0A,#6H3X4x0!r+0A:2]",
     "[#6H2X3,#6H3X4;!r;+0;x0;A:1]-,=;!@[#6H2X3,#6H3X4;!r;+0;x0;A:2]"
     ),
    # 2
    (make_cluster_graph(['CC', 'C=C', 'C1CC1']),
     "[#6H2X3x0!r+0A,#6H2X4x2r3+0A,#6H3X4x0!r+0A:1]-,=[#6H2X3x0!r+0A,#6H2X4x2r3+0A,#6H3X4x0!r+0A:2]",
     "[#6!rH2X3x0,#6!rH3X4x0,#6H2X4r3x2;+0;A:1]-,=[#6!rH2X3x0,#6!rH3X4x0,#6H2X4r3x2;+0;A:2]"
     ),
    # 3
    (make_cluster_graph(['CC'], layers=1),
     "[#6H3X4x0!r+0A:1](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])" \
     "-;!@[#6H3X4x0!r+0A:2](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]",
     "[#6H3X4x0!r+0A:1](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])" \
     "-;!@[#6H3X4x0!r+0A:2](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]"
     ),
    # 4
    (make_cluster_graph(['CC'], layers='all'),
     "[#6H3X4x0!r+0A:1](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])" \
     "-;!@[#6H3X4x0!r+0A:2](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]",
     "[#6H3X4x0!r+0A:1](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])" \
     "-;!@[#6H3X4x0!r+0A:2](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]"
     ),
    # 5
    (make_cluster_graph(['C#C'], 1), # one layer
     '[#6H1X2x0!r+0A:1](-;!@[#1H0X1x0!r+0A])#;!@[#6H1X2x0!r+0A:2]-;!@[#1H0X1x0!r+0A]',
     '[#6H1X2x0!r+0A:1](-;!@[#1H0X1x0!r+0A])#;!@[#6H1X2x0!r+0A:2]-;!@[#1H0X1x0!r+0A]'
     ),
    # 6
    (make_cluster_graph(['CO'], 'all'), # infinite layers
     '[#6H3X4x0!r+0A:1](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])' \
     '-;!@[#8H1X2x0!r+0A:2]-;!@[#1H0X1x0!r+0A]',
     '[#6H3X4x0!r+0A:1](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])' \
     '-;!@[#8H1X2x0!r+0A:2]-;!@[#1H0X1x0!r+0A]'
     ),
    # 7
    (make_cluster_graph(['C#CC'], 3), # one layer
     '[#6H1X2x0!r+0A:1](-;!@[#1H0X1x0!r+0A])#;!@[#6H0X2x0!r+0A:2]-;!@[#6H3X4x0!r+0A]' \
     "(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]",
     "[#6H1X2x0!r+0A:1](-;!@[#1H0X1x0!r+0A])#;!@[#6H0X2x0!r+0A:2]-;!@[#6H3X4x0!r+0A]" \
     "(-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]"
     ),
    # 8
    (make_cluster_graph(['CC', 'CCC'], layers=1),
     '[#6H2X4x0!r+0A,#6H3X4x0!r+0A:1](-;!@[#1H0X1x0!r+0A,#6H3X4x0!r+0A])(-;!@[#1H0X1x0!r+0A])' \
     '(-;!@[#1H0X1x0!r+0A])-;!@[#6H3X4x0!r+0A:2](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]',
     '[#6H2,#6H3;!r;+0;X4;x0;A:1](-;!@[#1H0X1,#6H3X4;!r;+0;x0;A])(-;!@[#1H0X1x0!r+0A])' \
     '(-;!@[#1H0X1x0!r+0A])-;!@[#6H3X4x0!r+0A:2](-;!@[#1H0X1x0!r+0A])(-;!@[#1H0X1x0!r+0A])-;!@[#1H0X1x0!r+0A]'
     ),
    # Make single molecule ChemPerGraphs
    # 9
    (make_frag_graph('C', 0), '[#6AH4X4x0!r+0:1]-!@[#1AH0X1x0!r+0:2]', "[#6:1]-!@[#1:2]"), # no layers
    # 10
    (make_frag_graph('C#C', 1), # one layer
     '[#6AH1X2x0!r+0:1](-!@[#1AH0X1x0!r+0])#!@[#6AH1X2x0!r+0:2]-!@[#1AH0X1x0!r+0]',
     "[#6:1](-!@[#1])#!@[#6:2]-!@[#1]"
     ),
    # 11
    (make_frag_graph('C#CC', 3), # three layers
     '[#6AH1X2x0!r+0:1](-!@[#1AH0X1x0!r+0])#!@[#6AH0X2x0!r+0:2]-!@[#6AH3X4x0!r+0]' \
     "(-!@[#1AH0X1x0!r+0])(-!@[#1AH0X1x0!r+0])-!@[#1AH0X1x0!r+0]",
     "[#6:1](-!@[#1])#!@[#6:2]-!@[#6](-!@[#1])(-!@[#1])-!@[#1]"
     ),
    # 12
    (make_frag_graph('CO', 'all'), # infinite layers
     '[#6AH3X4x0!r+0:1](-!@[#1AH0X1x0!r+0])(-!@[#1AH0X1x0!r+0])(-!@[#1AH0X1x0!r+0])' \
     '-!@[#8AH1X2x0!r+0:2]-!@[#1AH0X1x0!r+0]',
     "[#6:1](-!@[#1])(-!@[#1])(-!@[#1])-!@[#8:2]-!@[#1]"
     ),
]

graphs = [(g,g) for g, e1, e2 in graph_data]


@pytest.mark.parametrize('graph1, graph2', graphs )
def test_other_cluster_graph(graph1, graph2):
    """
    Check graphs build successfully by looking at minimal examples
    """
    # check that a duplicate of the graph is treated as identical
    assert graph1 == graph2

    # check that some of the extraction functions work as expected
    bonds = graph1.get_bonds()
    atom = graph1.get_atoms()[0]
    neighbor = graph1.get_neighbors(atom)[0]
    bond = graph1.get_connecting_bond(atom, neighbor)
    assert bond is not None

    # check graph makes a valid smirks
    assert is_valid_smirks(graph1.as_smirks())


@pytest.mark.parametrize('graph,expected,expected_compressed', graph_data)
def test_smirks_frag_graph(graph, expected, expected_compressed):
    """
    Checking the smirks pattern is the easiest way to check that the
    graphs were built correctly and that new changes have not affected the output
    """
    smirks = graph.as_smirks()
    print(smirks)
    assert smirks == expected

    compressed_smirks = graph.as_smirks(True)
    print(compressed_smirks)
    assert compressed_smirks == expected_compressed

