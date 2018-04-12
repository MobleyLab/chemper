"""
This test is created to check how a single molecule SMIRKS graph works
"""

from chemper.graphs.fragment_graph import ChemPerGraphFromMol
from chemper.mol_toolkits import mol_toolkit
import pytest


graph_data = [
    ('C', 0, '[#6AH4X4x0r0+0:1]-!@[#1AH0X1x0r0+0:2]'), # no layers
    ('C#C', 1, # one layer
     '[#6AH1X2x0r0+0:1](-!@[#1AH0X1x0r0+0])#!@[#6AH1X2x0r0+0:2]-!@[#1AH0X1x0r0+0]'),
    ('CO', 'all', # infinite layers
     '[#6AH3X4x0r0+0:1](-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])'\
     '-!@[#8AH1X2x0r0+0:2]-!@[#1AH0X1x0r0+0]'
     )
]

@pytest.mark.parametrize('smiles,layers,expected', graph_data)
def test_single_molecule_graph(smiles, layers, expected):
    mol = mol_toolkit.MolFromSmiles(smiles)
    smirks_dict = {1:0, 2:1}
    graph = ChemPerGraphFromMol(mol, smirks_dict, layers)
    smirks = graph.as_smirks()
    print(smirks)
    assert smirks == expected
