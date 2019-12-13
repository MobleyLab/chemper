from chemper.mol_toolkits.mol_toolkit import Mol
from chemper.graphs.single_graph import SingleGraph

# make molecule from SMILES
smiles = 'CCO'
mol = Mol.from_smiles(smiles)

tagged = (0,1)   # atom in carbon-carbon bond
# try multiple options for layers
for layers in [0, 1, 'all']:
    # make graph and extract SMIRKS
    graph = SingleGraph(mol, tagged, layers)
    print(graph.as_smirks())   # complex SMIRKS with all decorators are the default
    print(graph.as_smirks(compress=True))   # compressed SMIRKS have only atomic numbers