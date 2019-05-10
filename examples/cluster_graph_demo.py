from chemper.mol_toolkits.mol_toolkit import Mol
from chemper.graphs.cluster_graph import ClusterGraph

# make molecules from smiles
mols = [
    Mol.from_smiles('CCO'),
    Mol.from_smiles('CC=C')
]
# identify atoms for tagging # one set of atoms in second molecule
tagged = [[ (0,1) ],  # one set of atoms in first molecule
          [ (0,1) ]   # one set of atoms in second molecule
          ]
# try multiple options for layers
for layers in [0,1,'all']:
    # make graph
    graph = ClusterGraph(mols, tagged, layers)
    print(graph.as_smirks())   # complex is the default output
    print(graph.as_smirks(compress=True))   # and's common decorators to the end of each atom