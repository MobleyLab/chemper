from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.cluster_graph import ClusterGraph

mol1 = mol_toolkit.mol_from_smiles('CCC')
mol2 = mol_toolkit.mol_from_smiles('CCCCC')
atoms1 = [(0,1)]
atoms2 = [(0,1), (1,2)]
graph = ClusterGraph([mol1, mol2], [atoms1, atoms2])
print(graph.as_smirks())
# "[#6AH2X4x0r0+0,#6AH3X4x0r0+0:1]-;!@[#6AH2X4x0r0+0:2]"
