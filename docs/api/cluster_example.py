from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.cluster_graph import ClusterGraph

mol1 = mol_toolkit.MolFromSmiles('CCC')
mol2 = mol_toolkit.MolFromSmiles('CCCCC')
smirks_dicts = [[{1:0, 2:1}], [{1:0,2:1}, {1:1, 2:2}]]
graph = ClusterGraph([mol1, mol2], smirks_dicts)
print(graph.as_smirks())
# "[#6AH2X4x0r0+0,#6AH3X4x0r0+0:1]-;!@[#6AH2X4x0r0+0:2]"
