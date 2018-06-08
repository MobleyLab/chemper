from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.fragment_graph import  ChemPerGraphFromMol

mol = mol_toolkit.MolFromSmiles('C=C') # note this adds explicit hydrogens to your molecule
smirks_dict = {1:0, 2: 1}
graph = ChemPerGraphFromMol(mol, smirks_dict, layers=1)
print(graph.as_smirks)
# [#6AH2X3x0r0+0:1](-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])=!@[#6AH2X3x0r0+0:2](-!@[#1AH0X1x0r0+0])-!@[#1AH0X1x0r0+0]
