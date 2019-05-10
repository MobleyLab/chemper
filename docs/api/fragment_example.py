from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.single_graph import  SingleGraph

mol = mol_toolkit.mol_from_smiles('C=C') # note this adds explicit hydrogens to your molecule
atoms = (0, 1)
graph = SingleGraph(mol, atoms, layers=1)
print(graph.as_smirks)
# [#6AH2X3x0r0+0:1](-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])=!@[#6AH2X3x0r0+0:2](-!@[#1AH0X1x0r0+0])-!@[#1AH0X1x0r0+0]
