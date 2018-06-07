
# coding: utf-8

# # Make SMIRKS from clustered subgraphs
# 
# This notebook will showcase how `chemper` `ClusterGraph`s store SMIRKS decorators for a set of molecular subgraphs. 
# Remember the ultimate goal here is to take clustered molecular subgraphs and create a set of SMIRKS patterns that would put those molecular subgraphs into the same groups. 
# 
# For example, if your initial clusters had 4 types of carbon-carbon bonds (single, aromatic, double, and triple), you would expect the final SMIRKS patterns to reflect those four categories. Keeping in mind that sorting by just bond order is unlikely to be detailed enough. 
# 
# The first step here is to store possible decorators for atoms and bonds in a given cluster. This notebook will attempt to walk through the code that currently exists and goals for how to move forward

# In[2]:


# import statements
from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.cluster_graph import ClusterGraph


# The following methods are useful for this example, but will likely be useful in general. 
# Therefore these will be moved to a utility script in chemper when I get a chance. 

# ## Clustering from other SMIRKS
# 
# This example attempts to show how `ClusterGraph` creates a SMIRKS for already clustered sub-graphs. 
# In order to create an example set of clustered subgraphs, I am taking an existing set of SMIRKS patterns
# and clustering molecular sub-graphs. 
# 
# The function `get_smirks_dict` makes these clusters for each molecule creating a dictionary in the form
# ```python
# {
#     label: [
#         {smirks_index: atom_index}
#     ]
# }
# ```

# In[3]:


def get_smirks_dict(mol, smirks_list):
    """
    mol - chemper Mol object
    smirks_list - list of tuples (SMIRKS, label)
    
    Returns a dictionary of listes
    {label: [ {smirks_index: atom_index} ] }
    """
    temp_dict = dict()
    for smirks, label in smirks_list:
        for dic in mol.smirks_search(smirks):
            atom_tuple = tuple([dic[i+1].get_index() for i in range(len(dic))])
            temp_dict[atom_tuple] = label

    label_dict = dict()
    for atom_tuple, label in temp_dict.items():
        if label not in label_dict:
            label_dict[label] = list()

        label_dict[label].append({i+1: atom_idx for i, atom_idx in enumerate(atom_tuple) })

    return label_dict


# ## Start with a single molecule
# 
# For the first example, we will start with just one molecule (ethane) and extract the clusters matching the SMIRKS patterns in the example smarts file `angles.smarts` in this directory.
# 
# Start by loading the smriks in angles.smarts, shown below

# In[4]:


smirks_list = [
    ("[*:1]~[#6X4:2]-[*:3]", "c1"),
    ("[#1:1]-[#6X4:2]-[#1:3]", "c2"),
    ("[*;r3:1]1~;@[*;r3:2]~;@[*;r3:3]1", "c3"),
    ("[*;r3:1]~;@[*;r3:2]~;!@[*:3]", "c4"),
    ("[*:1]~;!@[*;r3:2]~;!@[*:3]", "c5"),
    ("[#1:1]-[*;r3:2]~;!@[*:3]", "c6"),
    ("[#6r4:1]-;@[#6r4:2]-;@[#6r4:3]", "c7"),
    ("[!#1:1]-[#6r4:2]-;!@[!#1:3]", "c8"),
    ("[!#1:1]-[#6r4:2]-;!@[#1:3]", "c9"),
    ("[*:1]~[#6X3:2]~[*:3]", "c10"),
    ("[#1:1]-[#6X3:2]~[*:3]", "c11"),
    ("[#1:1]-[#6X3:2]-[#1:3]", "c12"),
    ("[*;r6:1]~;@[*;r5:2]~;@[*;r5;x2:3]", "c13"),
    ("[*:1]~;!@[*;r5:2]~;@[*;r5:3]", "c14"),
    ("[#8X1:1]~[#6X3:2]~[#8:3]", "c15"),
    ("[*:1]~[#6X2:2]~[*:3]", "c16"),
]
for smirks, label in smirks_list:
    print(label,'\t',smirks)


# Next we need to create the ethan molecule and type it with the SMIRKS from angles.smarts

# In[5]:


mol = mol_toolkit.MolFromSmiles('CC')
atom_index_list = get_smirks_dict(mol, smirks_list)
# ethane only has two matching SMIRKS patterns
print(atom_index_list.keys())


# Then we will look at the `ClusterGraph` for the set of atoms matching the Angle parameter `c1` (`[*:1]~[#6X4:2]-[*:3]`). That is the angle centered on a tetrahedral carbon with any two outer atoms. For ethane, this parameter is chosen when at least one of the outer atoms is not hydrogen, so H-C-C, C-C-H, or C-C-C are all assigned `c1`
# 
# This ethan molecule has 12 sets of atoms matching this label.

# In[6]:


c1_atoms = atom_index_list['c1']
for idx, dic in enumerate(c1_atoms):
    print(dic)
print(idx)


# However, `ClusterGraph` stores only stores unique options for each atom. 
# So atom `:2` is always a tetrahedral carbon, while atoms `:1` and `:3` can be hydrogens or carbons

# In[7]:


graph = ClusterGraph([mol], [c1_atoms])
print(graph.as_smirks())


# ## Multiple molecules at once
# 
# Now that you have the general idea, lets consider a more complex case,
# Lets create a `ClusterGraph` for every label the `smirks_list` from above for
# a more diverse set of molecules (from `carbon.smi`). 
# 
# To do this I created `make_cluster_graphs` which takes a list of chemper Mol ojects and a `smirks_list` like the one created above having SMIRKS and labels. It returns a dictionary storing a `ClusterGraph` for every label in smirks_list that matched any molecule provided.

# In[8]:


def make_cluster_graphs(molecules, smirks_list):
    """
    molecules - list of chemper mols
    smirks_list - list of tuples (SMIRKS, label)
    
    returns a dictionary of chemper ClusterGraphs:
    {label: ClusterGraph} object
    """
    graph_dict = dict()
    for mol in molecules:
        label_dict = get_smirks_dict(mol, smirks_list)
        for label, atom_list in label_dict.items():
            if label not in graph_dict:
                graph_dict[label] = ClusterGraph()
            graph_dict[label].add_mol(mol, atom_list)
    
    print("Created %i ClusterGraphs" % len(graph_dict))
    return graph_dict


# In[16]:


smiles = ['CC', 'CCC', 'C1CC1', 'CCCC', 'CC(C)C', 'C1CCC1', 'CCCCC', 'CC(C)(C)C', 'C1CCCC1', 'C=C', 
          'CC=C', 'CC(C)=C', 'C=C=C', 'C=CC=C', 'C\\C=C\\C', 'C/C=C\\C', 'CC(C)=C(C)C', 'C1=CCCC1', 'C1=CCC1', 
          'C1=CC=CC1', 'CC(CC)=C(CC)C', 'C#C', 'C#CC', 'CC#CC', 'C#CCC', 'C#CC#C', 'C#CCC#C', 'C1C#CCC1', 'C#CC=C', 
          'C1C#CCCC1', 'CC#CCCC#C', 'c1ccccc1', 'c1ccccc1C', 'c1cccc(C)c1C', 'c1ccc(C)cc1C', 'c1cc(C)ccc1C', 
          'c1ccc(C)c(C)c1C', 'c1cc(C)cc(C)c1C', 'c1c(C)cc(C)cc1C', 'c1ccccc1c2ccccc2', 'c1ccccc1Cc2ccccc2', 
          'c1ccccc1CCc2ccccc2', 'c1ccccc1CCCc2ccccc2']
mols = [mol_toolkit.MolFromSmiles(s) for s in smiles]
clusters = make_cluster_graphs(mols, smirks_list)


# Lets start by just looking at the SMIRKS written by ClusterGraph first

# In[17]:


for smirks, label in smirks_list: 
    if label in clusters:
        print(label,'\t',smirks)
        print(clusters[label].as_smirks())
        print()


# ## Where do you go from here
# 
# As you see above, the `ClusterGraph` SMIRKS are significantly more complicated and specific than the input SMIRKS. 
# For example, the input SMIRKS for `c1` is `[*:1]~[#6X4:2]-[*:3]`, 
# however `ClusterGraph` creates this monstrosity:
# 
# ```
# [#1AH0X1x0r0+0,#6AH0X2x0r0+0,#6AH0X2x2r5+0,#6AH0X2x2r6+0,#6AH0X3x0r0+0,#6AH0X4x0r0+0,#6AH1X3x0r0+0,#6AH1X3x2r5+0,#6AH1X4x0r0+0,#6AH2X4x0r0+0,#6AH2X4x2r5+0,#6AH2X4x2r6+0,#6AH3X4x0r0+0,#6aH0X3x2r6+0:1]-[#6AH0X4x0r0+0,#6AH1X4x0r0+0,#6AH2X4x0r0+0,#6AH2X4x2r3+0,#6AH2X4x2r4+0,#6AH2X4x2r5+0,#6AH2X4x2r6+0,#6AH3X4x0r0+0:2]-[#1AH0X1x0r0+0,#6AH0X2x0r0+0,#6AH0X2x2r5+0,#6AH0X2x2r6+0,#6AH0X3x0r0+0,#6AH0X4x0r0+0,#6AH1X3x0r0+0,#6AH1X3x2r4+0,#6AH1X3x2r5+0,#6AH1X4x0r0+0,#6AH2X4x0r0+0,#6AH2X4x2r3+0,#6AH2X4x2r4+0,#6AH2X4x2r5+0,#6AH2X4x2r6+0,#6AH3X4x0r0+0,#6aH0X3x2r6+0:3]
# ```
# 
# The next step of this project is to take the information stored in a set of `ClusterGraphs` and determine the most simple combination that will retain the matching provided by the input sub-structures. 
# There are two possible options here
# 
# 1. Use a MC type algorithm like in SMARTY/SMIRKY, but use `ClusterGraphs` to come up with better moves in SMIRKS space
# 2. Use a less random approach by choosing an order the graph objects should go in and then removing unnecessary details to replicate a generic to most detailed pattern. 
# 
