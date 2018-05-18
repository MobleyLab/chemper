[![Build Status](https://travis-ci.org/MobleyLab/chemper.svg?branch=master)](https://travis-ci.org/MobleyLab/chemper) [![codecov](https://codecov.io/gh/MobleyLab/chemper/branch/master/graph/badge.svg)](https://codecov.io/gh/MobleyLab/chemper) [![Documentation Status](https://readthedocs.org/projects/chemper/badge/?version=latest)](http://chemper.readthedocs.io/en/latest/?badge=latest)

# chemper

This repository contains a variety of tools that will be useful in automating the process
of chemical perception for the new SMIRKS Native Open Force Field (SMIRNOFF) format
as a part of the [Open Force Field Consortium](http://openforcefield.org) [1].

This idea originated from the tools [SMARTY and SMIRKY](https://github.com/openforcefield/smarty) which
were designed to use an automated monte carlo algorithm to sample the chemical perception used
in existing force fields. SMARTY samples SMARTS patterns corresponding to traditional atom types and was
tested compared to the parm99/parm@frosst force field. SMIRKY is an extension of SMARTY created to sample SMIRKS
patterns corresponding to SMIRNOFF parameter types (nonbonded, bond, angle, and proper and improper torsions).

One of the most important lessons learned while testing SMARTY and SMIRKY is that the combinatorial problem
in SMIRKS space is very large. These tools currently use very naive moves in SMIRKS space chosing atom or
bond decorators to add to a pattern at complete random. This wastes a signficant amount of time making
chemically insensible moves. One of the take away conclusions on that project was that future chemical perception
sampling tools would need to take atom and bond information from input molecules in order to be feasible [2].

We developed `chemper` based on the knowledge of the SMARTY project outcomes.
The goal here is to take clustered molecular subgraphs and generate SMIRKS patterns.
These tools will use information stored in the atoms and bonds of a molecule to drive
choices about SMIRKS decorators. Then will automatically generate reasonable SMIRKS patterns
matching clustering of molecular subgraphs.

For example, if you know you want to assign certain group of angles (sets of three atoms)
the same equilibrium bond angle and force constant,
then chemper should generate SMIRKS patterns that maintain that clustering.

## Prerequisites

We currently test with Python 3.5, though we expect anything 3.5+ should work.

This is a python tool kit with a few dependencies. We recommend installing
[miniconda](http://conda.pydata.org/miniconda.html). Then you can create an
environment with the following commands:

```bash
conda create -n [my env name] python=3.5 numpy networkx pytest
source activate [my env name]
```

This command will install all dependencies besides a toolkit for cheminformatics or storing of molecule
information. We seek to keep this tool independent of cheminformatics toolkit, but currently only support
[RDKit](http://www.rdkit.org/docs/index.html) and [OpenEye Toolkits](https://www.eyesopen.com/).
If you wish to add support please feel free to submit a pull request.
Make sure one of these tool kits is installed in your environment before installing chemper.

#### RDKit environment

Conda installation according to [RDKit documentation](http://www.rdkit.org/docs/Install.html):
```bash
conda install -c rdkit rdkit
```

#### OpenEye environment
Conda installation according to [OpenEye documentation](https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html)
```bash
cond install -c openeye openeye-toolkits
```

## Installation

Hopefully chemper will be conda installable in the near future, but for now the best option
is to download or clone this repository and then install `chemper` from inside the `chemper` directory with:
```bash
pip install -e .
```

# Documentation

Below are some details on the tools provided in `chemper` see
[examples](https://github.com/MobleyLab/chemper/tree/master/examples) for more detailed usage examples

### mol_toolkits

As noted [above](#installation), we seek to keep `chemper` independent of the cheminformatics toolkit.
`mol_toolkits` is created to keep all code dependent on the toolkit installed. It can create molecules from
an RDK or OE molecule object or from a SMILES string. It includes a variety of functions for extracting information
about atoms, bonds, and molecules. Also included here are SMIRKS pattern searches.

### ChemPerGraph

The goal of this tool was to create an example of how you could create a SMIRKS pattern from a
molecule and set of atom indices.
While this isn't ultimately useful in sampling chemical perception as they
only work for a single molecule, however it is a tool that did not exist to the best of the authors knowledge before.
For a detailed example see the [single_mol_smirks](examples/using_fragment_graph/single_mol_smirks.ipynb)
jupyter notebook.

Here is a brief usage example for creating the SMIRKS pattern for the bond between the two carbon
atoms in ethene including atoms one bond away from the indexed atoms. The indexed atoms are the two carbon
atoms at indices 0 and 1 in the molecule are assigned to SMIRKS indices `:1` and `:2` respectively

```python
from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.fragment_graph import  ChemPerGraphFromMol

mol = mol_toolkit.MolFromSmiles('C=C') # note this adds explicit hydrogens to your molecule
smirks_dict = {1:0, 2: 1}
graph = ChemPerGraphFromMol(mol, smirks_dict, layers=1)
print(graph.as_smirks)
# [#6AH2X3x0r0+0:1](-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])=!@[#6AH2X3x0r0+0:2](-!@[#1AH0X1x0r0+0])-!@[#1AH0X1x0r0+0]
```

### ClusterGraph

The goal of this tool is to store all information about the atoms and bonds that could be in a SMIRKS pattern.
These are created assuming you already have a clustered set of molecular subgraphs. As our primary goal is to
determine chemical perception for force field parameterization we image the input data being clustered subgraphs
based on what parameter we wish to assign those atoms, such as equilibrium
bond lengths and force constants. However, you could imagine other reasons for wanting to store how you clustered groups
of atoms.

For more detailed examples and illustration of how this works see [SMIRKS_from_molecules](examples/using_cluster_graph/SMIRKS_from_molecules.ipynb).
Below is a brief example showing the SMIRKS for the bond between two carbon atoms in propane and pentane.

```python
from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.cluster_graph import ClusterGraph

mol1 = mol_toolkit.MolFromSmiles('CCC')
mol2 = mol_toolkit.MolFromSmiles('CCCCC')
smirks_dicts = [[{1:0, 2:1}], [{1:0,2:1}, {1:1, 2:2}]]
graph = ClusterGraph([mol1, mol2], smirks_dicts)
print(graph.as_smirks())
# '[#6AH2X4x0r0+0,#6AH3X4x0r0+0:1]-;!@[#6AH2X4x0r0+0:2]'
```

The idea with ClusterGraph objects is that they store all possible decorator information for each atom.
In this case the SMIRKS indexed atoms for propane (mol1) are one of the terminal and the middle carbons.
In pentane (mol2) however atom1 can be a terminal or middle of the chain carbon atom. This changes the number of
hydrogen atoms (`Hn` decorator) on the carbon, thus there are two possible SMIRKS patterns for atom `:1`
`#6AH2X4x0r0+0` or (indicated by the "`,`") `#6AH3X4x0r0+0`. But, atom `:2` only has one possibility `#6AH2X4x0r0+0`.


# What's coming next?

The ClusterGraph code can accurately create a SMIRKS pattern for a group of clustered molecular subgraphes.
However as you can see in the [SMIRKS_from_molecules](examples/using_cluster_graph/SMIRKS_from_molecules.ipynb)
the SMIRKS created by `ClusterGraph` are highly specific.
Our final goal here is to maintain a given set of clustering, but to generate relatively general SMIRKS patterns
that can then be used to assign force field parameters. This use of relatively generic SMIRKS patterns is part of what
we believe makes the SMIRNOFF force field format so powerful.

One option here would be to go back to the MC sampling used in [SMARTY/SMIRKY](https://github.com/openforcefield/smarty).
Where the information stored in `ClusterGraph` could be used to make more intelligent/efficient moves.
However, we believe there is a yet more efficient option where the differences and similarities in `ClusterGraph` objects
could be used to determine the correct SMIRKS patterns. Thus, the next step is to essentially find the similarities
and differences in `ClusterGraph`s so that the most general SMIRKS can be used to maintain clustering as the user
inputs, but not specify more information than necessary.

## Contributors

* [Caitlin Bannan (UCI)](https://github.com/bannanc)
* [David L. Mobley (UCI)](https://github.com/davidlmobley)

## Acknowledgments

CCB is funded by a fellowship from [The Molecular Sciences Software Institute](http://molssi.org/) under NSF grant ACI-1547580.

## References

1. D. Mobley et al. bioRxiv 2018, 286542. [doi.org/10.1101/286542](http://doi.org/10.1101/286542)
2. C. Zanette and C.C. Bannan et al. chemRxiv 2018, [doi.org/10.26434/chemrxiv.6230627.v1](https://doi.org/10.26434/chemrxiv.6230627.v1)
