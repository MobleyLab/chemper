# <img src="https://github.com/mobleylab/chemper/blob/master/chemper_logo.svg" height=150>

[![Documentation Status](https://readthedocs.org/projects/chemper/badge/?version=latest)](http://chemper.readthedocs.io/en/latest/?badge=latest)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/MobleyLab/chemper.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/MobleyLab/chemper/context:python)
[![codecov](https://codecov.io/gh/MobleyLab/chemper/branch/master/graph/badge.svg)](https://codecov.io/gh/MobleyLab/chemper)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![DOI](https://zenodo.org/badge/120546688.svg)](https://zenodo.org/badge/latestdoi/120546688)

This repository contains a variety of tools that will be useful in automating the process
of chemical perception for the new SMIRKS Native Open Force Field (SMIRNOFF) format
as a part of the [Open Force Field Initiative](http://openforcefield.org) [1].

`ChemPer` can be used to automatically generate SMIRKS patterns to match clustered molecular fragments.
For example, you may have calculated bond lengths and force constants for a variety of bonds in one group of molecules.
You could use that data to cluster those bonds and then use `ChemPer` to generate SMIRKS patterns which would allow
you to apply those lengths and force constants to a new set of molecules.
The algorithms implemented here were inspired by
[SMARTY and SMIRKY](https://github.com/openforcefield/smarty) which were proven to be too inefficient for
practical use in force field parameterization [2].

For a more extensive history and explanation, see our [preprint](http://doi.org/10.26434/chemrxiv.8304578.v1) [3].

## Installation

Chemper is available via `conda-forge`:

```shell
conda install -c conda-forge chemper
```

This command will install all dependencies besides a toolkit for cheminformatics or storing of molecule
information. Also install [RDKit](http://www.rdkit.org/docs/Install.html) and/or
[OpenEye toolkits](https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html) by
running:

```bash
conda install -c conda-forge rdkit
```

and/or

```bash
conda install -c openeye openeye-toolkits
```


## Supported Python versions

We test with whatever Python versions are found in `.github/workflows/ci.yaml`. Chemper may function
on some older and/or newer versions as well.

## Supported chemiformatics toolkits

We seek to keep this tool independent of cheminformatics toolkit, but currently only support
[RDKit](http://www.rdkit.org/docs/index.html) and [OpenEye Toolkits](https://www.eyesopen.com/).
If you wish to add support please feel free to submit a pull request.
Make sure one of these toolkits is installed in your environment before installing chemper.

# Documentation

Below are some details on the tools provided in `chemper` see
[examples](https://github.com/MobleyLab/chemper/tree/master/examples)
and [documentation](https://chemper.readthedocs.io/en/latest/)
for more detailed usage examples

### SMIRKSifier

This is `chempers` main function.
It takes groups of molecular fragments which should be typed together and generates a heirarchical list
of SMIRKS patterns which maintains this typing.
`chemper`'s `SMIRKSifier` takes a list of molecules and groups of atoms based on index and generates
a hierarchical list of SMIRKS in just a few lines of code.
In the example, [general_smirks_for_clusters](https://chemper.readthedocs.io/en/latest/examples/general_smirks_for_clusters.html)
we cluster bonds in a set of simple hydrocarbons based on order. Then `SMIRKSifer` turns these clusters into a list of SMIRKS patterns.
The following functionalities are used to make the `SMIRKSifier` possible, but may be useful on their own.

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

mol1 = mol_toolkit.Mol.from_smiles('CCC')
mol2 = mol_toolkit.Mol.from_smiles('CCCCC')
smirks_atom_lists = [[(0,1)], [(0,1), (1,2)]]
graph = ClusterGraph([mol1, mol2], smirks_atom_lists)
print(graph.as_smirks())
# '[#6AH2X4x0r0+0,#6AH3X4x0r0+0:1]-;!@[#6AH2X4x0r0+0:2]'
```

The idea with ClusterGraph objects is that they store all possible decorator information for each atom.
In this case the SMIRKS indexed atoms for propane (mol1) are one of the terminal and the middle carbons.
In pentane (mol2) however atom1 can be a terminal or middle of the chain carbon atom. This changes the number of
hydrogen atoms (`Hn` decorator) on the carbon, thus there are two possible SMIRKS patterns for atom `:1`
`#6AH2X4x0r0+0` or (indicated by the "`,`") `#6AH3X4x0r0+0`. But, atom `:2` only has one possibility `#6AH2X4x0r0+0`.

### SingleGraph

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
from chemper.graphs.single_graph import  SingleGraph

mol = mol_toolkit.Mol.from_smiles('C=C') # note this adds explicit hydrogens to your molecule
smirks_atoms = (0,1)
graph = SingleGraph(mol, smirks_atoms, layers=1)
print(graph.as_smirks())
# [#6AH2X3x0r0+0:1](-!@[#1AH0X1x0r0+0])(-!@[#1AH0X1x0r0+0])=!@[#6AH2X3x0r0+0:2](-!@[#1AH0X1x0r0+0])-!@[#1AH0X1x0r0+0]
```

### mol\_toolkits

As noted [above](#installation), we seek to keep `chemper` independent of the underlying cheminformatics toolkits.
`mol_toolkits` was created to keep all code dependent on the toolkit isolated. It can create molecules from
an RDK or OE molecule object or from a SMILES string. It includes a variety of functions for extracting information
about atoms, bonds, and molecules. Also included here are subsearchs using indexed SMARTS (or SMIRKS) patterns.

## Versions

### 0.1.0 Alpha Release
This is a first release of the Alpha testing version of `chemper`. As you can follow in the [issue tracker](https://github.com/MobleyLab/chemper/issues) there are still
on going problems to resolve. This first release will allow for reference to the concepts and algorithms included here
for automated chemical perception. However, the API is still in flux and nothing should be considered permanent at this time.

### Version 1.0.0
This release matches the work published in our [preprint](http://doi.org/10.26434/chemrxiv.8304578.v1).
While the code is stable and there are tests showing how it should work the science it represents is still in the early stages and big changes to the algorithms and API should be expected in future releases.

### Version 1.0.1
This release includes non-behavior-breaking changes to support distribution on Python 3.10.

## Contributors

* [Caitlin C. Bannan (UCI)](https://github.com/bannanc)
* [Jessica Maat (UCI)](https://github.com/jmaat)
* [David L. Mobley (UCI)](https://github.com/davidlmobley)

## Acknowledgments

CCB is funded by a fellowship from [The Molecular Sciences Software Institute](http://molssi.org/) under NSF grant ACI-1547580.

## References

1. D.L. Mobley et al. _JCTC,_ **2018**, _14_(11), pp 6076-6092. ([JCTC](http://doi.org/10.1021/acs.jctc.8b00640) or [bioRxiv](http://doi.org/10.1101/286542))
2. C. Zanette and C.C. Bannan et al. _JCTC_ **2019** _15_(1), pp 402-423. ([JCTC](https://doi.org/10.1021/acs.jctc.8b00821) or [ChemRxiv](https://doi.org/10.26434/chemrxiv.6230627.v1))
3. C.C. Bannan and D.L. Mobley _ChemRxiv_ **2019** [doi:10.26434/chemrxiv.8304578.v1](http://doi.org/10.26434/chemrxiv.8304578.v1)
