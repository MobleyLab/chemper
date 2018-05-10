[![Build Status](https://travis-ci.org/MobleyLab/chemper.svg?branch=master)](https://travis-ci.org/MobleyLab/chemper) [![codecov](https://codecov.io/gh/MobleyLab/chemical_perception/branch/master/graph/badge.svg)](https://codecov.io/gh/MobleyLab/chemical_perception)
# ChemPer

This repository contains a variety of tools that will be useful in automating the process
of chemical perception for the new SMIRKS Native Open Force Field (SMIRNOFF) format 
as a part of the [Open Force Field Consortium](http://openforcefield.org) [1]. 

This idea originated from the tools [SMARTY and SMIRKY](https://github.com/openforcefield/smarty) which 
were designed to use an automated monte carlo algorithm to sample the chemical perception used 
in existing force fields. SMARTY samples SMARTS patterns corresponding to traditional atom types and was 
tested compared to the parm99/parm@frosst force field. SMIRKY is an extension of SMARTY created to sample SMIRKS 
patterns corresponding to SMIRNOFF parameter types (nonbonded, bond, angle, and proper and improper torsions) [2].

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

## Contributors

* [Caitlin Bannan (UCI)](https://github.com/bannanc)
* [David L. Mobley (UCI)](https://github.com/davidlmobley)

## Acknowledgments

CCB is funded by a fellowship from [The Molecular Sciences Software Institute](http://molssi.org/) under NSF grant ACI-1547580.

## References

1. D. Mobley et al. bioRxiv 2018, 286542. [doi.org/10.1101/286542](www.doi.org/10.1101/286542)
2. C. Zanette and C.C. Bannan et al. chemRxiv 2018, [doi.org/10.26434/chemrxiv.6230627.v1](www.doi.org/10.26434/chemrxiv.6230627.v1)