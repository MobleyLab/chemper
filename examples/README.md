# Examples

This is directory stores examples for how to use and work with ChemPer

### Using Cluster Graph

This is a starting example for how to use ChemPer's `ClusterGraph` class 
to create SMIRKS patterns from clusters of molecular graphs. 

The primary file is a jupyter notebook, `SMIRKS_from_molecules.ipynb`. Other files are there for support

Note - for now this notebook has some functions that will probably be moved 
into the `chemper.utils` script for general use, but part of building this 
notebook was thinking about what other methods we will need or need for a
temporary amount of time.  

### Conda environments

These are files can be used to create a replica of a conda environments developers have been working in. This is a temporary measure until ChemPer can be installed using conda. The command to create a new conda environment from a `yml` file are:
 
```
conda env create -f [environment.yml]
source activate [environment name]
```

For more details on conda environments see [Managing Environments](https://conda.io/docs/user-guide/tasks/manage-environments.html) in the anaconda documentation. 

* `oe_environment` to use openeye toolkits
and 
* `oe_rdk` to use RDKit
