# Examples

This is directory stores examples for how to use and work with ChemPer

### Single Mol SMIRKS

`ChemPerGraph`s convert a molecule into a SMIRKS pattern by 
converting information about the atoms and bonds. 
The `ChemPerGraph` stores only one molecule. It can convert only certain atoms or
the *entire* molecule into a SMIRKS pattern. 
This direction of molecule to SMIRKS string is new as far as the authors know with chemper. 

### SMIRKS from Molecules 

`ClusterGraph`s are an expansion of the initial `ChemPerGraph`. 
They can store information about multiple atoms and bonds simultaneously.
This is a starting example for how to use ChemPer's `ClusterGraph` class 
to create SMIRKS patterns from clusters of molecular graphs. 

