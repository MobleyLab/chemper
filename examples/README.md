# Examples

This is directory stores examples for how to use and work with ChemPer

### Generating reasonable SMIRKS patterns

`SMIRKSifier` is used to create a hierarchical list of SMIRKS patterns for 
bonds which are clustered by order. This is done by first grouping the bonds in 
a small set of molecules based on their order. Then, the `SMIRKSifier` takes those
molecules and groups of bonds and create SMIRKS patterns that are specific to these molecules.
Lastly, these SMIRKS patterns are made more general by removing unncessary decorators. 

### SMIRKS from Molecules 

`ClusterGraph`s are an expansion of the initial `SingleGraph`. 
They can store information about multiple atoms and bonds simultaneously.
This is a starting example for how to use ChemPer's `ClusterGraph` class 
to create SMIRKS patterns from clusters of molecular graphs. 

### Single Mol SMIRKS

`SingleGraph`s convert a molecule into a SMIRKS pattern by 
converting information about the atoms and bonds. 
The `SingleGraph` stores only one molecule. It can convert only certain atoms or
the *entire* molecule into a SMIRKS pattern. 
This direction of molecule to SMIRKS string is new as far as the authors know with chemper. 
