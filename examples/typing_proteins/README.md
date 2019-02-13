# SMIRNOFF Protein Force Field

In this example, we will ideally show how 
chemper can be used to convert an existing force field
to the SMIRNOFF format. 
Specifically, to start we will attempt to convert the
amber ff14SB force field using the following general steps.

1. Determine a list of small di- and/or tri- peptides
2. Type the small amino acids tleap
3. Using OpenMM and oemmtools identify the force field parameters associated with each atom, bond, angle, torsion, and improper.
4. Group fragments by the parameters they are assigned. 
5. Use ChemPer to generate hierarchical lists of SMIRKS patterns for each parameter type.
6. Use the utilities script `convert_frcmod` from `openforcefield` to convert these SMIRKS and parameters into a SMIRNOFF `offxml` file. 

Below is more details about each of these steps. 

## 1. Determine which small amino acids to type

Below is a list of files with small amino acids
listed in increasing complexity, starting with a 
couple of very simple dipeptides to the process working 
and then increasing to cover all standard amino acids. 

* `test_set.smi` - 3 example dipeptides test the examples
* `dipeptides.smi` - all dipeptides possible with standard amino acids. 

## 2. Type the small amino acids with tleap

Sukanya gave me this command:
```
tleap -f leap.in
```
where `leap.in` is 
```
source leaprc.protein.ff14SB
protein = loadpdb BACE_18-oe.pdb
saveamberparm protein protein_nowat.prmtop protein_nowat.inpcrd
```
In order to automate this however, I'm going to look at using openmoltools functions the way we would to type a mol2 file with GAFF. 

## 3. Using OpenMM and oemmtools identify the force field parameters associated with each atom, bond angle, proper torsion, and improper.


## 4. Group fragments by parameters


## 5. Use ChemPer to generate hierarchical list of SMIRKS patterns for each parameter type


## 6. Conver to SMIRNOFF format
