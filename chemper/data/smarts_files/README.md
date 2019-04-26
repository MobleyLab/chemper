# SMARTS Files

SMARTS files (`*.smarts`) are similar to smiles (`*.smi`) files.
They simply contain a list of SMARTS strings, but they can also include
a label one word for each SMARTS/SMIRKS with a space separating the
label and the pattern.

### From smirnoff99Frosst

The SMARTS files:
* `nonbond_smirks.smarts`
* `bond_smirks.smarts`
* `angle_smirks.smarts`
* `proper_torsion_smirks.smarts`
* `improper_torsion_smirks.smarts`
in this directory were created using the intermediate file
[`smirnoffishFrcmod.parm99Frosst.txt`](https://github.com/openforcefield/openforcefield/blob/master/utilities/convert_frosst/smirnoffishFrcmod.parm99Frosst.txt).
In each there are the SMARTS patterns from SMIRNOFF99Frosst for the given category. 
These patterns are used while testing and validating ChemPer.
