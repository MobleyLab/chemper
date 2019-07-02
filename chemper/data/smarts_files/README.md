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

These files were created using `smirnoffishFrcmod.parm99Frosst.txt` for 
`smirnoff99Frosst` version 1.0.7. That file is included here for 
posterity. 
The patterns were split up by fragment type 
(nonbonded, bond, angle, proper torsion, and improper torsion).
These files were used for testing and validating `ChemPer`, but
they also make writing examples simpler since the user does not 
necessarily need to have their own fragment clusters to try the tool.

The python script below was used to make the `*.smarts` files from 
`smirnoffishFrcmod.parm99Frosst.txt`:

```python
import os
from chemper import chemper_utils

# get smirnoffish file out of chemper/data/smarts_files
fn = os.path.join('smarts_files', 'smirnoffishFrcmod.parm99Frosst.txt')
smirnoffish = chemper_utils.get_data_path(fn)
f = open(smirnoffish, 'r')
lines = f.readlines()
f.close()

# use the categories in smirnoffish file to identify 
# which parameter you're working with
cats = ["NONBON", "BOND", "ANGL", "DIHE", "IMPR"]
smirks_by_cat = {c:list() for c in cats}
for l in lines:
    strip = l.strip()
    # if its a category label update the current category
    if strip in cats:
        cat = strip
    # skip lines that don't start with a SMIRKS pattern
    if len(strip) == 0 or strip[0] != '[':
        continue
    smirks_by_cat[cat].append(strip.split()[0])
    
cat_file = {"NONBON": "nonbond_smirks.smarts",
            "BOND": "bond_smirks.smarts",
            "ANGL": "angle_smirks.smarts",
            "DIHE": "proper_torsion_smirks.smarts",
            "IMPR": "improper_torsion_smirks.smarts"}

# save SMIRKS into a *.smarts file by category
for cat, smirks in smirks_by_cat.items():
    f = open(cat_file[cat], 'w')
    written_smirks = set()
    for s in smirks:
        # for proper torsions there are multiple lines for
        # the same SMIRKS so we only want to save one of each
        if s not in written_smirks:
            f.write(s+'\n')
            written_smirks.add(s)
    f.close()
```
