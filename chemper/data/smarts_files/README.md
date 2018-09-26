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
in this directory were created with the following code snippet
using an intermediate file from the Open Force Field Initiative
for smirnoff99frosst:
[`smirnoffishFrcmod.parm99Frosst.txt`](https://github.com/openforcefield/openforcefield/blob/master/utilities/convert_frosst/smirnoffishFrcmod.parm99Frosst.txt)

```python

f = open('smirnoffishFrcmod.parm99Frosst.txt', 'r')
lines = f.readlines()
f.close()

cats = ["NONBON", "BOND", "ANGL", "DIHE", "IMPR"]
smirks_by_cat = {c:list() for c in cats}

for l in lines:
    strip = l.strip()
    if strip in cats:
        cat = strip
    if len(strip) == 0 or strip[0] != '[':
        continue

    smirks_by_cat[cat].append(strip.split()[0])

cat_file = {"NONBON": "nonbond_smirks.smarts",
            "BOND": "bond_smirks.smarts",
            "ANGL": "angle_smirks.smarts",
            "DIHE": "proper_torsion_smirks.smarts",
            "IMPR": "improper_torsion_smirks.smarts"}
for cat, smirks in smirks_by_cat.items():
    f = open(cat_file[cat], 'w')
    for s in smirks:
        f.write(s+'\n')
    f.close()
```

