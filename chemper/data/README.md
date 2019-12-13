# Data

This directory contains a series of files separated into subdirectories
by category which are useful for testing chemper or allowing users a 
set of real examples.

The files are imported with the chemper package and can be accessed with:
```python
import os
from chemper import chemper_utils

relative_path = os.path.join('molecules', 'MiniDrugBank_tripos.mol2')
path_in_chemper = chemper_utils.get_data_path(relative_path)
```
