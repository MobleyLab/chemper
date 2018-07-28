# Data

This directory contains a series of files separated into subdirectories
by category which are useful for testing chemper or allowing users a 
set of real examples.

The files are imported with the chemper package and can be accessed with:
```python
from chemper import chemper_utils

relative_path = 'molecules/MiniDrugBank.mol2'
chemper_utils.get_data_path(relative_path)
```
