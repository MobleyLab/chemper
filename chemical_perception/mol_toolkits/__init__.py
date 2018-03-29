from . import adapters

from pkgutil import iter_modules
modules = (name for _, name, _ in iter_modules())

if 'openeye' in modules:
    from . import cp_openeye as mol_toolkit
elif 'rdkit' in modules:
    from . import cp_rdk as mol_toolkit
else:
    print("No cheminformatics package (openeye, rdkit) was found")
