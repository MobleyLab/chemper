from . import adapters
from . import mol_toolkit

if mol_toolkit.HAS_OE:
    from . import cp_openeye

if mol_toolkit.HAS_RDK:
    from . import cp_rdk

if not mol_toolkit.HAS_OE and not mol_toolkit.HAS_RDK:
    raise Exception("Neither OpenEye or RDKit is installed"\
"ChemPer requires at least one of these toolkits")
