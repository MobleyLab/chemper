from . import adapters
try:
    from . import cp_openeye as mol_toolkit
except ImportError:
    try:
        from . import cp_rdk as mol_toolkit
    except Exception as e:
        print(e)
        print("Warning: No cheminformatics package (openeye, rdkit) was found")
