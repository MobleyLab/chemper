from pkgutil import iter_modules

def toolkit_factory():
    modules = (name for _, name, _ in iter_modules())
    support = ['openeye', 'rdkit']
    has_package = {s:s in modules for s in support}

    if has_package['openeye']:
        from chemical_perception.mol_toolkits import cp_openeye
        return cp_openeye
    if has_package['rdkit']:
        from chemical_perception.mol_toolkits import cp_rdk
        return cp_rdk

    raise Exception("No molecule toolkit was found (%s)" % ', '.join(support))


toolkit_pkg = toolkit_factory()

class Mol(toolkit_pkg.Mol): pass

class Atom(toolkit_pkg.Atom): pass

class Bond(toolkit_pkg.Bond): pass
