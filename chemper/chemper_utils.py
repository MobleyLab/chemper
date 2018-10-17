"""
utils.py

This file provides simple functions that might be of use to people using the chemper package

"""

import os

def get_data_path(relative_path, package='chemper'):
    """
    Returns the absolute path to a specified relative path inside data directory of a given package
    there for this will find a file at:

    [absolute path to Python packages]/[package(chemper]/data/[filename]

    Parameters
    ----------
    relative_path: str
              filename or relative path to a file stored in the data/ directory
              of a given package
    package: str
             name of the Python package that should contain the specified file
             default: "chemper"

    Returns
    -------
    path: str
          absolute path to filename in specified package
    """
    from pkg_resources import resource_filename
    fn = resource_filename(package, os.path.join('data', relative_path))

    if not os.path.exists(fn):
        raise IOError("Sorry! The absolute path %s was not found, that is %s is not in the data directory for %s" \
                          % (fn, relative_path, package))

    return fn

def get_full_path(relative_path, package="chemper"):
    """
    This function checks to see if a file/path at this relative path is available locally
    if it is it returns the absolute path to that file.
    Otherwise it checks if that relative path is located at [package location]/data/relative path
    it will raise an IOError if neither location exists.

    Parameters
    ----------
    relative_path: str
                   relative path to a file available locally or in package/data
    package: str
             package name with a data directory the file might be located in

    Returns
    -------
    abs_path: str
              The absolute path to a local file if available, otherwise the abolsute path
              to [package]/data/[relative_path]
    """
    if os.path.exists(relative_path):
        return os.path.abspath(relative_path)
    return get_data_path(relative_path, package)

# =======================================
# Check SMIRKS validity
# =======================================

def is_valid_smirks(smirks):
    """
    This function checks if a given smirks is valid.

    Parameters
    ----------
    smirks: SMIRKS (or SMARTS) string

    Returns
    -------
    is_valid: boolean
              is the provided SMIRKS a valid pattern
    """
    from chemper.mol_toolkits import mol_toolkit

    mol = mol_toolkit.MolFromSmiles('C')
    try:
        mol.smirks_search(smirks)
        return True
    except ValueError:
        return False

# ===================================================================
# custom classes and functions for matching sets of SMIRKS to molecules
# ===================================================================
"""
custom_dicts

This custom dictionaries are taken from
openforcefield.typing.engines.smirnoff.forcefield
Like other openforcefield typing tools, I wanted to be
able to test ideas for chemper without adding dependencies
assuming things work and openforcefield all gets
merged to support OE/RDK I will just import these from there.
"""


import collections
class TransformedDict(collections.MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

    def __str__(self): return str(self.store)


class ValenceDict(TransformedDict):
    """Enforce uniqueness in atom indices"""
    def __keytransform__(self, key):
        """Reverse tuple if first element is larger than last element."""
        # Ensure key is a tuple.
        key = tuple(key)
        # Reverse the key if the first element is bigger than the last.
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key


class ImproperDict(TransformedDict):
    """Symmetrize improper torsions"""
    def __keytransform__(self,key):
        """Reorder tuple in numerical order except for element[1] which is the central atom; it retains its position."""
        # Ensure key is a tuple
        key = tuple(key)
        # Retrieve connected atoms
        connectedatoms = [key[0], key[2], key[3]]
        # Sort connected atoms
        connectedatoms.sort()
        # Re-store connected atoms
        key = tuple( [connectedatoms[0], key[1], connectedatoms[1], connectedatoms[2]])
        return key


def get_typed_molecules(smirks_list, molecules):
    """
    Creates a dictionary assigning a typename
    for each set of atom indices in each molecule

    Parameters
    ----------
    smirks_list: list of tuples in the form (label, smirks)
    molecules: list of chemper Mols

    Returns
    -------
    typeDict: embedded dictionary
        keys: SMILES string for each molecule
            keys: tuple of indices assigned a parameter type
    """
    type_dict = dict()
    for mol_idx, mol in enumerate(molecules):
        type_dict[mol_idx] = {}
        for [label, smirks] in smirks_list:
            matches = get_smirks_matches(mol, smirks)
            for match in matches:
                type_dict[mol_idx][match] = label

    return type_dict


def create_tuples_for_clusters(smirks_list, molecules):
    """
    This function is used to get clusters of molecular fragments based
    on input SMIRKS.

    For example, let's assume you wanted to find all of the
    atoms that match this SMIRKS list
    'any', '[*:1]~[*:2]'
    'single', '[*:1]-[*:2]'

    In this case, the "any" bond would match all bonds, but then
    the "single" would match all single bonds.
    If you were typing Ethene (C=C) then you expect the double bond
    between atoms 0 and 1 to match any bond and all C-H bonds to match single.

    The output in this case would be:
    [ ('any', [[ (0, 1) ]] ),
      ('single', [[ (0, 2), (0, 3), (1,4), (1,5) ]] )
    ]

    Parameters
    ----------
    smirks_list: list of tuples
        This is a list of tuples with the form (label, SMIRKS)
    molecules: list of ChemPer Mols
        This is a list of ChemPer molecules you want to type with this list of SMIRKS

    Returns
    -------
    type_list: list of atom index tuples for each molecule for each type
        These are a list of tuples with the form (label, [ [ (atom indices by molecule)], ]),...

    """
    ordered_labels = [l for l, s in smirks_list]
    type_dict = get_typed_molecules(smirks_list, molecules)
    # start by getting things in the form
    # {label: {mol_idx: list of dictionary of indices} }
    label_dict = dict()
    for mol_idx, mol_dict in type_dict.items():
        for match, label in mol_dict.items():
            if label not in label_dict:
                label_dict[label] = [list() for i in range(len(molecules))]

            label_dict[label][mol_idx].append(match)

    # now we need to resort the mol_lists
    final_list = [ (l, label_dict[l]) for l in ordered_labels if l in label_dict]
    return final_list


def get_smirks_matches(mol, smirks):
    """
    Gets atom indices for a smirks string in a given molecule

    Parameters
    ----------
    mol : a chemper Mol
    smirks : str
             SMIRKS pattern being matched to the molecule

    Returns
    --------
    matches: list of tuples
        atom indices for labeled atom in the smirks
    """
    from chemper.graphs.environment import ChemicalEnvironment

    env = ChemicalEnvironment(smirks)
    if env.getType().lower() == 'impropertorsion':
        matches = ImproperDict()
    else:
        matches = ValenceDict()

    for match in mol.smirks_search(smirks):
        smirks_indices = sorted(list(match.keys()))
        atom_indices = tuple([match[s].get_index() for s in smirks_indices])
        matches[atom_indices] = ''

    return matches.keys()



