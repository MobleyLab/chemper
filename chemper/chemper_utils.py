"""
utils.py

This file provides simple functions that might be of use to people using the chemper package

"""

import os
from chemper.mol_toolkits import mol_toolkit

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
    mol = mol_toolkit.MolFromSmiles('C')
    try:
        mol.smirks_search(smirks)
        return True
    except ValueError:
        return False

# =======================================
# Functions for parsing molecule files
# =======================================

def mols_fom_mol2(mol2_file, toolkit=None):
    """
    Creates a list of chemper Mols from the provided mol2 file
    using a specified or default toolkit

    Parameters
    ----------
    mol2_file: str
               path to mol2 file, this can be a relative or absolute path locally
               or the path to a molecules file stored in chemper at chemper/data/molecules/
    toolkit: None or str
             "oe" for openeye, "rdk" for RDKit
             if None then the first package available (in the order listed here)
             will be used

    Returns
    -------
    mol2s: list of chemper Mol
           list of molecules in the provided mol2 file
    """

    if os.path.exists(mol2_file):
        mol2_path = os.path.abspath(mol2_file)
    else:
        mol2_path = get_data_path('molecules/%s' % mol2_file)

    if not os.path.exists(mol2_path):
        raise IOError("Mol2 file (%s) was not found locally or in chemper/data/molecules" % mol2_file)

    if toolkit is None:
        try:
            from openeye import oechem
            toolkit = 'oe'
        except:
            try:
                from rdkit import Chem
                toolkit = 'rdk'
            except:
                raise ImportError("Could not find OpenEye or RDKit toolkits in this pythong environment")

    if toolkit.lower() == 'oe':
        return oe_mols_from_file(mol2_path)

    elif toolkit.lower() == 'rdk':
        return rdk_mols_from_mol2(mol2_path)

def oe_mols_from_file(mol_file):
    """
    Parses a standard molecule file into chemper molecules using OpenEye toolkits

    Parameters
    ----------
    mol_file: str
              relative or full path to molecule containing the molecule file
              that is accessible from the current working directory


    Returns
    -------
    mols: list of chemper Mols
          list of molecules in the mol2 file as chemper Mols
    """
    try:
        from openeye import oechem
        from chemper.mol_toolkits.cp_openeye import Mol
    except:
        raise ImportError("Could not find OpenEye Toolkits, try the rdk_mols_from_mol2 method instead")

    if not os.path.exists(mol_file):
        raise IOError("File '%s' not found." % mol_file)

    molecules = list()

    # make Openeye input file stream
    ifs = oechem.oemolistream(mol_file)

    oemol = oechem.OECreateOEGraphMol()
    while oechem.OEReadMolecule(ifs, oemol):
        # if an SD file, the molecule name may be in the SD tags
        if oemol.GetTitle() == '':
            name = oechem.OEGetSDData(oemol, 'name').strip()
            oemol.SetTitle(name)
        # Append to list.
        molecules.append(Mol(oechem.OEMol(oemol)))
    ifs.close()

    return molecules


def rdk_mols_from_mol2(mol2_file):
    """
    Parses a mol2 file into chemper molecules using RDKit

    This is a hack for separating mol2 files taken from a Source Forge discussion here:
    https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01510.html
    It splits up a mol2 file into blocks and then uses RDKit to parse those blocks

    Parameters
    ----------
    mol2_file: str
               relative or absolute path to a mol2 file you want to parse
               accessible form the current directory

    Returns
    -------
    mols: list of chemper Mols
          list of molecules in the mol2 file as chemper molecules
    """
    # TODO: check that this works with mol2 files with a single molecule
    # TODO: figure out if @<TRIPOS>MOLECULE is the only delimiter acceptable in this file type
    try:
        from rdkit import Chem
        from chemper.mol_toolkits.cp_rdk import Mol
    except:
        raise ImportError("Could not find RDKit toolkit, try the oe_mols_from_mol2 method instead")

    delimiter="@<TRIPOS>MOLECULE"

    if mol2_file.split('.')[-1] != "mol2":
        raise IOError("File '%s' is not a mol2 file" % mol2_file)

    if not os.path.exists(mol2_file):
        raise IOError("File '%s' not found." % mol2_file)

    molecules = list()
    mol2_block = list()

    file_open = open(mol2_file, 'r')

    for line in file_open:
        if line.startswith(delimiter) and mol2_block:
            rdmol = Chem.MolFromMol2Block("".join(mol2_block))
            if rdmol is not None:
                molecules.append(Mol(rdmol))
            mol2_block = []
        mol2_block.append(line)
    if mol2_block:
        rdmol = Chem.MolFromMol2Block("".join(mol2_block))
        if rdmol is not None:
            molecules.append(Mol(rdmol))

    file_open.close()
    return molecules

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

    For example, lets assume you wanted to find all of the
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
    from chemper.optimize_smirks.environment import ChemicalEnvironment

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



