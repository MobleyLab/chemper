"""
utils.py

This file provides simple functions that might be of use to people using the chemper package

"""

import os
from chemper.mol_toolkits import mol_toolkit

def get_filename(relative_path, package='chemper'):
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
        return oe_mols_from_file(mol2_file)

    elif toolkit.lower() == 'rdk':
        return rdk_mols_from_mol2(mol2_file)

def oe_mols_from_file(mol_file):
    """
    Parses a standard molecule file into chemper molecules using OpenEye toolkits

    Parameters
    ----------
    mol_file: str
               relative or full path to molecule containing file

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
        raise IOError("File '%s' not found." % mol2_file)

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
               relative or full path to a mol2 file you want to parse

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


