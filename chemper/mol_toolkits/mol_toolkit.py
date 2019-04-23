"""
mol_toolkit is a set of factory methods which
allow the user to interact with ChemPer Mol, Atom, and Bond
objects without having to explicitly define which toolkit they
are using.

For example one user could be loading their molecules through openeye
and the others could be using RDKit to create what we would call a "user mol"

In this case they can call:
from chemper.mol_toolkits import mol_toolkit
mol = mol_toolkit.Mol(user_mol)
without needing to tell us what toolkit they are using or they
could specify
mol = mol_toolkit.Mol(user_mol, toolkit='openeye')

"""
import os
from chemper.chemper_utils import get_data_path
from chemper.mol_toolkits.adapters import MolAdapter, BondAdapter, AtomAdapter
# define default_toolkit:
try:
    from openeye.oechem import OEChemIsLicensed
    if not OEChemIsLicensed():
        HAS_OE = False
    else:
        from chemper.mol_toolkits import cp_openeye
        HAS_OE = True
except ImportError:
    HAS_OE = False

try:
    from rdkit import Chem
    from chemper.mol_toolkits import cp_rdk
    HAS_RDK = True
except ImportError:
    HAS_RDK = False

if not HAS_OE and not HAS_RDK:
    raise ImportWarning("No Cheminformatics toolkit was found."\
                        " currently chemper supports OpenEye and RDKit")


def Mol(mol):
    """

    Parameters
    ----------
    mol - a molecule object from any supported toolkit

    Returns
    -------
    mol - a chemper Mol

    """
    # if it is already a chemper molecule return as is
    if isinstance(mol, MolAdapter):
        return mol

    # check if this is an Openeye molecule
    if HAS_OE and isinstance(mol, oechem.OEMolBase):
        return cp_openeye.Mol(mol)

    # check if it is an RDK molecule
    if HAS_RDK and isinstance(mol, Chem.rdchem.Mol):
        return cp_rdk.Mol(mol)

    err_msg = """
    Your molecule has the type %s.
    Currently chemper only supports OpenEye and RDKit.
    To add support to a new toolkit submit an issue on GitHub at
    github.com/MobleyLab/chemper
    """
    raise TypeError(err_msg % type(mol))


def Atom(atom):
    """

    Parameters
    ----------
    atom - Atom object from any supported toolkit

    Returns
    -------
    atom - a chemper Atom object

    """
    if isinstance(atom, AtomAdapter):
        return atom

    if HAS_OE and isinstance(atom, oechem.OEAtomBase):
        return cp_openeye.Atom(atom)

    if HAS_RDK and isinstance(atom, Chem.rdchem.Atom):
        return cp_rdk.Atom(atom)

    err_msg = """
    Your atom has the type %s.
    Currently chemper only supports OpenEye and RDKit.
    To add support to a new toolkit submit an issue on GitHub at
    github.com/MobleyLab/chemper
    """
    raise TypeError(err_msg % type(atom))


def Bond(bond):
    """

    Parameters
    ----------
    bond - Bond object from any supported toolkit

    Returns
    -------
    bond - a chemper Bond object

    """
    if isinstance(bond, BondAdapter):
        return bond

    if HAS_OE and isinstance(bond, oechem.OEBondBase):
        return cp_openeye.Bond(bond)

    if HAS_RDK and isinstance(bond, Chem.rdchem.Bond):
        return cp_rdk.Bond(bond)

    err_msg = """
    Your bond has the type %s.
    Currently chemper only supports OpenEye and RDKit.
    To add support to a new toolkit submit an issue on GitHub at
    github.com/MobleyLab/chemper
    """
    raise TypeError(err_msg % type(bond))

# =======================================
# check user specifications
# =======================================

def check_toolkit(toolkit=None):
    """

    Parameters
    ----------
    toolkit - str
        'openeye', 'rdkit', or None
        if None then the toolkit will be picked automatically

    Returns
    -------
    toolkit - str
    returns the name of the toolkit to be used.
    If the package isn't available for the specified toolkit
    then an error is raised instead
    """
    # check for a stable
    if toolkit is None:
        if HAS_OE:
            return 'openeye'
        elif HAS_RDK:
            return 'rdkit'

    if toolkit.lower() == 'openeye' and HAS_OE:
        return 'openeye'

    if toolkit.lower() == 'rdkit' and HAS_RDK:
        return 'rdkit'

    if toolkit.lower() == 'openeye' or toolkit.lower() == 'rdkit':
        raise ImportError("Toolkit (%s) was not importable" % toolkit)

    else:
        raise ImportError("The provided toolkit (%s) isn't supported,"\
                          " ChemPer only supports 'openeye' and 'rdkit'" \
                          % toolkit)

def check_mol_file(file_name):
    """

    Parameters
    ----------
    file_name - str
        path to a molecule file

    Returns
    -------
    path - str
        absolute path to a molecule file
        raises error if file isn't available
    """
    # is it a local file?
    if os.path.exists(file_name):
        return os.path.abspath(file_name)

    path = get_data_path(os.path.join('molecules', file_name))

    if not os.path.exists(path):
        raise IOError("Molecule file (%s) was not found locally or in chemper/data/molecules" % file_name)

    return path


# ================================================================
# get molecule from SMILES
# ================================================================

def MolFromSmiles(smiles, toolkit=None):
    """

    Parameters
    ----------
    smiles - str
        SMILES string
    toolkit - str or None
        'openeye' or 'rdkit' or None to let chemper pick

    Returns
    -------
    mol - ChemPer Mol
    """
    toolkit = check_toolkit(toolkit)
    if toolkit.lower() == 'openeye':
        return cp_openeye.MolFromSmiles(smiles)

    return cp_rdk.MolFromSmiles(smiles)

# =======================================
# get molecules from files
# =======================================

def mols_from_mol2(mol2_file, toolkit=None):
    """
    Creates a list of chemper Mols from the provided mol2 file
    using a specified or default toolkit

    Parameters
    ----------
    mol2_file: str
               path to mol2 file, this can be a relative or absolute path locally
               or the path to a molecules file stored in chemper at chemper/data/molecules/
    toolkit: None or str
             'openeye' or 'rdkit' are the two supported toolkits
             if None then the first package available (in the order listed here)
             will be used

    Returns
    -------
    mol2s: list of chemper Mol
           list of molecules in the provided mol2 file
    """
    toolkit = check_toolkit(toolkit)
    mol2_path = check_mol_file(mol2_file)

    if toolkit.lower() == 'openeye':
        return cp_openeye.mols_from_mol2(mol2_path)

    return cp_rdk.mols_from_mol2(mol2_path)



