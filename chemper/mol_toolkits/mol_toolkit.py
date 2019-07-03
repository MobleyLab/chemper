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

# identify which toolkits are available
HAS_OE = False
HAS_RDK = False
try:
    from openeye import oechem
    if oechem.OEChemIsLicensed():
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


# ======================================================================
# Find super Mol/Atom/Bond classes
# ======================================================================

if not HAS_OE and not HAS_RDK:
    raise ImportWarning("No Cheminformatics toolkit was found." \
                        " currently ChemPer supports OpenEye and RDKit")


class Mol:
    def __init__(self, mol):
        # check if its a ChemPer Mol with OE wrapper
        if HAS_OE and isinstance(mol, cp_openeye.Mol):
            self.mol = mol.mol
            self.__class__ = cp_openeye.Mol

        # check if this is an Openeye molecule
        elif HAS_OE and isinstance(mol, oechem.OEMolBase):
            self.__class__ = cp_openeye.Mol
            self.__class__.__init__(self,mol)

        # check if its a ChemPer Mol with RDK wrapper
        elif HAS_RDK and isinstance(mol, cp_rdk.Mol):
            self.mol = mol.mol
            self.__class__ = cp_rdk.Mol

        # check if it is an RDK molecule
        elif HAS_RDK and isinstance(mol, Chem.rdchem.Mol):
            self.__class__ = cp_rdk.Mol
            self.__class__.__init__(self, mol)

        else:
            err_msg = """
            Your molecule has the type %s.
            Currently ChemPer only supports OpenEye and RDKit.
            To add support to a new toolkit submit an issue on GitHub at
            github.com/MobleyLab/chemper
            """
            raise TypeError(err_msg % type(mol))

    @staticmethod
    def from_smiles(smiles):
        if HAS_OE:
            return cp_openeye.Mol.from_smiles(smiles)
        return cp_rdk.Mol.from_smiles(smiles)


class Atom:
    def __init__(self, atom):

        if HAS_OE and isinstance(atom, cp_openeye.Atom):
            self.atom = atom.atom
            self.__class__ = cp_openeye.Atom

        elif HAS_OE and isinstance(atom, oechem.OEAtomBase):
            self.__class__ = cp_openeye.Atom
            self.__class__.__init__(self, atom)

        elif HAS_RDK and isinstance(atom, cp_rdk.Atom):
            self.atom = atom.atom
            self.__class__ = cp_rdk.Atom

        elif HAS_RDK and isinstance(atom, Chem.rdchem.Atom):
            self.__class__ = cp_rdk.Atom
            self.__class__.__init__(self, atom)

        else:
            err_msg = """
            Your atom has the type %s.
            Currently ChemPer only supports OpenEye and RDKit.
            To add support to a new toolkit submit an issue on GitHub at
            github.com/MobleyLab/chemper
            """
            raise TypeError(err_msg % type(atom))


class Bond:
    def __init__(self, bond):
        if HAS_OE and isinstance(bond, cp_openeye.Bond):
            self.bond = bond.bond
            self.__class__ = cp_openeye.Bond

        elif HAS_OE and isinstance(bond, oechem.OEBondBase):
            self.__class__ = cp_openeye.Bond
            self.__class__.__init__(self,bond)

        elif HAS_RDK and isinstance(bond, cp_rdk.Bond):
            self.__class__ = cp_rdk.Bond
            self.bond = bond.bond

        elif HAS_RDK and isinstance(bond, Chem.rdchem.Bond):
            self.__class__ = cp_rdk.Bond
            self.__class__.__init__(self,bond)

        else:
            err_msg = """
            Your bond has the type %s.
            Currently ChemPer only supports OpenEye and RDKit.
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
    toolkit : str or None
              'openeye', 'rdkit', or None
              if None then the toolkit will be picked automatically

    Returns
    -------
    toolkit : str
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
        raise ImportError("The provided toolkit (%s) is not supported,"\
                          " ChemPer only supports 'openeye' and 'rdkit'" \
                          % toolkit)


def check_mol_file(file_name):
    """

    Parameters
    ----------
    file_name : str
                path to a molecule file

    Returns
    -------
    path : str
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


# =======================================
# get molecules from files
# =======================================

def mols_from_mol2(mol2_file, toolkit=None):
    """
    Creates a list of ChemPer Mols from the provided mol2 file
    using a specified or default toolkit

    Parameters
    ----------
    mol2_file : str
                path to mol2 file, this can be a relative or absolute path locally
                or the path to a molecules file stored in ChemPer at chemper/data/molecules/
    toolkit : None or str
              'openeye' or 'rdkit' are the two supported toolkits
              if None then the first package available (in the order listed here)
              will be used

    Returns
    -------
    mol2s : list[ChemPer Mol]
            List of molecules in the provided mol2 file
    """
    toolkit = check_toolkit(toolkit)
    mol2_path = check_mol_file(mol2_file)

    if toolkit.lower() == 'openeye':
        return cp_openeye.mols_from_mol2(mol2_path)

    return cp_rdk.mols_from_mol2(mol2_path)
