"""
cp_rdk.py

Cheminformatics tools using RDKit

The classes provided here follow the structure in adapters.
This is a wrapper allowing our actual package to use RDKit

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

from chemper.mol_toolkits.adapters import MolAdapter, AtomAdapter, BondAdapter
from rdkit import Chem


# =======================================
# Molecule Class
# =======================================

class Mol(MolAdapter):
    def __init__(self, mol):
        """
        Create a ChemPer Mol from an RDMol

        Parameters
        ----------
        mol : openeye RDKMol object
              openeye molecule to convert to ChemPer Mol object
        """
        if not isinstance(mol, Chem.rdchem.Mol):
            raise TypeError("Expecting an rdchem.Mol instead of %s" % type(mol))
        self.mol = mol

    def __str__(self): return self.get_smiles()

    @classmethod
    def from_smiles(cls, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError('Could not parse SMILES %s' % smiles)
        cls(Chem.AddHs(mol))

    def set_aromaticity_mdl(self):
        Chem.SanitizeMol(self.mol, Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(self.mol, Chem.AromaticityModel.AROMATICITY_MDL)

    def get_atoms(self):
        return [Atom(a) for a in self.mol.GetAtoms()]

    def get_atom_by_index(self, idx):
        return Atom(self.mol.GetAtomWithIdx(idx))

    def get_bonds(self):
        return [Bond(b) for b in self.mol.GetBonds()]

    def get_bond_by_index(self, idx):
        return Bond(self.mol.GetBondWithIdx(idx))

    def get_bond_by_atoms(self, atom1, atom2):
        if not atom1.is_connected_to(atom2):
            return None
        return Bond(self.mol.GetBondBetweenAtoms(atom1.get_index(), atom2.get_index()))

    def smirks_search(self, smirks):
        cmol = Chem.Mol(self.mol)

        matches = list()

        ss = Chem.MolFromSmarts(smirks)
        if ss is None:
            raise ValueError("Error parsing SMIRKS %s" % smirks)

        # get atoms in query mol with smirks index
        maps = dict()
        for qatom in ss.GetAtoms():
            smirks_idx = qatom.GetAtomMapNum()
            if smirks_idx != 0:
                maps[smirks_idx] = qatom.GetIdx()

        for match in cmol.GetSubstructMatches(ss, False):
            d = {k:self.get_atom_by_index(match[e]) for k,e in maps.items()}
            matches.append(d)

        return matches

    def get_smiles(self):
        smiles = Chem.MolToSmiles(Chem.RemoveHs(self.mol))
        return smiles


# =======================================
# Atom Class
# =======================================


class Atom(AtomAdapter):
    def __init__(self, atom):
        """
        Create a ChemPer Atom from an RDAtom

        Parameters
        ----------
        atom : RDKit Atom
               Atom object from an RDK molecule
        """
        if not isinstance(atom, Chem.rdchem.Atom):
            raise TypeError("Expecting a rdchem.Atom instead of %s" % type(atom))
        self.atom = atom
        self._idx = self.atom.GetIdx()

    def __str__(self): return '%i%s' % (self._idx, self.atom.GetSymbol())

    def atomic_number(self): return self.atom.GetAtomicNum()

    def degree(self): return self.atom.GetDegree()

    def connectivity(self): return self.atom.GetTotalDegree()

    def valence(self): return self.atom.GetTotalValence()

    def formal_charge(self): return self.atom.GetFormalCharge()

    def hydrogen_count(self):
        return self.atom.GetTotalNumHs(includeNeighbors=True)

    def ring_connectivity(self):
        return len([b for b in self.atom.GetBonds() if b.IsInRing()])

    def min_ring_size(self):
        if not self.atom.IsInRing():
            return 0
        for i in range(10000):
            if self.atom.IsInRingSize(i):
                return i
        # TODO: raise exception instead?
        return 10000

    def is_aromatic(self): return self.atom.GetIsAromatic()

    def get_index(self): return self._idx

    def is_connected_to(self, atom2):
        if not isinstance(atom2.atom, Chem.rdchem.Atom):
            return False
        neighbors = [a.GetIdx() for a in self.atom.GetNeighbors()]
        return atom2.get_index() in neighbors

    def get_neighbors(self):
        return [Atom(a) for a in self.atom.GetNeighbors()]

    def get_bonds(self):
        return [Bond(b) for b in self.atom.GetBonds()]

    def get_molecule(self):
        mol = Chem.Mol(self.atom.GetOwningMol())
        return Mol(mol)

# =======================================
# Bond Class
# =======================================


class Bond(BondAdapter):
    def __init__(self, bond):
        """
        Creates a ChemPer Bond from an RDK Bond

        Parameters
        ----------
        bond : RDK Bond
               Bond object from an RDK molecule
        """
        if not isinstance(bond, Chem.rdchem.Bond):
            raise TypeError("Expecting an rdchem.Bond instead of %s" % type(bond))
        self.bond = bond

        # save index
        self._idx = self.bond.GetIdx()

        # save order information
        self._order = self.bond.GetBondTypeAsDouble()
        orders = {1:'-', 2:'=', 3:'#', 1.5:':'}
        self._order_symbol = orders.get(self._order, '~')

        # save atoms in bond
        self._beginning = Atom(self.bond.GetBeginAtom())
        self._end = Atom(self.bond.GetEndAtom())

    def __str__(self):
        return "%i %s%s%s" % (self.get_index(), self._beginning,
                              self._order_symbol, self._end)

    def get_order(self): return self._order

    def get_atoms(self): return [self._beginning, self._end]

    def is_ring(self): return self.bond.IsInRing()

    def is_aromatic(self): return self.bond.GetIsAromatic()

    def is_single(self): return self._order == 1

    def is_double(self): return self._order == 2

    def is_triple(self): return self._order == 3

    def get_molecule(self):
        mol = Chem.Mol(self.bond.GetOwningMol())
        return Mol(mol)

    def get_index(self): return self._idx

# =====================================================================
# functions for importing molecules from files
# =====================================================================

def mols_from_mol2(mol2_file):
    """
    Parses a mol2 file into ChemPer molecules using RDKit

    This is a hack for separating mol2 files taken from a Source Forge discussion here:
    https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01510.html
    It splits up a mol2 file into blocks and then uses RDKit to parse those blocks

    Parameters
    ----------
    mol2_file : str
                relative or absolute path to a mol2 file you want to parse
                accessible form the current directory

    Returns
    -------
    mols : list[ChemPer Mol]
           list of molecules in the mol2 file as ChemPer molecules
    """
    # TODO: check that this works with mol2 files with a single molecule
    # TODO: figure out if @<TRIPOS>MOLECULE is the only delimiter acceptable in this file type
    import os

    if not os.path.exists(mol2_file):
        from chemper.chemper_utils import get_data_path
        mol_path = get_data_path(os.path.join('molecules', mol2_file))

        if not os.path.exists(mol_path):
            raise IOError("File '%s' not found locally or in chemper/data/molecules." % mol2_file)
        else:
            mol2_file = mol_path

    delimiter="@<TRIPOS>MOLECULE"

    if mol2_file.split('.')[-1] != "mol2":
        raise IOError("File '%s' is not a mol2 file" % mol2_file)

    if not os.path.exists(mol2_file):
        raise IOError("File '%s' not found." % mol2_file)

    molecules = list()
    mol2_block = list()

    file_open = open(mol2_file)

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

