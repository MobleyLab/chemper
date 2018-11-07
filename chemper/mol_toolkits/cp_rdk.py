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
    """
    Wrapper for RDKMol to create a chemper Mol
    """
    def __init__(self, mol):
        """
        Parameters
        ----------
        mol: openeye RDKMol object
            openeye molecule to convert to chemper Mol object
        """
        if type(mol) != Chem.rdchem.Mol:
            raise Exception("Expecting an rdchem.Mol instead of %s" % type(mol))
        self.mol = mol

    def __str__(self): return self.get_smiles()

    def set_aromaticity_mdl(self):
        """
        Sets the aromaticity flags in this molecule to use the MDL model
        """
        Chem.SanitizeMol(self.mol, Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(self.mol, Chem.AromaticityModel.AROMATICITY_MDL)

    def get_atoms(self):
        """
        Returns
        -------
        atom_list: list of chemper Atoms
            list of all atoms in the molecule
        """
        return [Atom(a) for a in self.mol.GetAtoms()]

    def get_atom_by_index(self, idx):
        """
        Parameters
        ----------
        idx: int
            atom index

        Returns
        -------
        atom: chemper Atom object
            atom with index idx
        """
        return Atom(self.mol.GetAtomWithIdx(idx))

    def get_bonds(self):
        """
        Returns
        -------
        bond_list: list of chemper Bonds
            list of all bonds in the molecule
        """
        return [Bond(b) for b in self.mol.GetBonds()]

    def get_bond_by_index(self, idx):
        """
        Parameters
        ----------
        idx: ing
            bond index

        Returns
        -------
        bond: chemper Bond object
            bond with index idx
        """
        return Bond(self.mol.GetBondWithIdx(idx))

    def get_bond_by_atoms(self, atom1, atom2):
        """
        Finds a bond between two atoms
        Parameters
        ----------
        atom1: chemper Atom object
        atom2: chemper Atom object

        Returns
        -------
        bond: chemper Bond object or None
            if atoms are connected returns bond otherwise None
        """
        if not atom1.is_connected_to(atom2):
            return None
        return Bond(self.mol.GetBondBetweenAtoms(atom1.get_index(), atom2.get_index()))

    def smirks_search(self, smirks, use_mdl=False):
        """
        Performs a substructure search on the molecule with the provided
        SMIRKS pattern. Note - this function expects SMIRKS patterns with indexed atoms
        that is with :n for at least some atoms.

        Parameters
        ----------
        smirks: str
            SMIRKS pattern with indexed atoms (:n)
        use_mdl: boolean (optional)
            Use MDL aromaticity model for this substructure search

        Returns
        -------
        matches: list of dictionaries
            dictionary for each match with the form {smirks index: atom index}
        """
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
        """
        Returns
        -------
        smiles: str
            SMILES string for the molecule
        """
        smiles = Chem.MolToSmiles(Chem.RemoveHs(self.mol))
        return smiles

class MolFromSmiles(Mol):
    """
    Creates a chemper Mol from a smiles string
    It automatically adds explicit hydrogens.
    """
    def __init__(self, smiles):
        """
        Parameters
        ----------
        smiles: str
            SMILES string for a molecule
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError('Could not parse SMILES %s' % smiles)
        Mol.__init__(self, Chem.AddHs(mol))

# =======================================
# Atom Class
# =======================================


class Atom(AtomAdapter):
    """
    Wrapper for RDKAtom to create a chemper Atom
    """
    def __init__(self, atom):
        """
        Parameters
        ----------
        atom: RDKAtom
            Atom object from an RDK molecule
        """
        if type(atom) != Chem.rdchem.Atom:
            raise Exception("Expecting an rdchem.Atom instead of %s" % type(atom))
        self.atom = atom

    def atomic_number(self):
        """
        Returns
        -------
        atomic_number: int
            atomic number for the atom
        """
        return self.atom.GetAtomicNum()

    def degree(self):
        """
        Returns
        -------
        degree: int
            degree or number of explicit bonds around the atom
        """
        return self.atom.GetDegree()

    def connectivity(self):
        """
        Returns
        -------
        connectivity: int
            connectivity or total number of bonds around the atom
        """
        return self.atom.GetTotalDegree()

    def valence(self):
        """
        Returns
        -------
        valence: int
            the atoms valence
        """
        return self.atom.GetTotalValence()

    def formal_charge(self):
        """
        Returns
        -------
        formal_charge: int
            the atom's formal charge
        """
        return self.atom.GetFormalCharge()

    def hydrogen_count(self):
        """
        Returns
        -------
        H_count: int
            total number of hydrogen atoms connected to this Atom
        """
        return self.atom.GetTotalNumHs(includeNeighbors=True)

    def ring_connectivity(self):
        """
        Returns
        -------
        ring_connectivity: int
            number of bonds on the atom that are a part of a ring
        """
        return len([b for b in self.atom.GetBonds() if b.IsInRing()])

    def min_ring_size(self):
        """
        Returns
        -------
        min_ring_size: int
            size of the smallest ring this atom is a part of
        """
        if not self.atom.IsInRing():
            return 0
        for i in range(10000):
            if self.atom.IsInRingSize(i):
                return i
        # TODO: raise exception instead?
        return 10000

    def is_aromatic(self):
        """
        Returns
        -------
        is_aromatic: boolean
            True if the atom is aromatic otherwise False
        """
        return self.atom.GetIsAromatic()

    def get_index(self):
        """
        Returns
        -------
        index: int
            atom index in its molecule
        """
        return self.atom.GetIdx()

    def is_connected_to(self, atom2):
        """
        Parameters
        ----------
        atom2: chemper Atom object
            atom to check if it is connected to this atom

        Returns
        -------
        connected: boolean
            True if atom2 is a direct neighbor or atom1
        """
        if not type(atom2.atom) is Chem.rdchem.Atom:
            # TODO: raise exception/return something else?
            return False
        neighbors = [a.GetIdx() for a in self.atom.GetNeighbors()]
        return atom2.get_index() in neighbors

    def get_neighbors(self):
        """
        Returns
        -------
        neighbors: list of chemper Atoms
            atoms that are one bond away from this atom
        """
        return [Atom(a) for a in self.atom.GetNeighbors()]

    def get_bonds(self):
        """
        Returns
        -------
        bonds: list of chemper Bonds
            bonds connected to this atom
        """
        return [Bond(b) for b in self.atom.GetBonds()]

    def get_molecule(self):
        """
        Extracts the parent molecule this atom is in

        Returns
        -------
        mol: chemper Mol
            molecule this atom is stored in
        """
        mol = Chem.Mol(self.atom.GetOwningMol())
        return Mol(mol)

# =======================================
# Bond Class
# =======================================


class Bond(BondAdapter):
    """
    Wrapper for RDKBond to create a chemper Bond
    """
    def __init__(self, bond):
        """
        Parameters
        ----------
        bond: RDKBond
            Bond object from an RDK molecule
        """
        if type(bond) != Chem.rdchem.Bond:
            raise Exception("Expecting an rdchem.Bond instead of %s" % type(bond))
        self.bond = bond
        self.order = self.bond.GetBondTypeAsDouble()
        self.beginning = Atom(self.bond.GetBeginAtom())
        self.end = Atom(self.bond.GetEndAtom())

    def get_order(self):
        """
        Returns
        -------
        order: int or float
            This is the absolute order, returns 1.5 if bond is aromatic
        """
        return self.order

    def get_atoms(self):
        """
        Returns
        -------
        atoms: list of chemper Atoms
            the two atoms connected by this bond
        """
        return [self.beginning, self.end]

    def is_ring(self):
        """
        Returns
        -------
        is_ring: boolean
            True if bond is a part of a ring, otherwise False
        """
        return self.bond.IsInRing()

    def is_aromatic(self):
        """
        Returns
        -------
        is_aromatic: boolean
            True if it is an aromatic bond
        """
        return self.bond.GetIsAromatic()

    def is_single(self):
        """
        Returns
        -------
        is_single: boolean
            True if it is a single bond
        """
        return self.order == 1

    def is_double(self):
        """
        Returns
        -------
        is_double: boolean
            True if it is a double bond
        """
        return self.order == 2

    def is_triple(self):
        """
        Returns
        -------
        is_triple: boolean
            True if it is a triple bond
        """
        return self.order == 3

    def get_molecule(self):
        """
        Extracts the parent molecule this bond is in

        Returns
        -------
        mol: chemper Mol
            molecule this bond is stored in
        """
        mol = Chem.Mol(self.bond.GetOwningMol())
        return Mol(mol)

    def get_index(self):
        """
        Returns
        -------
        index: int
            index of this bond in its parent molecule
        """
        return self.bond.GetIdx()

def mols_from_mol2(mol2_file):
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
    import os

    if not os.path.exists(mol2_file):
        from chemper.chemper_utils import get_data_path
        mol_path = get_data_path(os.path.join('molecules', mol2_file))

        if not os.path.exists(mol_path):
            raise IOError("File '%s' not found locally or in chemper/data/molecules." % mol_file)
        else:
            mol2_file = mol_path

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

