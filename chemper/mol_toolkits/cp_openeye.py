"""
cp_openeye.py

Cheminformatics tools using OpenEye Toolkits

The classes provided here follow the structure in adapters.
This is a wrapper allowing our actual package to use openeye toolkits

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

from chemper.mol_toolkits.adapters import MolAdapter, AtomAdapter, BondAdapter
from openeye import oechem


# =======================================
# Molecule Class
# =======================================

class Mol(MolAdapter):
    """
    Wrapper for OEMol to create a chemper Mol
    """
    def __init__(self, mol):
        """
        Parameters
        ----------
        mol: openeye OEMol object
            openeye molecule to convert to chemper Mol object
        """
        if type(mol) != oechem.OEMol:
            raise Exception("Expecting an OEMol object instead of %s" % type(mol))
        self.mol = mol

    def __str__(self): return self.get_smiles()

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
        return Atom(self.mol.GetAtom(oechem.OEHasAtomIdx(idx)))

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
        return Bond(self.mol.GetBond(oechem.OEHasBondIdx(idx)))

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
        return Bond(self.mol.GetBond(atom1.atom, atom2.atom))

    def smirks_search(self, smirks):
        """
        Performs a substructure search on the molecule with the provided
        SMIRKS pattern. Note - this function expects SMIRKS patterns with indexed atoms
        that is with :n for at least some atoms.

        Parameters
        ----------
        smirks: str
            SMIRKS pattern with indexed atoms (:n)

        Returns
        -------
        matches: list of dictionaries
            dictionary for each match with the form {smirks index: atom index}
        """
        matches = list()

        ss = oechem.OESubSearch()
        if not ss.Init(smirks):
            raise ValueError("Error parsing SMIRKS %s" % smirks)

        for match in ss.Match(self.mol, False):
            d = dict()
            for ma in match.GetAtoms():
                smirks_idx = ma.pattern.GetMapIdx()
                # if the map index is 0 then it isn't a "tagged" atom in the SMIRKS
                if smirks_idx !=0:
                    d[smirks_idx] = self.get_atom_by_index(ma.target.GetIdx())

            matches.append(d)

        return matches

    def get_smiles(self):
        """
        Returns
        -------
        smiles: str
            SMILES string for the molecule
        """
        smiles = oechem.OEMolToSmiles(self.mol)
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
        mol = oechem.OEMol()
        if not oechem.OESmilesToMol(mol, smiles):
            raise ValueError('Could not parse SMILES %s' % smiles)
        oechem.OEAddExplicitHydrogens(mol)
        Mol.__init__(self, mol)

# =======================================
# Atom Class
# =======================================

class Atom(AtomAdapter):
    """
    Wrapper for OEAtomBase to create a chemper Atom
    """
    def __init__(self, atom):
        """
        Parameters
        ----------
        atom: OEAtomBase
            Atom object from an OpenEye molecule
        """
        if type(atom) != oechem.OEAtomBase:
            raise Exception("Expecting an OEAtomBase object instead of %s" % type(atom))
        self.atom = atom
        self._index = self.atom.GetIdx()

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
        return len([b for b in self.atom.GetBonds()])

    def valence(self):
        """
        Returns
        -------
        valence: int
            the atoms valence
        """
        return self.atom.GetValence()

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
        return self.atom.GetTotalHCount()

    def ring_connectivity(self):
        """
        Returns
        -------
        ring_connectivity: int
            number of bonds on the atom that are a part of a ring
        """
        return len([b for b in self.atom.GetBonds(oechem.OEBondIsInRing())])

    def min_ring_size(self):
        """
        Returns
        -------
        min_ring_size: int
            size of the smallest ring this atom is a part of
        """
        return oechem.OEAtomGetSmallestRingSize(self.atom)

    def is_aromatic(self):
        """
        Returns
        -------
        is_aromatic: boolean
            True if the atom is aromatic otherwise False
        """
        return self.atom.IsAromatic()

    def get_index(self):
        """
        Returns
        -------
        index: int
            atom index in its molecule
        """
        return self._index

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
        return self.atom.IsConnected(atom2.atom)

    def get_neighbors(self):
        """
        Returns
        -------
        neighbors: list of chemper Atoms
            atoms that are one bond away from this atom
        """
        return [Atom(a) for a in self.atom.GetAtoms()]

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
        mol = oechem.OEMol(self.atom.GetParent())
        self.atom = mol.GetAtom(oechem.OEHasAtomIdx(self._index))
        return Mol(mol)

# =======================================
# Bond Class
# =======================================


class Bond(BondAdapter):
    """
    Wrapper for OEBondBase to create a chemper Bond
    """
    def __init__(self, bond):
        """
        Parameters
        ----------
        bond: OEBondBase
            Bond object from an OpenEye molecule
        """
        if type(bond) != oechem.OEBondBase:
            raise Exception("Expecting an OEBondBase object instead of %s" % type(bond))
        self.bond = bond
        self._order = self.bond.GetOrder()
        if self.is_aromatic():
            self._order = 1.5

        self._idx = self.bond.GetIdx()

    def get_order(self):
        """
        Returns
        -------
        order: int or float
            This is the absolute order, returns 1.5 if bond is aromatic
        """
        return self._order

    def get_atoms(self):
        """
        Returns
        -------
        atoms: list of chemper Atoms
            the two atoms connected by this bond
        """
        beginning = Atom(self.bond.GetBgn())
        end = Atom(self.bond.GetEnd())
        return [beginning, end]

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
        return self.bond.IsAromatic()

    def is_single(self):
        """
        Returns
        -------
        is_single: boolean
            True if it is a single bond
        """
        return self._order == 1

    def is_double(self):
        """
        Returns
        -------
        is_double: boolean
            True if it is a double bond
        """
        return self._order == 2

    def is_triple(self):
        """
        Returns
        -------
        is_triple: boolean
            True if it is a triple bond
        """
        return self._order == 3

    def get_molecule(self):
        """
        Extracts the parent molecule this bond is in

        Returns
        -------
        mol: chemper Mol
            molecule this bond is stored in
        """
        mol = oechem.OEMol(self.bond.GetParent())
        self.bond = mol.GetBond(oechem.OEHasBondIdx(self._idx))
        return Mol(mol)

    def get_index(self):
        """
        Returns
        -------
        index: int
            index of this bond in its parent molecule
        """
        return self._idx
