"""
adapters.py

This script contains adapters or the structure for atoms, molecules and
substructure searches.
Our chemical perception code is designed to be independent of the users
cheminformatics packages. For each cheminformatics package we support we
will provide classes following the structure in these adapters.

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

from abc import ABC, abstractmethod

# =======================================
# Molecule Class
# =======================================

class MolAdapter(ABC):
    """
    Template class for wrapping a molecule object
    from a given cheminformatics package into a chemper Mol.
    chemper `Mol` are initiated from the reference package molecule object.
    The class MolFromSmiles initiates a chemper Mol from a SMILES string.
    Currently we support OpenEye toolkits and RDKit
    """
    @abstractmethod
    def get_atoms(self):
        """
        Returns
        -------
        atom_list: list of chemper Atoms
            list of all atoms in the molecule
        """
        return

    @abstractmethod
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
        return

    @abstractmethod
    def get_bonds(self):
        """
        Returns
        -------
        bond_list: list of chemper Bonds
            list of all bonds in the molecule
        """
        return

    @abstractmethod
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
        return

    @abstractmethod
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
        return

    @abstractmethod
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
        return

    @abstractmethod
    def get_smiles(self):
        """
        Returns
        -------
        smiles: str
            SMILES string for the molecule
        """
        return

# =======================================
# Atom Class
# =======================================


class AtomAdapter(ABC):
    """
    Template class for wrapping an atom object
    from a given cheminformatics package into a chemper Atom.
    These are always initiated from a reference package Atom object.
    Currently we support OpenEye toolkits and RDKit
    """

    @abstractmethod
    def atomic_number(self):
        """
        Returns
        -------
        atomic_number: int
            atomic number for the atom
        """
        return

    @abstractmethod
    def degree(self):
        """
        Returns
        -------
        degree: int
            degree or number of explicit bonds around the atom
        """
        return

    @abstractmethod
    def connectivity(self):
        """
        Returns
        -------
        connectivity: int
            connectivity or total number of bonds around the atom
        """
        return

    @abstractmethod
    def valence(self):
        """
        Returns
        -------
        valence: int
            the atoms valence
        """
        return

    @abstractmethod
    def formal_charge(self):
        """
        Returns
        -------
        formal_charge: int
            the atom's formal charge
        """
        return

    @abstractmethod
    def hydrogen_count(self):
        """
        Returns
        -------
        H_count: int
            total number of hydrogen atoms connected to this Atom
        """
        return

    @abstractmethod
    def min_ring_size(self):
        """
        Returns
        -------
        min_ring_size: int
            size of the smallest ring this atom is a part of
        """
        return

    @abstractmethod
    def ring_connectivity(self):
        """
        Returns
        -------
        ring_connectivity: int
            number of bonds on the atom that are a part of a ring
        """
        return

    @abstractmethod
    def is_aromatic(self):
        """
        Returns
        -------
        is_aromatic: boolean
            True if the atom is aromatic otherwise False
        """
        return

    @abstractmethod
    def get_index(self):
        """
        Returns
        -------
        index: int
            atom index in its molecule
        """
        return

    @abstractmethod
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
        return

    @abstractmethod
    def get_neighbors(self):
        """
        Returns
        -------
        neighbors: list of chemper Atoms
            atoms that are one bond away from this atom
        """
        return

    @abstractmethod
    def get_bonds(self):
        """
        Returns
        -------
        bonds: list of chemper Bonds
            bonds connected to this atom
        """
        return

    @abstractmethod
    def get_molecule(self):
        """
        Extracts the parent molecule this atom is in

        Returns
        -------
        mol: chemper Mol
            molecule this atom is stored in
        """
        return


# =======================================
# Bond Class
# =======================================

class BondAdapter(ABC):
    """
    Template class for wrapping a bond object
    from a given cheminformatics package into a chemper Bond.
    These are always initiated from a reference package Bond object.
    Currently we support OpenEye toolkits and RDKit
    """

    @abstractmethod
    def get_order(self):
        """
        Returns
        -------
        order: int or float
            This is the absolute order, returns 1.5 if bond is aromatic
        """
        return

    @abstractmethod
    def get_atoms(self):
        """
        Returns
        -------
        atoms: list of chemper Atoms
            the two atoms connected by this bond
        """
        return

    @abstractmethod
    def is_ring(self):
        """
        Returns
        -------
        is_ring: boolean
            True if bond is a part of a ring, otherwise False
        """
        return

    @abstractmethod
    def is_aromatic(self):
        """
        Returns
        -------
        is_aromatic: boolean
            True if it is an aromatic bond
        """
        return

    @abstractmethod
    def is_single(self):
        """
        Returns
        -------
        is_single: boolean
            True if it is a single bond
        """
        return

    @abstractmethod
    def is_double(self):
        """
        Returns
        -------
        is_double: boolean
            True if it is a double bond
        """
        return

    @abstractmethod
    def is_triple(self):
        """
        Returns
        -------
        is_triple: boolean
            True if it is a triple bond
        """
        return

    @abstractmethod
    def get_molecule(self):
        """
        Extracts the parent molecule this bond is in

        Returns
        -------
        mol: chemper Mol
            molecule this bond is stored in
        """
        return

    @abstractmethod
    def get_index(self):
        """
        Returns
        -------
        index: int
            index of this bond in its parent molecule
        """
        return
