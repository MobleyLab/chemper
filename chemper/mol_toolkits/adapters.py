"""
adapters.py

This script contains adapters or the structure for
molecules, atoms, and bonds.
Our chemical perception code is designed to be independent of the users
cheminformatics packages. For each cheminformatics package we support we
will provide classes following the structure in these adapters.

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

from abc import ABC, abstractmethod, abstractclassmethod


# =======================================
# Molecule Class
# =======================================

class MolAdapter(ABC):
    """
    This is a ChemPer wrapper for a molecule from
    one of the cheminformatics toolkits.
    ChemPer molecules are initiated from the reference package molecule object.
    Currently we support OpenEye and RDKit toolkits.

    Attributes
    ----------
    mol : toolkit Mol
          Mol object from the reference cheminformatics toolkit
    """
    @abstractclassmethod
    def from_smiles(cls, smiles):
        """
        Creates a ChemPer Mol form a SMILES string

        Parameters
        ----------
        smiles : str
                 SMILES used to create molecule with wrapped toolkit

        Returns
        -------
        Mol : ChemPer Mol
        """
    @abstractmethod
    def set_aromaticity_mdl(self):
        """
        Sets the aromaticity flags in this molecule to use the MDL model
        """
        pass

    @abstractmethod
    def get_atoms(self):
        """
        Returns
        -------
        atom_list : list[ChemPer Atoms]
            list of all atoms in the molecule
        """
        pass

    @abstractmethod
    def get_atom_by_index(self, idx):
        """
        Parameters
        ----------
        idx : int
            atom index

        Returns
        -------
        atom : ChemPer Atom
            atom with index idx
        """
        pass

    @abstractmethod
    def get_bonds(self):
        """
        Returns
        -------
        bond_list : list[ChemPer Bonds]
            list of all bonds in the molecule
        """
        pass

    @abstractmethod
    def get_bond_by_index(self, idx):
        """
        Parameters
        ----------
        idx: int
            bond index

        Returns
        -------
        bond: ChemPer Bond
            bond with index idx
        """
        pass

    @abstractmethod
    def get_bond_by_atoms(self, atom1, atom2):
        """
        Finds a bond between two atoms

        Parameters
        ----------
        atom1 : ChemPer Atom
        atom2 : ChemPer Atom

        Returns
        -------
        bond : ChemPer Bond or None
            If atoms are connected returns bond otherwise None
        """
        pass

    @abstractmethod
    def smirks_search(self, smirks):
        """
        Performs a substructure search on the molecule with the provided
        SMIRKS pattern. Note - this function expects SMIRKS patterns with indexed atoms
        that is with :n for at least some atoms.

        Parameters
        ----------
        smirks : str
            SMIRKS pattern with indexed atoms (:n)

        Returns
        -------
        matches : list[match dictionary]
            match dictionaries have the form {smirks index: atom index}
        """
        pass

    @abstractmethod
    def get_smiles(self):
        """
        Returns
        -------
        smiles: str
            SMILES string for the molecule
        """
        pass

# =======================================
# Atom Class
# =======================================


class AtomAdapter(ABC):
    """
    This is a ChemPer wrapper for an atom from
    one of the cheminformatics toolkits.
    ChemPer Atoms are initiated from the reference package object.
    Currently we support OpenEye and RDKit toolkits.

    Attributes
    ----------
    atom : Atom from reference toolkit
    """
    @abstractmethod
    def atomic_number(self):
        """
        Returns
        -------
        atomic_number : int
            atomic number for the atom
        """
        pass

    @abstractmethod
    def degree(self):
        """
        Returns
        -------
        degree : int
            degree or number of explicit bond orders around the atom
        """
        pass

    @abstractmethod
    def connectivity(self):
        """
        Returns
        -------
        connectivity : int
            connectivity or total number of bonds (regardless of order) around the atom
        """
        pass

    @abstractmethod
    def valence(self):
        """
        Returns
        -------
        valence : int
            the atoms valence (equivalent to degree when all bonds are explicit)
        """
        pass

    @abstractmethod
    def formal_charge(self):
        """
        Returns
        -------
        formal_charge : int
            the atom's formal charge
        """
        pass

    @abstractmethod
    def hydrogen_count(self):
        """
        Returns
        -------
        H_count : int
            total number of hydrogen atoms connected to this Atom
        """
        pass

    @abstractmethod
    def min_ring_size(self):
        """
        Returns
        -------
        min_ring_size : int
            size of the smallest ring this atom is a part of
        """
        pass

    @abstractmethod
    def ring_connectivity(self):
        """
        Returns
        -------
        ring_connectivity : int
            number of bonds on the atom that are a part of a ring
        """
        pass

    @abstractmethod
    def is_aromatic(self):
        """
        Returns
        -------
        is_aromatic : boolean
            True if the atom is aromatic otherwise False
        """
        pass

    @abstractmethod
    def get_index(self):
        """
        Returns
        -------
        index : int
            atom index in its molecule
        """
        pass

    @abstractmethod
    def is_connected_to(self, atom2):
        """
        Parameters
        ----------
        atom2 : ChemPer Atom
            Atom to check if it is bonded to this atom

        Returns
        -------
        connected : boolean
            True if atom2 is a bonded to atom1
        """
        pass

    @abstractmethod
    def get_neighbors(self):
        """
        Returns
        -------
        neighbors : list[ChemPer Atoms]
            Atoms that are one bond away from this atom
        """
        pass

    @abstractmethod
    def get_bonds(self):
        """
        Returns
        -------
        bonds : list[ChemPer Bonds]
            Bonds connected to this atom
        """
        pass

    @abstractmethod
    def get_molecule(self):
        """
        Extracts the parent molecule this atom is from.

        Returns
        -------
        mol : ChemPer Mol
            Molecule this atom is stored in
        """
        pass


# =======================================
# Bond Class
# =======================================

class BondAdapter(ABC):
    """
    This is a ChemPer wrapper for a bond from
    one of the cheminformatics toolkits.
    ChemPer Bonds are initiated from the reference package object.
    Currently we support OpenEye and RDKit toolkits.

    Attributes
    ----------
    bond : Bond from reference class
    """
    @abstractmethod
    def get_order(self):
        """
        Returns
        -------
        order : int or float
            This is the absolute order, returns 1.5 if bond is aromatic
        """
        pass

    @abstractmethod
    def get_atoms(self):
        """
        Returns
        -------
        atoms : list[ChemPer Atoms]
            The two atoms connected by this bond
        """
        pass

    @abstractmethod
    def is_ring(self):
        """
        Returns
        -------
        is_ring : boolean
            True if bond is a part of a ring, otherwise False
        """
        pass

    @abstractmethod
    def is_aromatic(self):
        """
        Returns
        -------
        is_aromatic : boolean
            True if it is an aromatic bond
        """
        pass

    @abstractmethod
    def is_single(self):
        """
        Returns
        -------
        is_single : boolean
            True if it is a single bond
        """
        pass

    @abstractmethod
    def is_double(self):
        """
        Returns
        -------
        is_double : boolean
            True if it is a double bond
        """
        pass

    @abstractmethod
    def is_triple(self):
        """
        Returns
        -------
        is_triple : boolean
            True if it is a triple bond
        """
        pass

    @abstractmethod
    def get_molecule(self):
        """
        Extracts the parent molecule this bond is from

        Returns
        -------
        mol : ChemPer Mol
            Molecule this bond is stored in
        """
        pass

    @abstractmethod
    def get_index(self):
        """
        Returns
        -------
        index : int
            index of this bond in its parent molecule
        """
        pass
