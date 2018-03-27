"""
adapters

This script contains adapters or the structure for atoms, molecules and
substructure searches.
Our chemical perception code is designed to be independent of the users
cheminformatics packages. For each cheminformatics package we support we
will provide classes following the structure in these adapters.
Authors: Caitlin C. Bannan
"""

from abc import ABC, abstractmethod

# =======================================
# Molecule Class
# =======================================

class MolAdapter(ABC):

    @abstractmethod
    def get_atoms(self):
        return

    @abstractmethod
    def get_atom_by_index(self, idx):
        return

    @abstractmethod
    def get_bonds(self):
        return

    @abstractmethod
    def get_bond_by_index(self, idx):
        return

    @abstractmethod
    def smirks_search(self, smirks):
        return

    @abstractmethod
    def get_smiles(self):
        return

# =======================================
# Atom Class
# =======================================


class AtomAdapter(ABC):

    @abstractmethod
    def atomic_number(self):
        return

    @abstractmethod
    def degree(self):
        return

    @abstractmethod
    def connectivity(self):
        return

    @abstractmethod
    def valence(self):
        return

    @abstractmethod
    def formal_charge(self):
        return

    @abstractmethod
    def hydrogen_count(self):
        return

    @abstractmethod
    def min_ring_size(self):
        return

    @abstractmethod
    def ring_connectivity(self):
        return

    @abstractmethod
    def is_aromatic(self):
        return

    @abstractmethod
    def get_index(self):
        return

    @abstractmethod
    def is_connected_to(self, atom2):
        return

    @abstractmethod
    def get_neighbors(self):
        return

    @abstractmethod
    def get_bonds(self):
        return

    @abstractmethod
    def get_molecule(self):
        return


# =======================================
# Bond Class
# =======================================

class BondAdapter(ABC):

    @abstractmethod
    def get_order(self):
        return

    @abstractmethod
    def get_atoms(self):
        return

    @abstractmethod
    def is_ring(self):
        return

    @abstractmethod
    def is_aromatic(self):
        return

    @abstractmethod
    def is_single(self):
        return

    @abstractmethod
    def is_double(self):
        return

    @abstractmethod
    def is_triple(self):
        return

    @abstractmethod
    def get_molecule(self):
        return

    @abstractmethod
    def get_index(self):
        return
