"""
Cheminformatics tools using RDKit

The classes provided here follow the structure in adapters.
This is a wrapper allowing our actual package to use RDKit
Authors: Caitlin C. Bannan
"""

from chemical_perception.mol_toolkits.adapters import MolAdapter, AtomAdapter, BondAdapter
from rdkit import Chem

# =======================================
# Molecule Class
# =======================================

class MolRDK(MolAdapter):

    def get_atoms(self):
        return

    def get_atom_by_index(self, idx):
        return

    def get_bonds(self):
        return

    def get_bond_by_index(self, idx):
        return

    def smirks_search(self, smirks):
        return

    def get_smiles(self):
        return

# =======================================
# Atom Class
# =======================================


class AtomRDK(AtomAdapter):

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
