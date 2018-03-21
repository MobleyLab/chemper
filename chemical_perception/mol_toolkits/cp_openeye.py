"""
Cheminformatics tools using OpenEye Toolkits

The classes provided here follow the structure in adapters.
This is a wrapper allowing our actual package to use openeye toolkits
Authors: Caitlin C. Bannan
"""

from chemical_perception.mol_toolkits.adapters import MolAdapter, AtomAdapter, BondAdapter
from openeye import oechem


# =======================================
# Molecule Class
# =======================================

class MolOE(MolAdapter):
    def __init__(self, mol):
        # TODO: add checks that mol is an OEMol?
        self.mol = mol

    def get_atoms(self):
        return [AtomOE(a) for a in self.mol.GetAtoms()]

    def get_atom_by_index(self, idx):
        return AtomOE(self.mol.GetAtom(oechem.OEHasAtomIdx(idx)))

    def get_bonds(self):
        return [BondOE(b) for b in self.mol.GetBonds()]

    def get_bond_by_index(self, idx):
        return BondOE(self.mol.GetBond(oechem.OEHasBondIdx(idx)))

    def smirks_search(self, smirks):
        matches = list()

        ss = oechem.OESubSearch()
        if not ss.Init(smirks):
            # TODO: write custom exceptions?
            raise Exception("Error parsing SMIRKS %s" % smirks)

        for match in ss.Match(self.mol, False):
            d = dict()
            for ma in match.GetAtoms():
                smirks_idx = ma.pattern.GetMapIdx()
                if smirks_idx !=0:
                    d[smirks_idx] = self.mol.get_atom_by_index(
                        ma.target.GetIdx)

            matches.append(d)

        return matches

    def get_smiles(self):
        return oechem.OEMolToSmiles(self.mol)

# =======================================
# Atom Class
# =======================================


class AtomOE(AtomAdapter):
    def __init__(self, atom):
        self.atom = atom

    def atomic_number(self):
        return self.atom.GetAtomicNum()

    def degree(self):
        return self.atom.GetDegree()

    def connectivity(self):
        return len([b for b in self.atom.GetBonds()])

    def valence(self):
        return self.atom.GetValence()

    def formal_charge(self):
        return self.atom.GetFormalCharge()

    def hydrogen_count(self):
        return self.atom.GetTotalHCount()

    def ring_connectivity(self):
        return len([b for b in self.atom.GetBonds(oechem.OEBondIsInRing())])

    def min_ring_size(self):
        return oechem.OEAtomGetSmallestRingSize(self.atom)

    def is_aromatic(self):
        return a.IsAromatic()

    def get_index(self):
        return self.atom.GetIdx()

    def is_connected_to(self, atom2):
        return self.atom.IsConnected(atom2.atom)

    def get_neighbors(self):
        return [AtomOE(a) for a in self.atom.GetAtoms()]

    def get_bonds(self):
        return [BondOE(b) for b in self.atom.GetBonds()]

    def get_molecule(self):
        return MolOE(self.atom.GetParent())

# =======================================
# Bond Class
# =======================================


class BondOE(BondAdapter):
    def __init__(self, bond):
        # TODO: check bond is an OEBond object?
        self.bond = bond
        self.order = self.bond.GetOrder()

    def get_order(self):
        return self.order

    def get_atoms(self):
        return [AtomOE(a) for a in self.bond.GetAtoms()]

    def is_ring(self):
        return self.bond.IsInRing()

    def is_aromatic(self):
        return self.bond.IsAromatic()

    def is_single(self):
        return self.order == 1

    def is_double(self):
        return self.order == 2

    def is_triple(self):
        return self.order == 3

    def get_molecule(self):
        return MolOE(self.bond.GetParent())