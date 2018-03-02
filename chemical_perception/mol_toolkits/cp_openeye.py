"""
Cheminformatics tools using OpenEye Toolkits

The classes provided here follow the structure in adapters.
This is a wrapper allowing our actual package to use openeye toolkits
Authors: Caitlin C. Bannan
"""

from chemical_perception.mol_toolkits.adapters import MolAdapter, AtomAdapter
from openeye import oechem


# =======================================
# Molecule Class
# =======================================

class MolOE(MolAdapter):
    def __init__(self, mol):
        self.mol = mol
        return

    def get_atoms(self):
        return [AtomOE(a) for a in self.mol.GetAtoms()]

    def get_atom_by_index(self, idx):
        return AtomOE(mol.GetAtom(oechem.OEHasAtomIdx(idx)))

    def smirks_search(self, smirks):
        matches = list()

        ss = oechem.OESubSearch()
        if not ss.Init(smirks):
            #TODO: write custom exceptions?
            raise Exception("Error parsing SMIRKS %s" % smirks)

        for match in ss.Match(self.mol, False):
            d = dict()
            for ma in match.GetAtoms():
                atom_idx = ma.pattern.GetMapIdx()
                if atom_idx !=0:
                    d[atom_idx] = AtomOE(ma.target)
            matches.append(d)

        return matches


# =======================================
# Atom Class
# =======================================


def AtomOE(AtomAdapter):
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
