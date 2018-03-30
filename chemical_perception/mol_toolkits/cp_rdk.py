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

class Mol(MolAdapter):
    def __init__(self, mol):
        if type(mol) != Chem.rdchem.Mol:
            raise Exception("Expecting an rdchem.Mol instead of %s" % type(mol))
        self.mol = mol

    def get_atoms(self):
        return [Atom(a) for a in self.mol.GetAtoms()]

    def get_atom_by_index(self, idx):
        return Atom(self.mol.GetAtomWithIdx(idx))

    def get_bonds(self):
        return [Bond(b) for b in self.mol.GetBonds()]

    def get_bond_by_index(self, idx):
        return Bond(self.mol.GetBondWithIdx(idx))

    def smirks_search(self, smirks):
        matches = list()

        ss = Chem.MolFromSmarts(smirks)
        if ss is None:
            # TODO: write custom exceptions?
            raise ValueError("Error parsing SMIRKS %s" % smirks)

        # get atoms in query mol with smirks index
        maps = dict()
        for qatom in ss.GetAtoms():
            smirks_idx = qatom.GetAtomMapNum()
            if smirks_idx != 0:
                maps[smirks_idx] = qatom.GetIdx()

        for match in self.mol.GetSubstructMatches(ss, False):
            d = {k:self.get_atom_by_index(match[e]) for k,e in maps.items()}
            matches.append(d)

        return matches

    def get_smiles(self):
        smiles = Chem.MolToSmiles(Chem.RemoveHs(self.mol))
        return smiles

class MolFromSmiles(Mol):
    def __init__(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError('Could not parse SMILES %s' % smiles)
        Mol.__init__(self, Chem.AddHs(mol))

# =======================================
# Atom Class
# =======================================


class Atom(AtomAdapter):
    def __init__(self, atom):
        if type(atom) != Chem.rdchem.Atom:
            raise Exception("Expecting an rdchem.Atom instead of %s" % type(atom))
        self.atom = atom

    def atomic_number(self):
        return self.atom.GetAtomicNum()

    def degree(self):
        return self.atom.GetDegree()

    def connectivity(self):
        return self.atom.GetTotalDegree()

    def valence(self):
        return self.atom.GetTotalValence()

    def formal_charge(self):
        return self.atom.GetFormalCharge()

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

    def is_aromatic(self):
        return self.atom.GetIsAromatic()

    def get_index(self):
        return self.atom.GetIdx()

    def is_connected_to(self, atom2):
        if not type(atom2.atom) is Chem.rdchem.Atom:
            # TODO: raise exception/return something else?
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
        if type(bond) != Chem.rdchem.Mol:
            raise Exception("Expecting an rdchem.Bond instead of %s" % type(bond))
        self.bond = bond
        self.order = self.bond.GetBondTypeAsDouble()
        self.beginning = Atom(self.bond.GetBeginAtom())
        self.end = Atom(self.bond.GetEndAtom())

    def get_order(self):
        return self.order

    def get_atoms(self):
        return [self.beginning, self.end]

    def is_ring(self):
        return self.bond.IsInRing()

    def is_aromatic(self):
        return self.bond.GetIsAromatic()

    def is_single(self):
        return self.order == 1

    def is_double(self):
        return self.order == 2

    def is_triple(self):
        return self.order == 3

    def get_molecule(self):
        mol = Chem.Mol(self.bond.GetOwningMol())
        return Mol(mol)

    def get_index(self):
        return self.bond.GetIdx()