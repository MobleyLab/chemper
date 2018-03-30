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

class Mol(MolAdapter):
    def __init__(self, mol):
        if type(mol) != oechem.OEMol:
            raise Exception("Expecting an OEMol object instead of %s" % type(mol))
        self.mol = mol

    def get_atoms(self):
        return [Atom(a) for a in self.mol.GetAtoms()]

    def get_atom_by_index(self, idx):
        return Atom(self.mol.GetAtom(oechem.OEHasAtomIdx(idx)))

    def get_bonds(self):
        return [Bond(b) for b in self.mol.GetBonds()]

    def get_bond_by_index(self, idx):
        return Bond(self.mol.GetBond(oechem.OEHasBondIdx(idx)))

    def smirks_search(self, smirks):
        matches = list()

        ss = oechem.OESubSearch()
        if not ss.Init(smirks):
            # TODO: write custom exceptions?
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
        smiles = oechem.OEMolToSmiles(self.mol)
        return smiles

class MolFromSmiles(Mol):
    def __init__(self, smiles):
        mol = oechem.OEMol()
        if not oechem.OESmilesToMol(mol, smiles):
            raise ValueError('Could not parse SMILES %s' % smiles)
        oechem.OEAddExplicitHydrogens(mol)
        Mol.__init__(self, mol)

# =======================================
# Atom Class
# =======================================


class Atom(AtomAdapter):
    def __init__(self, atom):
        if type(atom) != oechem.OEAtomBase:
            raise Exception("Expecting an OEAtomBase object instead of %s" % type(atom))
        self.atom = atom
        self._index = self.atom.GetIdx()

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
        return self.atom.IsAromatic()

    def get_index(self):
        return self._index

    def is_connected_to(self, atom2):
        return self.atom.IsConnected(atom2.atom)

    def get_neighbors(self):
        return [Atom(a) for a in self.atom.GetAtoms()]

    def get_bonds(self):
        return [Bond(b) for b in self.atom.GetBonds()]

    def get_molecule(self):
        mol = oechem.OEMol(self.atom.GetParent())
        self.atom = mol.GetAtom(oechem.OEHasAtomIdx(self._index))
        return Mol(mol)

# =======================================
# Bond Class
# =======================================


class Bond(BondAdapter):
    def __init__(self, bond):
        if type(bond) != oechem.OEBondBase:
            raise Exception("Expecting an OEBondBase object instead of %s" % type(bond))
        self.bond = bond
        self._order = self.bond.GetOrder()
        if self.is_aromatic():
            self._order = 1.5

        self._idx = self.bond.GetIdx()

    def get_order(self):
        return self._order

    def get_atoms(self):
        beginning = Atom(self.bond.GetBgn())
        end = Atom(self.bond.GetEnd())
        return [beginning, end]

    def is_ring(self):
        return self.bond.IsInRing()

    def is_aromatic(self):
        return self.bond.IsAromatic()

    def is_single(self):
        return self._order == 1

    def is_double(self):
        return self._order == 2

    def is_triple(self):
        return self._order == 3

    def get_molecule(self):
        mol = oechem.OEMol(self.bond.GetParent())
        self.bond = mol.GetBond(oechem.OEHasBondIdx(self._idx))
        return Mol(mol)

    def get_index(self):
        return self._idx