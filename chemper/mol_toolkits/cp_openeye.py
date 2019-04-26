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


# Note - doc strings on these functions are inherited from
#        there Adapters. To read these strings see adapters.py.

# =======================================
# Molecule Class
# =======================================

class Mol(MolAdapter):
    def __init__(self, mol):
        """
        ChemPer created from an OEMol

        Parameters
        ----------
        mol : openeye OEMol object
            openeye molecule to convert to ChemPer Mol object
        """
        if not isinstance(mol, oechem.OEMolBase):
            raise TypeError("Expecting an OEMol object instead of %s" % type(mol))
        self.mol = mol

    def __str__(self): return self.get_smiles()

    @classmethod
    def from_smiles(cls, smiles):
        mol = oechem.OEMol()
        if not oechem.OESmilesToMol(mol, smiles):
            raise ValueError('Could not parse SMILES %s' % smiles)
        oechem.OEAddExplicitHydrogens(mol)
        return cls(mol)

    def set_aromaticity_mdl(self):
        oechem.OEClearAromaticFlags(self.mol)
        oechem.OEAssignAromaticFlags(self.mol, oechem.OEAroModel_MDL)
        oechem.OEAssignHybridization(self.mol)

    def get_atoms(self):
        return [Atom(a) for a in self.mol.GetAtoms()]

    def get_atom_by_index(self, idx):
        return Atom(self.mol.GetAtom(oechem.OEHasAtomIdx(idx)))

    def get_bonds(self):
        return [Bond(b) for b in self.mol.GetBonds()]

    def get_bond_by_index(self, idx):
        return Bond(self.mol.GetBond(oechem.OEHasBondIdx(idx)))

    def get_bond_by_atoms(self, atom1, atom2):
        if not atom1.is_connected_to(atom2):
            return None
        return Bond(self.mol.GetBond(atom1.atom, atom2.atom))

    def smirks_search(self, smirks):
        cmol = oechem.OEMol(self.mol)

        matches = list()

        ss = oechem.OESubSearch()
        if not ss.Init(smirks):
            raise ValueError("Error parsing SMIRKS %s" % smirks)

        # set maximum matches in substructure search to infinite (0 in API)
        ss.SetMaxMatches(0)
        for match in ss.Match(cmol, False):
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

# =======================================
# Atom Class
# =======================================

class Atom(AtomAdapter):
    def __init__(self, atom):
        """
        ChemPer Atom created from an OEAtom

        Parameters
        ----------
        atom: OEAtomBase
            Atom object from an OpenEye molecule
        """
        if not isinstance(atom, oechem.OEAtomBase):
            raise TypeError("Expecting an OEAtomBase object instead of %s" % type(atom))
        self.atom = atom
        self._idx = self.atom.GetIdx()

    def __str__(self): return "%i%s" % (self._idx,
                                        oechem.OEGetAtomicSymbol(self.atomic_number()))

    def atomic_number(self): return self.atom.GetAtomicNum()

    def degree(self): return self.atom.GetDegree()

    def connectivity(self):
        return len([b for b in self.atom.GetBonds()])

    def valence(self): return self.atom.GetValence()

    def formal_charge(self): return self.atom.GetFormalCharge()

    def hydrogen_count(self): return self.atom.GetTotalHCount()

    def ring_connectivity(self):
        return len([b for b in self.atom.GetBonds(oechem.OEBondIsInRing())])

    def min_ring_size(self):
        return oechem.OEAtomGetSmallestRingSize(self.atom)

    def is_aromatic(self): return self.atom.IsAromatic()

    def get_index(self): return self._idx

    def is_connected_to(self, atom2):
        if not isinstance(atom2.atom, oechem.OEAtomBase):
            return False
        return self.atom.IsConnected(atom2.atom)

    def get_neighbors(self):
        return [Atom(a) for a in self.atom.GetAtoms()]

    def get_bonds(self):
        return [Bond(b) for b in self.atom.GetBonds()]

    def get_molecule(self):
        mol = oechem.OEMol(self.atom.GetParent())
        self.atom = mol.GetAtom(oechem.OEHasAtomIdx(self._idx))
        return Mol(mol)

# =======================================
# Bond Class
# =======================================


class Bond(BondAdapter):
    def __init__(self, bond):
        """
        ChemPer Bond created from an OEBond

        Parameters
        ----------
        bond: OEBondBase
            Bond object from an OpenEye molecule
        """
        if not isinstance(bond, oechem.OEBondBase):
            raise TypeError("Expecting an OEBondBase object instead of %s" % type(bond))
        self.bond = bond

        # save index
        self._idx = self.bond.GetIdx()

        # store order information
        self._order = self.bond.GetOrder()
        if self.is_aromatic():
            self._order = 1.5

        orders = {1:'-', 2:'=', 3:'#', 1.5:':'}
        self._order_symbol = orders.get(self._order, '~')

        # save atoms in bond
        self._beginning = Atom(self.bond.GetBgn())
        self._end = Atom(self.bond.GetEnd())

    def __str__(self):
        return "%i %s%s%s" % (self.get_index(), self._beginning,
                              self._order_symbol, self._end)

    def get_order(self): return self._order

    def get_atoms(self): return [self._beginning, self._end]

    def is_ring(self): return self.bond.IsInRing()

    def is_aromatic(self): return self.bond.IsAromatic()

    def is_single(self): return self._order == 1

    def is_double(self): return self._order == 2

    def is_triple(self): return self._order == 3

    def get_molecule(self):
        mol = oechem.OEMol(self.bond.GetParent())
        self.bond = mol.GetBond(oechem.OEHasBondIdx(self._idx))
        return Mol(mol)

    def get_index(self): return self._idx

# =====================================================================
# functions for importing molecules from files
# =====================================================================

def mols_from_mol2(mol2_file):
    return mols_from_file(mol2_file)

def mols_from_file(mol_file):
    """
    Parses a standard molecule file into ChemPer molecules using OpenEye toolkits

    Parameters
    ----------
    mol_file: str
              relative or full path to molecule containing the molecule file
              that is accessible from the current working directory

    Returns
    -------
    mols: list[ChemPer Mol]
          list of molecules in the mol2 file as ChemPer Mols
    """
    import os
    if not os.path.exists(mol_file):
        from chemper.chemper_utils import get_data_path
        mol_path = get_data_path(os.path.join('molecules', mol_file))

        if not os.path.exists(mol_path):
            raise IOError("File '%s' not found locally or in chemper/data/molecules." % mol_file)
        else:
            mol_file = mol_path

    molecules = list()

    # make Openeye input file stream
    ifs = oechem.oemolistream(mol_file)

    oemol = oechem.OECreateOEGraphMol()
    while oechem.OEReadMolecule(ifs, oemol):
        # if an SD file, the molecule name may be in the SD tags
        if oemol.GetTitle() == '':
            name = oechem.OEGetSDData(oemol, 'name').strip()
            oemol.SetTitle(name)
        # Append to list.
        molecules.append(Mol(oechem.OEMol(oemol)))
    ifs.close()

    return molecules

