"""
fragment_graph.py

ChemPerGraph are a class for tracking molecular fragments based on information about the atoms and bonds they contain.
You can combine them or take a difference in order to find distinguishing characteristics in a set of clustered
molecular sub-graphs.

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

import networkx as nx
from chemper.graphs.fragment_graph import ChemPerGraph


class ClusterGraph(ChemPerGraph):
    """
    ChemPerGraphs are a graph based class for storing atom and bond information.
    They use the chemper.mol_toolkits Atoms, Bonds, and Mols
    """
    class AtomStorage(object):
        """
        AtomStorage tracks information about an atom
        """
        def __init__(self, atoms=None, smirks_index=None):
            """
            Parameters
            ----------
            atoms: list of chemper Atom objects
                this is one or more atoms whose information should be stored
            smirks_index: int
                SMIRKS index (:n) for writing SMIRKS
                if the value is less than zero it is used for storage purposes
                only as SMIRKS can only be written with positive integer indices
            """
            self.decorators = set()
            if atoms is not None:
                for atom in atoms:
                    self.decorators.add(self.make_atom_decorators(atom))
            self.smirks_index = smirks_index

        def make_atom_decorators(self, atom):
            """
            Extract information from a chemper Atom that would be useful in a SMIRKS

            Parameters
            ----------
            atom: chemper Atom object

            Returns
            -------
            decorators: tuple of str
                tuple of all possible decorators for this atom
            """
            aromatic = 'a' if atom.is_aromatic() else 'A'
            charge = atom.formal_charge()
            if charge >= 0:
                charge = '+%i' % charge
            else:
                charge = '%i' % charge

            return (
                '#%i' % atom.atomic_number(),
                aromatic,
                'H%i' % atom.hydrogen_count(),
                'X%i' % atom.connectivity(),
                'x%i' % atom.ring_connectivity(),
                'r%i' % atom.min_ring_size(),
                charge
                )

        def as_smirks(self):
            """
            Returns
            -------
            smirks: str
                how this atom would be represented in a SMIRKS string
                with the minimal combination of SMIRKS decorators
            """
            if len(self.decorators) == 0:
                if self.smirks_index is None or self.smirks_index <= 0:
                    return '[*]'
                return '[*:%i]' % self.smirks_index

            base_smirks = ','.join([''.join(l) for l in sorted(list(self.decorators))])
            if self.smirks_index is None or self.smirks_index <= 0:
                return '[%s]' % base_smirks

            return '[%s:%i]' % (base_smirks, self.smirks_index)

        def add_atom(self, atom):
            """
            Expand current AtomStorage by adding information about
            a new chemper Atom

            Parameters
            ----------
            atom: chemper Atom
            """
            self.decorators.add(self.make_atom_decorators(atom))

    class BondStorage(object):
        """
        BondStorage tracks information about a bond
        """
        def __init__(self, bonds=None, smirks_index=None):
            """
            Parameters
            ----------
            bonds: list of chemper Bond objects
                this is one or more bonds whose information should be stored
            smirks_index: int or float
                bonds do not have smirks indices so this is only used for internal storage
            """
            self.order = set()
            self.ring = set()
            self.order_dict = {1:'-', 1.5:':', 2:'=', 3:'#'}
            if bonds is not None:
                for bond in bonds:
                    self.order.add(self.order_dict.get(bond.get_order(), '~'))
                    self.ring.add(bond.is_ring())

            self.smirks_index = smirks_index

        def as_smirks(self):
            """
            Returns
            -------
            smirks: str
                how this bond would be represented in a SMIRKS string
                using only the required number of
            """
            if len(self.order) == 0:
                order = '~'
            else:
                order = ','.join(sorted(list(self.order)))

            if len(self.ring) == 1:
                if list(self.ring)[0]:
                    return order+';@'
                else:
                    return order+';!@'

            return order

        def add_bond(self, bond):
            """
            Expand current BondStorage by adding information about
            a new chemper Bond

            Parameters
            ----------
            bond: chemper Bond
            """
            self.order.add(self.order_dict.get(bond.get_order()))
            self.ring.add(bond.is_ring())

    # Initiate ClusterGraph
    def __init__(self, mols=None, smirks_atoms_lists=None, layers=0):
        # TODO: figure out layers, currently only supporting layers=0
        """
        Initialize a ChemPerGraph from a molecule and list of indexed atoms

        Parameters
        ----------
        mols: list of chemper Mols
            (optional) molecules to initiate ClusterGraph
        smirks_atoms_lists: list of list of dict
            (optional) atom indices by smirks index for each molecule
            required if a list of molecules is provided.
            This is a list of dictionaries of the form [{smirks_index: atom_index}]
            for each molecule provided
        layers: int
            (optional) currently only 0 is supported
            how many atoms out from the smirks indexed atoms do you wish save (default=0)
            'all' will lead to all atoms in the molecule being specified (not recommended)
            # TODO: figure out the best form for the smirks_atoms_lists
        """
        ChemPerGraph.__init__(self)

        # TODO: make sure to remove this when layers works
        if layers != 0:
            raise Exception("Currently we only support 0 layers for ClusterGraphs")

        self.mols = list()
        self.smirks_atoms_lists = list()

        if mols is not None:
            if len(mols) != len(smirks_atoms_lists):
                raise Exception('Number of molecules and smirks dictionaries should be equal')

            for idx, mol in enumerate(mols):
                self.add_mol(mol, smirks_atoms_lists[idx])

        # TODO: figure out how to extend beyond just "SMIRKS indexed" atoms
        # self.atom_by_index = dict()
        #for smirks_key, atom_storage in self.atom_by_smirks_index.items():
        #    self._add_layers(atom_storage, layers)

    def add_mol(self, mol, smirks_atoms_list):
        """
        Expand the information in this graph by adding a new molecule

        Parameters
        ----------
        mol: chemper Mol object
        smirks_atoms_list: list of dicts
            This is a list of dictionaries of the form [{smirks_index: atom_index}]
            each atom (by index) in the dictionary will be added the relevant
            AtomStorage by smirks index
        """
        if len(self.mols) == 0:
            self._add_first_smirks_atoms(mol, smirks_atoms_list[0])
            self._add_mol(mol, smirks_atoms_list[1:])
        else:
            self._add_mol(mol, smirks_atoms_list)

        self.mols.append(mol)
        self.smirks_atoms_lists.append(smirks_atoms_list)

    def _add_first_smirks_atoms(self, mol, smirks_atoms):
        """
        private function for adding the first molecule to an empty ClusterGraph
        add_mol calls this if the graph is empty

        Parameters
        ----------
        mol: chemper Mol
        smirks_atoms: dict
            dictionary for first atoms to add to the graph in the form {smirks_index: atom_index}]
        """
        sorted_keys = sorted(list(smirks_atoms.keys()))
        for key in sorted_keys:
            atom_index = smirks_atoms[key]

            atom1 = mol.get_atom_by_index(atom_index)
            new_atom_storage = self.AtomStorage([atom1], key)
            self._graph.add_node(new_atom_storage)
            self.atom_by_smirks_index[key] = new_atom_storage
            # self.atom_by_index[atom_index] = new_atom_storage

            # Check for bonded atoms already in the graph
            for neighbor_key in reversed(sorted_keys):
                if neighbor_key not in self.atom_by_smirks_index:
                    continue

                # check if atoms are already connected on the graph
                neighbor_storage = self.atom_by_smirks_index[neighbor_key]
                if nx.has_path(self._graph, new_atom_storage, neighbor_storage):
                    continue

                # check if atoms are connected in the molecule
                atom2 = mol.get_atom_by_index(smirks_atoms[neighbor_key])
                bond = mol.get_bond_by_atoms(atom1, atom2)

                if bond is not None: # Atoms are connected add edge
                    bond_smirks = (neighbor_key, key)
                    bond_storage = self.BondStorage([bond], bond_smirks)
                    self.bond_by_smirks_index[bond_smirks] = bond_storage
                    self._graph.add_edge(new_atom_storage,
                                         neighbor_storage,
                                         bond=bond_storage)

    def _add_mol(self, mol, smirks_atoms_list):
        """
        private function for adding a new molecule
        This is used by add_mol if the graph is not empty, allowing the user to
        not have to track if the graph already has information before adding molecules

        Parameters
        ----------
        mol: chemper Mol
        smirks_atoms_list: list of dicts
            This is a list of dictionaries of the form [{smirks_index: atom_index}]
            each atom (by index) in the dictionary will be added the relevant
            AtomStorage by smirks index
        """
        for smirks_atoms in smirks_atoms_list:
            for key, atom_index in smirks_atoms.items():
                atom1 = mol.get_atom_by_index(atom_index)
                self.atom_by_smirks_index[key].add_atom(atom1)

                for neighbor_key, neighbor_index in smirks_atoms.items():
                    # check for connecting bond
                    atom2 = mol.get_atom_by_index(neighbor_index)
                    bond = mol.get_bond_by_atoms(atom1, atom2)
                    if bond is not None and (neighbor_key, key) in self.bond_by_smirks_index:
                        bond_smirks = (neighbor_key, key)
                        self.bond_by_smirks_index[bond_smirks].add_bond(bond)
