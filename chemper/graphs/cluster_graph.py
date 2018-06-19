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
            atoms: chemper Atom or list of chemper Atoms
                this is one or more atoms whose information should be stored
            smirks_index: int
                SMIRKS index (:n) for writing SMIRKS
                if the value is less than zero it is used for storage purposes
                only as SMIRKS can only be written with positive integer indices
            """
            self.decorators = set()
            if atoms is not None:
                # check if this is a single atom
                if 'Atom' in str(type(atoms)):
                    atoms = [atoms]

                # otherwise it should be iterable
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

        def compare_atom(self, atom):
            """

            Parameters
            ----------
            atom: chemper Atom

            Returns
            -------
            score: int
                 how similar is atom to current storage, max of 7 for all decorators identical
                 0 if atom's atomic number not included in current set
            """
            score = 0
            decs = self.make_atom_decorators(atom)
            atomic_nums = [d[0] for d in self.decorators]

            for ref in self.decorators:
                # TODO: determine how to handle different atomic numbers?
                # right now they are treated as 0 agreement, but you could consider ring decorators as importantly similar
                if ref[0] != decs[0]:
                    continue

                # if the atomic numbers are the same check how many decorators they have in common
                current = len(set(ref) & set(decs))
                if current > score:
                    score = current

            return score

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
                if 'Bond' in str(type(bonds)):
                    bonds = [bonds]
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
        #if layers != 0:
        #    raise Exception("Currently we only support 0 layers for ClusterGraphs")

        self.mols = list()
        self.smirks_atoms_lists = list()
        self.layers = layers

        if mols is not None:
            if len(mols) != len(smirks_atoms_lists):
                raise Exception('Number of molecules and smirks dictionaries should be equal')

            for idx, mol in enumerate(mols):
                self.add_mol(mol, smirks_atoms_lists[idx])

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
        atom_dict = dict()
        sorted_keys = sorted(list(smirks_atoms.keys()))
        for key in sorted_keys:
            atom_index = smirks_atoms[key]
            atom_dict[atom_index] = key

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

        for smirks_index, atom_index in smirks_atoms.items():
            atom = mol.get_atom_by_index(atom_index)
            storage = self.atom_by_smirks_index[smirks_index]
            self._add_layers(mol, atom, storage, self.layers, atom_dict)

    def _add_layers(self, mol, atom, storage, layers, idx_dict):
        """

        Parameters
        ----------
        mol: chemper Mol
            molecule containing provided atom
        atom: chemper Atom
        storage: AtomStorage
            corresponding to the chemper Atom provided
        layers: int or 'all'
            number of layers left to add (or all)
        idx_dict: dict
            form {atom_index: smirks_index} for this molecule
        """
        if layers == 0:
            return

        atom_neighbors = [a for a in atom.get_neighbors() if a.get_index() not in idx_dict]

        smirks = [e for k,e in idx_dict.items()]
        storage_neighbors = [s for s in self.get_neighbors(storage) if s.smirks_index not in smirks]

        if storage_neighbors:
            min_smirks = min([s.smirks_index for s in storage_neighbors]) - 1
            pairs = self.find_pairs(atom_neighbors, storage_neighbors)

        else:
            min_smirks = storage.smirks_index * 10
            if min_smirks > 0:
                min_smirks = min_smirks * -1
            pairs = [(a, None) for a in atom_neighbors]

        new_pairs = list()

        for new_atom, new_storage in pairs:
            new_bond = mol.get_bond_by_atoms(atom, new_atom)
            if new_storage is None:
                new_bond_smirks = (storage.smirks_index, min_smirks)

                adding_new_storage = self.add_atom(new_atom,new_bond,storage,
                                            min_smirks, new_bond_smirks)

                idx_dict[new_atom.get_index()] = min_smirks
                self.atom_by_smirks_index[min_smirks] = adding_new_storage
                min_smirks -= 1
                new_pairs.append((new_atom, adding_new_storage))

            else:
                new_storage.add_atom(new_atom)

                bond_smirks = (storage.smirks_index, new_storage.smirks_index)
                self.bond_by_smirks_index[bond_smirks].add_bond(new_bond)
                new_pairs.append((new_atom, new_storage))
                idx_dict[new_atom.get_index()] = new_storage.smirks_index

        if layers == 'all':
            new_layers = 'all'
        else:
            new_layers = layers - 1

        for new_atom, new_storage in new_pairs:
            self._add_layers(mol, new_atom, new_storage, new_layers, idx_dict)

    def find_pairs(self, atoms, storages):
        """
        Determines the best pairing of atoms to an existing set of AtomStorages

        Parameters
        ----------
        atoms: list of chemper Atoms
        storages: list of AtomStorage

        Returns
        -------
        pairs: list of tuples
            pairs of atoms and storage objects that are most similar
        """
        pairs = list()

        g = nx.Graph()
        for a in atoms:
            g.add_node(a, bipartite=0)
        for s in storages:
            g.add_node(s, bipartite=1)

        for a in atoms:
            for s in storages:
                # TODO: scores only consider one atom and storage, should it look at neighbors?
                score = s.compare_atom(a)
                g.add_edge(a,s,weight=score)

        matching = nx.algorithms.max_weight_matching(g,maxcardinality=False)

        for a in atoms:
            if a in matching:
                s = matching[a]
                pairs.append((a,s))
            else:
                # All atoms need to be in final pairs, if not in matching then pair with None
                pairs.append((a,None))

        return pairs

    def _add_first_layers(self, mol, atom, storage, layers, idx_dict):
        """

        Parameters
        ----------
        mol: chemper Mol
            molecule atom is coming from (allows easier access to neighbors)
        atom: chemper Atom
        storage: AtomStorage
            corresponding to the chemper Atom
        layers: int or 'all'
            number of layers left to add (or all)
        idx_dict: dict
            form {atom_index: smirks_index in current graph}
        """
        if layers == 0:
            return

        top_smirks_idx = storage.smirks_index * 10
        if top_smirks_idx > 0:
            top_smirks_idx = -1 * top_smirks_idx

        atom_neighbors = [a for a in atom.get_neighbors() if a.get_index() not in idx_dict]

        for idx, new_atom in enumerate(atom_neighbors):
            new_bond = mol.get_bond_by_atoms(atom, new_atom)
            new_smirks_index = top_smirks_idx - idx - 1
            new_bond_smirks = (storage.smirks_index, new_smirks_index)
            new_storage = self.add_atom(new_atom,new_bond,storage,
                                        new_smirks_index,new_bond_smirks)

            idx_dict[new_atom.get_index()] = new_smirks_index
            self.atom_by_smirks_index[new_smirks_index] = new_storage

        go_again = False
        if layers == 'all':
            new_layers = 'all'
            go_again = True
        elif layers > 1:
            new_layers = layers - 1
            go_again = True

        if go_again:
            for new_atom in atom_neighbors:
                new_smirks_index = idx_dict[new_atom.get_index()]
                new_storage = self.atom_by_smirks_index[new_smirks_index]
                self._add_first_layers(mol, new_atom, new_storage, layers, idx_dict)

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
