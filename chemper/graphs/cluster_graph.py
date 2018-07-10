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

        def __str__(self):
            """
            Overrides the default implementation
            returns a string of this atom storage as a SMIRKS
            """
            return self.as_smirks()

        def __lt__(self, other):
            """
            Overrides the default implementation
            This method was primarily written for making SMIRKS patterns predictable.
            If atoms are sortable, then the SMIRKS patterns are always the same making
            tests easier to write. However, the specific sorting was created to also make SMIRKS
            output as human readable as possible, that is to at least make it easier for a
            human to see how the indexed atoms are related to each other.
            It is typically easier for humans to read SMILES/SMARTS/SMIRKS with less branching (indicated with ()).

            For example in:
            [C:1]([H])([H])~[N:2]([C])~[O:3]
            it is easier to see that the atoms C~N~O are connected in a "line" instead of:
            [C:1]([N:2]([O:3])[C])([H])[H]
            which is equivalent, but with all the () it is hard for a human to read the branching

            Parameters
            ----------
            other: AtomStorage

            Returns
            -------
            is_less_than: boolean
                self is less than other
            """
            # if either smirks index is None, then you can't directly compare
            # make a temporary index that is negative if it was None
            self_index = self.smirks_index if self.smirks_index is not None else -1000
            other_index = other.smirks_index if other.smirks_index is not None else -1000
            # if either index is greater than 0, the one that is largest should go at the end of the list
            if self_index > 0 or other_index > 0:
                return self_index < other_index

            # Both SMIRKS indices are not positive or None so compare the SMIRKS patterns instead
            return self.as_smirks() < other.as_smirks()

        def make_atom_decorators(self, atom):
            """
            extract information from a chemper atom that would be useful in a smirks

            parameters
            ----------
            atom: chemper atom object

            returns
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
            # TODO: add better description here
            Parameters
            ----------
            atom: chemper Atom

            Returns
            -------
            score: int
                 how similar is atom to current storage, max of 7 for all decorators identical
                 0 if atom's atomic number not included in current set
            """
            # If decorators is empty (no known atom information, return 7 (current max)
            if len(self.decorators) == 0:
                return 7

            score = 0
            decs = self.make_atom_decorators(atom)
            # extract atomic number for new atom
            test_atomic = int(decs[0][1:])

            for ref in self.decorators:
                # get atomic number for this set of decorators
                ref_atomic = int(ref[0][1:])
                current = len(set(ref) & set(decs))

                # if atomic numbers don't agree, get the number of common decorators / 10
                # if there are no matching atomic numbers, priority should still be given
                # when the current atom matches stored decorators most closely
                if ref[0] != decs[0]:
                    current = current / 10.0

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
                    self.order.add(bond.get_order())
                    self.ring.add(bond.is_ring())

            self.smirks_index = smirks_index

        def __str__(self):
            """
            Overrides the default implementation
            returns a string of this bond storage as a SMIRKS
            """
            return self.as_smirks()

        def __lt__(self, other):
            """
            Overrides the default implementation
            Used for sorting, while I don't know why you would want to sort Bond Storage
            it seemed like a good idea to add this one when I added the sorting for Atom Storage
            Parameters
            ----------
            other: BondStorage object

            Returns
            -------
            is_less_than: boolean
                self is less than other
            """
            return self.as_smirks() < other.as_smirks()

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
                order = ','.join([self.order_dict.get(o, '~') for o in sorted(list(self.order))])

            # the ring set has booleans, if the length of the set is 1 then only ring (@) or non-ring (!@)
            # bonds haven been added to this storage and we AND that decorator to the end of the bond
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
            self.order.add(bond.get_order())
            self.ring.add(bond.is_ring())

        def compare_bond(self, bond):
            """

            Parameters
            ----------
            bond: chemper Bond
                bond you want to compare to the current storage

            Returns
            -------
            score: int (0,1,2)
                1 for if the bond order is in storage plus
                1 base on if this is a ring bond
            """
            score = 0
            if bond.get_order() in self.order:
                score += 1

            # the ring set has booleans, if the length of the set is 1 then only ring or non-ring
            # bonds haven been added to this storage. That is the only time the ring contributes to the score
            if len(self.ring) == 1 and self.ring[0] == bond.is_ring():
                score += 1

            return score

    # Initiate ClusterGraph
    def __init__(self, mols=None, smirks_atoms_lists=None, layers=0):
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
        # if layers is 0 there are no more atoms to add so end the recursion
        if layers == 0:
            return

        # find atom neighbors that are not already included in SMIRKS indexed atoms
        atom_neighbors = [a for a in atom.get_neighbors() if a.get_index() not in idx_dict]

        # get the smirks indices already added to the storage
        # This includes all previous layers since the idx_dict is updated as you go
        smirks = [e for k,e in idx_dict.items()]

        # similar to atoms find neighbors already in the graph that haven't already been used
        storage_neighbors = [s for s in self.get_neighbors(storage) if s.smirks_index not in smirks]

        # If the storage neighbors isn't empty
        if storage_neighbors:
            # min_smirks in this case will be used for any new storages created
            # for example, if the current graph has 2 neighbors, but the new atom has 3
            min_smirks = min([s.smirks_index for s in storage_neighbors]) - 1
            # pair up the atoms with their storage
            pairs = self.find_pairs(atom_neighbors, storage_neighbors)

        # storage doesn't have any neighbors
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
            if new_layers == 0:
                return

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

        atom_dict = dict()
        storage_dict = dict()

        for idx, a in enumerate(atoms):
            g.add_node(idx+1, bipartite=0)
            atom_dict[idx+1] = a
        for idx, s in enumerate(storages):
            g.add_node((idx*-1)-1, bipartite=1)
            storage_dict[(idx*-1)-1] = s

        for a_idx, a in atom_dict.items():
            for s_idx, s in storage_dict.items():
                # TODO: scores only consider one atom and storage, should it look at neighbors?
                score = s.compare_atom(a)
                g.add_edge(a_idx,s_idx,weight=score+1)

        matching = nx.algorithms.max_weight_matching(g,maxcardinality=False)
        pair_set = set()

        for idx_1, idx_2 in matching:
            if idx_1 in atom_dict:
                a = atom_dict[idx_1]
                s = storage_dict[idx_2]
                pair_set.add(idx_1)
            else:
                a = atom_dict[idx_2]
                s = storage_dict[idx_1]
                pair_set.add(idx_2)
            pairs.append((a,s))

        for a_idx, a in atom_dict.items():
            if a_idx not in pair_set:
                pairs.append((a,None))

        return pairs

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
            atom_dict = dict()
            for key, atom_index in smirks_atoms.items():
                atom_dict[atom_index] = key
                atom1 = mol.get_atom_by_index(atom_index)
                self.atom_by_smirks_index[key].add_atom(atom1)

                for neighbor_key, neighbor_index in smirks_atoms.items():
                    # check for connecting bond
                    atom2 = mol.get_atom_by_index(neighbor_index)
                    bond = mol.get_bond_by_atoms(atom1, atom2)
                    if bond is not None and (neighbor_key, key) in self.bond_by_smirks_index:
                        bond_smirks = (neighbor_key, key)
                        self.bond_by_smirks_index[bond_smirks].add_bond(bond)

            for smirks_index, atom_index in smirks_atoms.items():
                atom = mol.get_atom_by_index(atom_index)
                storage = self.atom_by_smirks_index[smirks_index]
                self._add_layers(mol, atom, storage, self.layers, atom_dict)

