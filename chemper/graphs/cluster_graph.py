"""
cluster_graph.py

ClusterGraph are a class for tracking all possible smirks decorators in a group (or cluster)
of molecular fragments. Moving forward these will be used to find the minimum number of
smirks decorators that are required to have a set of smirks patterns that maintain
a given clustering of fragments.

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
        def __init__(self, atoms=None, label=None):
            """
            Parameters
            ----------
            atoms: chemper Atom or list of chemper Atoms
                this is one or more atoms whose information should be stored
            label: int
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
            self.label = label

        def __str__(self): return self.as_smirks()

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
            self_index = self.label if self.label is not None else -1000
            other_index = other.label if other.label is not None else -1000
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

        def as_smirks(self, compress=False):
            """
            Parameters
            ----------
            compress: boolean
                should decorators common to all sets be combined

            Returns
            -------
            smirks: str
                how this atom would be represented in a SMIRKS string
                with the minimal combination of SMIRKS decorators
            """
            if len(self.decorators) == 0:
                if self.label is None or self.label <= 0:
                    return '[*]'
                return '[*:%i]' % self.label

            if compress:
                base_smirks = self._compress_smirks()
            else:
                base_smirks = ','.join([''.join(l) for l in sorted(list(self.decorators))])

            if self.label is None or self.label <= 0:
                return '[%s]' % base_smirks

            return '[%s:%i]' % (base_smirks, self.label)

        def _compress_smirks(self):
            """
            Returns
            -------
            """
            set_decs = [set(d) for d in self.decorators]
            ands = set_decs[0]

            for d_set in set_decs:
                ands = ands & d_set

            ors = [''.join(list(d_set.difference(ands))) for d in set_decs]

            # add commas between ors
            base = ','.join(sorted(ors))
            # add and decorators
            if len(ands) > 0:
                base += ';'+''.join(ands)
            return base


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
            score: float
                A score describing how similar the input atom is to any set of decorators currently
                in this storage, based on its SMIRKS decorators.
                This score ranges from 0 to 7. 7 comes from the number of decorators
                on any atom, if this atom matches perfectly with one of the current decorator sets
                then 7 decorators agree.
                However, if the atomic number doesn't agree, then that set of decorators is considered
                less ideal, thus if the atomic numbers don't agree, then the score is given by
                the number other decorators divided by 10.
                If the current storage is empty, then the score is given as 7 any atom matches
                a wildcard atom.
            """
            # If decorators is empty (no known atom information, return 7 (current max)
            if len(self.decorators) == 0:
                return 7

            score = 0
            decs = self.make_atom_decorators(atom)

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
        def __init__(self, bonds=None, label=None):
            """
            Parameters
            ----------
            bonds: list of chemper Bond objects
                this is one or more bonds whose information should be stored
            label: a label for the object, it can be anything
                unlike atoms, bonds in smirks don't have labels
                so this is only used for labeling the object if wanted
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

            self.label = label

        def __str__(self): return self.as_smirks()

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
                A score describing how similar the input bond is to any set of decorators currently
                in this storage, based on its SMIRKS decorators.

                1 for the bond order +
                1 base on if this is a ring bond
            """
            score = 0
            if bond.get_order() in self.order or len(self.order) == 0:
                score += 1

            # the ring set has booleans, if the length of the set is 1 then only ring or non-ring
            # bonds haven been added to this storage. That is the only time the ring contributes to the score
            if len(self.ring) == 1 and list(self.ring)[0] == bond.is_ring():
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
            This is a list of dictionaries of the form [{smirks index: atom index}]
            for each molecule provided
        layers: int
            (optional) currently only 0 is supported
            how many atoms out from the smirks indexed atoms do you wish save (default=0)
            'all' will lead to all atoms in the molecule being specified (not recommended)
        """
        ChemPerGraph.__init__(self)

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
            This is a list of dictionaries of the form [{smirks index: atom index}]
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
            dictionary for first atoms to add to the graph in the form {smirks index: atom index}]
        """
        atom_dict = dict()
        sorted_keys = sorted(list(smirks_atoms.keys()))
        for key in sorted_keys:
            atom_index = smirks_atoms[key]
            atom_dict[atom_index] = key

            atom1 = mol.get_atom_by_index(atom_index)
            new_atom_storage = self.AtomStorage([atom1], key)
            self._graph.add_node(new_atom_storage)
            self.atom_by_label[key] = new_atom_storage

            # Check for bonded atoms already in the graph
            for neighbor_key in reversed(sorted_keys):
                if neighbor_key not in self.atom_by_label:
                    continue

                # check if atoms are already connected on the graph
                neighbor_storage = self.atom_by_label[neighbor_key]
                if nx.has_path(self._graph, new_atom_storage, neighbor_storage):
                    continue

                # check if atoms are connected in the molecule
                atom2 = mol.get_atom_by_index(smirks_atoms[neighbor_key])
                bond = mol.get_bond_by_atoms(atom1, atom2)

                if bond is not None: # Atoms are connected add edge
                    bond_smirks = (neighbor_key, key)
                    bond_storage = self.BondStorage([bond], bond_smirks)
                    self.bond_by_label[bond_smirks] = bond_storage
                    self._graph.add_edge(new_atom_storage,
                                         neighbor_storage,
                                         bond=bond_storage)

        for atom_label, atom_index in smirks_atoms.items():
            atom = mol.get_atom_by_index(atom_index)
            storage = self.atom_by_label[atom_label]
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
            form {atom index: label} for this molecule
        """
        # if layers is 0 there are no more atoms to add so end the recursion
        if layers == 0:
            return

        # find atom neighbors that are not already included in SMIRKS indexed atoms
        atom_neighbors = [(a, mol.get_bond_by_atoms(a,atom)) for a in atom.get_neighbors() \
                          if a.get_index() not in idx_dict]

        # get the smirks indices already added to the storage
        # This includes all previous layers since the idx_dict is updated as you go
        storage_labels = [e for k,e in idx_dict.items()]

        # similar to atoms find neighbors already in the graph that haven't already been used
        storage_neighbors = [(s, self.get_connecting_bond(s, storage)) for s in self.get_neighbors(storage) \
                             if s.label not in storage_labels]

        new_pairs = list()
        # If the storage doesn't have any neighbors, add storage
        # Make new storages for all neighbors
        if len(storage_neighbors) == 0:
            min_smirks = storage.label * 10
            if min_smirks > 0:
                min_smirks = min_smirks * -1

            for a, b in atom_neighbors:
                new_bond_smirks = (storage.label, min_smirks)

                adding_new_storage = self.add_atom(a,b,storage,
                                                   min_smirks, new_bond_smirks)

                idx_dict[a.get_index()] = min_smirks
                self.atom_by_label[min_smirks] = adding_new_storage
                min_smirks -= 1
                new_pairs.append((a, adding_new_storage))

        else:
            # pair up the atoms with their storage
            pairs = self.find_pairs(atom_neighbors, storage_neighbors)
            for new_atom, new_bond, new_storage_atom, new_storage_bond in pairs:
                if new_storage_atom is None:
                    continue
                # add atom and bond information to the storage
                new_storage_atom.add_atom(new_atom)
                new_storage_bond.add_bond(new_bond)
                new_pairs.append((new_atom, new_storage_atom))
                idx_dict[new_atom.get_index()] = new_storage_atom.label

        # Repeat for the extra layers
        if layers == 'all':
            new_layers = 'all'
        else:
            new_layers = layers - 1
            if new_layers == 0:
                return

        for new_atom, new_storage in new_pairs:
            self._add_layers(mol, new_atom, new_storage, new_layers, idx_dict)

    def find_pairs(self, atoms_and_bonds, storages):
        """
        Find pairs is used to determine which current AtomStorage from storages
        atoms should be paired with.
        This function takes advantage of the maximum scoring function in networkx
        to find the pairing with the highest "score".
        Scores are determined using functions in the atom and bond storage objects
        that compare those storages to the new atom or bond.

        If there are less atoms than storages then the atoms with the lowest pair are
        assigned a None pairing.

        Parameters
        ----------
        atoms_and_bonds: list of tuples in form (chemper Atom, chemper Bond)
        storages: list of tuples in form (AtomStorage, BondStorage)

        Returns
        -------
        pairs: list of tuples
            pairs of atoms and storage objects that are most similar,
            note these are only the Atom and AtomStorage objects since the bonds and
        """
        pairs = list()

        g = nx.Graph()

        atom_dict = dict()
        storage_dict = dict()

        # create a bipartite graph with atoms/bonds on one side
        for idx, (a, b) in enumerate(atoms_and_bonds):
            g.add_node(idx+1, bipartite=0)
            atom_dict[idx+1] = (a,b)
        # and atom/bond storage objects on the other
        for idx, (s,sb) in enumerate(storages):
            g.add_node((idx*-1)-1, bipartite=1)
            storage_dict[(idx*-1)-1] = (s,sb)

        # Fill in the weight on each edge of the graph using the compare_atom/bond functions
        for a_idx, (a, b) in atom_dict.items():
            for s_idx, (sa, sb) in storage_dict.items():
                score = sa.compare_atom(a)
                score += sb.compare_bond(b) / 10.
                print(sa.as_smirks(), sb.as_smirks(), a.atomic_number(), score)
                g.add_edge(a_idx,s_idx,weight=score+1)

        # calculate maximum matching, that is the pairing of atoms/bonds to
        # storage objects that leads the the highest overall score
        matching = nx.algorithms.max_weight_matching(g,maxcardinality=False)
        # track the atoms assigned a paired storage object
        pair_set = set()

        for idx_1, idx_2 in matching:
            if idx_1 in atom_dict:
                (a,b) = atom_dict[idx_1]
                (sa,sb) = storage_dict[idx_2]
                pair_set.add(idx_1)
            else:
                (a,b) = atom_dict[idx_2]
                (sa,sb) = storage_dict[idx_1]
                pair_set.add(idx_2)
            pairs.append((a, b, sa, sb))

        for a_idx, (a,b) in atom_dict.items():
            if a_idx not in pair_set:
                pairs.append((a, b, None, None))

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
            This is a list of dictionaries of the form [{smirks index: atom index}]
            each atom (by index) in the dictionary will be added the relevant
            AtomStorage by smirks index
        """
        for smirks_atoms in smirks_atoms_list:
            atom_dict = dict()
            for key, atom_index in smirks_atoms.items():
                atom_dict[atom_index] = key
                atom1 = mol.get_atom_by_index(atom_index)
                self.atom_by_label[key].add_atom(atom1)

                for neighbor_key, neighbor_index in smirks_atoms.items():
                    # check for connecting bond
                    atom2 = mol.get_atom_by_index(neighbor_index)
                    bond = mol.get_bond_by_atoms(atom1, atom2)
                    if bond is not None and (neighbor_key, key) in self.bond_by_label:
                        bond_smirks = (neighbor_key, key)
                        self.bond_by_label[bond_smirks].add_bond(bond)

            for atom_label, atom_index in smirks_atoms.items():
                atom = mol.get_atom_by_index(atom_index)
                storage = self.atom_by_label[atom_label]
                self._add_layers(mol, atom, storage, self.layers, atom_dict)

