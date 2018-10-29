"""
cluster_graph.py

ClusterGraph are a class for tracking all possible smirks decorators in a group (or cluster)
of molecular fragments. Moving forward these will be used to find the minimum number of
smirks decorators that are required to have a set of smirks patterns that maintain
a given clustering of fragments.

# TODO: add specific example like the one at the top of fragment_graph

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

import networkx as nx
from chemper.graphs.fragment_graph import ChemPerGraph
from chemper.graphs.environment import ChemicalEnvironment as CE
from chemper.mol_toolkits import mol_toolkit


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
            min_ring_size = atom.min_ring_size()
            if min_ring_size == 0:
                ring = '!r'
            else:
                ring = 'r%i' % min_ring_size

            return (
                '#%i' % atom.atomic_number(),
                aromatic,
                'H%i' % atom.hydrogen_count(),
                'X%i' % atom.connectivity(),
                'x%i' % atom.ring_connectivity(),
                ring,
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

            if compress and len(self.decorators) > 1:
                base_smirks = self._compress_smirks()
            else:
                base_smirks = ','.join(sorted([''.join(l) for l in self.decorators]))

            if self.label is None or self.label <= 0:
                return '[%s]' % base_smirks

            return '[%s:%i]' % (base_smirks, self.label)

        def _sort_decs(self, dec_set, wild=True):
            """
            Parameters
            ----------
            dec_set: list like
                single set of atom decorators

            Returns
            -------
            sorted_dec_set: list
                same set of decorators sorted with atomic number or * first
            """
            atom_num = [i for i in dec_set if '#' in i]
            if len(atom_num) == 0 and wild:
                atom_num = ["*"]

            return atom_num + sorted(list(set(dec_set) - set(atom_num)))

        def _compress_smirks(self):
            """
            Returns
            -------
            smirks: str
                This SMIRKS is compressed with all common decorators and'd to
                the end of the pattern
            """
            set_decs = [set(d) for d in self.decorators]
            ands = set_decs[0]

            for d_set in set_decs:
                ands = ands & d_set

            or_sets = [self._sort_decs(d.difference(ands)) for d in set_decs]
            ors = [''.join(o) for o in or_sets]

            # add commas between ors
            base = ','.join(sorted(ors))
            # add and decorators
            if len(ands) > 0:
                base += ';'+ ';'.join(self._sort_decs(ands, wild=False))
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

        For the example, imagine we wanted to get a SMIRKS that
        would match the carbon-carbon bonds in ethane and propane.
        The carbon atoms are have indices (0,1) in ethane and (0,1) and (1,2)
        in propane. For this example, we will assume we also want to include
        the atoms one bond away from the indexed atoms (1 layer away).

        Parameters
        ----------
        mols: list of molecules (optional)
            these can be ChemPer mol's or molecule objects from
            any supported toolkit (current OpenEye or RDKit)

        smirks_atoms_lists: list of list of tuples (optional)
            There is a list of tuples for each molecule, where each tuple specifies
            a molecular fragment using the atoms' indices.
            In the ethane and propane example, the `smirks_atoms_lists` would be
                [ [ (0,1) ], [ (0,1), (1,2) ] ]
            with one carbon-carbon bond in ethane and two carbon-carbon bonds in propane

        layers: int (optional, default=0)
            layers specifies how many bonds away from the indexed atoms should be included in the
            the SMIRKS patterns.
            Instead of an int, the string 'all' would lead to all atoms in the molecules
            being included in the SMIRKS (not recommended)
        """
        ChemPerGraph.__init__(self)

        self.mols = list()
        self.smirks_atoms_lists = list()
        self.layers = layers
        self._symmetry_funct = self._no_symmetry

        if mols is not None:
            temp_mols = [mol_toolkit.Mol(m) for m in mols]
            if len(temp_mols) != len(smirks_atoms_lists):
                raise Exception('Number of molecules and smirks dictionaries should be equal')

            for idx, mol in enumerate(temp_mols):
                self.add_mol(mol, smirks_atoms_lists[idx])

    def as_smirks(self, compress=False):
        """
        Parameters
        ----------
        compress: boolean
                  returns the shorter version of atom SMIRKS patterns
                  that is atoms have decorators "anded" to the end rather than listed
                  in each set that are OR'd together.
                  For example "[#6AH2X3x0r0+0,#6AH1X3x0r0+0:1]-;!@[#1AH0X1x0r0+0]"
                  compresses to: "[#6H2,#6H1;AX3x0r0+0:1]-;!@[#1AH0X1x0r0+0]"

        Returns
        -------
        SMIRKS: str
            a SMIRKS string matching the exact atom and bond information stored
        """
        return ChemPerGraph.as_smirks(self, compress)

    def get_symmetry_funct(self, sym_label):
        """
        Parameters
        ----------
        sym_label: str or None
            type of symmetry, options which will change the way symmetry is
            handled in the graph are "bond", "angle", "ProperTorsion", and "ImproperTorsion"

        Returns
        -------
        symmetry_funct: function
            returns the function that should be used to handle the appropriate symmetry
        """
        if sym_label is None:
            return self._no_symmetry
        if sym_label.lower() == 'bond':
            return self._bond_symmetry
        if sym_label.lower() == 'angle':
            return self._angle_symmetry
        if sym_label.lower() == 'propertorsion':
            return self._proper_torsion_symmetry
        if sym_label.lower() == 'impropertorsion':
            return self._improper_torsion_symmetry
        return self._no_symmetry

    def add_mol(self, input_mol, smirks_atoms_list):
        """
        Expand the information in this graph by adding a new molecule

        Parameters
        ----------
        input_mol: chemper Mol object
        smirks_atoms_list: list of tuples
            This is a list of tuples with atom indices [ (indices), ... ]
        """
        mol = mol_toolkit.Mol(input_mol)

        if len(smirks_atoms_list) == 0:
            return

        if len(self.mols) == 0:
            self._add_first_smirks_atoms(mol, smirks_atoms_list[0])
            self._symmetry_funct = self.get_symmetry_funct(CE(self.as_smirks()).getType())
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
        smirks_atoms: tuple
            tuple of atom indices for the first atoms to add to the graph. i.e. (0, 1)
        """
        atom_dict = dict()
        for key, atom_index in enumerate(smirks_atoms, 1):
            atom_dict[atom_index] = key

            atom1 = mol.get_atom_by_index(atom_index)
            new_atom_storage = self.AtomStorage([atom1], key)
            self._graph.add_node(new_atom_storage)
            self.atom_by_label[key] = new_atom_storage

            # Check for bonded atoms already in the graph
            for neighbor_key in range(len(smirks_atoms), 0, -1):
                if neighbor_key not in self.atom_by_label:
                    continue

                # check if atoms are already connected on the graph
                neighbor_storage = self.atom_by_label[neighbor_key]
                if nx.has_path(self._graph, new_atom_storage, neighbor_storage):
                    continue

                # check if atoms are connected in the molecule
                atom2 = mol.get_atom_by_index(smirks_atoms[neighbor_key-1])
                bond = mol.get_bond_by_atoms(atom1, atom2)

                if bond is not None: # Atoms are connected add edge
                    bond_smirks = tuple(sorted([neighbor_key, key]))
                    bond_storage = self.BondStorage([bond], bond_smirks)
                    self.bond_by_label[bond_smirks] = bond_storage
                    self._graph.add_edge(new_atom_storage,
                                         neighbor_storage,
                                         bond=bond_storage)

        for atom_label, atom_index in enumerate(smirks_atoms, 1):
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
                new_bond_smirks = tuple(sorted([storage.label, min_smirks]))

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
            with tuples in the form (Atom, Bond, AtomStorage, BondStorage)
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
        mol: any Mol
        smirks_atoms_list: list of dicts
            This is a list of dictionaries of the form [{smirks index: atom index}]
            each atom (by index) in the dictionary will be added the relevant
            AtomStorage by smirks index
        """
        for smirks_atoms in smirks_atoms_list:
            atom_dict = dict()
            sorted_smirks_atoms = self._symmetry_funct(mol, smirks_atoms)
            for key, atom_index in enumerate(sorted_smirks_atoms, 1):
                atom_dict[atom_index] = key
                atom1 = mol.get_atom_by_index(atom_index)
                self.atom_by_label[key].add_atom(atom1)

                for neighbor_key, neighbor_index in enumerate(sorted_smirks_atoms, 1):
                    # check for connecting bond
                    atom2 = mol.get_atom_by_index(neighbor_index)
                    bond = mol.get_bond_by_atoms(atom1, atom2)
                    if bond is not None and (neighbor_key, key) in self.bond_by_label:
                        bond_smirks = tuple(sorted([neighbor_key, key]))
                        self.bond_by_label[bond_smirks].add_bond(bond)

            for atom_label, atom_index in enumerate(sorted_smirks_atoms, 1):
                atom = mol.get_atom_by_index(atom_index)
                storage = self.atom_by_label[atom_label]
                self._add_layers(mol, atom, storage, self.layers, atom_dict)

    def _no_symmetry(self, mol, smirks_atoms):
        return smirks_atoms

    def _bond_symmetry(self, mol, smirks_atoms):
        """
        Returns a tuple of two atom indices in the order that
        leads to the atoms that match with previously stored atoms.
        """
        # pair atoms and bonds
        atom1 = mol.get_atom_by_index(smirks_atoms[0])
        atom2 = mol.get_atom_by_index(smirks_atoms[1])
        bond = mol.get_bond_by_atoms(atom1, atom2)
        if bond is None:
            return smirks_atoms
        # Find potential storages for those atoms and bonds
        atoms_and_bonds = [(atom1, bond), (atom2, bond)]
        bond_storage = self.bond_by_label[(1,2)]
        storages = [
            (self.atom_by_label[1], bond_storage),
            (self.atom_by_label[2], bond_storage)
        ]
        pairs = self.find_pairs(atoms_and_bonds, storages)
        ordered_smirks_atoms = [p[0].get_index() for p in sorted(pairs, key=lambda x: x[2].label)]
        return tuple(ordered_smirks_atoms)

    def _angle_symmetry(self, mol, smirks_atoms):
        """
        Returns a tuple of three atom indices in the order that
        leads to the atoms that match with previously stored atoms.
        """
        # get all three atoms
        atom1 = mol.get_atom_by_index(smirks_atoms[0])
        atom2 = mol.get_atom_by_index(smirks_atoms[1])
        atom3 = mol.get_atom_by_index(smirks_atoms[2])
        # get both bonds
        bond1 = mol.get_bond_by_atoms(atom1, atom2)
        bond2 = mol.get_bond_by_atoms(atom2, atom3)
        if None in (bond1, bond2):
            return smirks_atoms
        # save atom and bond pairs that could be reordered
        atoms_and_bonds = [(atom1, bond1), (atom3, bond2)]
        # find current atom and bond storage
        storages = [
            (self.atom_by_label[1], self.bond_by_label[(1,2)]),
            (self.atom_by_label[3], self.bond_by_label[(2,3)])
        ]
        pairs = self.find_pairs(atoms_and_bonds, storages)
        order = [p[0].get_index() for p in sorted(pairs, key=lambda x: x[2].label)]
        return tuple((order[0], smirks_atoms[1], order[1]))

    def _proper_torsion_symmetry(self, mol, smirks_atoms):
        """
        Returns a tuple of four atom indices for a proper torsion
        reordered to match with previously stored atoms.
        """
        # get all four atoms
        atom1 = mol.get_atom_by_index(smirks_atoms[0])
        atom2 = mol.get_atom_by_index(smirks_atoms[1])
        atom3 = mol.get_atom_by_index(smirks_atoms[2])
        atom4 = mol.get_atom_by_index(smirks_atoms[3])
        # get two relevant bonds
        bond1 = mol.get_bond_by_atoms(atom1, atom2)
        bond3 = mol.get_bond_by_atoms(atom3, atom4)
        if None in (bond1, bond3):
            return smirks_atoms
        # make pairs
        atoms_and_bonds = [ (atom2, bond1), (atom3, bond3) ]
        # get atom and bond storages
        storages = [
            (self.atom_by_label[2], self.bond_by_label[(1,2)]),
            (self.atom_by_label[3], self.bond_by_label[(3,4)])
        ]
        pairs = self.find_pairs(atoms_and_bonds, storages)
        order = [p[0].get_index() for p in sorted(pairs, key=lambda x: x[2].label)]
        if order[0] == smirks_atoms[1]:
            return smirks_atoms
        temp = list(smirks_atoms)
        temp.reverse()
        return tuple(temp)

    def _improper_torsion_symmetry(self, mol, smirks_atoms):
        """
        Returns a tuple of four atom indices for an improper torsion
        reordered to match with previously stored atoms.
        """
        # get all four atoms
        atom1 = mol.get_atom_by_index(smirks_atoms[0])
        atom2 = mol.get_atom_by_index(smirks_atoms[1])
        atom3 = mol.get_atom_by_index(smirks_atoms[2])
        atom4 = mol.get_atom_by_index(smirks_atoms[3])
        # get all three bonds
        bond1 = mol.get_bond_by_atoms(atom1, atom2)
        bond2 = mol.get_bond_by_atoms(atom2, atom3)
        bond3 = mol.get_bond_by_atoms(atom2, atom4)
        if None in (bond1, bond2, bond3):
            return smirks_atoms
        # make pairs of atoms and bonds to be reordered
        atoms_and_bonds = [
            (atom1, bond1), (atom3, bond2), (atom4, bond3)
        ]
        # find current atom and bond storages
        storages = [
            (self.atom_by_label[1], self.bond_by_label[(1,2)]),
            (self.atom_by_label[3], self.bond_by_label[(2,3)]),
            (self.atom_by_label[4], self.bond_by_label[(2,4)])
        ]
        pairs = self.find_pairs(atoms_and_bonds, storages)
        order = [p[0].get_index() for p in sorted(pairs, key=lambda x: x[2].label)]
        return tuple((order[0], smirks_atoms[1], order[1], order[2]))