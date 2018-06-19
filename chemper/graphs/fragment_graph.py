"""
fragment_graph.py

ChemPerGraph are a class for tracking molecular fragments based on information about the atoms and bonds they contain.
You can combine them or take a difference in order to find distinguishing characteristics in a set of clustered
molecular sub-graphs.

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

import networkx as nx


class ChemPerGraph(object):
    """
    ChemPerGraphs are a graph based class for storing atom and bond information.
    They use the chemper.mol_toolkits Atoms, Bonds, and Mols
    """
    class AtomStorage(object):
        """
        AtomStorage tracks information about an atom
        """
        def __init__(self, atom=None, smirks_index=None):
            """
            Initializes AtomStorage based on a provided atom

            Parameters
            ----------
            atom: chemper Atom object
            smirks_index: int
                integer for labeling this atom in a SMIRKS
                or if negative number just used to track the atom locally
            """
            self.atom = atom

            if atom is None:
                self.atomic_number = None
                self.aromatic = None
                self.charge = None
                self.hydrogen_count = None
                self.connectivity = None
                self.ring_connectivity = None
                self.min_ring_size = None
                self.atom_index = None

            else:
                self.atomic_number = atom.atomic_number()
                self.aromatic = atom.is_aromatic()
                self.charge = atom.formal_charge()
                self.hydrogen_count = atom.hydrogen_count()
                self.connectivity = atom.connectivity()
                self.ring_connectivity = atom.ring_connectivity()
                self.min_ring_size = atom.min_ring_size()
                self.atom_index = atom.get_index()

            self.smirks_index = smirks_index

        def as_smirks(self):
            """
            Returns
            -------
            smirks: str
                how this atom would be represented in a SMIRKS string
            """
            if self.atom is None:
                if self.smirks_index is None or self.smirks_index <= 0:
                    return '[*]'
                return '[*:%i]' % self.smirks_index

            aromatic = 'a' if self.aromatic else 'A'
            if self.charge >= 0:
                charge = '+%i' % self.charge
            else:
                charge = '%i' % self.charge

            base_smirks = '#%i%sH%iX%ix%ir%i%s' % (self.atomic_number,
                                                   aromatic,
                                                   self.hydrogen_count,
                                                   self.connectivity,
                                                   self.ring_connectivity,
                                                   self.min_ring_size,
                                                   charge)

            if self.smirks_index is None or self.smirks_index <= 0:
                return '[%s]' % base_smirks

            return '[%s:%i]' % (base_smirks, self.smirks_index)

    class BondStorage(object):
        """
        BondStorage tracks information about a bond
        """
        def __init__(self, bond=None, smirks_index=None):
            """
            Parameters
            ----------
            bond: chemper Bond object
            smirks_index: int or float
                Bonds don't have SMIRKS indices so this is only used for internal
                tracking of the object.
            """
            if bond is None:
                self.order = None
                self.ring = None
                self.bond_index = None
            else:
                self.order = bond.get_order()
                self.ring = bond.is_ring()
                self.bond_index = bond.get_index()

            self._bond = bond
            self.smirks_index = smirks_index

        def as_smirks(self):
            """
            Returns
            -------
            SMIRKS: str
                how this bond should appear in a SMIRKS string
            """
            if self.ring is None:
                ring = ''
            elif self.ring:
                ring = '@'
            else:
                ring = '!@'

            order = {1:'-', 1.5:':', 2:'=', 3:'#', None:'~'}.get(self.order)

            return order+ring

    def __init__(self):
        """
        Initialize empty ChemPerGraph
        """
        self._graph = nx.Graph()
        self.atom_by_smirks_index = dict() # stores a dictionary of atoms with smirks_index
        self.bond_by_smirks_index = dict() # stores a dictionary of bonds with smirks_index

    def as_smirks(self):
        """
        Returns
        -------
        SMIRKS: str
            a SMIRKS string matching the exact atom and bond information stored
        """

        # If no atoms have been added
        if len(self._graph.nodes()) == 0:
            return None

        if self.atom_by_smirks_index:
            # sometimes we use negative numbers for internal indexing
            # the first atom in a smirks pattern should be based on actual smirks indices (positive)
            min_smirks = min([k for k in self.atom_by_smirks_index.keys() if k > 0])
            init_atom = self.atom_by_smirks_index[min_smirks]
        else:
            init_atom = self.get_atoms()[0]

        # sort neighboring atoms to keep consist output
        neighbors = sorted(self.get_neighbors(init_atom), key=lambda a: a.smirks_index)
        return self._as_smirks(init_atom, neighbors)

    def _as_smirks(self, init_atom, neighbors):
        """
        This is an internal/private method used to add all AtomStorage to the SMIRKS pattern

        Parameters
        ----------
        init_atom: AtomStorage object
            current atom
        neighbors: list of AtomStorage objects
            list of neighbor atoms you wanted added to the SMIRKS pattern

        Returns
        -------
        SMIRKS: str
            This graph as a SMIRKS string
        """

        smirks = init_atom.as_smirks()

        for idx, neighbor in enumerate(neighbors):
            bond = self.get_connecting_bond(init_atom, neighbor)
            bond_smirks = bond.as_smirks()

            new_neighbors = sorted(self.get_neighbors(neighbor), key=lambda a: a.smirks_index)
            new_neighbors.remove(init_atom)

            atom_smirks = self._as_smirks(neighbor, new_neighbors)

            if idx < len(neighbors) - 1:
                smirks += '(' + bond_smirks + atom_smirks + ')'
            else:
                smirks += bond_smirks + atom_smirks

        return smirks

    def get_atoms(self):
        """
        Returns
        -------
        atoms: list of AtomStorage objects
            all atoms stored in the graph
        """
        return list(self._graph.nodes())

    def get_connecting_bond(self, atom1, atom2):
        """
        Parameters
        ----------
        atom1: AtomStorage object
        atom2: AtomStorage object

        Returns
        -------
        bond: BondStorage object
            bond between the two given atoms or None if not connected
        """
        bond = self._graph.get_edge_data(atom1, atom2)
        if bond is not None:
            return bond['bond']
        return None

    def get_bonds(self):
        """
        Returns
        -------
        bonds: list of BondStorage objects
            all bonds stored as edges in this graph
        """
        return [data['bond'] for a1, a2, data in self._graph.edges(data=True)]

    def get_neighbors(self, atom):
        """
        Parameters
        ----------
        atom: an AtomStorage object

        Returns
        -------
        atoms: list of AtomStorage objects
            list of atoms one bond (edge) away from the given atom
        """
        return list(self._graph.neighbors(atom))

    def add_atom(self, new_atom, new_bond=None, bond_to_atom=None,
                 new_smirks_index=None, new_bond_index=None):
        """
        Expand the graph by adding one new atom including relevant bond

        Parameters
        ----------
        new_atom: a chemper Atom object
        new_bond: a chemper Bond object
        bond_to_atom: AtomStorage object
            This is where you want to connect the new atom, required if the graph isn't empty
        new_smirks_index: int
            (optional) index for SMIRKS or internal storage if less than zero
        new_bond_index: int or float
            (optional) index used to track bond storage

        Returns
        -------
        AtomStorage: AtomStorage object or None
            If the atom was successfully added then the AtomStorage object is returned
            None is returned if the atom wasn't able to be added
        """
        if bond_to_atom is None and len(self.get_atoms()) > 0:
            return None

        new_atom_storage = self.AtomStorage(new_atom, smirks_index=new_smirks_index)
        self._graph.add_node(new_atom_storage)
        if new_smirks_index is not None and new_smirks_index > 0:
            self.atom_by_smirks_index[new_smirks_index] = new_atom_storage

        # This is the first atom added to the graph
        if bond_to_atom is None:
            return new_atom_storage

        new_bond_storage = self.BondStorage(new_bond, new_bond_index)
        self.bond_by_smirks_index[new_bond_index] = new_bond_storage

        self._graph.add_edge(bond_to_atom, new_atom_storage, bond = new_bond_storage)
        return new_atom_storage


class ChemPerGraphFromMol(ChemPerGraph):
    """
    Creates a ChemPerGraph from a chemper Mol object
    """
    def __init__(self, mol, smirks_atoms, layers=0):
        """
        Parameters
        ----------
        mol: chemper Mol
        smirks_atoms: dict
            dictionary of the form {smirks_index: atom_index}
        layers: int or 'all'
            how many atoms out from the smirks indexed atoms do you wish save (default=0)
            'all' will lead to all atoms in the molecule being specified
        """
        ChemPerGraph.__init__(self)

        self.mol = mol
        self.atom_by_index = dict()
        self._add_smirks_atoms(smirks_atoms)
        for smirks_key, atom_storage in self.atom_by_smirks_index.items():
            self._add_layers(atom_storage, layers)

    def _add_smirks_atoms(self, smirks_atoms):
        """
        private function for adding atoms to the graph

        Parameters
        ----------
        smirks_atoms: dict
            dictionary of the form {smirks_index: atom_index}
        """
        # add all smirks atoms to the graph
        for key, atom_index in smirks_atoms.items():
            atom1 = self.mol.get_atom_by_index(atom_index)
            new_atom_storage = self.AtomStorage(atom1, key)
            self._graph.add_node(new_atom_storage)
            self.atom_by_smirks_index[key] = new_atom_storage
            self.atom_by_index[atom_index] = new_atom_storage
            # Check for bonded atoms already in the graph
            for neighbor_key, neighbor_index in smirks_atoms.items():
                if not neighbor_key in self.atom_by_smirks_index:
                    continue

                # check if atoms are already connected on the graph
                neighbor_storage = self.atom_by_smirks_index[neighbor_key]
                if nx.has_path(self._graph, new_atom_storage, neighbor_storage):
                    continue

                # check if atoms are connected in the molecule
                atom2 = self.mol.get_atom_by_index(neighbor_index)
                bond = self.mol.get_bond_by_atoms(atom1, atom2)

                if bond is not None: # Atoms are connected add edge
                    bond_index = max(neighbor_key, key)-1
                    bond_storage = self.BondStorage(bond, bond_index)
                    self.bond_by_smirks_index[bond_index] = bond_storage
                    self._graph.add_edge(new_atom_storage,
                                         self.atom_by_smirks_index[neighbor_key],
                                         bond=bond_storage)

    # TODO: I could probably do this with a while loop, is that better?
    def _add_layers(self, atom_storage, add_layer):
        """
        private function for expanding beyond the initial SMIRKS atoms.
        For now this is recursive so the input is:

        Parameters
        ----------
        atom_storage: AtomStorage object
            atom whose's neighbors you currently need to add
        add_layer: int
            how many more layers need to be added
        """
        if add_layer == 0:
            return

        new_smirks_index = min(1, atom_storage.smirks_index) - 1

        for new_atom in atom_storage.atom.get_neighbors():
            if new_atom.get_index() in self.atom_by_index:
                continue

            new_bond = self.mol.get_bond_by_atoms(atom_storage.atom, new_atom)
            new_storage = self.add_atom(new_atom, new_bond, atom_storage,
                                        new_smirks_index, new_smirks_index)
            self.atom_by_index[new_atom.get_index()] = new_storage
            if add_layer == 'all':
                self._add_layers(new_storage, add_layer)
            elif add_layer > 1:
                self._add_layers(new_storage, add_layer-1)
