"""
fragment_graph.py

ChemPerGraph are a class for tracking molecular fragments based on information about the atoms and bonds they contain.
You can combine them or take a difference in order to find distinguishing characteristics in a set of clustered
molecular sub-graphs.

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

# TODO: import mol_toolkit?
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

            :param atom: chemper atom object
            :param smirks_index: integer for labeling this atom in a SMIRKS
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
            :return: smirks pattern for this atom
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
            :return: string of how this bond should appear in a SMIRKS string
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
        Initialize ChemPerGraph
        """
        self._graph = nx.Graph()
        self.atom_by_smirks_index = dict() # stores a dictionary of atoms with smirks_index
        self.bond_by_smirks_index = dict() # stores a dictionary of bonds with smirks_index

    def as_smirks(self):
        """
        :return: a SMIRKS string matching the exact atom and bond information stored
        """

        # If no atoms have been added
        if len(self._graph.nodes()) == 0:
            return None

        if self.atom_by_smirks_index:
            min_smirks = min(self.atom_by_smirks_index)
            init_atom = self.atom_by_smirks_index[min_smirks]
        else:
            init_atom = self.get_atoms()[0]

        # sort neighboring atoms to keep consist output
        neighbors = sorted(self.get_neighbors(init_atom), key=lambda a: a.as_smirks())
        return self._as_smirks(init_atom, neighbors)

    def _as_smirks(self, init_atom, neighbors):

        smirks = init_atom.as_smirks()

        for idx, neighbor in enumerate(neighbors):
            bond = self.get_connecting_bond(init_atom, neighbor)
            bond_smirks = bond.as_smirks()

            new_neighbors = self.get_neighbors(neighbor)
            new_neighbors.remove(init_atom)

            atom_smirks = self._as_smirks(neighbor, new_neighbors)

            if idx < len(neighbors) - 1:
                smirks += '(' + bond_smirks + atom_smirks + ')'
            else:
                smirks += bond_smirks + atom_smirks

        return smirks

    def get_atoms(self):
        """
        :return: list of all AtomStorage objects in graph
        """
        return list(self._graph.nodes())

    def get_connecting_bond(self, atom1, atom2):
        """
        :param atom1: AtomStorage object in this graph
        :param atom2: AtomStorage object in this graph
        :return: bond connecting them
        """
        return self._graph.get_edge_data(atom1, atom2)['bond']

    def get_bonds(self):
        """
        :return: list of all BondStorage objects in graph
        """
        return [data['bond'] for a1, a2, data in self._graph.edges(data=True)]

    def get_neighbors(self, atom):
        """
        :param atom: an AtomStorage object
        :return: neighboring AtomStorage objects
        """
        return list(self._graph.neighbors(atom))

    def add_atom(self, new_atom, new_bond=None, bond_to_atom=None,
                 new_smirks_index=None, new_bond_index=None):
        """
        :param new_atom: a chemper mol_toolkit atom object
        :param new_bond: a chemper mol_toolkit bond object
        :param bond_to_atom: AtomStorage object to connect this bond to
        :param new_smirks_index: integer for SMIRKS indexing the new atom
        :param new_bond_index: integer for 'smirks_index' on bond storage object
        :return: new AtomStorage object or None if atom wasn't added
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
    def __init__(self, mol, smirks_atoms, layers=0):
        """
        Initialize a ChemPerGraph from a molecule and a set of indexed atoms

        :param mol: chemper mol
        :param smirks_atoms: dictionary of the form {smirks_index: atom_index}
        :param layers: how many atoms out from the smirks indexed atoms do you wish save (default=0)
                       'all' will lead to all atoms in the molecule being specified (not recommended)
        """
        ChemPerGraph.__init__(self)

        self.mol = mol
        self.atom_by_index = dict()
        self._add_smirks_atoms(smirks_atoms)
        for smirks_key, atom_storage in self.atom_by_smirks_index.items():
            self._add_layers(atom_storage, layers)

    def _add_smirks_atoms(self, smirks_atoms):
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
