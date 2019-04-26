"""
fragment_graph.py

ChemPerGraph is a class for storing smirks decorators for a molecular fragment.
These can be used to convert a molecular sub-graph or an entire molecule into a SMIRKS
pattern with all decorators specified.

For example, imagine you want a SMIRKS for the carbon in methane, it would become:

"[#6AH4X4x0!r+0:1]"

with decorators:
#6: atomic number 6 for carbon
A: aliphatic (a would be aromatic)
H4: a total hydrogen count of 4, 4 neighbors are hydrogen
X4: connectivity of 4, that is number of neighbors, not valence or sum of bond orders
x0: ring connectivity of 0, no ring bonds
!r: not in a ring, for atoms in a ring this decorator is `rn` where n is the size of the smallest ring
+0: 0 formal charge

To the best of the authors knowledge, this is the first open source tool capable
of converting a molecule (or sub-graph) into a detailed SMIRKS pattern.

AUTHORS:

Caitlin C. Bannan <bannanc@uci.edu>, Mobley Group, University of California Irvine
"""

import networkx as nx
from functools import total_ordering
from chemper.mol_toolkits import mol_toolkit


@total_ordering
class ChemPerGraph(object):
    """
    ChemPerGraphs are a graph based class for storing atom and bond information.
    They use the chemper.mol_toolkits Atoms, Bonds, and Mols
    """
    @total_ordering
    class AtomStorage(object):
        """
        AtomStorage tracks information about an atom
        """
        def __init__(self, atom=None, label=None):
            """
            Initializes AtomStorage based on a provided atom

            Parameters
            ----------
            atom: chemper Atom object
            label: int
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

            self.label = label

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

        def __eq__(self, other): return self.as_smirks() == other.as_smirks() and self.label == other.label

        def __hash__(self): return id(self)

        def __str__(self): return self.as_smirks()

        def as_smirks(self, compress=False):
            """
            Returns
            -------
            smirks: str
                how this atom would be represented in a SMIRKS string
            """
            if self.atom is None:
                if self.label is None or self.label <= 0:
                    return '[*]'
                return '[*:%i]' % self.label

            aromatic = 'a' if self.aromatic else 'A'
            if self.charge >= 0:
                charge = '+%i' % self.charge
            else:
                charge = '%i' % self.charge
            if self.min_ring_size == 0:
                ring = '!r'
            else:
                ring = 'r%i' % self.min_ring_size

            if compress:
                base_smirks = "#%i" % self.atomic_number
            else:
                base_smirks = '#%i%sH%iX%ix%i%s%s' % (self.atomic_number,
                                                      aromatic,
                                                      self.hydrogen_count,
                                                      self.connectivity,
                                                      self.ring_connectivity,
                                                      ring,
                                                      charge)

            if self.label is None or self.label <= 0:
                return '[%s]' % base_smirks

            return '[%s:%i]' % (base_smirks, self.label)

    @total_ordering
    class BondStorage(object):
        """
        BondStorage tracks information about a bond
        """
        def __init__(self, bond=None, label=None):
            """
            Parameters
            ----------
            bond: chemper Bond object
            label: int or float
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
            self.label = label

        def __str__(self): return self.as_smirks()

        def __lt__(self, other):
            if self.as_smirks() == other.as_smirks():
                return self.label < other.label
            return self.as_smirks() < other.as_smirks()

        def __eq__(self, other):
            return self.label == other.label and self.as_smirks() == other.as__smirks()

        def __hash__(self): return id(self)

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
        self.atom_by_label = dict() # stores a dictionary of atoms by label
        self.bond_by_label = dict() # stores a dictionary of bonds by label

    def __str__(self): return self.as_smirks()

    def __lt__(self, other): return self.as_smirks() < other.as_smirks()

    def __eq__(self, other): return self.as_smirks() == self.as_smirks()

    def __hash__(self): return id(self)

    def as_smirks(self, compress=False):
        """
        Parameters
        ----------
        compress: boolean
                  returns the shorter version of atom SMIRKS patterns
                  that is the atoms only include atomic numbers rather
                  than the full list of decorators
        Returns
        -------
        SMIRKS: str
            a SMIRKS string matching the exact atom and bond information stored
        """

        # If no atoms have been added
        if len(self._graph.nodes()) == 0:
            return None

        if self.atom_by_label:
            # sometimes we use negative numbers for internal indexing
            # the first atom in a smirks pattern should be based on actual smirks indices (positive)
            smirks_indices = [k for k in self.atom_by_label.keys() if k > 0]
            if len(smirks_indices) != 0:
                min_smirks = min(smirks_indices)
            else:
                min_smirks = min([k for k in self.atom_by_label.keys()])
            init_atom = self.atom_by_label[min_smirks]
        else:
            init_atom = self.get_atoms()[0]

        # sort neighboring atoms to keep consist output
        neighbors = sorted(self.get_neighbors(init_atom))
        return self._as_smirks(init_atom, neighbors, compress)

    def _as_smirks(self, init_atom, neighbors, compress=False):
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
        smirks = init_atom.as_smirks(compress)
        for idx, neighbor in enumerate(neighbors):
            bond = self.get_connecting_bond(init_atom, neighbor)
            bond_smirks = bond.as_smirks()

            new_neighbors = sorted(self.get_neighbors(neighbor))
            new_neighbors.remove(init_atom)

            atom_smirks = self._as_smirks(neighbor, new_neighbors,compress)

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

    def remove_atom(self, atom):
        """
        Removes the provided atom and all connected atoms
        """
        # if atom isn't in the graph, it can't be removed
        if atom not in self._graph.nodes():
            return False
        # if atom is "indexed" that is has a SMIRKS index > 0 it can't be removed
        if atom.label > 0:
            return False
        # remove specified atom
        self._graph.remove_node(atom)
        # find atoms on that "branch" of the molecule
        # we do this by looking for atoms that are no longer connected to
        # the base of the graph, where we consider the base a positively indexed atom
        ref_atom = [n for n in self._graph.nodes if n.label > 0][0]
        remove_atoms_list = list()
        for n in self._graph.nodes:
            if not nx.has_path(self._graph, n, ref_atom):
                remove_atoms_list.append(n)
        # remove the disconnected atoms
        self._graph.remove_nodes_from(remove_atoms_list)
        return True

    def add_atom(self, new_atom, new_bond=None, bond_to_atom=None,
                 new_label=None, new_bond_label=None):
        """
        Expand the graph by adding one new atom including relevant bond

        Parameters
        ----------
        new_atom: a chemper Atom object
        new_bond: a chemper Bond object
        bond_to_atom: AtomStorage object
            This is where you want to connect the new atom, required if the graph isn't empty
        new_label: int
            (optional) index for SMIRKS or internal storage if less than zero
        new_bond_label: int or float
            (optional) index used to track bond storage

        Returns
        -------
        AtomStorage: AtomStorage object or None
            If the atom was successfully added then the AtomStorage object is returned
            None is returned if the atom wasn't able to be added
        """
        if bond_to_atom is None and len(self.get_atoms()) > 0:
            return None

        new_atom_storage = self.AtomStorage(new_atom, label=new_label)
        self._graph.add_node(new_atom_storage)
        if new_label is not None:
            self.atom_by_label[new_label] = new_atom_storage

        # This is the first atom added to the graph
        if bond_to_atom is None:
            return new_atom_storage

        new_bond_storage = self.BondStorage(new_bond, new_bond_label)
        self.bond_by_label[new_bond_label] = new_bond_storage

        self._graph.add_edge(bond_to_atom, new_atom_storage, bond = new_bond_storage)
        return new_atom_storage


# ==============================================================================
# TODO: Isn't this the same thing as starting with a ChemPerGraph with mols=None
# and smirks_atoms=None as the default?
# ==============================================================================
class ChemPerGraphFromMol(ChemPerGraph):
    """
    Creates a ChemPerGraph from a chemper Mol object
    """
    def __init__(self, mol, smirks_atoms, layers=0):
        """
        Parameters
        ----------
        mol: Mol
            this can be a chemper mol or a molecule from any supported toolkit
            (currently OpenEye or RDKit)
        smirks_atoms: tuple of integers
            This is a tuple of the atom indices which will have SMIRKS indices.
            For example, if (1,2) is provided then the atom in molecule with indices
            1 and 2 will be used to create a SMIRKS with two indexed atoms.
        layers: int or 'all'
            how many atoms out from the smirks indexed atoms do you wish save (default=0)
            'all' will lead to all atoms in the molecule being specified
        """
        ChemPerGraph.__init__(self)

        self.mol = mol_toolkit.Mol(mol)
        self.atom_by_index = dict()
        self._add_smirks_atoms(smirks_atoms)
        keys = list(self.atom_by_label.keys())
        for smirks_key in keys:
            atom_storage = self.atom_by_label[smirks_key]
            self._add_layers(atom_storage, layers)

    def _add_smirks_atoms(self, smirks_atoms):
        """
        private function for adding atoms to the graph

        Parameters
        ----------
        smirks_atoms: tuple of integers
            This is a tuple of the atom indices which will have SMIRKS indices.
        """
        # add all smirks atoms to the graph
        for key, atom_index in enumerate(smirks_atoms, 1):
            atom1 = self.mol.get_atom_by_index(atom_index)
            new_atom_storage = self.AtomStorage(atom1, key)
            self._graph.add_node(new_atom_storage)
            self.atom_by_label[key] = new_atom_storage
            self.atom_by_index[atom_index] = new_atom_storage
            # Check for bonded atoms already in the graph
            for neighbor_key, neighbor_index in enumerate(smirks_atoms, 1):
                if not neighbor_key in self.atom_by_label:
                    continue

                # check if atoms are already connected on the graph
                neighbor_storage = self.atom_by_label[neighbor_key]
                if nx.has_path(self._graph, new_atom_storage, neighbor_storage):
                    continue

                # check if atoms are connected in the molecule
                atom2 = self.mol.get_atom_by_index(neighbor_index)
                bond = self.mol.get_bond_by_atoms(atom1, atom2)

                if bond is not None: # Atoms are connected add edge
                    bond_index = max(neighbor_key, key)-1
                    bond_storage = self.BondStorage(bond, bond_index)
                    self.bond_by_label[bond_index] = bond_storage
                    self._graph.add_edge(new_atom_storage,
                                         self.atom_by_label[neighbor_key],
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

        new_label = min(1, atom_storage.label) - 1

        for new_atom in atom_storage.atom.get_neighbors():
            if new_atom.get_index() in self.atom_by_index:
                continue

            new_bond = self.mol.get_bond_by_atoms(atom_storage.atom, new_atom)
            new_storage = self.add_atom(new_atom, new_bond, atom_storage,
                                        new_label, new_label)
            self.atom_by_index[new_atom.get_index()] = new_storage
            if add_layer == 'all':
                self._add_layers(new_storage, add_layer)
            elif add_layer > 1:
                self._add_layers(new_storage, add_layer-1)
