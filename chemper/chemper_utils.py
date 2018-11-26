"""
utils.py

This file provides simple functions that might be of use to people using the chemper package

"""

import os

def get_data_path(relative_path, package='chemper'):
    """
    Returns the absolute path to a specified relative path inside data directory of a given package
    there for this will find a file at:

    [absolute path to Python packages]/[package(chemper]/data/[filename]

    Parameters
    ----------
    relative_path: str
              filename or relative path to a file stored in the data/ directory
              of a given package
    package: str
             name of the Python package that should contain the specified file
             default: "chemper"

    Returns
    -------
    path: str
          absolute path to filename in specified package
    """
    from pkg_resources import resource_filename
    fn = resource_filename(package, os.path.join('data', relative_path))

    if not os.path.exists(fn):
        raise IOError("Sorry! The absolute path %s was not found, that is %s is not in the data directory for %s" \
                          % (fn, relative_path, package))

    return fn

def get_full_path(relative_path, package="chemper"):
    """
    This function checks to see if a file/path at this relative path is available locally
    if it is it returns the absolute path to that file.
    Otherwise it checks if that relative path is located at [package location]/data/relative path
    it will raise an IOError if neither location exists.

    Parameters
    ----------
    relative_path: str
                   relative path to a file available locally or in package/data
    package: str
             package name with a data directory the file might be located in

    Returns
    -------
    abs_path: str
              The absolute path to a local file if available, otherwise the abolsute path
              to [package]/data/[relative_path]
    """
    if os.path.exists(relative_path):
        return os.path.abspath(relative_path)
    return get_data_path(relative_path, package)

# =======================================
# Check SMIRKS validity
# =======================================

def is_valid_smirks(smirks):
    """
    This function checks if a given smirks is valid.

    Parameters
    ----------
    smirks: SMIRKS (or SMARTS) string

    Returns
    -------
    is_valid: boolean
              is the provided SMIRKS a valid pattern
    """
    from chemper.mol_toolkits.mol_toolkit import MolFromSmiles
    mol = MolFromSmiles('C')
    try:
        mol.smirks_search(smirks)
        return True
    except ValueError:
        return False

# ===================================================================
# custom classes and functions for matching sets of SMIRKS to molecules
# ===================================================================
"""
custom_dicts

This custom dictionaries are taken from
openforcefield.typing.engines.smirnoff.forcefield
Like other openforcefield typing tools, I wanted to be
able to test ideas for chemper without adding dependencies
assuming things work and openforcefield all gets
merged to support OE/RDK I will just import these from there.
"""


import collections
class TransformedDict(collections.MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

    def __str__(self): return str(self.store)


class ValenceDict(TransformedDict):
    """Enforce uniqueness in atom indices"""
    def __keytransform__(self, key):
        """Reverse tuple if first element is larger than last element."""
        # Ensure key is a tuple.
        key = tuple(key)
        # Reverse the key if the first element is bigger than the last.
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key


class ImproperDict(TransformedDict):
    """Symmetrize improper torsions"""
    def __keytransform__(self,key):
        """Reorder tuple in numerical order except for element[1] which is the central atom; it retains its position."""
        # Ensure key is a tuple
        key = tuple(key)
        # Retrieve connected atoms
        connectedatoms = [key[0], key[2], key[3]]
        # Sort connected atoms
        connectedatoms.sort()
        # Re-store connected atoms
        key = tuple( [connectedatoms[0], key[1], connectedatoms[1], connectedatoms[2]])
        return key


def get_typed_molecules(smirks_list, input_molecules):
    """
    Creates a dictionary assigning a typename
    for each set of atom indices in each molecule

    Parameters
    ----------
    smirks_list: list of tuples in the form (label, smirks)
    input_molecules: list of Mols
        These can be molecules from OpenEye, RDKit, or ChemPer

    Returns
    -------
    typeDict: embedded dictionary
        keys: SMILES string for each molecule
            keys: tuple of indices assigned a parameter type
    """
    from chemper.mol_toolkits.mol_toolkit import Mol
    molecules = [Mol(m) for m in input_molecules]
    type_dict = dict()
    for mol_idx, mol in enumerate(molecules):
        type_dict[mol_idx] = {}
        for [label, smirks] in smirks_list:
            matches = get_smirks_matches(mol, smirks)
            for match in matches:
                type_dict[mol_idx][match] = label

    return type_dict

def create_tuples_for_clusters(smirks_list, molecules):
    """
    This function is used to get clusters of molecular fragments based
    on input SMIRKS.

    For example, let's assume you wanted to find all of the
    atoms that match this SMIRKS list
    'any', '[*:1]~[*:2]'
    'single', '[*:1]-[*:2]'

    In this case, the "any" bond would match all bonds, but then
    the "single" would match all single bonds.
    If you were typing Ethene (C=C) then you expect the double bond
    between atoms 0 and 1 to match any bond and all C-H bonds to match single.

    The output in this case would be:
    [ ('any', [[ (0, 1) ]] ),
      ('single', [[ (0, 2), (0, 3), (1,4), (1,5) ]] )
    ]

    Parameters
    ----------
    smirks_list: list of tuples
        This is a list of tuples with the form (label, SMIRKS)
    molecules: list of Mols
        This is a list of molecules you want to type with this list of SMIRKS
        They can be from any toolkit chemper supports (currently OpenEye and RDKit)

    Returns
    -------
    type_list: list of atom index tuples for each molecule for each type
        These are a list of tuples with the form (label, [ [ (atom indices by molecule)], ]),...

    """
    ordered_labels = [l for l, s in smirks_list]
    type_dict = get_typed_molecules(smirks_list, molecules)
    # start by getting things in the form
    # {label: {mol_idx: list of dictionary of indices} }
    label_dict = dict()
    for mol_idx, mol_dict in type_dict.items():
        for match, label in mol_dict.items():
            if label not in label_dict:
                label_dict[label] = [list() for i in range(len(molecules))]

            label_dict[label][mol_idx].append(match)

    # now we need to resort the mol_lists
    final_list = [ (l, label_dict[l]) for l in ordered_labels if l in label_dict]
    return final_list


def get_smirks_matches(mol, smirks):
    """
    Gets atom indices for a smirks string in a given molecule

    Parameters
    ----------
    mol : a chemper Mol
    smirks : str
             SMIRKS pattern being matched to the molecule

    Returns
    --------
    matches: list of tuples
        atom indices for labeled atom in the smirks
    """
    from chemper.graphs.environment import ChemicalEnvironment

    env = ChemicalEnvironment(smirks)
    if env.getType().lower() == 'impropertorsion':
        matches = ImproperDict()
    else:
        matches = ValenceDict()

    for match in mol.smirks_search(smirks):
        smirks_indices = sorted(list(match.keys()))
        atom_indices = tuple([match[s].get_index() for s in smirks_indices])
        matches[atom_indices] = ''

    return matches.keys()


# ===================================================================
# classes for evaluating lists of SMIRKS patterns
# ===================================================================

def score_match_reference(current_assignments, ref_assignments):
    """

    Acknowledgement:
    Josh Fass <josh.fass@choderalab.org>
    created this scoring algorithm using the bipartite graph for smarty,
    which has been published (doi: ...) # TODO: add DOI
    """
    import networkx as nx
    total_counts = dict()
    total = 0
    for mol_idx, ref_dict in ref_assignments.items():
        for indices, label in ref_dict.items():
            if label not in total_counts:
                total_counts[label] = 0
            total += 1
            total_counts[label] +=1

    cur_labels = set()
    ref_labels = set()
    # check for missing tuples in dictionaries
    for mol_idx, current_dict in current_assignments.items():
        cur_keys = set(current_dict.keys())
        ref_keys = set(ref_assignments[mol_idx].keys())
        # check if there are atom index sets in references not in current
        if ref_keys - cur_keys:
            # TODO:
            return None, -1

        # store a set of labels
        c_labs = [e for k,e in current_dict.items()]
        r_labs = [e for k,e in ref_assignments[mol_idx].items()]
        cur_labels = cur_labels.union(set(c_labs))
        ref_labels = ref_labels.union(set(r_labs))

    # Create bipartite graph (U,V,E) matching current types U with
    # reference types V via edges E with weights equal to number of types in common.
    # create a dictionary for each possible pair of current and reference labels
    types_in_common = dict()
    for c_lab in cur_labels:
        for r_lab in ref_labels:
            types_in_common[(c_lab, r_lab)] = 0

    # up the count by +1 for each time a current and reference type matches the same set of atoms
    for mol_idx, index_dict in ref_assignments.items():
        for indices, r_lab in index_dict.items():
            c_lab = current_assignments[mol_idx][indices]
            types_in_common[(c_lab, r_lab)] += 1

    # Actually make the graph
    graph = nx.Graph()
    # Add current types
    for c_lab in cur_labels:
        graph.add_node(c_lab, bipartite=0)
    # add reference types
    for r_lab in ref_labels:
        graph.add_node(r_lab, bipartite=1)
    # add edges
    for c_lab in cur_labels:
        for r_lab in ref_labels:
            weight = types_in_common[(c_lab, r_lab)]
            graph.add_edge(c_lab, r_lab, weight=weight)

    mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)

    # Compute match dictionary and total number of matches.
    type_matches = list()
    total_type_matches = 0
    for lab1, lab2 in mate:
        counts = graph[lab1][lab2]['weight']
        total_type_matches += counts
        # labels come back in an arbitrary order, so determine if which is current/reference
        if lab1 in cur_labels:
            type_matches.append(( lab1, lab2, counts, total_counts[lab2]))
        else:
            type_matches.append( (lab2, lab1, counts, total_counts[lab1]))

    # compute fractional score:
    score = total_type_matches / total

    return type_matches, score

def match_reference(current_assignments, ref_assignments):
    """
    Determine if two sets of labels agree.
    This could be used to test if two sets of SMIRKS type a set of molecules
    in exactly the same way.

    This starts with two sets of "assignments" which are dictionaries
    like those created with get_typed_molecules.

    Let's imagine you wanted to type the molecule CF and CCl
    where most bonds are assigned X, but one bond in each molecule is
    assigned 'Y'. Then the reference_assignments would be:
    { 0: (0,1): 'X', (0,2): 'X', (0,3): 'X', (0,4): 'Y'},
      1: (0,1): 'X', (0,2): 'X', (0,3): 'X', (0,4): 'Y'},
    }

    Then you discover a set of SMIRKS patterns which match
    C-H bonds with one type and C-halogen have another.
    Then your current_assignments would be:
    { 0: (0,1): 'C-H', (0,2): 'C-H', (0,3): 'C-H', (0,4): 'C-hal.'},
      1: (0,1): 'C-H', (0,2): 'C-H', (0,3): 'C-H', (0,4): 'C-hal.'},
    }
    This would return a set tuples for the corresponding labels
    [('C-H', 'X'), ('C-hal.', 'Y')] and True for the two assigments matching.

    However, if your current set assigned different types to the two
    different carbon halogen bonds (C-F and C-Cl) then these two
    sets of assignments would not match and the result would be an empty set and False

    Parameters
    ----------
    current_assignments : dictionary of indices with current types
    ref_assignments : dictionary of indices with reference types

    Returns
    -------
    type_matches : set of tuples (current_label, reference_label)
        pair of current and reference labels, None if the sets don't match
    is_match : boolean
        Does the current assignment exactly match the reference?
    """
    import networkx as nx

    graph = nx.Graph()
    matches = set()
    r_labs = set()
    c_labs = set()

    for mol_idx, index_dict in ref_assignments.items():
        # check that mol_idx is in current_assignments
        if mol_idx not in current_assignments:
            return set(), False

        for indices, r_label in index_dict.items():
            # check this set of indices is in the current set
            if indices not in current_assignments[mol_idx]:
                return set(), False

            c_label = current_assignments[mol_idx][indices]
            # check if nodes exist for the two labels (and make them)
            if c_label not in graph:
                graph.add_node(c_label)
                c_labs.add('current_'+c_label)

            if r_label not in graph:
                graph.add_node(r_label)
                r_labs.add(r_label)

            # create an edge connecting the reference and current type
            graph.add_edge('current_'+c_label, r_label)
            matches.add((c_label, r_label))

    for c_label in c_labs:
        if len(graph.edges(c_label)) != 1:
            return set(), False

    for r_label in r_labs:
        if len(graph.edges(r_label)) != 1:
            return set(), False

    return matches, True


def check_smirks_to_reference(current_types, reference_assignments, molecules):
    """

    Parameters
    ----------
    current_types: list of tuples
        SMIRKS types in the form (label, smirks)
    reference_assignments: dictionary of tuples and labels
        This could be the output from get_typed_molecules
        the dictionary has the form:
        {mol_idx: {(atom indices tuple): label}, mol_idx2: {} }
    molecules: list of molecules
        These can be OpenEye, RDKit, or ChemPer molecules

    Returns
    -------
    agree: boolean
        Returns True if the SMIRKS type the molecules in the same way
        as the reference assignments
    """
    r_labs = [l for m,d in reference_assignments.items() for a,l in d.items()]
    if len(current_types) != len(set(r_labs)):
        return False

    current_assignments = get_typed_molecules(current_types, molecules)
    type_matches, matched = match_reference(current_assignments, reference_assignments)
    return matched


def check_smirks_agree(current_types, reference_types, molecules):
    """
    Checks if two lists of SMIRKS patterns type a list of molecules in the same way.

    Parameters
    ----------
    current_types: list of tuples
        SMIRKS types in the form (label, smirks)
    reference_types: list of tuples
        SMIRKS types in the form (label, smirks)
    molecules: list of molecules
        molecules to be typed with this patterns
        Can be any molecule type chemper supports (OpenEye, RDKit, or Chemper)

    Returns
    -------
    match: boolean
        True if the two sets of SMIRKS match the set of molecules in the exact same way
    """
    if len(current_types) != len(reference_types):
        return False

    current_assignments = get_typed_molecules(current_types, molecules)
    reference_assignments = get_typed_molecules(reference_types, molecules)

    # check if they agree
    type_matches, matched = match_reference(current_assignments, reference_assignments)
    return matched

