#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
smirksify.py

In this script, we start with a set of clustered molecular fragments with specified
indexed atoms as those you would use to build a ClusterGraph.
We then build cluster Graphs to create the initial SMIRKS patterns and check
that the generated SMIRKS patterns retain the typing from the input cluster.
Next we run a series of iterations removing SMIRKS decorators.
If this "move" doesn't change the way the molecules are typed then the change is accepted.

This class takes inspiration from the tool SMIRKY previously published by the
Open Force Field Initiative:
github.com/openforcefield/smarty

In theory, it is possible this process of removing decorators could be more
systematic/deterministic, however this is a first approach to see
if extracted SMIRKS patterns can do better than SMIRKY.
Also, this approach will be more general since the input clusters do not
rely on a reference force field.
"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy

from chemper.graphs.environment import ChemicalEnvironment as CE
from chemper.mol_toolkits import mol_toolkit
from chemper.chemper_utils import ImproperDict, ValenceDict, \
    get_typed_molecules, match_reference
from chemper.graphs.cluster_graph import ClusterGraph

import numpy as np


# =============================================================================================
# private subroutines
# =============================================================================================

def print_smirks(smirks_list):
    """
    Prints out the provided smirks list
    in a table like format with label | SMIRKS

    Parameters
    -----------
    smirks_list : list of tuples
        list in the form [ (label, SMIRKS), ...]
    """
    str_form = " {0:<20} | {1:} "

    print()
    print(str_form.format("Label", "SMIRKS"))
    print('='*80)
    for label, smirks in smirks_list:
        print(str_form.format(label, smirks))
        print('-'*80)
    print()


class ClusteringError(Exception):
    """
    Exception for when the SMIRKSifier is unable to create
    a list of SMIRKS to maintain the input clusters.
    """
    def __init__(self, msg):
        Exception.__init__(self, msg)
        self.msg = msg


# =============================================================================================
# SMIRKSifier
# =============================================================================================

class SMIRKSifier(object):
    """
    Generates complex SMIRKS for a given cluster of substructures
    and then reduces the decorators in those smirks
    """
    def __init__(self, molecules, cluster_list,
                 max_layers=5, verbose=True, strict_smirks=True):
        """
        Parameters
        ----------
        molecules : list of Mols
            These can be chemper Mols or molecules from any supported toolkit
            (currently OpenEye or RDKit)

        cluster_list : list of labels and smirks_atom_lists
            For each label the user should provide a list tuples for atom indices
            in each molecule you want included in that cluster.

            For example, if you wanted all atoms with indices (0,1) and (1,2) to be in cluster 'c1'
            and atoms (2,3) in cluster 'c2' for each of two molecules then cluster_list would be

            [ ('c1', [ (0,1), (1,2) ], [ (0,1), (1,2) ]),
              ('c2', [ (2,3)        ], [ (2,3)        ]) ]

            To see an example of this in action checkout
            https://github.com/MobleyLab/chemper/tree/master/examples

        max_layers : int (optional)
            default = 5
            how many atoms away from the indexed atoms should
            we consider at the maximum

        verbose : boolean (optional)
            default = True
            If true information is printed to the command line during reducing

        strict_smirks : boolean (optional)
            default = True
            If False it will not raise an error when incapable of making SMIRKS
            This setting is not recommended unless you are a master user
            or developer trying to test current behavior.
            The variable SMIRKSifier.checks will tell you if the SMIRKS
            generation failed when strict_smirks = False
        """
        self.molecules = [mol_toolkit.Mol(m) for m in molecules]
        self.intermediate_smirks = dict()
        self.cluster_list = cluster_list
        self.verbose = verbose
        self.max_layers = max_layers
        self.strict_smirks = strict_smirks

        # determine the type of SMIRKS for symmetry in indices purposes
        # This is done by making a test SMIRKS
        graph = ClusterGraph(self.molecules, cluster_list[0][1], 0)
        test_smirks = graph.as_smirks(compress=True)
        env = CE(test_smirks)
        if env.get_type() is None:
            # corresponds to an unknown chemical pattern
            self.dict_type = dict
        elif env.get_type().lower() == 'impropertorsion':
            self.dict_type = ImproperDict
        else:
            self.dict_type = ValenceDict

        # Convert input "smirks_atom_list" into a dictionary with the form:
        # {mol_idx: {(atom indices): label, ...}, ... }
        self.cluster_dict = dict()
        self.ref_labels = set()
        self.total = 0
        # form of cluster_list is [(label, [for each mol [ (tuples of atom indices)] ) ]
        for label, mol_list in self.cluster_list:
            self.ref_labels.add(label)
            # [for each mol [ (tuples of atom indices)]
            for mol_idx, atom_indice_tuples in enumerate(mol_list):
                if mol_idx not in self.cluster_dict:
                    self.cluster_dict[mol_idx] = self.dict_type()
                for atom_tuple in atom_indice_tuples:
                    self.total += 1
                    self.cluster_dict[mol_idx][atom_tuple] = label

        # make SMIRKS patterns for input clusters
        self.current_smirks, self.layers = self.make_smirks()
        if self.verbose: print_smirks(self.current_smirks)
        # check SMIRKS and save the matches to input clusters
        self.type_matches, self.checks = self.types_match_reference()

        if not self.checks:
            msg = """
                      SMIRKSifier was not able to create SMIRKS for the provided
                      clusters with %i layers. Try increasing the number of layers
                      or changing your clusters
                      """ % self.max_layers
            if self.strict_smirks:
                raise ClusteringError(msg)
            else:
                print("WARNING!", msg)

    def make_smirks(self):
        """
        Create a list of SMIRKS patterns for the input clusters.
        This includes a determining how far away from the indexed atom should
        be included in the SMIRKS (or the number of max_layers is reached)

        Returns
        -------
        smirks_list : list of tuples
            list of tuples in the form (label, smirks)
        layers : int
            number of layers actually used to specify the set clusters
        """
        layers = 0
        # try generating smirks with no layers
        smirks_list = self._make_cluster_graphs(layers)
        # store intermediate smirks and check them
        self.intermediate_smirks[layers] = smirks_list
        _, checks = self.types_match_reference(current_types=smirks_list)

        while not checks and (layers < self.max_layers):
            layers += 1
            smirks_list = self._make_cluster_graphs(layers)
            # store intermediate smirks
            self.intermediate_smirks[layers] = smirks_list
            # check current smirks patterns
            _, checks = self.types_match_reference(current_types=smirks_list)

        return smirks_list, layers

    def _make_cluster_graphs(self, layers):
        """
        Creates a list of SMIRKS using the stored
        molecules and clusters with the specified number
        of layers (atoms away from the indexed atoms)

        Parameters
        -----------
        layers : int
            number of layers (atoms away from indexed atoms) to
            include in this round of graphs

        Returns
        --------
        smirks_list : list of two tuples
            SMIRKS list in the form [ (label: SMIRKS), ...]
        """
        smirks_list = list()

        # loop through the list of fragment clusters
        for label, smirks_atom_list in self.cluster_list:
            # make a ClusterGraph for that label
            graph = ClusterGraph(self.molecules, smirks_atom_list, layers)

            # extract and save the SMIRKS for the cluster
            smirks = graph.as_smirks(compress=True)
            smirks_list.append(('zz_'+str(label), smirks))

        return smirks_list

    def types_match_reference(self, current_types=None):
        """
        Determine best match for each parameter with reference types

        Parameters
        ----------
        current_types : list of tuples with form [ (label, smirks), ]

        Returns
        -------
        type_matches : list of tuples (current_label, reference_label, counts)
            pair of current and reference labels with the number of fragments that match it
        """
        if current_types is None:
            current_types = self.current_smirks

        if len(current_types) != len(self.ref_labels):
            return set(), False

        current_assignments = get_typed_molecules(current_types, self.molecules)

        type_matches, checks_pass = match_reference(current_assignments, self.cluster_dict)
        return type_matches, checks_pass

    def reduce(self, max_its=1000, verbose=None):
        """
        Reduce the SMIRKS decorators for a number of iterations

        Parameters
        ----------
        max_its : int
            default = 1000
            The specified number of iterations
        verbose : boolean
            default = None
            will set the verboseness while running
            (if None, the current verbose variable will be used)

        Returns
        ----------
        final_smirks : list of tuples
            list of final smirks patterns after reducing in the form
            [(label, smirks)]
        """
        if not self.checks:
            print("Cannot reduce since unable to create reliable SMIRKS")
            return self.current_smirks

        red = Reducer(self.current_smirks, self.molecules, verbose)
        self.current_smirks = red.run(max_its)
        return self.current_smirks


# ==================================================================================
# Reducer class
# ==================================================================================

class Reducer():
    """
    Reducer starts with any list of SMIRKS and removes unnecessary decorators
    while maintaining typing on input molecules.
    This was created to be used as a part of the SMIRKSifier.reduce function.
    However, if you have complex SMIRKS and a list of molecules you can
    also reduce those patterns independently.

    Attributes
    ----------
    current_smirks : list of tuples
                    current SMIRKS patterns in the form (label, smirks)
    mols : list of chemper molecules
          molecules being used to reduce the input SMIRKS
    cluster_dict : dictionary
                  Dictionary specifying typing using current SMIRKS in the form:
                  {mol_idx:
                        { (tuple of atom indices): label } }
    """
    def __init__(self, smirks_list, mols, verbose=False):
        """
        Parameters
        ----------
        smirks_list : list of tuples
            set of SMIRKS patterns in the form (label, smirks)
        mols : list of molecules
            Any chemper compatible molecules accepted (ChemPer, OpenEye, or RDKit)
        """
        self.verbose = verbose
        self.current_smirks = copy.deepcopy(smirks_list)

        # make reference clusters
        ref_smirks = [('ref_%s' % l, smirks) for l, smirks in smirks_list]
        self.molecules = [mol_toolkit.Mol(m) for m in mols]
        self.cluster_dict = get_typed_molecules(ref_smirks, self.molecules)

    def remove_one_sub_dec(self, input_ors, ref_idx):
        """
        Remove one OR decorator from the specified index
        # i.e. [(#6, [X4, +0]), (#7, [X3]) ] --> [(#6, [+0]), (#7, [X3]) ]

        Parameters
        -----------
        input_ors : list of two tuples
            OR decorators in the form from ChemicalEnvironments
            that is [ (base, [decorators, ]), ... ]
        ref_idx : int
            The index from this list to use when removing one sub-decorator

        Returns
        --------
        new_ors : list of two tuples
            New OR decorators
        """
        new_ors = copy.deepcopy(input_ors)
        # remove one or decorator from a single type
        ref_or = new_ors[ref_idx]
        decs = ref_or[1]
        decs.remove(np.random.choice(decs))
        new_ors[ref_idx] = (ref_or[0], decs)
        return new_ors

    def remove_ref_sub_decs(self, input_ors, ref_idx):
        """
        Remove all of the ORdecorators at the specified index
        i.e. [(#6, [X4, +0]), (#7, [X3]) ] --> [(#6, []), (#7, [X3]) ]

        Parameters
        -----------
        input_ors : list of two tuples
            OR decorators in the form from ChemicalEnvironments
            that is [ (base, [decorators, ]), ... ]
        ref_idx : int
            The index from this list to use when removing one set of sub-decorators

        Returns
        --------
        new_ors : list of two tuples
            New OR decorators
        """
        new_ors = copy.deepcopy(input_ors)
        # remove all decorators on one OR type
        ref_or = new_ors[ref_idx]
        new_ors[ref_idx] = (ref_or[0], list())
        return new_ors

    def remove_ref(self, input_ors, ref_idx):
        """
        Remove the decorator at the referenced index
        i.e. [(#6, [X4, +0]), (#7, [X3]) ] --> [(#7, [X3])]

        Parameters
        -----------
        input_ors : list of two tuples
            OR decorators in the form from ChemicalEnvironments
            that is [ (base, [decorators, ]), ... ]
        ref_idx : int
            The OR decorators at ref_idx will be removed entirely

        Returns
        --------
        new_ors : list of two tuples
            New OR decorators
        """
        new_ors = copy.deepcopy(input_ors)
        # remove the single OR type at or_idx
        ref = new_ors[ref_idx]
        new_ors.remove(ref)
        return new_ors

    def remove_all_bases(self, input_ors):
        """
        convert all bases to [*]
        i.e. [(#6, [X4, +0]), (#7, [X3]) ] --> [(*, [X4, +0]), (*, [X3]) ]

        Parameters
        -----------
        input_ors : list of two tuples
            OR decorators in the form from ChemicalEnvironments
            that is [ (base, [decorators, ]), ... ]

        Returns
        --------
        new_ors : list of two tuples
            New OR decorators
        """
        new_all_ors = [('*', d) for b, d in input_ors]
        return new_all_ors

    def remove_all_dec_type(self, input_ors):
        """
        remove all decorators of the same type, like all 'X' decorators
        i.e. [(#6, [X4, +0]), (#7, [X3]) ] --> [(#6, [+0]), (#7, []) ]

        Parameters
        -----------
        input_ors : list of two tuples
            OR decorators in the form from ChemicalEnvironments
            that is [ (base, [decorators, ]), ... ]

        Returns
        --------
        new_ors : list of two tuples
            New OR decorators
        """
        all_decs = set([d for b,decs in input_ors for d in decs])
        remove_dec = np.random.choice(list(all_decs))

        # we start by defining a criteria for when the decorator
        # should not be removed.
        def criteria(x): return remove_dec[0] not in x

        if '+' in remove_dec or '-' in remove_dec:
            def criteria(x): return not ('+' in x or '-' in x)
        elif remove_dec in ['a', 'A']:
            def criteria(x): return x.lower() != 'a'
        elif 'r' in remove_dec:
            def criteria(x): return 'r' not in x

        new_ors = list()
        for b, decs in input_ors:
            new_decs = [d for d in decs if criteria(d)]
            new_ors.append((b, new_decs))
        return new_ors

    def remove_or_atom(self, input_all_ors, or_idx):
        """
        makes specific types of changes based on atom OR decorators

        Parameters
        ----------
        input_all_ors: list of OR decorators
            [ (base, [decs]), ...]
        or_idx: index that should be used to guide changes

        Returns
        --------
        new_ors : list of two tuples
            new or decorators
        """
        ref_or = input_all_ors[or_idx]

        # Start by checking what removal choices are available with the
        # current list of OR decorators
        # ---------------------------------------------------------------------
        # one choice is to remove all decorators (0)
        # it is always an option
        choices = ['all']
        # are there non-generic bases? If so removing all bases is an option
        if len(set([b for b,ds in input_all_ors])) > 1:
            choices.append('gen_base')
        # there are two options if the ref type has OR decorators
        if len(ref_or[1]) > 0:
            # remove all OR decorators (#6, [X4, +0]) --> (#6, [])
            choices.append('all_ref_decs')
            if len(ref_or[1]) > 1:
                # remove one OR decorator (#6, [X4, +0]) --> (#6, [+0])
                choices.append('one_dec')
        # if there are more than one or type
        if len(input_all_ors) > 1:
            # you can remove just one ORtype (the reference)
            choices.append('remove_ref')
            # check that there are decorators
            all_decs = set([d for b, decs in input_all_ors for d in decs])
            if len(all_decs) > 0:
                # you can remove 1 type of decorator, i.e. all 'Xn' decorators
                choices.append('all_one_dec')
        # ---------------------------------------------------------------------

        # Make a random choice from the available change options
        change = np.random.choice(choices)

        # Based on the option chosen call the appropriate method
        if change == 'all':
            # remove ALL OR decorators
            # i.e.[(  # 6, [X4, +0]), (#7, [X3]) ] --> []
            return list()

        if change == 'remove_ref':
            return self.remove_ref(input_all_ors, or_idx)

        if change == 'all_ref_decs':
            return self.remove_ref_sub_decs(input_all_ors, or_idx)

        if change == 'one_dec':
            return self.remove_one_sub_dec(input_all_ors, or_idx)

        if change == 'gen_base':
            return self.remove_all_bases(input_all_ors)

        # change = 'all_one_dec'
        return self.remove_all_dec_type(input_all_ors)

    def remove_or(self, input_all_ors, is_bond=False):
        """
        Changes the OR decorators by removing some of them

        Parameters
        -----------
        input_all_ors : list of tuples
            these are the OR decorators for an atom or bond
            from a ChemicalEnvironment
        is_bond : boolean
            are these decorators from from a bond (False for atom)

        Returns
        --------
        new_ors : list of two tuples
            new OR decorators
        """
        if len(input_all_ors) == 0:
            return input_all_ors, False

        # chose a set of or decorators to change
        # these come in the form [base, (decorators)]
        all_ors = copy.deepcopy(input_all_ors)
        or_idx = np.random.randint(len(all_ors))

        # atoms have more OR decorators and therefore more options for
        # how they can be removed
        if not is_bond:
            return self.remove_or_atom(all_ors, or_idx), True

        # For a bond, either one OR type is removed
        # or they are all removed.
        if np.random.rand() > 0.5:
            # remove just one or type
            return self.remove_ref(all_ors, or_idx), True
        # remove all ors
        return list(), True

    def remove_and(self, input_all_ands):
        """
        removes a decorator that is AND'd in the original SMIRKS

        Parameters
        -----------
        input_all_ands : list
            List of AND decorators

        Returns
        --------
        new_ands : list
            List of new AND decorators
        """
        # if there are no ands return with no changes
        if len(input_all_ands) == 0:
            return input_all_ands, False

        if np.random.rand() > 0.5:
            # otherwise randomly remove an AND decorator
            all_ands = copy.deepcopy(input_all_ands)
            all_ands.remove(np.random.choice(all_ands))
        else:
            all_ands = list()
        return all_ands, True

    def _get_item_and_remove_options(self, env):
        """
        This function chooses and Atom or Bond
        which will then have some of its decorators removed.
        It also determines what part of the component can be
        removed: OR decorators(0), AND decorators(1), the whole atom(2)

        Parameters
        ----------
        env: ChemicalEnvironment

        Returns
        -------
        item: an Atom or Bond from the Chemical Environment
        dec_opts: list
            0 = has OR decorators to remove
            1 = has AND decorators to remove
            2 = is an atom that could be removed
        """
        items = list()
        odds = list()
        for a_b in env.get_atoms() + env.get_bonds():
            count = len(a_b.and_types)
            for o in a_b.or_types:
                if o[0] != '*':
                    count += 1
                count += len(o[1])

            # a wild card, atom with no index should be considered ([*])
            # should also be on this list so it can be removed
            if isinstance(a_b, CE.Atom) and not isinstance(a_b, CE.Bond) and env.is_unindexed(a_b):
                count += 1

            items.append(a_b)
            odds.append(count)

        # normalize odds to weights
        odds = np.array(odds)

        # if the odds for changing all atoms and bonds is zero return None
        if odds.sum() == 0:
            return None, list()

        weights = odds / odds.sum()

        # choose an atom or bond with the probabilities:
        item = np.random.choice(items, p=weights)
        dec_opts = list()
        if len(item.or_types) > 0:
            dec_opts.append('remove_ors')
        if len(item.and_types) > 0:
            dec_opts.append('remove_ands')

        if not isinstance(item, CE.Bond): # then it is an atom
            if env.get_valence(item) == 1 and env.is_unindexed(item):
                dec_opts.append('remove_atom')

        return item, dec_opts

    def remove_decorator(self, smirks):
        """
        Chose an atom or bond in the input smirks pattern
        and then remove one decorator from it.

        Parameters
        -----------
        smirks : str
            A SMIRKS string which should be reduced

        Returns
        --------
        new_smirks : str
            A new SMIRKS pattern
        is_changed : bool
            True if some of the decorators were successfully removed
        """
        env = CE(smirks)
        sub, dec_opts = self._get_item_and_remove_options(env)

        # note this should be impossible
        if sub is None or len(dec_opts) == 0:
            return smirks, False

        change = np.random.choice(dec_opts)
        if change == 'remove_ors':
            new_or_types, changed = self.remove_or(sub.or_types, isinstance(sub, CE.Bond))
            if not changed:
                return smirks, False
            sub.or_types = new_or_types
        elif change == 'remove_ands':
            new_and_types, changed = self.remove_and(sub.and_types)
            if not changed:
                return smirks, False
            sub.and_types = new_and_types
        else: # change == 'remove_atom'
            remove = env.remove_atom(sub)
            if not remove:
                return smirks, False

        return env.as_smirks(), True

    def run(self, max_its=1000, verbose=None):
        """
        Reduce the SMIRKS decorators for a number of iterations

        Parameters
        ----------
        max_its : int
            The specified number of iterations
        verbose : boolean
            will set the verboseness while running
            (if None, the current verbose variable will be used)

        Returns
        ----------
        final_smirks : list of tuples
            list of final smirks patterns after reducing in the form
            [(label, smirks)]
        """
        temp_verbose = self.verbose
        if verbose is not None:
            self.verbose = verbose

        for it in range(max_its):
            if self.verbose: print("Iteration: ", it)

            proposed_smirks = copy.deepcopy(self.current_smirks)
            # chose a SMIRKS to change
            change_idx = np.random.randint(len(proposed_smirks))
            change_entry = proposed_smirks[change_idx]

            # generate new smirks
            new_smirks, changed = self.remove_decorator(change_entry[1])
            if self.verbose: print("Attempting to change SMIRKS #%i\n%s  -->  %s"\
                                   % (change_idx, change_entry[1], new_smirks))
            if not changed:
                if self.verbose:
                    print("Rejected!\nNo change made to SMIRKS\n%s" % change_entry[1])
                    print('-'*90)
                continue

            # update proposed list and check if it still matches the reference
            proposed_smirks[change_idx] = (change_entry[0], new_smirks)

            current_assignments = get_typed_molecules(proposed_smirks, self.molecules)
            _, proposed_checks = match_reference(current_assignments, self.cluster_dict)

            if proposed_checks:
                if self.verbose: print("Accepted! ")
                self.current_smirks = copy.deepcopy(proposed_smirks)
            else:
                if self.verbose: print("Rejected!\n proposed SMIRKS changed the way fragments are clustered")
            if self.verbose: print('-'*90)

        if self.verbose: print_smirks(self.current_smirks)
        self.verbose = temp_verbose

        return self.current_smirks




