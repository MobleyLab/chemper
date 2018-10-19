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

AUTHORS

Caitlin Bannan <bannanc@uci.edu>, UC Irvine
Additional contributions from the Mobley lab, UC Irvine,
including David Mobley, and Camila Zanette
and from the Chodera lab, John Chodera and Josh Fass

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy

import networkx as nx
import time
from chemper.graphs.environment import ChemicalEnvironment
from chemper.mol_toolkits import mol_toolkit
from chemper.chemper_utils import ImproperDict, ValenceDict, \
    get_typed_molecules, is_valid_smirks, match_reference

import numpy as np
from numpy import random

# =============================================================================================
# SMIRKSifier
# =============================================================================================

class SMIRKSifier(object):
    """
    Generates complex SMIRKS for a given cluster of substructures
    and then reduces the decorators in those smirks
    """
    def __init__(self, molecules, cluster_list,
                 layers=2, verbose=True):
        """
        Parameters
        ----------
        molecules: list of Mols
            These can be chemper Mols or molecules from any supported toolkit
            (currently OpenEye or RDKit)

        cluster_list: list of labels and smirks_atom_lists
            For each label the user should provide a list tuples for atom indices
            in each molecule you want included in that cluster.

            For example, if you wanted all atoms with indices (0,1) and (1,2) to be in cluster 'c1'
            and atoms (2,3) in cluster 'c2' for each of two molecules then cluster_list would be

            [ ('c1', [ (0,1), (1,2) ], [ (0,1), (1,2) ]),
              ('c2', [ (2,3)        ], [ (2,3)        ]) ]

            To see an example of this in action checkout
            https://github.com/MobleyLab/chemper/tree/master/examples

        layers: int (optional)
            how many atoms away from the indexed atoms should we consider
            default = 2

        verbose: boolean (optional)
            If true information is printed to the command line during reducing
            default = True
        """
        self.molecules = [mol_toolkit.Mol(m) for m in molecules]
        self.cluster_list = cluster_list
        self.layers = layers
        self.verbose = verbose

        # TODO: internally determine the number of layers!!!
        self.current_smirks = self.make_cluster_graphs()
        self.print_smirks(self.current_smirks)

        # determine the type of SMIRKS for symmetry in indices purposes
        test_smirks = self.current_smirks[0][1]
        env = ChemicalEnvironment(test_smirks)
        if env.getType().lower() == 'impropertorsion':
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

        # save matches and score
        self.type_matches, checks = self.types_match_reference()
        # TODO: figure out how to test score if less than 100% raise error? or maintain that score?
        # TODO: we will include this in the PR handling internal layer management

    def make_cluster_graphs(self):
        """
        Creates a dictionary of SMIRKS with the form
        {label: SMIRKS}
        using the stored molecules and cluster_list
        """
        from chemper.graphs.cluster_graph import ClusterGraph
        smirks_list = list()

        # loop through the list of fragment clusters
        for label, smirks_atom_list in self.cluster_list:
            # make a ClusterGraph for that label
            graph = ClusterGraph(self.molecules, smirks_atom_list, self.layers)

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

        current_assignments = get_typed_molecules(current_types, self.molecules)

        type_matches, checks_pass = match_reference(current_assignments, self.cluster_dict)
        return type_matches, checks_pass

    def print_smirks(self, smirks_list=None):
        """
        Prints out the current or provided smirks list
        in a table like format
        """
        if smirks_list is None:
            smirks_list = self.current_smirks

        str_form = " {0:<20} | {1:} "

        print()
        print(str_form.format("Label", "SMIRKS"))
        print('='*80)
        for label, smirks in smirks_list:
            print(str_form.format(label, smirks))
            print('-'*80)
        print()

    def remove_or(self, input_all_ors):
        """
        removes a decorator that is OR'd in the original SMIRKS
        """
        if len(input_all_ors) == 0:
            return input_all_ors, False

        # chose a set of or decorators to change
        # these come in the form [base, (decorators)]
        all_ors = copy.deepcopy(input_all_ors)
        or_idx = random.randint(len(all_ors))
        change_or = all_ors[or_idx]
        # temporarily remove the change OR from the list
        all_ors.remove(change_or)

        # if it has no decorators, just remove the whole ORtype
        decs = change_or[1]
        if len(decs) == 0:
            return all_ors, True

        # otherwise remove one decorator from this ORtype
        base = change_or[0]
        decs.remove(random.choice(decs))
        all_ors.append((base, decs))
        return all_ors, True

    def remove_and(self, input_all_ands):
        """
        removes a decorated that is AND'd in the original SMIRKS
        """
        # if there are no ands return with no changes
        if len(input_all_ands) == 0:
            return input_all_ands, False

        # otherwise randomly remove an AND decorator
        all_ands = copy.deepcopy(input_all_ands)
        all_ands.remove(random.choice(all_ands))
        return all_ands, True

    def _get_item_and_remove_options(self, env):
        """
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
        for a_b in env.getAtoms() + env.getBonds():
            count = len(a_b.getANDtypes())
            for o in a_b.getORtypes():
                count += len(o[1])

            # a wild card, unidexed atom ([*])
            # should also be on this list so it can be removed
            if isinstance(a_b, ChemicalEnvironment.Atom) and env.isUnindexed(a_b):
                count += 1

            items.append(a_b)
            odds.append(count)

        # normalize odds to weights
        odds = np.array(odds)
        weights = odds / np.sum(odds)

        # choose an atom or bond with the probabilities:
        item = random.choice(items, p=weights)
        dec_opts = list()
        if len(item.getORtypes()) > 0:
            dec_opts.append(0)
        if len(item.getANDtypes()) > 0:
            dec_opts.append(1)

        if isinstance(item, ChemicalEnvironment.Atom):
            if env.getValence(item) == 1 and env.isUnindexed(item):
                dec_opts.append(2)

        return item, dec_opts

    def remove_decorator(self, smirks):
        """
        Chose an atom or bond in the input smirks pattern
        and then remove one decorator from it.
        """
        env = ChemicalEnvironment(smirks)
        sub, dec_opts = self._get_item_and_remove_options(env)

        # note this should be impossible
        if len(dec_opts) == 0:
            return smirks, False

        change = random.choice(dec_opts)
        if change == 0:
            new_or_types, changed = self.remove_or(sub.getORtypes())
            if not changed:
                return smirks, False
            sub.setORtypes(new_or_types)
        elif change == 1:
            new_and_types, changed = self.remove_and(sub.getANDtypes())
            if not changed:
                return smirks, False
            sub.setANDtypes(new_and_types)
        else:
            remove = env.removeAtom(sub)
            if not remove:
                return smirks, False

        return env.asSMIRKS(), True

    def reduce(self, max_its=1000, verbose=None):
        """
        Reduce the SMIRKS decorators for a number of iterations

        Parameters
        ----------
        max_its : int
            The specified number of iterations
        verbose: boolean
            will set the verboseness while running
            (if None, the current verbose variable will be used)

        Returns
        ----------
        final_smirks: list of tuples
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
            change_idx = random.randint(len(proposed_smirks))
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
            proposed_type_matches, proposed_checks = self.types_match_reference(proposed_smirks)

            if proposed_checks:
                if self.verbose: print("Accepted! ")
                self.current_smirks = copy.deepcopy(proposed_smirks)
                self.type_matches = copy.deepcopy(proposed_type_matches)
            else:
                if self.verbose: print("Rejected!\n proposed SMIRKS changed the way fragments are clustered")
            if self.verbose: print('-'*90)

        if self.verbose: self.print_smirks()
        self.verbose = temp_verbose

        return self.current_smirks




