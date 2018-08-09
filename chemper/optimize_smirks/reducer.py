#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
reducer.py

In this script, we start with a set of clustered molecular fragments with specified
indexed atoms as those you would use to build a ClusterGraph.
We then build cluster Graphs to create the initial SMIRKS patterns and check
that the generated SMIRKS patterns retain the typing from the input cluster.
Next we run a series of iterations removing SMIRKS decorators.
If this "move" of removing a decorator doesn't change the "score" then the move
is retained.

This class takes inspiration from the tool smirky previously published by the
Open Force Field Initiative:
github.com/openforcefield/smarty

Some authors theorize this process of removing decorators could be more
systematic/deterministic, however this is a first approach to see
if extracted SMIRKS patterns can do better than smirky.
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
from chemper.optimize_smirks.environment import ChemicalEnvironment
from chemper.mol_toolkits import mol_toolkit
from chemper.chemper_utils import ImproperDict, ValenceDict, \
    get_typed_molecules, is_valid_smirks
from chemper.optimize_smirks.environment import ChemicalEnvironment

import numpy
from numpy import random


# ==============================================================================
# PRIVATE SUBROUTINES
# TODO: determine which private subroutines are necessary in the new approach
# ==============================================================================

# =============================================================================================
# SMIRKS reducer
# =============================================================================================

class Reducer(object):
    """
    Generates complex SMIRKS for a given cluster of substructures 
    and then reduces the decorators in those smirks
    """
    def __init__(self, molecules, cluster_list,
                 layers=2, verbose=True):
        """
        Parameters
        ----------
        mols: list of chemper Mols
        cluster_list: list of labels and smirks_atom_lists
            with the form [(label, smirks_atom_list)] where the
            smirks_atoms_lists is a list of list of dict
            atom indices by smirks index for each molecule
            required if a list of molecules is provided.
            This is a list of dictionaries of the form [{smirks index: atom index}]
            for each molecule provided
            # TODO: this is complicated, but I'm not going to do anything else
            # until I can show this works, then we can brain storm a
            # better way to format inputs.
        layers: int (optional)
            how many atoms away from the indexed atoms should we consider
            default = 2
        max_its: int (optional)
            maximum number of iterations
            default = 1000
        verbose: boolean (optional)
            If true information is printed to the command line during reducing
            default = True
        """
        self.molecules = molecules
        self.cluster_list = cluster_list
        self.layers = layers
        self.verbose = verbose
        #TODO: figure out how to test score if less than 100% raise error? exit?

        # TODO: figure out how to handle layers (self determine or user set)
        self.current_smirks = self.make_cluster_graphs()

        # TODO: I want to change ClusterGraph to take a "cluster_dict"
        # instead of the weird list of list of dictionaries it uses now...
        test_smirks = self.current_smirks[0][1]
        env = ChemicalEnvironment(test_smirks)
        if env.getType().lower() == 'impropertorsion':
            self.dict_type = ImproperDict
        else:
            self.dict_type = ValenceDict

        self.cluster_dict = dict()
        self.ref_labels = set()
        self.total = 0
        for label, smirks_atom_list in self.cluster_list:
            self.ref_labels.add(label)
            for mol_idx, smirks_sets in enumerate(smirks_atom_list):
                self.cluster_dict[mol_idx] = self.dict_type()
                for smirks_dict in smirks_sets:
                    sorted_keys = sorted(list(smirks_dict.keys()))
                    atom_indices = tuple([smirks_dict[k] for k in sorted_keys])
                    self.total += 1
                    self.cluster_dict[mol_idx][atom_indices] = label

        self.type_matches, self.score = self.best_match_reference(self.current_smirks)

    def make_cluster_graphs(self):
        """
        Creates a dictionary of SMIRKS with the form
        {label: SMIRKS}
        using the stored molecules and cluster_list
        """
        from chemper.graphs.cluster_graph import ClusterGraph
        smirks_list = list()
        for label, smirks_atom_list in self.cluster_list:
            print("making graph: ", len(smirks_atom_list))
            graph = ClusterGraph(self.molecules, smirks_atom_list, self.layers)
            smirks = graph.as_smirks(True)
            smirks_list.append(('zz_'+str(label), smirks))
        return smirks_list

    def best_match_reference(self, current_types=None):
        """
        Determine best match for each parameter with reference types

        Parameters
        ----------
        current_types : list of list with form [label, smirks]

        Returns
        -------
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference types, along with number of current types equivalently typed in reference molecule set.
        total_type_matches : int
            The total number of corresponding types in the reference molecule set.

        Contributor:
        * Josh Fass <josh.fass@choderalab.org> contributed this algorithm.

        """
        if current_types is None:
            current_types = self.current_smirks

        current_assignments = get_typed_molecules(current_types, self.molecules)

        # check for missing tuples in dictionaries
        for mol_idx, current_dict in current_assignments.items():
            print("Molecule ", mol_idx)
            cur_keys = set(current_dict.keys())
            ref_keys = set(self.cluster_dict[mol_idx].keys())
            # check if there are indice sets in references not in current
            if ref_keys - cur_keys:
                print("missing keys for mol ", self.molecules[mol_idx].get_smiles())
                print('reference: ', ref_keys)
                print('current: ', cur_keys)
                print('difference: ', ref_keys-cur_keys)
                #return None, None

        # Create bipartite graph (U,V,E) matching current types U with
        # reference types V via edges E with weights equal to number of types in common.
        if self.verbose: print('Creating graph matching current types with reference types...\n')
        initial_time = time.time()
        graph = nx.Graph()

        # Get current types and reference types
        cur_labels = [ lab for (lab, smirks) in current_types ]
        # ref_names = self.ref_lables
        # check that current types are not in reference types
        if set(cur_labels) & set(self.ref_labels):
            raise Exception("Current and reference type names must be unique")

        # Add current types
        for c_lab in cur_labels:
            graph.add_node(c_lab, bipartite=0)
        # add reference types
        for r_lab in self.ref_labels:
            graph.add_node(r_lab, bipartite=1)
        # Add edges.
        types_in_common = dict()
        for c_lab in cur_labels:
            for r_lab in self.ref_labels:
                types_in_common[(c_lab, r_lab)] = 0

        print(self.cluster_dict)
        for mol_idx, index_dict in self.cluster_dict.items():
            for indices, r_lab in index_dict.items():
                c_lab = current_assignments[mol_idx][indices]
                print(c_lab)
                types_in_common[(c_lab, r_lab)] += 1

        for c_lab in cur_labels:
            for r_lab in self.ref_labels:
                weight = types_in_common[(c_lab, r_lab)]
                graph.add_edge(c_lab, r_lab, weight=weight)

        elapsed_time = time.time() - initial_time
        if self.verbose:
            print('Graph creation took %.3f s\n' % elapsed_time)
            print('Computing maximum weight match...\n')

        initial_time = time.time()
        for e in graph.edges(data=True):
            print(e)
        mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)
        print('matching: ', mate)
        elapsed_time = time.time() - initial_time

        if self.verbose: print('Maximum weight match took %.3f s\n' % elapsed_time)

        # Compute match dictionary and total number of matches.
        type_matches = list()
        total_type_matches = 0
        for lab1, lab2 in mate:
            counts = graph[lab1][lab2]['weight']
            total_type_matches += counts
            if lab1 in cur_labels:
                type_matches.append(( lab1, lab2, counts))
            else:
                type_matches.append( (lab2, lab1, counts))

            # else:
            #    type_matches.append( (c_lab, None, None ))
            # figure out how to know what is missing

        # TODO: determine how/if we want to print matching types

        return type_matches, total_type_matches

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
        for label, smirks in smirks_list:
            print(str_form.format(label, smirks))
        print()

    def remove_or(self, input_all_ors):
        """
        removes a decorator that is OR'd in the original SMIRKS
        """
        if len(input_all_ors) == 0:
            return input_all_ors

        all_ors = copy.deepcopy(input_all_ors)
        change_or = random.choice(all_ors)
        all_ors.remove(change_or)

        decs = change_or[1]
        if len(decs) == 0:
            # just remove the whole ORtype
            return all_ors

        # otherwise remove one decorator from this ORtype
        base = change_or[0]
        decs.remove(random.choice(decs))
        all_ors.append((base, decs))
        return all_ors

    def remove_and(self, input_all_ands):
        """
        removes a decorated that is AND'd in the original SMIRKS
        """
        if len(input_all_ands) == 0:
            return input_all_ands
        all_ands = copy.deepcopy(input_all_ands)
        all_ands.remove(random.choice(all_ands))
        return all_ands

    def remove_decorator(self, smirks):
        """
        TODO: write docs
        """
        env = ChemicalEnvironment(smirks)
        # change atom or bond with equal probability
        # also chose type of decorator removal
        if random.random() > 0.5:
            sub = random.choice(env.getAtoms())
        else: # chose a bond
            sub = random.choice(env.getBonds())

        no_ors = True
        no_ands = True
        dec_opts = list()
        if len(sub.GetORtypes()) > 0:
            dec_opts.append(0)
            no_ors = False
        if len(sub.GetANDtypes()) > 0:
            dec_opts.append(1)
            no_ands = False

        if no_ands and no_ors:
            return smirks, False

        if random.choice(dec_opts) == 0:
            new_or_types = self.remove_or(sub.GetORtypes())
            sub.SetORtypes(new_or_types)
        else:
            new_and_types = self.remove_and(sub.GetANDtypes)
            sub.SetANDtypes(new_and_types)

        return env.asSMIRKS(), True

    def reduce_smirks(self, smirks_list):
        """
        TODO:write docs
        """
        proposed_smirks = copy.deepcopy(smirks_list)
        changed = False
        idx = 0
        while not changed and idx < 1000:
            proposed_smirks = copy.deepcopy(smirks_list)
            change_entry = random.choice(proposed_smirks)
            change_idx = proposed_smirks.index(change_entry)
            new_smirks, changed = self.remove_decorator(change_entry[1])
            proposed_smirks[change_idx] = (change_entry[0], new_smirks)
            idx += 1

        return proposed_smirks, changed

    def run(self, max_its=1000):
        """
        Run sampler for the specified number of iterations.

        Parameters
        ----------
        max_its : int
            The specified number of iterations

        Returns
        ----------
        final_smirks: list of tuples
            list of final smirks patterns after reducing in the form
            [(label, smirks)]
        """
        for it in range(max_its):
            if self.verbose: print("Iteration: ", it)

            proposed_smirks = self.reduce_smirks(copy.deepcopy(self.current_smirks))
            proposed_type_matches, proposed_score = self.best_match_reference(proposed_smirks)

            if proposed_score == self.total:
                if self.verbose: print("Accepted! ")
                self.current_smirks = copy.deepcopy(proposed_smirks)
                self.type_matches = copy.deepcopy(proposed_type_matches)
                self.score = proposed_score
            else:
                if self.verbose: print("Rejected! ")

        if self.verbose: self.print_smirks()
        return self.current_smirks




