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
Open Force Field Initative:
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
from chemper.optimize_smirks import custom_dicts

import numpy
from numpy import random


# ==============================================================================
# PRIVATE SUBROUTINES
# TODO: determine which private subroutines are necessary in the new approach
# ==============================================================================

# =============================================================================================
# SMIRKS reducer
# =============================================================================================

class Reducer(object)
    """
    Generates complex SMIRKS for a given cluster of substructures 
    and then reduces the decorators in those smirks
    """
    def __init__(self, molecules, cluster_list, layers=2,
                 max_its=1000, verbose=True):
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
        # TODO: figure out how to handle layers (self determine or user set)
        self.current_smirks, self.current_score = self.make_cluster_graphs()

        # get smirks type
        test_smirks = self.smirks[0][1]
        env = ChemicalEnvironment(test_smirks)
        if env.getType().lower() != 'impropertorsion':
            self.dict_type = custom_dicts.ImproperDict
        else:
            self.dict_type = custom_dicts.ValenceDict

        self.cluster_dict = dict()
        self.ref_labels = set()
        for label, smirks_atom_list in self.cluster_list:
            self.ref_labels.add(label)
            for mol_idx, smirks_sets in enumerate(smirks_atom_list):
                self.cluster_dict[mol_idx] = self.dict_type()
                for smirks_dict in smirks_sets:
                    atom_indices = tuple([a for s,a in smirks_dict.items()])
                    self.cluster_dict[mol_idx][atom_indices] = label

    def make_cluster_graphs(self):
        """
        Creates a dictionary of SMIRKS with the form
        {label: SMIRKS}
        using the stored molecules and cluster_list
        """
        from chemper.graphs.cluster_graph import ClusterGraph
        smirks_list = list()
        for label, smirks_atom_list in self.cluster_list:
            graph = ClusterGraph(self.molecules, smirks_atom_list, self.layers)
            smirks = graph.as_smirks(True)
            smirks_list.append(('zz_'+str(label), smirks))
        matches, score = self.best_match_reference(smirks_list)
        return smirks_list, score

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

        current_assignments = self.get_typed_molecules(current_types)

        # check for missing tuples in dictionaries
        for mol_idx, current_dict in current_assignments.items():
            cur_keys = set(current_dict.keys())
            ref_keys = set(self.reference[mol_idx].keys())
            # check if there are indice sets in references not in current
            if ref_keys - cur_keys:
                raise Exception("current and reference types do not match all molecules")

        # Create bipartite graph (U,V,E) matching current types U with
        # reference types V via edges E with weights equal to number of types in common.
        if self.verbose: print('Creating graph matching current types with reference types...\n')
        initial_time = time.time()
        graph = nx.Graph()

        # Get current types and reference types
        cur_labels = [ typename for (smirks, typename) in typelist ]
        # ref_names = self.ref_lables
        # check that current types are not in reference types
        if set(cur_labels) & set(self.ref_labels):
            raise Exception("Current and reference type names must be unique")

        # Add current types
        for clab in cur_labels:
            graph.add_node(clab, bipartite=0)
        # add reference types
        for rlab in self.ref_labels:
            graph.add_node(rlab, bipartite=1)
        # Add edges.
        types_in_common = dict()
        for clab in cur_labels:
            for rlab in self.ref_labels:
                types_in_common[(clab, rlab)] = 0

        current_typed_molecules = self.get_typed_molecules(typelist)

        for smile, indexDict in current_typed_molecules.items():
            for indices, current_typename in indexDict.items():
                reference_typename = self.reference_typed_molecules[smile][indices]
                types_in_common[(current_typename, reference_typename)] += 1

        for current_typename in current_typenames:
            for reference_typename in reference_typenames:
                weight = types_in_common[(current_typename,reference_typename)]
                graph.add_edge(current_typename, reference_typename, weight=weight)
        elapsed_time = time.time() - initial_time
        self.log.write('Graph creation took %.3f s\n' % elapsed_time)

        # Compute maximum match
        self.log.write('Computing maximum weight match...\n')
        initial_time = time.time()
        mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)
        elapsed_time = time.time() - initial_time
        self.log.write('Maximum weight match took %.3f s\n' % elapsed_time)

        # Compute match dictionary and total number of matches.
        type_matches = list()
        total_type_matches = 0
        for current_typename in current_typenames:
            if current_typename in mate:
                reference_typename = mate[current_typename]
                counts = graph[current_typename][reference_typename]['weight']
                total_type_matches += counts
                type_matches.append( (current_typename, reference_typename, counts) )
            else:
                type_matches.append( (current_typename, None, None) )

        self.log.write("PROPOSED:\n")
        self.write_type_matches(typelist, type_matches)

        return (type_matches, total_type_matches)


class FragmentSampler(object):
    """
    SMIRKS sampler for atoms, bonds, angles, torsions, and impropers.
    """
    def __init__(self, molecules, typetag, AtomORbases, AtomORdecorators,
            AtomANDdecorators, BondORbases, BondANDdecorators,
            AtomIndexOdds = None, BondIndexOdds = None,
            replacements = None,  initialtypes = None,
            SMIRFF = None, temperature = 0.1, outputFile = None):
        """
        Initialize a fragment type sampler
        For VdW, Bond, Angle, Torsion, or Improper

        Parameters
        -----------
        molecules : list of OEmol objects, required
            List of molecules used for sampling
        typetag : string, required
            Must 'Bond', 'Angle', 'Torsion', 'Improper', or 'VdW'
            'VdW' is for single labeled atom

        The following parameters come in the form of tuples of two
        with ( [features], [odds or probabilities])
        if [odds or probabilities] = None all are treated equally
        odds like [1,4] will be converted to probabilities [.2, .8]
        ------------------------------------------------
        AtomORbases : list of strings and their probabilities, required
            each element can be an atomic number ('#1') or
            shorthand name that is in the replacements list
        AtomORdecorators: list of strings and their probabilities, required
            List of decorators that can be combined directly with an atom
            for example: for [#6X4, #8X2] 'X4' and 'X2' are ORdecorators
        AtomANDdecorators: list of strings and their probabilities, required
            List of decorators that are AND'd to the end of an atom
            for example: in [#6,#7,#8;H0;+0] 'H0' and '+0' are ANDdecorators
        BondORbases : list of strings and their probabilities, required
            each bond type and the probability of using it
        BondANDdecorators : of strings and their probabilitis, required
            bond decorator options such as [ ('@', 1), ('@', 1) ]
        AtomIndexOdds and BondIndexOdds have list of atom and bond indices
            options for property lists: atom indices (integer),
            or string descriptor from : indexed, unindexed, alpha, beta
            if None all atoms and bonds will be treated equally
        ------------------------------------------------
        replacements: list of the form [short hand, smarts], optional
        initialtypes: initial typelist in form [smirks, typename], optional
            if None, the typetag is used to make an empty environment, such as [*:1]~[*:2] for a bond
        SMIRFF: string, optional
            file with the SMIRFF you wish to compare fragment typing with
        temperature : float, optional, default=0.1
            Temperature for Monte Carlo acceptance/rejection
        outputFile: string, optional, default = typetag_temperature
            output name base for log files and trajectory files

        Notes
        -----

        """
        # Save properties that remain unchanged
        self.AtomORdecorators = _OddsToWeights(AtomORdecorators)
        self.AtomANDdecorators = _OddsToWeights(AtomANDdecorators)
        self.BondORbases = _OddsToWeights(BondORbases)
        self.BondORdecorators = ( [''], [1.])
        self.BondANDdecorators = _OddsToWeights(BondANDdecorators)
        self.temperature = temperature
        self.SMIRFF = SMIRFF
        self.replacements = replacements
        self.types_with_no_matches = []

        if outputFile is None:
            self.output = "%s_%.2e" % (typetag, temperature)
        else:
            self.output = outputFile

        self.log = open("%s.log" % self.output, 'w')

        if AtomIndexOdds is None:
            AtomIndexOdds = [['all'], [1]]
        if BondIndexOdds is None:
            BondIndexOdds = [ ['all'], [1]]

        self.AtomIndexOdds = _OddsToWeights(AtomIndexOdds)
        self.BondIndexOdds = _OddsToWeights(BondIndexOdds)

        if not typetag.lower() in ['bond', 'angle', 'torsion','improper','vdw']:
            raise Exception("Error typetag %s is not recognized, please use 'Bond', 'Angle', 'Torsion', 'Improper', or 'VdW' ")
        self.typetag = typetag

        # get information based on typetag
        self.forcetype, self.edges, self.sym_odds = self.get_type_info(self.typetag)

        # get molecules and add explicit hydrogens
        self.molecules = list()
        smiles = set()

        # loop through input molecules, remove repeats
        for mol in molecules:
            smile = OECreateIsoSmiString(mol)
            if not smile in smiles:
                self.molecules.append(OEMol(mol))
                smiles.add(smile)

        # if no initialtypes specified start with empty type
        self.emptyEnv = self._makeEnvironments(self.typetag, None)[0]

        # Compute total types being sampled
        empty_typelist = [[self.emptyEnv.asSMIRKS(), 'empty']]
        [empty_counts, empty_molecule_counts] = self.compute_type_statistics(empty_typelist)
        self.total_types = empty_counts['empty']
        self.IndexDict = self.get_typed_molecules(empty_typelist)

        temp_envList = self._makeEnvironments(self.typetag, initialtypes)
        self.envList = list()
        # Make typelist to fit method set up
        typelist = [[env.asSMIRKS(), env.label] for env in temp_envList]

        # check that current smirks match all types in self.molecules
        if not self.check_typed_molecules(typelist):
            raise Exception("Initial types do not type all %s in the molecules" % self.typetag)

        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        self.write_type_statistics(typelist, typecounts, molecule_typecounts)
        # Only keep fragments that are initially populated
        for env in temp_envList:
            if typecounts[env.label] > 0:
                self.envList.append(copy.deepcopy(env))
        self.baseTypes = copy.deepcopy(self.envList)

        # This might be better as a graph with unlabeled nodes where the nodes are environment objects?
        self.parents = dict()
        for env in self.envList:
            self.parents[env.label] = dict()
            self.parents[env.label]['children'] = list()
            self.parents[env.label]['SMIRKS'] = env.asSMIRKS()
            self.parents[env.label]['parent'] = None

        # Store elements without the [ ]
        self.AtomORbases = ( [], [])
        elementList = _OddsToWeights(AtomORbases)
        elementList = zip(elementList[0], elementList[1])
        for (element, prob) in elementList:
            e = element.replace('[','')
            e = e.replace(']','')
            # Check if element is used in any molecule
            for mol in self.molecules:
                useE = False
                matches = self.get_SMIRKS_matches(mol, '[%s:1]' % e)
                if len(matches) > 0:
                    useE = True
                    break
            # Save the used ones only
            if useE:
                self.AtomORbases[0].append(e)
                self.AtomORbases[1].append(prob)
            else:
                self.log.write("removing unused element (%s) from list\n" % element)

        # Store reference molecules
        self.reference_typed_molecules = dict()
        self.type_matches = None
        self.total_type_matches = None
        self.reference_typename_dict = dict()
        if self.SMIRFF is not None:
            self.log.write("Creating labeler from %s...\n" % self.SMIRFF)
            # get labeler for specified SMIRFF
            self.labeler = forcefield.ForceField(self.SMIRFF)
            labels = self.labeler.labelMolecules(self.molecules, verbose=False)

            # save the type we are considering
            self.ref_labels = [l[self.forcetype] for l in labels]

            # get smiles to key reference typed molecule dictionary
            smiles = [OECreateIsoSmiString(mol) for mol in molecules]
            # Extract list of reference SMIRKS types present in molecules
            for idx, label_set in enumerate(self.ref_labels):
                smile = smiles[idx]
                self.reference_typed_molecules[smile] = {}
                for (indices, pid, smirks) in label_set:
                    self.reference_typename_dict[pid] = smirks
                    self.reference_typed_molecules[smile][tuple(indices)] = pid

            # Check that all indices are typed for the reference
            for smile, indicesList in self.IndexDict.items():
                for indices in indicesList:
                    if indices not in self.reference_typed_molecules[smile].keys():
                        raise Exception("Reference types in SMIRFF (%s) do not \
                                type all %ss in the molecules" % (self.SMIRFF, self.typetag))

            # Compute current type matches
            [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
            # Count  types.
            self.reference_type_counts = { pid : 0 for pid, smirks in \
                    self.reference_typename_dict.items() }
            for label_set in self.ref_labels:
                for (atom_indices, pid, smirks) in label_set:
                    self.reference_type_counts[pid] += 1
            self.write_type_statistics(typelist, typecounts, molecule_typecounts, self.type_matches)
        return

    def _makeEnvironments(self, typetag, smirksList):
        """
        Given a typetag and list of SMIRKS strings
        returns a list of chemical environment objects
        of the correct type
        """
        if typetag.lower() == 'vdw':
            chemEnv = AtomChemicalEnvironment
        elif typetag.lower() == 'bond':
            chemEnv = BondChemicalEnvironment
        elif typetag.lower() == 'angle':
            chemEnv = AngleChemicalEnvironment
        elif typetag.lower() == 'torsion':
            chemEnv = TorsionChemicalEnvironment
        elif typetag.lower() == 'improper':
            chemEnv = ImproperChemicalEnvironment
        else:
            return None

        if smirksList is None:
            return [chemEnv(label=0, replacements=self.replacements)]

        envList = list()
        for smirks, typename in smirksList:
            envList.append(chemEnv(smirks, typename, self.replacements))

        return envList

    def get_type_info(self, typetag):
        """
        Uses typetag to get the Force type key word used to read the SMIRFF file

        Parameters
        ----------
        typetag: string, required
            'vdw', 'bond', 'angle', 'torsion', 'improper'
            indicates the type of system being sampled

        Returns
        ------------
        force - term used in forcefield objects
        edges - edges of fragments, for example [1,3] for angle
        sym_odds - odds of the edges being identical, based on smirff99Frosst version 1.0.1
        """
        force = None
        edges = []
        sym_odds = ([True, False], [1,1])
        if typetag.lower() == 'vdw':
            force = 'NonbondedGenerator'
        if typetag.lower() == 'bond':
            force = 'HarmonicBondGenerator'
            edges = [1,2]
            sym_odds = ([True, False], [1,5])
        if typetag.lower() == 'angle':
            force = 'HarmonicAngleGenerator'
            edges = [1,3]
            sym_odds = ([True, False], [3,2])
        if typetag.lower() == 'torsion' or typetag.lower() == 'improper':
            force = 'PeriodicTorsionGenerator'
            if typetag.lower() == 'torsion':
                edges = [1,4]
                sym_odds = ([True, False], [15, 100])
        return force, edges, sym_odds

    def get_SMIRKS_matches(self, mol, smirks):
        """
        Gets atom indices for a smirks string in a given molecule

        Parameters
        ----------
        mol : an OpenEye OEMol object
        smirks : a string for the SMIRKS string being parsed

        Returns
        --------
        matches: list of tuples
            atom indices for labeled atom in the smirks
        """
        if self.replacements is not None:
            smirks = OESmartsLexReplace(smirks, self.replacements)

        qmol = OEQMol()
        if not OEParseSmarts(qmol, smirks):
            raise Exception("Error parsing SMIRKS %s" % smirks)

        # Impropers have different symmetry from other fragments
        if self.typetag.lower() == 'improper':
            matches = dict()
        else:
            # ValenceDict allow for symmetric fragments
            # for example bond (1,2) is identical to (2,1)
            matches = forcefield.ValenceDict()

        # then require non-unique matches
        unique = False
        ss = OESubSearch(qmol)

        for match in ss.Match(mol, unique):
            indices = dict()
            for ma in match.GetAtoms():
                patMap = ma.pattern.GetMapIdx()
                # if patMap == 0, then it's an unidexed atom
                if patMap != 0:
                    indices[patMap-1] = ma.target.GetIdx()

            indices = [indices[idx] for idx in range(len(indices))]
            matches[indices] = ''

        return matches.keys()

    def get_typed_molecules(self, typelist):
        """
        Creates a dictionary assigning a typename
        for each set of atom indices in each molecule

        Parameters
        ----------
        typelist: list of tuples in the form (smirks, typename)

        Returns
        -------
        typeDict: embedded dictionary
            keys: SMILES string for each molecule
                keys: tuple of indices assigned a parameter type
        """
        typeDict = dict()
        for mol in self.molecules:
            smiles = OECreateIsoSmiString(mol)
            typeDict[smiles] = {}
            for [smirks, typename] in typelist:
                matches = self.get_SMIRKS_matches(mol, smirks)
                for match in matches:
                    typeDict[smiles][match] = typename

        return typeDict

    def check_SMIRK_matches(self, smirks):
        """
        checks that a give smirks matches at least one set of indices in the molecules
        """
        for mol in self.molecules:
            matches = self.get_SMIRKS_matches(mol, smirks)
            if len(matches) > 0:
                return True

        return False

    def check_typed_molecules(self, typelist):
        """
        given a typelist of [smirks, typename] it check that
        all types in each molecule can be typed
        """
        typed_molecules = self.get_typed_molecules(typelist)
        for smile, indicesList in self.IndexDict.items():
            for indices in indicesList:
                if indices not in typed_molecules[smile].keys():
                    return False
        return True

    def best_match_reference_types(self, typelist):
        """
        Determine best match for each parameter with reference types

        Parameters
        ----------
        typelist : list of list with form [smarts, typename]

        Returns
        -------
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference types, along with number of current types equivalently typed in reference molecule set.
        total_type_matches : int
            The total number of corresponding types in the reference molecule set.

        Contributor:
        * Josh Fass <josh.fass@choderalab.org> contributed this algorithm.

        """
        if self.SMIRFF is None:
            self.log.write('No reference SMIRFF specified, so skipping likelihood calculation.\n')
            return None, None # None for type_matches and total_type_matches

        # Create bipartite graph (U,V,E) matching current types U with reference types V via edges E with weights equal to number of types in common.
        self.log.write('Creating graph matching current types with reference types...\n')
        initial_time = time.time()
        graph = nx.Graph()

        # Get current types and reference types
        current_typenames = [ typename for (smirks, typename) in typelist ]
        reference_typenames = [ typename for typename, smirks in \
                self.reference_typename_dict.items()]
        # check that current types are not in reference types
        if set(current_typenames) & set(reference_typenames):
            raise Exception("Current and reference type names must be unique")
        # Add current types
        for typename in current_typenames:
            graph.add_node(typename, bipartite=0)
        # add reference types
        for typename in reference_typenames:
            graph.add_node(typename, bipartite=1)
        # Add edges.
        types_in_common = dict()
        for current_typename in current_typenames:
            for reference_typename in reference_typenames:
                types_in_common[(current_typename,reference_typename)] = 0

        current_typed_molecules = self.get_typed_molecules(typelist)

        for smile, indexDict in current_typed_molecules.items():
            for indices, current_typename in indexDict.items():
                reference_typename = self.reference_typed_molecules[smile][indices]
                types_in_common[(current_typename, reference_typename)] += 1

        for current_typename in current_typenames:
            for reference_typename in reference_typenames:
                weight = types_in_common[(current_typename,reference_typename)]
                graph.add_edge(current_typename, reference_typename, weight=weight)
        elapsed_time = time.time() - initial_time
        self.log.write('Graph creation took %.3f s\n' % elapsed_time)

        # Compute maximum match
        self.log.write('Computing maximum weight match...\n')
        initial_time = time.time()
        mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)
        elapsed_time = time.time() - initial_time
        self.log.write('Maximum weight match took %.3f s\n' % elapsed_time)

        # Compute match dictionary and total number of matches.
        type_matches = list()
        total_type_matches = 0
        for current_typename in current_typenames:
            if current_typename in mate:
                reference_typename = mate[current_typename]
                counts = graph[current_typename][reference_typename]['weight']
                total_type_matches += counts
                type_matches.append( (current_typename, reference_typename, counts) )
            else:
                type_matches.append( (current_typename, None, None) )

        self.log.write("PROPOSED:\n")
        self.write_type_matches(typelist, type_matches)

        return (type_matches, total_type_matches)

    def write_type_matches(self, typelist, type_matches):
        """
        Show pairing of current to reference types.

        Parameters
        ----------
        typelist: list of [smirks, typenames]
        type_matches : list of (current_typename, reference_typename, counts)
            List of type matches.

        Returns
        --------
        fraction_matched_types, the fractional count of matched types

        """
        self.log.write('%s type matches:\n' % self.typetag)
        if type_matches is None:
            return 0
        total_type_matches = 0
        current_dict = dict()
        for (current_typename, reference_typename, counts) in type_matches:
            current_dict[current_typename] = (reference_typename, counts)

        for [current_smirks, current_typename] in typelist:
            (reference_typename, counts) = current_dict[current_typename]
            current_combo = "%s: %s" % (current_typename, current_smirks)
            if reference_typename is not None:
                reference_smirks = self.reference_typename_dict[reference_typename]
                reference_combo = "%s: %s" % (reference_typename, reference_smirks)
                self.log.write('%-64s matches %64s: %8d %-10s types matched\n' % (current_combo, reference_combo, counts, self.typetag))
                total_type_matches += counts
            else:
                self.log.write('%-64s no match\n' % (current_combo))

        fraction_matched = float(total_type_matches) / float(self.total_types)
        self.log.write('%d / %d total %ss match (%.3f %%)\n' % (total_type_matches, self.total_types, self.typetag, fraction_matched * 100))

        return fraction_matched

    def pick_an_atom(self, env):
        """
        Uses AtomIndexOdds to chose an atom
        returns an Atom object and its associated probability
        """
        choices = copy.deepcopy(self.AtomIndexOdds)
        atom = None
        while atom is None:
            descriptor,prob = _PickFromWeightedChoices(choices)
            # TODO: selectAtom has random built in, account for those odds?
            atom = env.selectAtom(descriptor)
            choices[1][choices[0].index(descriptor)] = 0
        return atom, prob

    def pick_a_bond(self, env):
        """
        Uses BondIndexOdds to chose a bond
        returns a Bond object and its associated probability
        """
        choices = copy.deepcopy(self.BondIndexOdds)
        bond = None
        while bond is None:
            descriptor, prob = _PickFromWeightedChoices(choices)
            # TODO: selectBond also uses random, not accounting for those odds
            bond = env.selectBond(descriptor)
            choices[1][choices[0].index(descriptor)] = 0
        return bond, prob

    def add_swap_delete(self, current, new, new_prob, probabilities = None):
        """
        Makes a change to a current list of properties

        Parameters
        -----------
        current: current list for that property
        new: new property for the list
        new_prob: probability for picking that new property
        probabilities: probability or Odds list for
            [add property, swap property, delete a property]

        Returns
        --------
        newList, probability
            of creating it
            returns None, 0 if no change was made
        """
        opts = [1,2,3]
        if probabilities is None:
            probabilities = [1,1,1]

        if len(current) == 0:
            # if nothing in the list then adding is the only option
            move = 1
            move_prob = 1
        elif new in current: # add and swap are not an option
            move = 3
            move_prob = 1
        else:
            move, move_prob = _PickFromWeightedChoices( (opts, probabilities))

        newList = copy.deepcopy(current)

        if move == 2 or move == 3: # swap or delete
            remove, remove_prob = _PickFromWeightedChoices((current, None))
            newList.remove(remove)
            move_prob *= remove_prob

        if move == 1 or move == 2: # add or swap property
            newList.append(new)
            move_prob *= new_prob

        return newList, move_prob

    def add_atom(self, env, atom):
        """
        Adds an atom bonded to the current atom
        New atoms have 1 ORbase and no ORdecorators or ANDdecorators

        Returns probability of creating this atom
        """
        #choose ORbase
        atom_or_options = _RemoveBlankOdds(copy.deepcopy(self.AtomORbases))
        base, prob = _PickFromWeightedChoices( atom_or_options )
        new_atom = env.addAtom(atom, newORtypes = [(base, [])])
        return prob

    def isremoveable(self,env, atom):
        """
        Returns true if the atom can be removed
        """
        if env.isIndexed(atom):
            return False # indexed atom isn't removeable
        if len(atom.getANDtypes() ) > 0:
            return False # AND decorators not added to new atoms
        ORs = atom.getORtypes()
        if len(ORs) > 1: # only 1 OR type on new atoms
            return False
        if len(ORs) == 1 and len(ORs[0][1]) > 0:
            return False # only 1 OR type with no OR decorators

        bonds = env.getBonds(atom)
        if len(bonds) > 1:
            return False # atom has more than 1 neighbor
        bond = bonds[0]
        if len(bond.getANDtypes() ) > 0:
            return False
        if len(bond.getORtypes() ) > 0:
            return False

        return True

    def _is_symmetric(self, env):
        """
        Returns true if environment is currently symmetric,
        meaning changes can be made to the outer parts of the
        fragment identically
        """
        # if there are unindexed atoms we will assume the environment is no longer symmetric
        if len(env.getUnindexedAtoms()) > 0:
            return False

        # check current atom edges
        atom1 = env.selectAtom(self.edges[0])
        atom2 = env.selectAtom(self.edges[1])
        if atom1.asSMARTS() != atom2.asSMARTS():
            return False

        # if looking at Bonds stop here
        if self.typetag.lower() == 'bond':
            return True

        # otherwise, check bonds are also symmetric
        bond1 = env.selectBond(self.edges[0])
        bond2 = env.selectBond(self.edges[1] - 1)
        if bond1.asSMARTS() != bond2.asSMARTS():
            return False

        # passes all checks
        return True

    def change_atom(self, env, atom):
        """
        Makes changes to the provided Atom object
            in the provided environment
        returns the probability of making this change
            if 0 no change was made
        """
        # Check if symmetric
        symmetric_change = False
        sym_change_odds = 1
        if atom.index in self.edges:
            if self._is_symmetric(env):
                symmetric_change, sym_change_odds = _PickFromWeightedChoices(self.sym_odds)

        # assign options and probabilities:
        # start with OR base change opt = 3
        opts = [3]
        probs = [5]
        # In order to have add_atom as an allowed option,
        # The current atom must have at least one decorator and are not Beta
        decorated = (len(atom.getORtypes()) > 0 or len(atom.getANDtypes()) > 0)
        if (not env.isBeta(atom)) and decorated:
            opts.append(1)
            probs.append(1)

        # check if removeable
        if self.isremoveable(env, atom):
            opts.append(0)
            probs.append(1)

        # AND decorators are optional, if they were provided check the only 1 isn't ''
        if not (len(self.AtomANDdecorators[0]) <=1 and self.AtomANDdecorators[0][0] ==''):
            opts.append(2)
            probs.append(5)
        # check for existing OR types, for OR decorator changes
        if len(atom.getORtypes()) > 0:
            # Check that or decorators were provided
            if not (len(self.AtomORdecorators[0]) <= 1 and self.AtomORdecorators[0][0] == ''):
                opts.append(4)
                probs.append(5)

        # Choose a move
        move, move_prob = _PickFromWeightedChoices( (opts, probs))
        if move == 0: # remove Atom
            removed = env.removeAtom(atom, False)
            change_prob = 1.
            return move_prob * change_prob
        elif move == 1: # add Atom
            change_prob = self.add_atom(env, atom)
            return move_prob * change_prob
        elif move == 2: # Change ANDdecorators
            change_prob = self.change_ANDdecorators(atom, self.AtomANDdecorators)
        elif move == 3: # Change ORbases
            change_prob = self.change_ORbase(atom, self.AtomORbases, self.AtomORdecorators)
        else: # move == 4 Change ORdecorators
            change_prob = self.change_ORdecorator(atom, self.AtomORdecorators)

        if symmetric_change: # make changes to mirrored atom
            temp_edges = copy.deepcopy(self.edges)
            temp_edges.remove(atom.index)
            atom2_index = temp_edges[0]
            atom2 = env.selectAtom(atom2_index)
            atom2.setORtypes(atom.getORtypes())
            atom2.setANDtypes(atom.getANDtypes())

        return move_prob * change_prob * sym_change_odds

    def change_ORdecorator(self, component, input_decorators):
        """
        Makes changes to the decorators associated with 1 ORbase for
        a given component (just atoms in this case)
        returns the probability of making this change
        """
        # Remove '' from choices, if we're changing
        # OR decorator we don't want that to count
        decorators = _RemoveBlankOdds(copy.deepcopy(input_decorators))
        # pick new decorator with the probability
        new_decor, decor_prob = _PickFromWeightedChoices(decorators)
        currentORs = component.getORtypes()
        if len(currentORs) > 0:
            # Pick an ORtype (base, decorators) to make changes to
            change_OR, base_prob = _PickFromWeightedChoices((currentORs, None))

            # Remove from currentORs to make changes
            currentORs.remove(change_OR)

        # get new decorators and add back into the ORtypes list
        new_decs, new_prob = self.add_swap_delete(change_OR[1], new_decor, decor_prob*base_prob, None)

        currentORs.append( (change_OR[0], new_decs))
        component.setORtypes(currentORs)
        return new_prob

    def change_ORbase(self, component, bases, decorators):
        """
        Changes a component (atom or bond)'s ORtype list from
            a given set of bases and decorators
        Returns probability of making the move
        """
        new_base, base_prob = _PickFromWeightedChoices(bases)
        new_decor, decor_prob = _PickFromWeightedChoices(decorators)
        current = component.getORtypes()
        new_OR = (new_base, [new_decor])
        new_list, new_prob = self.add_swap_delete(current, new_OR, base_prob*decor_prob, None)
        component.setORtypes(new_list)
        return new_prob

    def change_ANDdecorators(self, component, decorators):
        """
        Changs a component's (atom or bond) ANDtype list from
            a given set of decorators
        Returns probability of making the change
        """
        new_decor, decor_prob = _PickFromWeightedChoices(decorators)
        current = component.getANDtypes()
        new_list, new_prob = self.add_swap_delete(current, new_decor, decor_prob, None)
        component.setANDtypes(new_list)
        return new_prob

    def change_bond(self,env,bond):
        """
        Makes changes to the Bond object
        returns probability of making change
        """
        #TODO: if symmetry in change_atom works for AlkEthOH, figure out how to handle bonds
        # Can only make changes to the bond OR or AND types
        changeOR = random.choice([True, False], p = [0.7, 0.3])
        if changeOR:
            # Bonds only have OR bases (no ORdecorators)
            new_prob = self.change_ORbase(bond, self.BondORbases, self.BondORdecorators)
            return 0.7 * new_prob

        else: # change AND type
            new_prob = self.change_ANDdecorators(bond, self.BondANDdecorators)
            return new_prob * 0.3

    def create_new_environment(self, env):
        """
        Given a parent environment type it creates a new child environment
        returns child environment type and probability of creating it
        """
        new_env = copy.deepcopy(env)
        new_env.label = _get_new_label([e.label for e in self.envList])

        if len(env.getBonds()) == 0:
            choices = ([True], None)
        else:
            choices = ([True, False], [0.9, 0.1])

        # pick to make changes to an atom or bond
        change_atom, choice_prob = _PickFromWeightedChoices(choices)

        if change_atom:
            atom, atom_prob = self.pick_an_atom(new_env)
            change_prob = self.change_atom(new_env, atom)
            prob = choice_prob * atom_prob * change_prob

        else: # Change Bond
            bond, bond_prob = self.pick_a_bond(new_env)
            change_prob = self.change_bond(new_env, bond)
            prob = choice_prob * bond_prob * change_prob

        return new_env, prob

    def sample_types(self):
        """
        Perform one step of sampling the current set of chemical environments

        """
        # Copy current sets for proposal.
        proposed_envList = copy.deepcopy(self.envList)
        proposed_parents = copy.deepcopy(self.parents)
        ntypes = len(proposed_envList)

        # chose an environment from the list to focus on:
        env = random.choice(proposed_envList)

        # determine create or destroy:
        if random.random() < 0.2:
            # TODO: determine how frequently to destroy entire types
            self.log.write("Attempting to destroy type %s : %s...\n" % (env.label, env.asSMIRKS()))

            # Reject deletion of (populated) base types as we want to retain
            # generics even if empty
            if env.label in [e.label for e in self.baseTypes]:
                self.log.write("Destruction rejected for type %s because this is a generic type which was initially populated.\n" % env.label)
                return False

            # Delete the type.
            proposed_envList.remove(env)

            # Try to type all molecules.
            typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]
            if not self.check_typed_molecules(typelist):
                self.log.write("Typing failed; rejecting.\n")
                return False

            # save parent and children for this type
            parent = proposed_parents[env.label]['parent']
            children = proposed_parents[env.label]['children']
            # update parent dicitonary
            proposed_parents[parent]['children'] += children
            proposed_parents[parent]['children'].remove(env.label)
            for c in children:
                proposed_parents[c]['parent'] = parent
            del proposed_parents[env.label]

        else: # create new type from chosen environment
            new_env, prob = self.create_new_environment(env)

            self.log.write("Attempting to create new subtype: '%s' (%s) from parent type '%s' (%s)\n" % (new_env.label, new_env.asSMIRKS(), env.label, env.asSMIRKS()))
            self.log.write("\tProbability of making this environment is %.3f %%" % prob)

            # Check the SMIRKS for new_env is valid
            qmol = OEQMol()
            smirks = new_env.asSMIRKS()
            if self.replacements is not None:
                smirks = OESmartsLexReplace(smirks, self.replacements)
            # check the SMIRKS pattern is valid
            if not OEParseSmarts(qmol, smirks):
                self.log.write("Type '%s' (%s) is invalid; rejecting.\n" % (new_env.label, new_env.asSMIRKS()))
                return False

            # Check if new_env is already in types with no matches
            if new_env.asSMIRKS() in self.types_with_no_matches:
                self.log.write("Type '%s' (%s) unused in dataset; rejecting.\n" % (new_env.label, new_env.asSMIRKS()))
                return False

            # Check if proposed type is already in set.
            if new_env.asSMIRKS() in [e.asSMIRKS() for e in self.envList]:
                self.log.write("Type '%s' (%s) already exists; rejecting to avoid duplication.\n" % (new_env.label, new_env.asSMIRKS()))
                return False

            # add new type to proposed list
            proposed_envList.append(new_env)
            proposed_typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]
            # Compute updated statistics
            [proposed_typecounts, proposed_molecule_typecounts] = self.compute_type_statistics(proposed_typelist)

            # Reject if new type matches nothing
            if proposed_typecounts[new_env.label] == 0:
                self.log.write("Type '%s' (%s) unused in dataset; rejecting.\n" % (new_env.label, new_env.asSMIRKS()))
                self.types_with_no_matches.append(new_env.asSMIRKS())
                return False

            # Reject if any type is emptied (UNLESS it is a basetype
            base_labels = [e.label for e in self.baseTypes]
            for label, count in proposed_typecounts.items():
                if not label in base_labels:
                    if count == 0:
                        self.log.write("Fragment with typename %s now has no matches, rejecting.\n" % label)
                        return False

            # updated proposed parent dictionary
            proposed_parents[env.label]['children'].append(new_env.label)
            proposed_parents[new_env.label] = {}
            proposed_parents[new_env.label]['children'] = list()
            proposed_parents[new_env.label]['parent'] = env.label
            proposed_parents[new_env.label]['SMIRKS'] = new_env.asSMIRKS()

        self.log.write('Proposal is valid...\n')

        # Accept automatically if no set was provided
        if self.SMIRFF is None:
            self.envList = proposed_envList
            self.parents = proposed_parents
            return True

        # Determine number of matches
        accept = False
        # Compute change in score
        typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]
        (proposed_type_matches, proposed_total_type_matches) = self.best_match_reference_types(typelist)
        score_dif = (proposed_total_type_matches - self.total_type_matches)

        # If temeperature is 0.0 only accept improved scores
        if self.temperature == 0.0:
            self.log.write('Proposal score: %d >> %d \n' % (self.total_type_matches, proposed_total_type_matches))
            accept = score_dif > 0.0

        # Finite temperature compute likelihood function
        else:
            # Compute effective temperature and likelihood
            effective_temperature = (self.total_types * self.temperature)
            log_P_accept = score_dif / effective_temperature
            self.log.write('Proposal score: %d >> %d : log_P_accept = %.5e\n' % (self.total_type_matches, proposed_total_type_matches, log_P_accept))
            accept = (log_P_accept > 0.0) or (numpy.random.uniform() < numpy.exp(log_P_accept))

        if accept:
            # Change accepted
            self.envList = proposed_envList
            self.parents = proposed_parents
            self.type_matches = proposed_type_matches
            self.total_type_matches = proposed_total_type_matches
            return True

        else: # Change not accpted
            return False

    def compute_type_statistics(self, typelist):
        """
        Compute statistics for numnber of molecules assigned each type.

        Parameters
        ----------
        typelist: list of lists with the form [smarts, typename]

        Returns
        -------
        typecounts (dict) - number of matches for each fragment type
        molecule_typecounts (dict) - number of molecules that contain each fragment type

        """
        # Zero type counts by typename and molecule.
        typecounts = dict()
        molecule_typecounts = dict()
        for [smarts, typename] in typelist:
            typecounts[typename] = 0
            molecule_typecounts[typename] = 0

        typed_molecules = self.get_typed_molecules(typelist)
        # Count number of indice sets with each type.
        for molecule, indexDict in typed_molecules.items():
            typenames_in_molecule = set()
            for indices, typename in indexDict.items():
                typecounts[typename] += 1
                typenames_in_molecule.add(typename)
            for typename in typenames_in_molecule:
                molecule_typecounts[typename] += 1

        return (typecounts, molecule_typecounts)

    def write_type_statistics(self, typelist, typecounts, molecule_typecounts, type_matches=None):
        """
        Print type statistics.

        Parameters
        -----------
        typelist - list of lists with form [smirks, typename]
        typecounts (dict) - number of matches for each fragment type
        molecule_typecounds (dict) - number of molecules that contain each fragment type
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference types, along with number of fragment types equivalently typed in reference molecule set.
        """
        index = 1
        ntypes = 0

        if type_matches is not None:
            reference_type_info = dict()
            for (current_typename, reference_typename, count) in type_matches:
                reference_type_info[current_typename] = (reference_typename, count)

        # Print header
        if type_matches is not None:
            self.log.write("%5s   %10sS %10s   %-50s %-50s %30s\n" % ('INDEX', self.typetag.upper(), 'MOLECULES', 'TYPE NAME: SMIRKS', 'REF TYPE: SMIRKS', 'FRACTION OF REF TYPED MOLECULES MATCHED'))
        else:
            self.log.write("%5s   %10sS %10s   %-50s\n" % ('INDEX', self.typetag.upper(), 'MOLECULES', 'TYPE NAME: SMIRKS'))

        # Print counts
        for [smarts, typename] in typelist:
            current_combo = "%s: %s" % (typename, smarts)
            if type_matches is not None:
                (reference_typename, reference_count) = reference_type_info[typename]
                if reference_typename is not None:
                    reference_total = self.reference_type_counts[reference_typename]
                    reference_fraction = float(reference_count) / float(reference_total)
                    reference_combo = "%s: %s" % (reference_typename, self.reference_typename_dict[reference_typename])
                    self.log.write("%5d : %10d %10d | %-50s %-50s %7d / %7d (%7.3f%%)\n" % (index, typecounts[typename], molecule_typecounts[typename], current_combo, reference_combo, reference_count, reference_total, reference_fraction*100))
                else:
                    self.log.write("%5d : %10d %10d | %-50s\n" % (index, typecounts[typename], molecule_typecounts[typename], current_combo))
            else:
                self.log.write("%5d : %10d %10d | %-50s\n" % (index, typecounts[typename], molecule_typecounts[typename], current_combo))

            ntypes += typecounts[typename]
            index += 1

        nmolecules = len(self.molecules)

        if type_matches is not None:
            frac_match = (float(self.total_type_matches) / float(self.total_types))
            self.log.write("%5s : %10d %10d |  %15s %32s %8d / %8d match (%.3f %%)\n" % ('TOTAL', ntypes, nmolecules, '', '', self.total_type_matches, self.total_types, frac_match * 100))
            return frac_match
        else:
            self.log.write("%5s : %10d %10d\n" % ('TOTAL', ntypes, nmolecules))
        return 0.0

    def save_type_statistics(self, typelist, typecounts, molecule_typecounts, type_matches=None):
        """
        Collects typecount information for a csv "trajectory" output file

        Parameters
        -----------
        typelist - list of lists with form [smirks, typename]
        typecounts (dict) - number of matches for each fragment type
        molecule_typecounts (dict) - number of molecules that contain each fragment type
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference types, along with number of fragment types equivalently typed in reference molecule set.
        """
        if type_matches is not None:
            reference_type_info = dict()
            for (current_typename, reference_typename, count) in type_matches:
                reference_type_info[current_typename] = (reference_typename, count)

        index = 1
        output = []
        # Print counts
        # INDEX, SMARTS, PARENT INDEX, REF TYPE, MATCHES, MOLECULES, FRACTION, OUT of, PERCENTAGE
        ntypes = 0
        for [smarts, typename] in typelist:
            parent = str(self.parents[typename]['parent'])
            if type_matches is not None:
                (reference_typename, reference_count) = reference_type_info[typename]
                if reference_typename is not None:
                    reference_total = self.reference_type_counts[reference_typename]
                    # Save output
                    output.append("%i,'%s','%s','%s','%s',%i,%i,%i,%i" % (index, smarts, typename, parent, reference_typename, typecounts[typename], molecule_typecounts[typename], reference_count, reference_total))
                else:
                    output.append("%i,'%s','%s','%s','%s',%i,%i,%i,%i" % (index, smarts, typename, parent, 'NONE', typecounts[typename], molecule_typecounts[typename], 0, 0))

            else:
                output.append("%i,'%s','%s','%s','%s',%i,%i,%i,%i" % (index, smarts, typename, parent, 'NONE', typecounts[typename], molecule_typecounts[typename], 0, 0))
            ntypes += typecounts[typename]
            index += 1
        # add totals to the reference count
        nmolecules = len(self.molecules)
        if type_matches is None:
            output.append("-1,'total','all','None','all',%i,%i,0,0" % (ntypes, nmolecules))
        else:
            output.append("-1,'total','all','None','all',%i,%i,%i,%i" % (ntypes,nmolecules,self.total_type_matches, self.total_types))

        return output

    def write_parent_tree(self, roots, start='', verbose=False):
        """
        Recursively writes out the parent tree.

        Parameters
        ----------
        roots = list of typenames to print with this start
        strart = string to print at beginning of current line
        verbose = boolean, print to commandline
        """
        for r in roots:
            branch_string = "%s%s (%s)" % (start,r,self.parents[r]['SMIRKS'])
            self.log.write("%s\n" % (branch_string))
            if verbose: print(branch_string)
            new_roots = self.parents[r]['children']
            self.write_parent_tree(new_roots, start+'\t', verbose)

    def run(self, niterations, verbose = False):
        """
        Run sampler for the specified number of iterations.

        Parameters
        ----------
        niterations : int
            The specified number of iterations
        verbose : boolean, optional, default = False
            if True prints information for each iteration

        Returns
        ----------
        fraction_matched : float
            fraction of total types matched successfully at end of run

        """
        self.traj = open("%s.csv" % self.output,'w')
        self.traj.write('Iteration,Index,Smarts,Typename,ParentTypename,RefType,Matches,Molecules,FractionMatched,Denominator\n')
        for iteration in range(niterations):
            itinfo = "Iteration %d / %d" % (iteration, niterations)
            self.log.write(itinfo+'\n')
            if verbose: print(itinfo)

            accepted = self.sample_types()
            typelist = [[env.asSMIRKS(), env.label] for env in self.envList]
            [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)

            # Get data as list of csv strings
            lines = self.save_type_statistics(typelist, typecounts, molecule_typecounts, type_matches=self.type_matches)
            # Add lines to trajectory with iteration number:
            for l in lines:
                self.traj.write('%i,%s\n' % (iteration, l))

            if accepted:
                self.log.write('Accepted.\n')
            else:
                self.log.write("Rejected.\n")

            # Compute type statistics on molecules.
            frac_match = self.write_type_statistics(typelist, typecounts, molecule_typecounts, type_matches=self.type_matches)
            if verbose: print("Current Score: %.2f %%" % (frac_match * 100.0))
            self.log.write('\n')

            # Print parent tree as it is now.
            self.log.write("%s type hierarchy: \n" % self.typetag)
            if verbose: print("Current %s type hierarchy:" % self.typetag)
            roots = [e.label for e in self.baseTypes]
            self.write_parent_tree(roots, '\t', verbose)
            self.log.write('\n\n')
            if verbose: print('')


        #Compute final type stats
        typelist = [ [env.asSMIRKS(), env.label] for env in self.envList]
        [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        fraction_matched = self.write_type_matches(typelist, self.type_matches)

        self.log.write("%s type hierarchy: \n" % self.typetag)
        roots = [e.label for e in self.baseTypes]
        self.write_parent_tree(roots, '\t', verbose)

        # close files
        self.log.close()
        self.traj.close()
        return fraction_matched

    def write_results_smarts_file(self):
        """
        Creates a smarts file with the most recent SMIRKS patterns compared
        to references SMIRKS patterns the format follows smarts files
        lines beginning with % are comments and the final format has

        % Results for sampling (typetag) at # temperature
        SMIRKS      label
        % matched reference SMIRKS      reference label
        ...
        % Final score was # %%
        file name is *_results.smarts
            where * is the provided output base
        returns results file name
        """
        # open results file
        smarts_file = self.output+"_results.smarts"
        smarts = open(smarts_file, 'w')

        # header
        smarts.write("%% Results for sampling %ss at %.2e\n" % (self.typetag, self.temperature))
        smarts.write("%% SMIRKS patterns for final results are below\n")

        # if there is a reference SMIRFF
        if self.SMIRFF is not None:
            smarts.write("%% followed by a their matched reference SMIRKS from %s\n" % (self.SMIRFF))
            score = 100.0 * (float(self.total_type_matches) / float(self.total_types))
            smarts.write("%%Final Score was %.3f %%\n" % (score))
            smarts.write("%%\n")
            # Make dictionary to retrieve SMIRKS
            current_dict = dict()
            for env in self.envList:
                current_dict[env.label] = env.asSMIRKS()
            # loop through current type_matches
            for (current, ref, count) in self.type_matches:
                smarts.write("%-50s %-20s\n" % (current_dict[current], current))
                if ref is not None:
                    smarts.write("%% %-48s %-20s\n" % (self.reference_typename_dict[ref], ref))

        else: # no reference SMIRFF, just print current SMIRKS
            smarts.write("%% No reference SMIRFF provided\n")
            for env in self.envList:
                smarts.write("%-50s %-20s\n" % (env.asSMIRKS(), env.label))

        smarts.close()
        return smarts_file

