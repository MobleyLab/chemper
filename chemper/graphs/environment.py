#!/usr/bin/env python

#==============================================================================
# MODULE DOCSTRING
#==============================================================================

"""
environment.py
this is adapted from the openforcefield ChemicalEnvironments class.

There have been some on going debates over where environment.py should live,
here or in the base openforcefield toolkit. Due to the want to update
this module for use in chemper amidst the openforcefield API overall,
this environment.py has been updated independently in this repository.
These updates have been fairly substantial, specifically the getter
and setter functions for decorators were removed and replaced with more
pythonic @property and @[property].setter functions instead.
All methods have also been renamed to use snake case.

The only function of real use in openforcefield.py is the `get_type` function,
which has also been updated here.
"""

#==============================================================================
# GLOBAL IMPORTS
#==============================================================================

import networkx as nx
import re
import copy

from numpy import random


#==============================================================================
# Functions
#==============================================================================

def _find_embedded_brackets(string, in_char, out_char):
    """
    Finds the substring surrounded by the in_char and out_char
    intended use is to identify embedded bracketed sequences

    For example, if you have the input
    string = "[#1$(*-C(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]):1]"
    sub_string, in_idx, out_idx = _find_embedded_brackets(string, '\(','\)')
    # sub_string = (*-C(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br])  in_idx = 4, out_idx = 50

    Parameters
    -----------
    string : str
        a string you want separated
    in_char : str
        regular expression for the character you're looking for '\(' for '('
    out_char : str
        regular expression for the closing character such as '\)' for ')'

    Returns
    --------
    substring : str
        string between the first occurances of the in_char and out_char
    in_idx : int
        index from initial string with the first in_char
    out_idx : int
        index from initial string with the first out_char
    """
    in_list = [m.start() for m in re.finditer(in_char, string)]
    out_list = [m.start() for m in re.finditer(out_char, string)]
    # If no occurance of the in_char return an empty string
    if len(in_list) == 0:
        return "", -1, -1
    # If no out_char returns the first in_char to the end
    if len(out_list) == 0:
        return string[in_list[0]:], in_list[0], -1

    # Otherwise find closure from the first in_char
    list_idx = 0
    while list_idx < len(in_list) - 1:
        if in_list[list_idx+1] > out_list[list_idx]:
            break
        list_idx+=1
    in_idx = in_list[0]
    out_idx = out_list[list_idx]
    return string[in_idx:out_idx+1], in_idx, out_idx


def _convert_embedded_smirks(smirks):
    """
    Converts a SMIRKS string with the $(...) in an atom to the
    form expected by the environment parser

    For example, if you provide initial_smirks = "[#1$(*~[#6]):1]"
    then new_smirks = _convert_embedded_smirks(initial_smirks)
    will return new_smirks = [#1:1]~[#6]

    Parameters
    -----------
    smirks : str
        any smirks string, if no $(...) then the original smirks is returned

    Returns
    --------
    updated_smirks: str
        smirks string with no recursive smirks
    """
    a_out = 0
    while smirks.find('$(') != -1:
        # Find first atom
        atom, a_in, a_out = _find_embedded_brackets(smirks, r'\[', r'\]')
        d = atom.find('$(')
        # Find atom with the $ string embedded
        while d == -1:
            atom, temp_in, temp_out = _find_embedded_brackets(smirks[a_out+1:], r'\[', r'\]')
            a_in = a_out + temp_in + 1
            a_out += temp_out + 1
            d = atom.find('$(')

        # Store the smirks pattern before and after the relevant atom
        pre_smirks = smirks[:a_in]
        post_smirks = smirks[a_out+1:]

        # Check for ring index, i.e. the 1s in "[#6:1]1-CCCCC1"
        match = re.match(r'(\d+)',post_smirks)
        if match is not None: # leftover starts with int
            ring_out = re.findall(r'(\d+)',post_smirks)[0]
            # update post_smirks
            post_smirks = post_smirks[match.end():]
        else:
            ring_out = ''

        embedded, p_in, p_out = _find_embedded_brackets(atom, r'\(', r'\)')
        # two forms of embedded strings $(*~stuff) or $([..]~stuff)
        # in the latter case the first atom refers the current atom
        if embedded[1] == '[':
            first, f_in, f_out = _find_embedded_brackets(embedded, r'\[',r'\]')
            first = _convert_embedded_smirks(first)
            new_atom = atom[:d]+first[1:-1]+atom[p_out+1:]
            embedded = embedded[f_out+1:]
            # if embedded is empty between brackets, remove it
            if embedded.replace('(','').replace(')','') == '':
                embedded = ''

        elif embedded[1] == '*': # embedded[1] = *
            new_atom = atom[:d]+atom[p_out+1:]
            embedded = embedded[2:]

        else: # embedded starts with a "no bracket" atom such as 'C'
            embedded = embedded[1:] # remove leading '('
            # atoms by symbol don't need brackets, this covers atomic symbols and aromatic atoms
            no_bracket = r'(!?[A-Z][a-z]?|!?[cnops])'
            match = re.match(no_bracket, embedded)
            if match is not None:
                new_atom = atom[:d]+embedded[:match.end()]+atom[p_out+1:]
                embedded = embedded[match.end():]
            else:
                new_atom = atom[:d]+atom[p_out+1]

        # Look for ring inside embedded SMIRKS "[#6$(*1CCC1)]"
        match = re.match(r'(\d+)', embedded)
        if match is not None: # embedded starts with an int
            ring_in = re.findall(r'(\d+)', embedded)[0]
            embedded = '(' + embedded[match.end():]
        else:
            ring_in = ''
            if embedded != '':
                embedded = '(' + embedded

        # Make new smirks
        smirks = pre_smirks+new_atom+ring_out+ring_in+embedded+post_smirks

    return smirks


def _remove_blanks_repeats(init_list, remove_list = ['']):
    """
    Returns the input list 'init_list'
    without any repeating entries or blank strings ''

    Parameters
    -----------
    init_list : list
        This is a list of anything, but intended for decorator strings
    remove_list : list
        List of things you want removed from the init_list
        For decorators we don't need empty strings

    Returns
    --------
    final_list : list
        The init_list with no duplicates and nothing from the remove_list

    TODO: this changes the order of inputs potentially so this function
          will need to be updated if order of init_list is important.
    """
    final_list = [item for item in init_list if item not in remove_list]
    return list( set(final_list) )


class SMIRKSMismatchError(Exception):
    """
    Exception for cases where smirks are inappropriate
    for the environment type they are being parsed into
    """
    def __init__(self, msg):
        super(SMIRKSMismatchError, self).__init__(self,msg)
        self.msg = msg


class SMIRKSParsingError(Exception):
    """
    Exception for when SMIRKS are not parseable for any environment
    """
    def __init__(self, msg):
        super(SMIRKSParsingError, self).__init__(self, msg)
        self.msg = msg


class ChemicalEnvironment(object):
    """Chemical environment abstract base class that matches an atom, bond, angle, etc.
    """
    class Atom(object):
        """Atom representation, which may have some or_types and ANDtypes properties.

        Attributes
        ----------
        or_types : list of tuples in the form (base, [list of decorators])
            where bases and decorators are both strings
            The descriptor types that will be combined with logical OR
        and_types : list of string
            The descriptor types  that will be combined with logical AND
        """
        def __init__(self, or_types = None, and_types = None, index = 0, ring = None):
            """Initialize an Atom object with optional descriptors.

            Parameters
            -----------
            or_types : list of tuples for ORbases and ORdecorators,
                in the form (base, [list of decorators])
                optional, default = []
            and_types : list of str
                strings that will be AND'd together in a SMARTS
                optional, default = None
            index : int
                If greater than zero,
                the specified index will be attached as a SMIRKS index (e.g. '[#6:1]')
                otherwise, it is only used for accessing atom information
            ring : int, optional, default = None
                If not None, the specified ring index will be attached at the end of the atom i.e. '[#6:1]1'
            """
            # List of 2 tuples in the form [ (ORbase, ORdecorator), ...]
            if or_types is not None:
                self._or_types = copy.deepcopy(or_types)
            else:
                self._or_types = list()

            # Set of strings that will be AND'd to the the end
            if and_types is not None:
                self._and_types = list(copy.deepcopy(and_types))
            else:
                self._and_types = list()

            self.index = index
            self.ring = ring
            self.is_atom = True

        def is_generic(self):
            """
            returns True if there are no decorators on this atom
            (IMPORTANT: this is newly added and in chemper only as of 8/9/18)
            """
            if not self._or_types:
                if not self._and_types:
                    return True
            return False

        def as_smarts(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
            The SMARTS string for this atom, meaning it has no :n index
            """

            smarts = '['

            # Add the OR'd features
            if self._or_types:
                or_list = list()
                for (base, or_decorators) in self._or_types:
                    if len(base) > 0 and base[0] == '$':
                        # after a $base an explicit '&' is necessary
                        if or_decorators:
                            or_bit = base + '&' + ''.join(or_decorators)
                        else:
                            or_bit = base
                    else: # base doesn't start with $
                        or_bit = base + ''.join(or_decorators)
                    or_list.append(or_bit)
                smarts += ','.join(or_list)
            else:
                smarts += '*'

            if len(self._and_types) > 0:
                smarts += ';' + ';'.join(self._and_types)

            if self.ring is not None:
                return smarts + ']' + str(self.ring)
            else:
                return smarts + ']'

        def as_smirks(self):
            """Return the atom representation as SMIRKS.

            Returns
            --------
            smirks : str
            The SMIRKS string for this atom, same as SMARTS, but with :n index
            """
            smirks = self.as_smarts()

            # No index specified so SMIRKS = SMARTS
            if self.index <= 0:
                return smirks

            # Add label to the end of SMARTS
            sub_string, start, end = _find_embedded_brackets(smirks, r'\[', r'\]')
            end_string = smirks[end:]
            return sub_string[:-1] + ':' + str(self.index) + end_string

        def add_or_type(self, or_base, or_decorators):
            """
            Adds ORtype to the set for this atom.

            Parameters
            -----------
            or_base : string, such as '#6'
            or_decorators : list of strings, such as ['X4','+0']
            """
            or_decorators = _remove_blanks_repeats(or_decorators, ['', or_base])
            self._or_types.append((or_base, or_decorators))

        def add_and_type(self, and_type):
            """
            Adds ANDtype to the set for this atom.

            Parameters
            --------
            and_type : string
                added to the list of and_types for this atom
            """
            self._and_types.append(and_type)
            self._and_types = _remove_blanks_repeats(self._and_types)

        @property
        def or_types(self):
            """Provides the or_types in this atom"""
            return self._or_types

        @or_types.setter
        def or_types(self, new_or_types):
            """
            sets new or_types for this atom

            Parameters
            ----------
            new_or_types : list of tuples in the form (base, [ORdecorators])
                for example : ('#6', ['X4','H0','+0']) --> '#6X4H0+0'
            """
            self._or_types = list()
            if new_or_types is not None:
                for (base, decs) in new_or_types:
                    adjusted_decs = _remove_blanks_repeats(decs, ['', base])
                    self._or_types.append((base, adjusted_decs))

        @property
        def and_types(self):
            """
            returns a copy of the list of and_types for this atom
            """
            return list(copy.deepcopy(self._and_types))

        @and_types.setter
        def and_types(self, new_and_types):
            """
            sets new and_types for this atom

            Parameters
            ----------
            new_and_types : list of strings
                strings that will be AND'd together in a SMARTS
            """
            if new_and_types is None:
                self._and_types = list()
            else:
                self._and_types = _remove_blanks_repeats(new_and_types)

    class Bond(Atom):
        """Bond representation, which may have ORtype and ANDtype descriptors.

        Attributes
        ----------
        or_types : list of tuples of ORbases and ORdecorators
            in form (base: [list of decorators])
            The ORtype types that will be combined with logical OR
        and_types : list of string
            The and_types that will be combined with logical AND

        """
        # Implementation identical to atoms apart from what is put in the asSMARTS/asSMIRKS strings

        def __init__(self, or_types = None, and_types = None, index = 0):
            """
            Parameters
            -----------
            or_types : list of tuples, optional, default = None
                tuples have form (base, [ORdecorators])
                bond descriptors that will be OR'd together in a SMARTS
            and_types : list of str, optional, default = None
                strings that will be AND'd together in a SMARTS
            index : integer, default = 0
                This is for book keeping inside environments and will not be shown in SMARTS or SMIRKS
                example: bond1 in a Bond is the bond between atom1 and atom2
            """
            super(ChemicalEnvironment.Bond,self).__init__(or_types, and_types, index)
            self.is_atom = False
            return

        def as_smarts(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
                The SMARTS string for just this atom
            """
            if self._or_types:
                or_combos = list()
                for (OR_base, OR_decorators) in self._or_types:
                    or_combos.append(OR_base + ''.join(OR_decorators))
                smarts = ','.join(or_combos)
            else:
                smarts = '~'

            if len(self._and_types) > 0:
                smarts += ';' + ';'.join(self._and_types)

            return smarts

        def as_smirks(self):
            """
            Returns
            --------
            smarts : str
                The SMIRKS string for just this bond
            """
            # the same as as_smarts()
            #    for consistency as_smarts() or as_smirks() can be called
            #    for all environment objects
            return self.as_smarts()

        def get_order(self):
            """
            Returns a float for the order of this bond
            for multiple or_types or ~ it returns the minimum possible order
            the intended application is for checking valence around a given atom

            Returns
            --------
            min_order : float
                minimum order for this bond (i.e. 1 for a '-' decorator)
            """
            # Minimum order for empty or_types is 1:
            if not self._or_types:
                return 1

            order_dict = {'~':1.,
                    '-':1., ':': 1.5, '=':2., '#':3.,
                    '!-':1.5, '!:':1., '!=':1., '!#':1.}
            order_list = [order_dict.get(base,1) for (base, decor) in self._or_types]
            return min(order_list)

    def __init__(self, smirks = None, label = None, replacements = None):
        """Initialize a chemical environment abstract base class.

        Parameters
        -----------
        smirks : string, optional
            if smirks is not None, a chemical environment is built
            from the provided SMIRKS string
        label : anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        replacements : list of lists, optional,
            [substitution, smarts] form for parsing SMIRKS
        """
        # Define the regular expressions used for all SMIRKS decorators
        # There are a limited number of descriptors for smirks string they are:
        # That is a # followed by one or more ints w/or w/o at ! in front '!#16'
        element_num = "!?[#]\d+"
        # covers element symbols, i.e. N,C,O,Br not followed by a number
        element_sym = "!?[A-Z][a-z]?"
        # covers element symbols that are aromatic:
        aro_sym = "!?[cnops]"
        # replacement strings
        replace_str = "\$\w+"
        # a or A w/ or w/o a ! in front 'A'
        aro_ali = "!?[aA]"
        # the decorators (D,H,j,r,V,X,^) followed by one or more integers
        needs_int = "!?[DjVX^]\d+"
        # R(x), +, - do not need to be followed by a integer w/ or w/o a ! 'R2'
        optional_int = "!?[RHhrx+-]\d*"
        # chirality options, "@", "@@", "@int" w/ or w/o a ! in front
        chirality = "!?[@]\d+|!?[@]@?"

        # Generate RegEx string for decorators:
        self.no_bracket_atom_reg = r'('+'|'.join([element_sym, aro_sym, replace_str])+')'
        self.atom_reg = '|'.join([element_num, aro_ali, needs_int,
                                  optional_int, chirality, replace_str,
                                  element_sym, aro_sym])
        self.atom_reg = r'('+self.atom_reg+')'

        # Define bond regular expression options below in order:
        # single, double, triple, aromatic, directional up bond, directional down bond
        # Each can have ! in from and directional can have ? after
        # up and down bonds have lots of \ to fit the python requirements
        self.bond_regs = ['!?[-]', '!?[=]', '!?[#]', '!?[:]', '!?[@]', '!?[\\\\]\\??', '!?[\\/]\\??']
        self.bond_regs = r'('+'|'.join(self.bond_regs)+')'
        # Note, not looking for ~ because that is used for empty bonds

        # Create an empty graph which will store Atom objects.
        self._graph = nx.Graph()
        self.label = label
        self.replacements = replacements

        if smirks is not None:
            # Check that it is a valid SMIRKS
            if not self.is_valid(smirks):
                raise SMIRKSParsingError("Error Provided SMIRKS ('%s') was \
not parseable with current toolkit" % smirks)

            # Check for SMIRKS not supported by Chemical Environments
            if smirks.find('.') != -1:
                raise SMIRKSParsingError("Error: Provided SMIRKS ('%s') \
contains a '.' indicating multiple molecules in the same pattern. This type \
of pattern is not parseable into ChemicalEnvironments" % smirks)
            if smirks.find('>') != -1:
                raise SMIRKSParsingError("Error: Provided SMIRKS ('%s') \
contains a '>' indicating a reaction. This type of pattern is not parseable \
into ChemicalEnvironments." % smirks)

            # try parsing into environment object
            try:
                self._parse_smirks(smirks)
            except:
                raise SMIRKSParsingError("Error SMIRKS (%s) was not parseable\
                        into a ChemicalEnvironment" % smirks)

        # Check that the created Environment is valid
        if not self.is_valid():
            raise SMIRKSParsingError("Input SMIRKS (%s), converted to %s \
                    is now invalid" % (smirks, self.as_smirks()))

        return

    def _graph_remove_node(self, node):
        """
        removes a node from the graph, kept separate from other
        functions so if (when) networkx has an API change we only
        have to change one place.

        Parameters
        -----------
        node : node in self._graph

        Returns
        --------
        node_removed : bool
        """
        if node not in self._graph:
            return False
        self._graph.remove_node(node)
        return True

    def _graph_nodes(self, data=False):
        """
        When data is False returns a list of nodes in graph
        otherwise returns a dictionary in the form {node: data}

        Parameters
        -----------
        data : bool
            include data for each node

        Returns
        --------
        nodes : list or dict
            if data is False, returns a list in the form:
                [node1, node2, ...]
            if data is True, returns a dictionary in the form:
                {node: {data_key: data, ...}, ... }
        """
        if data:
            return dict(self._graph.nodes(data=True))
        return list(self._graph.nodes())

    def _graph_edges(self, data=False, node=None):
        """
        Returns all edges (node=None) or edges associated
        with a specific node. We use a custom internal function
        so that if (when) networkx changes their API we only
        have to change one place in the script.

        Parameters
        -----------
        data : bool
            include data on edges (bonds)?
        node : graph node
            get only edges connected to that edge

        Returns
        --------
        edges : list of edges
            Returns all edges (node=None)
            or edges connected to the specified node.
            If data is False then the list has the form:
                [ (node1, node2), ... ]
            otherwise, if data is True is has the form:
                [ (node1, node2, {dictionary of data}), ...]
        """
        if node is None:
            return list(self._graph.edges(data=data))
        return list(self._graph.edges(node, data=data))

    def _graph_neighbors(self, node):
        """
        Returns a list of neighbors for the given node.
        This is done in a custom function so we have only
        one place to change if (when) networkx changes the API.

        Parameters
        -----------
        node : graph node

        Returns
        --------
        neighbors : list of neighboring nodes
        """
        return list(self._graph.neighbors(node))

    def _graph_get_edge_data(self, node1, node2):
        """
        Returns a dictionary for the data at the edged connecting
        node1 and node2 in graph. We set this in a custom function
        so we only have to change one place if (when) networkx
        changes their API.

        Parameters
        -----------
        node1 : graph node
        node2 : a different graph node

        Returns
        --------
        data_dict : dict
            dictionary of the data stored on the edge between node1 and node2
        """
        return self._graph.get_edge_data(node1, node2)

    def is_valid(self, smirks=None):
        """
        Checks if the provided SMIRKS or the one created
        by the environment is valid according to ChemPer rules.

        Parameters
        -----------
        smirks : str or None
            if None then we call self.as_smirks()

        Returns
        --------
        is_valid : bool
            True if this is a valid ChemPer SMIRKS
        """
        if smirks is None:
            smirks = self._as_smirks()
        from chemper.chemper_utils import is_valid_smirks
        return is_valid_smirks(smirks)

    def _parse_smirks(self,input_smirks):
        """
        This function converts a smirks string to a Chemical Environment
        """
        smirks = _convert_embedded_smirks(input_smirks)
        atoms = dict() # store created atom
        idx = 1 # current atom being created
        store = list() # to store indices while branching
        bonding_to = idx # which atom are we going to bond to

        atom_string, start, end = _find_embedded_brackets(smirks, r'\[', r'\]')

        if start != 0: # first atom is not in square brackets
            if start != -1:
                start_string = smirks[:start]
            else:
                start_string = smirks

            # Check for atoms not between square brackets
            split = re.split(self.no_bracket_atom_reg, start_string)
            atom_string = split[1]

            # update leftover for this condition
            if start != -1: # there is at least 1 more bracketed atom
                leftover = ''.join(split[2:])+smirks[start:]
            else:
                leftover = ''.join(split[2:])

        else: # First atom is in square brackets
            leftover = smirks[end+1:]
            # remove square brackets for parsing
            atom_string = atom_string[1:-1]

        # Check for ring index, i.e. the 1s in "[#6:1]1-CCCCC1"
        match = re.match(r'(\d+)',leftover)
        if match is not None:  # leftover starts with int
            ring = re.findall(r'(\d+)',leftover)[0]
            leftover = leftover[match.end():]
        else:
            ring = None

        # Get atom information and create first atom
        ors, ands, index = self._get_atom_info(atom_string)
        new_atom = self.add_atom(None, new_or_types= ors, new_and_types= ands,
                                 new_atom_index= index, new_atom_ring= ring, beyond_beta= True)
        atoms[idx] = new_atom

        while len(leftover) > 0:
            idx += 1
            # Check for branching
            if leftover[0] == ')':
                bonding_to = store.pop()
                leftover = leftover[1:]
                continue

            if leftover[0] == '(':
                store.append(bonding_to)
                leftover = leftover[1:]
                continue

            # find beginning and end of next [atom]
            atom_string, start, end = _find_embedded_brackets(leftover, r'\[', r'\]')

            if start != -1: # no more square brackets
                bond_string = leftover[:start]
            else:
                bond_string = leftover

            # Check for atoms not between square brackets
            bond_split = re.split(self.no_bracket_atom_reg, bond_string)
            # Next atom is not in brackets for example C in "[#7:1]-C"
            if len(bond_split) > 1:
                bond_string = bond_split[0]
                atom_string = '['+bond_split[1]+']'
                # update leftover for this condition
                if start != -1: # ther is at least 1 more bracketed atom
                    leftover = ''.join(bond_split[2:])+leftover[start:]
                else:
                    leftover = ''.join(bond_split[2:])

            else: # next atom is in the brackets [atom]
                # bond and atom string stay the same, update leftover
                leftover = leftover[end+1:]

            # Get bond and atom info
            b_or, b_and = self._get_bond_info(bond_string)
            a_or, a_and, index = self._get_atom_info(atom_string[1:-1])

            # Check for ring index, i.e. the 1s in "[#6:1]1-CCCCC1"
            match = re.match(r'(\d+)',leftover)
            if match is not None: # leftover starts with int
                ring = re.findall(r'(\d+)',leftover)[0]
                leftover = leftover[match.end():]
            else:
                ring = None

            # create new atom
            new_atom = self.add_atom(atoms[bonding_to], bond_or_types=b_or,
                                     bond_and_types=b_and, new_or_types=a_or, new_and_types=a_and,
                                     new_atom_index=index, new_atom_ring=ring, beyond_beta=True)

            # update state
            atoms[idx] = new_atom
            bonding_to = idx
        return

    def _get_atom_info(self, atom):
        """
        Parses string for one atom

        Parameters
        -----------
        atom : str
            string for one atom (the part between brackets)

        Returns
        --------
        or_types : list of tuples
            OR decorators are in the form [ (base, [decorators]), ...]
        and_types : list
        index : int
        """
        # Find atom index
        colon = atom.find(':')
        if colon == -1:
            index = None
        else:
            index = int(atom[colon+1:])
            atom = atom[:colon]

        split = atom.split(';')

        # Get and_types (and split them if they don't use ;)
        and_types = list()
        for a in split[1:]:
            and_types += re.findall(self.atom_reg, a)

        # Get or_types
        or_list = split[0].split(',')
        or_types = list()
        # Separate or_types into bases and decorators
        for OR in or_list:
            or_base, or_decors = self._separate_or_types(OR)
            if or_base is not None:
                or_types.append((or_base, or_decors))

        return or_types, and_types, index

    def _separate_or_types(self, or_type):
        """
        Separates ORtype (i.e. "#6X4R+0") into
        a base and decorators (i.e. '#6', ['X4','R','+0'] )

        Parameters
        -----------
        or_type : str
            string for one or_type

        Returns
        --------
        base : str
            #n, element symbol, or *
        decs : list
            list of decorators
        """
        # special case 1: wild card
        if or_type == '*':
            return None, []

        # if OR base is a wildcard
        if or_type[0] == '*':
            return '*', re.findall(self.atom_reg, or_type[1:])

        # Split up decorators by RegEx strings for atoms
        split = re.findall(self.atom_reg, or_type)
        if len(split) == 0:
            return None, []

        base = split[0]
        decs = _remove_blanks_repeats(split[1:], ['',base])
        return base, decs

    def _get_bond_info(self, bond):
        """
        Given bond strings returns or_types and and_types

        Parameters
        -----------
        bond : str
            string for one bond (i.e. '-,:;!@')

        Returns
        --------
        or_types : list
            list of or_type decorators, following atom tuple format
            in the '-,:;!@' example you get [ ('-', []), (':', []) ]
        and_types : list
            list of and_type decorators
            in this example you get ['!@']
        """
        # blank bond string is single or aromatic
        # empty or_types in Chemical Environments are treated as ~ bonds
        if bond == "":
            and_types = list()
            or_types = [('-', []), (':', [])]
            return or_types, and_types

        # AND types indicated by ; at the end
        split = bond.split(';')
        and_types = list()
        for a in split[1:]:
            and_types += re.findall(self.bond_regs, a)

        # or_types are divided by ,
        or_list = split[0].split(',')
        or_types = list()
        for OR in or_list:
            if OR == '~':
                continue
            or_divide = re.findall(self.bond_regs, OR)
            if len(or_divide) > 0:
                or_types.append((or_divide[0], or_divide[1:]))

        return or_types, and_types

    def as_smirks(self, smarts = False):
        """
        Returns a SMIRKS representation of the chemical environment

        Parameters
        -----------
        smarts : optional, boolean
            if True, returns a SMARTS instead of SMIRKS without index labels

        Returns
        --------
        smirks : str
            SMIRKS string for this environment
        """
        init_atom = self.select_atom(1)
        return self._as_smirks(init_atom, None, smarts)

    def _as_smirks(self, initial_atom = None, neighbors = None, smarts = False):
        """
        Return a SMIRKS representation of the chemical environment.
        This uses a recursive structure to combine SMIRKS for every
        atom in this environment.
        TODO: figure out if this can be done with a while loop instead

        Parameters
        -----------
        inital_atom : optional, atom object
            This is randomly selected if not chosen.
        neighbors : optional, list of atom objects
            This is all of the initalAtom neighbors if not specified
            generally this is used only for the recursive calls
            so initial atoms are not reprinted
        smarts : optional, boolean
            if True, returns a SMARTS string instead of SMIRKS
        """
        # If empty chemical environment
        if len(self._graph_nodes()) == 0:
            return ""

        if initial_atom is None:
            initial_atom = self.get_atoms()[0]

        if neighbors is None:
            neighbors = self._graph_neighbors(initial_atom)

        # sort neighbors to guarantee order is constant
        neighbors = sorted(neighbors, key=lambda atom: atom.as_smirks())

        # initialize smirks for starting atom
        if smarts:
            smirks = initial_atom.as_smarts()
        else:
            smirks = initial_atom.as_smirks()

        # loop through neighbors
        for idx, neighbor in enumerate(neighbors):
            # get the SMIRKS for the bond between these atoms
            # bonds are the same if smarts or smirks
            bond_edge = self._graph_get_edge_data(initial_atom, neighbor)
            bond_smirks = bond_edge['bond'].as_smirks()

            # Get the neighbors for this neighbor
            new_neighbors = self._graph_neighbors(neighbor)
            # Remove initialAtom so it doesn't get reprinted
            new_neighbors.remove(initial_atom)

            # Call asSMIRKS again to get the details for that atom
            atom_smirks = self._as_smirks(neighbor, new_neighbors, smarts)

            # Use ( ) for branch atoms (all but last)
            if idx < len(neighbors) - 1:
                smirks += '(' + bond_smirks + atom_smirks + ')'
            # This is for the atoms that are a part of the main chain
            else:
                smirks += bond_smirks + atom_smirks

        return smirks

    def select_atom(self, descriptor=None):
        """
        Select a random atom fitting the provided descriptor.

        Parameters
        ----------
        descriptor : optional, None
            describe what type of atom you want with the
            follow options:
                None - returns any atom with equal probability
                int - will return an atom with that index
                'Indexed' - returns a random indexed atom
                'Unindexed' - returns a random unindexed atom
                'Alpha' - returns a random alpha atom
                'Beta' - returns a random beta atom

        Returns
        --------
        atom : Atom
            one atom from this chemical environment which
            fits the provided description. If no atom matched
            the description then None is returned.
        """
        if descriptor is None or isinstance(descriptor,str):
            atoms = self.get_component_list('atom', descriptor)
            if len(atoms) == 0:
                return None
            return random.choice(atoms)

        if not isinstance(descriptor, int):
            return None

        for atom in self.get_atoms():
            if atom.index == descriptor:
                return atom
        return None

    def get_component_list(self, component_type, descriptor = None):
        """
        Returns a list of atoms or bonds matching the descriptor

        Parameters
        -----------
        component_type : string: 'atom' or 'bond'
        descriptor : string, optional
            'all', 'Indexed', 'Unindexed', 'Alpha', 'Beta'

        Returns
        -------
        component_list : list
            list of atoms or bonds (from component type)
            which match the provided descriptor
        """
        des_list = ['indexed', 'unindexed', 'alpha', 'beta', 'all']
        if descriptor is not None:
            d = descriptor.lower()
            if d not in des_list:
                raise LookupError("Error: descriptor must be in the list [%s]" %
                                ', '.join(des_list))
        else:
            d = None

        if not component_type.lower() in ['atom', 'bond']:
            raise LookupError("Error: component_type must be 'atom' or 'bond'")

        if component_type.lower() == 'atom':
            if d == 'indexed':
                return self.get_indexed_atoms()
            elif d == 'unindexed':
                return self.get_unindexed_atoms()
            elif d == 'alpha':
                return self.get_alpha_atoms()
            elif d == 'beta':
                return self.get_beta_atoms()
            else:
                return self.get_atoms()

        elif component_type.lower() == 'bond':
            if d == 'indexed':
                return self.get_indexed_bonds()
            elif d == 'unindexed':
                return self.get_unindexed_bonds()
            elif d == 'alpha':
                return self.get_alpha_bonds()
            elif d == 'beta':
                return self.get_beta_bonds()

            return self.get_bonds()

        return None

    def select_bond(self, descriptor = None):
        """Select a random bond fitting the descriptor.

        Parameters
        ----------
        descriptor : optional, None
            describe what type of atom you want with the
            follow options:
                None - returns any bond with equal probability
                int - will return an bond with that index
                'Indexed' - returns a random indexed bond
                'Unindexed' - returns a random unindexed bond
                'Alpha' - returns a random alpha bond
                'Beta' - returns a random beta bond

        Returns
        --------
        bond : Bond
            one bond from this chemical environment which
            fits the provided description. If no bond matched
            the description then None is returned.
        """
        if descriptor is None or isinstance(descriptor,str):
            bonds = self.get_component_list('bond', descriptor)
            if len(bonds) == 0:
                return None
            return random.choice(bonds)

        if not isinstance(descriptor, int):
            return None
        for bond in self.get_bonds():
            if bond.index == descriptor:
                return bond

        return None

    def add_atom(self, bond_to_atom, bond_or_types = None, bond_and_types = None,
                 new_or_types = None, new_and_types = None, new_atom_index = None,
                 new_atom_ring = None, beyond_beta = False):
        """Add an atom to the specified target atom.

        Parameters
        -----------
        bond_to_atom : atom object, required
            atom the new atom will be bound to
        bond_or_types : list of tuples, optional
            strings that will be used for the or_types for the new bond
        bond_and_types : list of strings, optional
            strings that will be used for the and_types for the new bond
        new_or_types : list of strings, optional
            strings that will be used for the or_types for the new atom
        new_and_types : list of strings, optional
            strings that will be used for the and_types for the new atom
        new_atom_index : int, optional
            integer label that could be used to index the atom in a SMIRKS string
        new_atom_ring : int, optional
            integer used to track the openning and closing of rings in SMARTS/SMIRKS patterns
        beyond_beta : boolean, optional
            if True, allows bonding beyond beta position

        Returns
        --------
        new_atom : Atom
            atom object for the newly created atom
        """
        if bond_to_atom is None:
            if len(self._graph_nodes()) > 0:
                return None

            if new_atom_index is None:
                new_atom_index = 0

            new_atom = self.Atom(new_or_types, new_and_types, new_atom_index, new_atom_ring)
            self._graph.add_node(new_atom)
            return new_atom

        # Check if we can get past beta position
        bond_to_index = bond_to_atom.index
        if bond_to_index < 0 and not beyond_beta:
            return None

        # determine the type integer for the new atom and bond
        if new_atom_index is None:
            if bond_to_index > 0:
                new_atom_index = 0
            else:
                new_atom_index = bond_to_index - 1

        if new_atom_index > 0 and bond_to_index > 0:
            bond_index = max(new_atom_index, bond_to_index) - 1
        else:
            bond_index = new_atom_index

        # create new bond
        new_bond = self.Bond(bond_or_types, bond_and_types, bond_index)

        # create new atom
        new_atom = self.Atom(new_or_types, new_and_types, new_atom_index, new_atom_ring)

        # Add node for new_atom
        self._graph.add_node(new_atom)

        # Connect original atom and new atom
        self._graph.add_edge(bond_to_atom, new_atom, bond = new_bond)

        return new_atom

    def remove_atom(self, atom, only_empty=False):
        """Remove the specified atom from the chemical environment.
        if the atom is not indexed for the SMIRKS string or
        used to connect two other atoms.

        Parameters
        ----------
        atom : atom object, required
            atom to be removed if it meets the conditions.
        only_empty : boolean, optional
            True only an atom with no and_types and 1 ORtype can be removed

        Returns
        --------
        removed : bool
            atom was removed, False: atom was not removed
        """
        # labeled atoms can't be removed
        if atom.index > 0:
            return False

        # Atom connected to more than one other atom cannot be removed
        if len(self._graph_neighbors(atom)) > 1:
            return False

        # if you can remove "decorated atoms" remove it
        if not only_empty:
            self._graph_remove_node(atom)
            return True

        if len(atom.ANDtypes) > 0:
            return False
        elif len(atom.ORtypes) > 1:
            return False

        self._graph_remove_node(atom)
        return True

    def get_atoms(self):
        """
        Returns
        -------
        atoms : list
            list of Atoms in the environment
        """
        return self._graph_nodes()

    def get_bonds(self, atom=None):
        """
        Parameters
        ----------
        atom : Atom, optional
            returns bonds in the environment, if atom is not None
            then it returns bonds connected to that atom.

        Returns
        --------
        bonds : list
            Bond objects in this environment or connected to
            the specified atom (if atom is not None)
        """
        if atom is None:
            edge_list = self._graph_edges(data=True)
            bonds = [data['bond'] for a1, a2, data in edge_list]
        else:
            bonds = []
            for (a1, a2, info) in self._graph_edges(data=True, node=atom):
                bonds.append(info['bond'])

        return bonds

    def get_bond(self, atom1, atom2):
        """
        Get bond between two atoms

        Parameters
        -----------
        atom1 : Atom
        atom2 : Atom

        Returns
        --------
        bond : Bond
            Bond between atoms 1 and 2, None if there
            is no bond connecting these atoms.
        """
        if atom2 in self._graph_neighbors(atom1):
            return self._graph_get_edge_data(atom1, atom2)['bond']
        else:
            return None

    def get_indexed_atoms(self):
        """
        Returns
        --------
        indexed_atoms : list
            The list of indexed Atom objects in this environment
        """
        index_atoms = []
        for atom in self.get_atoms():
            if atom.index > 0:
                index_atoms.append(atom)
        return index_atoms

    def get_unindexed_atoms(self):
        """
        Returns
        --------
        unindexed_atoms : list
            A list of unindexed Atom objects
        """
        unindexed_atoms = []
        for atom in self.get_atoms():
            if atom.index < 1:
                unindexed_atoms.append(atom)
        return unindexed_atoms

    def get_alpha_atoms(self):
        """
        Returns
        --------
        alpha_atoms : list
            A list atoms which are alpha to any indexed atom
        """
        alpha_atoms = []
        for atom in self.get_atoms():
            if atom.index == 0:
                alpha_atoms.append(atom)

        return alpha_atoms

    def get_beta_atoms(self):
        """
        Returns
        --------
        beta_atoms : list
            A list of atoms which are beta to any indexed atom
        """
        beta_atoms = []
        for atom in self.get_atoms():
            if atom.index == -1:
                beta_atoms.append(atom)
        return beta_atoms

    def get_indexed_bonds(self):
        """
        Returns
        --------
        indexed_bonds : list
            A list of Bond objects that connect two indexed atoms
        """
        indexed_bonds = []
        for bond in self.get_bonds():
            if bond.index > 0:
                indexed_bonds.append(bond)
        return indexed_bonds

    def get_unindexed_bonds(self):
        """
        Returns
        --------
        unindexed_bonds : list
            A list of Bond objects that connect
            an indexed atom to an unindexed atom
            two unindexed atoms
        """
        unindexed_bonds = []
        for bond in self.get_bonds():
            if bond.index < 1:
                unindexed_bonds.append(bond)
        return unindexed_bonds

    def get_alpha_bonds(self):
        """
        Returns
        --------
        alpha_bonds : list
            A list of Bond objects that connect
            an indexed atom to an alpha atom
        """
        alpha_bonds = []
        for bond in self.get_bonds():
            if bond.index == 0:
                alpha_bonds.append(bond)
        return alpha_bonds

    def get_beta_bonds(self):
        """
        Returns
        --------
        beta_bonds : list
            A list of Bond objects that connect
            an alpha atom to a beta atom
        """
        beta_bonds = []
        for bond in self.get_bonds():
            if bond.index == -1:
                beta_bonds.append(bond)
        return beta_bonds

    def is_alpha(self, component):
        """
        Checks if the component is alpha

        Parameters
        -----------
        component : Atom or Bond

        Returns
        --------
        is_alpha : bool
            True if this is an alpha atom or
            a bond which connects an indexed atom and alpha atom
        """
        return component.index == 0

    def is_unindexed(self, component):
        """
        Checks if component is not indexed

        Parameters
        -----------
        component : Atom or Bond

        Returns
        --------
        is_unindexed : bool
            True if the atom is not indexed or it is a bond that does not
            connect two indexed atoms.
        """
        return component.index < 1

    def is_indexed(self, component):
        """
        Checks if component is indexed

        Parameters
        -----------
        component : Atom or Bond

        Returns
        --------
        is_indexed : bool
            True if the atom is indexed or the bond connects
            two indexed atoms.
        """
        return component.index > 0

    def is_beta(self, component):
        """
        Checks if the component is beta

        Parameters
        -----------
        component : Atom or Bond

        Returns
        --------
        is_beta : bool
            True if the atom is beta to an indexed atom or
            the bond connects an alpha and beta atom
        """
        return component.index == -1

    def get_type(self):
        """
        Uses number of indexed atoms and bond connectivity
        to determine the type of chemical environment

        Returns
        -------
        chemical_environment_type : str
            'Atom', 'Bond', 'Angle', 'ProperTorsion', 'ImproperTorsion'
            None if number of indexed atoms is 0 or > 4
        """
        index_atoms = self.get_indexed_atoms()
        natoms = len(index_atoms)

        if natoms == 0 or natoms > 4:
            return None

        # one indexed atom can only be atom
        if natoms == 1:
            return "Atom"

        # if 2 indexed atoms check they are connected
        bond12 = self.get_bond(self.select_atom(1), self.select_atom(2))
        if bond12 is None:
            return None
        if natoms == 2:
            return "Bond"

        # if 3 indexed atoms check connected in an angle
        bond23 = self.get_bond(self.select_atom(2), self.select_atom(3))
        if bond23 is None:
            return None
        if natoms == 3:
            return "Angle"

        # if 4 indexed atoms check for proper or improper torsion connections
        bond34 = self.get_bond(self.select_atom(3), self.select_atom(4))
        bond24 = self.get_bond(self.select_atom(2), self.select_atom(4))
        if bond24 is not None:
            return "ImproperTorsion"
        if bond34 is not None:
            return "ProperTorsion"

        # must have 4 indexed atoms with incorrect connectivity
        return None

    def get_neighbors(self, atom):
        """
        Finds neighbors for the specified atom

        Parameters
        -----------
        atom : Atom

        Returns
        --------
        neighbors : list
            A list of Atoms that are directly connected to the specified atom
        """
        return self._graph_neighbors(atom)

    def get_valence(self, atom):
        """
        Calculated the valence (number of neighboring atoms)

        Parameters
        -----------
        atom : Atom

        Returns
        --------
        valence : int
            Number of Atoms neighboring the specified atom
        """
        return len(self._graph_neighbors(atom))

    def get_bond_order(self, atom):
        """
        Calculates the minimum bond order around a given atom.
        This could be used to check if the atom has too many bonds.
        For example, for a carbon atom you cannot have a
        total bond order greater than 4, so if the minimum bond order
        is greater than 4 you know there is a problem.

        Parameters
        -----------
        atom : Atom

        Returns
        --------
        bond_order : float
            Minimum total bond order around the specified atom
            0 if atom has no neighbors
            aromatic bonds count as 1.5
            any bond counts as 1.0
        """
        order = 0.
        for a1, a2, info in self._graph_edges(data=True, node=atom):
            order += info['bond'].get_order()
        return order

