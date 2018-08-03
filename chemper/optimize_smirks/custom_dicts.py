#!/usr/bin/env python

# ===================================================================
# MODULE DOCSTRING
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
        return(key)

