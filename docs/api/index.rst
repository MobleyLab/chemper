chemper API
==============

``chemper`` is open source. If you have suggestions or concerns
please add them to our `issue tracker <https://github.com/MobleyLab/chemper/issues>`_.
Documentation for contributing to this work will be available soon.
Below is an outline of tools provided in  ``chemper`` see :doc:`examples <../examples>`
for more details.

mol_toolkits
------------

As noted in :doc:`installation <../installation>`, we seek to
keep ``chemper`` independent of the cheminformatics toolkit.
:doc:`mol_toolkits <mol>` is created to keep all code dependent
on the toolkit installed. It can create molecules from
an RDK or OE molecule object or from a SMILES string.
It includes a variety of functions for extracting information
about atoms, bonds, and molecules.
Also included here are SMIRKS pattern searches.

.. toctree::
    :maxdepth: 2

    mol.rst


chemperGraph
------------

The goal of :doc:`chemperGraph <single_graph>` was to create an example of how you could
create a SMIRKS pattern from a molecule and set of atom indices.
Creating SMIRKS for one molecule may not be useful for sampling chemical perception
in the long run. However, it is a tool that did not previously exist to the best
of our knowledge. For a detailed example, see
single_mol_smirks_.

.. _single_mol_smirks: ../../examples/single_mol_smirks.ipynb

Here is a brief usage example for creating the SMIRKS pattern
for the bond between the two carbon atoms in ethene including
atoms one bond away from the indexed atoms. The indexed atoms
are the two carbon atoms at indices 0 and 1 in the molecule, which are
specified with the tuple `(0,1)`.
These atoms are assigned to SMIRKS indices ``:1`` and ``:2`` respectively.
The variable `layers` tells `ChemPerGraph` to include atoms up to 1 bond away
from the indexed atoms.

.. literalinclude:: fragment_example.py
    :language: python

.. toctree::
    :maxdepth: 2

    single_graph.rst

ClusterGraph
------------

The goal of :doc:`ClusterGraph <cluster_graph>` is to store all
information about the atoms and bonds that could be in a SMIRKS
pattern. These are created assuming you already have a clustered
set of molecular fragments. Our primary goal is to determine
chemical perception for force field parameterization. We imagine
parameters for each molecule (for example bond lengths and force constants)
could be clustered by fragment. Then we could generate a hierarchical
list of SMIRKS patterns that maintain those clusters for typing purposes.
However, you could imagine other reasons
for wanting to store how you clustered groups of atoms.
For example, using atom or bond types in a machine learning model.

For more detailed examples and illustration of how this works see
smirks_from_molecules_ example.
Below is a brief example showing the SMIRKS for the bond
between two carbon atoms in propane and pentane.

.. _smirks_from_molecules: ../../examples/smirks_from_molecules.ipynb

.. literalinclude:: cluster_example.py
    :language: python

The idea with ClusterGraph objects is that they store all
possible decorator information for each atom. In this case the
SMIRKS indexed atoms for propane (mol1) are one of the terminal
and the middle carbons. In pentane (mol2) however the first atom can
be a terminal or middle of the chain carbon atom. This changes
the number of hydrogen atoms (``Hn`` decorator) on the carbon,
thus there are two possible SMIRKS patterns for atom ``:1``
``#6AH2X4x0r0+0`` or (indicated by the "`,`") ``#6AH3X4x0r0+0``. But, atom `:2` only has one possibility `#6AH2X4x0r0+0`.

.. toctree::
    :maxdepth: 2

    cluster_graph.rst

SMIRKSifier
------------

Lets assume you have a few clusters of fragments that you want
assigned the same force field parameter. For example, you could
have clusters of carbon-carbon bonds based on the type of bond
between them (single, double, etc).
In this case ChemPer would use the SMIRKSifier to generate
a hierarchical list of SMIRKS patterns for those clusters.
This process creates SMIRKS using ClusterGraph and then takes
a stochastic approach to removing unnecessary decorators.
See the general_smirks_for_clusters_ example for how this process could
be applied to different bonding parameters.

.. _general_smirks_for_clusters: ../../examples/general_smirks_for_clusters.ipynb

.. toctree::
    :maxdepth: 2

    smirksify.rst