ChemPer API
==============

``Chemper`` is open source if you have suggestions or concerns
please add them to our `issue tracker <https://github.com/MobleyLab/chemper/issues>`_.
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


ChemPerGraph
------------

The goal of :doc:`ChemPerGraph <single_graph>` was to create an example of how you could
create a SMIRKS pattern from a molecule and set of atom indices.
While this isn't ultimately useful in sampling chemical
perception as they only work for a single molecule, however it
is a tool that did not exist to the best of the authors
knowledge before. For a detailed example see the
single_mol_smirks_ example.

.. _single_mol_smirks: ../../examples/single_mol_smirks.ipynb

Here is a brief usage example for creating the SMIRKS pattern
for the bond between the two carbon atoms in ethene including
atoms one bond away from the indexed atoms. The indexed atoms
are the two carbon atoms at indices 0 and 1 in the molecule
are assigned to SMIRKS indices ``:1`` and ``:2`` respectively

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
set of molecular subgraphs. As our primary goal is to determine
chemical perception for force field parameterization we imagine
the input data being clustered subgraphs based on what parameter
we wish to assign those atoms, such as equilibrium bond lengths
and force constants. However, you could imagine other reasons
for wanting to store how you clustered groups of atoms.

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
and the middle carbons. In pentane (mol2) however atom1 can
be a terminal or middle of the chain carbon atom. This changes
the number of hydrogen atoms (``Hn`` decorator) on the carbon,
thus there are two possible SMIRKS patterns for atom ``:1``
``#6AH2X4x0r0+0`` or (indicated by the "``,``") ``#6AH3X4x0r0+0``. But, atom `:2` only has one possibility `#6AH2X4x0r0+0`.

.. toctree::
    :maxdepth: 2

    cluster_graph.rst

