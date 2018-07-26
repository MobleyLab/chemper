.. chemper documentation master file, created by
   sphinx-quickstart on Tue May 15 13:36:15 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to chemper
==================

``chemper`` is python package that works to automatically learn chemical perception
-- that is, the way force field parameters are assigned to a molecule based on chemical environment.
Traditional force fields use atom types to assign parameters.
This effort is in conjunction with the Open Force Field Initiative
-- an academic and industrial collaboration focused on developing a machinery that could automatically
generate a classical all-atom force field given reference data and a functional form [:ref:`offinitiative`].

Molecular mechanics force fields enable a wide variety of computational chemistry calculations
by giving the energy and forces of an atomistic system as a function of the coordinates.
The literature is littered with concerns about force field limitations, but due to the considerable number of
human hours it takes to build a new force field, quantifying these limitations is nearly impossible.
In an effort to help improve these essential models, the Open Force Field Initiative is working
to build a machinery that could automatically generate a classical all-atom force field given
reference data and a functional form. This set of open source parameterization
tools will increase our ability to do force field science. In the traditional form, new general all atom
force fields take years to generate. With an automated tool, it will be easier to test choices in functional form.

Our new SMIRKS Native Open Force Field (SMIRNOFF) format replaces atom types with direct chemical perception [:ref:`escaping`].
Most current force fields use atom types for chemical perception.
Atom types are considered indirect chemical perception because first a molecule is labeled with these types;
then all other chemical information (i.e. bond orders) are removed; and parameters are then assigned using only
the atom types and their connectivity.
This means the chemistry required for all force field parameters -- Lennard-Jones interactions, bond-stretching,
angle-bending, and torsions -- must be included in those those atom types.
SMIRNOFF uses direct chemical perception in the form of SMIRKS strings.
SMIRKS are a cheminformatics language for describing molecular fragments with decorators to
distinguish characteristics of the atoms and bonds in that fragment [:ref:`daylight`,:ref:`smarts_lit`].
This allows each parameter type to be specified independently.
We have shown that the first SMIRNOFF implementation, smirnoff99Frosst [:ref:`smirnoff99frosst`],
performs comparably to GAFF and is capable of
typing significantly more molecules in both the DrugBank and eMolecules databases [:ref:`escaping`].
However, like traditional atom types, the SMIRKS patterns in smirnoff99Frosst were determined by a chemical expert.
In order to automate force field parameterization we need a way to automatically learn these SMIRKS patterns based on
reference quantum mechanical or experimental reference data, which is why we introduce ``chemper``.

``chemper`` provides

Contributors
------------

* `Caitlin C. Bannan (UCI) <https://github.com/bannanc>`_
* `Jessica Maat (UCI) <https://github.com/jmaat>`_
* `David L. Mobley (UCI) <https://github.com/davidlmobley>`_

.. toctree::
      :maxdepth: 1
      :caption: Contents:

      docs/intro
      docs/installation
      docs/api/index
      docs/examples
      docs/future
      docs/references

Acknowledgments
---------------

CCB is funded by a fellowship from
`The Molecular Sciences Software Institute <http://molssi.org/>`_
under NSF grant ACI-1547580.
JM and DLM appreciate appreciate the financial support from the
National Science Foundation (CHE 1352608) and the National Institutes of Health (1R01GM108889-01).

