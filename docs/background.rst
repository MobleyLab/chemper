Background
============

``ChemPer`` is python package that works to automatically learn chemical perception
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
distinguish characteristics of the atoms and bonds in that fragment [:ref:`daylight`, :ref:`smartslit`].
This allows each parameter type to be specified independently.
We have shown that the first SMIRNOFF implementation, smirnoff99Frosst [:ref:`smirnoff99frosst`],
performs comparably to GAFF and is capable of
typing significantly more molecules in both the DrugBank and eMolecules databases [:ref:`escaping`].
However, like traditional atom types, the SMIRKS patterns in smirnoff99Frosst were determined by a chemical expert.
In order to automate force field parameterization we need a way to automatically learn these SMIRKS patterns based on
reference quantum mechanical or experimental reference data, which is why we introduce ``ChemPer``.



``ChemPer`` contains a variety of tools that will be useful in
automating the process of chemical perception for the new
SMIRKS Native Open Force Field (`SMIRNOFF <https://github.com/openforcefield/openforcefield>`_)
format as a part of the `Open Force Field Consortium <http://openforcefield.org>`_ :ref:`[1] <refs>`.

This idea originated from the tools `SMARTY and SMIRKY <https://github.com/openforcefield/smarty>`_
which were designed to use an automated monte carlo algorithm to
sample the chemical perception used in existing force fields.
SMARTY samples SMARTS patterns corresponding to traditional atom
types and was tested compared to the `parm99/parm@frosst <http://www.ccl.net/cca/data/parm_at_Frosst/>`_
force field. SMIRKY is an extension of SMARTY created to sample SMIRKS
patterns corresponding to SMIRNOFF parameter types
(nonbonded, bond, angle, and proper and improper torsions).
SMIRKY was tested against our own `smirnoff99Frosst <https://github.com/openforcefield/smirnoff99Frosst>`_

One of the most important lessons learned while testing SMARTY
and SMIRKY is that the combinatorial problem in SMIRKS space is
very large. These tools currently use very naive moves in SMIRKS
space chosing atom or bond decorators to add to a pattern at
complete random. This wastes a signficant amount of time making
chemically insensible moves. One of the take away conclusions
on that project was that future chemical perception sampling
tools would need to take atom and bond information from input
molecules in order to be feasible :ref:`[2] <refs>`.

We developed ``ChemPer`` based on the knowledge of the SMARTY
project outcomes. The goal here is to take clustered molecular
subgraphs and generate SMIRKS patterns. These tools will use
information stored in the atoms and bonds of a molecule to drive
choices about SMIRKS decorators. Then will automatically
generate reasonable SMIRKS patterns matching clustering of
molecular subgraphs.

