.. chemper documentation master file, created by
   sphinx-quickstart on Tue May 15 13:36:15 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ChemPer
==================

`ChemPer` is a Python module for creating custom chemical perception based on
user provided reference data or specified molecular examples.
Chemical perception is typically defined in the context of molecular mechanics force fields
and refers to the way parameters are assigned to a molecule based on chemical environment.

For example, consider the highlighted angles in the molecules above.
You could assume these types of assignment should be taken generally.
That is angles around carbon could always be categorized as
tetrahedral ([insert color] above) and trivalent ([insert color] above).
Using the molecules above as training you would want to recover the same typing in these molecules:

`ChemPer` can automatically generate a hierarchical list of
general SMIRKS patterns using the input molecules which will retain these types.
SMIRKS are a cheminformatics language for describing molecular fragments [:ref:`daylight`, :ref:`smartslit`].
A hierarchical list of SMIRKS can be used to recover typing rules by storing the last SMIRKS to match
a set of atoms as the type to assign that set.

Available Now
--------------
You can install `ChemPer` now from our `GitHub <https://github.com/Mobleylab/chemper/>`_
page. See :doc:`docs/installation` for more details.

A good place after installing is to checkout the :doc:`docs/examples` which are also included in the GitHub repository.

`ChemPer`'s key ingrediant is `smirksify` which takes input molecules.
`smirksify` takes user specified clustering (such as the colored angles above)
and generates a general list of SMIRKS patterns to act as typing rules.
Other modules include `ChemPerGraph` which will generate a SMIRKS pattern for a single molecule and
`ClusterGraph` which generates a very specific SMIRKS pattern for one group of molecular fragments.
Also included are a number of utility functions for using SMIRKS patterns with molecules.

Contributors
------------

* `Caitlin C. Bannan (UCI) <https://github.com/bannanc>`_
* `Jessica Maat (UCI) <https://github.com/jmaat>`_
* `David L. Mobley (UCI) <https://github.com/davidlmobley>`_

.. toctree::
      :maxdepth: 1
      :caption: Contents:

      docs/installation
      docs/api/index
      docs/examples
      docs/future
      docs/background
      docs/references

Coming Soon
------------
Currently, `ChemPer` requires the user to know how they want to type their molecules.
To implement the example above, the user would need to know they want to
group carbons by tetrahedral and trivalent.
Soon `ChemPer` will also provide examples for how to cluster those angles based on
quantitative data, such as angle in a minimized conformer.

Acknowledgments
---------------

This project is a part of the larger `Open Force Field Initiative <http://openforcefield.org>`_
We are thankful for the useful conversations about chemical perception with many Initiative members.
CCB is funded by a fellowship from
`The Molecular Sciences Software Institute <http://molssi.org/>`_
(MolSSI) under NSF grant ACI-1547580.
We would also like to thank Jessica Nash and Daniel G. A. Smith
from MolSSI for their help on this project from the beginning.
JM and DLM appreciate financial support from the
National Science Foundation (CHE 1352608) and the National Institutes of Health (1R01GM108889-01).

