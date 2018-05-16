Mol Toolkits
============

``ChemPer`` wraps around existing cheminformatics packages
to extract information about atoms and bonds stored in molecules.
These are accessed through `chemper.mol_toolks.mol_toolkit` which
accesses the source code for only the installed cheminformatics package.
Right now we support OpenEye Toolkits

.. py:module::chemper.mol_toolkits.cp_openeye

and RDKit

.. py:module::chemper.mol_toolkits.cp_rdk

These and any future packge support follow the template laid
out as adapters.

.. automodule:: chemper.mol_toolkits.adapters
    :members:
