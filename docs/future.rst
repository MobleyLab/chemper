What's Next?
============

The ClusterGraph code can accurately create a SMIRKS pattern
for a group of clustered molecular subgraphes.
However as you can see in the :dod:`examples <../examples>`
the SMIRKS created by ``ClusterGraph`` are highly specific.
Our final goal here is to maintain a given set of clustering,
but to generate relatively general SMIRKS patterns
that can then be used to assign force field parameters.
This use of relatively generic SMIRKS patterns is part of
what we believe makes the SMIRNOFF force field format so
powerful.


One option here would be to go back to the MC sampling used
in `SMARTY/SMIRKY <https://github.com/openforcefield/smarty>`_.
Where the information stored in ``ClusterGraph`` could be
used to make more intelligent/efficient moves.
However, we believe there is a yet more efficient option
where the differences and similarities in ``ClusterGraph``
objects could be used to determine the correct SMIRKS
patterns. Thus, the next step is to essentially find the
similarities and differences in ``ClusterGraph`` so that the
most general SMIRKS can be used to maintain clustering as
the user innputs, but not specify more information than
necessary.

