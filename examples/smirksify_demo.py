from chemper.mol_toolkits.mol_toolkit import Mol
from chemper.smirksify import SMIRKSifier, print_smirks

# make molecules from smiles
mols = [
    Mol.from_smiles('CCO'),
    Mol.from_smiles('CC=C')
]
# make clusters for each of the 6 bond types:
# carbon-carbon single bond 1 ethanol, 1 propene
cc_single = ('cc_single', [ [(0,1)],  [ (0,1) ] ]  )
# carbon-carbon double bond 0 ethanol, 1 propene
cc_double = ('cc_double', [ [ ], [ (1,2) ] ]  )
# carbon-oxygen bond 1 ethanol, 0 propene
co = ('co', [ [ (1,2) ], [ ] ]  )
# hydrogen-tetrahedral carbon 5 ethanol, 3 propene
hc_tet = ('hc_tet', [ [ (0,3), (0,4), (0,5), (1,6), (1,7) ],  [ (0,3), (0,4), (0,5)] ]  )
# hydrogen-planar carbon bond 0 ethanol, 3 propene
hc_plan = ('hc_plan', [ [ ], [ (1,6), (2,7), (2,8)] ])
# hydrogen-oxygen bond 1 ethanol, 0 propene
ho = ('ho', [ [ (2,8) ],  [ ] ] )

# initiate SMIRKSifier with default max_layers = 5 and verbose = True
fier = SMIRKSifier(mols, [cc_single, cc_double, co, hc_tet, hc_plan, ho])
# print initial SMIRKS
print_smirks(fier.current_smirks)
# Reduce SMIRKS with default 1000 iterations
fier.reduce()
# print final SMIRKS
print_smirks(fier.current_smirks)
