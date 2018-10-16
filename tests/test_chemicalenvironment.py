"""
These tests were taken from openforcefield/tests/test_chemicalenvironments.py
on GitHub at:
https://github.com/openforcefield/openforcefield/blob/master/openforcefield/tests/test_chemicalenvironment.py

Functions were updated to work with chemper's current implementation
As with the environment code itself, if/when openforcefield sufficiently
supports RDKit we will make it a dependency and then these tests will be removed.

"""
from functools import partial
import pytest
from chemper.graphs.environment import ChemicalEnvironment
from chemper.chemper_utils import is_valid_smirks

input_SMIRKS = [
    ('[#6X4:1]','Atom'),
    ('[#6X4:1]-[#6X4:2]', 'Bond'),
    ('[#6X4:1]-[#6X4:2]-[#6X4:3]', 'Angle'),
    ('[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]', 'ProperTorsion'),
    ('[#6X4:1]-[#6X4:2](-[#6X4:3])-[#6X4:4]', 'ImproperTorsion'),
    ("[#6](-[#1])-[#8]", None),
    ("[#6&X4&H0:1](-[#1])-[#6&X4]", 'Atom'),
    ("[#6&X4&H0:1](-[#1])-[#6&X4:2]", 'Bond'),
    ("[*:1]-[*:2](-[#6&X4])-[*:3]", 'Angle'),
    ("[#6&X4&H0:1](-[#1])-[#6&X4:2]-[#6&X4&H0:3](-[#1])-[#6&X4:4]", 'ProperTorsion'),
    ("[#1:1]-[#6&X4:2](-[#8:3])-[#1:4]", 'ImproperTorsion'),
    ("[#1:1]-[#6&X4:2](-[#8:3])-[*:4](-[#6&H1])-[#8:5]", None),
    ("CCC", None),
    ("[#6:1]1(-;!@[#1,#6])=;@[#6]-;@[#6]1", 'Atom'),
    ("C(O-[#7,#8])CC=[*]", None),
    ("[*:1]-[#7X3:2](-[#6a$(*1ccc(-[#8-1X1])cc1):3])-[*:4]", 'ImproperTorsion'),
    ("[#6X4:1]1~[*:2]~[*$(*~[#1]):3]1", 'Angle'),
    ("[$([#7]1~[#6]-CC1)]", None),
    ("[$(c1ccccc1)]", None),
    ("[#1$(*-[#6](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]):1]~[$([#1]~[#6])]", 'Atom')
]

@pytest.mark.parametrize('smirks, frag_type', input_SMIRKS)
def test_create_environments(smirks, frag_type):
    """
    Test all types of ChemicalEnvironment objects with defined atoms and bonds
    Each will be tetrahedral carbons connected by ring single bonds
    """
    env = ChemicalEnvironment(smirks)
    output_type = env.getType()
    assert output_type == frag_type

def test_complicated_torsion():
    """
    Test ChemicalEnvironment objects with complicated torsion
    test methods that add atoms, remove atoms
    add ORtypes and ANDtypes to existing atoms

    This is the SMIRK for the final torsion
    "[*:1] - [#6:2](=[#8,#7;H0]) - [#6:3](-[#7X3,#8X2;+0]-[#1])(-[#1]) - [*:4]"
    """
    torsion_smirks = "[*:1]-[#6:2]-[#6:3]-[*:4]"
    torsion = ChemicalEnvironment(torsion_smirks)
    # save atoms (use selectAtom)
    atom1 = torsion.selectAtom(1)
    atom2 = torsion.selectAtom(2)
    atom3 = torsion.selectAtom(3)

    # Add atoms with names so I can try to remove them
    atom2alpha = torsion.addAtom(atom2, [('=',[])], None, [('#8',[]),('#7',[])], ['H0'])
    atom3alpha1 = torsion.addAtom(atom3)
    atom3beta1 = torsion.addAtom(atom3alpha1, [('-',[])], None, [('#1',[])])
    atom3alpha2 = torsion.addAtom(atom3, [('-',[])], None, [('#1',[])])

    # Get bond for atom3 and alpha and add ANDtype
    bond = torsion.getBond(atom3, atom3alpha1)
    assert bond is not None
    bond.addORtype('-', [])

    # Add ORtypes and ANDtypes to atom3 alpha atom
    atom3alpha1.addORtype('#7', ['X3'])
    atom3alpha1.addORtype('#8', ['X2'])
    atom3alpha1.addANDtype('+0')

    # Call getAtoms and getBonds just to make sure they work
    torsion.getAtoms()
    torsion.getBonds()

    # get smarts and smirks for the large torsion
    smarts = torsion.asSMIRKS(smarts=True)
    assert is_valid_smirks(smarts)
    smirks = torsion.asSMIRKS()
    assert is_valid_smirks(smirks)


    # Try removing atoms
    # if it was labeled:
    removed = torsion.removeAtom(atom1)
    assert not removed
    removed = torsion.removeAtom(atom3alpha1)
    assert not removed
    removed = torsion.removeAtom(atom3beta1)
    assert removed


selections = [(4,True), ('Beta', True), (2, False), ('Indexed', False),
              ('Unindexed', False), ('Alpha', False), (None, False)]
@pytest.mark.parametrize("descriptor, is_none", selections)
def test_selection_by_descriptor(descriptor, is_none):
    """
    test selection by description works
    """
    angle_smirks = "[#6X3;R1:1]=,:;@[#6X3;R1;a:2](-,:;@[#7])-;!@[#8X2H1;!R:3]"
    angle = ChemicalEnvironment(angle_smirks)
    atom = angle.selectAtom(descriptor)
    bond = angle.selectBond(descriptor)

    assert (atom is None) == is_none
    assert (bond is None) == is_none


comp_dict = {'atom': [ ('all', 4), ('Indexed', 3), ('Unindexed', 1), ('Alpha', 1), ('Beta', 0)],
               'bond': [ ('all', 3), ('Indexed', 2), ('Unindexed', 1), ('Alpha', 1), ('Beta', 0)]}
comp_list = [(k, opt[0], opt[1]) for k,e in comp_dict.items() for opt in e]
@pytest.mark.parametrize('comp, option, expected_len', comp_list)
def test_get_component_list(comp, option, expected_len):
    """
    Test getting full component list works correctly
    """
    angle_smirks = "[#6X3;R1:1]=,:;@[#6X3;R1;a:2](-,:;@[#7])-;!@[#8X2H1;!R:3]"
    angle = ChemicalEnvironment(angle_smirks)
    components = angle.getComponentList(comp, option)
    assert len(components) == expected_len

def test_other_env_methods():
    """
    Test the other minor class functions for ChemicalEnvironments
    """
    angle_smirks = "[#6X3;R1:1]=,:;@[#6X3;R1;a:2](-,:;@[#7])-;!@[#8X2H1;!R:3]"
    angle = ChemicalEnvironment(angle_smirks)
    # Check is__ descriptors
    atom2 = angle.selectAtom(2)
    bond1 = angle.selectBond(1)
    alpha_atom = angle.selectAtom('Alpha')
    beta_atom = angle.addAtom(alpha_atom)
    alpha_bond = angle.getBond(atom2, alpha_atom)
    beta_bond = angle.getBond(alpha_atom, beta_atom)

    # list of lists:
    # [ [[components], [(method, expected)]], [...]]
    check_is_methods = [
            [[atom2,bond1], [(angle.isAlpha, False), (angle.isBeta, False),
                (angle.isIndexed, True), (angle.isUnindexed, False)]],
            [[alpha_atom, alpha_bond], [(angle.isAlpha, True), (angle.isIndexed, False),
                (angle.isUnindexed, True)]],
            [[beta_atom, beta_bond], [(angle.isBeta, True)]]]

    for compSet, methodList in check_is_methods:
        for comp in compSet:
            for (method, expected) in methodList:
                # same message
                classify = method(comp)
                assert classify == expected

    # Check getBond when atoms aren't bonded
    atom1 = angle.selectAtom(1)
    beta_to_atom1 = angle.getBond(beta_atom, atom1)
    assert beta_to_atom1 is None

    # Check valence: should be 3 for atom2
    val = angle.getValence(atom2)
    assert val == 3

    # Check bond order
    # For bond1 =,:;@ it should be 1.5 because order returns lowest possible
    order = bond1.getOrder()
    assert order == 1.5

    # For atom
    order = angle.getBondOrder(atom2)
    assert order == 3.5

def test_wrong_smirks_error():
    """
    Check that unparseable SMIRKS raises errors
    """
    smirks = "[*;m:1]"
    msg = "SMIRKS (%s) should not be parseable, but an environment was successfully created"
    from chemper.graphs.environment import SMIRKSParsingError
    with pytest.raises(SMIRKSParsingError):
        env = ChemicalEnvironment(smirks)
