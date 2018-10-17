"""
This script provides tests for all functions in chemper_utils
"""
from chemper import chemper_utils
from chemper.mol_toolkits import mol_toolkit
import pytest

# -----------------------
# Check SMIRKS validity
# -----------------------
smirks_checks = [
    ("[#6AH3X4x0!r+0:1]-;!@[#6AH3X4x0!r+0:2]", True),
    ("C1CC1", True),
    ("[#6AH3X4x0!r+0]-;!@[#6AH3X4x0!r+0]", True),
    ("#6X4H3", False),
    ("]-;=[", False),
]

@pytest.mark.parametrize('smirks, is_valid', smirks_checks)
def test_smirks_validity(smirks, is_valid):
    assert chemper_utils.is_valid_smirks(smirks) == is_valid


# -------------------------------
# check file locator functions
# -------------------------------

chemper_data = [
    'molecules/MiniDrugBank_tripos.mol2',
    'smarts_files/angle_smirks.smarts',
    'smarts_files/bond_smirks.smarts',
    'smarts_files/improper_torsion_smirks.smarts',
    'smarts_files/proper_torsion_smirks.smarts',
    'smarts_files/nonbond_smirks.smarts'
]
@pytest.mark.parametrize('fn', chemper_data)
def test_valid_files(fn):
    import chemper
    import os
    ref_path = os.path.join(chemper.__path__[0], 'data')
    ref_path = os.path.join(ref_path, fn)

    data_path = chemper_utils.get_data_path(fn)
    assert data_path == ref_path

    full_path = chemper_utils.get_full_path(fn)
    assert full_path == data_path
    assert full_path == ref_path

def test_failing_files():
    fn = 'zzz.zzz'
    with pytest.raises(IOError):
        chemper_utils.get_full_path(fn)

    with pytest.raises(IOError):
        chemper_utils.get_data_path(fn)

smirks1 = [
    ('any', "[*:1]~[*:2]"),
    ('single', "[*:1]-[*:2]"),
    ('double', "[*:1]=[*:2]"),
    ('aromatic', "[*:1]:[*:2]"),
    ('triple', "[*:1]#[*:2]")
]
smirks2 = [
    ('any2', "[a,A:1]~[a,A:2]"),
    ('single2', "[a,A:1]-[a,A:2]"),
    ('double2', "[a,A:1]=[a,A:2]"),
    ('aromatic2', "[a:1]:[a:2]"),
    ('triple2', "[!#5:1]#[!#5:2]")
]
smirks3 = [
    ('any3', "[a,A:1]~[a,A:2]"),
    ('double3', "[a,A:1]=[a,A:2]"),
    ('aromatic3', "[a:1]:[a:2]"),
]
smirks_match_sets = [(smirks1, smirks2, True),
                     (smirks1, smirks3, False),
                     (smirks3, smirks2, False)
                     ]

@pytest.mark.parametrize('smirks1,smirks2,checks',smirks_match_sets)
def test_matching_smirks(smirks1, smirks2, checks):
    path = chemper_utils.get_data_path('molecules/MiniDrugBank_tripos.mol2')
    mols = mol_toolkit.mols_from_mol2(path)
    did_it_work = chemper_utils.check_smirks_agree(smirks1, smirks2, mols)
    assert did_it_work == checks


# test scoring function
d1 = {
    0: {(1,):'a', (2,):'b'},
    1: {(1,):'b', (2,):'c'}
}
d2 = {
    0: {(1,):'a2', (2,):'b2'},
    1: {(1,):'b2', (2,):'c2'}

}
d3 = {
    0: {(1,):'a3', (2,):'b3'},
    1: {(1,):'b3' }
}

# options:
pairs = [
    # d1 and d2 should both be able to be references to each other
    (d1, d2, {('a', 'a2'), ('b', 'b2'), ('c', 'c2')}, True),
    (d2, d1, {('a2', 'a'), ('b2', 'b'), ('c2', 'c')}, True),
    (d3, d1, set(), False),
]

@pytest.mark.parametrize('cur,ref,matches,checks',pairs)
def test_match_reference(cur,ref,matches,checks):
    type_matches, did_it_work = chemper_utils.match_reference(cur,ref)
    assert type_matches == matches
    assert did_it_work == checks