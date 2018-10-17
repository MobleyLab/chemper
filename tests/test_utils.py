"""
This script provides tests for all functions in chemper_utils
"""
from chemper import chemper_utils
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


