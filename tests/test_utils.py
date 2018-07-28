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


# -------------------------------
# check molecule file parsers
# -------------------------------

# set up tests based on which package is expected
try:
    from openeye import oechem
    OE = True
except:
    OE = False

try:
    from rdkit import Chem
    RDK = True
except:
    RDK = False

mol2_abs_file = chemper_utils.get_data_path('molecules/MiniDrugBank_tripos.mol2')
mol2_rel_path = 'MiniDrugBank_tripos.mol2'
paths = [mol2_abs_file, mol2_rel_path]

# For the following functions, we will test default behavior and
# look for exceptions based on the available mol toolkit

@pytest.mark.parametrize('path', paths)
def test_mols_from_mol2(path):
    mols = chemper_utils.mols_fom_mol2(path)
    assert len(mols) == 363

    if OE:
        mols = chemper_utils.mols_fom_mol2(path, toolkit='oe')
        assert len(mols) == 363
    else:
        with pytest.raises(ImportError):
            mols = chemper_utils.mols_fom_mol2(path, toolkit='oe')

    if RDK:
        mols = chemper_utils.mols_fom_mol2(path, toolkit='rdk')
        assert len(mols) == 363
    else:
        with pytest.raises(ImportError):
            mols = chemper_utils.mols_fom_mol2(path, toolkit='rdk')


function_pairs = [
    (OE, chemper_utils.oe_mols_from_file),
    (RDK, chemper_utils.rdk_mols_from_mol2)
]

@pytest.mark.parametrize('toolkit,file_function', function_pairs)
def test_specific_mols_from_mol2(toolkit, file_function):
    if toolkit:
        mols = file_function(mol2_abs_file)
        assert len(mols) == 363

        with pytest.raises(IOError):
            file_function(mol2_rel_path)

    else:
        with pytest.raises(ImportError):
            file_function(mol2_abs_file)

