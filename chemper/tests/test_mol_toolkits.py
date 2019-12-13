"""
This is a test for the molecule, atom, and bond objects.
First, it checks which toolkits (OpenEye or RDKit) are available
It will run tests for each object for all toolkits that are available.
"""

from chemper.mol_toolkits import mol_toolkit
from chemper import chemper_utils
import pytest
import itertools

# make a list of mol toolkits
mts = list()
if mol_toolkit.HAS_RDK:
    from chemper.mol_toolkits import cp_rdk
    mts.append(cp_rdk)
if mol_toolkit.HAS_OE:
    from chemper.mol_toolkits import cp_openeye
    mts.append(cp_openeye)

if mol_toolkit.HAS_RDK or mol_toolkit.HAS_OE:
    mts.append(mol_toolkit)


@pytest.mark.parametrize('toolkit', mts)
def test_molecule(toolkit):
    """
    Test MolOE functions
    """
    mol = toolkit.Mol.from_smiles('C')

    atoms = 0
    for a in mol.get_atoms():
        atoms += 1
    assert atoms == 5

    bonds = 0
    for b in mol.get_bonds():
        bonds += 1
    assert bonds == 4

    carbon = mol.get_atom_by_index(0)
    bond = mol.get_bond_by_index(0)

    smiles = mol.get_smiles()
    assert smiles == "C"


@pytest.mark.parametrize('toolkit', mts)
def test_smirks_search(toolkit):
    """
    test SMIRKS searching
    """
    mol = toolkit.Mol.from_smiles('C')

    # smirks for C-H bond
    smirks = "[#6:1]-[#1:2]"

    matches = mol.smirks_search(smirks)

    assert len(matches) == 4

    for match in matches:
        assert 1 in match
        assert 2 in match


@pytest.mark.parametrize('toolkit', mts)
def test_bad_smirks(toolkit):
    """
    Check a ValueError is raised with improper SMIRKS
    """
    mol = toolkit.Mol.from_smiles('C')

    with pytest.raises(ValueError):
        mol.smirks_search(']X[')


@pytest.mark.parametrize('toolkit', mts)
def test_bad_smiles(toolkit):
    """
    Check a ValueError is raised with a bad SMILES
    """
    with pytest.raises(ValueError):
        mol = toolkit.Mol.from_smiles('ZZZ')


@pytest.mark.parametrize('toolkit', mts)
def test_bond(toolkit):
    mol = toolkit.Mol.from_smiles('C')
    print('made molecule')
    bond = mol.get_bond_by_index(0)

    assert bond.get_order() == 1

    atoms = bond.get_atoms()
    assert len(atoms) == 2

    assert not bond.is_ring()

    assert not bond.is_aromatic()

    assert bond.is_single()

    assert not bond.is_double()

    assert not bond.is_triple()

    mol = bond.get_molecule()
    smiles = mol.get_smiles()
    assert smiles == "C"

    print('trying to get index')
    assert bond.get_index() == 0
    print('past bond index')


@pytest.mark.parametrize('toolkit', mts)
def test_atom(toolkit):
    mol = toolkit.Mol.from_smiles('C')
    atom = mol.get_atom_by_index(0)

    assert atom.atomic_number() == 6

    assert atom.degree() == 4

    assert atom.connectivity() == 4

    assert atom.valence() == 4

    assert atom.formal_charge() == 0

    assert atom.hydrogen_count() == 4

    assert atom.ring_connectivity() == 0

    assert atom.min_ring_size() == 0

    assert not atom.is_aromatic()

    assert atom.get_index() == 0

    neighbors = atom.get_neighbors()
    assert len(neighbors) == 4

    atom2 = neighbors[0]
    assert atom.is_connected_to(atom2)

    assert len(atom.get_bonds()) == 4

    # at least run get_molecule function, not sure how to check this
    mol = atom.get_molecule()
    smiles = mol.get_smiles()
    assert smiles == "C"


@pytest.mark.parametrize('toolkit', mts)
def test_atom_exception(toolkit):
    with pytest.raises(TypeError):
        toolkit.Atom(None)


@pytest.mark.parametrize('toolkit', mts)
def test_bond_exception(toolkit):
    with pytest.raises(TypeError):
        toolkit.Bond(None)


@pytest.mark.parametrize('toolkit', mts)
def test_mol_exception(toolkit):
    with pytest.raises(TypeError):
        toolkit.Mol(None)


# -------------------------------
# check molecule file parsers
# -------------------------------

mol2_abs_file = chemper_utils.get_data_path('molecules/MiniDrugBank_tripos.mol2')
mol2_rel_path = 'MiniDrugBank_tripos.mol2'
paths = [mol2_abs_file, mol2_rel_path]


# For the following functions, we will test default behavior and
# look for exceptions based on the available mol toolkit

@pytest.mark.parametrize('toolkit,path', itertools.product(mts, paths))
def test_file_parsing(toolkit,path):
    mols = toolkit.mols_from_mol2(path)
    assert len(mols) == 363


@pytest.mark.parametrize('path', paths)
def test_mols_specified_toolkit(path):
    if mol_toolkit.HAS_OE:
        mols = mol_toolkit.mols_from_mol2(path, toolkit='openeye')
        assert len(mols) == 363
    else:
        with pytest.raises(ImportError):
            mols = mol_toolkit.mols_from_mol2(path, toolkit='openeye')

    if mol_toolkit.HAS_RDK:
        mols = mol_toolkit.mols_from_mol2(path, toolkit='rdkit')
        assert len(mols) == 363
    else:
        with pytest.raises(ImportError):
            mols = mol_toolkit.mols_from_mol2(path, toolkit='rdkit')


# -------------------------------
# check functions for checking
# -------------------------------
tk_fails = ['babel']
if not mol_toolkit.HAS_OE:
    tk_fails.append('openeye')
if not mol_toolkit.HAS_RDK:
    tk_fails.append('rdkit')


@pytest.mark.parametrize('tks', tk_fails)
def test_fake_toolkit(tks):
    with pytest.raises(ImportError):
        mol_toolkit.check_toolkit(tks)


def test_default_toolkit():
    tk = mol_toolkit.check_toolkit(None)
    if mol_toolkit.HAS_OE:
        assert tk == 'openeye'
    else:
        assert tk == 'rdkit'
