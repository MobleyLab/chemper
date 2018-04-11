"""
This is a general test for importing the tool for now
"""

from chemper.mol_toolkits import mol_toolkit
print('testing mol_toolkits %s' % mol_toolkit.__name__)
import pytest

def test_molecule():
    """
    Test MolOE functions
    """
    mol = mol_toolkit.MolFromSmiles('C')

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


def test_smirks_search():
    """
    test SMIRKS searching
    """
    mol = mol_toolkit.MolFromSmiles('C')

    # smirks for C-H bond
    smirks = "[#6:1]-[#1:2]"

    matches = mol.smirks_search(smirks)

    assert len(matches) == 4

    for match in matches:
        assert 1 in match
        assert 2 in match

def test_bad_smirks():
    """
    Check a ValueError is raised with improper SMIRKS
    """
    mol = mol_toolkit.MolFromSmiles('C')
    with pytest.raises(ValueError):
        mol.smirks_search(']X[')

def test_bad_smiles():
    """
    Check a ValueError is raised with a bad SMILES
    """
    with pytest.raises(ValueError):
        mol = mol_toolkit.MolFromSmiles('ZZZ')

def test_bond():
    mol = mol_toolkit.MolFromSmiles('C')
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


def test_atom():
    mol = mol_toolkit.MolFromSmiles('C')
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

