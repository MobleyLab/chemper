"""
Tests for the script chemper.clusters
"""
import pytest
from chemper.clusters import cluster_fragments, cluster_fragments_to_smirks
from chemper.mol_toolkits.mol_toolkit import MolFromSmiles

data_dicts2 = [{
    (0, 1): (0, 0.001412823150423126,), (0, 2): (0, -0.42223584033038275,),
    (0, 3): (0, -0.42221759356346955,), (0, 4): (0, -0.4222139322856049,),
    (1, 5): (0, -0.4222363689959121,), (1, 6): (0, -0.422217217432993,),
    (1, 7): (0, -0.4222137474302068,)}
]
data_dicts1 = [
    {atoms: (d[1],) for atoms, d in data_dicts2[0].items()}
]

list_data_dicts = [ data_dicts1, data_dicts2 ]


@pytest.mark.parametrize('data_dicts',list_data_dicts)
def test_cluster_fragments(data_dicts):
    output = cluster_fragments(data_dicts)
    assert len(output) == 2
    o1 = output[0][1][0]
    assert len(o1) == 6
    for atoms in [(1, 5), (1, 6), (0, 4), (0, 3), (1, 7), (0, 2)]:
        assert atoms in o1

    o2 = output[1][1][0]
    assert len(o2) == 1
    assert o2[0] == (0,1)


smirks_input = [
    None,
    [('c1', '[*:1]~[*:2]')],
    [('c1', '[*:1]~[*:2]'), ('c2', '[#1:1]~[#6:2]'), ('c1', '[#6:1]~[#6:2]')]
]
test_list = [
    (s, d) for s in smirks_input for d in list_data_dicts
]


@pytest.mark.parametrize('test_smirks, data_dicts', test_list)
def test_cluster_fragments_to_smirks(test_smirks, data_dicts):
    mol = MolFromSmiles('CC')
    d = cluster_fragments_to_smirks(data_dicts, [mol], test_smirks=test_smirks, max_its=100)
    output = dict()
    i = 0
    for smirks, mol_info in d.items():
        output[i] = mol_info[0]
        i += 1
    assert len(output[0]) == 6
    assert output[1][0] == (0,1)

@pytest.mark.parametrize('data_dicts', list_data_dicts)
def test_give_cluster(data_dicts):
    mol = MolFromSmiles('CC')
    ts = [('c1', '[*:1]~[*:2]'), ('c2', '[*X4:1]~[*X4:2]')]
    d = cluster_fragments_to_smirks(data_dicts, [mol], test_smirks=ts)
    for l, s in ts:
        assert s in d
