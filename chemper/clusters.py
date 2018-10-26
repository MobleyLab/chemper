"""
clusters.py

This script contains functions for clustering data
and uses the SMIRKSifier class to then generate SMIRKS patterns for
those clusters.
"""

import numpy as np
from collections import OrderedDict
from chemper.mol_toolkits.mol_toolkit import Mol
from chemper.smirksify import SMIRKSifier
from chemper.chemper_utils import get_typed_molecules


def cluster_fragments(data_dicts, cluster_class=None):
    """

    Parameters
    ----------
    data_dicts: list of dictionaries
        This is a list of dictionaries for data for each input molecule
        For example, if you were looking to cluster on two pieces of
        data for a bond, you would provide a dictionary in the form
        { (0,1): (x,y), (1,2): (x2, y2) ..} for each molecule being
        considered
    cluster_class: (optional) clustering function
        This is a tool for clustering the provided data,
        such as a prepped BayesianGaussianMixture from SciKit-Learn.
        This class should have the functions .fit and .predict
        If this is set to None, then a BayesianGaussianMixture with
        default settings will be used.

    Returns
    -------
    assignments: list of tuples
        These are set up to be the input for SMIRKSifier
        So for each cluster you would have a list of atom indices
        from each molecule. For example, if the clustered
        ended up putting the bonds (0,1) and (1,2) for two molecules into
        cluster 'c1' and the bond (2,3) into a second cluster the
        output would look like this:

        [ ('c1', [ (0,1), (1,2) ], [ (0,1), (1,2) ]),
          ('c2', [ (2,3)        ], [ (2,3)        ]) ]
    """
    # organize the data into an array that can be fit
    # store the atom indices in a list of the same shape
    # to be referenced later
    atom_list = list()
    data_list = list()
    num_mols = len(data_dicts)

    for mol_idx, mol_dict in enumerate(data_dicts):
        for atoms, data in mol_dict.items():
            atom_list.append((mol_idx, atoms))
            data_list.append(data)

    x = np.array(data_list)
    if cluster_class is None:
        from sklearn.mixture import BayesianGaussianMixture
        n_components = min(x.shape[0], 20)
        cluster_class = BayesianGaussianMixture(n_components=n_components,
                                                covariance_type='full',
                                                max_iter=100)
    # fit data
    gmm = cluster_class.fit(x)

    # organize output
    output = dict()
    for s, mol_info in zip(gmm.predict(x), atom_list):
        if s not in output:
            output[s] = [list() for i in range(num_mols)]
        output[s][mol_info[0]].append(mol_info[1])

    output = [(k, e) for k, e in output.items()]
    return sorted(output, key=lambda x: sum([len(i) for i in x[1]]), reverse=True)


def cluster_fragments_to_smirks(data_dicts, mols,
                                test_smirks=None, max_its=2000,
                                cluster_class=None):
    """
    Parameters
    ----------
    data_dicts: list of dictionaries
        This is a list of dictionaries for data for each input molecule
        For example, if you were looking to cluster on two pieces of
        data for a bond, you would provide a dictionary in the form
        { (0,1): (x,y), (1,2): (x2, y2) ..} for each molecule being
        considered
    mols: list of molecules
        These can be from any toolkit supported by chemper,
        OpenEye, RDKit, or ChemPer
    test_smirks: (optional) list of SMIRKS
        default = None
        This is a list of SMIRKS in the form (label, smirks)
        if the clustered fragments and these SMIRKS agree, then
        this function will just return how those SMIRKS type a list of molecules
    max_its: (optional) maximum iterations
        default = 2,000
        How many iterations to use while reducing molecules
    cluster_class: (optional) clustering function
        This is a tool for clustering the provided data,
        such as a prepped BayesianGaussianMixture from SciKit-Learn.
        This class should have the functions .fit and .predict
        If this is set to None, then a BayesianGaussianMixture with
        default settings will be used.
    """
    clusters = cluster_fragments(data_dicts, cluster_class)
    red = SMIRKSifier(mols, clusters, max_layers=10, verbose=False)
    matches = False
    if test_smirks is not None:
        _, matches = red.types_match_reference(current_types=test_smirks)

    if matches:
        smirks_list = test_smirks
    else:
        smirks_list = red.reduce(max_its)

    # return OrderedDict{ smirks: [list of tuples for each molecule]
    d = OrderedDict()
    lab_to_smirks = dict()
    for lab, smirks in smirks_list:
        d[smirks] = [list() for i in range(len(mols))]
        lab_to_smirks[lab] = smirks

    typed_mols = get_typed_molecules(smirks_list, mols)
    for mol_idx, match_dict in typed_mols.items():
        for atom_indices, label in match_dict.items():
            smirks = lab_to_smirks[label]
            d[smirks][mol_idx].append(atom_indices)

    return d
