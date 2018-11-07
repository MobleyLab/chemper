"""
type_minidrugbank.py

Author: Caitlin C. Bannan U.C. Irvine Mobley Group

This script shows an example of chemper behavior by attempting to generate
SMIRKS patterns to mirror those in smirnoff99Frosst.

The final script will take the following steps

1. Load molecules
2. For each type of SMIRKS (bonds, angles, proper torsions, non-bonds):

    a. Type molecules with the SMIRKS patterns from smirnoff99Frosst version 1.0.7.
       These SMIRKS patterns are available at chemper/data/smarts_files.
    b. Use the clusters based on smirnoff99Frosst to generate custom SMIRKS with chemper's SMIRKSifier
    c. Test that the generate SMIRKS and the SMIRKS from smirnoff99Frosst both type MiniDrugBank in the same way
    d. store these SMIRKS

5. After it has been created I'll add an additional test on outside molecules
(Zinc set or more of DrugBank)
"""

from copy import deepcopy
import time
from chemper.chemper_utils import get_full_path, create_tuples_for_clusters
from chemper.mol_toolkits.mol_toolkit import mols_from_mol2, Mol
from chemper.smirksify import SMIRKSifier, print_smirks

def parse_smarts_file(file_path):
    """

    Parameters
    ----------
    file_path: str
        relative path in chemper/data or absolute path

    Returns
    -------
    smirks_list: list of tuples (label, smirks)
        This is the ordered list of SMIRKS from the smarts file
        if a label is provided in the file then it is assigned, otherwise
        the indices from the file is used
    """
    fn = get_full_path(file_path)
    f = open(fn)
    lines = f.readlines()
    f.close()

    lines = [l.split() for l in lines]
    if len(lines[0]) > 1:
        return [(label, smirks) for smirks, label in lines]

    return [(str(i), l[0]) for i, l in enumerate(lines)]


def make_smarts_file(smirks_list, output_file):
    """
    Parameters
    ----------
    smirks_list: list of tuples
        smirks tuples have the form (label, smirks
    output_file: relative or absolute path
    """
    f = open(output_file, 'w')
    for label, smirks in smirks_list:
        f.write("%s %s\n" % (smirks, label))
    f.close()


def run_fragement(frag, mols, steps=1000):
    """

    Parameters
    ----------
    frag: str
        angle, bond, nonbond, proper_torsion
    mols: list of molecules
    steps: number of reducing steps

    Returns
    -------
    new_smirks_list: list of tuples
        Final list of tuples
    """
    print("Searching %s" % frag)

    # open and parse SMIRKS file
    smirks_list = parse_smarts_file('smarts_files/%s_smirks.smarts' % frag)

    # copy the molecules:
    current_mols = deepcopy(mols)

    # type molecules with these smirks:
    print("clustering molecules:")
    init = time.time()
    clusters = create_tuples_for_clusters(smirks_list, current_mols)
    end = time.time()
    print("Took %.3f minutes to type %i mols\n" % ((end-init)/60., len(mols)))

    # Try to make SMIRKS for these clusters:
    print("Creating initial SMIRKS...")
    init = time.time()
    create_frag = SMIRKSifier(current_mols, clusters, verbose=False)
    end = time.time()
    print(end, init, type(end), type(init))
    print("generating smirksifier took %.3f minutes\n" % ((end-init)/60.))

    print_smirks(create_frag.current_smirks)

    print("Reducing SMIRKS...")
    init = time.time()
    new_smirks_list = create_frag.reduce(steps)
    end = time.time()
    print("reducing SMIRKS took %.3f minutes\n" % ((end-init)/60.))

    return new_smirks_list


# 1. load molecules
mol_file = 'MiniDrugBank_tripos.mol2'
mol_file = '/Users/bannanc/Desktop/baby_minidrugbank.mol2'
mol_file = "/Users/bannanc/gitHub/openforcefield/openforcefield/data/molecules/AlkEthOH_test_filt1_tripos.mol2"
mols = list()
for m in mols_from_mol2(mol_file):
    new_m = Mol(m)
    m.set_aromaticity_mdl()
    mols.append(m)

# Loop over all fragment types
fragments = ['nonbond', 'bond', 'angle', 'torsion']

for fragment in fragments:
    try:
        frag_smirks_list = run_fragement(fragment, mols, steps=10)
        print_smirks(frag_smirks_list)
        frag_file = "%s_output_smirks.smarts" % fragment
        make_smarts_file(frag_smirks_list, frag_file)
    except:
        print("unable to make smirks for %s" % fragment)