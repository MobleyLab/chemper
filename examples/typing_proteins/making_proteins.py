"""
making_proteins.py


"""

import os
import parmed
from simtk.openmm import app
from simtk import unit
from oeommtools import utils as oeo_utils
from chemper.smirksify import SMIRKSifier
from chemper.graphs.cluster_graph import ClusterGraph
import copy
import random
from openeye import oechem


class ParameterDict:

    def __init__(self):
        self.d = dict()

    def items(self):
        return self.d.items()

    def add_key(self, key):
        if key not in self.d:
            self.d[key] = {'atom_indices': dict(), 'parameters': set(), 'units': None}

    def add_atoms(self, key, mol_id, atom_tuple):
        self.add_key(key)
        if mol_id not in self.d[key]['atom_indices']:
            self.d[key]['atom_indices'][mol_id] = list()
        self.d[key]['atom_indices'][mol_id].append(tuple(atom_tuple))

    def add_param(self, key, params):
        self.add_key(key)
        new_tuple = [x._value for x in params]
        self.d[key]['parameters'].add(tuple(new_tuple))
        self.d[key]['units'] = tuple([x.unit for x in params])


class ParameterSystem:

    def __init__(self, openmm_xml='amber14-all.xml'):
        self.openmm_xml = openmm_xml
        self.lj_dict = ParameterDict()
        self.charge_dict = ParameterDict()
        self.bond_dict = ParameterDict()
        self.angle_dict = ParameterDict()
        self.proper_dict = ParameterDict()
        self.improper_dict = ParameterDict()
        self.mol_dict = dict()

    def add_system_from_fasta(self, fasta):
        base = os.path.abspath(fasta).split('.')[0]
        mol_id = base.split('/')[-1]

        oemol = oechem.OEMol()
        ifs = oechem.oemolistream(fasta)

        oechem.OEReadFASTAFile(ifs, oemol)
        oechem.OEAddExplicitHydrogens(oemol)

        oechem.OEPerceiveResidues(oemol)
        oechem.OEPDBOrderAtoms(oemol)

        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat_PDB)
        ofs.openstring()
        oechem.OEWriteMolecule(ofs, oemol)

        ifs = oechem.oemolistream()
        ifs.openstring(ofs.GetString())
        m = oechem.OEMol()
        oechem.OEReadPDBFile(ifs, m)
        m.SetTitle(mol_id)
        top = oeo_utils.oemol_to_openmmTop(m)[0]
        ff = app.ForceField(self.openmm_xml)
        protein_sys = ff.createSystem(top)
        parm = parmed.openmm.load_topology(top, protein_sys)

        oechem.OEAssignFormalCharges(m)
        self.mol_dict[mol_id] = {
            'parmed': parm,
            'oemol': oechem.OEMol(m)
        }
        self._add_parameters_from_system(parm, mol_id)

        return parm, m

    def _add_parameters_from_system(self, sys, mol_id):
        self.add_nonbonds(sys, mol_id)
        self.add_bonds(sys, mol_id)
        self.add_angles(sys, mol_id)
        self.add_torsions(sys, mol_id)

    def add_nonbonds(self, sys, mol_id):
        """
        Cluster atoms based on their partial charge

        Parameters
        ----------
        sys: list like of parmed system
        charge_dict: dictionary to store data that will be updated in this function
        lj_dict: dictionary to store LJ parameters for this molecule
        mol_id: key for this system to store data in the dictionaries

        Returns
        -------
        clusters: dictionary with the form
                  {string parameter: {'atom_idices': {}}
        """
        # TODO raise error if mol_id not in mol_dict
        for a in sys.atoms:
            # Update charge dictionary:
            charge_str = "%.5f" % a.charge
            charge_param = [a.ucharge]
            self.charge_dict.add_param(charge_str, charge_param)
            self.charge_dict.add_atoms(charge_str, mol_id, [a.idx])

            # Update LJ dictionary
            lj_str = "%.3f\t%.3f" % (a.epsilon, a.rmin)
            lj_params = [a.uepsilon, a.urmin]
            self.lj_dict.add_param(lj_str, lj_params)
            self.lj_dict.add_atoms(lj_str, mol_id, [a.idx])

    def add_bonds(self, sys, mol_id):
        # TODO raise error if mol_id not in mol_dict
        for b in sys.bonds:
            bond_str = "%.3f\t%.3f" % (b.type.k, b.type.req)
            bond_params = [b.type.uk, b.type.ureq]
            self.bond_dict.add_param(bond_str, bond_params)
            self.bond_dict.add_atoms(bond_str, mol_id, [b.atom1.idx, b.atom2.idx])

    def add_angles(self, sys, mol_id):
        # TODO raise error if mol_id not in mol_dict
        for an in sys.angles:
            angle_str = "%.3f\t%.3f" % (an.type.k, an.type.theteq)
            angle_params = [an.type.uk, an.type.utheteq]
            self.angle_dict.add_param(angle_str, angle_params)
            self.angle_dict.add_atoms(angle_str, mol_id, [an.atom1.idx, an.atom2.idx, an.atom3.idx])

    def convert_for_smirksifying(self, param_type=None):
        """
        param_type: string specifying the parameter you want clusters for
        must chose from ['lj', 'charge', 'proper_torsion', 'improper_torsion', 'angle', 'bond']

        Returns
        -------
        - list of molecules
        - either dictionary or list of clustered atomic indices

        """
        idx_list = list()
        mol_list = list()
        cluster_types = dict()

        dictionaries = {
            'lj': self.lj_dict,
            'charge': self.charge_dict,
            'proper_torsion': self.proper_dict,
            'improper_torsion': self.improper_dict,
            'angle': self.angle_dict,
            'bond': self.bond_dict,
        }

        if param_type is not None:
            if param_type.lower() not in dictionaries.keys():
                return cluster_types
            dictionaries = {param_type.lower(): dictionaries[param_type.lower()]}

        for idx, me in self.mol_dict.items():
            idx_list.append(idx)
            mol_list.append(me['oemol'])

        for label, par_dict in dictionaries.items():
            cluster_types[label] = list()
            for cluster_label, entry in par_dict.items():
                atom_list = list()
                for idx in idx_list:
                    if idx in entry['atom_indices']:
                        atom_list.append(entry['atom_indices'][idx])
                    else:
                        atom_list.append(list())
                cluster_types[label].append((cluster_label, atom_list))

        if param_type is None:
            return mol_list, cluster_types

        return mol_list, cluster_types[param_type.lower()]

    def add_torsions(self, sys, mol_id):
        # TODO raise error if mol_id not in mol_dict
        temp_dict = dict()
        for d in sys.dihedrals:
            if d.improper:
                imp_str = "%.3f\t%.3f\t%.3f" % (d.type.phi_k, d.type.phase, d.type.per)
                imp_params = [d.type.uphi_k, d.type.uphase, unit.Quantity(d.type.per)]
                self.improper_dict.add_param(imp_str, imp_params)
                self.improper_dict.add_atoms(imp_str, mol_id, [d.atom1.idx, d.atom3.idx, d.atom2.idx, d.atom4.idx])
            else:
                atoms = tuple([d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx])
                params = (d.type.uphi_k, d.type.uphase, unit.Quantity(d.type.per))
                if atoms not in temp_dict:
                    temp_dict[atoms] = list()
                temp_dict[atoms].append(params)

        for atoms, param_list in temp_dict.items():
            new_params = [p for t in param_list for p in t]
            prop_str = '\t'.join(['%.3f' % p._value for p in new_params])
            self.proper_dict.add_param(prop_str, new_params)
            self.proper_dict.add_atoms(prop_str, mol_id, atoms)


def change_order_smirksified(mols, cluster_types, order_type_names=None, smirks_verbose=False,
                             include_params=None):
    """
    Returns smirs_order_types object
    """

    def reverse_clusters(clusters):
        return list(reversed(clusters))

    def shuffle(clusters):
        temp_c = copy.deepcopy(clusters)
        random.shuffle(temp_c)
        return temp_c

    def by_smallest_size(clusters):
        return sorted(clusters, key=lambda x: len([a for l in x[1] for a in l]))

    def by_smallest_num_molecule(clusters):
        temp_c = sorted(clusters, key=lambda x: len([1 for l in x[1] if len(l) > 0]))
        return temp_c

    def by_biggest_size(clusters):
        return reverse_clusters(by_smallest_size(clusters))

    def by_biggest_num_molecule(clusters):
        return reverse_clusters(by_smallest_num_molecule(clusters))

    def by_smallest_smirks(clusters):
        temp_c = sorted(clusters, key=lambda x: len(ClusterGraph(mols, x[1]).as_smirks()))
        return temp_c

    def by_biggest_smirks(clusters):
        return reverse_clusters(by_smallest_smirks(clusters))

    smirs_order_types = dict()

    order_types_dict = {
        'original': None,
        'reversed': reverse_clusters,
        'shuffle': shuffle,
        'small_size': by_smallest_size,
        'biggest_size': by_biggest_size,
        'fewest_mols': by_smallest_num_molecule,
        'most_mols': by_biggest_num_molecule,
        'small_smirks': by_smallest_smirks,
        'big_smirks': by_biggest_smirks}

    if order_type_names is None:
        order_type_names = ['shuffle']

    order_types = list()
    for n in set(order_type_names):
        o_funct = order_types_dict.get(n, None)
        order_types.append((n, o_funct))
        for i in range(1, order_type_names.count(n)):
            temp_n = '%s_%i' % (n, i)
            order_types.append((temp_n, o_funct))

    if include_params is None:
        include_params = list(cluster_types.keys())

    for o_type, o_funct in order_types:
        print(o_type)
        smirs_order_types[o_type] = dict()

        for label, clusters in cluster_types.items():
            if label not in include_params:
                continue

            if o_funct is None:
                o_clusters = clusters
            else:
                o_clusters = o_funct(clusters)

            smirs_order_types[o_type][label] = SMIRKSifier(mols, o_clusters,
                                                           max_layers=10,
                                                           strict_smirks=False,
                                                           verbose=smirks_verbose)

    return smirs_order_types


def print_order_type_data(smirks_order_types, print_all=True):
    for o_type, order_smirks in smirks_order_types.items():
        final_print = ''
        all_passed = True
        for label, output in order_smirks.items():
            if not output.checks:
                final_print += '%-23s FAILED to make smirks\n' % label
                all_passed = False
            else:
                final_print += '%-23s PASSED\n' % label

        if all_passed:
            print('-' * 80)
            print('%-23s ALL PASSED' % o_type)
            print('-' * 80)
            print(final_print)
        elif print_all:
            print('-' * 80)
            print(o_type)
            print('-' * 80)
            print(final_print)


def at_least_one_passed(smirks_order_types):
    passed = True
    for o_type, order_smirks in smirks_order_types.items():
        passed = True
        for label, output in order_smirks.items():
            if not output.checks:
                passed = False
        if passed:
            return True
    return passed


def everything_from_fastas(list_fastas,
                           protein_xml='amber14-all.xml',
                           order_type_names=None,
                           verbose=True,
                           include_params=None):
    """
    Parameters
    ----------
    list_fastas: list of str
                 list of full paths to .fasta files with peptide sequences

    protein_xml: str
                 file name or complete path to a openMM .xml file for assign protein parameters
    order_type_names: list of str

    Returns
    -------
    data_storage: ParameterSystem object
                  ParameterSystem storing the parameterized molecules and clusters
                  by parameter type from the provided fasta files
    smirs_order_types: dictionary with the form...
            {
            order type (ie 'shuffle' or 'smallest_size): {
                param_type (ie 'proper_torsion'): SMIRKSifier object (even if it failed)
                }
            }
    """
    if order_type_names is None:
        order_type_names = ['shuffle']

    store_data = ParameterSystem(openmm_xml=protein_xml)
    for fasta in list_fastas:
        parm, oemol = store_data.add_system_from_fasta(fasta)
    mols, cluster_types = store_data.convert_for_smirksifying()

    if verbose:
        table_form = "%-20s %-10s %-10s %s"
        print('=' * 80)
        print(table_form % ('parameter', 'mols', 'clusters', 'mols in clusters'))
        print('-' * 80)
        for label, clusters in cluster_types.items():
            print(table_form % (label, len(mols), len(clusters), len(clusters[0][1])))
        print('=' * 80)

    smirs_order_types = change_order_smirksified(mols, cluster_types,
                                                 order_type_names=order_type_names,
                                                 include_params=include_params)
    if verbose: print_order_type_data(smirs_order_types)

    return store_data, smirs_order_types, mols, cluster_types

