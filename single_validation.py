#!/usr/bin/python3

import numpy as np
import pandas as pd
from d_utils import runner_utils as ru, config, eval_utils as eu
from evaluation import comparator as comp
import random
from mappers.mapper import Mapper, FileMapper


def single_validation(tar, tar_id, mode, ref=None, ref_id=None, enriched: bool = False,
                      mapper: Mapper = FileMapper(), out_dir=""):
    """
    Single validation of a set, cluster a id versus set and set versus set.

    :param tar:
    :param tar_id:
    :param mode:
    :param ref:
    :param ref_id:
    :param enriched:
    :param mapper:
    :param out_dir:
    :return:
    """
    # ===== Comparison with a set =====
    ru.print_current_usage('Starting validation ...')
    ru.print_current_usage('Load mappings for input into cache ...')
    mapper.load_mappings()

    if mode in ["set", "set-set", "id-set"]:
        if mode == "set-set":
            comparator = comp.SetSetComparator(mapper=mapper, enriched=enriched)
            comparator.load_reference(ref=pd.read_csv(ref, header=None, sep="\t")[0], ref_id_type=ref_id)
        elif mode == "id-set":
            comparator = comp.IDSetComparator(mapper=mapper)
            comparator.load_reference(ref=ref, ref_id_type=ref_id)
        else:  # mode == "set"
            comparator = comp.SetComparator(mapper=mapper)
            ru.print_current_usage('Load distances for input into cache ...')
            mapper.load_distances(set_type=tar_id)
        comparator.load_target(id_set=pd.read_csv(tar, header=None, sep="\t")[0], id_type=tar_id)

        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...')
        my_value = comparator.compare()
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(comparator=comparator, mode=mode, mapper=mapper, tar_id=tar_id)
        ru.print_current_usage('Calculating p-values ...')
        if mode == "set":
            p_values = eu.calc_pvalue(test_value=my_value, value_df=pd.DataFrame(comp_values), maximize=False)
        else:
            p_values = eu.calc_pvalue(test_value=my_value, value_df=pd.DataFrame(comp_values))

    # ===== Special case cluster =====
    elif args.mode == "cluster":
        ru.print_current_usage('Load distances for input into cache ...')
        mapper.load_distances(set_type=tar_id)
        ru.print_current_usage('Load input data ...')
        comparator = comp.ClusterComparator(mapper=mapper)
        comparator.load_target(id_set=pd.read_csv(tar, header=None, sep="\t"), id_type=tar_id)
        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...')
        my_value_di, my_value_ss = comparator.compare()
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(comparator=comparator, mode=mode, mapper=mapper, tar_id=tar_id)
        p_values_di = eu.calc_pvalue(test_value=my_value_di, value_df=pd.DataFrame(comp_values[0]))
        p_values_ss = eu.calc_pvalue(test_value=my_value_ss, value_df=pd.DataFrame(comp_values[1]))
        p_values = {'di': p_values_di, 'ss': p_values_ss}
    else:
        return None
    # mapper.save_mappings()
    # mapper.save_distances()
    ru.print_current_usage('Finished validation')
    print(p_values)


def get_random_runs_values(comparator: comp.Comparator, mode, mapper: Mapper, tar_id: str) -> list:
    # ===== On the fly =====
    results = list()
    if not mode == "cluster":
        if tar_id in config.SUPPORTED_DISEASE_IDS:
            full_id_list = mapper.get_full_set(id_type=tar_id, mapping_name='disorder_ids')
        else:  # if tar_id in config.SUPPORTED_GENE_IDS:
            full_id_list = mapper.get_full_set(id_type=config.ID_TYPE_KEY[tar_id], mapping_name='gene_ids')
        # ===== Remove empty values =====
        full_id_list = list(filter(None, full_id_list))
        # ===== Calculate values =====
        for run in range(0, config.NUMBER_OF_RANDOM_RUNS):
            comparator.load_target(id_set=random.sample(full_id_list, len(comparator.id_set)), id_type=tar_id)
            results.append(comparator.compare())

    # ===== Special case cluster =====
    else:
        # ===== Calculate values =====
        results.extend((list(), list()))
        for run in range(0, config.NUMBER_OF_RANDOM_RUNS):
            comparator.clustering['cluster_index'] = np.random.permutation(comparator.clustering['cluster_index'])
            value_di, value_ss = comparator.compare()
            results[0].append(value_di)
            results[1].append(value_ss)

    return results


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    args = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o', 'e'))
    single_validation(tar=args.target, tar_id=args.target_id_type,
                      mode=args.mode, ref=args.reference, ref_id=args.reference_id_type,
                      enriched=args.enriched, out_dir=args.out_dir)
