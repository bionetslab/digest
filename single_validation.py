#!/usr/bin/python3

import pandas as pd
from d_utils import runner_utils as ru, config, eval_utils as eu
from evaluation import comparator as comp
import random
import numpy as np
from mappers.mapper import Mapper, FileMapper


def single_validation(tar, tar_id, mode, ref=None, ref_id=None, mapper: Mapper = FileMapper(), out_dir=""):
    """
    Single validation of a set, cluster a id versus set and set versus set.

    :param tar:
    :param tar_id:
    :param mode:
    :param ref:
    :param ref_id:
    :param mapper:
    :param out_dir:
    :return:
    """
    reference, target = None, None
    # ===== Comparison with a set =====
    ru.print_current_usage('Starting validation ...')
    ru.print_current_usage('Load mappings for input into cache ...')
    mapper.load_mappings(set_type=ref_id)
    mapper.load_mappings(set_type=tar_id)
    if mode in ["set", "set-set", "id-set"]:
        target = pd.read_csv(tar, header=None, sep="\t")[0]
        if mode == "set-set":
            reference = pd.read_csv(ref, header=None, sep="\t")[0]
        elif mode == "id-set":
            reference = ref
        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...')
        my_value = get_validation(ref=reference, ref_id=ref_id,
                                  tar=target, tar_id=tar_id, mode=mode, mapper=mapper)
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(ref=reference, ref_id=ref_id,
                                             tar=target, tar_id=tar_id, mode=mode, mapper=mapper)
        ru.print_current_usage('Calculating p-values ...')
        if mode == "set":
            p_values = eu.calc_pvalue(test_value=my_value, value_df=pd.DataFrame(comp_values), maximize=False)
        else:
            p_values = eu.calc_pvalue(test_value=my_value, value_df=pd.DataFrame(comp_values))

    # ===== Special case cluster =====
    elif args.mode == "cluster":
        target = pd.read_csv(tar, header=None, sep="\t")
        # ===== Get validation values of input =====
        my_value_di, my_value_ss = get_validation(ref=reference, ref_id=ref_id,
                                                  tar=target, tar_id=tar_id, mode=mode, mapper=mapper)
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(ref=reference, ref_id=ref_id,
                                             tar=target, tar_id=tar_id, mode=mode, mapper=mapper)
        p_values_di = eu.calc_pvalue(test_value=my_value_di, value_df=pd.DataFrame(comp_values[0]))
        p_values_ss = eu.calc_pvalue(test_value=my_value_ss, value_df=pd.DataFrame(comp_values[1]))
        p_values = {'di': p_values_di, 'ss': p_values_ss}
    else:
        return None
    mapper.save_mappings()
    ru.print_current_usage('Finished validation')
    print(p_values)


def get_validation(ref, ref_id, tar, tar_id, mode, mapper: Mapper, verbose=True):
    if mode == "set":
        return comp.compare_set(id_set=tar, id_type=tar_id, mapper=mapper, verbose=verbose)
    elif mode == "set-set":
        return comp.compare_set_to_set(ref=ref, ref_id_type=ref_id, tar=tar, tar_id_type=tar_id, mapper=mapper)
    elif mode == "id-set":
        return comp.compare_id_to_set(ref_id=ref, ref_id_type=ref_id, tar=tar, tar_id_type=tar_id, mapper=mapper)
    elif mode == "cluster":
        return comp.compare_clusters(clusters=tar, id_type=tar_id, mapper=mapper, verbose=verbose)
    else:
        return None


def get_random_runs_values(ref, ref_id, tar, tar_id, mode, mapper: Mapper) -> list:
    # ===== On the fly =====
    results = list()
    if not mode == "cluster":
        if tar_id in config.SUPPORTED_DISEASE_IDS:
            full_id_list = mapper.get_full_set(id_type=config.ID_TYPE_KEY[tar_id], mapping_name='disorder_ids')
        else:  # if tar_id in config.SUPPORTED_GENE_IDS:
            full_id_list = mapper.get_full_set(id_type=config.ID_TYPE_KEY[tar_id], mapping_name='gene_ids')
        # ===== Remove empty values =====
        full_id_list = list(filter(None, full_id_list))
        # ===== Calculate values =====
        for run in range(0, config.NUMBER_OF_RANDOM_RUNS):
            results.append(get_validation(ref=ref, ref_id=ref_id, tar=random.sample(full_id_list, len(tar)),
                                          tar_id=tar_id, mode=mode, mapper=mapper, verbose=False))
    else:
        # ===== Calculate values =====
        results.extend((list(), list()))
        for run in range(0, config.NUMBER_OF_RANDOM_RUNS):
            tar[1] = np.random.permutation(tar[1].values)
            value_di, value_ss = get_validation(ref=ref, ref_id=ref_id, tar=tar, tar_id=tar_id, mode=mode,
                                                mapper=mapper, verbose=False)
            results[0].append(value_di)
            results[1].append(value_ss)

    return results


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    args = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o'))
    single_validation(tar=args.target, tar_id=args.target_id_type,
                      mode=args.mode, ref=args.reference, ref_id=args.reference_id_type, out_dir=args.out_dir)
