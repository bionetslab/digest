#!/usr/bin/python3

import pandas as pd
from d_utils import runner_utils as ru, config, mapping_utils as mu, eval_utils as eu
from evaluation import comparator as comp
import random
import numpy as np
from mappers.mapper import Mapper, FileMapper


def single_validation(args, mapper: Mapper = FileMapper()):
    """
    Single validation of a set, cluster a id versus set and set versus set.

    :param mapper:
    :param args: user given parameters
    :return:
    """
    reference, target = None, None
    # ===== Comparison with a set =====
    ru.print_current_usage('Starting validation ...')
    ru.print_current_usage('Load mappings for input into cache ...')
    mapper.load_mappings(set_type=args.reference_id_type)
    mapper.load_mappings(set_type=args.target_id_type)
    if args.mode in ["set", "set-set", "id-set"]:
        target = pd.read_csv(args.target, header=None, sep="\t")[0]
        if args.mode == "set-set":
            reference = pd.read_csv(args.reference, header=None, sep="\t")[0]
        elif args.mode == "id-set":
            reference = args.reference
        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...')
        my_value = get_validation(ref=reference, ref_id=args.reference_id_type,
                                  tar=target, tar_id=args.target_id_type, mode=args.mode, mapper=mapper)
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(ref=reference, ref_id=args.reference_id_type,
                                             tar=target, tar_id=args.target_id_type, mode=args.mode, mapper=mapper)
        p_values = eu.calc_pvalue(test_value=my_value, value_df=pd.DataFrame(comp_values))

    # ===== Special case cluster =====
    elif args.mode == "cluster":
        target = pd.read_csv(args.target, header=None, sep="\t")
        # ===== Get validation values of input =====
        my_value_di, my_value_ss = get_validation(ref=reference, ref_id=args.reference_id_type,
                                                  tar=target, tar_id=args.target_id_type, mode=args.mode, mapper=mapper)
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(ref=reference, ref_id=args.reference_id_type,
                                             tar=target, tar_id=args.target_id_type, mode=args.mode,
                                             mapper=mapper)
        p_values_di = eu.calc_pvalue(test_value=my_value_di, value_df=pd.DataFrame(comp_values[0]))
        p_values_ss = eu.calc_pvalue(test_value=my_value_ss, value_df=pd.DataFrame(comp_values[1]))
        p_values = {'di': p_values_di, 'ss': p_values_ss}
    else:
        return None
    mapper.save_mappings()
    ru.print_current_usage('Finished validation')
    return p_values


def get_validation(ref, ref_id, tar, tar_id, mode, mapper: Mapper, verbose=True):
    if mode == "set":
        return comp.compare_set(id_set=tar, id_type=tar_id, mapper=mapper)
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
            disorders_map = pd.read_csv(config.FILES_DIR + 'disorders.map', sep="\t", dtype=str)
            if tar_id == "ICD-10":
                disorders_map = mu.split_and_expand_column(data=disorders_map, split_string=",", column_name="ICD-10")
            full_id_list = set(disorders_map[tar_id])
        else:  # if tar_id in config.SUPPORTED_GENE_IDS:
            full_id_list = set(pd.read_csv(config.FILES_DIR + 'gene.list', sep="\t", dtype=str)['entrez'])
            tar_id = 'entrez'
        # ===== Remove empty values =====
        full_id_list = list(filter(None, full_id_list))
        # ===== Calculate values =====

        for run in range(0, config.NUMBER_OF_RANDOM_RUNS):
            results.append(
                get_validation(ref=ref, ref_id=ref_id, tar=random.sample(full_id_list, len(tar)), tar_id=tar_id,
                               mode=mode, mapper=mapper))
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
    arguments = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o'))
    print(single_validation(args=arguments))
