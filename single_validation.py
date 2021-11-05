#!/usr/bin/python3

import pandas as pd
from d_utils import runner_utils as ru, config, mapping_utils as mu
from evaluation import comparator as comp
import random


def single_validation(args):
    """
    Single validation of a set, cluster a id versus set and set versus set.

    :param args: user given parameters
    :return:
    """
    reference, target = None, None
    # ===== Comparison with a set =====
    if args.mode in ["set", "set-set", "id-set"]:
        target = pd.read_csv(args.target, header=None, sep="\t")[0]
        if args.mode == "set-set":
            reference = pd.read_csv(args.reference, header=None, sep="\t")[0]
        elif args.mode == "id-set":
            reference = args.reference
        # ===== Get validation values of input =====
        my_value = get_validation(ref=reference, ref_id=args.reference_id_type,
                                  tar=target, tar_id=args.target_id_type, mode=args.mode)
        comp_values = get_random_runs_values(ref=reference, ref_id=args.reference_id_type,
                                             tar=target, tar_id=args.target_id_type, mode=args.mode)

    # ===== Special case cluster =====
    elif args.mode == "cluster":
        target = pd.read_csv(args.target, header=None, sep="\t")
        # ===== Get validation values of input =====
        my_value = get_validation(ref=reference, ref_id=args.reference_id_type,
                                  tar=target, tar_id=args.target_id_type, mode=args.mode)
    else:
        return None
    return comp_values


def get_validation(ref, ref_id, tar, tar_id, mode):
    if mode == "set":
        return comp.compare_set(id_set=tar, id_type=tar_id)
    elif mode == "set-set":
        return comp.compare_set_to_set(ref=ref, ref_id_type=ref_id, tar=tar, tar_id_type=tar_id)
    elif mode == "id-set":
        return comp.compare_id_to_set(ref_id=ref, ref_id_type=ref_id, tar=tar, tar_id_type=tar_id)
    elif mode == "cluster":
        return comp.compare_clusters(clusters=tar, id_type=tar_id)
    else:
        return None


def get_random_runs_values(ref, ref_id, tar, tar_id, mode):
    # pd.read_csv(file=config.FILES_DIR + 'random_runs.csv', sep=",")
    # ===== On the fly =====
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
    results = list()
    for run in range(0, config.NUMBER_OF_RANDOM_RUNS):
        results.append(get_validation(ref=ref, ref_id=ref_id,
                                      tar=random.sample(full_id_list, len(tar)), tar_id=tar_id, mode=mode))
    return results


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    arguments = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o'))
    print(single_validation(args=arguments))
