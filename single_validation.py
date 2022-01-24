#!/usr/bin/python3

import numpy as np
import pandas as pd
from d_utils import runner_utils as ru, config, eval_utils as eu
from evaluation import comparator as comp
import random
from mappers.mapper import Mapper, FileMapper
import json


def single_validation(tar: str, tar_id: str, mode: str, ref: str = None, ref_id: str = None, enriched: bool = False,
                      mapper: Mapper = FileMapper(), out_dir: str = "", runs: int = config.NUMBER_OF_RANDOM_RUNS,
                      verbose: bool = False):
    """
    Single validation of a set, cluster a id versus set and set versus set.

    :param tar: path to the file with the target input
    :param tar_id: id type of target input
    :param mode: comparison mode [set, id-set, set-set, cluster]
    :param ref: path to the file with the reference input or string with reference id [Default=None]
    :param ref_id: id type of reference input
    :param enriched: bool setting if values of reference set should be filtered for enriched values [Default=False]
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :param out_dir: output directory for results
    :param runs: number of random runs to create p-values [Default=1000]
    :param verbose: bool if additional info like ids without assigned attributes should be printed
    """
    # ===== Comparison with a set =====
    ru.print_current_usage('Starting validation ...')
    ru.print_current_usage('Load mappings for input into cache ...')
    mapper.load_mappings()

    if mode in ["set", "set-set", "id-set"]:
        if mode == "set-set":
            comparator = comp.SetSetComparator(mapper=mapper, enriched=enriched, verbose=verbose)
            comparator.load_reference(ref=pd.read_csv(ref, header=None, sep="\t")[0], ref_id_type=ref_id)
        elif mode == "id-set":
            comparator = comp.IDSetComparator(mapper=mapper, verbose=verbose)
            comparator.load_reference(ref=ref, ref_id_type=ref_id)
        else:  # mode == "set"
            comparator = comp.SetComparator(mapper=mapper, verbose=verbose)
            ru.print_current_usage('Load distances for input into cache ...')
            mapper.load_distances(set_type=tar_id)
        comparator.load_target(id_set=pd.read_csv(tar, header=None, sep="\t")[0], id_type=tar_id)

        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...')
        my_value = comparator.compare()
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(comparator=comparator, mode=mode, mapper=mapper, tar_id=tar_id,
                                             runs=runs)
        ru.print_current_usage('Calculating p-values ...')
        if mode == "set":
            p_values = eu.calc_pvalue(test_value=my_value, random_values=pd.DataFrame(comp_values), maximize=False)
        else:
            p_values = eu.calc_pvalue(test_value=my_value, random_values=pd.DataFrame(comp_values))
        result = {'input_values': my_value, 'p_values:': p_values}

    # ===== Special case cluster =====
    elif args.mode == "cluster":
        ru.print_current_usage('Load distances for input into cache ...')
        mapper.load_distances(set_type=tar_id)
        ru.print_current_usage('Load input data ...')
        comparator = comp.ClusterComparator(mapper=mapper, verbose=verbose)
        comparator.load_target(id_set=pd.read_csv(tar, header=None, sep="\t"), id_type=tar_id)
        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...')
        my_value_di, my_value_ss = comparator.compare()
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(comparator=comparator, mode=mode, mapper=mapper, tar_id=tar_id,
                                             runs=runs)
        p_values_di = eu.calc_pvalue(test_value=my_value_di, random_values=pd.DataFrame(comp_values[0]))
        p_values_ss = eu.calc_pvalue(test_value=my_value_ss, random_values=pd.DataFrame(comp_values[1]))
        p_values = {'di': p_values_di, 'ss': p_values_ss}
        result = {'input_values': {'di': my_value_di, 'ss': my_value_ss}, 'p_values:': p_values}
    else:
        result = {None}

    # ===== Saving final files and results =====
    ru.print_current_usage('Save files')
    mapper.save_mappings()
    mapper.save_distances()
    ru.print_current_usage('Finished validation')
    with open(out_dir+"digest_"+args.mode+"_result.json", "w") as outfile:
        json.dump(result, outfile)


def get_random_runs_values(comparator: comp.Comparator, mode: str, mapper: Mapper, tar_id: str, runs: int) -> list:
    """
    Pick random ids to recreate a target input and run the comparison against reference or itself.
    The random ids are of the same id type of the original target input.

    :param comparator: comparator with defined way how to make the comparison
    :param mode: comparison mode [set, id-set, set-set, cluster]
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :param tar_id: id type of target input
    :param runs: number of random runs to create p-values [Default=1000]
    :return: comparison
    """
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
        for run in range(0, runs):
            comparator.load_target(id_set=random.sample(full_id_list, len(comparator.id_set)), id_type=tar_id)
            results.append(comparator.compare())

    # ===== Special case cluster =====
    else:
        # ===== Calculate values =====
        results.extend([list(), list()])
        for run in range(0, runs):
            comparator.clustering['cluster_index'] = np.random.permutation(comparator.clustering['cluster_index'])
            value_di, value_ss = comparator.compare()
            results[0].append(value_di)
            results[1].append(value_ss)

    return results


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    args = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o', 'e', 'c', 'v'))
    single_validation(tar=args.target, tar_id=args.target_id_type,
                      mode=args.mode, ref=args.reference, ref_id=args.reference_id_type,
                      enriched=args.enriched, out_dir=args.out_dir, runs=args.runs)
