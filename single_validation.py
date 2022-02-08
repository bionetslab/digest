#!/usr/bin/python3

import numpy as np
import pandas as pd
from d_utils import runner_utils as ru, config, eval_utils as eu, mapping_utils as mu
from evaluation import comparator as comp
import random
from mappers.mapper import Mapper, FileMapper
import json
import time


def single_validation(tar: str, tar_id: str, mode: str, ref: str = None, ref_id: str = None, enriched: bool = False,
                      mapper: Mapper = FileMapper(), out_dir: str = "", runs: int = config.NUMBER_OF_RANDOM_RUNS,
                      background_model: str = "complete", replace=100, verbose: bool = False):
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
    :param background_model
    :param replace
    :param verbose: bool if additional info like ids without assigned attributes should be printed
    """
    ru.start_time = time.time()
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
                                             runs=runs, background_model=background_model, replace=replace)
        ru.print_current_usage('Calculating p-values ...')
        if mode == "set":
            p_values = eu.calc_pvalue(test_value=my_value, random_values=pd.DataFrame(comp_values), maximize=False)
        else:
            p_values = eu.calc_pvalue(test_value=my_value, random_values=pd.DataFrame(comp_values))
        result = {'input_values': my_value, 'p_values': p_values}

    # ===== Special case cluster =====
    elif mode == "cluster":
        ru.print_current_usage('Load distances for input into cache ...')
        mapper.load_distances(set_type=tar_id)
        ru.print_current_usage('Load input data ...')
        comparator = comp.ClusterComparator(mapper=mapper, verbose=verbose)
        comparator.load_target(id_set=pd.read_csv(tar, header=None, sep="\t"), id_type=tar_id)
        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...')
        my_value_di, my_value_ss, my_value_dbi = comparator.compare()
        ru.print_current_usage('Validation of random runs ...')
        comp_values = get_random_runs_values(comparator=comparator, mode=mode, mapper=mapper, tar_id=tar_id,
                                             runs=runs)
        p_values_di = eu.calc_pvalue(test_value=my_value_di, random_values=pd.DataFrame(comp_values[0]))
        p_values_ss = eu.calc_pvalue(test_value=my_value_ss, random_values=pd.DataFrame(comp_values[1]))
        p_values_dbi = eu.calc_pvalue(test_value=my_value_dbi, random_values=pd.DataFrame(comp_values[2]))
        p_values = {'di': p_values_di, 'ss': p_values_ss, 'dbi': p_values_dbi}
        result = {'input_values': {'di': my_value_di, 'ss': my_value_ss, 'dbi': my_value_dbi}, 'p_values': p_values}
    else:
        result = {None}

    # ===== Saving final files and results =====
    ru.print_current_usage('Save files')
    # mapper.save_mappings()
    # mapper.save_distances()
    ru.print_current_usage('Finished validation')
    with open(out_dir + "digest_" + mode + "_result.json", "w") as outfile:
        json.dump(result, outfile)


def get_random_runs_values(comparator: comp.Comparator, mode: str, mapper: Mapper, tar_id: str, runs: int,
                           background_model: str = "complete", replace=100, term: str = "sum") -> list:
    """
    Pick random ids to recreate a target input and run the comparison against reference or itself.
    The random ids are of the same id type of the original target input.

    :param comparator: comparator with defined way how to make the comparison
    :param mode: comparison mode [set, id-set, set-set, cluster]
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :param tar_id: id type of target input
    :param runs: number of random runs to create p-values [Default=1000]
    :param background_model
    :param replace
    :param term
    :return: comparison
    """
    results = list()
    if not mode == "cluster":
        # ===== Get full id mapping =====
        if tar_id in config.SUPPORTED_DISEASE_IDS:
            full_id_map = mapper.get_full_set(id_type=tar_id, mapping_name='disorder_ids')
            new_id_type, map_id_type = "mondo", "disorder_atts"
        else:  # if tar_id in config.SUPPORTED_GENE_IDS:
            full_id_map = mapper.get_full_set(id_type=config.ID_TYPE_KEY[tar_id], mapping_name='gene_ids')
            new_id_type, map_id_type = "entrez", "gene_atts"
        # ===== Precalculate sizes =====
        orig_ids = set(comparator.id_set)
        size = len(orig_ids)
        random_size = int((size / 100) * replace)
        # ===== Precalculate attribute sizes for term-pres =====
        if background_model == "term-pres":
            att_size = atts_to_size(pd_map=mapper.loaded_mappings[map_id_type])
            att_dict = size_mapping_to_dict(pd_size_map=att_size, id_col=config.ID_TYPE_KEY[new_id_type], term_col=term,
                                            threshold=100)
        for run in range(0, runs):
            # ===== Pick new samples =====
            old_sample = set(random.sample(orig_ids, (size - random_size)))
            if background_model == "complete":
                full_id_list = list(
                    set(full_id_map[full_id_map[config.ID_TYPE_KEY[tar_id]] != ""][config.ID_TYPE_KEY[tar_id]]))
                random_sample = set(random.sample(full_id_list, random_size))
            elif background_model == "term-pres":
                to_replace = orig_ids.difference(old_sample)
                to_replace = full_id_map[full_id_map[config.ID_TYPE_KEY[tar_id]].isin(to_replace)][comparator.att_id]
                random_sample = set()
                for replace_id in to_replace:
                    random_sample.add(
                        att_size[att_size[term].isin(att_dict[replace_id])][config.ID_TYPE_KEY[new_id_type]].sample(
                            n=1).values[0])
            # ===== Get corresponding id set =====
            id_set = full_id_map[full_id_map[config.ID_TYPE_KEY[tar_id]].isin(random_sample.union(old_sample))][
                comparator.att_id]
            # ===== Calculate values =====
            comparator.load_target(id_set=set(id_set), id_type=new_id_type)
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


def atts_to_size(pd_map: pd.DataFrame) -> pd.DataFrame:
    att_len = pd_map.copy()
    att_len[att_len.columns[1:]] = att_len[att_len.columns[1:]].applymap(mu.set_to_len)
    att_len['sum'] = att_len[att_len.columns[1:]].sum(axis=1)
    return att_len


def size_mapping_to_dict(pd_size_map: pd.DataFrame, id_col: str, term_col: str, threshold: int = 100):
    size_to_occ = pd.DataFrame(pd_size_map[term_col].value_counts()).sort_index().to_dict()
    pd_size_map = pd_size_map.sort_values(by=[id_col]).reset_index(drop=True)
    new_dict = dict()
    term_sizes = pd_size_map[term_col].unique().tolist()
    for index, key in enumerate(term_sizes):
        curr_keys = [key]
        if size_to_occ[key] < threshold:
            sum_tmp, add_top, add_bottom = size_to_occ[key], index, index
            while sum_tmp < threshold:
                if add_top - 1 >= 0:
                    add_top = add_top - 1
                    sum_tmp = sum_tmp + size_to_occ[term_sizes[add_top]]
                    curr_keys.append(term_sizes[add_top])
                if add_bottom + 1 < len(term_sizes):
                    add_bottom = add_bottom + 1
                    sum_tmp = sum_tmp + size_to_occ[term_sizes[add_bottom]]
                    curr_keys.append(term_sizes[add_bottom])
        for cur_id in pd_size_map[pd_size_map[term_col] == key][id_col]:
            new_dict[cur_id] = curr_keys
    return new_dict


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    args = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o', 'e', 'c', 'v', 'b', 'pr'))
    single_validation(tar=args.target, tar_id=args.target_id_type,
                      mode=args.mode, ref=args.reference, ref_id=args.reference_id_type,
                      enriched=args.enriched, out_dir=args.out_dir, runs=args.runs,
                      background_model=args.background_model, replace=args.replace)
