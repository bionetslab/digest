#!/usr/bin/python3

import os
import numpy as np
import pandas as pd
from evaluation.d_utils import runner_utils as ru, eval_utils as eu, plotting_utils as pu
from evaluation.mappers.mapper import Mapper, FileMapper
from evaluation import config, comparator as comp, background_models as bm
import random
import json
import time
from pathlib import Path
from typing import Union, Callable


def single_validation(tar: Union[pd.DataFrame, set], tar_id: str, mode: str, distance: str = "jaccard",
                      ref: set = None, ref_id: str = None, enriched: bool = False,
                      mapper: Mapper = FileMapper(), runs: int = config.NUMBER_OF_RANDOM_RUNS,
                      background_model: str = "complete", replace=100, verbose: bool = False,
                      network_data: dict = None,
                      progress: Callable[[float, str], None] = None):
    """
    Single validation of a set, cluster a id versus set and set versus set.

    :param tar: cluster as type dataframe with columns=["id","cluster"] or set as type set
    :param tar_id: id type of target input
    :param mode: comparison mode [set, set-set, clustering, subnetwork, subnetwork-set]
    :param distance: distance measure used for comparison. [Default=jaccard]
    :param ref: set of reference ids [Default=None]
    :param ref_id: id type of reference input [Default=None]
    :param enriched: bool setting if values of reference set should be filtered for enriched values [Default=False]
    :param mapper: mapper from type Mapper defining where the precalculated information comes from [Default=FileMapper]
    :param runs: number of random runs to create p-values [Default=1000]
    :param background_model: which background model to use for random picks [Default="complete"]
    :param replace: how many % of target input should be replaced by random picks [Default=100]
    :param verbose: bool if additional info like ids without assigned attributes should be printed [Default=False]
    :param progress: method that will get a float [0,1] and message indicating the current progress
    """
    ru.start_time = time.time()
    ru.print_current_usage('Check for proper setup ...')
    mapper.check_for_setup_sources()
    # ===== Comparison with a set =====
    ru.print_current_usage('Starting validation ...')
    progress(0.01, "Load mappings...") if progress is not None else None
    if mapper.load:
        ru.print_current_usage('Load mappings for input into cache ...') if verbose else None
        mapper.load_mappings()
    if mode in ["set", "set-set", "subnetwork", "subnetwork-set"]:
        error_mappings = []
        if mode in ["set-set", "subnetwork-set"]:
            comparator = comp.SetSetComparator(mapper=mapper, enriched=enriched, verbose=verbose,
                                               distance_measure=distance)
            comparator.load_reference(ref=ref, ref_id_type=ref_id, tar_id_type=tar_id)
            if not comparator.ref_dict:
                error_mappings.append("reference")
        else:  # mode == "set"
            comparator = comp.SetComparator(mapper=mapper, verbose=verbose, distance_measure=distance)
            if mapper.load:
                ru.print_current_usage('Load distances for input into cache ...')
                progress(0.05, "Load distances...") if progress is not None else None
                mapper.load_distances(set_type=tar_id, distance_measure=distance)
        comparator.load_target(id_set=tar, id_type=tar_id)
        # ===== Check if mappings possible =====
        if comparator.mapping.empty:
            error_mappings.append("target set")
        if len(error_mappings) > 0:
            ru.print_current_usage(
                'Exit validation: No mapping found for ' + ' and '.join(error_mappings)) if verbose else None
            return {'status': 'No mapping found for ' + ' and '.join(error_mappings),
                    'input_values': {'values': None, 'mapped_ids': []}, 'random_values': None,
                    'p_values': {'values': None}}
        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...') if verbose else None
        progress(0.1, "Validation of input...") if progress is not None else None
        my_value, mapped = comparator.compare()
        ru.print_current_usage('Validation of random runs ...') if verbose else None
        progress(0.1 + (0.9 / min(runs + 1, 100)),
                 "Validation with background model...") if progress is not None else None
        comparator.verbose = False
        comparator.input_run = False
        comp_values = get_random_runs_values(comparator=comparator, mode=mode, mapper=mapper, tar_id=tar_id,
                                             runs=runs, background_model=background_model, replace=replace,
                                             network_data=network_data, progress=progress)
        # ===== Statistical analysis =====
        ru.print_current_usage('Calculating p-values ...') if verbose else None
        set_value = eu.calc_pvalue(test_value=my_value, random_values=comp_values[0], maximize=True)
        measure_short = {"jaccard": "JI-based", "overlap": "OC-based"}
        p_values = {measure_short[distance]: set_value}
        results = {'status': 'ok',
                   'input_values': {'values': {measure_short[distance]: my_value}, 'mapped_ids': mapped},
                   'random_values': {measure_short[distance]: comp_values[0].to_dict('list')},
                   'p_values': {'values': p_values}}

    # ===== Special case cluster =====
    elif mode == "clustering":
        if mapper.load:
            ru.print_current_usage('Load distances for input into cache ...') if verbose else None
            progress(0.05, "Load distances...") if progress is not None else None
            mapper.load_distances(set_type=tar_id, distance_measure=distance)
        ru.print_current_usage('Load input data ...') if verbose else None
        comparator = comp.ClusterComparator(mapper=mapper, verbose=verbose, distance_measure=distance)
        comparator.load_target(id_set=tar, id_type=tar_id)
        # ===== Check if mappings possible =====
        if comparator.mapping.empty:
            ru.print_current_usage(
                'Exit validation: No mapping found for target cluster') if verbose else None
            return {'status': 'No mapping found for target cluster',
                    'input_values': {'values': None, 'mapped_ids': []},
                    'random_values': None, 'p_values': {'values': None}}
        # ===== Get validation values of input =====
        ru.print_current_usage('Validation of input ...') if verbose else None
        progress(0.1, "Validation of input...") if progress is not None else None
        my_value_di, my_value_ss, my_value_dbi, my_value_ss_inter, mapped = comparator.compare()
        ru.print_current_usage('Validation of random runs ...') if verbose else None
        progress(0.1 + (0.9 / min(runs + 1, 100)),
                 "Validation with background model...") if progress is not None else None
        comparator.verbose = False
        comparator.input_run = False
        comp_values = get_random_runs_values(comparator=comparator, mode=mode, mapper=mapper, tar_id=tar_id,
                                             runs=runs, progress=progress)
        # ===== Statistical analysis =====
        ru.print_current_usage('Calculating p-values ...') if verbose else None
        p_values_di = eu.calc_pvalue(test_value=my_value_di, random_values=comp_values[0], maximize=False)
        p_values_ss = eu.calc_pvalue(test_value=my_value_ss, random_values=comp_values[1], maximize=True)
        p_values_dbi = eu.calc_pvalue(test_value=my_value_dbi, random_values=comp_values[2],
                                      maximize=False)
        p_values = {'DI-based': p_values_di, 'SS-based': p_values_ss, 'DBI-based': p_values_dbi}
        results = {'status': 'ok',
                   'input_values': {'values': {'DI-based': my_value_di, 'SS-based': my_value_ss,
                                               'DBI-based': my_value_dbi},
                                    'values_inter': my_value_ss_inter,
                                    'mapped_ids': mapped},
                   'random_values': {'DI-based': comp_values[0].to_dict('list'),
                                     'SS-based': comp_values[1].to_dict('list'),
                                     'DBI-based': comp_values[2].to_dict('list')},
                   'p_values': {'values': p_values}}
    else:
        results = {None}
    ru.print_current_usage('Finished validation')
    return results


def get_random_runs_values(comparator: comp.Comparator, mode: str, mapper: Mapper, tar_id: str, runs: int,
                           background_model: str = "complete", replace=100, term: str = "sum",
                           network_data: dict = None, progress: Callable[[float, str], None] = None) -> list:
    """
    Pick random ids to recreate a target input and run the comparison against reference or itself.
    The random ids are of the same id type of the original target input.

    :param comparator: comparator with defined way how to make the comparison
    :param mode: comparison mode [set, id-set, set-set, cluster]
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :param tar_id: id type of target input
    :param runs: number of random runs to create p-values [Default=1000]
    :param background_model: which background model to use for random picks [Default="complete"]
    :param replace: how many % of target input should be replaced by random picks [Default=100]
    :param term: on what term the term preserving background model should calculate [Default="sum"]
    :param network_data: ={"network_file": network_file, "prop_name": "name", "id_type": "symbol"},
    :param progress: method that will get a float [0,1] and message indicating the current progress
    :return: comparison
    """
    results = list()
    threshold = min(runs + 1, 100)
    limit, counter = int((runs + 1) / threshold), 1
    if not mode == "clustering":
        results.extend([list()])
        # ===== Get full id mapping =====
        if tar_id in config.SUPPORTED_DISEASE_IDS:
            full_id_map = mapper.get_full_set(id_type=tar_id, mapping_name='disorder_ids')
            new_id_type, map_id_type, map_att_type = "mondo", "disorder_ids", "disorder_atts"
        else:  # if tar_id in config.SUPPORTED_GENE_IDS:
            full_id_map = mapper.get_full_set(id_type=config.ID_TYPE_KEY[tar_id], mapping_name='gene_ids')
            new_id_type, map_id_type, map_att_type = "entrez", "gene_ids", "gene_atts"
        # ===== Precalculate sizes =====
        orig_ids = set(comparator.id_set)
        size = len(orig_ids)
        random_size = int((size / 100) * replace)

        if background_model == "complete":
            background = bm.CompleteModel(prev_id_type=tar_id, full_id_map=full_id_map)
        elif background_model == "term-pres":
            background = bm.TermPresModel(mapper=mapper, prev_id_type=tar_id, new_id_type=new_id_type,
                                          map_id_type=map_id_type, map_att_type=map_att_type, term=term)
        elif background_model == "network":
            if network_data is None:  # load default network from config
                print("load default network from config")
            if network_data["id_type"] != tar_id:  # remap ids
                input_ids = set(full_id_map[full_id_map[config.ID_TYPE_KEY[tar_id]].isin(orig_ids)][
                    config.ID_TYPE_KEY[network_data["id_type"]]])
                background = bm.NetworkModel(network_data=network_data, to_replace=input_ids, N=runs)
            else:
                background = bm.NetworkModel(network_data=network_data, to_replace=orig_ids, N=runs)

        else:
            return list()
        for run in range(1, runs + 1):
            # ===== Update progress =====
            if run % limit == 0:
                progress(0.1 + (counter * (0.9 / threshold)),
                         str(
                             limit * counter) + " run(s) with background model finished...") if progress is not None else None
                counter += 1
            # ===== Pick new samples =====
            old_sample = set(random.sample(list(orig_ids), (size - random_size)))  # ignore when network model
            to_replace = orig_ids.difference(old_sample)

            if background_model == "complete":
                random_sample = background.get_module(to_replace=to_replace)
            elif background_model == "term-pres":
                random_sample = background.get_module(to_replace=to_replace, term=term, prev_id_type=tar_id)
            elif background_model == "network":
                random_sample = background.get_module(index=run - 1)
                if network_data["id_type"] != tar_id:  # remap ids
                    random_sample = \
                    set(full_id_map[full_id_map[config.ID_TYPE_KEY[network_data["id_type"]]].isin(random_sample)][
                        config.ID_TYPE_KEY[tar_id]])
            else:
                return list()

            # ===== Get corresponding id set =====
            id_set = full_id_map[full_id_map[config.ID_TYPE_KEY[tar_id]].isin(random_sample.union(old_sample))][
                comparator.att_id]
            # ===== Calculate values =====
            comparator.load_target(id_set=set(id_set), id_type=new_id_type)
            result, _ = comparator.compare()
            results[0].append(result)

    # ===== Special case cluster =====
    else:
        # ===== Calculate values =====
        results.extend([list(), list(), list()])
        # ===== Precalculate sizes =====
        orig_ids = set(comparator.clustering['id'])
        orig_clusters = comparator.clustering
        size = len(orig_ids)
        random_size = int((size / 100) * replace)
        for run in range(0, runs):
            # ===== Update progress =====
            if run % limit == 0:
                progress(0.1 + (counter * (0.9 / threshold)),
                         str(
                             limit * counter) + " run(s) with background model finished...") if progress is not None else None
                counter += 1
            old_sample = set(random.sample(orig_ids, (size - random_size)))
            # ===== Shuffle subset of clusters =====
            subset = orig_clusters[~orig_clusters['id'].isin(old_sample)]
            subset["cluster_index"] = np.random.permutation(subset["cluster_index"])
            comparator.clustering = pd.concat([orig_clusters[orig_clusters['id'].isin(old_sample)], subset])
            # ===== Start validating =====
            value_di, value_ss, value_dbi, value_ss_inter, mapped = comparator.compare()
            results[0].append(value_di)
            results[1].append(value_ss)
            results[2].append(value_dbi)
    # convert to dataframes
    final_results = list()
    for result in results:
        final_results.append(pd.DataFrame(result))
    return final_results


def save_results(results: dict, prefix: str, out_dir):
    # ===== Save complete output =====
    with open(os.path.join(out_dir, prefix + "_result.json"), "w") as outfile:
        json.dump(results, outfile)
    # ===== Save output to tables =====
    if results["status"] == "ok":
        pd.DataFrame(results["input_values"]['values']).to_csv(os.path.join(out_dir, prefix + "_input_validation.csv"))
        pd.DataFrame(results["p_values"]['values']).to_csv(os.path.join(out_dir, prefix + "_p-value_validation.csv"))


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    args = ru.save_parameters(script_desc=desc,
                              arguments=('r', 'ri', 't', 'ti', 'm', 'o', 'e', 'c', 'v', 'b', 'pr', 'p', 'dg'))
    res = single_validation(tar=args.target, tar_id=args.target_id_type, verbose=args.verbose,
                            mode=args.mode, ref=args.reference, ref_id=args.reference_id_type,
                            enriched=args.enriched, runs=args.runs, distance=args.distance_measure,
                            background_model=args.background_model, replace=args.replace)
    # ===== Saving final files and results =====
    ru.print_current_usage('Save files') if args.verbose else None
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)  # make sure output dir exists
    pref = "_".join([args.mode, args.target_id_type, args.background_model])
    save_results(results=res, prefix=pref, out_dir=args.out_dir)
    if args.plot and res["status"] == "ok":
        pu.create_plots(results=res, mode=args.mode, tar=args.target, tar_id=args.target_id_type,
                        out_dir=args.out_dir, prefix=pref)
