#!/usr/bin/python3

import pandas as pd
import numpy as np
from d_utils import runner_utils as ru, eval_utils as eu
from mappers import gene_mapper as gm, disease_mapper as dm
from stats import score_calculator as sc


def single_validation():
    return


def compare_set_to_set(ref, ref_id_type, targets, targets_id_type, threshold=0.0):
    reference_mapping = gm.get_gene_to_attributes(ref, ref_id_type)
    target_mapping = gm.get_gene_to_attributes(targets, targets_id_type)
    ref_dict = eu.create_ref_dict(mapping=reference_mapping, keys=reference_mapping.columns[1:])
    result = eu.evaluate_values(mapping=target_mapping, ref_dict=ref_dict, threshold=threshold,
                                keys=target_mapping.columns[1:])
    return result


def compare_id_to_set(ref_id, ref_id_type, targets, targets_id_type, threshold=0.0):
    disease_id_atts = dm.get_disease_mapping({ref_id}, ref_id_type)
    disease_id_atts = disease_id_atts.rename(columns={'ctd.pathway_related_to_disease': 'pathway.kegg'})
    target_mapping = gm.get_gene_to_attributes(targets, targets_id_type)
    ref_dict = eu.create_ref_dict(mapping=disease_id_atts, keys=['pathway.kegg'])
    result = eu.evaluate_values(mapping=target_mapping, ref_dict=ref_dict, threshold=threshold, keys=['pathway.kegg'])
    return result


def compare_ids(disease_id_set, id_type):
    disease_df = dm.get_disease_mapping(disease_set=disease_id_set, id_type=id_type)
    result = list()
    for attribute in disease_df.columns[1:]:
        subset_df = disease_df[disease_df[attribute].str.len() > 0].reset_index()
        missing_values = len(disease_id_set) - len(subset_df)
        print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(len(disease_id_set)))
        comp_mat = eu.get_distance_matrix(eval_set=subset_df[attribute])
        comp_mean = (comp_mat.sum() - np.diag(comp_mat).sum()) / (
                    len(subset_df[attribute]) * (len(subset_df[attribute]) - 1))
        result.append([attribute, comp_mean])
    return result


def compare_id_clusters(disease_clusters, id_type):
    disease_clusters['cluster_index'] = disease_clusters.groupby(1).ngroup()
    disease_clusters_df = dm.get_disease_mapping(disease_set=disease_clusters[0], id_type=id_type)
    result = list()
    for attribute in disease_clusters_df.columns[1:]:
        subset_df = disease_clusters_df[disease_clusters_df[attribute].str.len() > 0].reset_index(drop=True)
        subset_clusters = disease_clusters[disease_clusters[0].isin(subset_df[id_type])]
        missing_values = len(disease_clusters) - len(subset_df)
        print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(len(disease_clusters_df)))
        dist_mat = eu.get_distance_matrix(subset_df[attribute])
        dist_df = pd.DataFrame(dist_mat, columns=subset_df[id_type], index=subset_df[id_type])
        ss_score = sc.silhouette_score(distance_matrix=dist_df, ids_cluster=subset_clusters)
        di_score = sc.dunn_index(distance_matrix=dist_df, ids_cluster=subset_clusters)
        result.append([attribute, di_score, ss_score[0], ss_score[1]])
    return result


if __name__ == "__main__":
    desc = "Evaluation of gene sets against seed gene sets or disease id; disease id sets and disease clusters."
    args = ru.save_parameters(script_desc=desc, arguments=('m', 'o', 'r'))
    #single_validation(mapping_file=args.disease_mapping, out_dir=args.out_dir)
