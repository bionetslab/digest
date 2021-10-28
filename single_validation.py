#!/usr/bin/python3

import pandas as pd
import numpy as np
from d_utils import runner_utils as ru, eval_utils as eu, config as c
from mappers import gene_mapper as gm, disease_mapper as dm
from stats import score_calculator as sc


def single_validation(args):
    if args.mode == "set":
        print(compare_set(id_set=pd.read_csv(args.target, header=None, sep="\t")[0], id_type=args.target_id_type))
    if args.mode == "set-set":
        print(compare_set_to_set(ref=pd.read_csv(args.reference, header=None, sep="\t")[0],
                                 ref_id_type=args.reference_id_type,
                                 tar=pd.read_csv(args.target, header=None, sep="\t")[0],
                                 tar_id_type=args.target_id_type))
    if args.mode == "id-set":
        print(compare_id_to_set(ref_id=args.reference, ref_id_type=args.reference_id_type,
                                tar=pd.read_csv(args.target, header=None, sep="\t")[0],
                                tar_id_type=args.target_id_type))
    if args.mode == "cluster":
        print(
            compare_clusters(clusters=pd.read_csv(args.target, header=None, sep="\t"), id_type=args.target_id_type))
    return


def compare_set(id_set, id_type):
    if id_type in c.SUPPORTED_DISEASE_IDS:
        mapping = dm.get_disease_to_attributes(disease_set=id_set, id_type=id_type)
    else:  # if id_type in c.SUPPORTED_GENE_IDS:
        mapping = gm.get_gene_to_attributes(gene_set=id_set, id_type=id_type)
    result = list()
    for attribute in mapping.columns[1:]:
        subset_df = mapping[mapping[attribute].str.len() > 0].reset_index()
        missing_values = len(mapping) - len(subset_df)
        print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(len(id_set)))
        comp_mat = eu.get_distance_matrix(eval_set=subset_df[attribute])
        comp_mean = (comp_mat.sum() - np.diag(comp_mat).sum()) / (
                len(subset_df[attribute]) * (len(subset_df[attribute]) - 1))
        result.append([attribute, comp_mean])
    return result


def compare_set_to_set(ref, ref_id_type, tar, tar_id_type, threshold=0.0):
    if ref_id_type in c.SUPPORTED_DISEASE_IDS and tar_id_type in c.SUPPORTED_DISEASE_IDS:
        reference_mapping = dm.get_disease_to_attributes(disease_set=ref, id_type=ref_id_type)
        target_mapping = dm.get_disease_to_attributes(disease_set=tar, id_type=tar_id_type)
    else:  # if ref_id_type in c.SUPPORTED_GENE_IDS and targets_id_type in c.SUPPORTED_GENE_IDS:
        reference_mapping = gm.get_gene_to_attributes(gene_set=ref, id_type=ref_id_type)
        target_mapping = gm.get_gene_to_attributes(gene_set=tar, id_type=tar_id_type)
    ref_dict = eu.create_ref_dict(mapping=reference_mapping, keys=reference_mapping.columns[1:])
    result = eu.evaluate_values(mapping=target_mapping, ref_dict=ref_dict, threshold=threshold,
                                keys=target_mapping.columns[1:])
    return result


def compare_id_to_set(ref_id, ref_id_type, tar, tar_id_type, threshold=0.0):
    id_mapping = dm.get_disease_to_attributes({ref_id}, ref_id_type)
    if tar_id_type in c.SUPPORTED_DISEASE_IDS:
        target_mapping = dm.get_disease_to_attributes(disease_set=tar, id_type=tar_id_type)
        keys = id_mapping.columns[1:]
    else:  # if targets_id_type in c.SUPPORTED_GENE_IDS:
        id_mapping = id_mapping.rename(columns={'ctd.pathway_related_to_disease': 'pathway.kegg'})
        target_mapping = gm.get_gene_to_attributes(gene_set=tar, id_type=tar_id_type)
        keys = ['pathway.kegg']
    ref_dict = eu.create_ref_dict(mapping=id_mapping, keys=keys)
    result = eu.evaluate_values(mapping=target_mapping, ref_dict=ref_dict, threshold=threshold, keys=keys)
    return result


def compare_clusters(clusters, id_type):
    clusters['cluster_index'] = clusters.groupby(1).ngroup()
    if id_type in c.SUPPORTED_DISEASE_IDS:
        clusters_mapping = dm.get_disease_to_attributes(disease_set=clusters[0], id_type=id_type)
    else:  # if id_type in c.SUPPORTED_GENE_IDS:
        clusters_mapping = gm.get_gene_to_attributes(gene_set=clusters[0], id_type=id_type)
    result = list()
    for attribute in clusters_mapping.columns[1:]:
        subset_df = clusters_mapping[clusters_mapping[attribute].str.len() > 0].reset_index(drop=True)
        subset_clusters = clusters[clusters[0].isin(subset_df[id_type])]
        missing_values = len(clusters) - len(subset_df)
        print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(len(clusters_mapping)))
        dist_mat = eu.get_distance_matrix(subset_df[attribute])
        dist_df = pd.DataFrame(dist_mat, columns=subset_df[id_type], index=subset_df[id_type])
        ss_score = sc.silhouette_score(distance_matrix=dist_df, ids_cluster=subset_clusters)
        di_score = sc.dunn_index(distance_matrix=dist_df, ids_cluster=subset_clusters)
        result.append([attribute, di_score, ss_score[0], ss_score[1]])
    return result


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    arguments = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o'))
    single_validation(args=arguments)
