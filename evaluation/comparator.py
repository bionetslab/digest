#!/usr/bin/python3

import pandas as pd
import numpy as np
from d_utils import runner_utils as ru, eval_utils as eu, config as c
from mappers import gene_mapper as gm, disease_mapper as dm
from evaluation import score_calculator as sc


def compare_set(id_set, id_type):
    """
    Compare the set on itself based on connected attributes. See config for more info.

    :param id_set: set of genes or diseases
    :param id_type: id type of the set
    :return: evaluation value for each connected attribute
    """
    if id_type in c.SUPPORTED_DISEASE_IDS:
        mapping = dm.get_disease_to_attributes(disease_set=id_set, id_type=id_type)
    else:  # if id_type in c.SUPPORTED_GENE_IDS:
        mapping = gm.get_gene_to_attributes(gene_set=id_set, id_type=id_type)
    result = dict()
    for attribute in mapping.columns[1:]:
        subset_df = mapping[mapping[attribute].str.len() > 0].reset_index()
        missing_values = len(mapping) - len(subset_df)
        print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(len(id_set)))
        comp_mat = eu.get_distance_matrix(eval_df=subset_df[attribute])
        comp_mean = (comp_mat.sum() - np.diag(comp_mat).sum()) / (
                len(subset_df[attribute]) * (len(subset_df[attribute]) - 1))
        result[attribute] = comp_mean
    return result


def compare_set_to_set(ref, ref_id_type, tar, tar_id_type, threshold=0.0):
    """
    Compare two sets of the same type (either genes or diseases) with each other.
    The tar set is evaluated how good it matches to the ref set.

    :param ref: reference set of genes or diseases
    :param ref_id_type: id type of reference set
    :param tar: target set of genes or diseases
    :param tar_id_type: id type of target set
    :param threshold: minimum similarity of each element in tar to reference set [0,1]
    :return: evaluation value for each connected attribute
    """
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
    """
    Compare disease id to either gene or disease set.
    The tar set is evaluated how good it matches to the ref id.

    :param ref_id: reference disease id
    :param ref_id_type: id type of reference set
    :param tar: target set of genes or diseases
    :param tar_id_type: id type of target set
    :param threshold: minimum similarity of each element in tar to reference set [0,1]
    :return: evaluation value for each connected attribute
    """
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
    """
    Evaluate the quality of clustering of given set with diseases or genes and
    assigned clusters. Additionally calculate statistical values with
    silhouette score and dunn index.

    :param clusters: set of genes or diseases and assigned clusters
    :param id_type: id type of the set
    :return: evaluation value for each connected attribute; silhouette score and dunn index
    """
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
        dist_mat = eu.get_distance_matrix(eval_df=subset_df[attribute])
        dist_df = pd.DataFrame(dist_mat, columns=subset_df[id_type], index=subset_df[id_type])
        ss_score = sc.silhouette_score(distance_matrix=dist_df, ids_cluster=subset_clusters)
        di_score = sc.dunn_index(distance_matrix=dist_df, ids_cluster=subset_clusters)
        result.append([attribute, di_score, ss_score[0], ss_score[1]])
    return result