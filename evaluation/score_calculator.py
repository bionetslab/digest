#!/usr/bin/python3
import pandas as pd
from mappers.mapper import Mapper


def min_inter_distance_to_entity(entity: pd.Series, distances: dict, current_cluster, ids_cluster: pd.DataFrame,
                                 ids_mapping: pd.DataFrame, ids: dict):
    min_distance = None
    for cluster in ids_cluster['cluster_index'].unique():
        if cluster != current_cluster:
            cluster_member = ids_cluster[ids_cluster['cluster_index'] == cluster][0]
            cluster_member = ids_mapping[ids_mapping[ids['id_type']].isin(set(cluster_member))][ids['att_id']]
            distance = inter_distance(cluster_one=entity, cluster_two=cluster_member, distances=distances,
                                      linkage='average')
            if min_distance is None or min_distance > distance:
                min_distance = distance
    return min_distance


def intra_distance(cluster: pd.Series, mapper: Mapper, linkage, ids: dict, entity: pd.Series = None):
    if len(cluster) == 1:
        return 0
    if entity is not None:
        sub_mat_1 = mapper.get_loaded_distances(in_series=cluster, id_type=ids['sparse_key'],
                                                key=ids['attribute'], to_series=entity)
        sub_mat_2 = mapper.get_loaded_distances(in_series=entity, id_type=ids['sparse_key'],
                                                key=ids['attribute'], to_series=cluster)
        if linkage == 'average':
            return (sub_mat_1.sum() + sub_mat_2.sum()) / (len(cluster)*len(entity))
        if linkage == 'complete':
            return max(sub_mat_1.max(), sub_mat_2.max())
    else:
        sub_mat = mapper.get_loaded_distances(in_series=cluster, id_type=ids['sparse_key'],
                                              key=ids['attribute'])
        if linkage == 'average':
            return sub_mat.sum() / (len(cluster)*(len(cluster)-1))/2
        if linkage == 'complete':
            return sub_mat.max()


def inter_distance(cluster_one: pd.Series, cluster_two: pd.Series, distances: dict, linkage):
    res = []
    for one in cluster_one:
        for two in cluster_two:
            if (one, two) in distances:
                res.append(distances[(one, two)])
    if linkage == 'average':
        return sum(res) / ((len(cluster_one)*len(cluster_two))/2)
    if linkage == 'complete':
        return max(res)
    if linkage == 'single':
        return min(res)


def silhouette_score(ids_cluster: pd.DataFrame, ids_mapping: pd.DataFrame, distances: dict, mapper: Mapper, ids: dict):
    cluster_sizes = ids_cluster['cluster_index'].value_counts().to_dict()
    # ==== calculate inter and intra distance for every entity
    dist = {}
    for index, row in ids_cluster.iterrows():
        id_series = ids_cluster[ids_cluster['cluster_index'] == row['cluster_index']][0]
        id_series = ids_mapping[ids_mapping[ids['id_type']].isin(set(id_series))][ids['att_id']]
        entity = ids_mapping[ids_mapping[ids['id_type']].isin({row[0]})][ids['att_id']]
        dist[row[0]] = {
            'intra_dist': intra_distance(cluster=id_series, mapper=mapper, linkage='average', ids=ids, entity=entity),
            'inter_dist': min_inter_distance_to_entity(entity=entity, distances=distances,
                                                       current_cluster=row['cluster_index'], ids_cluster=ids_cluster,
                                                       ids_mapping=ids_mapping, ids=ids)}
    # ==== calculate score
    s_score = 0
    intra_s_scores = dict()
    for entity in dist:
        current_cluster = ids_cluster[ids_cluster[0] == entity]['cluster_index'].iloc[0]
        if cluster_sizes[current_cluster] > 1 and max(dist[entity]['inter_dist'],
                                                      dist[entity]['intra_dist']) > 0.0:
            score = ((dist[entity]['inter_dist'] - dist[entity]['intra_dist']) /
                     max(dist[entity]['inter_dist'], dist[entity]['intra_dist']))
        else:
            score = 0.0
        # ==== save score for every cluster separately
        if current_cluster not in intra_s_scores:
            intra_s_scores[current_cluster] = 0
        intra_s_scores[current_cluster] += score
        # ==== save for total score
        s_score += score
    for cluster in cluster_sizes:
        intra_s_scores[ids_cluster[ids_cluster['cluster_index'] == cluster][0].iloc[0]] = \
            intra_s_scores[cluster] / cluster_sizes[cluster]
        del intra_s_scores[cluster]
    return s_score / len(distances), intra_s_scores


def dunn_index(ids_cluster: pd.DataFrame, ids_mapping: pd.DataFrame, distances: dict, mapper: Mapper, ids: dict):
    clusters = ids_cluster['cluster_index'].unique()
    max_intra_dist = None
    min_inter_dist = None
    for i in range(len(clusters)):
        id_series = ids_cluster[ids_cluster['cluster_index'] == clusters[i]][0]
        id_series = ids_mapping[ids_mapping[ids['id_type']].isin(set(id_series))][ids['att_id']]
        distance = intra_distance(cluster=id_series, mapper=mapper, linkage='average', ids=ids)
        if max_intra_dist is None or max_intra_dist < distance:
            max_intra_dist = distance
        for j in range(i + 1, len(clusters)):
            id_series_to = ids_cluster[ids_cluster['cluster_index'] == clusters[j]][0]
            id_series_to = ids_mapping[ids_mapping[ids['id_type']].isin(set(id_series_to))][ids['att_id']]
            distance = inter_distance(cluster_one=id_series, cluster_two=id_series_to, distances=distances,
                                      linkage='average')
            if min_inter_dist is None or min_inter_dist > distance:
                min_inter_dist = distance
    return min_inter_dist / max_intra_dist
