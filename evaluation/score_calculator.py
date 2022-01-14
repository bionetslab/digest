#!/usr/bin/python3
import pandas as pd


def min_inter_distance_to_entity(entity: pd.Series, comparator: dict, current_cluster, ids_cluster: pd.DataFrame,
                                 ids_mapping: pd.DataFrame = None):
    min_distance = None
    for cluster in ids_cluster['cluster_index'].unique():
        if cluster != current_cluster:
            cluster_member = ids_cluster[ids_cluster['cluster_index'] == cluster][0]
            if ids_mapping is not None:  # remap ids for use of distance sparse
                cluster_member = ids_mapping[ids_mapping[0].isin(set(cluster_member))][1]
            distance = inter_distance(cluster_one=entity, cluster_two=cluster_member, comparator=comparator,
                                      linkage='average')
            if min_distance is None or min_distance > distance:
                min_distance = distance
    return min_distance


def intra_distance(cluster: pd.Series, comparator: dict, linkage, entity: set = None):
    distances = []
    if len(cluster) == 1:
        return 0
    if entity is not None:
        for element in cluster:
            if element not in entity:
                for entity_el in entity:
                    distances.append(
                        comparator['mapper'].get_loaded_distance(id1=element, id2=entity_el,
                                                                 id_type=comparator['sparse_key'],
                                                                 key=comparator['att_id']))
    else:
        distances, _ = comparator['mapper'].get_loaded_distances(in_series=cluster, id_type=comparator['sparse_key'],
                                                                 key=comparator['att_id'])
    if linkage == 'average':
        return sum(distances) / len(distances)
    if linkage == 'complete':
        return max(distances)


def inter_distance(cluster_one: pd.Series, cluster_two: pd.Series, comparator: dict, linkage):
    distances = []
    for one in cluster_one:
        for two in cluster_two:
            distances.append(
                comparator['mapper'].get_loaded_distance(id1=one, id2=two, id_type=comparator['sparse_key'],
                                                         key=comparator['att_id']))
    if linkage == 'average':
        return sum(distances) / len(distances)
    if linkage == 'complete':
        return max(distances)
    if linkage == 'single':
        return min(distances)


def silhouette_score(ids_cluster, ids_mapping: pd.DataFrame, comparator: dict):
    cluster_sizes = ids_cluster['cluster_index'].value_counts().to_dict()
    # ==== calculate inter and intra distance for every entity
    distances = {}
    for index, row in ids_cluster.iterrows():
        id_series = ids_cluster[ids_cluster['cluster_index'] == row['cluster_index']][0]
        entity = {row[0]}
        if ids_mapping is not None:  # remap ids for use of distance sparse
            id_series = ids_mapping[ids_mapping[0].isin(set(id_series))][1]
            entity = ids_mapping[ids_mapping[0].isin(entity)][1]
        distances[row[0]] = {
            'intra_dist': intra_distance(cluster=id_series, comparator=comparator, linkage='average', entity=entity),
            'inter_dist': min_inter_distance_to_entity(entity=entity, comparator=comparator,
                                                       current_cluster=row['cluster_index'], ids_cluster=ids_cluster,
                                                       ids_mapping=ids_mapping)}
    # ==== calculate score
    s_score = 0
    intra_s_scores = dict()
    for entity in distances:
        current_cluster = ids_cluster[ids_cluster[0] == entity]['cluster_index'].iloc[0]
        if cluster_sizes[current_cluster] > 1 and max(distances[entity]['inter_dist'],
                                                      distances[entity]['intra_dist']) > 0.0:
            score = ((distances[entity]['inter_dist'] - distances[entity]['intra_dist']) /
                     max(distances[entity]['inter_dist'], distances[entity]['intra_dist']))
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


def dunn_index(ids_cluster: pd.DataFrame, ids_mapping: pd.DataFrame, comparator: dict):
    clusters = ids_cluster['cluster_index'].unique()
    max_intra_dist = None
    min_inter_dist = None
    for i in range(len(clusters)):
        id_series = ids_cluster[ids_cluster['cluster_index'] == clusters[i]][0]
        if ids_mapping is not None:  # remap ids for use of distance sparse
            id_series = ids_mapping[ids_mapping[0].isin(set(id_series))][1]
        distance = intra_distance(cluster=id_series, comparator=comparator, linkage='average')
        if max_intra_dist is None or max_intra_dist < distance:
            max_intra_dist = distance
        for j in range(i + 1, len(clusters)):
            id_series_to = ids_cluster[ids_cluster['cluster_index'] == clusters[j]][0]
            if ids_mapping is not None:  # remap ids for use of distance sparse
                id_series_to = ids_mapping[ids_mapping[0].isin(set(id_series_to))][1]
            distance = inter_distance(cluster_one=id_series, cluster_two=id_series_to, comparator=comparator,
                                      linkage='average')
            if min_inter_dist is None or min_inter_dist > distance:
                min_inter_dist = distance
    return min_inter_dist / max_intra_dist
