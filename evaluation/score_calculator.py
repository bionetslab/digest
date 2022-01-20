#!/usr/bin/python3
import pandas as pd
from mappers.mapper import Mapper
from collections import defaultdict


def precalc_distance_dicts(ids_cluster: pd.DataFrame, ids_mapping: pd.DataFrame, distances: dict, index_to_id: dict,
                           ids: dict):
    entity_intra = defaultdict(lambda: list())
    entity_inter = defaultdict(lambda: defaultdict(lambda: list()))
    intra = defaultdict(lambda: list())
    inter = defaultdict(lambda: defaultdict(lambda: list()))
    # map ids to cluster
    id_to_cluster = ids_cluster.set_index(0).to_dict()['cluster_index']
    # map att ids to ids
    att_id_to_id = dict()
    for att_id in ids_mapping[ids['att_id']].unique():
        att_id_to_id[att_id] = ids_mapping[ids_mapping[ids['att_id']] == att_id][ids['id_type']]
    # assign distances > 0.0  to dicts
    for index1, index2 in distances:
        att_id1 = index_to_id[index1]
        att_id2 = index_to_id[index2]
        for id1 in att_id_to_id[att_id1]:
            for id2 in att_id_to_id[att_id2]:
                id1_cluster = id_to_cluster[id1]
                id2_cluster = id_to_cluster[id2]
                if id1_cluster == id2_cluster:
                    # if id1 not in entity_intra:
                    #     entity_intra[id1] = list()
                    entity_intra[id1].append(distances[(index1, index2)])
                    # if id2 not in entity_intra:
                    #     entity_intra[id2] = list()
                    entity_intra[id2].append(distances[(index1, index2)])
                    # if id1_cluster not in intra:
                    #     intra[id1_cluster] = list()
                    intra[id1_cluster].append(distances[(index1, index2)])
                else:
                    # if id1 not in entity_inter:
                    #     entity_inter[id1] = dict()
                    # if id2_cluster not in entity_inter[id1]:
                    #     entity_inter[id1][id2_cluster] = list()
                    entity_inter[id1][id2_cluster].append(distances[(index1, index2)])
                    # if id2 not in entity_inter:
                    #     entity_inter[id2] = dict()
                    # if id1_cluster not in entity_inter[id2]:
                    #     entity_inter[id2][id1_cluster] = list()
                    entity_inter[id2][id1_cluster].append(distances[(index1, index2)])
                    # if id1_cluster not in inter:
                    #     inter[id1_cluster] = dict()
                    # if id2_cluster not in inter[id1_cluster]:
                    #     inter[id1_cluster][id2_cluster] = list()
                    inter[id1_cluster][id2_cluster].append(distances[(index1, index2)])
                    # if id2_cluster not in inter:
                    #     inter[id2_cluster] = dict()
                    # if id1_cluster not in inter[id2_cluster]:
                    #     inter[id2_cluster][id1_cluster] = list()
                    inter[id2_cluster][id1_cluster].append(distances[(index1, index2)])

    # # calculate average
    # def calc_linkage(value_list, size, linkage="average"):
    #     if linkage == "average":
    #         return sum(value_list) / size
    #     if linkage == "complete":
    #         return max(value_list)
    #     if linkage == "single":
    #         return min(value_list)
    #     else:
    #         return None
    #
    # entity_intra = {k: calc_linkage(value_list=v, size=cluster_to_size[id_to_cluster[k]], linkage="average") for k, v in
    #                 entity_intra.items()}
    # intra = {k: calc_linkage(value_list=v, size=cluster_to_size[k], linkage="average") for k, v in intra.items()}
    # for entity in entity_inter:
    #     entity_inter[entity] = {k: calc_linkage(value_list=v, size=cluster_to_size[k], linkage="average") for k, v in entity.items()}
    # for cluster in inter:
    #     inter[cluster] = {k: calc_linkage(value_list=v, size=cluster_to_size[k], linkage="average") for k, v in cluster.items()}
    return {'entity_intra': entity_intra, 'entity_inter': entity_inter, 'intra': intra, 'inter': inter}


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
            return (sub_mat_1.sum() + sub_mat_2.sum()) / (len(cluster) * len(entity))
        if linkage == 'complete':
            return max(sub_mat_1.max(), sub_mat_2.max())
    else:
        sub_mat = mapper.get_loaded_distances(in_series=cluster, id_type=ids['sparse_key'],
                                              key=ids['attribute'])
        if linkage == 'average':
            return sub_mat.sum() / (len(cluster) * (len(cluster) - 1)) / 2
        if linkage == 'complete':
            return sub_mat.max()


def inter_distance(cluster_one: pd.Series, cluster_two: pd.Series, distances: dict, linkage):
    res = []
    for one in cluster_one:
        for two in cluster_two:
            if (one, two) in distances:
                res.append(distances[(one, two)])
    if linkage == 'average':
        return sum(res) / ((len(cluster_one) * len(cluster_two)) / 2)
    if linkage == 'complete':
        return max(res)
    if linkage == 'single':
        return min(res)


def calc_linkage(value_list, size, linkage="average"):
    if linkage == "average":
        return sum(value_list) / size
    if linkage == "complete":
        return max(value_list)
    if linkage == "single":
        return min(value_list)
    else:
        return None


def silhouette_score(ids_cluster: pd.DataFrame, ids_mapping: pd.DataFrame, distances: dict, linkage="average"):
    cluster_sizes = ids_cluster['cluster_index'].value_counts().to_dict()
    # ==== calculate inter and intra distance for every entity
    # dist = {}
    # for index, row in ids_cluster.iterrows():
    #     id_series = ids_cluster[ids_cluster['cluster_index'] == row['cluster_index']][0]
    #     id_series = ids_mapping[ids_mapping[ids['id_type']].isin(set(id_series))][ids['att_id']]
    #     entity = ids_mapping[ids_mapping[ids['id_type']].isin({row[0]})][ids['att_id']]
    #     dist[row[0]] = {
    #         'intra_dist': intra_distance(cluster=id_series, mapper=mapper, linkage='average', ids=ids, entity=entity),
    #         'inter_dist': min_inter_distance_to_entity(entity=entity, distances=distances,
    #                                                    current_cluster=row['cluster_index'], ids_cluster=ids_cluster,
    #                                                    ids_mapping=ids_mapping, ids=ids)}
    # map ids to cluster
    id_to_cluster = ids_cluster.set_index(0).to_dict()['cluster_index']
    # ==== calculate score
    s_score = 0
    intra_s_scores = dict()
    for entity in set(distances['entity_intra'].keys()).union(set(distances['entity_inter'].keys())):
        current_cluster = id_to_cluster[entity]
        # calc intra distance
        entity_intra = calc_linkage(value_list=distances['entity_intra'][entity], size=cluster_sizes[current_cluster], linkage=linkage) if entity in distances['entity_intra'] else 0
        # calc min inter distance
        min_entity_inter = None
        if entity in distances['entity_inter'] and len(distances['entity_inter'][entity]) < (len(cluster_sizes) - 1):
            for cluster in distances['entity_inter'][entity]:
                distance = calc_linkage(value_list=distances['entity_inter'][entity][cluster], size=cluster_sizes[cluster], linkage=linkage)
                if min_entity_inter is None or min_entity_inter > distance:
                    min_entity_inter = distance
        else:
            min_entity_inter = 0

        if cluster_sizes[current_cluster] > 1 and max(min_entity_inter, entity_intra) > 0.0:
            score = ((min_entity_inter - entity_intra) / max(min_entity_inter, entity_intra))
        else:
            score = 0.0
        # ==== save score for every cluster separately
        if current_cluster not in intra_s_scores:
            intra_s_scores[current_cluster] = 0
        intra_s_scores[current_cluster] += score
        # ==== save for total score
        s_score += score
    for cluster in cluster_sizes:
        if cluster in intra_s_scores:
            intra_s_scores[cluster] = intra_s_scores[cluster] / cluster_sizes[cluster]
        else:
            intra_s_scores[cluster] = 0
    return s_score / len(ids_cluster[0]), intra_s_scores


def dunn_index(ids_cluster: pd.DataFrame, distances: dict, linkage="average"):
    max_intra_dist = 0
    min_inter_dist = None
    # for i in range(len(clusters)):
    #     id_series = ids_cluster[ids_cluster['cluster_index'] == clusters[i]][0]
    #     id_series = ids_mapping[ids_mapping[ids['id_type']].isin(set(id_series))][ids['att_id']]
    #     distance = intra_distance(cluster=id_series, mapper=mapper, linkage='average', ids=ids)
    #     if max_intra_dist is None or max_intra_dist < distance:
    #         max_intra_dist = distance
    #     for j in range(i + 1, len(clusters)):
    #         id_series_to = ids_cluster[ids_cluster['cluster_index'] == clusters[j]][0]
    #         id_series_to = ids_mapping[ids_mapping[ids['id_type']].isin(set(id_series_to))][ids['att_id']]
    #         distance = inter_distance(cluster_one=id_series, cluster_two=id_series_to, distances=distances,
    #                                   linkage='average')
    #         if min_inter_dist is None or min_inter_dist > distance:
    #             min_inter_dist = distance

    # count cluster size
    cluster_to_size = ids_cluster['cluster_index'].value_counts().to_dict()
    for cluster in cluster_to_size:
        if cluster in distances['intra']:
            distance = calc_linkage(value_list=distances['intra'][cluster], size=cluster_to_size[cluster],
                                    linkage=linkage)
            if max_intra_dist is None or max_intra_dist < distance:
                max_intra_dist = distance
        if cluster in distances['inter']:
            # distance to at least one cluster == 0 -> not in dict
            if len(distances['inter'][cluster]) < (len(cluster_to_size) - 1):
                min_inter_dist = 0
            # has distance to all other clusters > 0
            else:
                for to_cluster in distances['inter'][cluster]:
                    distance = calc_linkage(value_list=distances['inter'][cluster][to_cluster],
                                            size=cluster_to_size[to_cluster], linkage=linkage)
                    if min_inter_dist is None or min_inter_dist > distance:
                        min_inter_dist = distance
        else:
            min_inter_dist = 0
    return min_inter_dist / max_intra_dist
