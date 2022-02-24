#!/usr/bin/python3

import pandas as pd
from collections import defaultdict


def precalc_distance_dicts(ids_cluster: pd.DataFrame, ids_mapping: pd.DataFrame, distances: dict, index_to_id: dict,
                           ids: dict):
    """
    Precalculate intra and inter distances based on pairwise distances calculated beforehand.

    :param ids_cluster: mapping of ids to corresponding cluster
    :param ids_mapping: mapping of different id types to original id type
    :param distances: all pairwise distances previously calculated
    :param index_to_id: dictionary with mapping of distance matrix indices to ids
    :param ids: dictionary with mapping of keys to mapping information in mapper
    :return: 4 dictionaries consisting of entity intra, entity inter, overall intra and overall inter distances
    """
    # ===== create empty default dicts =====
    entity_intra = defaultdict(lambda: {'max': None, 'sum': 0, 'min': None, 'count': 0})
    entity_inter = defaultdict(lambda: defaultdict(lambda: {'max': None, 'sum': 0, 'min': None, 'count': 0}))
    intra = defaultdict(lambda: {'max': None, 'sum': 0, 'min': None, 'count': 0})
    inter = defaultdict(lambda: defaultdict(lambda: {'max': None, 'sum': 0, 'min': None, 'count': 0}))
    # ===== map ids to cluster =====
    id_to_cluster = ids_cluster.set_index('id').to_dict()['cluster_index']
    # ===== map att ids to ids =====
    if ids['att_id'] != ids['id_type']:
        att_id_to_id = \
            ids_mapping[[ids['att_id'], ids['id_type']]].groupby(ids['att_id']).agg(lambda g: set(g)).to_dict()[
                ids['id_type']]
    else:
        att_id_to_id = None

    # ===== method to add values =====
    def add_value(destination, distance):
        dist = 1 - distance
        if destination['max'] is None or destination['max'] < dist:
            destination['max'] = dist
        if destination['min'] is None or destination['min'] > dist:
            destination['min'] = dist
        destination['sum'] = destination['sum'] + dist
        destination['count'] = destination['count'] + 1
        return destination

    # ===== assign distances > 0.0  to dicts =====
    for index1, index2 in distances:
        att_id1 = index_to_id[index1]
        att_id2 = index_to_id[index2]
        att_id_to_id1 = att_id_to_id[att_id1] if att_id_to_id is not None else {att_id1}
        att_id_to_id2 = att_id_to_id[att_id2] if att_id_to_id is not None else {att_id2}
        for id1 in att_id_to_id1:
            for id2 in att_id_to_id2:
                id1_cluster = id_to_cluster[id1]
                id2_cluster = id_to_cluster[id2]
                if id1_cluster == id2_cluster:
                    entity_intra[id1] = add_value(destination=entity_intra[id1], distance=distances[(index1, index2)])
                    entity_intra[id2] = add_value(destination=entity_intra[id2], distance=distances[(index1, index2)])
                    intra[id1_cluster] = add_value(destination=intra[id1_cluster], distance=distances[(index1, index2)])
                else:
                    entity_inter[id1][id2_cluster] = add_value(destination=entity_inter[id1][id2_cluster],
                                                               distance=distances[(index1, index2)])
                    entity_inter[id2][id1_cluster] = add_value(destination=entity_inter[id2][id1_cluster],
                                                               distance=distances[(index1, index2)])
                    inter[id1_cluster][id2_cluster] = add_value(destination=inter[id1_cluster][id2_cluster],
                                                                distance=distances[(index1, index2)])
                    inter[id2_cluster][id1_cluster] = add_value(destination=inter[id2_cluster][id1_cluster],
                                                                distance=distances[(index1, index2)])
    return {'entity_intra': entity_intra, 'entity_inter': entity_inter, 'intra': intra, 'inter': inter}


def calc_linkage(value_dict: dict, size: int, linkage="average"):
    """
    Calculate desired linkage based on precalculated values.

    :param value_dict: dictionary with precalculated values
    :param size: size of cluster
    :param linkage: linkage to be calculated [Default="average"]
    :return: linkage value
    """
    if linkage == "average":
        return (value_dict['sum'] + (size - value_dict['count'])) / size
    if linkage == "complete":
        return value_dict['max']
    if linkage == "single":
        return value_dict['min']
    else:
        return None


def silhouette_score(ids_cluster: pd.DataFrame, distances: dict, linkage="average"):
    """
    Calculate the silhouette score for the given cluster.

    :param ids_cluster: mapping of ids to corresponding cluster
    :param distances: all pairwise distances previously calculated
    :param linkage: linkage type for intra and inter distance [Default="average"]
    :return: sillhouette score for the whole clustering and for each cluster separately
    """
    cluster_sizes = ids_cluster['cluster_index'].value_counts().to_dict()
    # ===== map ids to cluster =====
    id_to_cluster = ids_cluster.set_index('id').to_dict()['cluster_index']
    # ===== calculate score =====
    s_score = 0
    intra_s_scores = dict()
    for entity in set(distances['entity_intra'].keys()).union(set(distances['entity_inter'].keys())):
        current_cluster = id_to_cluster[entity]
        # ===== calc intra distance =====
        if entity in distances['entity_intra']:
            entity_intra = calc_linkage(value_dict=distances['entity_intra'][entity],
                                        size=cluster_sizes[current_cluster], linkage=linkage)
        else:
            entity_intra = 1
        # ===== calc min inter distance =====
        min_entity_inter = None
        if entity in distances['entity_inter'] and len(distances['entity_inter'][entity]) < (len(cluster_sizes) - 1):
            for cluster in distances['entity_inter'][entity]:
                distance = calc_linkage(value_dict=distances['entity_inter'][entity][cluster],
                                        size=cluster_sizes[cluster], linkage=linkage)
                if min_entity_inter is None or min_entity_inter > distance:
                    min_entity_inter = distance
        if min_entity_inter is None:
            min_entity_inter = 1

        if cluster_sizes[current_cluster] > 1 and max(min_entity_inter, entity_intra) > 0.0:
            score = ((min_entity_inter - entity_intra) / max(min_entity_inter, entity_intra))
        else:
            score = 0.0
        # ===== save score for every cluster separately =====
        if current_cluster not in intra_s_scores:
            intra_s_scores[current_cluster] = 0
        intra_s_scores[current_cluster] += score
        # ===== save for total score =====
        s_score += score
    for cluster in cluster_sizes:
        if cluster in intra_s_scores:
            intra_s_scores[cluster] = intra_s_scores[cluster] / cluster_sizes[cluster]
        else:
            intra_s_scores[cluster] = 0
    return s_score / len(ids_cluster['id']), intra_s_scores


def dunn_index(ids_cluster: pd.DataFrame, distances: dict, linkage="average") -> float:
    """
    Calculate the dunn index for the given cluster.

    :param ids_cluster: mapping of ids to corresponding cluster
    :param distances: all pairwise distances previously calculated
    :param linkage: linkage type for intra and inter distance [Default="average"]
    :return: dunn index score as float
    """
    max_intra_dist = 0
    min_inter_dist = None

    # ===== count cluster size =====
    cluster_to_size = ids_cluster['cluster_index'].value_counts().to_dict()
    for cluster in cluster_to_size:
        # ===== calc intra distance =====
        if cluster in distances['intra']:
            distance = calc_linkage(value_dict=distances['intra'][cluster],
                                    size=(cluster_to_size[cluster] * cluster_to_size[cluster]) / 2, linkage=linkage)
        else:  # if no distances saved for intra cluster: all pairwise have distance 1
            distance = 1
        if max_intra_dist is None or max_intra_dist < distance:
            max_intra_dist = distance
        # ===== calc min inter distance =====
        if cluster in distances['inter']:
            # ===== distance to at least one cluster == 1 -> not in dict =====
            if len(distances['inter'][cluster]) < (len(cluster_to_size) - 1):
                min_inter_dist = 1
            # ===== calc distance to all clusters =====
            for to_cluster in distances['inter'][cluster]:
                distance = calc_linkage(value_dict=distances['inter'][cluster][to_cluster],
                                        size=cluster_to_size[cluster] * cluster_to_size[to_cluster], linkage=linkage)
                if min_inter_dist is None or min_inter_dist > distance:
                    min_inter_dist = distance
    if max_intra_dist == 0 or min_inter_dist is None:
        return 0.0
    return min_inter_dist / max_intra_dist


def davies_bouldin_index(ids_cluster: pd.DataFrame, distances: dict, linkage="average") -> float:
    """
    Calculate the davies bouldin index for the given cluster.

    :param ids_cluster: mapping of ids to corresponding cluster
    :param distances: all pairwise distances previously calculated
    :param linkage: linkage type for intra and inter distance [Default="average"]
    :return: davies bouldin index score as float
    """
    max_values = 0
    # ===== count cluster size =====
    cluster_to_size = ids_cluster['cluster_index'].value_counts().to_dict()
    for cluster_i in cluster_to_size:
        cur_max = 0
        for cluster_j in cluster_to_size:
            if cluster_i != cluster_j:
                # ===== calc intra distances =====
                distance_i = calc_linkage(value_dict=distances['intra'][cluster_i],
                                          size=(cluster_to_size[cluster_i] * cluster_to_size[cluster_i]) / 2,
                                          linkage=linkage) if cluster_i in distances['intra'] else 1
                distance_j = calc_linkage(value_dict=distances['intra'][cluster_j],
                                          size=(cluster_to_size[cluster_j] * cluster_to_size[cluster_j]) / 2,
                                          linkage=linkage) if cluster_j in distances['intra'] else 1
                # ===== calc inter distance =====
                distance_ij = calc_linkage(value_dict=distances['inter'][cluster_i][cluster_j],
                                           size=cluster_to_size[cluster_i] * cluster_to_size[cluster_j],
                                           linkage=linkage)
                # ===== calc (q(ci)*q(cj))/p(ci,cj) =====
                value = (distance_i * distance_j) / distance_ij if distance_ij != 0 else 0
                if value > cur_max:
                    cur_max = value
        max_values = max_values + cur_max
    return max_values / len(cluster_to_size)
