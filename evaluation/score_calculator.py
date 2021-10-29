#!/usr/bin/python3


def min_inter_distance_to_entity(entity, distance_matrix, current_cluster, ids_cluster):
    min_distance = None
    for cluster in ids_cluster['cluster_index'].unique():
        if cluster != current_cluster:
            cluster_member = ids_cluster[ids_cluster['cluster_index'] == cluster][0]
            distance = inter_distance([entity], cluster_member, distance_matrix, 'average')
            if min_distance is None or min_distance > distance:
                min_distance = distance
    return min_distance


def intra_distance(cluster, distance_matrix, linkage, entity=None):
    distances = []
    if len(cluster) == 1:
        return 0
    for i in range(len(cluster)):
        if entity is not None:
            if cluster[i] != entity:
                if entity in distance_matrix.columns and cluster[i] in distance_matrix.columns:
                    distances.append(distance_matrix[entity][cluster[i]])
        else:
            for j in range(i + 1, len(cluster)):
                if cluster[j] in distance_matrix.columns and cluster[i] in distance_matrix.columns:
                    distances.append(distance_matrix[cluster[j]][cluster[i]])
    calculations = len(cluster) * (len(cluster) - 1) if entity is None else len(cluster) - 1
    if linkage == 'average':
        return sum(distances) / calculations
    if linkage == 'complete':
        return max(distances)


def inter_distance(cluster_one, cluster_two, distance_matrix, linkage):
    distances = []
    for one in cluster_one:
        for two in cluster_two:
            if one in distance_matrix.columns and two in distance_matrix.columns:
                distances.append(distance_matrix[one][two])
    if linkage == 'average':
        return sum(distances) / (len(cluster_one) * len(cluster_two))
    if linkage == 'complete':
        return max(distances)
    if linkage == 'single':
        return min(distances)


def silhouette_score(distance_matrix, ids_cluster):
    cluster_sizes = ids_cluster['cluster_index'].value_counts().to_dict()
    # ==== calculate inter and intra distance for every entity
    distances = {}
    for entity in ids_cluster[0]:
        current_cluster = ids_cluster[ids_cluster[0] == entity]['cluster_index'].iloc[0]
        distances[entity] = {
            'intra_dist': intra_distance(list(ids_cluster[ids_cluster['cluster_index'] == current_cluster][0]),
                                         distance_matrix, 'average', entity),
            'inter_dist': min_inter_distance_to_entity(entity, distance_matrix, current_cluster, ids_cluster)}
    s_score = 0
    intra_s_scores = dict()
    for entity in distances:
        current_cluster = ids_cluster[ids_cluster[0] == entity]['cluster_index'].iloc[0]
        score = ((distances[entity]['inter_dist'] - distances[entity]['intra_dist']) /
                 max(distances[entity]['inter_dist'], distances[entity]['intra_dist'])) if cluster_sizes[
                                                                                               current_cluster] > 1 else 0
        # ==== save score for every cluster separately
        if current_cluster not in intra_s_scores:
            intra_s_scores[current_cluster] = 0
        intra_s_scores[current_cluster] += score
        # ==== save for total score
        s_score += score
    for cluster in cluster_sizes:
        intra_s_scores[ids_cluster[ids_cluster['cluster_index'] == cluster][1].iloc[0]] = \
        intra_s_scores[cluster] / cluster_sizes[cluster]
        del intra_s_scores[cluster]
    return s_score / len(distances), intra_s_scores


def dunn_index(distance_matrix, ids_cluster):
    clusters = ids_cluster['cluster_index'].unique()
    max_intra_dist = None
    min_inter_dist = None
    for i in range(len(clusters)):
        distance = intra_distance(list(ids_cluster[ids_cluster['cluster_index'] == clusters[i]][0]),
                                  distance_matrix, 'average')
        if max_intra_dist is None or max_intra_dist < distance:
            max_intra_dist = distance

        for j in range(i + 1, len(clusters)):
            distance = inter_distance(ids_cluster[ids_cluster['cluster_index'] == clusters[i]][0],
                                      ids_cluster[ids_cluster['cluster_index'] == clusters[j]][0],
                                      distance_matrix, 'average')
            if min_inter_dist is None or min_inter_dist > distance:
                min_inter_dist = distance
    return min_inter_dist / max_intra_dist
