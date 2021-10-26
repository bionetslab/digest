#!/usr/bin/python3

import numpy as np


def get_distance_matrix(eval_set):
    dis_mat = np.zeros((len(eval_set), len(eval_set)))
    for index1 in range(0, len(eval_set)):
        for index2 in range(index1, len(eval_set)):
            calc_dis = 1-(len(eval_set[index1] & eval_set[index2]) / min(len(eval_set[index1]), len(eval_set[index2])))
            # assign to matrix
            dis_mat[index1][index2] = calc_dis
            dis_mat[index2][index1] = calc_dis
    return dis_mat


def create_ref_dict(mapping, keys):
    reference_dict = dict()
    for att_type in keys:
        reference_dict[att_type] = set.union(*mapping[att_type])
    return reference_dict


def get_intersection(values_set, ref_set):
    if len(values_set) == 0:
        return 0.0
    return len(values_set & ref_set) / len(values_set)


def evaluate_values(mapping, ref_dict, threshold, keys):
    evaluation = list()
    for attribute in keys:
        evaluated_series = mapping[attribute].apply(get_intersection, ref_set=ref_dict[attribute])
        evaluation.append([attribute, str(len(evaluated_series[evaluated_series > threshold]) / len(evaluated_series))])
    return evaluation
