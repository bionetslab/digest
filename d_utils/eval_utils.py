#!/usr/bin/python3

import numpy as np


def get_distance_matrix(eval_df):
    """
    Calculating the distance of each element in id column to each other based on
    the number of shared attribute values divided by the smaller set of values.

    :param eval_df: dataframe with 2 columns (id, attribute values)
    :return: distance matrix
    """
    dis_mat = np.zeros((len(eval_df), len(eval_df)))
    for index1 in range(0, len(eval_df)):
        for index2 in range(index1, len(eval_df)):
            calc_dis = 1 - (len(eval_df[index1] & eval_df[index2]) / min(len(eval_df[index1]), len(eval_df[index2])))
            # assign to matrix
            dis_mat[index1][index2] = calc_dis
            dis_mat[index2][index1] = calc_dis
    return dis_mat


def create_ref_dict(mapping, keys):
    """
    Create reference dictionary with each attribute type as key
    and the union of all attribute values in the set.

    :param mapping: mapping of reference to attributes
    :param keys: attribute names
    :return: reference dictionary with unified values
    """
    reference_dict = dict()
    for att_type in keys:
        reference_dict[att_type] = set.union(*mapping[att_type])
    return reference_dict


def get_intersection(tar_att_set, ref_att_set):
    """
    Calculate special intersection value by dividing the
    length of overlapping elements of two sets by the length of the target set.

    :param tar_att_set: target set of attribute values
    :param ref_att_set: reference set of attribute values
    :return: intersection value
    """
    if len(tar_att_set) == 0:
        return 0.0
    return len(tar_att_set & ref_att_set) / len(tar_att_set)


def evaluate_values(mapping, ref_dict, threshold, keys):
    """
    Evaluate mapped attribute values of target set with the unified
    attribute values of the reference based on a threshold.

    :param mapping: mapping of target set to attributes
    :param ref_dict: reference dictionary with unified values
    :param threshold: minimum similarity of each element in tar to reference set [0,1]
    :param keys: attribute names
    :return:
    """
    evaluation = dict()
    for attribute in keys:
        evaluated_series = mapping[attribute].apply(get_intersection, ref_att_set=ref_dict[attribute])
        evaluation[attribute] = str(len(evaluated_series[evaluated_series > threshold]) / len(evaluated_series))
    return evaluation


def calc_pvalue(test_value, value_df, maximize=True):
    pvalue = dict()
    for keys in test_value:
        pvalue[keys] = 1 - (sum(value_df[keys] < test_value[keys]) / len(value_df.index)) if maximize else 1 - (
                    sum(value_df[keys] > test_value[keys]) / len(value_df.index))
    return pvalue
