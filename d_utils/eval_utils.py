#!/usr/bin/python3

import numpy as np
import pandas as pd
import config
import scipy.sparse as sp


def get_distance_matrix(full_att_series: pd.Series, from_ids: pd.Series, id_to_index: dict, to_ids: pd.Series = None,
                        coefficient='jaccard') -> sp.coo_matrix:
    """
    Calculating the distance of each element in from_ids to to_ids based on coefficient.

    :param full_att_series: dataframe with 2 columns (id, attribute values)
    :param from_ids: dataframe with 2 columns (id, attribute values)
    :param to_ids: dataframe with 2 columns (id, attribute values)
    :param coefficient: coefficient type for the distance. Possible: jaccard or overlap [Default="jaccard"]
    :return: distance sparse matrix
    """
    def get_distance(index1: int, index2: int):
        if coefficient == "jaccard":
            return jaccard_coefficient(tar_att_set=full_att_series[index1], ref_att_set=full_att_series[index2])
        else:  # coefficient == "overlap"
            return overlap_coefficient(tar_att_set=full_att_series[index1], ref_att_set=full_att_series[index2])

    row, col, data = list(), list(), list()
    if to_ids is None:
        from_ids = from_ids.to_numpy()
        for id1_index in range(0, len(from_ids) - 1):
            for id2_index in range(id1_index + 1, len(from_ids)):
                calc_dis = get_distance(index1=id_to_index[from_ids[id1_index]],
                                        index2=id_to_index[from_ids[id2_index]])
                # assign to matrix
                if calc_dis > 0.0:
                    row.append(id_to_index[from_ids[id1_index]])
                    col.append(id_to_index[from_ids[id2_index]])
                    data.append(calc_dis)

    else:
        for id1 in from_ids:
            for id2 in to_ids:
                calc_dis = get_distance(index1=id_to_index[id1], index2=id_to_index[id2])
                # assign to matrix
                if calc_dis > 0.0:
                    row.append(id_to_index[id1])
                    col.append(id_to_index[id2])
                    data.append(calc_dis)
    return sp.coo_matrix((np.array(data), (np.array(row), np.array(col))),
                         shape=(len(full_att_series), len(full_att_series)))


def create_ref_dict(mapping, keys: set, enriched=False):
    """
    Create reference dictionary with each attribute type as key
    and the union of all attribute values in the set.

    :param mapping: mapping of reference to attributes
    :param keys: attribute names
    :return: reference dictionary with unified values
    """
    reference_dict = dict()
    if enriched:
        for key in keys:
            if key in mapping:
                reference_dict[config.ENRICH_KEY[key]] = set(mapping[key].dropna())
            else:
                reference_dict[config.ENRICH_KEY[key]] = set()
    else:
        for att_type in keys:
            reference_dict[att_type] = set.union(*mapping[att_type])
    return reference_dict


def overlap_coefficient(tar_att_set, ref_att_set):
    """
    Calculate overlap coefficient by dividing the length of overlapping elements
    of two sets by the minimum length of the two sets.

    :param tar_att_set: target set of attribute values
    :param ref_att_set: reference set of attribute values
    :return: overlap coefficient
    """
    if len(tar_att_set) == 0 & len(ref_att_set) == 0:
        return 0.0
    intersection = len(tar_att_set.intersection(ref_att_set))
    if intersection == 0:
        return 0.0
    return intersection / min(len(tar_att_set), len(ref_att_set))


def jaccard_coefficient(tar_att_set, ref_att_set):
    """
    Calculate jaccard coefficient by dividing the length of overlapping elements
    of two sets by the combined length of the two sets.

    :param tar_att_set: target set of attribute values
    :param ref_att_set: reference set of attribute values
    :return: jaccard coefficient
    """
    if len(tar_att_set) == 0 & len(ref_att_set) == 0:
        return 0.0
    intersection = len(tar_att_set.intersection(ref_att_set))
    if intersection == 0:
        return 0.0
    return intersection / len(tar_att_set.union(ref_att_set))


def evaluate_values(mapping, ref_dict, threshold, keys, coefficient="jaccard"):
    """
    Evaluate mapped attribute values of target set with the unified
    attribute values of the reference based on a threshold.

    :param mapping: mapping of target set to attributes
    :param ref_dict: reference dictionary with unified values
    :param threshold: minimum similarity of each element in tar to reference set [0,1]
    :param keys: attribute names
    :param coefficient: type of coefficient to validate two sets [Default="jaccard"]
    :return:
    """
    evaluation = dict()
    for attribute in keys:
        if coefficient == "jaccard":
            evaluated_series = mapping[attribute].apply(jaccard_coefficient, ref_att_set=ref_dict[attribute])
        else:  # == "overlap_coefficient"
            evaluated_series = mapping[attribute].apply(overlap_coefficient, ref_att_set=ref_dict[attribute])
        evaluation[attribute] = str(len(evaluated_series[evaluated_series > threshold]) / len(evaluated_series))
    return evaluation


def calc_pvalue(test_value, value_df, maximize=True):
    pvalue = dict()
    print(test_value)
    print(value_df)
    for keys in test_value:
        pvalue[keys] = (1 + sum(value_df[keys] <= test_value[keys])) / (len(value_df.index) + 1) if maximize else \
            (1 + sum(value_df[keys] >= test_value[keys])) / (len(value_df.index) + 1)
    return pvalue
