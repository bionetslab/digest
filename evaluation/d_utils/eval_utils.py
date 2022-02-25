#!/usr/bin/python3

import numpy as np
import pandas as pd
from .. import config
import scipy.sparse as sp


def get_distance_matrix(full_att_series: pd.Series, from_ids: pd.Series, id_to_index: dict, to_ids: pd.Series = None,
                        coefficient='jaccard') -> sp.coo_matrix:
    """
    Calculating the distance of each element in from_ids to elements in to_ids, if provided, or each
    element in from_ids based on coefficient.

    :param full_att_series: dataframe with 2 columns (id, attribute values)
    :param from_ids: dataframe with 2 columns (id, attribute values)
    :param id_to_index: dict with mapping of id to index inside sparse matrix
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
    # ===== from_ids against from_ids =====
    if to_ids is None:
        from_ids = from_ids.to_numpy()
        for id1_index in range(0, len(from_ids) - 1):
            for id2_index in range(id1_index + 1, len(from_ids)):
                calc_dis = get_distance(index1=id_to_index[from_ids[id1_index]],
                                        index2=id_to_index[from_ids[id2_index]])
                # assign to matrix
                if calc_dis > 0.0:
                    if id_to_index[from_ids[id1_index]] < id_to_index[from_ids[id2_index]]:
                        row.append(id_to_index[from_ids[id1_index]])
                        col.append(id_to_index[from_ids[id2_index]])
                    else:
                        row.append(id_to_index[from_ids[id2_index]])
                        col.append(id_to_index[from_ids[id1_index]])
                    data.append(calc_dis)

    # ===== from_ids against to_ids =====
    else:
        for id1 in from_ids:
            for id2 in to_ids:
                calc_dis = get_distance(index1=id_to_index[id1], index2=id_to_index[id2])
                # assign to matrix
                if calc_dis > 0.0:
                    if id_to_index[id1] < id_to_index[id2]:
                        row.append(id_to_index[id1])
                        col.append(id_to_index[id2])
                    else:
                        row.append(id_to_index[id2])
                        col.append(id_to_index[id1])
                    data.append(calc_dis)
    return sp.coo_matrix((np.array(data), (np.array(row), np.array(col))),
                         shape=(len(full_att_series), len(full_att_series)))


def create_ref_dict(mapping: pd.DataFrame, keys: set, enriched: bool = False):
    """
    Create reference dictionary with each attribute type as key
    and the union of all attribute values in the set.

    :param mapping: mapping of reference to attributes
    :param keys: attribute names
    :param enriched: bool if set resulted from enrichment analysis or not
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


def overlap_coefficient(tar_att_set: set, ref_att_set: set):
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


def jaccard_coefficient(tar_att_set: set, ref_att_set: set):
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


def calc_pvalue(test_value: dict, random_values: pd.DataFrame, maximize=True):
    """
    Calculate pvalue based on the values of the original input and the values of random runs.

    :param test_value: values from the original input as dict with attributes as key
    :param random_values: values from random runs as dataframe and attributes as columns
    :param maximize: true if the goal is to have a high test value or false if low test value [Default=True]
    :return: pvalue
    """
    pvalue = dict()
    for keys in test_value:
        pvalue[keys] = (1 + sum(test_value[keys] <= random_values[keys])) / (
                    len(random_values.index) + 1) if maximize else (1 + sum(
            test_value[keys] >= random_values[keys])) / (len(random_values.index) + 1)
    return pvalue
