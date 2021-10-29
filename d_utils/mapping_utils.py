#!/usr/bin/python3

import pandas as pd
import config


def preprocess_results(mapping, multicol, singlecol, key, explode=False):
    """
    Depending on the input id the mapping can map either one value of an attribute
    or multiple values resulting in two columns for an attribute. If a mapping
    has two columns for an attribute, both columns will be combined.

    :param mapping: mapping for an attribute
    :param multicol: name of the column with multiple mapped values
    :param singlecol: name of the column with single mapped values
    :param key: key of the dict in the multiple mapped values
    :param explode: (optional) can expand dataframe by giving each value a separate row
    :return: transformed dataframe
    """
    def convert_to_string(cell, key):
        if str(cell) != 'nan':
            extracted_ids = [val.get(key) for val in cell]
            return ';'.join(str(e) for e in list(set(extracted_ids)))
        return cell

    mapping[multicol] = mapping[multicol].apply(lambda x: convert_to_string(x, key)) if multicol in mapping else np.nan
    if singlecol in mapping:
        mapping[multicol].fillna(mapping[singlecol], inplace=True)
        mapping = mapping.drop(columns=[singlecol])
    if explode:
        mapping = mapping[multicol].split(';').explode(multicol)
        mapping.rename(columns={multicol: singlecol}, inplace=True)
    return mapping


def get_prev_mapping(in_set, id_type, file, sep):
    """
    Get previous mapping from file sub setting it to the
    element in in_set. Also returning the list of elements
    from the set that are not present in the previous mapping.

    :param in_set: set of genes or diseases
    :param id_type: id type of set
    :param file: file with previous mapping
    :param sep: separator of values in file
    :return: mapped_set: dataframe with previous mapping of input set or subset;
             missing: list of missing values of the set in the previous mapping;
             mapping: complete mapping of file
    """
    # ===== Get mapping from local mapping file =====
    mapping = pd.read_csv(file, sep=sep, header=0, dtype=str)
    if id_type == "ICD-10":
        mapping = split_and_expand_column(data=mapping, split_string=",", column_name="ICD-10")
    # ==== Map given disease set ====
    id_type = config.ID_TYPE_KEY[id_type] if id_type in config.ID_TYPE_KEY else id_type
    mapped_set = mapping[mapping[id_type].isin(in_set)].copy()
    # ===== Get missing values =====
    missing = list(set(in_set) - set(mapping[id_type]))
    return mapped_set, missing, mapping


def split_and_expand_column(data, split_string, column_name):
    """
    Split column value in data by split_string and expand the dataframe
    to have a separate row for each value in split set.

    :param data: dataframe with data
    :param split_string: separator of values in cell
    :param column_name: column to split each cell of
    :return: expanded data dataframe
    """
    s = data[column_name].str.split(split_string, expand=True).stack()
    i = s.index.get_level_values(0)
    df2 = data.loc[i].copy()
    df2[column_name] = s.values
    return df2


def combine_rows(x):
    return set(filter(None, ';'.join(x).split(';')))


