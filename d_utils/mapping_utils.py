#!/usr/bin/python3

import pandas as pd
import numpy as np


def preprocess_results(mapping: pd.DataFrame, multicol: str, singlecol: str, key: str, explode=False):
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

    def convert_to_string(cell):
        if str(cell) != 'nan':
            extracted_ids = [val.get(key) for val in cell if key in val]
            return ';'.join(str(e) for e in list(set(extracted_ids)))
        return cell

    mapping[multicol] = mapping[multicol].apply(lambda x: convert_to_string(x)) if multicol in mapping else np.nan
    if singlecol in mapping:
        mapping[multicol].fillna(mapping[singlecol], inplace=True)
        mapping = mapping.drop(columns=[singlecol])
    if explode:
        mapping[multicol] = mapping[multicol].str.split(';').explode(multicol)
        mapping.rename(columns={multicol: singlecol}, inplace=True)
    return mapping


def split_and_expand_column(data: pd.DataFrame, split_string: str, column_name: str):
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


def combine_rows(x: str):
    return set(filter(None, ';'.join(x).split(';')))


def combine_rowsets(x: list):
    return set().union(*x)


def combine_rowsets_list(x: set):
    return set().union(x)


def string_to_set(x: str, sep: str = ';'):
    return set(filter(None, x.split(sep)))


def set_to_string(x: set, sep: str = ';'):
    return sep.join(x)


def transform_disgenet_mapping(mapping: pd.DataFrame, file: str, col_old, col_new):
    """
    Transform mapping from disgenet database to create one combined dataframe of attributes to mondo id.

    :param mapping:
    :param file: path to file with disease mapping from disgenet
    :param col_old: attribute column name in raw disgenet file
    :param col_new: desired new column name for col_old
    :return: combined dataframe
    """
    disease_mapping = pd.read_csv(file, compression='gzip', sep='\t', dtype=str)
    df = pd.merge(mapping[['diseaseId', 'mondo']], disease_mapping[['diseaseId', col_old]],
                  on="diseaseId", how="left")
    df = df.rename(columns={col_old: col_new})
    df[col_new] = df[col_new].str.strip()
    df = df[['mondo', col_new]].fillna('').groupby(['mondo'], as_index=False).agg(combine_rows)
    return df
