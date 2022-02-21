#!/usr/bin/python3

import pandas as pd
import numpy as np
from .. import config as c


def preprocess_results(mapping: pd.DataFrame, multicol: str, singlecol: str, key: str):
    """
    Depending on the input id the mapping can map either one value of an attribute
    or multiple values resulting in two columns for an attribute. If a mapping
    has two columns for an attribute, both columns will be combined.

    :param mapping: mapping for an attribute
    :param multicol: name of the column with multiple mapped values
    :param singlecol: name of the column with single mapped values
    :param key: key of the dict in the multiple mapped values
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


def combine_rows_to_set(x):
    if isinstance(x, list):
        return combine_rowsets_list(x)
    elif isinstance(x, set):
        return combine_rowsets_set(x)
    elif isinstance(x, pd.Series):
        return combine_rowsets_series_to_set(x)
    elif isinstance(x, str):
        return string_to_set(x, sep=";")
    return None


def combine_rows_to_string(x):
    if isinstance(x, list):
        return list_to_string(x)
    elif isinstance(x, set):
        return set_to_string(x)
    elif isinstance(x, pd.Series):
        return combine_rowsets_series_to_string(x)
    elif isinstance(x, str):
        return ';'.join(x).split(';')
    return None


def combine_rowsets_list(x: list):
    return set().union(*x)


def combine_rowsets_series_to_set(x: pd.Series):
    if isinstance(x.iloc[0], set):
        return set().union(*x)
    else:
        return set(filter(None, ';'.join(x).split(';')))


def combine_rowsets_series_to_string(x: pd.Series):
    if isinstance(x.iloc[0], set):
        return ";".join(*x)
    else:
        return ';'.join(x)


def combine_rowsets_set(x: set):
    return set().union(x)


def string_to_set(x: str, sep: str = ';'):
    return set(filter(None, x.split(sep)))


def set_to_string(x: set, sep: str = ';'):
    return sep.join(x)


def list_to_string(x, sep: str = ';'):
    if isinstance(x, list):
        return sep.join(x)
    else:
        return x


def set_to_len(x: set):
    return len(x)


def map_to_prev_id(main_id_type: str, id_type: str, id_mapping: pd.DataFrame, att_mapping: pd.DataFrame):
    """
    Map attribute mapping back to original id.

    :param main_id_type: main id type, here mondo for diseases or entrez for genes
    :param id_type: target id type of user input
    :param id_mapping: full mapping of all id types of either genes or diseases
    :param att_mapping: full attribute mapping from terms to main id
    :return: return attribute mapping mapped to target id type
    """
    if id_type not in ['entrezgene', 'mondo']:
        id_mapping = id_mapping.explode(id_type)
    attributes = c.DISEASE_ATTRIBUTES_KEY if main_id_type == "mondo" else c.GENE_ATTRIBUTES_KEY
    columns = [main_id_type, id_type] if id_type != main_id_type else [main_id_type]
    mapping_subset = id_mapping[columns].drop_duplicates()
    hit_mapping = pd.merge(mapping_subset, att_mapping, on=[main_id_type], how='outer')
    hit_mapping = hit_mapping.drop(columns=[main_id_type]) if id_type != main_id_type else hit_mapping
    hit_mapping = hit_mapping.fillna('')
    hit_mapping = hit_mapping[hit_mapping[id_type] != ""]
    hit_mapping = hit_mapping.groupby(id_type, as_index=False).agg(
        {x: combine_rows_to_set for x in attributes})
    return hit_mapping


def transform_disgenet_mapping(mapping: pd.DataFrame, file: str, col_old, col_new):
    """
    Transform mapping from disgenet database to create one combined dataframe of attributes to mondo id.

    :param mapping: disease id mapping from disgenet
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
    df = df[['mondo', col_new]].fillna('').groupby(['mondo'], as_index=False).agg(combine_rowsets_series_to_set)
    return df
