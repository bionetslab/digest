#!/usr/bin/python3

import pandas as pd
import config


def preprocess_results(mapping, multicol, singlecol, key, explode=False):
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
    s = data[column_name].str.split(split_string, expand=True).stack()
    i = s.index.get_level_values(0)
    df2 = data.loc[i].copy()
    df2[column_name] = s.values
    return df2


def combine_rows(x):
    return set(filter(None, ';'.join(x).split(';')))


