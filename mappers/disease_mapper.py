#!/usr/bin/python3

import pandas as pd
from d_utils import config, mapping_utils as mu
from biothings_client import get_client
from mapper import Mapper


def get_disease_to_attributes(disease_set, id_type, mapper: Mapper):
    """
    Simple attribute mapper using a local mapping file
    and the myDisease.info database.
    Mapped attributes:
    DisGeNET related genes, DisGeNET related variants, ctd + kegg pathways

    :param disease_set: set of disease ids
    :param id_type: id type of set
    :return:
    """
    # ==== Get Mondo IDs ====
    disorder_mapping, _ = mapper.get_loaded_mapping(in_set=disease_set, id_type=id_type, key='disorder_ids')
    # ===== Get mapping from previous mappings =====
    hit_mapping, missing_hits = mapper.get_loaded_mapping(in_set=set(disorder_mapping['mondo']), id_type='mondo',
                                                          key='disorder_atts')
    missing_hits = ['MONDO:' + x for x in missing_hits]
    # ==== Get att for missing values ====
    if len(missing_hits) > 0:
        mapping = get_attributes_from_database(missing=missing_hits)
        mapping = mapping.fillna('').groupby(id_type, as_index=False).agg(
            {x: mu.combine_rows for x in config.DISEASE_ATTRIBUTES_KEY})
        # ===== Add results from missing values =====
        mapper.update_mappings(in_df=mapping, key='disorder_atts')
        hit_mapping = pd.concat([hit_mapping, mapping]) if not hit_mapping.empty else mapping
    # ==== Map back to previous ids ====
    columns = ['mondo', id_type] if id_type != 'mondo' else ['mondo']
    mapping_subset = disorder_mapping[columns].drop_duplicates()
    hit_mapping = pd.merge(mapping_subset, hit_mapping, on=['mondo'], how='outer')
    hit_mapping = hit_mapping.drop(columns=['mondo']) if id_type != 'mondo' else hit_mapping
    hit_mapping = hit_mapping.fillna('').groupby(id_type, as_index=False).agg(
        {x: mu.combine_rowsets for x in config.DISEASE_ATTRIBUTES_KEY})
    return hit_mapping


def get_attributes_from_database(missing, attributes: list = config.DISEASE_ATTRIBUTES_KEY.keys()):
    md = get_client("disease")
    mapping = md.getdiseases(missing,
                             fields=','.join(attributes),
                             species='human', returnall=False, as_dataframe=True, df_index=False)
    mapping.rename(columns={'query': 'mondo'}, inplace=True)
    # transform dataframe to combine single and multiple results
    for attribute in attributes:
        mapping = mu.preprocess_results(mapping=mapping, multicol=attribute,
                                        singlecol=attribute + '.' + config.DISEASE_ATTRIBUTES_KEY[attribute],
                                        key=config.DISEASE_ATTRIBUTES_KEY[attribute])
    mapping.drop(columns=list(set(mapping.columns[1:]) - set(attributes)), inplace=True)
    mapping = mapping.fillna('')
    mapping = mapping.astype(str)
    mapping["mondo"] = mapping["mondo"].str.replace('MONDO:', '')
    return mapping
