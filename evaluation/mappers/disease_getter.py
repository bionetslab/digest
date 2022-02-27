#!/usr/bin/python3

import pandas as pd
from .. import config
from . import mapping_utils as mu
from biothings_client import get_client
from .mapper import Mapper


def get_disease_to_attributes(disease_set, id_type, mapper: Mapper):
    """
    Attribute mapper using local mapping files generated during setup
    and the myDisease.info database for missing values.
    Mapped attributes:
    DisGeNET related genes, DisGeNET related variants, ctd + kegg pathways

    :param disease_set: set of disease ids
    :param id_type: id type of set
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :return: disease to attribute mapping as dataframe
    """
    # ==== Get Mondo IDs ====
    disorder_mapping, _ = mapper.get_loaded_mapping(in_set=disease_set, id_type=id_type, key='disorder_ids')
    # ===== Return empty Dataframe if IDs were not mappable =====
    if disorder_mapping.empty:
        return disorder_mapping
    # ===== Get mapping from previous mappings =====
    hit_mapping, missing_hits = mapper.get_loaded_mapping(in_set=set(disorder_mapping['mondo']), id_type='mondo',
                                                          key='disorder_atts')
    missing_hits = ['MONDO:' + x for x in missing_hits]
    # ===== Get att for missing values =====
    if len(missing_hits) > 0:
        mapping = get_attributes_from_database(missing=missing_hits)
        mapping = mapping.fillna('').groupby(id_type, as_index=False).agg(
            {x: mu.combine_rows_to_set for x in config.DISEASE_ATTRIBUTES_KEY})
        # ===== Add results from missing values =====
        mapper.update_mappings(in_df=mapping, key='disorder_atts')
        hit_mapping = pd.concat([hit_mapping, mapping]) if not hit_mapping.empty else mapping
    # ===== Map back to previous ids =====
    hit_mapping = mu.map_to_prev_id(main_id_type="mondo", id_type=id_type,
                                    id_mapping=disorder_mapping, att_mapping=hit_mapping)
    return hit_mapping


def get_attributes_from_database(missing: list, attributes: list = config.DISEASE_ATTRIBUTES_KEY.keys()):
    """
    Get mapping from myDisease.info.

    :param missing: list of missing values that should be mapped
    :param attributes: attributes that should be mapped to the missing values
    :return: retrieved mapping as dataframe
    """
    md = get_client("disease")
    mapping = md.getdiseases(missing, fields=','.join(attributes),
                             species='human', returnall=False, as_dataframe=True, df_index=False)
    mapping.rename(columns={'query': 'mondo'}, inplace=True)
    if 'notfound' in mapping:
        mapping = mapping[mapping['notfound'] != True]
    if not mapping.empty:
        # ===== transform dataframe to combine single and multiple results =====
        for attribute in attributes:
            mapping = mu.preprocess_results(mapping=mapping, multicol=attribute,
                                            singlecol=attribute + '.' + config.DISEASE_ATTRIBUTES_KEY[attribute],
                                            key=config.DISEASE_ATTRIBUTES_KEY[attribute])
        mapping.drop(columns=list(set(mapping.columns[1:]) - set(attributes)), inplace=True)
        mapping = mapping.fillna('')
        mapping = mapping.astype(str)
        mapping["mondo"] = mapping["mondo"].str.replace('MONDO:', '')
        mapping[mapping.columns[1:]] = mapping[mapping.columns[1:]].fillna('').applymap(mu.string_to_set)
    return mapping
