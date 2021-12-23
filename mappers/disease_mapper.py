#!/usr/bin/python3

import pandas as pd
from d_utils import config, mapping_utils as mu
from biothings_client import get_client
from mapper import Mapper


def get_disease_to_attributes(disease_set, id_type, mapper:Mapper):
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
    mondo_set = list(set('MONDO:' + disorder_mapping['mondo']))
    # ===== Get mapping from previous mappings =====
    hit_mapping, missing_hits = mapper.get_loaded_mapping(in_set=mondo_set, id_type='mondo', key='disorder_atts')
    # ==== Get disgenet values ====
    if len(missing_hits) > 0:
        md = get_client("disease")
        mapping = md.getdiseases(missing_hits,
                                 fields=','.join(config.DISEASE_ATTRIBUTES),
                                 species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping.rename(columns={'query': 'mondo'}, inplace=True)
        # transform dataframe to combine single and multiple results
        for attribute in config.DISEASE_ATTRIBUTES_KEY:
            mapping = mu.preprocess_results(mapping=mapping, multicol=attribute,
                                            singlecol=attribute + '.' + config.DISEASE_ATTRIBUTES_KEY[attribute],
                                            key=config.DISEASE_ATTRIBUTES_KEY[attribute])
        mapping.drop(columns=['_id', '_version', 'disgenet._license'], inplace=True)
        mapping.drop(columns=['notfound'], inplace=True) if 'notfound' in mapping.columns else None
        # ==== Get additional pathways from file ====
        mondo_to_pathway = pd.read_csv(config.FILES_DIR + 'mondo_to_pathways.csv')
        mapping = mapping.merge(mondo_to_pathway, on='mondo', how='left')
        #  work with nan float values
        mapping = mapping.fillna('')
        mapping = mapping.astype(str)
        # combine with ctd pathway mapping
        mapping['ctd.pathway_related_to_disease'] = (
                mapping['ctd.pathway_related_to_disease'] + ";" + mapping['pathways']).str.strip(';')
        mapping = mapping.drop(columns=['pathways'])
        mapping = mapping.drop_duplicates()
        # ===== Add results from missing values =====
        mapper.update_mappings(in_df=mapping, key='disorder_atts')
        hit_mapping = pd.concat([hit_mapping, mapping]) if not hit_mapping.empty else mapping
    # ==== Map back to previous ids ====
    hit_mapping["mondo"] = hit_mapping["mondo"].str.replace('MONDO:', '')
    columns = ['mondo', id_type] if id_type != 'mondo' else ['mondo']
    mapping_subset = disorder_mapping[columns].drop_duplicates()
    df = pd.merge(mapping_subset, hit_mapping, on=['mondo'], how='outer')
    df = df.drop(columns=['mondo']) if id_type != 'mondo' else df
    df = df.fillna('').groupby(id_type, as_index=False).agg(
        {x: mu.combine_rows for x in config.DISEASE_ATTRIBUTES_KEY})
    return df
