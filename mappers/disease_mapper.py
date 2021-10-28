#!/usr/bin/python3

import pandas as pd
from d_utils import config, mapping_utils as mu
from biothings_client import get_client


def get_disease_to_attributes(disease_set, id_type):
    # ==== Get Mondo IDs ====
    disease_id_set, _, _ = mu.get_prev_mapping(in_set=disease_set, id_type=id_type,
                                               file=config.FILES_DIR + "new_disorders.map", sep="\t")
    mondo_set = list(set('MONDO:' + disease_id_set['mondo']))
    # ===== Get mapping from previous mappings =====
    df, missing, prev_mapping = mu.get_prev_mapping(in_set=mondo_set, id_type='mondo',
                                                    file=config.FILES_DIR + 'disease_disgenet_mapping.csv', sep=",")
    # ==== Get disgenet values ====
    if len(missing) > 0:
        md = get_client("disease")
        mapping = md.getdiseases(missing,
                                 fields=','.join(config.DISEASE_ATTRIBUTES),
                                 species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping.rename(columns={'query': 'mondo'}, inplace=True)
        # transform dataframe to combine single and multiple results
        for attribute in config.DISEASE_ATTRIBUTES_KEY:
            mapping = mu.preprocess_results(mapping=mapping, multicol=attribute,
                                            singlecol=attribute + '.' + config.DISEASE_ATTRIBUTES_KEY[attribute],
                                            key=config.DISEASE_ATTRIBUTES_KEY[attribute])
        mapping = mapping.drop(columns=['_id', '_version', 'disgenet._license'])
        # ==== Get additional pathways from file ====
        mondo_to_pathway = pd.read_csv(config.FILES_DIR + 'mondo_to_pathways.csv')
        mapping = mapping.merge(mondo_to_pathway, on='mondo', how='left')
        #  work with nan float values
        mapping = mapping.fillna('')
        mapping = mapping.astype(str)
        # combine with ctd pathway mapping
        mapping.loc[:, 'ctd.pathway_related_to_disease'] = (
                mapping.loc[:, 'ctd.pathway_related_to_disease'] + ";" + mapping.loc[:, 'pathways']).str.strip(';')
        mapping = mapping.drop(columns=['pathways'])
        mapping = mapping.drop_duplicates()
        # ===== Add results from missing values =====
        pd.concat([prev_mapping, mapping]).to_csv(config.FILES_DIR + 'disease_disgenet_mapping.csv', index=False)
        df = pd.concat([df, mapping]).reset_index(drop=True)
    # ==== Map back to previous ids ====
    df.loc[:, "mondo"] = df.loc[:, "mondo"].str.split(':').str[1]
    columns = ['mondo', id_type] if id_type != 'mondo' else ['mondo']
    mapping_subset = disease_id_set[columns].drop_duplicates()
    df = pd.merge(mapping_subset, df, on=['mondo'], how='outer')
    df = df.drop(columns=['mondo']) if id_type != 'mondo' else df
    df = df.fillna('').groupby(id_type, as_index=False).agg(
        {x: mu.combine_rows for x in config.DISEASE_ATTRIBUTES_KEY})
    return df
