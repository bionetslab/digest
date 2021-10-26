#!/usr/bin/python3

import pandas as pd
from digest.d_utils import config, mapping_utils as mu
from biothings_client import get_client


def get_gene_mapping(gene_set, id_type):
    """
    Simple gene ID mapper using a local mapping file
    and the myGene.info database.
    Supported ID types:
    uniprot, ensembl, entrezgene, genesymbol

    :param gene_set: Set of gene ids
    :param id_type: Gene ID type of input gene set
    :return: Dataframe
    """
    # ===== Get mapping from previous mappings =====
    df, missing, prev_mapping = mu.get_prev_mapping(in_set=gene_set, id_type=id_type,
                                                    file=config.FILES_DIR + 'gene_id_mapping.csv',
                                                    sep=",")
    # ===== Get mapping for missing values =====
    if len(missing) > 0:
        mg = get_client("gene")
        mapping = mg.querymany(missing, scopes=config.ID_TYPE_KEY[id_type], fields=','.join(config.GENE_IDS),
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping = mapping.drop(columns=[config.ID_TYPE_KEY[id_type]])
        mapping.rename(columns={'query': config.ID_TYPE_KEY[id_type]}, inplace=True)
        # ===== Split if there are multiple ensembl ids =====
        if 'ensembl' in mapping:
            mapping = mu.preprocess_results(mapping=mapping, multicol='ensembl', singlecol='ensembl.gene', key='gene',
                                            explode=True)
        mapping = mapping.drop(columns=['_id', '_score'])
        # ===== Add results from missing values =====
        pd.concat([prev_mapping, mapping]).to_csv(config.FILES_DIR + 'gene_id_mapping.csv', index=False)
        df = pd.concat([df, mapping]).reset_index(drop=True)
    return df


def get_gene_to_attributes(gene_set, id_type):
    """
    Simple attribute mapper using a local mapping file
    and the myGene.info database.
    Mapped attributes:
    KEGG pathway, GO biological process, GO molecular function, GO cellular component

    :param gene_set: Set of gene ids
    :param id_type: Gene ID type of input gene set
    :return: Dataframe
    """
    # ===== Get gene ID mappings =====
    gene_mapping, _, _ = mu.get_prev_mapping(in_set=gene_set, id_type=id_type,
                                             file=config.FILES_DIR + 'gene_id_mapping.csv', sep=",")
    df, missing, prev_mapping = mu.get_prev_mapping(in_set=set(gene_mapping['entrezgene']), id_type='entrez',
                                                    file=config.FILES_DIR + 'gene_att_mapping.csv', sep=",")
    if len(missing) > 0:
        mg = get_client("gene")
        mapping = mg.querymany(missing, scopes=','.join(config.GENE_IDS),
                               fields=','.join(config.GENE_ATTRIBUTES),
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping.rename(columns={'query': 'entrezgene'}, inplace=True)
        for attribute in config.GENE_ATTRIBUTES_KEY:
            mapping = mu.preprocess_results(mapping=mapping, multicol=attribute,
                                            singlecol=attribute + '.' + config.GENE_ATTRIBUTES_KEY[attribute],
                                            key=config.GENE_ATTRIBUTES_KEY[attribute])
        mapping = mapping.drop(columns=['_id', '_score'])
        # ===== Add results from missing values =====
        pd.concat([prev_mapping, mapping]).to_csv(config.FILES_DIR + 'gene_att_mapping.csv', index=False)
        df = pd.concat([df, mapping]).reset_index(drop=True)
    # work with not unique values...
    mapping_subset = gene_mapping[['entrezgene', config.ID_TYPE_KEY[id_type]]].drop_duplicates()
    df = pd.merge(mapping_subset, df, on=['entrezgene'], how='outer')
    df = df.drop(columns=['entrezgene'])
    df = df.fillna('').groupby([config.ID_TYPE_KEY[id_type]], as_index=False).agg(
        {x: mu.combine_rows for x in config.GENE_ATTRIBUTES_KEY})
    return df
