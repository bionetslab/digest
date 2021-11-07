#!/usr/bin/python3

import pandas as pd
from d_utils import config, mapping_utils as mu
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
    hit_mapping, missing_hits, full_mapping = mu.get_prev_mapping(in_set=gene_set, id_type=id_type,
                                                                  file=config.FILES_DIR + 'gene_id_mapping.csv',
                                                                  sep=",")
    # ===== Get mapping for missing values =====
    if len(missing_hits) > 0:
        mg = get_client("gene")
        mapping = mg.querymany(missing_hits, scopes=config.ID_TYPE_KEY[id_type], fields=','.join(config.GENE_IDS),
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping = mapping.drop(columns=[config.ID_TYPE_KEY[id_type]])
        mapping.rename(columns={'query': config.ID_TYPE_KEY[id_type]}, inplace=True)
        # ===== Split if there are multiple ensembl ids =====
        if 'ensembl' in mapping:
            mapping = mu.preprocess_results(mapping=mapping, multicol='ensembl', singlecol='ensembl.gene', key='gene',
                                            explode=True)
        mapping = mapping.drop(columns=['_id', '_score'])
        # ===== Add results from missing values =====
        hit_mapping = mu.save_mapping_in_file(full_map=full_mapping, prev_map=hit_mapping, new_map=mapping,
                                       file=config.FILES_DIR + 'gene_id_mapping.csv')
    return hit_mapping


def get_gene_to_attributes(gene_set, id_type):
    """
    Simple attribute mapper using a local mapping file
    and the myGene.info database.
    Mapped attributes:
    KEGG pathway, GO biological process, GO molecular function, GO cellular component

    :param gene_set: set of gene ids
    :param id_type: id type of set
    :return: Dataframe
    """
    # ===== Get gene ID mappings =====
    gene_mapping = get_gene_mapping(gene_set=gene_set, id_type=id_type)
    hit_mapping, missing_hits, full_mapping = mu.get_prev_mapping(in_set=set(gene_mapping['entrezgene']),
                                                                  id_type='entrez',
                                                                  file=config.FILES_DIR + 'gene_att_mapping.csv',
                                                                  sep=",")
    if len(missing_hits) > 0:
        mg = get_client("gene")
        mapping = mg.querymany(missing_hits, scopes=','.join(config.GENE_IDS),
                               fields=','.join(config.GENE_ATTRIBUTES),
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping.rename(columns={'query': 'entrezgene'}, inplace=True)
        for attribute in config.GENE_ATTRIBUTES_KEY:
            mapping = mu.preprocess_results(mapping=mapping, multicol=attribute,
                                            singlecol=attribute + '.' + config.GENE_ATTRIBUTES_KEY[attribute],
                                            key=config.GENE_ATTRIBUTES_KEY[attribute])
        mapping = mapping.drop(columns=['_id', '_score'])
        # ===== Add results from missing values =====
        print(mapping)
        hit_mapping = mu.save_mapping_in_file(full_map=full_mapping, prev_map=hit_mapping, new_map=mapping,
                                              file=config.FILES_DIR + 'gene_att_mapping.csv')
    # work with not unique values...
    columns = ['entrezgene', config.ID_TYPE_KEY[id_type]] if id_type != 'entrez' else ['entrezgene']
    mapping_subset = gene_mapping[columns].drop_duplicates()
    hit_mapping = pd.merge(mapping_subset, hit_mapping, on=['entrezgene'], how='outer')
    hit_mapping = hit_mapping.drop(columns=['entrezgene']) if id_type != 'entrez' else hit_mapping
    hit_mapping = hit_mapping.fillna('').groupby([config.ID_TYPE_KEY[id_type]], as_index=False).agg(
        {x: mu.combine_rows for x in config.GENE_ATTRIBUTES_KEY})
    return hit_mapping
