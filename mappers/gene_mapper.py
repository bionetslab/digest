#!/usr/bin/python3

import pandas as pd
from d_utils import config, mapping_utils as mu
from mapper import Mapper
from biothings_client import get_client
import gseapy


def get_gene_mapping(gene_set, id_type, mapper: Mapper):
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
    hit_mapping, missing_hits = mapper.get_loaded_mapping(in_set=gene_set, id_type=config.ID_TYPE_KEY[id_type], key='gene_ids')
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
        drop_cols = ['_id', '_score', 'notfound'] if 'notfound' in mapping.columns else ['_id', '_score']
        mapping = mapping.drop(columns=drop_cols)
        mapping.dropna(subset=config.GENE_IDS[2:], inplace=True)
        # ===== Add results from missing values =====
        mapper.update_mappings(in_df=mapping, key='gene_ids')
        hit_mapping = pd.concat([hit_mapping, mapping])
    return hit_mapping


def get_gene_to_attributes(gene_set, id_type, mapper: Mapper):
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
    gene_mapping = get_gene_mapping(gene_set=gene_set, id_type=id_type, mapper=mapper)
    hit_mapping, missing_hits = mapper.get_loaded_mapping(in_set=set(gene_mapping['entrezgene']), id_type='entrezgene',
                                                          key='gene_atts')
    #print(missing_hits)
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
        mapping[mapping.columns[1:]] = mapping[mapping.columns[1:]].fillna('').applymap(mu.string_to_set)
        # ===== Add results from missing values =====
        mapper.update_mappings(in_df=mapping, key='gene_atts')
        hit_mapping = pd.concat([hit_mapping, mapping])
    # work with not unique values...
    columns = ['entrezgene', config.ID_TYPE_KEY[id_type]] if id_type != 'entrez' else ['entrezgene']
    mapping_subset = gene_mapping[columns].drop_duplicates()
    hit_mapping = pd.merge(mapping_subset, hit_mapping, on=['entrezgene'], how='outer')
    hit_mapping = hit_mapping.drop(columns=['entrezgene']) if id_type != 'entrez' else hit_mapping
    hit_mapping = hit_mapping.fillna('').groupby([config.ID_TYPE_KEY[id_type]], as_index=False).agg(
        {x: mu.combine_rowsets for x in config.GENE_ATTRIBUTES_KEY})
    return hit_mapping


def get_enriched_attributes(gene_set, id_type, mapper: Mapper):
    gene_mapping = get_gene_mapping(gene_set=gene_set, id_type=id_type, mapper=mapper)
    enrich_df = gseapy.enrichr(
        gene_list=list(gene_mapping['symbol']),
        description='atts',
        gene_sets=list(config.ENRICH_KEY.keys()),
        cutoff=0.05).results
    enrich_df = enrich_df[enrich_df['Adjusted P-value'] < 0.05]
    if len(enrich_df) > 0:
        enrich_df.insert(2, 'Term_ID', enrich_df['Term'].str.extract(r'(GO:[0-9]*|hsa[0-9]*)')[0])
    enrich_df = enrich_df[['Gene_set', 'Term_ID']].pivot(columns='Gene_set')['Term_ID']
    return enrich_df
