#!/usr/bin/python3

import pandas as pd
from .. import config
from . import mapping_utils as mu
from .mapper import Mapper
from biothings_client import get_client
import gseapy


def get_gene_mapping(gene_set: set, id_type: str, mapper: Mapper):
    """
    Simple gene ID mapper using local mapping files generated during setup
    and the myGene.info database.
    Supported ID types:
    uniprot, ensembl, entrezgene, genesymbol

    :param gene_set: set of gene ids
    :param id_type: gene ID type of input gene set
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :return: gene to id mapping as dataframe
    """
    # ===== Get mapping from previous mappings =====
    hit_mapping, missing_hits = mapper.get_loaded_mapping(in_set=gene_set, id_type=config.ID_TYPE_KEY[id_type], key='gene_ids')
    # ===== Get mapping for missing values =====
    if len(missing_hits) > 0:
        mg = get_client("gene")
        mapping = mg.querymany(missing_hits, scopes=config.ID_TYPE_KEY[id_type], fields=','.join(config.GENE_IDS),
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        if 'notfound' in mapping:
            mapping = mapping[mapping['notfound']!=True]
        # ===== Only when missing values found a mapping =====
        if not mapping.empty:
            mapping = mapping.drop(columns=[config.ID_TYPE_KEY[id_type]])
            mapping.rename(columns={'query': config.ID_TYPE_KEY[id_type]}, inplace=True)
            # ===== Fill if ID types are missing =====
            for id_col in config.GENE_IDS:
                if id_col not in mapping:
                    mapping[id_col] = ""
            # ===== Split if there are multiple ensembl ids =====
            if 'ensembl' in mapping:
                mapping = mu.preprocess_results(mapping=mapping, multicol='ensembl', singlecol='ensembl.gene', key='gene')
                mapping.rename(columns={'ensembl': 'ensembl.gene'}, inplace=True)
            mapping['uniprot.Swiss-Prot'] = mapping['uniprot.Swiss-Prot'].fillna("").apply(mu.list_to_string)
            drop_cols = ['_id', '_score', 'notfound'] if 'notfound' in mapping.columns else ['_id', '_score']
            mapping = mapping.drop(columns=drop_cols)
            mapping = mapping.fillna('').groupby(['entrezgene','symbol'], as_index=False).agg(
                {x: mu.combine_rows_to_string for x in config.GENE_IDS[2:]})
            # ===== Add results from missing values =====
            mapper.update_mappings(in_df=mapping, key='gene_ids')
            hit_mapping = pd.concat([hit_mapping, mapping])
    return hit_mapping


def get_gene_to_attributes(gene_set: set, id_type: str, mapper: Mapper):
    """
    Simple attribute mapper using local mapping files generated during setup
    and the myGene.info database.
    Mapped attributes:
    KEGG pathway, GO biological process, GO molecular function, GO cellular component

    :param gene_set: set of gene ids
    :param id_type: id type of set
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :return: gene to attribute mapping as dataframe
    """
    # ===== Get gene ID mappings =====
    gene_mapping = get_gene_mapping(gene_set=gene_set, id_type=id_type, mapper=mapper)
    # ===== Return empty Dataframe if IDs were not mappable =====
    if gene_mapping.empty:
        return gene_mapping
    hit_mapping, missing_hits = mapper.get_loaded_mapping(in_set=set(gene_mapping['entrezgene']), id_type='entrezgene',
                                                          key='gene_atts')
    # ===== Get mapping for missing values =====
    if len(missing_hits) > 0:
        mg = get_client("gene")
        mapping = mg.querymany(missing_hits, scopes=','.join(config.GENE_IDS),
                               fields=','.join(config.GENE_ATTRIBUTES),
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping.rename(columns={'query': 'entrezgene'}, inplace=True)
        if 'notfound' in mapping:
            mapping = mapping[mapping['notfound']!=True]
        if not mapping.empty:
            for attribute in config.GENE_ATTRIBUTES_KEY:
                mapping = mu.preprocess_results(mapping=mapping, multicol=attribute,
                                                singlecol=attribute + '.' + config.GENE_ATTRIBUTES_KEY[attribute],
                                                key=config.GENE_ATTRIBUTES_KEY[attribute])
            drop_cols = ['_id', '_score', 'notfound'] if 'notfound' in mapping.columns else ['_id', '_score']
            mapping = mapping.drop(columns=drop_cols)
            mapping[mapping.columns[1:]] = mapping[mapping.columns[1:]].fillna('').applymap(mu.combine_rows_to_set)
            # ===== Add results from missing values =====
            mapper.update_mappings(in_df=mapping, key='gene_atts')
            hit_mapping = pd.concat([hit_mapping, mapping])
    # ===== work with not unique values =====
    hit_mapping = mu.map_to_prev_id(main_id_type="entrezgene", id_type=config.ID_TYPE_KEY[id_type],
                                    id_mapping=gene_mapping, att_mapping=hit_mapping)
    return hit_mapping


def get_enriched_attributes(gene_set: set, id_type: str, mapper: Mapper):
    """
    Get enriched attribute values for given gene set.

    :param gene_set: gene set to get enriched values for
    :param id_type: id type of gene ids in set
    :param mapper: mapper from type Mapper defining where the precalculated information comes from
    :return: enriched attributes as dataframe
    """
    gene_mapping = get_gene_mapping(gene_set=gene_set, id_type=id_type, mapper=mapper)
    # ===== Return empty Dataframe if IDs were not mappable =====
    if gene_mapping.empty:
        return gene_mapping
    enrich_df = gseapy.enrichr(
        gene_list=list(mu.combine_rows_to_set(gene_mapping['symbol'])),
        description='atts',
        gene_sets=list(config.ENRICH_KEY.keys()),
        cutoff=0.05).results
    enrich_df = enrich_df[enrich_df['Adjusted P-value'] < 0.05]
    if enrich_df.empty:
        return pd.DataFrame()
    enrich_df.insert(2, 'Term_ID', enrich_df['Term'].str.extract(r'(GO:[0-9]*|hsa[0-9]*)')[0])
    enrich_df = enrich_df[['Gene_set', 'Term_ID']].pivot(columns='Gene_set')['Term_ID']
    return enrich_df
