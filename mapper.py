#!/usr/bin/python3

import pandas as pd
import numpy as np
from biothings_client import get_client

ID_TYPE_KEY = {'entrez': 'entrezgene', 'ensembl': 'ensembl.gene', 'symbol': 'symbol', 'uniprot': 'uniprot.Swiss-Prot'}
GENE_IDS = ['uniprot.Swiss-Prot', 'symbol', 'ensembl.gene', 'entrezgene']


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
    df, missing, prev_mapping = _get_prev_mapping(gene_set=gene_set, id_type=id_type, file='gene_id_mapping.csv')
    # ===== Get mapping for missing values =====
    if len(missing) > 0:
        mg = get_client("gene")
        mapping = mg.querymany(missing, scopes=ID_TYPE_KEY[id_type], fields=','.join(GENE_IDS),
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping = mapping.drop(columns=[ID_TYPE_KEY[id_type]])
        mapping.rename(columns={'query': ID_TYPE_KEY[id_type]}, inplace=True)
        # ===== Split if there are multiple ensembl ids =====
        if 'ensembl' in mapping:
            mapping = _preprocess_results(mapping=mapping, multicol='ensembl', singlecol='ensembl.gene', key='gene',
                                          explode=True)
        mapping = mapping.drop(columns=['_id', '_score'])
        # ===== Add results from missing values =====
        pd.concat([prev_mapping, mapping]).to_csv('gene_id_mapping.csv', index=False)
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
    gene_mapping, _, _ = _get_prev_mapping(gene_set=gene_set, id_type=id_type, file='gene_id_mapping.csv')
    # ===== Get mapping from previous mappings =====
    df, missing, prev_mapping = _get_prev_mapping(gene_set=set(gene_mapping['entrezgene']),
                                                  id_type='entrez', file='gene_att_mapping.csv')
    # ===== Get mapping for missing values =====
    if len(missing) > 0:
        mg = get_client("gene")
        gene_ids = ['uniprot.Swiss-Prot', 'symbol', 'ensembl.gene', 'entrezgene']
        mapping = mg.querymany(missing, scopes=','.join(gene_ids),
                               fields='pathway.kegg.id, go.BP.id, go.CC.id, go.MF.id',
                               species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping.rename(columns={'query': 'entrezgene'}, inplace=True)
        for column in ['go.BP', 'go.CC', 'go.MF', 'pathway.kegg']:
            mapping = _preprocess_results(mapping=mapping, multicol=column, singlecol=column+'.id', key='id')
        mapping = mapping.drop(columns=['_id', '_score'])
        # ===== Add results from missing values =====
        pd.concat([prev_mapping, mapping]).to_csv('gene_att_mapping.csv', index=False)
        df = pd.concat([df, mapping]).reset_index(drop=True)
    # ===== Reverse map to input gene ids =====
    mapping_subset = gene_mapping[['entrezgene', ID_TYPE_KEY[id_type]]].drop_duplicates()
    df = pd.merge(mapping_subset, df, on=['entrezgene'], how='outer')
    df = df.drop(columns=['entrezgene'])
    # ===== Combine values if original gene id appears multiple times due to mapping =====
    df = df.fillna('').groupby(['uniprot.Swiss-Prot'], as_index=False).agg({'go.BP': ';'.join, 'go.CC': ';'.join,
                                                                            'go.MF': ';'.join, 'pathway.kegg': ';'.join})
    return df


def _get_prev_mapping(gene_set, id_type, file):
    """
    Get previous mappings from local file and filter which genes are
    still missing and also returned.
    :param gene_set: Set of genes to map
    :param id_type: Gene ID type of the genes
    :param file: Filename of the pre-mapped values
    :returns:
        - df: Dataframe of previously mapped genes of the list
        - missing: List of not previously mapped genes of the set
        - prev_mapping: Full dataframe of all previously mapped genes
    """
    # ===== Get mapping from local mapping file =====
    prev_mapping = pd.read_csv(file, header=0, dtype=str)
    df = prev_mapping[prev_mapping[ID_TYPE_KEY[id_type]].isin(gene_set)]
    # ===== Get missing values =====
    missing = list(set(gene_set)-set(prev_mapping[ID_TYPE_KEY[id_type]]))
    return df, missing, prev_mapping


def _preprocess_results(mapping, multicol, singlecol, key, explode=False):
    """
    Mapping may result in two columns depending if there are multiple values for
    one id or only one value. This method combines the resulting two columns,
    if they exist, into one.
    :param mapping: Dataframe with mapping
    :param multicol: Name of the column with multiple values
    :param singlecol: Name of column with single values
    :param key: Name of key in dict of multiple values
    :param explode: If True, each value from multiple values is
                    separated into separate rows. [Default=False]
    :return: Dataframe
    """
    def _convert_to_string(cell):
        if str(cell) != 'nan':
            extracted_ids = [val.get(key) for val in cell]
            return ';'.join(extracted_ids)
        return cell

    mapping[multicol] = mapping[multicol].apply(lambda x: _convert_to_string(x)) if multicol in mapping else np.nan
    if singlecol in mapping:
        mapping[multicol].fillna(mapping[singlecol], inplace=True)
        mapping = mapping.drop(columns=[singlecol])
    if explode:
        mapping = mapping[multicol].split(';').explode(multicol)
        mapping.rename(columns={multicol: singlecol}, inplace=True)
    return mapping
