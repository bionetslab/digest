#!/usr/bin/python3

import pandas as pd
import numpy as np
from biothings_client import get_client

ID_TYPE_KEY = {'entrez': 'entrezgene', 'ensembl': 'ensembl.gene', 'symbol': 'symbol', 'uniprot': 'uniprot.Swiss-Prot'}
GENE_IDS = ['uniprot.Swiss-Prot', 'symbol', 'ensembl.gene', 'entrezgene']
IN_DIR = 'Input/'


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
    df, missing, prev_mapping = _get_prev_mapping(in_set=gene_set, id_type=id_type, file='mapping_files/gene_id_mapping.csv', sep=",")
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
    gene_mapping, _, _ = _get_prev_mapping(in_set=gene_set, id_type=id_type, file='mapping_files/gene_id_mapping.csv', sep=",")
    # ===== Get mapping from previous mappings =====
    df, missing, prev_mapping = _get_prev_mapping(in_set=set(gene_mapping['entrezgene']),
                                                  id_type='entrez', file='mapping_files/gene_att_mapping.csv', sep=",")
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
    # ===== Reverse map to input gene ids if not entrez id type =====
    if ID_TYPE_KEY[id_type] != 'entrezgene':
        mapping_subset = gene_mapping[['entrezgene', ID_TYPE_KEY[id_type]]].drop_duplicates()
        df = pd.merge(mapping_subset, df, on=['entrezgene'], how='outer')
        df = df.drop(columns=['entrezgene'])
    # ===== Combine values if original gene id appears multiple times due to mapping =====
    df = df.fillna('').groupby([ID_TYPE_KEY[id_type]], as_index=False).agg({'go.BP': _combine_rows, 'go.CC': _combine_rows,
                                                                            'go.MF': _combine_rows, 'pathway.kegg': _combine_rows})
    return df


def get_disease_mapping(disease_set, id_type):
    # ==== Get Mondo IDs ====
    disease_id_set, _, _ = _get_prev_mapping(in_set=disease_set, id_type=id_type, file="mapping_files/disorders.map", sep="\t")
    mondo_set = list(set('MONDO:'+disease_id_set['mondo']))
    # ===== Get mapping from previous mappings =====
    df, missing, prev_mapping = _get_prev_mapping(in_set=mondo_set, id_type='mondo', file='mapping_files/disease_disgenet_mapping.csv', sep=",")
    # ==== Get disgenet values ====
    if len(missing) > 0:
        md = get_client("disease")
        mapping = md.getdiseases(missing,
                                 fields='disgenet.genes_related_to_disease.gene_id,disgenet.variants_related_to_disease.rsid',
                                 species='human', returnall=False, as_dataframe=True, df_index=False)
        mapping.rename(columns={'query': 'mondo'}, inplace=True)
        # transform dataframe to combine single and multiple results
        mapping = _preprocess_results(mapping=mapping, multicol='disgenet.genes_related_to_disease',
                                      singlecol='disgenet.genes_related_to_disease.gene_id', key='gene_id')
        mapping = _preprocess_results(mapping=mapping, multicol='disgenet.variants_related_to_disease',
                                      singlecol='disgenet.variants_related_to_disease.rsid', key='rsid')
        mapping = mapping.drop(columns=['_id', '_version', 'disgenet._license'])
        # ===== Add results from missing values =====
        pd.concat([prev_mapping, mapping]).to_csv('disease_disgenet_mapping.csv', index=False)
        df = pd.concat([df, mapping]).reset_index(drop=True)
    # ==== Map back to previous ids ====
    df["mondo"] = df["mondo"].str.replace("MONDO:", "")
    # work with not unique values...
    mapping_subset = disease_id_set[['mondo', id_type]].drop_duplicates()
    df = pd.merge(mapping_subset, df, on=['mondo'], how='outer')
    df = df.drop(columns=['mondo'])
    df = df.fillna('').groupby(id_type, as_index=False).agg({'disgenet.genes_related_to_disease': _combine_rows, 'disgenet.variants_related_to_disease': _combine_rows})
    return df


def _split_and_expand_column(data, split_string, column_name):
    s = data[column_name].str.split(split_string, expand=True).stack()
    i = s.index.get_level_values(0)
    df2 = data.loc[i].copy()
    df2[column_name] = s.values
    return df2


def _get_prev_mapping(in_set, id_type, file, sep):
    """
    Get previous mappings from local file and filter which elements of the set
    are still missing and also returned.
    :param in_set: Set of genes to map
    :param id_type: Gene ID type of the genes
    :param file: Filename of the pre-mapped values
    :returns:
        - mapped_set: Dataframe of previously mapped elements of the list
        - missing: List of not previously mapped elements of the set
        - prev_mapping: Full dataframe of all previously mapped genes
    """
    # ===== Get mapping from local mapping file =====
    mapping = pd.read_csv(file, sep=sep, header=0, dtype=str)
    if id_type == "ICD-10":
        mapping = _split_and_expand_column(data=mapping, split_string=",", column_name="ICD-10")
        mapping_copy = mapping.copy()
        mapping_copy['ICD-10'] = mapping_copy['ICD-10'].str.split('.', expand=True)[0]
        mapping = pd.concat([mapping, mapping_copy], ignore_index=True)
    # ==== Map given disease set ====
    set_id_type = ID_TYPE_KEY[id_type] if id_type in ID_TYPE_KEY else id_type
    mapped_set = mapping[mapping[set_id_type].isin(in_set)]
    # ===== Get missing values =====
    missing = list(set(in_set) - set(mapping[set_id_type]))
    return mapped_set, missing, mapping


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


def _combine_rows(x):
    return set(filter(None,';'.join(x).split(';')))