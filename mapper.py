#!/usr/bin/python3

import pandas as pd
from biothings_client import get_client


def get_gene_mapping(gene_set):
    """
    Simple gene ID mapper using a local mapping file
    and the myGene.info database.
    Supported ID types:
    uniprot, ensembl, entrezgene, genesymbol

    :param gene_set: Set of gene ids
    :return: Dataframe
    """
    gene_ids = ['uniprot.Swiss-Prot', 'symbol', 'ensembl.gene', 'entrezgene']
    mg = get_client("gene")
    mg.querymany(gene_set, scopes=','.join(gene_ids),
                 fields=','.join(gene_ids),
                 species='human', returnall=False, as_dataframe=True, df_index=False)