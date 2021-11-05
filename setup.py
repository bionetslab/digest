#!/usr/bin/python3
import pandas as pd

from d_utils import config as c
from mappers import mapping_transformer as mt


def load_files():
    # ---- Load disease ids ----
    disorder_ids = pd.read_csv(c.NEDREX_DISORDER_IDS, sep="\t")
    icd10_ids = pd.read_csv(c.NEDREX_ICD10_IDS, sep="\t")
    # ---- Tranform mappings ----
    disease_mapping = mt.transform_id_mapping(disorder_ids.merge(icd10_ids, on='primaryDomainId', how='outer'))
    # ---- Save mapping ----
    disease_mapping.to_csv(c.FILES_DIR+"disorders.map", sep="\t", index=False)

    # ---- Load gene ids ----
    gene_ids = pd.read_csv(c.NEDREX_GENE_IDS, sep="\t")
    # ---- Tranform list ----
    gene_ids = gene_ids.rename(columns={'primaryDomainId': "entrez"})
    gene_ids.entrez = gene_ids.entrez.str.replace(r'entrez.', '', regex=True)
    # ---- Save mapping ----
    gene_ids.to_csv(c.FILES_DIR + "gene.list", sep="\t", index=False)
    return


if __name__ == "__main__":
    load_files()
