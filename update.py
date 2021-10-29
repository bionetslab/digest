#!/usr/bin/python3
import pandas as pd

from d_utils import config as c
from mappers import mapping_transformer as mt


def update_files():
    # ---- Load files ----
    disorder_ids = pd.read_csv(c.NEDREX_DISORDER_IDS, sep="\t")
    icd10_ids = pd.read_csv(c.NEDREX_ICD10_IDS, sep="\t")
    # ---- Tranform mappings ----
    disease_mapping = mt.transform_id_mapping(disorder_ids.merge(icd10_ids, on='primaryDomainId', how='outer'))

    # ---- Save mapping ----
    disease_mapping.to_csv(c.FILES_DIR+"disorders.map", sep="\t", index=False)
    return


if __name__ == "__main__":
    update_files()
