#!/usr/bin/python3
import pandas as pd

from d_utils import config as c
from mappers import mapping_transformer as mt


def update_files():
    load_nedrex_files()
    return


def load_nedrex_files():
    disorder_ids = pd.read_csv(c.NEDREX_DISORDER_IDS, sep="\t")
    icd10_ids = pd.read_csv(c.NEDREX_ICD10_IDS, sep="\t")
    mt.transform_id_mapping(disorder_ids.merge(icd10_ids, on='primaryDomainId', how='outer'))


if __name__ == "__main__":
    update_files()
