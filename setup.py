#!/usr/bin/python3
import pandas as pd
from d_utils import config as c, runner_utils as ru
from mappers import mapping_transformer as mt, gene_mapper as gm, disease_mapper as dm
from mappers.mapper import FileMapper


def load_files():
    # ---- Load disorder ids ----
    ru.print_current_usage('Starting Setup ...')
    mapper = FileMapper()
    ru.print_current_usage('Load NeDrEx disorder ids ...')
    disorder_ids = pd.read_csv(c.NEDREX_DISORDER_IDS, sep="\t")
    icd10_ids = pd.read_csv(c.NEDREX_ICD10_IDS, sep="\t")
    # ---- Tranform mappings ----
    ru.print_current_usage('Transform loaded mappings ...')
    disease_mapping = mt.transform_id_mapping(disorder_ids.merge(icd10_ids, on='primaryDomainId', how='outer'))
    # ---- Save mapping ----
    mapper.update_mappings(in_df=disease_mapping, key='disorder_ids')
    disease_mapping.to_csv(c.FILES_DIR + "disorders.map", sep="\t", index=False)
    # ---- Load gene ids ----
    ru.print_current_usage('Load NeDrEx gene ids ...')
    gene_ids = pd.read_csv(c.NEDREX_GENE_IDS, sep="\t")
    # ---- Tranform list ----
    ru.print_current_usage('Transform loaded mappings ...')
    gene_ids = gene_ids.rename(columns={'primaryDomainId': "entrez"})
    gene_ids.entrez = gene_ids.entrez.str.replace(r'entrez.', '', regex=True)
    # ---- Save mapping ----
    gene_ids.to_csv(c.FILES_DIR + "gene.list", sep="\t", index=False)

    # ---- Get attributes ----
    ru.print_current_usage('Get gene attributes ...')
    gm.get_gene_to_attributes(gene_set=gene_ids.entrez, id_type="entrez", mapper=mapper)
    dm.get_disease_to_attributes(disease_set=disease_mapping.mondo, id_type="mondo", mapper=mapper)
    mapper.save_mappings()
    ru.print_current_usage('Finished Setup ...')


if __name__ == "__main__":
    load_files()
