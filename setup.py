#!/usr/bin/python3
import pandas as pd
from d_utils import config as c, runner_utils as ru
from mappers import mapping_transformer as mt, gene_mapper as gm, disease_mapper as dm
from mappers.mapper import Mapper, FileMapper


def load_files(mapper: Mapper):
    # ---- Load disorder ids ----
    ru.print_current_usage('Starting Setup ...')
    ru.print_current_usage('Load NeDrEx disorder ids ...')
    disorder_ids = pd.read_csv(c.NEDREX_DISORDER_IDS, sep="\t")
    icd10_ids = pd.read_csv(c.NEDREX_ICD10_IDS, sep="\t")

    # ---- Tranform mappings ----
    ru.print_current_usage('Transform loaded mappings ...')
    disease_mapping = mt.transform_id_mapping(disorder_ids.merge(icd10_ids, on='primaryDomainId', how='outer'))

    # ---- Save mapping ----
    mapper.update_mappings(in_df=disease_mapping, key='disorder_ids')

    # ---- Load gene ids ----
    ru.print_current_usage('Load NeDrEx gene ids ...')
    gene_ids = pd.read_csv(c.NEDREX_GENE_IDS, sep="\t")

    # ---- Tranform list ----
    ru.print_current_usage('Transform loaded mappings ...')
    gene_ids = gene_ids.rename(columns={'primaryDomainId': "entrez"})
    gene_ids.entrez = gene_ids.entrez.str.replace(r'entrez.', '', regex=True)

    # ---- Get attributes ----
    ru.print_current_usage('Get gene attributes ...')
    gm.get_gene_to_attributes(gene_set=gene_ids.entrez, id_type="entrez", mapper=mapper)
    dm.get_disease_to_attributes(disease_set=disease_mapping.mondo, id_type="mondo", mapper=mapper)
    mapper.save_mappings()

    # ---- Get PPI network ----
    ru.print_current_usage('Precalculate pairwise distances ...')

    # ---- Get PPI network ----
    # ru.print_current_usage('Get PPI network ...')

    ru.print_current_usage('Finished Setup ...')


if __name__ == "__main__":
    load_files(mapper=FileMapper())
