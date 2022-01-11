#!/usr/bin/python3
import os
import pandas as pd
from d_utils import config as c, runner_utils as ru, eval_utils as eu, mapping_utils as mu
from mappers import mapping_transformer as mt, gene_mapper as gm, disease_mapper as dm
from mappers.mapper import Mapper, FileMapper


def load_files(mapper: Mapper):
    ru.print_current_usage('Starting Setup ...')

    # ---- Load disorder ids ----
    ru.print_current_usage('Load NeDrEx disorder ids ...')
    disorder_ids = pd.read_csv(c.NEDREX_DISORDER_IDS, sep="\t")
    icd10_ids = pd.read_csv(c.NEDREX_ICD10_IDS, sep="\t")

    # ---- Transform disease mapping ----
    ru.print_current_usage('Transform loaded mappings ...')
    disease_mapping = mt.transform_id_mapping(disorder_ids.merge(icd10_ids, on='primaryDomainId', how='outer'))

    # ---- Save disease mapping ----
    mapper.update_mappings(in_df=disease_mapping, key='disorder_ids')

    # ---- Load gene ids ----
    ru.print_current_usage('Load NeDrEx gene ids ...')
    gene_ids = pd.read_csv(c.NEDREX_GENE_IDS, sep="\t")

    # ---- Transform gene mapping ----
    ru.print_current_usage('Transform loaded mappings ...')
    gene_ids = gene_ids.rename(columns={'primaryDomainId': "entrez"})
    gene_ids.entrez = gene_ids.entrez.str.replace(r'entrez.', '', regex=True)

    # ---- Get attributes ----
    ru.print_current_usage('Get gene attributes ...')
    gene_att_mapping = gm.get_gene_to_attributes(gene_set=gene_ids.entrez, id_type="entrez", mapper=mapper)

    ru.print_current_usage('Get disease attributes ...')
    # ---- Get DisGeNet attributes ----
    disgenet_mapping = pd.read_csv(c.DISGENET_DIS_MAP, compression='gzip', sep='\t', dtype=str)
    disgenet_mapping = disgenet_mapping[disgenet_mapping['vocabulary'] == 'MONDO'][['diseaseId', 'code']]
    disgenet_mapping = disgenet_mapping.rename(columns={'code': 'mondo'})

    var_mapping = mu.transform_disgenet_mapping(mapping=disgenet_mapping, file=c.DISGENET_REL_VARS, col_old='snpId',
                                                col_new='disgenet.variants_related_to_disease')
    gene_mapping = mu.transform_disgenet_mapping(mapping=disgenet_mapping, file=c.DISGENET_REL_GENES, col_old='geneId',
                                                 col_new='disgenet.genes_related_to_disease')
    disease_att_mapping = pd.merge(var_mapping[['mondo', 'disgenet.variants_related_to_disease']],
                                   gene_mapping[['mondo', 'disgenet.genes_related_to_disease']],
                                   on="mondo", how="outer")

    # ---- Get KEGG attributes ----
    omim_to_hsa = pd.read_csv(c.KEGG_OMIM_TO_HSA, names=['hsa', 'omim', 'dir'], sep="\t", dtype=str)
    hsa_to_pathway = pd.read_csv(c.KEGG_HSA_TO_PATH, names=['hsa', 'pathway'], sep="\t", dtype=str)
    omim_to_pathway = pd.merge(omim_to_hsa[['hsa', 'omim']], hsa_to_pathway[['hsa', 'pathway']], on="hsa", how="inner")
    omim_to_pathway.omim = omim_to_pathway.omim.str.replace('omim:', '')
    omim_to_pathway.pathway = omim_to_pathway.pathway.str.replace('path:', '')
    omim_to_pathway = pd.merge(disease_mapping[['mondo', 'omim']], omim_to_pathway[['omim', 'pathway']], on="omim",
                               how="inner")
    omim_to_pathway = omim_to_pathway[['mondo', 'pathway']].fillna('').groupby(['mondo'], as_index=False).agg(
        mu.combine_rows)
    omim_to_pathway.rename(columns={'pathway': 'ctd.pathway_related_to_disease'}, inplace=True)
    mapping = dm.get_attributes_from_database(
        missing=['MONDO:' + x for x in set(disease_att_mapping.mondo) - set(omim_to_pathway.mondo)],
        attributes=['ctd.pathway_related_to_disease'])
    mapping = pd.concat([omim_to_pathway, mapping])
    disease_att_mapping = pd.merge(disease_att_mapping, mapping[['mondo', 'ctd.pathway_related_to_disease']],
                                   on="mondo", how="left")

    # ---- Save mapping ----
    mapper.update_mappings(in_df=disease_att_mapping, key='disorder_atts')
    disease_att_mapping = dm.get_disease_to_attributes(disease_set=disease_mapping.mondo, id_type="mondo",
                                                       mapper=mapper)
    mapper.save_mappings()
    mapper.drop_mappings()

    # ---- Calculate pairwise comparisons ----
    ru.print_current_usage('Precalculate pairwise distances ...')

    ru.print_current_usage('Precalculate pairwise distances for genes ...')
    mapper.update_distance_ids(in_series=gene_att_mapping[c.ID_TYPE_KEY['entrez']].tolist(), key='gene_mat_ids')
    for attribute in gene_att_mapping.columns[1:]:
        ru.print_current_usage('Precalculate pairwise distances for '+attribute)
        comp_mat = eu.get_distance_matrix(,
        mapper.update_distances(in_mat=comp_mat, key=c.GENE_DISTANCES[attribute])

    ru.print_current_usage('Precalculate pairwise distances for diseases ...')
    mapper.update_distance_ids(in_series=disease_att_mapping['mondo'].tolist(), key='disease_mat_ids')
    for attribute in disease_att_mapping.columns[1:]:
        ru.print_current_usage('Precalculate pairwise distances for ' + attribute)
        comp_mat = eu.get_distance_matrix(,
        mapper.update_distances(in_mat=comp_mat, key=c.DISEASE_DISTANCES[attribute])

    mapper.save_distances()

    # ---- Get PPI network ----
    # ru.print_current_usage('Get PPI network ...')

    ru.print_current_usage('Finished Setup ...')


if __name__ == "__main__":
    load_files(mapper=FileMapper())
