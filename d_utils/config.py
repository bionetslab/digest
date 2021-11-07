#!/usr/bin/python3

import sys
import os
import inspect

# =============================================================================
# SETUP
# ============================================================================
MAIN_GENE_ID = 'entrezgene'
SUPPORTED_GENE_IDS = ['entrez', 'ensembl', 'symbol', 'uniprot']
SUPPORTED_DISEASE_IDS = ['mondo', 'omim', 'snomedct', 'umls', 'orpha', 'mesh', 'doid', 'ICD-10']
NUMBER_OF_RANDOM_RUNS = 1000

# =============================================================================
# Load all paths of required folders with scripts
# =============================================================================
UTILS_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
BASE_DIR = os.path.dirname(UTILS_DIR)
sys.path.insert(0, BASE_DIR + '/setup')
sys.path.insert(0, BASE_DIR + '/d_utils')
sys.path.insert(0, BASE_DIR + '/mappers')

# =============================================================================
# Set directories
# ============================================================================
FILES_DIR = BASE_DIR + '/mapping_files/'

# =============================================================================
# Set ID names mapping
# ============================================================================
ID_TYPE_KEY = {'entrez': 'entrezgene', 'ensembl': 'ensembl.gene', 'symbol': 'symbol', 'uniprot': 'uniprot.Swiss-Prot'}

# =============================================================================
# Set mapping attributes
# ============================================================================
GENE_IDS = ['uniprot.Swiss-Prot', 'symbol', 'ensembl.gene', 'entrezgene']

GENE_ATTRIBUTES = ['go.BP.id', 'go.CC.id', 'go.MF.id', 'pathway.kegg.id']
GENE_ATTRIBUTES_KEY = {'go.BP': 'id', 'go.CC': 'id', 'go.MF': 'id', 'pathway.kegg': 'id'}

DISEASE_ATTRIBUTES = ['disgenet.genes_related_to_disease.gene_id', 'disgenet.variants_related_to_disease.rsid',
                      'ctd.pathway_related_to_disease.kegg_pathway_id']
DISEASE_ATTRIBUTES_KEY = {'disgenet.genes_related_to_disease': 'gene_id',
                          'disgenet.variants_related_to_disease': 'rsid',
                          'ctd.pathway_related_to_disease': 'kegg_pathway_id'}

# =============================================================================
# Set API paths to nedrex data
# ============================================================================
NEDREX_DISORDER_IDS = "https://api.nedrex.net/disorder/attributes/domainIds/tsv"
NEDREX_ICD10_IDS = "https://api.nedrex.net/disorder/attributes/icd10/tsv"
NEDREX_GENE_IDS = "https://api.nedrex.net/gene/attributes/primaryDomainId/tsv"


