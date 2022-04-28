#!/usr/bin/python3

import os
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# SETUP
# ============================================================================
MAIN_GENE_ID = 'entrezgene'
SUPPORTED_GENE_IDS = ['entrez', 'ensembl', 'symbol', 'uniprot']
SUPPORTED_DISEASE_IDS = ['mondo', 'omim', 'snomedct', 'umls', 'orpha', 'mesh', 'doid', 'ICD-10']
NUMBER_OF_RANDOM_RUNS = 1000

# =============================================================================
# Set directories
# ============================================================================
FILES_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/mapping_files/'

# =============================================================================
# Set ID names mapping
# ============================================================================
ID_TYPE_KEY = {'entrez': 'entrezgene', 'ensembl': 'ensembl.gene', 'symbol': 'symbol', 'uniprot': 'uniprot.Swiss-Prot',
               'mondo': 'mondo', 'omim': 'omim', 'snomedct': 'snomedct', 'umls': 'umls', 'orpha': 'orpha',
               'mesh': 'mesh', 'doid': 'doid', 'ICD-10': 'ICD-10'}
ENRICH_KEY = {'GO_Molecular_Function_2015': 'go.MF', 'GO_Biological_Process_2015': 'go.BP',
              'GO_Cellular_Component_2015': 'go.CC', 'KEGG_2016': 'pathway.kegg'}

# =============================================================================
# Set mapping attributes
# ============================================================================
GENE_IDS = ['entrezgene', 'symbol', 'ensembl.gene', 'uniprot.Swiss-Prot']

GENE_ATTRIBUTES = ['go.BP.id', 'go.CC.id', 'go.MF.id', 'pathway.kegg.id']
GENE_ATTRIBUTES_KEY = {'go.BP': 'id', 'go.CC': 'id', 'go.MF': 'id', 'pathway.kegg': 'id'}

DISEASE_ATTRIBUTES = ['disgenet.genes_related_to_disease.gene_id', 'disgenet.variants_related_to_disease.rsid',
                      'ctd.pathway_related_to_disease.kegg_pathway_id']
DISEASE_ATTRIBUTES_KEY = {'disgenet.genes_related_to_disease': 'gene_id',
                          'disgenet.variants_related_to_disease': 'rsid',
                          'ctd.pathway_related_to_disease': 'kegg_pathway_id'}

# =============================================================================
# Set distance attributes
# ============================================================================
DISTANCES = {  # GENES
             'go.BP': 'go_BP',
             'go.CC': 'go_CC',
             'go.MF': 'go_MF',
             'pathway.kegg': 'pathway_kegg',
             # DISEASES
             'disgenet.genes_related_to_disease': 'related_genes',
             'disgenet.variants_related_to_disease': 'related_variants',
             'ctd.pathway_related_to_disease': 'related_pathways'}

# =============================================================================
# Set naming replacements
# ============================================================================
replacements = {"disgenet.genes_related_to_disease": "related_genes",
                "disgenet.variants_related_to_disease": "related_variants",
                "ctd.pathway_related_to_disease": "KEGG",
                "go.BP": "GO.BP", "go.CC": "GO.CC", "go.MF": "GO.MF", "pathway.kegg": "KEGG"}

# =============================================================================
# Set API paths to digest data
# ============================================================================
DIGEST = "https://api.digest-validation.net/files?"

# =============================================================================
# Set API paths to nedrex data
# ============================================================================
NEDREX_DISORDER_IDS = "https://api.nedrex.net/disorder/attributes/domainIds/tsv"
NEDREX_ICD10_IDS = "https://api.nedrex.net/disorder/attributes/icd10/tsv"
NEDREX_GENE_IDS = "https://api.nedrex.net/gene/attributes/primaryDomainId/tsv"

# =============================================================================
# Set API paths to DisGeNET data
# ============================================================================
DISGENET_REL_GENES = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz"
DISGENET_REL_VARS = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_variant_disease_associations.tsv.gz"
DISGENET_DIS_MAP = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/disease_mappings.tsv.gz"

# =============================================================================
# Set API paths to KEGG data
# ============================================================================
KEGG_OMIM_TO_HSA = "http://rest.genome.jp/link/omim/hsa"
KEGG_HSA_TO_PATH = "http://rest.kegg.jp/link/pathway/hsa"
