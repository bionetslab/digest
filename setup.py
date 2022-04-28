#!/usr/bin/python3

import os
import pandas as pd
import requests
import scipy.sparse as sp
from evaluation.d_utils import runner_utils as ru, eval_utils as eu
from evaluation import config as c
from evaluation.mappers import mapping_transformer as mt, gene_getter as gm, mapping_utils as mu
from evaluation.mappers import disease_getter as dm
from evaluation.mappers.mapper import Mapper, FileMapper
# only in full biodigest
import graph_tool as gt
import graph_tool.topology as gtt


def load_files(mapper: Mapper):
    """
    Run setup to load all needed files from api.

    :param mapper: object of type Mapper defining where and how to save the generated data
    """
    ru.print_current_usage('Starting Setup ...')


    ru.print_current_usage('Get id and attribute mappings ...')
    for file_id in mapper.loaded_mappings:
        r = requests.get(c.DIGEST+"name="+mapper.file_names[file_id])
        with open(os.path.join(mapper.files_dir, mapper.file_names[file_id]), 'wb') as f:
            f.write(r.content)

    ru.print_current_usage('Get distance mappings ...')
    for distance_measure in ["jaccard", "overlap"]:
        ru.print_current_usage('Get distance mappings for '+distance_measure+' ...')
        os.system("mkdir -p " + os.path.join(mapper.files_dir, distance_measure, ""))
        for file_id in mapper.loaded_distance_ids[distance_measure]:
            r = requests.get(c.DIGEST + "name=" + mapper.file_names[file_id] + "&measure="+distance_measure)
            with open(os.path.join(mapper.files_dir, distance_measure, mapper.file_names[file_id]), 'wb') as f:
                f.write(r.content)

        for file_id in mapper.loaded_distances[distance_measure]:
            r = requests.get(c.DIGEST + "name=" + mapper.file_names[file_id] + "&measure=" + distance_measure)
            with open(os.path.join(mapper.files_dir, distance_measure, mapper.file_names[file_id]), 'wb') as f:
                f.write(r.content)

    ru.print_current_usage('Finished Setup ...')


def create_files(mapper: Mapper):
    """
    Run setup to generate all files from scratch and gather data from databases.

    :param mapper: object of type Mapper defining where and how to save the generated data
    """
    ru.print_current_usage('Starting Setup ...')

    # ===== Start graph builder =====
    ggi = requests.post(url=c.NEDREX_GRAPH_BUILDER, json=c.NEDREX_GGI_POST)
    ddi = requests.post(url=c.NEDREX_GRAPH_BUILDER, json=c.NEDREX_DDI_POST)

    # ===== Load disorder ids =====
    ru.print_current_usage('Load NeDrEx disorder ids ...')
    disorder_ids = pd.read_csv(c.NEDREX_DISORDER_IDS, sep="\t")
    icd10_ids = pd.read_csv(c.NEDREX_ICD10_IDS, sep="\t")

    # ===== Transform disease mapping =====
    ru.print_current_usage('Transform loaded mappings ...')
    disease_mapping = mt.transform_id_mapping(disorder_ids.merge(icd10_ids, on='primaryDomainId', how='outer'))
    disease_mapping["ICD-10"] = disease_mapping["ICD-10"].apply(mt.reduce_to_parent)

    # ===== Save disease mapping =====
    mapper.update_mappings(in_df=disease_mapping, key='disorder_ids')

    # ===== Load gene ids =====
    ru.print_current_usage('Load NeDrEx gene ids ...')
    gene_ids = pd.read_csv(c.NEDREX_GENE_IDS, sep="\t")

    # ===== Transform gene mapping =====
    ru.print_current_usage('Transform loaded mappings ...')
    gene_ids = gene_ids.rename(columns={'primaryDomainId': "entrez"})
    gene_ids.entrez = gene_ids.entrez.str.replace(r'entrez.', '', regex=True)

    # ===== Get attributes =====
    ru.print_current_usage('Get gene attributes ...')
    gene_att_mapping = gm.get_gene_to_attributes(gene_set=gene_ids.entrez, id_type="entrez", mapper=mapper)

    ru.print_current_usage('Get disease attributes ...')
    # ===== Get DisGeNet attributes =====
    disgenet_mapping = pd.read_csv(c.DISGENET_DIS_MAP, compression='gzip', sep='\t', dtype=str)
    disgenet_mapping = disgenet_mapping[disgenet_mapping['vocabulary'] == 'MONDO'][['diseaseId', 'code']]
    disgenet_mapping = disgenet_mapping.rename(columns={'code': 'mondo'})

    var_mapping = mu.transform_disgenet_mapping(mapping=disgenet_mapping, file=c.DISGENET_REL_VARS, col_old='snpId',
                                                col_new='disgenet.variants_related_to_disease')
    gene_mapping = mu.transform_disgenet_mapping(mapping=disgenet_mapping, file=c.DISGENET_REL_GENES,
                                                 col_old='geneId',
                                                 col_new='disgenet.genes_related_to_disease')
    disease_att_mapping = pd.merge(var_mapping[['mondo', 'disgenet.variants_related_to_disease']],
                                   gene_mapping[['mondo', 'disgenet.genes_related_to_disease']],
                                   on="mondo", how="outer")

    # ===== Get KEGG attributes =====
    omim_to_hsa = pd.read_csv(c.KEGG_OMIM_TO_HSA, names=['hsa', 'omim', 'dir'], sep="\t", dtype=str)
    hsa_to_pathway = pd.read_csv(c.KEGG_HSA_TO_PATH, names=['hsa', 'pathway'], sep="\t", dtype=str)
    omim_to_pathway = pd.merge(omim_to_hsa[['hsa', 'omim']], hsa_to_pathway[['hsa', 'pathway']], on="hsa",
                               how="inner")
    omim_to_pathway.omim = omim_to_pathway.omim.str.replace('omim:', '')
    omim_to_pathway.pathway = omim_to_pathway.pathway.str.replace('path:', '')
    omim_to_pathway = pd.merge(disease_mapping[['mondo', 'omim']], omim_to_pathway[['omim', 'pathway']], on="omim",
                               how="inner")
    omim_to_pathway = omim_to_pathway[['mondo', 'pathway']].fillna('').groupby(['mondo'], as_index=False).agg(
        mu.combine_rows_to_set)
    omim_to_pathway.rename(columns={'pathway': 'ctd.pathway_related_to_disease'}, inplace=True)
    mapping = dm.get_attributes_from_database(
        missing=['MONDO:' + x for x in set(disease_att_mapping.mondo) - set(omim_to_pathway.mondo)],
        attributes=['ctd.pathway_related_to_disease'])
    mapping = pd.concat([omim_to_pathway, mapping])
    disease_att_mapping = pd.merge(disease_att_mapping, mapping[['mondo', 'ctd.pathway_related_to_disease']],
                                   on="mondo", how="left")

    # ===== Save mapping =====
    mapper.update_mappings(in_df=disease_att_mapping, key='disorder_atts')
    disease_att_mapping = dm.get_disease_to_attributes(disease_set=disease_mapping.mondo, id_type="mondo",
                                                       mapper=mapper)

    mapper.save_mappings()
    mapper.drop_mappings()

    # ===== Calculate pairwise comparisons =====
    ru.print_current_usage('Precalculate pairwise distances ...')

    for distance_measure in ["jaccard", "overlap"]:
        os.system("mkdir -p " + os.path.join(mapper.files_dir, distance_measure, ""))
        ru.print_current_usage('Precalculate pairwise distances for genes [' + distance_measure + '] ...')
        mapper.update_distance_ids(in_series=gene_att_mapping[c.ID_TYPE_KEY['entrez']], key='gene_mat_ids',
                                   distance_measure=distance_measure)
        for attribute in gene_att_mapping.columns[1:]:
            ru.print_current_usage('Precalculate pairwise distances for ' + attribute)
            subset_df = gene_att_mapping[gene_att_mapping[attribute].str.len() > 0]
            comp_mat = eu.get_distance_matrix(full_att_series=gene_att_mapping[attribute],
                                              from_ids=subset_df[c.ID_TYPE_KEY['entrez']],
                                              coefficient=distance_measure,
                                              id_to_index=mapper.loaded_distance_ids[distance_measure][
                                                  'gene_mat_ids'])
            sp.save_npz(os.path.join(mapper.files_dir, distance_measure, mapper.file_names[c.DISTANCES[attribute]]),
                        comp_mat)

        ru.print_current_usage('Precalculate pairwise distances for diseases [' + distance_measure + '] ...')
        mapper.update_distance_ids(in_series=disease_att_mapping['mondo'], key='disease_mat_ids',
                                   distance_measure=distance_measure)
        for attribute in disease_att_mapping.columns[1:]:
            ru.print_current_usage('Precalculate pairwise distances for ' + attribute)
            subset_df = disease_att_mapping[disease_att_mapping[attribute].str.len() > 0]
            comp_mat = eu.get_distance_matrix(full_att_series=disease_att_mapping[attribute],
                                              from_ids=subset_df['mondo'],
                                              coefficient=distance_measure,
                                              id_to_index=mapper.loaded_distance_ids[distance_measure][
                                                  'disease_mat_ids'])
            sp.save_npz(os.path.join(mapper.files_dir, distance_measure, mapper.file_names[c.DISTANCES[attribute]]),
                        comp_mat)

    mapper.save_distances()

    # ===== Calculate pairwise comparisons =====
    ru.print_current_usage('Create default networks ...')

    def recreate_network(graph_id, id_type, file_name, edges_num):
        response = requests.get(c.NEDREX_GRAPH_DOWNLOADER + graph_id.json() + ".graphml")
        open(os.path.join(mapper.files_dir, graph_id.json() + ".graphml"), "wb").write(response.content)
        graph = gt.load_graph(os.path.join(mapper.files_dir, graph_id.json() + ".graphml"))

        vertices = list()
        for v in graph.iter_vertices():
            if graph.vp['primaryDomainId'][v].startswith(id_type):
                vertices.append(graph.vertex(v))

        edges = list()
        for v in vertices:
            res = gtt.shortest_distance(graph, source=graph.vertex(v), max_dist=edges_num, target=vertices, directed=False)
            for index in range(len(res)):
                if res[index] == edges_num:
                    edges.append([graph.vp['primaryDomainId'][v],
                                  graph.vp['primaryDomainId'][vertices[index]]])

        edges = pd.DataFrame(edges, columns=["from", "to"])
        edges["from"] = edges["from"].str[len(id_type)+1:]
        edges["to"] = edges["to"].str[len(id_type)+1:]

        G = gt.Graph(directed=False)
        v_ids = G.add_edge_list(edges[["from", "to"]].values, hashed=True)
        G.vertex_properties['id'] = v_ids
        gt.stats.remove_parallel_edges(G)
        v_ids_set = set(v_ids)
        for v in vertices:
            cur_v = graph.vp['primaryDomainId'][v][len(id_type)+1:]
            if cur_v not in v_ids_set:
                new_v = G.add_vertex()
                G.vertex_properties['id'][new_v] = cur_v

        G.save(os.path.join(mapper.files_dir, file_name), fmt="graphml")
        os.system("rm " + os.path.join(mapper.files_dir, graph_id.json() + ".graphml"))

    recreate_network(graph_id=ggi, id_type="entrez", file_name="ggi_graph.graphml", edges_num=3)
    recreate_network(graph_id=ddi, id_type="mondo", file_name="ddi_graph.graphml", edges_num=2)

    ru.print_current_usage('Finished Setup ...')


def main(setup_type: str, replace: bool=True):
    os.system("mkdir -p " + c.FILES_DIR + "tmp/")
    if setup_type == "create":
        create_files(mapper=FileMapper(files_dir=os.path.join(c.FILES_DIR, "tmp", "")))
    elif setup_type == "api":
        load_files(mapper=FileMapper(files_dir=os.path.join(c.FILES_DIR, "tmp", "")))
    if replace:
        os.system("cp -r " + os.path.join(c.FILES_DIR, "tmp", "") + "* " + c.FILES_DIR)
        os.system("rm -rf " + os.path.join(c.FILES_DIR, "tmp", ""))


if __name__ == "__main__":
    desc = "     Run setup to create/load precalculated files."
    args = ru.save_parameters(script_desc=desc, arguments=('s', 'o'))
    main(setup_type=args.setup_type)
