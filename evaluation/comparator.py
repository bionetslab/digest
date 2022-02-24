#!/usr/bin/python3

from .d_utils import eval_utils as eu
from . import config as c
from .mappers import gene_getter as gg, disease_getter as dg
from .mappers.mapper import Mapper
from . import score_calculator as sc
from abc import abstractmethod


class Comparator:
    def __init__(self, mapper: Mapper, distance_measure: str, verbose: bool = False):
        self.mapper = mapper
        self.verbose = verbose
        self.distance_measure = distance_measure
        self.mapping = None
        self.id_set, self.id_type = None, None
        self.att_id, self.att_key, self.sparse_key = None, None, None

    def load_target(self, id_set, id_type):
        self.id_set = id_set
        self.id_type = id_type
        if id_type in c.SUPPORTED_DISEASE_IDS:
            self.mapping = dg.get_disease_to_attributes(disease_set=id_set, id_type=id_type, mapper=self.mapper)
            self.sparse_key, self.att_key, self.att_id = 'disease_mat_ids', 'disorder_atts', 'mondo'
        else:  # if id_type in c.SUPPORTED_GENE_IDS:
            self.mapping = gg.get_gene_to_attributes(gene_set=id_set, id_type=id_type, mapper=self.mapper)
            self.sparse_key, self.att_key, self.att_id = 'gene_mat_ids', 'gene_atts', 'entrezgene'

    @abstractmethod
    def compare(self, threshold: float = 0.0):
        pass


class SetComparator(Comparator):
    """
    Compare the set on itself based on connected attributes. See config for more info.
    """

    def compare(self, threshold: float = 0.0):
        result, mapped = dict(), dict()
        new_ids = self.mapper.update_distance_ids(in_series=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  key=self.sparse_key, distance_measure=self.distance_measure)
        for attribute in self.mapping.columns[1:]:
            subset_df = self.mapping[self.mapping[attribute].str.len() > 0]
            missing_values = len(self.mapping) - len(subset_df)
            if missing_values > 0:
                print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(
                    len(self.id_set))) if self.verbose else None
            if len(new_ids) > 0:
                comp_mat = eu.get_distance_matrix(full_att_series=self.mapper.loaded_mappings[self.att_key][attribute],
                                                  from_ids=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  id_to_index=self.mapper.loaded_distance_ids[self.sparse_key],
                                                  to_ids=new_ids)
                self.mapper.update_distances(in_mat=comp_mat, key=c.DISTANCES[attribute], id_type=self.sparse_key,
                                             distance_measure=self.distance_measure)
            if subset_df.empty:
                result[attribute] = 0
            else:
                ids = self.mapper.get_loaded_mapping_ids(in_ids=set(subset_df[subset_df.columns[0]]),
                                                         id_type=self.id_type)
                sub_mat = self.mapper.get_loaded_distances(in_series=ids[self.att_id], id_type=self.sparse_key,
                                                           key=c.DISTANCES[attribute],
                                                           distance_measure=self.distance_measure)
                missing_distances = ((len(self.mapping) * (len(self.mapping) - 1)) / 2) - sub_mat.getnnz()
                result[attribute] = ((sub_mat.getnnz()-sub_mat.sum()) + missing_distances) / \
                                    ((len(self.mapping) * (len(self.mapping) - 1)) / 2)
                mapped[attribute] = list(subset_df[c.ID_TYPE_KEY[self.id_type]])
        return result, mapped


class SetSetComparator(Comparator):
    """
    Compare two sets of the same type (either genes or diseases) with each other.
    The tar set is evaluated how good it matches to the ref set.
    """

    def __init__(self, mapper: Mapper, distance_measure: str, enriched: bool = False, verbose: bool = False):
        super().__init__(mapper=mapper, verbose=verbose, distance_measure=distance_measure)
        self.enriched = enriched
        self.ref_dict = None

    def load_reference(self, ref, ref_id_type):
        if ref_id_type in c.SUPPORTED_DISEASE_IDS:
            reference_mapping = dg.get_disease_to_attributes(disease_set=ref, id_type=ref_id_type, mapper=self.mapper)
        else:  # if ref_id_type in c.SUPPORTED_GENE_IDS:
            if self.enriched:
                reference_mapping = gg.get_enriched_attributes(gene_set=ref, id_type=ref_id_type, mapper=self.mapper)
            else:
                reference_mapping = gg.get_gene_to_attributes(gene_set=ref, id_type=ref_id_type, mapper=self.mapper)
        if self.enriched:
            self.ref_dict = eu.create_ref_dict(mapping=reference_mapping, keys=c.ENRICH_KEY.keys(), enriched=True)
        else:
            self.ref_dict = eu.create_ref_dict(mapping=reference_mapping, keys=reference_mapping.columns[1:],
                                               enriched=False)

    def compare(self, threshold: float = 0.0):
        return eu.evaluate_values(mapping=self.mapping, ref_dict=self.ref_dict, threshold=threshold,
                                  keys=self.mapping.columns[1:])


class IDSetComparator(Comparator):
    """
    Compare single id to either gene or disease set.
    The tar set is evaluated how good it matches to the ref id.
    """

    def __init__(self, mapper: Mapper, distance_measure: str, verbose: bool = False):
        super().__init__(mapper=mapper, verbose=verbose, distance_measure=distance_measure)
        self.ref_dict = None

    def load_reference(self, ref, ref_id_type, tar_id_type):
        if ref_id_type in c.SUPPORTED_DISEASE_IDS:
            id_mapping = dg.get_disease_to_attributes(disease_set={ref}, id_type=ref_id_type, mapper=self.mapper)
            if tar_id_type in c.SUPPORTED_DISEASE_IDS:
                self.ref_dict = eu.create_ref_dict(mapping=id_mapping, keys=id_mapping.columns[1:])
            else:  # if targets_id_type in c.SUPPORTED_GENE_IDS:
                id_mapping = id_mapping.rename(columns={'ctd.pathway_related_to_disease': 'pathway.kegg'})
                self.ref_dict = eu.create_ref_dict(mapping=id_mapping, keys={'pathway.kegg'})
        else:  # if targets_id_type in c.SUPPORTED_GENE_IDS:
            id_mapping = gg.get_gene_to_attributes(gene_set={ref}, id_type=ref_id_type, mapper=self.mapper)
            if tar_id_type in c.SUPPORTED_DISEASE_IDS:
                id_mapping = id_mapping.rename(columns={'pathway.kegg': 'ctd.pathway_related_to_disease'})
                self.ref_dict = eu.create_ref_dict(mapping=id_mapping, keys={'ctd.pathway_related_to_disease'})
            else:  # if targets_id_type in c.SUPPORTED_GENE_IDS:
                self.ref_dict = eu.create_ref_dict(mapping=id_mapping, keys=id_mapping.columns[1:])

    def compare(self, threshold: float = 0.0):
        return eu.evaluate_values(mapping=self.mapping, ref_dict=self.ref_dict, threshold=threshold,
                                  keys=self.ref_dict.keys())


class ClusterComparator(Comparator):
    """
    Evaluate the quality of clustering of given set with diseases or genes and
    assigned clusters. Additionally calculate statistical values with
    silhouette score and dunn index.
    """

    def __init__(self, mapper: Mapper, distance_measure: str, verbose: bool = False):
        super().__init__(mapper=mapper, verbose=verbose, distance_measure=distance_measure)
        self.clustering = None

    def load_target(self, id_set, id_type):
        id_set['cluster_index'] = id_set.groupby('cluster').ngroup()
        super().load_target(id_set=id_set["id"], id_type=id_type)
        self.clustering = id_set[['id', 'cluster', 'cluster_index']]

    def compare(self, threshold: float = 0.0):
        result_di, result_ss, result_ss_intermediate, result_dbi, mapped = dict(), dict(), dict(), dict(), dict()
        new_ids = self.mapper.update_distance_ids(in_series=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  key=self.sparse_key, distance_measure=self.distance_measure)
        for attribute in self.mapping.columns[1:]:
            subset_df = self.mapping[self.mapping[attribute].str.len() > 0]
            subset_clusters = self.clustering[self.clustering['id'].isin(subset_df[c.ID_TYPE_KEY[self.id_type]])][
                ['id', 'cluster_index']]
            missing_values = len(self.mapping) - len(subset_df)
            if missing_values > 0:
                print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(
                    len(self.mapping)) + "") if self.verbose else None

            if len(new_ids) > 0:
                comp_mat = eu.get_distance_matrix(full_att_series=self.mapper.loaded_mappings[self.att_key][attribute],
                                                  from_ids=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  id_to_index=self.mapper.loaded_distance_ids[self.distance_measure][
                                                      self.sparse_key],
                                                  to_ids=new_ids)
                self.mapper.update_distances(in_mat=comp_mat, key=c.DISTANCES[attribute], id_type=self.sparse_key,
                                             distance_measure=self.distance_measure)

            if subset_df.empty:
                result_di[attribute], result_ss[attribute], result_ss_intermediate[attribute] = None, None, None
                result_dbi[attribute], mapped[attribute] = None, []
            else:
                ids = self.mapper.get_loaded_mapping_ids(in_ids=set(subset_df[subset_df.columns[0]]),
                                                         id_type=self.id_type)
                distances = self.mapper.get_loaded_distances(in_series=ids[self.att_id], id_type=self.sparse_key,
                                                             key=c.DISTANCES[attribute],
                                                             distance_measure=self.distance_measure)
                distances = dict(distances.todok().items())
                inv_index_to_id = {index: value for index, value in ids[self.att_id].reset_index(drop=True).items()}
                precalc_dist = sc.precalc_distance_dicts(ids_cluster=subset_clusters, ids_mapping=ids,
                                                         distances=distances,
                                                         index_to_id=inv_index_to_id,
                                                         ids={'id_type': c.ID_TYPE_KEY[self.id_type],
                                                              'sparse_key': self.sparse_key, 'att_id': self.att_id,
                                                              'attribute': c.DISTANCES[attribute]})
                ss_score = sc.silhouette_score(ids_cluster=subset_clusters, distances=precalc_dist, linkage="average")
                di_score = sc.dunn_index(ids_cluster=subset_clusters, distances=precalc_dist, linkage="average")
                dbi_score = sc.davies_bouldin_index(ids_cluster=subset_clusters, distances=precalc_dist,
                                                    linkage="average")
                result_di[attribute] = di_score
                result_ss[attribute] = ss_score[0]
                result_ss_intermediate[attribute] = ss_score[1]
                result_dbi[attribute] = dbi_score
                mapped[attribute] = list(subset_df[c.ID_TYPE_KEY[self.id_type]])
        return result_di, result_ss, result_dbi, result_ss_intermediate, mapped
